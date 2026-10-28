#!/usr/bin/env Rscript
rm(list = ls())

suppressPackageStartupMessages({
  library(yaml)
  library(Matrix)
  library(edgeR)
  library(arrow)
  library(data.table)
  library(fs)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: 03_pseudobulk.R <config.yaml> <cellset_id>")

config_path <- args[1]
cellset_id  <- args[2]

cat("=== S03_build_pseudobulk ===\n")
cat("Cell set:", cellset_id, "\n")

config <- yaml::read_yaml(config_path)
base_dir <- config$input$base_dir

pb_cfg <- config$pseudobulk

unit_col      <- pb_cfg$aggregation$unit_col
condition_col <- pb_cfg$aggregation$condition_col
batch_col     <- pb_cfg$aggregation$batch_col
donor_col     <- pb_cfg$aggregation$donor_col

min_cells_per_unit        <- pb_cfg$filtering$units$min_cells_per_unit
min_total_counts_per_unit <- pb_cfg$filtering$units$min_total_counts_per_unit
use_group_for_filter      <- pb_cfg$filtering$genes$use_group_for_filter

downsample_enabled <- pb_cfg$downsampling$enabled
downsample_max     <- pb_cfg$downsampling$max_cells_per_unit
downsample_seed    <- pb_cfg$downsampling$seed

set.seed(downsample_seed)

merged_dir      <- config$input$merged
cellset_dir     <- file.path(base_dir, "outputs", "cell_sets", cellset_id)
out_dir         <- file.path(base_dir, "outputs", "pseudobulk", cellset_id)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cells_path      <- file.path(merged_dir, "cells.parquet")
genes_path      <- file.path(merged_dir, "genes.tsv")
membership_path <- file.path(cellset_dir, "membership.parquet")

cat("Loading merged metadata and membership...\n")
merged_cells <- as.data.table(read_parquet(cells_path))
membership <- as.data.table(read_parquet(membership_path))
membership_ids <- membership$cell_id

required_cols <- c("cell_id", unit_col, donor_col, condition_col, batch_col)
missing_cols <- setdiff(required_cols, colnames(merged_cells))
if (length(missing_cols) > 0) {
  stop("Missing required columns in merged cells: ", paste(missing_cols, collapse = ", "))
}

cat("Loading genes...\n")
genes <- fread(genes_path)
gene_col <- if ("gene_id" %in% names(genes)) "gene_id" else names(genes)[1]

data_dir <- config$input$data_dir
runs <- dir_ls(data_dir, type = "directory")

pb_counts_list <- list()
unit_list <- list()

aggregate_one_run <- function(run_path) {
  run_name <- path_file(run_path)
  s01_dir <- file.path("/mnt/odap-beegfs/a008/DGE/data/ingested", run_name, "S01_ingest")
  
  if (!dir_exists(s01_dir)) {
    warning("Skipping run without S01_ingest: ", run_name)
    return(NULL)
  }
  
  cat("Processing run:", run_name, "\n")
  
  run_cells <- as.data.table(read_parquet(file.path(s01_dir, "cells.parquet")))
  run_genes <- fread(file.path(s01_dir, "genes.tsv"))
  
  if (!identical(run_genes[[gene_col]], genes[[gene_col]])) {
    stop("Gene order/content differs in run: ", run_name)
  }
  
  missing <- setdiff(required_cols, names(run_cells))
  if (length(missing)) {
    stop("Run ", run_name, " missing columns: ", paste(missing, collapse = ", "))
  }
  
  run_cells[, original_cell_id := cell_id]
  run_cells[, source_run := run_name]
  run_cells[, cell_id := paste0(run_name, "_", cell_id)]
  
  keep <- run_cells$cell_id %in% membership_ids
  if (!any(keep)) {
    cat("  No cells for this cell set in run:", run_name, "\n")
    return(NULL)
  }
  
  counts <- Matrix::readMM(file.path(s01_dir, "counts.mtx.gz"))
  counts <- as(counts, "dgCMatrix")
  
  if (ncol(counts) != nrow(run_cells)) {
    stop("Counts columns do not match cells rows for run: ", run_name)
  }
  
  counts <- counts[, keep, drop = FALSE]
  run_cells <- run_cells[keep]
  
  if (downsample_enabled) {
    cat("  Applying downsampling...\n")
    keep_cells <- unlist(lapply(split(seq_len(nrow(run_cells)), run_cells[[unit_col]]), function(idx) {
      if (length(idx) > downsample_max) sample(idx, downsample_max) else idx
    }), use.names = FALSE)
    
    run_cells <- run_cells[keep_cells]
    counts <- counts[, keep_cells, drop = FALSE]
  }
  
  if (ncol(counts) == 0) return(NULL)
  
  group <- factor(run_cells[[unit_col]])
  group_levels <- levels(group)
  
  design_sparse <- sparseMatrix(
    i = seq_along(group),
    j = as.integer(group),
    x = 1,
    dims = c(length(group), length(group_levels)),
    dimnames = list(NULL, group_levels)
  )
  
  pb_counts <- counts %*% design_sparse
  pb_counts <- as(pb_counts, "dgCMatrix")
  
  meta_cols <- unique(c(unit_col, donor_col, condition_col, batch_col))
  unit_df <- run_cells[, .SD[1], by = unit_col, .SDcols = setdiff(meta_cols, unit_col)]
  unit_df <- unit_df[match(group_levels, unit_df[[unit_col]])]
  
  unit_df[, n_cells := as.integer(table(group)[group_levels])]
  unit_df[, library_size := as.numeric(Matrix::colSums(pb_counts))]
  unit_df[, source_run := run_name]
  
  # Make matrix column names unique even if sample_id values are reused across plates
  unit_df[, pseudobulk_id := paste(source_run, get(unit_col), sep = "_")]
  colnames(pb_counts) <- unit_df$pseudobulk_id
  
  list(counts = pb_counts, units = unit_df)
}

for (run in runs) {
  res <- aggregate_one_run(run)
  if (!is.null(res)) {
    pb_counts_list[[length(pb_counts_list) + 1]] <- res$counts
    unit_list[[length(unit_list) + 1]] <- res$units
  }
}

if (length(pb_counts_list) == 0) {
  stop("No pseudobulk counts generated for cell set: ", cellset_id)
}

cat("Combining pseudobulk matrices across runs...\n")
pb_counts <- pb_counts_list[[1]]
if (length(pb_counts_list) > 1) {
  for (i in 2:length(pb_counts_list)) {
    pb_counts <- Matrix::cbind2(pb_counts, pb_counts_list[[i]])
    pb_counts <- as(pb_counts, "dgCMatrix")
  }
}

unit_df <- rbindlist(unit_list, use.names = TRUE, fill = TRUE)

if (ncol(pb_counts) != nrow(unit_df)) {
  stop("Combined pseudobulk counts columns do not match unit metadata rows")
}

cat("Filtering units...\n")
keep_units <- unit_df$n_cells >= min_cells_per_unit &
  unit_df$library_size >= min_total_counts_per_unit

if (sum(keep_units) == 0) stop("All units filtered out for cell set: ", cellset_id)

pb_counts <- pb_counts[, keep_units, drop = FALSE]
unit_df <- unit_df[keep_units]

cat("Filtering genes with edgeR::filterByExpr...\n")
dge <- DGEList(counts = pb_counts)

if (use_group_for_filter) {
  group <- factor(unit_df[[condition_col]])
  keep_genes <- filterByExpr(dge, group = group)
} else {
  keep_genes <- filterByExpr(dge)
}

pb_counts <- pb_counts[keep_genes, , drop = FALSE]
genes_used <- genes[keep_genes]

qc_df <- unit_df[, .SD, .SDcols = intersect(
  c("pseudobulk_id", unit_col, donor_col, condition_col, batch_col,
    "source_run", "n_cells", "library_size"),
  names(unit_df)
)]

cat("Writing outputs...\n")
tmp <- tempfile(fileext = ".mtx")
pb_counts <- as(pb_counts, "dgCMatrix")
Matrix::writeMM(pb_counts, tmp)
system2("gzip", c("-c", tmp), stdout = file.path(out_dir, "counts.mtx.gz"))
unlink(tmp)

fwrite(unit_df, file.path(out_dir, "units.tsv"), sep = "\t")
fwrite(qc_df, file.path(out_dir, "qc.tsv"), sep = "\t")
fwrite(genes_used, file.path(out_dir, "genes_used.tsv"), sep = "\t")

cat("S03 complete for:", cellset_id, "\n")

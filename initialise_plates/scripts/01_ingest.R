#!/usr/bin/env Rscript
rm(list = ls())
# ----------------------------
# Packages
# ----------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(arrow)
  library(Matrix)
  library(yaml)
  library(SeuratObject)
  library(Seurat)
})

# ----------------------------
# Arguments
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: 01_ingest.R <config.yaml> <sample_name>")
}

config <- yaml::read_yaml(args[1])
sample_name <- args[2]

# ----------------------------
# Config & paths
# ----------------------------
base <- config$input$base_dir
if (is.null(base)) stop("config$input$base_dir is required")

data_dir<-config$input$data_dir

ingest_cfg <- config$ingest
if (is.null(ingest_cfg$default)) stop("ingest$default section is required")

# Input folder for this sample
input_dir <- file.path(data_dir, "raw", sample_name)
if (!dir.exists(input_dir)) stop("Input directory does not exist: ", input_dir)

# Pick mapping: per-sample override or default
mapping <- ingest_cfg$samples[[sample_name]]
if (is.null(mapping)) mapping <- ingest_cfg$default

required_fields <- c("sample_id_col", "donor_id_col", "condition_col",
                     "batch_col", "raw_cell_id_cols")
missing <- setdiff(required_fields, names(mapping))
if (length(missing)) stop("Missing ingest mapping fields: ", paste(missing, collapse=", "))

# Output dir
out_dir <- file.path(data_dir, "ingested", sample_name, "S01_ingest")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Base directory: ", base)
message("Input directory: ", input_dir)
message("Sample name: ", sample_name)
message("Output directory: ", out_dir)

# ----------------------------
# Source helper functions
# ----------------------------
rfunc_dir <- file.path(base, "scripts", "R_functions")
source(file.path(rfunc_dir, "detect_input_type.R"))
source(file.path(rfunc_dir, "read_parse_inputs.R"))
source(file.path(rfunc_dir, "read_10x_inputs.R"))
source(file.path(rfunc_dir, "read_zarr_inputs.R"))
source(file.path(rfunc_dir, "read_rds_inputs.R"))
source(file.path(rfunc_dir, "read_h5ad_inputs.R"))

# ----------------------------
# Harmonisation
# ----------------------------
harmonise_cells <- function(cells, mapping, raw_cell_id_col) {
  setDT(cells)
  
  required <- c(mapping$sample_id_col,
                mapping$donor_id_col,
                mapping$condition_col,
                mapping$batch_col)
  missing <- setdiff(required, colnames(cells))
  if (length(missing)) stop("Missing metadata columns: ", paste(missing, collapse=", "))
  
  # Resolve raw cell IDs
  if (raw_cell_id_col == "rownames") {
    cells[, .raw_cell_id := cell_id]
  } else {
    if (!raw_cell_id_col %in% colnames(cells)) {
      stop("Raw cell ID column not found: ", raw_cell_id_col)
    }
    cells[, .raw_cell_id := get(raw_cell_id_col)]
  }
  
  # Construct cell_id
  if(raw_cell_id_col == "rownames" &&
     "cell_id" %in% names(cells)) {
    
    cells[, cell_id := .raw_cell_id]
    } else {
      cells[, cell_id := paste(
        get(mapping$sample_id_col),
        .raw_cell_id,
        sep = "_"
      )]
      }
  # Harmonise standard columns
  cells[, `:=`(
    sample_id = get(mapping$sample_id_col),
    donor_id  = get(mapping$donor_id_col),
    condition = get(mapping$condition_col),
    batch     = get(mapping$batch_col)
  )]
  
  if (anyDuplicated(cells$cell_id)) stop("cell_id is not unique after harmonisation")
  
  cells[, .raw_cell_id := NULL]
  cells[]
}

# ----------------------------
# Validation
# ----------------------------
validate_inputs <- function(counts, genes, cells) {
  if (!inherits(counts, "dgCMatrix")) stop("Counts must be a dgCMatrix")
  if (ncol(counts) != nrow(cells)) stop("Counts columns != number of cells")
  if (nrow(counts) != nrow(genes)) stop("Counts rows != number of genes")
  if (anyDuplicated(cells$cell_id)) stop("cell_id not unique")
  invisible(TRUE)
}

# ----------------------------
# Detect & read input
# ----------------------------
input_type <- detect_input_type(input_dir)
message("Detected input type: ", input_type)

dat <- switch(
  input_type,
  parse_splitpipe = read_parse_inputs(input_dir),
  standard_matrix = read_10x_inputs(input_dir),
  rds             = read_rds_input(input_dir),
  h5ad            = read_h5ad_input(input_dir),
  zarr            = read_zarr_input(input_dir),
  stop("Unsupported input type: ", input_type)
)

# ----------------------------
# Extract data
# ----------------------------

counts<-dat$counts
count_mats<-dat$count_mats
cells<-dat$cells
genes<-dat$genes

if(is.null(count_mats)) {
  count_mats <- list(counts = counts)
}

# ----------------------------
# Harmonise metadata
# ----------------------------

raw_cell_id_col <- mapping$raw_cell_id_cols[[input_type]]

if(is.null(raw_cell_id_col)) {
  stop("No raw_cell_id_col defined for: ", input_type)
}

cells<-harmonise_cells(
  cells,
  mapping,
  raw_cell_id_col
)

# ----------------------------
# Match each matrix to harmonised metadata
# ----------------------------

raw_ids<-unlist(
  lapply(count_mats, colnames),
  use.names = F
)

if(anyDuplicated(raw_ids)) {
  stop("Duplicate raw cell IDs across matrices")
}

if(raw_cell_id_col == "rownames") {
  metadata_raw_ids <- as.character(dat$cells$cell_id)
} else {
  metadata_raw_ids <- as.character(
    dat$cells[[raw_cell_id_col]]
  )
}

if(anyDuplicated(metadata_raw_ids)) {
  stop("Duplicate raw cell IDs in metadata")
}

if(length(raw_ids) != length(metadata_raw_ids)) {
  stop(
    "Count matrices countain ", length(raw_ids),
    " cells, but metadata contains ",
    length(metadata_raw_ids), " rows"
  )
}

match_idx<-match(raw_ids, metadata_raw_ids)

if(anyNA(match_idx)) {
  stop("Some count matric cell IDs are missing from metadata")
}

cells<-cells[match_idx, ]

if(!identical(
  as.character(raw_ids),
  as.character(metadata_raw_ids[match_idx])
)) {
  stop("Count matrix and metadata alignment failed")
}

message("Successfully aligned ", length(raw_ids), " cells with count matrices")
# ----------------------------
# Rename matrix columns and calculate UMI totals
# ----------------------------

cell_ids <- cells$cell_id

offset <- 0L
umi_totals <- numeric(length(cell_ids))

for(j in seq_along(count_mats)) {
  m<-count_mats[[j]]
  
  idx<-seq.int(
    offset + 1L,
    length.out = ncol(m)
  )
  
  colnames(m) <- cell_ids[idx]
  
  umi_totals[idx] <- Matrix::colSums(m)
  
  count_mats[[j]] <- m
  
  offset <- offset + ncol(m)
}

if(offset != length(cell_ids)) {
  stop("Not all cells were assigned to count matrices")
}

cells[, n_umis := umi_totals]

# ----------------------------
# Validate
# ----------------------------

if(anyDuplicated(cells$cell_id)) {
  stop("Duplicate harmonised cell IDs")
}

if(anyDuplicated(genes$gene_id)) {
  stop("Duplicated gene IDs")
}

for (m in count_mats) {
  
  if(!inherits(m, "dgCMatrix")) {
    stop("Counts must be dgCMatrix")
  }
  
  if(nrow(m) != nrow(genes)) {
    stop("Counts/genes dimension mismatch")
  }
  
  if(!identical(rownames(m), as.character(genes$gene_id))) {
    stop("Counts/genes ordering mismatch")
  }
}



matrix_cells<-unlist(
  lapply(count_mats, colnames),
  use.names = F
)

if(!identical(matrix_cells, cells$cell_id)) {
  stop("Matrix column order does not match cells.parquet")
}

# ----------------------------
# Write outputs
# ----------------------------

dir.create(
  out_dir,
  recursive = T,
  showWarnings = F
)

fwrite(
  genes,
  file.path(out_dir, "genes.tsv"),
  sep = "\t"
)

write_parquet(
  cells,
  file.path(out_dir, "cells.parquet")
)

manifest <- vector(
  "list",
  length(count_mats)
)

cell_offset<-0L

for(j in seq_along(count_mats)) {
  m<-count_mats[[j]]
  
  filename<-paste0(
    "counts_",
    j,
    ".mtx.gz"
  )
  
  output_path<-file.path(
    out_dir,
    filename
  )
  
  temporary_path<-sub(
    "\\.gz$",
    "",
    output_path
  )
  
  Matrix::writeMM(
    m,
    temporary_path
  )
  
  status<-system2(
    "gzip",
    c("-f", shQuote(temporary_path))
  )
  
  if(status != 0){
    stop("Failed to compress: ", temporary_path)
  }
  
  cell_start<-cell_offset + 1L
  cell_end <- cell_offset + ncol(m)
  
  manifest[[j]]<-data.table(
    layer = names(count_mats)[j],
    counts_file = output_path,
    n_genes = nrow(m),
    n_cells = ncol(m),
    cell_start = cell_start,
    cell_end = cell_end
  )
  
  cell_offset <- cell_end
  
  message(
    "Written layer ",
    j,
    "/",
    length(count_mats),
    ": ",
    filename
  )
}

if(cell_offset != nrow(cells)) {
  stop("Manifest cell ranges do not matcdh metadata")
}

manifest<-rbindlist(manifest)

fwrite(
  manifest,
  file.path(out_dir, "counts_manifest.tsv"),
  sep = "\t"
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

message(
  "S01_ingest completed successfully for sample: ",
  sample_name
)
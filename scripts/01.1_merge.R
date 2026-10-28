#!/usr/bin/env Rscript
rm(list = ls())

library(data.table)
library(arrow)
library(fs)
library(yaml)

# ----------------------------
# Arguments
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)


# ----------------------------
# Config
# ----------------------------

config <- yaml::read_yaml(args[1])
ingest_cfg <- config$ingest
base<-config[["input"]][["base_dir"]]

required_cols <- c("cell_id", "donor_id", "condition", "batch")

# ----------------------------
# Sample directories
# ----------------------------
data_dir <- file.path(base, "data")
samples <- dir_ls(data_dir, type = "directory")

out_dir <- file.path(base, "outputs", "merged")
dir_create(out_dir, recurse = TRUE)

all_cells <- list()
all_genes <- list()
all_counts_files <- list()

read_mtx_dims <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  
  repeat {
    line <- readLines(con, n = 1)
    
    if (length(line) == 0) {
      stop("Could not find dimensions in Matrix Market file: ", path)
    }
    
    if (!startsWith(line, "%")) {
      break
    }
  }
  
  dims <- as.numeric(strsplit(trimws(line), "\\s+")[[1]])
  
  c(n_genes = dims[1], n_cells = dims[2])
}

# ----------------------------
# Read S01 outputs per sample
# ----------------------------
for (sample in samples) {
  sample_name <- path_file(sample)
  s01_dir <- file.path(base, "outputs", sample_name, "S01_ingest")
  
  if (!dir_exists(s01_dir)) {
    warning("S01_ingest folder not found for sample: ", sample_name)
    next
  }
  
  # Read outputs
  cells <- read_parquet(file.path(s01_dir, "cells.parquet"))
  genes <- fread(file.path(s01_dir, "genes.tsv"))
  count_file <- file.path(s01_dir, "counts.mtx.gz")
  count_dims <- read_mtx_dims(count_file)
  
  # Prefix cell_id with sample to ensure uniqueness
  cells$cell_id <- paste0(sample_name, "_", cells$cell_id)
  cells$source_dataset <- sample_name

  # Check for missing required columns
  missing_cols <- setdiff(required_cols, colnames(cells))
  if (length(missing_cols)) {
    stop("Sample ", sample_name, " is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  
  if (count_dims["n_cells"] != nrow(cells)) {
    stop(
      "Sample ", sample_name,
      ": count matrix has ", count_dims["n_cells"],
      " cells but cells.parquet has ", nrow(cells)
    )
  }
  
  if (count_dims["n_genes"] != nrow(genes)) {
    stop(
      "Sample ", sample_name,
      ": count matrix has ", count_dims["n_genes"],
      " genes but genes.tsv has ", nrow(genes)
    )
  }
  # Track all extra columns
  extra_cols <- setdiff(colnames(cells), required_cols)
  if (length(extra_cols)) {
    message("Sample ", sample_name, " has extra metadata columns: ",
            paste(extra_cols, collapse = ", "))
  }
  
  all_cells[[sample_name]] <- as.data.table(cells)
  all_genes[[sample_name]] <- genes
  all_counts_files[[sample_name]] <- data.table(
    source_dataset = sample_name,
    counts_file = count_file,
    n_genes = unname(count_dims["n_genes"]),
    n_cells = unname(count_dims["n_cells"])
  )
}

# ----------------------------
# Merge outputs
# ----------------------------
# Cells
merged_cells <- rbindlist(all_cells, use.names = TRUE, fill = TRUE)

# Genes
reference_genes <- all_genes[[1]]

for (nm in names(all_genes)[-1]) {
  
  if (!identical(
    reference_genes$gene_id,
    all_genes[[nm]]$gene_id
  )) {
    stop(
      "Gene IDs or gene order differ between datasets: ",
      names(all_genes)[1],
      " and ",
      nm
    )
  }
}

merged_genes <- reference_genes

count_manifest <- rbindlist(all_counts_files)

# ----------------------------
# Validation
# ----------------------------
if (!all(required_cols %in% colnames(merged_cells))) {
  stop("Merged cells table is missing required columns")
}

if (sum(count_manifest$n_cells) != nrow(merged_cells)) {
  stop(
    "Total number of cells in count matrices does not match ",
    "number of rows in merged metadata"
  )
}

# ----------------------------
# Save merged outputs
# ----------------------------
fwrite(merged_genes, file.path(out_dir, "genes.tsv"), sep = "\t")
write_parquet(merged_cells, file.path(out_dir, "cells.parquet"))

fwrite(
  count_manifest,
  file.path(out_dir, "counts_manifest.tsv"),
  sep = "\t"
)

message("Merged S01 outputs saved to: ", out_dir)


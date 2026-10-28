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
data<-config$input$data_dir

required_cols <- c("cell_id", "donor_id", "condition", "batch")

# ----------------------------
# Sample directories
# ----------------------------
data_dir <- file.path(data, "ingested")
samples <- dir_ls(data_dir, type = "directory")

out_dir <- file.path(data, "merged")
dir_create(out_dir, recursive = TRUE)

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
  
  s01_dir <- file.path(
    data,
    "ingested",
    sample_name,
    "S01_ingest"
  )
  
  if(!dir_exists(s01_dir)) {
    warning("S01_ingest folder not found: ", sample_name)
    next
  }
  
  cells<-as.data.table(
    read_parquet(file.path(s01_dir, "cells.parquet"))
  )
  
  genes<-fread(
    file.path(s01_dir, "genes.tsv")
  )
  
  manifest<-fread(
    file.path(s01_dir, "counts_manifest.tsv")
  )
  
  required_manifest_cols <-c(
    "counts_file",
    "n_genes",
    "n_cells",
    "cell_start",
    "cell_end"
  )
  
  setorder(manifest, cell_start)
  
  if(!identical(
    as.numeric(manifest$cell_start),
    c(1, head(manifest$cell_end, -1) + 1)
  )) {
    stop("Non-contiguous cell ranges: ", sample_name)
  }
  
  if(tail(manifest$cell_end, 1) != nrow(cells)) {
    stop("Manifest does not cover all cells: ", sample_name)
  }
  
  if(any(
    manifest$cell_end - manifest$cell_start + 1 != manifest$n_cells
  )) {
    stop("Incorrect cell ranges: ", sample_name)
  }
  
  cells[, cell_id := paste0(sample_name, "_", cell_id)]
  cells[, source_dataset := sample_name]
  
  missing_cols<-setdiff(required_cols, names(cells))
  
  if(length(missing_cols)) {
    stop(
      "Missing metadata columns for ", sample_name, ": ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  sample_counts<-vector("list", nrow(manifest))
  
  for(j in seq_len(nrow(manifest))) {
    count_file <- manifest$counts_file[j]
    
    if(!file.exists(count_file)) {
      stop("Count matrix not found: ", count_file)
    }
    
    dims<-read_mtx_dims(count_file)
    
    if(
      dims["n_genes"] != manifest$n_genes[j] ||
      dims["n_cells"] != manifest$n_cells[j]
    ) {
      stop("Matrix dime stions differ from manifest: ", count_file)
    }
  
    if(dims["n_genes"] != nrow(genes)) {
      stop("Gene count mismatch: ", count_file)
    }
    
    sample_counts[[j]]<-data.table(
      source_dataset = sample_name,
      counts_file = count_file,
      n_genes = unname(dims["n_genes"]),
      n_cells = unname(dims["n_cells"]),
      cell_start = manifest$cell_start[j],
      cell_end = manifest$cell_end[j]
    )
  }
  
  all_cells[[sample_name]]<-cells
  all_genes[[sample_name]]<-genes
  all_counts_files[[sample_name]]<-rbindlist(sample_counts)
  
  extra_cols<-setdiff(names(cells), required_cols)
  
  if(length(extra_cols)) {
    message("Sample ", sample_name,
            " has extra metadata columns: ",
            paste(extra_cols, collapse = ", "))
  }
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

offset<-0

for(nm in names(all_counts_files)) {
  
  all_counts_files[[nm]][
    ,
    `:=`(
      cell_start = cell_start + offset,
      cell_end = cell_end + offset
    )
  ]
  
  offset<-offset + nrow(all_cells[[nm]])
}

count_manifest <- rbindlist(all_counts_files)

if(anyDuplicated(merged_cells$cell_id)) {
  stop("Duplicated cell IDs in merged metadata")
}
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


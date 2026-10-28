read_rds_input <- function(input_dir, assay = "RNA") {
  library(Matrix)
  library(data.table)
  
  rds_files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0)
    stop("No .rds file found in ", input_dir)
  if (length(rds_files) > 1)
    stop("Multiple .rds files found: ", paste(rds_files, collapse = ", "))
  
  obj <- readRDS(rds_files[1])
  
  # Optional user-provided metadata
  metadata_file <- list.files(input_dir, pattern = "_metadata\\.(csv|tsv)$", full.names = TRUE)
  use_metadata <- length(metadata_file) == 1
  
  counts<-NULL
  count_mats<-NULL
  cells<-NULL
  genes<-NULL
  
  # -------------------------
  # SingleCellExperiment
  # -------------------------
  if (inherits(obj, "SingleCellExperiment")) {
    message("Detected SingleCellExperiment")
    
    counts <- SingleCellExperiment::counts(obj)
    counts <- as(counts, "dgCMatrix")
    
    # Ensure genes x cells
    if (nrow(counts) < ncol(counts)) {
      message("Counts appears transposed. Fixing orientation.")
      counts <- t(counts)
    }
    
    # Cells
    if (use_metadata) {
      cells <- fread(metadata_file[1])
      message("Using user-provided metadata: ", basename(metadata_file[1]))
    } else {
      cells <- as.data.table(SummarizedExperiment::colData(obj))
    }
    cells[, cell_id := colnames(counts)]
    
    # Genes
    genes <- as.data.table(SummarizedExperiment::rowData(obj))
    genes[, gene_id := rownames(counts)]
  }
  
  # -------------------------
  # Seurat
  # -------------------------
  else if (inherits(obj, "Seurat")) {
    message("Detected Seurat object")
    
    if (!assay %in% names(obj@assays))
      stop("Assay '", assay, "' not found in Seurat object")
    
    assay_obj<-obj[[assay]]
    
    count_layers<-SeuratObject::Layers(
      assay_obj,
      search = "^counts"
    )
    
    if(!length(count_layers)) {
      stop("No count layers found in assay: ", assay)
    }
    
    message("Found ", length(count_layers), " count layers")
      
      count_mats<-lapply(
        count_layers,
        function(layer_name) {
          m<-SeuratObject::LayerData(
            assay_obj,
            layer = layer_name
          )
          
          as(m, "dgCMatrix")
        }
      )
      
      names(count_mats)<-count_layers
      
      reference_genes<-rownames(count_mats[[1]])
      
      for (m in count_mats) {
        if(!identical(rownames(m), reference_genes)) {
          stop("Count laters have inconsistent genes/order")
        }
      }
      
      all_cells<-unlist(lapply(count_mats, colnames), use.names = F)
      
      if(anyDuplicated(all_cells)) {
        stop("Duplicates cell IDs across count layers")
      }
    
    
    # Cells
    if (use_metadata) {
      cells <- fread(metadata_file[1])
      message("Using user-provided metadata: ", basename(metadata_file[1]))
      
      if(!"cell_id"%in% names(cells)) {
        stop("External metadata must contain cell_id")
      }
      
      if(anyDuplicated(cells$cell_id)) {
        stop("Duplicate cell IDs in metadata")
      }
      
      if(!setequal(Cells$cell_id, all_cells)) {
        stop("External metadata cells do not match count laters")
      }
      
      cells<-cells[
        match(all_cells, cell_id)
      ]
    } else {
      cells <- as.data.table(obj@meta.data)
      
      cells[, cell_id := rownames(obj@meta.data)]
      
      if(!setequal(cells$cell_id, all_cells)) {
        stop("Seurat metadata cells do not match count laters")
      }
      
      cells<-cells[
        match(all_cells, cell_id)
      ]
    }
    
    # Genes
    genes <- data.table(
      gene_id = reference_genes
    )
    
    } else {
    stop("Unsupported RDS object class: ", paste(class(obj), collapse = ", "))
  }
  
  # -------------------------
  # Final sanity checks
  # -------------------------
  
  if(!is.null(count_mats)) {
    total_cells <- sum(
      vapply(count_mats, ncol, integer(1))
    )
    
    if(total_cells != nrow(cells)) {
      stop("Layer cell count does not match metadata")
    }
    
    if(anyDuplicated(genes$gene_id)) {
      stop("Duplicate gene names")
    }
  } else {
    if(ncol(counts) != nrow(cells)) {
      stop("Counts/metadata cell count mismatch")
    }
    
    if(nrow(counts) != nrow(genes)) {
      stop("Count/genes mismatch")
    }
    
    if(anyDuplicated(colnames(counts))) {
      stop("Duplicate cell names")
    }
    
    if(anyDuplicated(rownames(counts))) {
      stop("Duplicate gene names")
    }
  }
  
  # -------------------------
  # Return
  # -------------------------
  
  list(
    counts = counts,   
    count_mats = count_mats,
    cells  = cells,
    genes  = genes
  )
}



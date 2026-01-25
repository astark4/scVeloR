#' @title Data Conversion Functions
#' @description Functions for converting data between Seurat objects and 
#' scVeloR internal formats.
#' @name data_conversion
NULL

#' Import Velocity Data from Seurat
#'
#' @description Extract spliced and unspliced counts from Seurat object.
#' Supports both Seurat V4 and V5 object structures.
#'
#' @param object Seurat object
#' @param spliced_layer Name of spliced count layer
#' @param unspliced_layer Name of unspliced count layer
#'
#' @return List with spliced and unspliced matrices
#' @export
import_velocity_data <- function(object,
                                  spliced_layer = "spliced",
                                  unspliced_layer = "unspliced") {
  
  # Detect Seurat version
  version <- seurat_version(object)
  
  if (version == 5) {
    # Seurat V5
    spliced <- tryCatch({
      Seurat::LayerData(object, layer = spliced_layer)
    }, error = function(e) NULL)
    
    unspliced <- tryCatch({
      Seurat::LayerData(object, layer = unspliced_layer)
    }, error = function(e) NULL)
  } else {
    # Seurat V4
    assay <- Seurat::DefaultAssay(object)
    
    spliced <- tryCatch({
      object@assays[[assay]]@layers[[spliced_layer]]
    }, error = function(e) {
      tryCatch({
        Seurat::GetAssayData(object, slot = spliced_layer)
      }, error = function(e2) NULL)
    })
    
    unspliced <- tryCatch({
      object@assays[[assay]]@layers[[unspliced_layer]]
    }, error = function(e) {
      tryCatch({
        Seurat::GetAssayData(object, slot = unspliced_layer)
      }, error = function(e2) NULL)
    })
  }
  
  if (is.null(spliced) || is.null(unspliced)) {
    stop("Could not find spliced/unspliced layers. Ensure they are present in the Seurat object.")
  }
  
  # Convert to dense matrices (cells x genes)
  spliced <- as.matrix(spliced)
  unspliced <- as.matrix(unspliced)
  
  # Transpose if genes are in rows
  if (nrow(spliced) != ncol(object)) {
    spliced <- t(spliced)
    unspliced <- t(unspliced)
  }
  
  list(
    spliced = spliced,
    unspliced = unspliced
  )
}

#' Create Seurat Object with Velocity Layers
#'
#' @description Create a Seurat object with spliced and unspliced layers
#' from count matrices.
#'
#' @param spliced Spliced count matrix (cells x genes or genes x cells)
#' @param unspliced Unspliced count matrix
#' @param meta.data Cell metadata (optional)
#' @param project Project name
#'
#' @return Seurat object
#' @export
create_velocity_seurat <- function(spliced,
                                    unspliced,
                                    meta.data = NULL,
                                    project = "scVeloR") {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required")
  }
  
  # Ensure matrices are in genes x cells format
  if (nrow(spliced) == nrow(unspliced) && ncol(spliced) == ncol(unspliced)) {
    # Assume genes x cells if more rows than columns
    if (nrow(spliced) > ncol(spliced)) {
      # Already genes x cells
    } else {
      # Transpose to genes x cells
      spliced <- t(spliced)
      unspliced <- t(unspliced)
    }
  }
  
  # Create Seurat object with spliced as main counts
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = spliced,
    project = project,
    meta.data = meta.data
  )
  
  # Add unspliced layer
  seurat_obj[["RNA"]]@layers$unspliced <- unspliced
  seurat_obj[["RNA"]]@layers$spliced <- spliced
  
  seurat_obj
}

#' Export Velocity Results
#'
#' @description Export velocity results from Seurat object to list format.
#'
#' @param object Seurat object with velocity computed
#'
#' @return List with velocity results
#' @export
export_velocity_results <- function(object) {
  
  if (is.null(object@misc$scVeloR)) {
    stop("No scVeloR results found")
  }
  
  results <- list()
  
  # Velocity
  if (!is.null(object@misc$scVeloR$velocity)) {
    results$velocity <- object@misc$scVeloR$velocity$velocity_s
    results$gamma <- object@misc$scVeloR$velocity$gamma
    results$r2 <- object@misc$scVeloR$velocity$r2
    results$velocity_genes <- object@misc$scVeloR$velocity$velocity_genes
  }
  
  # Dynamics
  if (!is.null(object@misc$scVeloR$dynamics)) {
    results$dynamics <- object@misc$scVeloR$dynamics$gene_params
    results$fit_t <- object@misc$scVeloR$dynamics$fit_t
  }
  
  # Velocity embedding
  if (!is.null(object@misc$scVeloR$velocity_embedding)) {
    results$velocity_embedding <- object@misc$scVeloR$velocity_embedding
  }
  
  # Latent time
  if ("latent_time" %in% colnames(object@meta.data)) {
    results$latent_time <- object@meta.data$latent_time
  }
  
  results
}

#' Import from AnnData
#'
#' @description Import velocity data from AnnData h5ad file.
#' Requires the anndata package.
#'
#' @param h5ad_path Path to h5ad file
#'
#' @return Seurat object with velocity layers
#' @export
import_from_anndata <- function(h5ad_path) {
  
  if (!requireNamespace("anndata", quietly = TRUE)) {
    stop("anndata package required. Install with: pip install anndata, then reticulate::py_install('anndata')")
  }
  
  # Read AnnData
  adata <- anndata::read_h5ad(h5ad_path)
  
  # Extract matrices
  if ("spliced" %in% names(adata$layers)) {
    spliced <- adata$layers[["spliced"]]
  } else {
    spliced <- adata$X
  }
  
  if ("unspliced" %in% names(adata$layers)) {
    unspliced <- adata$layers[["unspliced"]]
  } else {
    stop("No unspliced layer found in AnnData")
  }
  
  # Convert to R matrices
  spliced <- as.matrix(spliced)
  unspliced <- as.matrix(unspliced)
  
  # Get metadata
  meta_data <- as.data.frame(adata$obs)
  
  # Create Seurat object
  seurat_obj <- create_velocity_seurat(
    spliced = t(spliced),  # Transpose to genes x cells
    unspliced = t(unspliced),
    meta.data = meta_data
  )
  
  # Import embeddings if available
  if ("X_umap" %in% names(adata$obsm)) {
    umap_coords <- adata$obsm[["X_umap"]]
    colnames(umap_coords) <- paste0("UMAP_", 1:ncol(umap_coords))
    rownames(umap_coords) <- colnames(seurat_obj)
    
    seurat_obj[["umap"]] <- Seurat::CreateDimReducObject(
      embeddings = umap_coords,
      key = "UMAP_",
      assay = "RNA"
    )
  }
  
  seurat_obj
}

#' Export to AnnData Format
#'
#' @description Export Seurat object with velocity to AnnData format.
#'
#' @param object Seurat object
#' @param filename Output h5ad filename
#'
#' @return Path to saved file
#' @export
export_to_anndata <- function(object, filename) {
  
  if (!requireNamespace("anndata", quietly = TRUE)) {
    stop("anndata package required")
  }
  
  # Get count matrix
  counts <- Seurat::GetAssayData(object, slot = "counts")
  
  # Get metadata
  obs <- object@meta.data
  
  # Create AnnData
  adata <- anndata::AnnData(
    X = t(as.matrix(counts)),
    obs = obs
  )
  
  # Add velocity layers if present
  if (!is.null(object@misc$scVeloR$Ms)) {
    adata$layers[["Ms"]] <- object@misc$scVeloR$Ms
    adata$layers[["Mu"]] <- object@misc$scVeloR$Mu
  }
  
  if (!is.null(object@misc$scVeloR$velocity$velocity_s)) {
    adata$layers[["velocity"]] <- object@misc$scVeloR$velocity$velocity_s
  }
  
  # Add embeddings
  if ("umap" %in% names(object@reductions)) {
    adata$obsm[["X_umap"]] <- Seurat::Embeddings(object, "umap")
  }
  
  # Save
  adata$write_h5ad(filename)
  
  filename
}

#' @title Data Conversion Functions
#' @description Functions for converting data between Seurat objects and 
#' scVeloR internal formats.
#' @name data_conversion
NULL

#' Import Velocity Data from Seurat
#'
#' @description Extract spliced and unspliced counts from Seurat object.
#' Supports both Seurat V4 and V5 object structures. Priority: V4 compatibility.
#'
#' @param object Seurat object
#' @param spliced_layer Name of spliced count layer (default: "spliced")
#' @param unspliced_layer Name of unspliced count layer (default: "unspliced")
#' @param assay Assay name (default: DefaultAssay)
#'
#' @return List with spliced and unspliced matrices (cells x genes)
#' @export
import_velocity_data <- function(object,
                                  spliced_layer = "spliced",
                                  unspliced_layer = "unspliced",
                                  assay = NULL) {
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
  }
  
  # Use get_layer_matrix which handles V4/V5 compatibility
  spliced <- tryCatch({
    get_layer_matrix(object, spliced_layer, assay = assay)
  }, error = function(e) {
    # Fallback: try "counts" as spliced
    message("Note: 'spliced' layer not found, trying 'counts' as spliced data")
    tryCatch({
      get_layer_matrix(object, "counts", assay = assay)
    }, error = function(e2) NULL)
  })
  
  unspliced <- tryCatch({
    get_layer_matrix(object, unspliced_layer, assay = assay)
  }, error = function(e) {
    NULL
  })
  
  # Validate data
  if (is.null(spliced)) {
    stop(sprintf(
      "Could not find spliced data.\n%s",
      "Ensure spliced counts are in the Seurat object as a layer or assay."
    ))
  }
  
  if (is.null(unspliced)) {
    stop(sprintf(
      "Could not find unspliced data in layer '%s'.\n%s\n%s",
      unspliced_layer,
      "For velocity analysis, unspliced counts are required.",
      "Load velocity data using velocyto, kallisto-bustools, or STARsolo output."
    ))
  }
  
  # Ensure consistent dimensions
  if (!all(dim(spliced) == dim(unspliced))) {
    stop("Spliced and unspliced matrices have different dimensions")
  }
  
  list(
    spliced = spliced,
    unspliced = unspliced,
    gene_names = colnames(spliced),
    cell_names = rownames(spliced)
  )
}

#' Create Seurat Object with Velocity Layers
#'
#' @description Create a Seurat object with spliced and unspliced layers
#' from count matrices. Compatible with both Seurat V4 and V5.
#'
#' @param spliced Spliced count matrix (cells x genes or genes x cells)
#' @param unspliced Unspliced count matrix
#' @param meta.data Cell metadata (optional)
#' @param project Project name
#'
#' @return Seurat object with velocity data stored appropriately for the Seurat version
#' @export
create_velocity_seurat <- function(spliced,
                                    unspliced,
                                    meta.data = NULL,
                                    project = "scVeloR") {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required")
  }
  
  # Ensure matrices are in genes x cells format for Seurat
  if (nrow(spliced) == nrow(unspliced) && ncol(spliced) == ncol(unspliced)) {
    # Assume genes x cells if more rows than columns
    if (nrow(spliced) <= ncol(spliced)) {
      # Transpose to genes x cells
      spliced <- t(spliced)
      unspliced <- t(unspliced)
    }
  }
  
  # Check Seurat package version
  pkg_version <- utils::packageVersion("Seurat")
  major_version <- as.integer(strsplit(as.character(pkg_version), "\\.")[[1]][1])
  
  if (major_version >= 5) {
    # Seurat V5: Create with layers
    seurat_obj <- Seurat::CreateSeuratObject(
      counts = spliced,
      project = project,
      meta.data = meta.data
    )
    
    # Add layers using V5 method
    tryCatch({
      seurat_obj[["RNA"]] <- Seurat::SetAssayData(
        seurat_obj[["RNA"]], 
        layer = "spliced", 
        new.data = spliced
      )
      seurat_obj[["RNA"]] <- Seurat::SetAssayData(
        seurat_obj[["RNA"]], 
        layer = "unspliced", 
        new.data = unspliced
      )
    }, error = function(e) {
      # Fallback: store in misc
      message("Note: Storing velocity data in misc slot (V5 layer creation failed)")
      seurat_obj@misc$spliced <<- spliced
      seurat_obj@misc$unspliced <<- unspliced
    })
    
  } else {
    # Seurat V4: Create object and store velocity data in misc
    seurat_obj <- Seurat::CreateSeuratObject(
      counts = spliced,
      project = project,
      meta.data = meta.data
    )
    
    # For V4, store spliced/unspliced in misc (standard slots only support counts/data/scale.data)
    seurat_obj@misc$spliced <- spliced
    seurat_obj@misc$unspliced <- unspliced
    
    # Also create separate assays for spliced and unspliced (alternative V4 approach)
    tryCatch({
      seurat_obj[["spliced"]] <- Seurat::CreateAssayObject(counts = spliced)
      seurat_obj[["unspliced"]] <- Seurat::CreateAssayObject(counts = unspliced)
    }, error = function(e) {
      message("Note: Could not create separate assays for spliced/unspliced")
    })
  }
  
  # Initialize scVeloR storage
  seurat_obj@misc$scVeloR <- list()
  
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
#' Compatible with both Seurat V4 and V5.
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
  
  # Get count matrix (compatible with V4 and V5)
  version <- seurat_version(object)
  
  if (version == 5) {
    counts <- tryCatch({
      Seurat::LayerData(object, layer = "counts")
    }, error = function(e) {
      # Try getting first available layer
      layers <- Seurat::Layers(object)
      if (length(layers) > 0) {
        Seurat::LayerData(object, layer = layers[1])
      } else {
        stop("No count data found")
      }
    })
  } else {
    counts <- Seurat::GetAssayData(object, slot = "counts")
  }
  
  # Get metadata
  obs <- object@meta.data
  
  # Create AnnData (needs genes x cells -> cells x genes)
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
  
  # Add spliced/unspliced if available
  tryCatch({
    spliced <- get_layer_matrix(object, "spliced")
    unspliced <- get_layer_matrix(object, "unspliced")
    adata$layers[["spliced"]] <- spliced
    adata$layers[["unspliced"]] <- unspliced
  }, error = function(e) {
    # Spliced/unspliced not available
  })
  
  # Add embeddings
  if ("umap" %in% names(object@reductions)) {
    adata$obsm[["X_umap"]] <- Seurat::Embeddings(object, "umap")
  }
  
  if ("pca" %in% names(object@reductions)) {
    adata$obsm[["X_pca"]] <- Seurat::Embeddings(object, "pca")
  }
  
  # Save
  adata$write_h5ad(filename)
  
  filename
}

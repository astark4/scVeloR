#' @title Data Conversion and Seurat Integration
#' @description Functions for extracting and setting velocity data in Seurat objects.
#' Supports both Seurat V4 and V5 data structures.
#' @name data_conversion
NULL

#' Detect Seurat Version
#'
#' @description Internal function to detect whether Seurat V4 or V5 is being used.
#'
#' @param seurat_obj A Seurat object
#' @return Character string "v4" or "v5"
#' @keywords internal
.get_seurat_version <- function(seurat_obj) {
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object", call. = FALSE)
  }
  
  pkg_version <- utils::packageVersion("Seurat")
  if (pkg_version >= "5.0.0") {
    return("v5")
  }
  return("v4")
}

#' Check if Layer Exists
#'
#' @description Check if a layer/assay exists in the Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param layer_name Name of the layer to check
#' @param assay_name Name of the assay (default: NULL uses default assay)
#' @return Logical
#' @keywords internal
.layer_exists <- function(seurat_obj, layer_name, assay_name = NULL) {
  if (is.null(assay_name)) {
    assay_name <- SeuratObject::DefaultAssay(seurat_obj)
  }
  
  version <- .get_seurat_version(seurat_obj)
  
  if (version == "v5") {
    # V5: Check layers
    tryCatch({
      layers <- SeuratObject::Layers(seurat_obj[[assay_name]])
      return(layer_name %in% layers)
    }, error = function(e) {
      return(FALSE)
    })
  } else {
    # V4: Check slots or separate assays
    if (layer_name %in% c("counts", "data", "scale.data")) {
      slot_data <- tryCatch(
        SeuratObject::GetAssayData(seurat_obj, slot = layer_name, assay = assay_name),
        error = function(e) NULL
      )
      return(!is.null(slot_data) && length(slot_data) > 0)
    } else {
      # Check if it's a separate assay
      return(layer_name %in% names(seurat_obj@assays))
    }
  }
}

#' Get Layer Data
#'
#' @description Get data from a layer in the Seurat object. Works with both V4 and V5.
#'
#' @param seurat_obj A Seurat object
#' @param layer_name Name of the layer
#' @param assay_name Name of the assay (default: NULL uses default assay)
#' @return A sparse matrix
#' @keywords internal
.get_layer_data <- function(seurat_obj, layer_name, assay_name = NULL) {
  if (is.null(assay_name)) {
    assay_name <- SeuratObject::DefaultAssay(seurat_obj)
  }
  
  version <- .get_seurat_version(seurat_obj)
  
  if (version == "v5") {
    # V5: Use LayerData
    tryCatch({
      data <- SeuratObject::LayerData(seurat_obj, layer = layer_name, assay = assay_name)
      return(as(data, "CsparseMatrix"))
    }, error = function(e) {
      # Try as separate assay
      if (layer_name %in% names(seurat_obj@assays)) {
        return(as(SeuratObject::GetAssayData(seurat_obj, assay = layer_name), "CsparseMatrix"))
      }
      stop("Layer '", layer_name, "' not found", call. = FALSE)
    })
  } else {
    # V4: Use GetAssayData or separate assay
    if (layer_name %in% c("counts", "data", "scale.data")) {
      return(as(SeuratObject::GetAssayData(seurat_obj, slot = layer_name, assay = assay_name), "CsparseMatrix"))
    } else if (layer_name %in% names(seurat_obj@assays)) {
      return(as(SeuratObject::GetAssayData(seurat_obj, assay = layer_name), "CsparseMatrix"))
    } else {
      stop("Layer '", layer_name, "' not found", call. = FALSE)
    }
  }
}

#' Set Layer Data
#'
#' @description Set data to a layer in the Seurat object. Works with both V4 and V5.
#'
#' @param seurat_obj A Seurat object
#' @param layer_name Name of the layer
#' @param data Matrix data to set
#' @param assay_name Name of the assay (default: NULL uses default assay)
#' @return Modified Seurat object
#' @keywords internal
.set_layer_data <- function(seurat_obj, layer_name, data, assay_name = NULL) {
  if (is.null(assay_name)) {
    assay_name <- SeuratObject::DefaultAssay(seurat_obj)
  }
  
  # Ensure sparse matrix
  if (!inherits(data, "sparseMatrix")) {
    data <- as(data, "CsparseMatrix")
  }
  
  version <- .get_seurat_version(seurat_obj)
  
  if (version == "v5") {
    # V5: Use LayerData<-
    tryCatch({
      seurat_obj[[assay_name]] <- SeuratObject::SetAssayData(
        seurat_obj[[assay_name]], 
        layer = layer_name, 
        new.data = data
      )
    }, error = function(e) {
      # Create new layer if needed
      seurat_obj[[assay_name]][[layer_name]] <- data
    })
  } else {
    # V4: Use SetAssayData for standard slots, or create new assay
    if (layer_name %in% c("counts", "data", "scale.data")) {
      seurat_obj <- SeuratObject::SetAssayData(
        seurat_obj, 
        slot = layer_name, 
        new.data = data, 
        assay = assay_name
      )
    } else {
      # Store in misc or as new assay
      seurat_obj@misc[[layer_name]] <- data
    }
  }
  
  return(seurat_obj)
}

#' Extract Velocity Data from Seurat Object
#'
#' @description Extract spliced and unspliced count matrices from a Seurat object.
#'
#' @param seurat_obj A Seurat object containing spliced and unspliced counts
#' @param spliced_layer Name of the spliced layer/assay (default: "spliced")
#' @param unspliced_layer Name of the unspliced layer/assay (default: "unspliced")
#' @param use_genes Character vector of genes to use, or NULL for all genes (default: NULL
#'
#' @return A list containing:
#' \describe{
#'   \item{spliced}{Sparse matrix of spliced counts (genes x cells)}
#'   \item{unspliced}{Sparse matrix of unspliced counts (genes x cells)}
#'   \item{genes}{Character vector of gene names}
#'   \item{cells}{Character vector of cell names}
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' data <- extract_velocity_data(seurat_obj)
#' }
extract_velocity_data <- function(seurat_obj,
                                   spliced_layer = "spliced",
                                   unspliced_layer = "unspliced",
                                   use_genes = NULL) {
  
  # Validate input
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object", call. = FALSE)
  }
  
  # Try to find spliced data
  spliced <- NULL
  for (layer in c(spliced_layer, "counts", "data")) {
    if (.layer_exists(seurat_obj, layer)) {
      spliced <- .get_layer_data(seurat_obj, layer)
      if (getOption("scVeloR.verbose", TRUE)) {
        message("Using '", layer, "' as spliced counts")
      }
      break
    }
  }
  
  if (is.null(spliced)) {
    stop("Could not find spliced counts. Please provide '", spliced_layer, "' layer.", call. = FALSE)
  }
  
  # Try to find unspliced data
  unspliced <- NULL
  if (.layer_exists(seurat_obj, unspliced_layer)) {
    unspliced <- .get_layer_data(seurat_obj, unspliced_layer)
  } else {
    stop("Could not find unspliced counts. Please provide '", unspliced_layer, "' layer.", call. = FALSE)
  }
  
  # Ensure dimensions match
  if (!all(dim(spliced) == dim(unspliced))) {
    stop("Spliced and unspliced matrices must have the same dimensions", call. = FALSE)
  }
  
  # Get gene and cell names
  genes <- rownames(spliced)
  cells <- colnames(spliced)
  
  # Subset genes if specified
  if (!is.null(use_genes)) {
    use_genes <- intersect(use_genes, genes)
    if (length(use_genes) == 0) {
      stop("No specified genes found in the data", call. = FALSE)
    }
    spliced <- spliced[use_genes, , drop = FALSE]
    unspliced <- unspliced[use_genes, , drop = FALSE]
    genes <- use_genes
  }
  
  # Return as list
  list(
    spliced = spliced,
    unspliced = unspliced,
    genes = genes,
    cells = cells
  )
}

#' Set Velocity Data to Seurat Object
#'
#' @description Store velocity-related data back into a Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param velocity Velocity matrix (genes x cells) or NULL
#' @param velocity_u Unspliced velocity matrix or NULL
#' @param Ms First-order moment of spliced counts or NULL
#' @param Mu First-order moment of unspliced counts or NULL
#' @param velocity_graph Velocity graph (sparse matrix) or NULL
#' @param velocity_graph_neg Negative velocity graph or NULL
#' @param velocity_embedding Velocity embedding coordinates or NULL
#' @param fit_params Data frame of fitted parameters or NULL
#'
#' @return Modified Seurat object
#'
#' @export
set_velocity_data <- function(seurat_obj,
                               velocity = NULL,
                               velocity_u = NULL,
                               Ms = NULL,
                               Mu = NULL,
                               velocity_graph = NULL,
                               velocity_graph_neg = NULL,
                               velocity_embedding = NULL,
                               fit_params = NULL) {
  
  # Store velocity
  if (!is.null(velocity)) {
    seurat_obj <- .set_layer_data(seurat_obj, "velocity", velocity)
  }
  
  if (!is.null(velocity_u)) {
    seurat_obj <- .set_layer_data(seurat_obj, "velocity_u", velocity_u)
  }
  
  # Store moments
  if (!is.null(Ms)) {
    seurat_obj <- .set_layer_data(seurat_obj, "Ms", Ms)
  }
  
  if (!is.null(Mu)) {
    seurat_obj <- .set_layer_data(seurat_obj, "Mu", Mu)
  }
  
  # Store velocity graph in misc
  if (!is.null(velocity_graph)) {
    seurat_obj@misc[["velocity_graph"]] <- velocity_graph
  }
  
  if (!is.null(velocity_graph_neg)) {
    seurat_obj@misc[["velocity_graph_neg"]] <- velocity_graph_neg
  }
  
  # Store velocity embedding
  if (!is.null(velocity_embedding)) {
    # Store as a reduction or in misc
    seurat_obj@misc[["velocity_embedding"]] <- velocity_embedding
  }
  
  # Store fit parameters
  if (!is.null(fit_params)) {
    # Add to gene metadata
    current_meta <- seurat_obj[[SeuratObject::DefaultAssay(seurat_obj)]]@meta.features
    for (col in colnames(fit_params)) {
      current_meta[[col]] <- fit_params[[col]][match(rownames(current_meta), rownames(fit_params))]
    }
    seurat_obj[[SeuratObject::DefaultAssay(seurat_obj)]]@meta.features <- current_meta
  }
  
  return(seurat_obj)
}

#' Get Embeddings from Seurat Object
#'
#' @description Extract embedding coordinates from a Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param basis Name of the dimensional reduction (e.g., "umap", "tsne", "pca")
#' @param dims Dimensions to use (default: NULL uses all)
#'
#' @return Matrix of embedding coordinates (cells x dimensions)
#'
#' @export
get_embedding <- function(seurat_obj, basis = "umap", dims = NULL) {
  
  # Handle basis name
  if (!grepl("^[A-Za-z]", basis)) {
    basis <- paste0("X_", basis)
  }
  basis <- tolower(gsub("^X_", "", basis))
  
  # Check if reduction exists
  reductions <- names(seurat_obj@reductions)
  if (!basis %in% reductions) {
    stop("Reduction '", basis, "' not found. Available: ", paste(reductions, collapse = ", "), call. = FALSE)
  }
  
  # Get embeddings
  emb <- Seurat::Embeddings(seurat_obj, reduction = basis)
  
  # Subset dimensions if specified
  if (!is.null(dims)) {
    if (max(dims) > ncol(emb)) {
      stop("Requested dimensions exceed available dimensions", call. = FALSE)
    }
    emb <- emb[, dims, drop = FALSE]
  }
  
  return(emb)
}

#' Get Variable Genes
#'
#' @description Get highly variable genes from Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param n_top Number of top variable genes to return (default: NULL returns all)
#'
#' @return Character vector of gene names
#'
#' @export
get_var_genes <- function(seurat_obj, n_top = NULL) {
  
  var_genes <- tryCatch({
    Seurat::VariableFeatures(seurat_obj)
  }, error = function(e) {
    character(0)
  })
  
  if (length(var_genes) == 0) {
    warning("No variable features found. Consider running FindVariableFeatures first.")
    return(rownames(seurat_obj))
  }
  
  if (!is.null(n_top) && length(var_genes) > n_top) {
    var_genes <- var_genes[seq_len(n_top)]
  }
  
  return(var_genes)
}

#' Show Proportions of Spliced/Unspliced Counts
#'
#' @description Display the proportions of spliced and unspliced counts in the data.
#'
#' @param seurat_obj A Seurat object
#' @param spliced_layer Name of the spliced layer (default: "spliced")
#' @param unspliced_layer Name of the unspliced layer (default: "unspliced")
#'
#' @return Invisibly returns a named vector of proportions
#'
#' @export
show_proportions <- function(seurat_obj,
                              spliced_layer = "spliced",
                              unspliced_layer = "unspliced") {
  
  # Extract data
  data <- extract_velocity_data(seurat_obj, spliced_layer, unspliced_layer)
  
  # Calculate sums per cell
  spliced_sum <- Matrix::colSums(data$spliced)
  unspliced_sum <- Matrix::colSums(data$unspliced)
  total_sum <- spliced_sum + unspliced_sum
  
  # Avoid division by zero
  total_sum[total_sum == 0] <- 1
  
  # Calculate proportions
  spliced_prop <- mean(spliced_sum / total_sum)
  unspliced_prop <- mean(unspliced_sum / total_sum)
  
  # Print results
  message("Abundance of counts:")
  message("  Spliced:   ", round(spliced_prop * 100, 2), "%")
  message("  Unspliced: ", round(unspliced_prop * 100, 2), "%")
  
  invisible(c(spliced = spliced_prop, unspliced = unspliced_prop))
}

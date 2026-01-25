#' @title Preprocessing Functions
#' @description Data preprocessing functions for RNA velocity analysis.
#' @name preprocessing
NULL

#' Filter Genes
#'
#' @description Filter genes based on expression criteria.
#'
#' @param object Seurat object
#' @param min_counts Minimum total counts per gene
#' @param min_cells Minimum cells expressing the gene
#' @param min_shared_counts Minimum shared counts (u and s > 0)
#' @param min_shared_cells Minimum cells with both u and s
#' @param spliced_layer Name of spliced layer
#' @param unspliced_layer Name of unspliced layer
#' @param verbose Print progress
#'
#' @return Seurat object with filtered gene metadata
#' @export
filter_genes <- function(object,
                         min_counts = 20,
                         min_cells = 10,
                         min_shared_counts = 10,
                         min_shared_cells = 5,
                         spliced_layer = "spliced",
                         unspliced_layer = "unspliced",
                         verbose = TRUE) {
  
  # Get spliced and unspliced matrices
  spliced <- get_layer_matrix(object, spliced_layer)
  unspliced <- get_layer_matrix(object, unspliced_layer)
  
  n_genes <- ncol(spliced)
  
  # Compute gene statistics
  s_counts <- colSums(spliced)
  u_counts <- colSums(unspliced)
  s_cells <- colSums(spliced > 0)
  u_cells <- colSums(unspliced > 0)
  
  # Shared expression (both u and s > 0)
  shared <- (spliced > 0) & (unspliced > 0)
  shared_counts <- colSums(spliced * shared + unspliced * shared)
  shared_cells <- colSums(shared)
  
  # Apply filters
  pass_filter <- (s_counts >= min_counts) &
                 (u_counts >= min_counts) &
                 (s_cells >= min_cells) &
                 (u_cells >= min_cells) &
                 (shared_counts >= min_shared_counts) &
                 (shared_cells >= min_shared_cells)
  
  n_passed <- sum(pass_filter)
  
  if (verbose) {
    message(sprintf("Filtered genes: %d/%d passed (%.1f%%)", 
                    n_passed, n_genes, 100 * n_passed / n_genes))
  }
  
  # Store filter results
  if (is.null(object@misc$scVeloR)) {
    object@misc$scVeloR <- list()
  }
  
  object@misc$scVeloR$gene_filter <- pass_filter
  object@misc$scVeloR$filtered_genes <- colnames(spliced)[pass_filter]
  
  object
}

#' Filter Cells
#'
#' @description Filter cells based on expression criteria.
#'
#' @param object Seurat object
#' @param min_counts Minimum total counts per cell
#' @param min_genes Minimum genes expressed per cell
#' @param spliced_layer Name of spliced layer
#' @param unspliced_layer Name of unspliced layer
#' @param verbose Print progress
#'
#' @return Seurat object (cells in metadata)
#' @export
filter_cells <- function(object,
                         min_counts = 100,
                         min_genes = 50,
                         spliced_layer = "spliced",
                         unspliced_layer = "unspliced",
                         verbose = TRUE) {
  
  spliced <- get_layer_matrix(object, spliced_layer)
  unspliced <- get_layer_matrix(object, unspliced_layer)
  
  n_cells <- nrow(spliced)
  
  # Compute cell statistics
  total_counts <- rowSums(spliced) + rowSums(unspliced)
  n_genes <- rowSums(spliced > 0) + rowSums(unspliced > 0)
  
  pass_filter <- (total_counts >= min_counts) & (n_genes >= min_genes)
  
  if (verbose) {
    message(sprintf("Cell filter: %d/%d passed (%.1f%%)", 
                    sum(pass_filter), n_cells, 100 * sum(pass_filter) / n_cells))
  }
  
  object@meta.data$pass_velocity_filter <- pass_filter
  
  object
}

#' Normalize Data
#'
#' @description Normalize spliced and unspliced counts.
#'
#' @param object Seurat object
#' @param method Normalization method: "log1p" or "total"
#' @param target_sum Target sum for total normalization
#' @param spliced_layer Name of spliced layer
#' @param unspliced_layer Name of unspliced layer
#'
#' @return Seurat object with normalized layers
#' @export
normalize_layers <- function(object,
                             method = "log1p",
                             target_sum = 10000,
                             spliced_layer = "spliced",
                             unspliced_layer = "unspliced") {
  
  spliced <- get_layer_matrix(object, spliced_layer)
  unspliced <- get_layer_matrix(object, unspliced_layer)
  
  if (method == "log1p") {
    # Size factor normalization + log transform
    s_size <- rowSums(spliced) / median(rowSums(spliced))
    u_size <- rowSums(unspliced) / median(rowSums(unspliced))
    
    s_size[s_size == 0] <- 1
    u_size[u_size == 0] <- 1
    
    spliced_norm <- log1p(spliced / s_size)
    unspliced_norm <- log1p(unspliced / u_size)
    
  } else if (method == "total") {
    s_sum <- rowSums(spliced)
    u_sum <- rowSums(unspliced)
    
    s_sum[s_sum == 0] <- 1
    u_sum[u_sum == 0] <- 1
    
    spliced_norm <- spliced / s_sum * target_sum
    unspliced_norm <- unspliced / u_sum * target_sum
  }
  
  # Store normalized data
  if (is.null(object@misc$scVeloR)) {
    object@misc$scVeloR <- list()
  }
  
  object@misc$scVeloR$spliced_norm <- spliced_norm
  object@misc$scVeloR$unspliced_norm <- unspliced_norm
  
  object
}

#' Compute First-Order Moments
#'
#' @description Compute neighbor-averaged expression (first moments).
#' Ms = smooth(spliced), Mu = smooth(unspliced)
#'
#' @param object Seurat object with neighbors computed
#' @param n_neighbors Number of neighbors (default uses existing neighbors)
#' @param use_normalized Use normalized data
#' @param spliced_layer Name of spliced layer
#' @param unspliced_layer Name of unspliced layer
#' @param verbose Print progress
#'
#' @return Seurat object with Ms and Mu in misc$scVeloR
#' @export
moments <- function(object,
                    n_neighbors = 30,
                    use_normalized = TRUE,
                    spliced_layer = "spliced",
                    unspliced_layer = "unspliced",
                    verbose = TRUE) {
  
  if (verbose) {
    message("Computing moments...")
  }
  
  # Ensure neighbors are computed
  if (is.null(object@misc$scVeloR$neighbors)) {
    object <- compute_neighbors(object, n_neighbors = n_neighbors, verbose = verbose)
  }
  
  conn <- object@misc$scVeloR$neighbors$connectivities
  
  # Get expression data
  if (use_normalized && !is.null(object@misc$scVeloR$spliced_norm)) {
    spliced <- object@misc$scVeloR$spliced_norm
    unspliced <- object@misc$scVeloR$unspliced_norm
  } else {
    spliced <- get_layer_matrix(object, spliced_layer)
    unspliced <- get_layer_matrix(object, unspliced_layer)
  }
  
  # Row normalize connectivities
  row_sums <- rowSums(conn)
  row_sums[row_sums == 0] <- 1
  conn_norm <- conn / row_sums
  
  # Compute moments (neighbor averaging)
  Ms <- as.matrix(conn_norm %*% spliced)
  Mu <- as.matrix(conn_norm %*% unspliced)
  
  # Ensure same dimensions and names
  colnames(Ms) <- colnames(spliced)
  colnames(Mu) <- colnames(unspliced)
  rownames(Ms) <- rownames(spliced)
  rownames(Mu) <- rownames(unspliced)
  
  # Store results
  object@misc$scVeloR$Ms <- Ms
  object@misc$scVeloR$Mu <- Mu
  
  if (verbose) {
    message(sprintf("  Computed moments: Ms %dx%d, Mu %dx%d", 
                    nrow(Ms), ncol(Ms), nrow(Mu), ncol(Mu)))
  }
  
  object
}

#' Compute Second-Order Moments
#'
#' @description Compute second-order moments for stochastic velocity.
#' Mss, Mus, Muu represent variances and covariances.
#'
#' @param object Seurat object with neighbors computed
#' @param use_normalized Use normalized data
#' @param spliced_layer Name of spliced layer
#' @param unspliced_layer Name of unspliced layer
#' @param verbose Print progress
#'
#' @return Seurat object with Mss, Mus, Muu in misc$scVeloR
#' @export
second_order_moments <- function(object,
                                  use_normalized = TRUE,
                                  spliced_layer = "spliced",
                                  unspliced_layer = "unspliced",
                                  verbose = TRUE) {
  
  if (verbose) {
    message("Computing second-order moments...")
  }
  
  if (is.null(object@misc$scVeloR$neighbors)) {
    stop("Run compute_neighbors() first")
  }
  
  conn <- object@misc$scVeloR$neighbors$connectivities
  
  # Get expression data
  if (use_normalized && !is.null(object@misc$scVeloR$spliced_norm)) {
    s <- object@misc$scVeloR$spliced_norm
    u <- object@misc$scVeloR$unspliced_norm
  } else {
    s <- get_layer_matrix(object, spliced_layer)
    u <- get_layer_matrix(object, unspliced_layer)
  }
  
  # Row normalize
  row_sums <- rowSums(conn)
  row_sums[row_sums == 0] <- 1
  conn_norm <- conn / row_sums
  
  # Compute <s^2>, <su>, <u^2>
  Mss <- as.matrix(conn_norm %*% (s^2))
  Mus <- as.matrix(conn_norm %*% (s * u))
  Muu <- as.matrix(conn_norm %*% (u^2))
  
  # Store
  object@misc$scVeloR$Mss <- Mss
  object@misc$scVeloR$Mus <- Mus
  object@misc$scVeloR$Muu <- Muu
  
  if (verbose) {
    message("  Computed second-order moments")
  }
  
  object
}

#' Get Layer Matrix
#'
#' @description Extract layer matrix from Seurat object.
#' Handles both Seurat V4 and V5 object structures.
#'
#' @param object Seurat object
#' @param layer_name Name of layer
#'
#' @return Dense matrix (cells x genes)
#' @keywords internal
get_layer_matrix <- function(object, layer_name) {
  
  # Try Seurat V5 structure first
  if ("Assay5" %in% class(object@assays[[Seurat::DefaultAssay(object)]])) {
    # V5: Use LayerData
    if (layer_name %in% Seurat::Layers(object)) {
      mat <- Seurat::LayerData(object, layer = layer_name)
    } else {
      stop(sprintf("Layer '%s' not found", layer_name))
    }
  } else {
    # V4: Use GetAssayData
    if (layer_name %in% names(object@assays[[Seurat::DefaultAssay(object)]]@layers)) {
      mat <- Seurat::GetAssayData(object, slot = layer_name)
    } else if (layer_name %in% slotNames(object@assays[[Seurat::DefaultAssay(object)]])) {
      mat <- slot(object@assays[[Seurat::DefaultAssay(object)]], layer_name)
    } else {
      stop(sprintf("Layer '%s' not found", layer_name))
    }
  }
  
  # Convert to dense matrix (cells x genes)
  if (inherits(mat, "sparseMatrix")) {
    mat <- as.matrix(mat)
  }
  
  if (nrow(mat) != ncol(object)) {
    mat <- t(mat)
  }
  
  mat
}

#' Get Velocity Data
#'
#' @description Get expression data for velocity computation.
#' Returns Ms/Mu if available, otherwise raw spliced/unspliced.
#'
#' @param object Seurat object
#' @param spliced_layer Name of spliced layer
#' @param unspliced_layer Name of unspliced layer
#' @param use_moments Use moment-averaged data if available
#'
#' @return List with Ms and Mu matrices
#' @keywords internal
get_velocity_data <- function(object,
                              spliced_layer = "spliced",
                              unspliced_layer = "unspliced",
                              use_moments = TRUE) {
  
  if (use_moments && !is.null(object@misc$scVeloR$Ms)) {
    Ms <- object@misc$scVeloR$Ms
    Mu <- object@misc$scVeloR$Mu
  } else {
    Ms <- get_layer_matrix(object, spliced_layer)
    Mu <- get_layer_matrix(object, unspliced_layer)
  }
  
  list(Ms = Ms, Mu = Mu)
}

#' Prepare Seurat Object for Velocity
#'
#' @description Convenience function to run all preprocessing steps.
#'
#' @param object Seurat object
#' @param min_counts Minimum counts for gene filtering
#' @param min_cells Minimum cells for gene filtering
#' @param n_neighbors Number of neighbors
#' @param n_pcs Number of PCs for neighbor computation
#' @param spliced_layer Name of spliced layer
#' @param unspliced_layer Name of unspliced layer
#' @param verbose Print progress
#'
#' @return Processed Seurat object ready for velocity computation
#' @export
prepare_velocity <- function(object,
                             min_counts = 20,
                             min_cells = 10,
                             n_neighbors = 30,
                             n_pcs = 30,
                             spliced_layer = "spliced",
                             unspliced_layer = "unspliced",
                             verbose = TRUE) {
  
  if (verbose) {
    message("Preparing Seurat object for velocity analysis...")
  }
  
  # Filter genes
  object <- filter_genes(object, 
                         min_counts = min_counts, 
                         min_cells = min_cells,
                         spliced_layer = spliced_layer,
                         unspliced_layer = unspliced_layer,
                         verbose = verbose)
  
  # Normalize
  object <- normalize_layers(object,
                             spliced_layer = spliced_layer,
                             unspliced_layer = unspliced_layer)
  
  if (verbose) message("  Normalized expression data")
  
  # Compute neighbors
  object <- compute_neighbors(object, 
                              n_neighbors = n_neighbors, 
                              n_pcs = n_pcs,
                              verbose = verbose)
  
  # Compute moments
  object <- moments(object, verbose = verbose)
  
  if (verbose) {
    message("Done. Object is ready for velocity computation.")
  }
  
  object
}

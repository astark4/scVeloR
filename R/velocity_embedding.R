#' @title Velocity Embedding Projection
#' @description Functions for projecting velocities onto low-dimensional embeddings.
#' @name velocity_embedding
NULL

#' Project Velocity to Embedding
#'
#' @description Project high-dimensional velocity vectors onto a 2D embedding space.
#'
#' @param seurat_obj A Seurat object with computed velocity graph
#' @param basis Name of the embedding to use (default: "umap")
#' @param vkey Key for velocity layer (default: "velocity")
#' @param scale Scale factor for velocity arrows (default: 10)
#' @param self_transitions Include self-transitions (default: TRUE)
#' @param use_negative Use negative correlations (default: FALSE)
#' @param autoscale Auto-scale velocities per cell (default: TRUE)
#' @param all_comps Use all embedding components (default: FALSE)
#' @param direct_pca Use direct PCA projection (default: FALSE)
#' @param n_pcs Number of PCs for direct projection (default: 30)
#'
#' @return Modified Seurat object with velocity_embedding in misc
#'
#' @export
project_velocity_embedding <- function(seurat_obj,
                                         basis = "umap",
                                         vkey = "velocity",
                                         scale = 10,
                                         self_transitions = TRUE,
                                         use_negative = FALSE,
                                         autoscale = TRUE,
                                         all_comps = FALSE,
                                         direct_pca = FALSE,
                                         n_pcs = 30) {
  
  .vmessage("Projecting velocity to ", basis, " embedding")
  
  # Get embedding
  emb <- get_embedding(seurat_obj, basis)
  n_cells <- nrow(emb)
  n_dims <- ifelse(all_comps, ncol(emb), min(2, ncol(emb)))
  
  # Initialize velocity embedding
  V_emb <- matrix(0, nrow = n_cells, ncol = n_dims)
  rownames(V_emb) <- rownames(emb)
  colnames(V_emb) <- paste0("V_", colnames(emb)[seq_len(n_dims)])
  
  if (direct_pca) {
    # Direct projection via PCA loadings
    V_emb <- project_velocity_pca(seurat_obj, basis, vkey, n_pcs, n_dims)
  } else {
    # Projection via transition matrix
    validate_seurat(seurat_obj, require_graph = TRUE)
    
    # Get transition matrix
    T_mat <- transition_matrix(seurat_obj, 
                                scale = scale, 
                                self_transitions = self_transitions,
                                use_negative = use_negative)
    
    # Project: V_emb = T * emb - emb
    emb_subset <- emb[, seq_len(n_dims), drop = FALSE]
    V_emb <- as.matrix(T_mat %*% emb_subset) - emb_subset
  }
  
  # Autoscale
  if (autoscale) {
    V_emb <- autoscale_velocity_embedding(V_emb, scale)
  }
  
  # Store result
  seurat_obj@misc$velocity_embedding <- V_emb
  seurat_obj@misc$velocity_embedding_params <- list(
    basis = basis,
    scale = scale,
    autoscale = autoscale,
    direct_pca = direct_pca
  )
  
  .vmessage("Velocity embedding computed")
  
  return(seurat_obj)
}

#' Get Velocity Embedding
#'
#' @description Extract velocity embedding from Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @return Matrix of velocity embedding coordinates
#'
#' @export
get_velocity_embedding <- function(seurat_obj) {
  V_emb <- seurat_obj@misc$velocity_embedding
  
  if (is.null(V_emb)) {
    stop("Velocity embedding not found. Run project_velocity_embedding first.", call. = FALSE)
  }
  
  return(V_emb)
}

#' Autoscale Velocity Embedding
#'
#' @description Scale velocity vectors to have consistent magnitudes.
#'
#' @param V_emb Velocity embedding matrix
#' @param scale Overall scale factor
#' @return Scaled velocity embedding
#' @keywords internal
autoscale_velocity_embedding <- function(V_emb, scale = 10) {
  # Compute magnitudes
  mags <- sqrt(rowSums(V_emb^2))
  
  # Scale to median magnitude
  median_mag <- stats::median(mags[mags > 0])
  if (median_mag > 0) {
    V_emb <- V_emb / median_mag
  }
  
  # Apply overall scale
  V_emb <- V_emb / scale
  
  return(V_emb)
}

#' Project Velocity via PCA
#'
#' @description Direct projection of velocity to embedding via PCA loadings.
#'
#' @param seurat_obj Seurat object
#' @param basis Embedding name
#' @param vkey Velocity key
#' @param n_pcs Number of PCs
#' @param n_dims Output dimensions
#' @return Velocity embedding matrix
#' @keywords internal
project_velocity_pca <- function(seurat_obj, basis, vkey, n_pcs, n_dims) {
  
  # Get velocity
  V <- get_velocity(seurat_obj, remove_na = TRUE)
  
  # Get PCA loadings
  if (!"pca" %in% names(seurat_obj@reductions)) {
    stop("PCA not found. Run RunPCA first.", call. = FALSE)
  }
  
  pca_loadings <- Seurat::Loadings(seurat_obj, "pca")
  
  # Subset to common genes
  common_genes <- intersect(rownames(V), rownames(pca_loadings))
  V <- V[common_genes, , drop = FALSE]
  pca_loadings <- pca_loadings[common_genes, , drop = FALSE]
  
  # Limit PCs
  n_pcs <- min(n_pcs, ncol(pca_loadings))
  pca_loadings <- pca_loadings[, seq_len(n_pcs), drop = FALSE]
  
  # Project velocity to PCA space: V_pca = t(V) %*% loadings
  V_pca <- t(make_dense(V)) %*% pca_loadings
  
  # Get embedding
  emb <- get_embedding(seurat_obj, basis)
  n_dims <- min(n_dims, ncol(emb))
  
  # If embedding is UMAP from PCA, we need to project further
  # For simplicity, we use the first n_dims of V_pca
  if (basis == "umap" && ncol(V_pca) >= n_dims) {
    # Use correlation with embedding to estimate projection
    V_emb <- matrix(0, nrow = nrow(V_pca), ncol = n_dims)
    
    for (i in seq_len(n_dims)) {
      # Find PC dimensions most correlated with this embedding dimension
      cors <- apply(V_pca, 2, function(x) stats::cor(x, emb[, i], use = "complete.obs"))
      best_pcs <- order(abs(cors), decreasing = TRUE)[1:min(5, n_pcs)]
      
      # Weighted combination
      weights <- cors[best_pcs]
      weights <- weights / sum(abs(weights))
      V_emb[, i] <- V_pca[, best_pcs, drop = FALSE] %*% weights
    }
  } else {
    V_emb <- V_pca[, seq_len(n_dims), drop = FALSE]
  }
  
  rownames(V_emb) <- colnames(seurat_obj)
  colnames(V_emb) <- paste0("V_", seq_len(n_dims))
  
  return(V_emb)
}

#' Compute Velocity Magnitude
#'
#' @description Compute the magnitude of velocity vectors.
#'
#' @param seurat_obj A Seurat object
#' @param use_embedding Use embedding-projected velocities (default: TRUE)
#'
#' @return Numeric vector of velocity magnitudes
#'
#' @export
velocity_magnitude <- function(seurat_obj, use_embedding = TRUE) {
  if (use_embedding) {
    V <- get_velocity_embedding(seurat_obj)
  } else {
    V <- get_velocity(seurat_obj)
    V <- t(make_dense(V))
  }
  
  mags <- sqrt(rowSums(V^2))
  names(mags) <- rownames(V)
  
  return(mags)
}

#' Compute Velocity Confidence
#'
#' @description Compute confidence scores for velocity estimates.
#'
#' @param seurat_obj A Seurat object with velocity graph
#' @param vkey Velocity key
#' @param n_neighbors Number of neighbors
#'
#' @return Seurat object with confidence scores in metadata
#'
#' @export
velocity_confidence <- function(seurat_obj, vkey = "velocity", n_neighbors = NULL) {
  
  .vmessage("Computing velocity confidence")
  
  # Get velocity graph
  graph <- get_velocity_graph(seurat_obj, negative = FALSE)
  graph_neg <- get_velocity_graph(seurat_obj, negative = TRUE)
  
  # Confidence based on positive vs negative correlations
  pos_sum <- Matrix::rowSums(graph)
  neg_sum <- Matrix::rowSums(graph_neg)
  total_sum <- pos_sum + neg_sum
  
  # Avoid division by zero
  total_sum[total_sum == 0] <- 1
  
  # Confidence = (positive - negative) / total
  confidence <- (pos_sum - neg_sum) / total_sum
  
  # Coherence = max correlation
  coherence <- apply(graph, 1, function(x) {
    if (all(x == 0)) return(0)
    max(x)
  })
  
  # Store in metadata
  seurat_obj@meta.data$velocity_confidence <- confidence
  seurat_obj@meta.data$velocity_coherence <- coherence
  
  .vmessage("Added velocity_confidence and velocity_coherence to metadata")
  
  return(seurat_obj)
}

#' Compute Cell Velocity Length
#'
#' @description Compute the expected velocity length (transition distance).
#'
#' @param seurat_obj A Seurat object
#' @param vkey Velocity key
#' @param basis Embedding basis
#'
#' @return Numeric vector of velocity lengths
#'
#' @export
velocity_length <- function(seurat_obj, vkey = "velocity", basis = "umap") {
  
  # Get velocity embedding
  V_emb <- get_velocity_embedding(seurat_obj)
  
  # Compute lengths
  lengths <- sqrt(rowSums(V_emb^2))
  names(lengths) <- rownames(V_emb)
  
  return(lengths)
}

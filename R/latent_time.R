#' @title Latent Time and Terminal States
#' @description Functions for computing latent time and identifying terminal states.
#' @name latent_time
NULL

#' Compute Shared Latent Time
#'
#' @description Compute a consensus latent time across all cells by aligning 
#' gene-specific times.
#'
#' @param seurat_obj A Seurat object with recovered dynamics
#' @param genes Genes to use (default: NULL uses velocity genes)
#' @param min_corr Minimum correlation for gene inclusion
#' @param normalize_time Whether to normalize to [0,1]
#'
#' @return Numeric vector of latent times
#'
#' @export
compute_shared_time <- function(seurat_obj,
                                  genes = NULL,
                                  min_corr = 0.3,
                                  normalize_time = TRUE) {
  
  # Get latent time matrix
  latent_time <- seurat_obj@misc$latent_time
  
  if (is.null(latent_time)) {
    stop("Latent time not computed. Run recover_dynamics first.", call. = FALSE)
  }
  
  # Select genes
  if (is.null(genes)) {
    fit_params <- get_fit_params(seurat_obj)
    genes <- rownames(fit_params)[!is.na(fit_params$fit_alpha)]
  }
  
  genes <- intersect(genes, rownames(latent_time))
  latent_time <- latent_time[genes, , drop = FALSE]
  
  # Remove genes with too many NAs
  na_frac <- rowMeans(is.na(latent_time))
  latent_time <- latent_time[na_frac < 0.5, , drop = FALSE]
  
  if (nrow(latent_time) < 3) {
    warning("Too few genes with valid latent times")
    return(rep(NA, ncol(latent_time)))
  }
  
  # Compute pairwise correlations to find consistent genes
  n_genes <- nrow(latent_time)
  
  # Reference time: first PCA component of latent times
  latent_time_clean <- latent_time
  latent_time_clean[is.na(latent_time_clean)] <- 0
  
  pca <- stats::prcomp(t(latent_time_clean), center = TRUE, scale. = TRUE)
  reference_time <- pca$x[, 1]
  
  # Flip if negatively correlated with median time
  median_time <- apply(latent_time, 2, stats::median, na.rm = TRUE)
  if (stats::cor(reference_time, median_time, use = "complete.obs") < 0) {
    reference_time <- -reference_time
  }
  
  # Weight genes by correlation with reference
  gene_weights <- apply(latent_time, 1, function(x) {
    stats::cor(x, reference_time, use = "complete.obs")
  })
  
  # Keep positively correlated genes
  good_genes <- gene_weights > min_corr
  
  if (sum(good_genes) < 3) {
    shared_time <- reference_time
  } else {
    # Weighted average
    weights <- pmax(0, gene_weights[good_genes])
    weights <- weights / sum(weights)
    
    shared_time <- apply(latent_time[good_genes, , drop = FALSE], 2, function(x) {
      valid <- !is.na(x)
      if (sum(valid) == 0) return(NA)
      sum(x[valid] * weights[valid]) / sum(weights[valid])
    })
  }
  
  # Normalize
  if (normalize_time) {
    shared_time <- (shared_time - min(shared_time, na.rm = TRUE)) /
                    (max(shared_time, na.rm = TRUE) - min(shared_time, na.rm = TRUE))
  }
  
  names(shared_time) <- colnames(latent_time)
  
  return(shared_time)
}

#' Root Time
#'
#' @description Determine the root time (starting point) for trajectory.
#'
#' @param seurat_obj A Seurat object with latent time
#' @param root_cells Root cell indices or names
#'
#' @return Adjusted latent time with specified root
#' @keywords internal
root_time <- function(seurat_obj, root_cells = NULL) {
  
  latent_time <- seurat_obj@meta.data$latent_time
  
  if (is.null(root_cells)) {
    # Auto-detect root from terminal states
    if (!is.null(seurat_obj@meta.data$root_cells_prob)) {
      root_cells <- which.max(seurat_obj@meta.data$root_cells_prob)
    } else {
      root_cells <- which.min(latent_time)
    }
  }
  
  if (is.character(root_cells)) {
    root_cells <- match(root_cells, colnames(seurat_obj))
  }
  
  # Adjust time so root has time 0
  root_time_val <- mean(latent_time[root_cells], na.rm = TRUE)
  
  # Flip direction if needed
  if (root_time_val > 0.5) {
    latent_time <- 1 - latent_time
  }
  
  # Shift to start at 0
  latent_time <- latent_time - min(latent_time, na.rm = TRUE)
  latent_time <- latent_time / max(latent_time, na.rm = TRUE)
  
  return(latent_time)
}

#' Compute Velocity Genes by Likelihood
#'
#' @description Select velocity genes based on dynamics likelihood.
#'
#' @param seurat_obj Seurat object with recovered dynamics
#' @param min_likelihood Minimum likelihood threshold
#' @param n_top Maximum number of genes
#'
#' @return Character vector of gene names
#' @keywords internal
get_likelihood_genes <- function(seurat_obj, 
                                   min_likelihood = 0.1,
                                   n_top = 1000) {
  
  fit_params <- get_fit_params(seurat_obj)
  
  if (!"fit_likelihood" %in% colnames(fit_params)) {
    return(rownames(fit_params)[!is.na(fit_params$fit_alpha)])
  }
  
  likelihood <- fit_params$fit_likelihood
  valid <- !is.na(likelihood) & (likelihood >= min_likelihood)
  
  genes <- rownames(fit_params)[valid]
  
  # Sort by likelihood
  if (length(genes) > n_top) {
    order_idx <- order(likelihood[valid], decreasing = TRUE)
    genes <- genes[order_idx[seq_len(n_top)]]
  }
  
  return(genes)
}

#' Compute Divergence
#'
#' @description Compute velocity divergence for each cell.
#'
#' @param seurat_obj A Seurat object with velocity
#'
#' @return Numeric vector of divergence values
#'
#' @export
velocity_divergence <- function(seurat_obj) {
  
  # Get velocity embedding
  V_emb <- get_velocity_embedding(seurat_obj)
  
  # Get embedding
  basis <- seurat_obj@misc$velocity_embedding_params$basis
  emb <- get_embedding(seurat_obj, basis)
  emb <- emb[, seq_len(ncol(V_emb)), drop = FALSE]
  
  # Get neighbors
  conn <- get_connectivities(seurat_obj)
  
  n_cells <- nrow(emb)
  divergence <- numeric(n_cells)
  
  for (i in seq_len(n_cells)) {
    # Get neighbors
    neighbors <- which(conn[i, ] > 0)
    if (length(neighbors) < 2) {
      divergence[i] <- 0
      next
    }
    
    # Compute divergence: how much do neighbors' velocities point away
    v_i <- V_emb[i, ]
    
    neighbor_divergence <- sapply(neighbors, function(j) {
      # Direction from i to j
      d_ij <- emb[j, ] - emb[i, ]
      d_norm <- sqrt(sum(d_ij^2))
      if (d_norm < 1e-10) return(0)
      d_ij <- d_ij / d_norm
      
      # Velocity at j
      v_j <- V_emb[j, ]
      
      # Divergence contribution: how aligned is v_j with outward direction
      sum(v_j * d_ij)
    })
    
    divergence[i] <- mean(neighbor_divergence)
  }
  
  names(divergence) <- rownames(emb)
  return(divergence)
}

#' Compute RNA Velocity Confidence Score
#'
#' @description Compute overall confidence scores for velocity estimates.
#'
#' @param seurat_obj A Seurat object
#' @param method Method for confidence: "graph" or "genes"
#'
#' @return Seurat object with confidence in metadata
#'
#' @export
compute_velocity_confidence <- function(seurat_obj, method = "graph") {
  
  if (method == "graph") {
    seurat_obj <- velocity_confidence(seurat_obj)
  } else {
    # Gene-based confidence
    fit_params <- get_fit_params(seurat_obj)
    velocity <- get_velocity(seurat_obj)
    
    # Confidence based on how many genes have good velocity
    n_good_genes <- colSums(!is.na(velocity))
    confidence <- n_good_genes / nrow(velocity)
    
    seurat_obj@meta.data$velocity_confidence <- confidence
  }
  
  return(seurat_obj)
}

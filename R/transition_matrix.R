#' @title Transition Matrix Computation
#' @description Functions for computing cell-to-cell transition probabilities.
#' @name transition_matrix
NULL

#' Compute Transition Matrix
#'
#' @description Compute the transition probability matrix from the velocity graph.
#' The transition matrix represents the probability of transitioning from one 
#' cell state to another based on velocity.
#'
#' @param seurat_obj A Seurat object with computed velocity graph
#' @param scale Scale parameter for Gaussian kernel (default: 10)
#' @param self_transitions Include self-transitions (default: TRUE)
#' @param use_negative Include negative correlations (default: FALSE)
#' @param backward Compute backward transitions (default: FALSE)
#' @param weight_diffusion Weight by diffusion distances (default: FALSE)
#'
#' @return Sparse transition matrix (row-stochastic)
#'
#' @export
transition_matrix <- function(seurat_obj,
                               scale = 10,
                               self_transitions = TRUE,
                               use_negative = FALSE,
                               backward = FALSE,
                               weight_diffusion = FALSE) {
  
  # Get velocity graph
  graph <- get_velocity_graph(seurat_obj, negative = FALSE)
  
  if (use_negative) {
    graph_neg <- get_velocity_graph(seurat_obj, negative = TRUE)
    graph <- graph - graph_neg
  }
  
  if (backward) {
    graph <- Matrix::t(graph)
  }
  
  n_cells <- nrow(graph)
  
  # Apply exponential kernel: T = exp(scale * graph)
  # For sparse matrix, apply to non-zero elements
  T_mat <- graph
  T_mat@x <- expm1(scale * T_mat@x)  # expm1(x) = exp(x) - 1
  T_mat@x <- T_mat@x + 1  # Add 1 back, so result is exp(scale * x)
  
  # Weight by diffusion distances if requested
  if (weight_diffusion) {
    distances <- get_distances(seurat_obj)
    if (!is.null(distances)) {
      # Gaussian kernel on distances
      sigma <- stats::median(distances@x[distances@x > 0])
      weights <- Matrix::sparseMatrix(
        i = distances@i + 1,
        j = rep(seq_len(n_cells), diff(distances@p)),
        x = exp(-distances@x^2 / (2 * sigma^2)),
        dims = c(n_cells, n_cells)
      )
      T_mat <- T_mat * weights
    }
  }
  
  # Self-transitions
  if (self_transitions) {
    # Set diagonal to 1
    diag(T_mat) <- 1
  } else {
    # Zero diagonal
    diag(T_mat) <- 0
  }
  
  # Row-normalize to get transition probabilities
  row_sums <- Matrix::rowSums(T_mat)
  row_sums[row_sums == 0] <- 1
  T_mat <- T_mat / row_sums
  
  return(as(T_mat, "CsparseMatrix"))
}

#' Get Transition Matrix
#'
#' @description Get or compute transition matrix from Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param recompute Whether to recompute even if cached
#' @param ... Additional arguments passed to transition_matrix
#'
#' @return Sparse transition matrix
#'
#' @export
get_transition_matrix <- function(seurat_obj, recompute = FALSE, ...) {
  
  if (!recompute && !is.null(seurat_obj@misc$transition_matrix)) {
    return(seurat_obj@misc$transition_matrix)
  }
  
  T_mat <- transition_matrix(seurat_obj, ...)
  
  return(T_mat)
}

#' Store Transition Matrix
#'
#' @description Compute and store transition matrix in Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param ... Arguments passed to transition_matrix
#'
#' @return Modified Seurat object
#'
#' @export
store_transition_matrix <- function(seurat_obj, ...) {
  T_mat <- transition_matrix(seurat_obj, ...)
  seurat_obj@misc$transition_matrix <- T_mat
  return(seurat_obj)
}

#' Velocity Pseudotime
#'
#' @description Compute pseudotime based on velocity-derived transition matrix.
#'
#' @param seurat_obj A Seurat object with velocity graph
#' @param root_cells Root cells (starting point for pseudotime)
#' @param end_points End points (optional terminal cells)
#' @param n_dcs Number of diffusion components (default: 10)
#' @param use_velocity_graph Use velocity graph instead of transition matrix
#'
#' @return Modified Seurat object with velocity_pseudotime in metadata
#'
#' @export
velocity_pseudotime <- function(seurat_obj,
                                 root_cells = NULL,
                                 end_points = NULL,
                                 n_dcs = 10,
                                 use_velocity_graph = TRUE) {
  
  .vmessage("Computing velocity pseudotime")
  
  n_cells <- ncol(seurat_obj)
  cell_names <- colnames(seurat_obj)
  
  # Get transition matrix
  T_mat <- transition_matrix(seurat_obj, self_transitions = TRUE)
  
  # Symmetrize for diffusion
  T_sym <- (T_mat + Matrix::t(T_mat)) / 2
  
  # Compute diffusion components
  # Using eigen decomposition of transition matrix
  n_dcs <- min(n_dcs, n_cells - 1)
  
  # Convert to dense for eigen (for small datasets) or use irlba
  if (n_cells < 2000) {
    T_dense <- as.matrix(T_sym)
    eig <- eigen(T_dense, symmetric = TRUE)
    dc <- eig$vectors[, 2:(n_dcs + 1), drop = FALSE]
  } else {
    # Use sparse eigenvalue decomposition
    eig <- irlba::partial_eigen(T_sym, n = n_dcs + 1, symmetric = TRUE)
    dc <- eig$vectors[, 2:(n_dcs + 1), drop = FALSE]
  }
  
  rownames(dc) <- cell_names
  
  # Compute pseudotime from first diffusion component
  if (is.null(root_cells)) {
    # Find root automatically based on velocity graph
    graph <- get_velocity_graph(seurat_obj)
    
    # Root cells have high outgoing and low incoming transitions
    outgoing <- Matrix::rowSums(graph)
    incoming <- Matrix::colSums(graph)
    
    root_score <- outgoing - incoming
    root_cells <- which.max(root_score)
  } else if (is.character(root_cells)) {
    root_cells <- match(root_cells, cell_names)
    root_cells <- root_cells[!is.na(root_cells)]
  }
  
  # Pseudotime: distance from root in diffusion space
  if (length(root_cells) == 1) {
    root_dc <- dc[root_cells, , drop = FALSE]
  } else {
    root_dc <- colMeans(dc[root_cells, , drop = FALSE])
  }
  
  # Euclidean distance in diffusion space
  distances <- sqrt(rowSums(sweep(dc, 2, root_dc, "-")^2))
  
  # Scale to [0, 1]
  pseudotime <- (distances - min(distances)) / (max(distances) - min(distances))
  
  # Store results
  seurat_obj@meta.data$velocity_pseudotime <- pseudotime
  seurat_obj@misc$diffusion_components <- dc
  seurat_obj@misc$pseudotime_root <- root_cells
  
  .vmessage("Added velocity_pseudotime to metadata")
  
  return(seurat_obj)
}

#' Terminal States
#'
#' @description Identify terminal states based on velocity graph.
#'
#' @param seurat_obj A Seurat object with velocity graph
#' @param groupby Column name for cell grouping (optional)
#' @param self_transitions Include self-transitions
#'
#' @return Modified Seurat object with terminal state information
#'
#' @export
terminal_states <- function(seurat_obj,
                             groupby = NULL,
                             self_transitions = FALSE) {
  
  .vmessage("Identifying terminal states")
  
  # Get velocity graph
  graph <- get_velocity_graph(seurat_obj)
  n_cells <- nrow(graph)
  cell_names <- colnames(seurat_obj)
  
  # Compute end point probability
  # Terminal cells have low outgoing velocity correlations
  
  # Forward transitions
  T_fwd <- transition_matrix(seurat_obj, self_transitions = self_transitions)
  
  # Backward transitions
  T_bwd <- transition_matrix(seurat_obj, self_transitions = self_transitions, backward = TRUE)
  
  # Absorption probabilities via iterating transition matrix
  # Initialize with uniform
  n_iter <- 100
  fwd_prob <- rep(1/n_cells, n_cells)
  bwd_prob <- rep(1/n_cells, n_cells)
  
  for (i in seq_len(n_iter)) {
    fwd_prob <- as.vector(T_fwd %*% fwd_prob)
    bwd_prob <- as.vector(T_bwd %*% bwd_prob)
  }
  
  # Normalize
  fwd_prob <- fwd_prob / max(fwd_prob)
  bwd_prob <- bwd_prob / max(bwd_prob)
  
  # End points have high forward probability (absorbing)
  # Root points have high backward probability
  
  end_points <- fwd_prob
  root_cells <- bwd_prob
  
  # Store
  seurat_obj@meta.data$end_points <- end_points
  seurat_obj@meta.data$root_cells_prob <- root_cells
  
  # Identify likely terminal clusters
  if (!is.null(groupby) && groupby %in% colnames(seurat_obj@meta.data)) {
    groups <- seurat_obj@meta.data[[groupby]]
    
    # Average terminal probability per group
    terminal_by_group <- tapply(end_points, groups, mean)
    root_by_group <- tapply(root_cells, groups, mean)
    
    seurat_obj@misc$terminal_groups <- sort(terminal_by_group, decreasing = TRUE)
    seurat_obj@misc$root_groups <- sort(root_by_group, decreasing = TRUE)
  }
  
  .vmessage("Added end_points and root_cells_prob to metadata")
  
  return(seurat_obj)
}

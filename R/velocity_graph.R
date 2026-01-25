#' @title Velocity Graph Construction
#' @description Build velocity graph and transition probability matrices.
#' The velocity graph represents the predicted direction of each cell based on
#' its RNA velocity vector.
#' @name velocity_graph
NULL

#' Build Velocity Graph
#'
#' @description Construct velocity graph using cosine similarity between 
#' velocity vectors and cell-cell displacement vectors.
#'
#' @param object Seurat object with velocity computed
#' @param n_neighbors Number of neighbors for graph construction
#' @param n_recurse_neighbors Number of times to recurse neighbors
#' @param mode Graph computation mode: "deterministic", "stochastic", or "dynamical"
#' @param approx Use approximate computation
#' @param sqrt_transform Square root transform of graph weights
#' @param n_cores Number of cores for parallel computation
#' @param verbose Print progress
#'
#' @return Seurat object with velocity_graph in misc$scVeloR$velocity_graph
#' @export
velocity_graph <- function(object,
                           n_neighbors = 30,
                           n_recurse_neighbors = 2,
                           mode = NULL,
                           approx = TRUE,
                           sqrt_transform = TRUE,
                           n_cores = 1,
                           verbose = TRUE) {
  
  if (is.null(object@misc$scVeloR$velocity$velocity_s)) {
    stop("Run velocity() first")
  }
  
  if (verbose) {
    message("Building velocity graph...")
  }
  
  # Get velocity matrix
  velocity <- object@misc$scVeloR$velocity$velocity_s
  velocity_genes <- object@misc$scVeloR$velocity$velocity_genes
  
  # Get spliced data
  data <- get_velocity_data(object)
  Ms <- data$Ms
  
  # Use velocity genes only
  velocity <- velocity[, velocity_genes, drop = FALSE]
  Ms <- Ms[, velocity_genes, drop = FALSE]
  
  n_cells <- nrow(Ms)
  
  # Get or compute neighbor graph
  if (!is.null(object@misc$scVeloR$neighbors)) {
    indices <- object@misc$scVeloR$neighbors$indices
    distances <- object@misc$scVeloR$neighbors$distances
    connectivities <- object@misc$scVeloR$neighbors$connectivities
    
    if (ncol(indices) < n_neighbors) {
      message("  Recomputing neighbors with more neighbors...")
      object <- compute_neighbors(object, n_neighbors = n_neighbors, verbose = FALSE)
      indices <- object@misc$scVeloR$neighbors$indices
      distances <- object@misc$scVeloR$neighbors$distances
      connectivities <- object@misc$scVeloR$neighbors$connectivities
    }
  } else {
    object <- compute_neighbors(object, n_neighbors = n_neighbors, verbose = FALSE)
    indices <- object@misc$scVeloR$neighbors$indices
    distances <- object@misc$scVeloR$neighbors$distances
    connectivities <- object@misc$scVeloR$neighbors$connectivities
  }
  
  # Extend neighbors by recursion
  if (n_recurse_neighbors > 0) {
    indices_ext <- recurse_neighbors(indices, n_recurse = n_recurse_neighbors)
  } else {
    indices_ext <- indices
  }
  
  # Compute cosine correlation between velocity and displacement
  if (verbose) {
    message("  Computing velocity correlations...")
  }
  
  # Use C++ implementation if available
  if (requireNamespace("Rcpp", quietly = TRUE)) {
    graph <- compute_cosine_correlation_graph(
      velocity = as.matrix(velocity),
      expression = as.matrix(Ms),
      indices = indices_ext,
      n_cores = n_cores
    )
  } else {
    graph <- compute_velocity_graph_r(
      velocity = velocity,
      expression = Ms,
      indices = indices_ext
    )
  }
  
  # Apply transformations
  if (sqrt_transform) {
    graph@x <- sqrt(abs(graph@x)) * sign(graph@x)
  }
  
  # Build transition matrix
  if (verbose) {
    message("  Computing transition matrix...")
  }
  
  transition_matrix <- compute_transition_matrix(
    graph = graph,
    self_transitions = TRUE,
    scale_by_distance = TRUE,
    distances = distances,
    indices = indices
  )
  
  # Store results
  if (is.null(object@misc$scVeloR$velocity_graph)) {
    object@misc$scVeloR$velocity_graph <- list()
  }
  
  object@misc$scVeloR$velocity_graph$graph <- graph
  object@misc$scVeloR$velocity_graph$transition_matrix <- transition_matrix
  object@misc$scVeloR$velocity_graph$n_neighbors <- n_neighbors
  
  # Also store in velocity for convenience
  object@misc$scVeloR$velocity$transition_matrix <- transition_matrix
  
  if (verbose) {
    message("Done. Velocity graph stored in object@misc$scVeloR$velocity_graph")
  }
  
  object
}

#' Compute Velocity Graph in R
#'
#' @description Pure R implementation of velocity graph computation.
#'
#' @param velocity Velocity matrix
#' @param expression Expression matrix
#' @param indices Neighbor indices matrix
#'
#' @return Sparse graph matrix
#' @keywords internal
compute_velocity_graph_r <- function(velocity, expression, indices) {
  
  n_cells <- nrow(velocity)
  n_neighbors <- ncol(indices)
  
  # Initialize sparse matrix components
  i_idx <- integer(0)
  j_idx <- integer(0)
  x_val <- numeric(0)
  
  for (cell in seq_len(n_cells)) {
    # Get velocity vector for this cell
    v <- velocity[cell, ]
    v_norm <- sqrt(sum(v^2))
    
    if (v_norm == 0 || !is.finite(v_norm)) next
    
    # Get neighbors
    neighbors <- indices[cell, ]
    neighbors <- neighbors[neighbors > 0 & neighbors <= n_cells]
    
    if (length(neighbors) == 0) next
    
    for (j in neighbors) {
      # Displacement from cell to neighbor
      dx <- expression[j, ] - expression[cell, ]
      dx_norm <- sqrt(sum(dx^2))
      
      if (dx_norm == 0 || !is.finite(dx_norm)) next
      
      # Cosine similarity
      cos_sim <- sum(v * dx) / (v_norm * dx_norm)
      
      # Store positive correlations
      if (cos_sim > 0) {
        i_idx <- c(i_idx, cell)
        j_idx <- c(j_idx, j)
        x_val <- c(x_val, cos_sim)
      }
    }
  }
  
  # Create sparse matrix
  sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = x_val,
    dims = c(n_cells, n_cells)
  )
}

#' Compute Cosine Correlation Graph Using C++
#'
#' @description Compute velocity graph using optimized C++ code.
#'
#' @param velocity Velocity matrix
#' @param expression Expression matrix
#' @param indices Neighbor indices matrix
#' @param n_cores Number of cores
#'
#' @return Sparse graph matrix
#' @keywords internal
compute_cosine_correlation_graph <- function(velocity, expression, indices, n_cores = 1) {
  
  n_cells <- nrow(velocity)
  n_neighbors <- ncol(indices)
  
  # Check if our C++ function is available
  if (exists("compute_cosine_similarity_cpp")) {
    result <- compute_cosine_similarity_cpp(velocity, expression, indices)
    
    # Convert to sparse matrix
    sparseMatrix(
      i = result$i + 1L,  # C++ uses 0-based indexing
      j = result$j + 1L,
      x = result$x,
      dims = c(n_cells, n_cells)
    )
  } else {
    # Fall back to R implementation
    compute_velocity_graph_r(velocity, expression, indices)
  }
}

#' Recurse Neighbors
#'
#' @description Extend neighbor graph by including neighbors of neighbors.
#'
#' @param indices Original neighbor indices
#' @param n_recurse Number of recursions
#'
#' @return Extended neighbor indices matrix
#' @keywords internal
recurse_neighbors <- function(indices, n_recurse = 1) {
  
  if (n_recurse <= 0) return(indices)
  
  n_cells <- nrow(indices)
  n_neighbors <- ncol(indices)
  
  # Convert to sets for easier manipulation
  neighbor_sets <- lapply(seq_len(n_cells), function(i) {
    unique(indices[i, ])
  })
  
  for (r in seq_len(n_recurse)) {
    new_neighbor_sets <- lapply(seq_len(n_cells), function(i) {
      current <- neighbor_sets[[i]]
      extended <- current
      
      for (j in current) {
        if (j > 0 && j <= n_cells) {
          extended <- union(extended, neighbor_sets[[j]])
        }
      }
      
      # Remove self
      extended <- setdiff(extended, i)
      extended
    })
    neighbor_sets <- new_neighbor_sets
  }
  
  # Convert back to matrix format
  max_neighbors <- max(sapply(neighbor_sets, length))
  
  indices_ext <- matrix(0L, n_cells, max_neighbors)
  
  for (i in seq_len(n_cells)) {
    neighbors <- neighbor_sets[[i]]
    if (length(neighbors) > 0) {
      indices_ext[i, seq_along(neighbors)] <- neighbors
    }
  }
  
  indices_ext
}

#' Compute Transition Matrix
#'
#' @description Compute transition probability matrix from velocity graph.
#'
#' @param graph Velocity graph (sparse matrix)
#' @param self_transitions Include self-transitions
#' @param scale_by_distance Scale by distance
#' @param distances Distance matrix
#' @param indices Neighbor indices
#'
#' @return Transition probability matrix
#' @keywords internal
compute_transition_matrix <- function(graph,
                                       self_transitions = TRUE,
                                       scale_by_distance = FALSE,
                                       distances = NULL,
                                       indices = NULL) {
  
  n_cells <- nrow(graph)
  
  # Make a copy
  T <- graph
  
  # Scale by distance if requested
  if (scale_by_distance && !is.null(distances) && !is.null(indices)) {
    # Create distance-based scaling
    scale_factor <- matrix(1, n_cells, n_cells)
    
    for (i in seq_len(n_cells)) {
      neighbors <- indices[i, ]
      dists <- distances[i, ]
      
      for (k in seq_along(neighbors)) {
        j <- neighbors[k]
        if (j > 0 && j <= n_cells) {
          # Gaussian kernel
          scale_factor[i, j] <- exp(-dists[k]^2 / (2 * median(dists[dists > 0])^2))
        }
      }
    }
    
    T <- T * scale_factor
  }
  
  # Add self-transitions based on velocity confidence
  if (self_transitions) {
    row_sums <- rowSums(T)
    max_sum <- max(row_sums[is.finite(row_sums)])
    
    self_trans <- 1 - row_sums / max_sum
    self_trans[!is.finite(self_trans)] <- 0
    self_trans <- pmax(self_trans, 0)
    
    diag(T) <- self_trans
  }
  
  # Row-normalize
  row_sums <- rowSums(T)
  row_sums[row_sums == 0] <- 1
  
  T <- T / row_sums
  
  T
}

#' Velocity Confidence
#'
#' @description Compute confidence scores for velocity estimates.
#'
#' @param object Seurat object with velocity graph
#'
#' @return Seurat object with velocity_confidence in metadata
#' @export
velocity_confidence <- function(object) {
  
  if (is.null(object@misc$scVeloR$velocity_graph$graph)) {
    stop("Run velocity_graph() first")
  }
  
  graph <- object@misc$scVeloR$velocity_graph$graph
  
  # Confidence = mean correlation with neighbors
  n_cells <- nrow(graph)
  
  confidence <- numeric(n_cells)
  
  for (i in seq_len(n_cells)) {
    row <- graph[i, ]
    nonzero <- which(row > 0)
    
    if (length(nonzero) > 0) {
      confidence[i] <- mean(row[nonzero])
    }
  }
  
  # Scale to [0, 1]
  confidence <- (confidence - min(confidence)) / (max(confidence) - min(confidence) + 1e-10)
  
  object@meta.data$velocity_confidence <- confidence
  
  object
}

#' Velocity Coherence
#'
#' @description Compute coherence of velocity with neighboring cells.
#'
#' @param object Seurat object with velocity
#' @param n_neighbors Number of neighbors
#'
#' @return Seurat object with velocity_coherence in metadata
#' @export
velocity_coherence <- function(object, n_neighbors = 30) {
  
  if (is.null(object@misc$scVeloR$velocity$velocity_s)) {
    stop("Run velocity() first")
  }
  
  velocity <- object@misc$scVeloR$velocity$velocity_s
  velocity_genes <- object@misc$scVeloR$velocity$velocity_genes
  velocity <- velocity[, velocity_genes, drop = FALSE]
  
  # Get neighbors
  if (is.null(object@misc$scVeloR$neighbors)) {
    object <- compute_neighbors(object, n_neighbors = n_neighbors, verbose = FALSE)
  }
  
  indices <- object@misc$scVeloR$neighbors$indices
  
  n_cells <- nrow(velocity)
  coherence <- numeric(n_cells)
  
  for (i in seq_len(n_cells)) {
    v_i <- velocity[i, ]
    v_i_norm <- sqrt(sum(v_i^2))
    
    if (v_i_norm == 0) next
    
    neighbors <- indices[i, ]
    neighbors <- neighbors[neighbors > 0 & neighbors <= n_cells]
    
    if (length(neighbors) == 0) next
    
    correlations <- numeric(length(neighbors))
    
    for (k in seq_along(neighbors)) {
      j <- neighbors[k]
      v_j <- velocity[j, ]
      v_j_norm <- sqrt(sum(v_j^2))
      
      if (v_j_norm > 0) {
        correlations[k] <- sum(v_i * v_j) / (v_i_norm * v_j_norm)
      }
    }
    
    coherence[i] <- mean(correlations)
  }
  
  # Scale
  coherence <- (coherence - min(coherence)) / (max(coherence) - min(coherence) + 1e-10)
  
  object@meta.data$velocity_coherence <- coherence
  
  object
}

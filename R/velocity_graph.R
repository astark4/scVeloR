#' @title Velocity Graph Computation
#' @description Functions for constructing the velocity graph based on cosine similarity.
#' @name velocity_graph
NULL

#' Compute Velocity Graph
#'
#' @description Compute the velocity graph, which represents the transition 
#' probabilities between cells based on velocity vectors.
#'
#' @param seurat_obj A Seurat object with computed velocity
#' @param vkey Key for velocity layer (default: "velocity")
#' @param xkey Key for expression layer (default: "Ms")
#' @param n_neighbors Number of neighbors (default: NULL uses stored value)
#' @param n_recurse_neighbors Recursion depth for neighbor expansion (default: 0)
#' @param sqrt_transform Apply sqrt transform to velocities (default: TRUE)
#' @param gene_subset Genes to use (default: NULL uses velocity genes)
#' @param n_jobs Number of parallel jobs (default: 1)
#'
#' @return Modified Seurat object with velocity_graph in misc
#'
#' @export
compute_velocity_graph <- function(seurat_obj,
                                    vkey = "velocity",
                                    xkey = "Ms",
                                    n_neighbors = NULL,
                                    n_recurse_neighbors = 0L,
                                    sqrt_transform = TRUE,
                                    gene_subset = NULL,
                                    n_jobs = 1L) {
  
  .vmessage("Computing velocity graph")
  
  # Validate
  validate_seurat(seurat_obj, require_velocity = TRUE)
  
  # Get velocity and expression
  V <- .get_layer_data(seurat_obj, vkey)
  X <- .get_layer_data(seurat_obj, xkey)
  
  # Ensure same dimensions
  common_genes <- intersect(rownames(V), rownames(X))
  V <- V[common_genes, , drop = FALSE]
  X <- X[common_genes, , drop = FALSE]
  
  # Subset genes if specified
  if (!is.null(gene_subset)) {
    gene_subset <- intersect(gene_subset, common_genes)
    V <- V[gene_subset, , drop = FALSE]
    X <- X[gene_subset, , drop = FALSE]
  }
  
  # Remove NA genes
  na_mask <- Matrix::rowSums(is.na(V)) < ncol(V)
  V <- V[na_mask, , drop = FALSE]
  X <- X[na_mask, , drop = FALSE]
  
  .vmessage("Using ", nrow(V), " genes for velocity graph")
  
  # Convert to dense and transpose (cells x genes)
  V_t <- t(make_dense(V))
  X_t <- t(make_dense(X))
  
  n_cells <- nrow(V_t)
  n_genes <- ncol(V_t)
  
  # Get neighbor information
  if (is.null(n_neighbors)) {
    n_neighbors <- get_n_neighbors(seurat_obj)
    if (n_neighbors == 0) {
      n_neighbors <- min(30, n_cells - 1)
    }
  }
  
  # Get neighbor indices
  neighbor_indices <- get_neighbor_indices(seurat_obj)
  
  # Expand neighbors if requested
  if (n_recurse_neighbors > 0) {
    indices_list <- lapply(seq_len(n_cells), function(i) {
      get_iterative_neighbors(seurat_obj, i, n_recurse = n_recurse_neighbors + 1)
    })
  } else {
    # Convert matrix to list of 0-based indices
    indices_list <- lapply(seq_len(n_cells), function(i) {
      as.integer(neighbor_indices[i, ] - 1L)
    })
  }
  
  # Compute velocity graph using C++
  .vmessage("Computing cosine correlations...")
  
  result <- compute_velocity_graph_batch_cpp(
    X_t, V_t, indices_list, sqrt_transform
  )
  
  # Create sparse matrix
  vals <- result$vals
  rows <- result$rows + 1L  # Convert to 1-based
  cols <- result$cols + 1L
  
  # Split into positive and negative
  pos_mask <- vals > 0
  neg_mask <- vals < 0
  
  # Positive graph
  if (sum(pos_mask) > 0) {
    graph_pos <- Matrix::sparseMatrix(
      i = rows[pos_mask],
      j = cols[pos_mask],
      x = vals[pos_mask],
      dims = c(n_cells, n_cells)
    )
  } else {
    graph_pos <- Matrix::sparseMatrix(i = integer(0), j = integer(0), 
                                      x = numeric(0), dims = c(n_cells, n_cells))
  }
  
  # Negative graph
  if (sum(neg_mask) > 0) {
    graph_neg <- Matrix::sparseMatrix(
      i = rows[neg_mask],
      j = cols[neg_mask],
      x = -vals[neg_mask],  # Store absolute values
      dims = c(n_cells, n_cells)
    )
  } else {
    graph_neg <- Matrix::sparseMatrix(i = integer(0), j = integer(0), 
                                      x = numeric(0), dims = c(n_cells, n_cells))
  }
  
  # Set names
  cell_names <- colnames(seurat_obj)
  rownames(graph_pos) <- cell_names
  colnames(graph_pos) <- cell_names
  rownames(graph_neg) <- cell_names
  colnames(graph_neg) <- cell_names
  
  # Store in Seurat object
  seurat_obj@misc$velocity_graph <- graph_pos
  seurat_obj@misc$velocity_graph_neg <- graph_neg
  
  seurat_obj@misc$velocity_graph_params <- list(
    n_neighbors = n_neighbors,
    sqrt_transform = sqrt_transform,
    n_genes = n_genes
  )
  
  .vmessage("Velocity graph computed with ", sum(graph_pos > 0), " positive edges")
  
  return(seurat_obj)
}

#' Get Velocity Graph
#'
#' @description Extract velocity graph from Seurat object.
#'
#' @param seurat_obj A Seurat object
#' @param negative Whether to return negative graph (default: FALSE)
#'
#' @return Sparse velocity graph matrix
#'
#' @export
get_velocity_graph <- function(seurat_obj, negative = FALSE) {
  if (negative) {
    graph <- seurat_obj@misc$velocity_graph_neg
  } else {
    graph <- seurat_obj@misc$velocity_graph
  }
  
  if (is.null(graph)) {
    stop("Velocity graph not found. Run compute_velocity_graph first.", call. = FALSE)
  }
  
  return(graph)
}

#' Cosine Similarity
#'
#' @description Compute cosine similarity between two vectors.
#'
#' @param a First vector
#' @param b Second vector
#' @return Cosine similarity value
#' @keywords internal
cosine_similarity <- function(a, b) {
  norm_a <- sqrt(sum(a^2))
  norm_b <- sqrt(sum(b^2))
  
  if (norm_a < 1e-10 || norm_b < 1e-10) {
    return(0)
  }
  
  sum(a * b) / (norm_a * norm_b)
}

#' Compute Velocity Correlation
#'
#' @description Compute correlation between velocity and cell-cell expression differences.
#'
#' @param dX Matrix of expression differences (n_neighbors x n_genes)
#' @param v Velocity vector (length n_genes)
#' @return Vector of correlations
#' @keywords internal
velocity_correlation <- function(dX, v) {
  # Use Rcpp function
  cosine_correlation_cpp(dX, v)
}

#' Subset Velocity Graph
#'
#' @description Subset velocity graph to specific cells.
#'
#' @param seurat_obj A Seurat object
#' @param cells Cell names or indices
#' @param negative Whether to subset negative graph
#' @return Subsetted sparse matrix
#' @keywords internal
subset_velocity_graph <- function(seurat_obj, cells, negative = FALSE) {
  graph <- get_velocity_graph(seurat_obj, negative)
  
  if (is.character(cells)) {
    cells <- match(cells, rownames(graph))
    cells <- cells[!is.na(cells)]
  }
  
  graph[cells, cells, drop = FALSE]
}

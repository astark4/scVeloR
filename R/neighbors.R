#' @title Neighbor Graph Computation
#' @description Functions for computing k-nearest neighbor graphs.
#' @name neighbors
NULL

#' Compute Neighbor Graph
#'
#' @description Compute k-nearest neighbors graph from PCA or other embeddings.
#'
#' @param seurat_obj A Seurat object with dimensional reduction
#' @param n_neighbors Number of neighbors (default: 30)
#' @param n_pcs Number of principal components to use (default: 30)
#' @param use_rep Representation to use ("pca", "umap", or key in reductions)
#' @param method KNN method: "fnn" (default), "annoy"
#' @param metric Distance metric: "euclidean" (default), "cosine", "manhattan"
#'
#' @return Modified Seurat object with neighbor information in misc
#'
#' @export
compute_neighbors <- function(seurat_obj,
                               n_neighbors = 30,
                               n_pcs = 30,
                               use_rep = "pca",
                               method = "fnn",
                               metric = "euclidean") {
  
  .vmessage("Computing neighbors with ", n_neighbors, " neighbors")
  
  # Get embedding
  if (use_rep == "pca") {
    if (!"pca" %in% names(seurat_obj@reductions)) {
      stop("PCA not found. Run Seurat::RunPCA first.", call. = FALSE)
    }
    X <- Seurat::Embeddings(seurat_obj, "pca")
    if (ncol(X) > n_pcs) {
      X <- X[, seq_len(n_pcs), drop = FALSE]
    }
  } else if (use_rep %in% names(seurat_obj@reductions)) {
    X <- Seurat::Embeddings(seurat_obj, use_rep)
  } else {
    stop("Representation '", use_rep, "' not found", call. = FALSE)
  }
  
  n_cells <- nrow(X)
  
  # Ensure we don't request more neighbors than cells
  n_neighbors <- min(n_neighbors, n_cells - 1)
  
  # Compute KNN
  if (method == "fnn") {
    knn_result <- FNN::get.knn(X, k = n_neighbors)
    indices <- knn_result$nn.index
    distances <- knn_result$nn.dist
  } else {
    stop("Method '", method, "' not yet implemented", call. = FALSE)
  }
  
  # Convert to 0-based for C++ compatibility
  indices_0based <- indices - 1L
  
  # Create sparse distance matrix
  dist_sparse <- indices_to_sparse_cpp(
    indices_0based,
    distances,
    n_cells,
    mode = "distances"
  )
  
  # Create sparse connectivity matrix
  conn_sparse <- indices_to_sparse_cpp(
    indices_0based,
    distances,
    n_cells,
    mode = "connectivities"
  )
  
  # Compute UMAP-style connectivities
  conn_umap <- compute_connectivities_umap_cpp(dist_sparse)
  
  # Set diagonal to 1
  diag(conn_umap) <- 1
  
  # Row-normalize
  conn_norm <- normalize_connectivity_cpp(conn_umap)
  
  # Store in Seurat object
  seurat_obj@misc$neighbors <- list(
    indices = indices,
    distances = distances,
    distances_sparse = dist_sparse,
    connectivities = conn_norm,
    connectivities_raw = conn_umap,
    params = list(
      n_neighbors = n_neighbors,
      n_pcs = n_pcs,
      use_rep = use_rep,
      method = method,
      metric = metric
    )
  )
  
  .vmessage("Added neighbor graph with ", n_neighbors, " neighbors")
  
  return(seurat_obj)
}

#' Get Connectivities Matrix
#'
#' @description Extract the connectivity matrix from a Seurat object.
#'
#' @param seurat_obj A Seurat object with computed neighbors
#' @param mode "connectivities" or "distances" (default: "connectivities")
#' @param n_neighbors Optional: select only top n neighbors
#' @param recurse_neighbors Whether to include neighbors of neighbors
#'
#' @return Sparse connectivity matrix
#'
#' @export
get_connectivities <- function(seurat_obj,
                                mode = "connectivities",
                                n_neighbors = NULL,
                                recurse_neighbors = FALSE) {
  
  if (is.null(seurat_obj@misc$neighbors)) {
    stop("Neighbors not computed. Run compute_neighbors first.", call. = FALSE)
  }
  
  neighbors <- seurat_obj@misc$neighbors
  
  if (mode == "connectivities") {
    conn <- neighbors$connectivities
  } else if (mode == "distances") {
    conn <- neighbors$distances_sparse
  } else {
    stop("Mode must be 'connectivities' or 'distances'", call. = FALSE)
  }
  
  # Select top n neighbors if specified
  if (!is.null(n_neighbors) && n_neighbors < neighbors$params$n_neighbors) {
    conn <- select_top_neighbors(conn, n_neighbors, mode)
  }
  
  # Recurse neighbors if requested
  if (recurse_neighbors) {
    conn2 <- conn %*% conn
    conn <- conn + 0.5 * conn2
    conn@x <- pmin(conn@x, 1)  # Clip to max 1
  }
  
  # Ensure diagonal is 1
  diag(conn) <- 1
  
  # Row-normalize
  row_sums <- Matrix::rowSums(conn)
  row_sums[row_sums == 0] <- 1
  conn <- conn / row_sums
  
  return(as(conn, "CsparseMatrix"))
}

#' Get Distances Matrix
#'
#' @description Extract the distance matrix from a Seurat object.
#'
#' @param seurat_obj A Seurat object with computed neighbors
#' @param n_neighbors Optional: select only top n neighbors
#'
#' @return Sparse distance matrix
#'
#' @export
get_distances <- function(seurat_obj, n_neighbors = NULL) {
  get_connectivities(seurat_obj, mode = "distances", n_neighbors = n_neighbors)
}

#' Get Number of Neighbors
#'
#' @description Get the number of neighbors used in computation.
#'
#' @param seurat_obj A Seurat object with computed neighbors
#' @return Integer number of neighbors
#' @keywords internal
get_n_neighbors <- function(seurat_obj) {
  if (is.null(seurat_obj@misc$neighbors)) {
    return(0L)
  }
  seurat_obj@misc$neighbors$params$n_neighbors
}

#' Select Top Neighbors
#'
#' @description Reduce connectivity matrix to top n neighbors per cell.
#'
#' @param conn Sparse connectivity matrix
#' @param n_neighbors Number of neighbors to keep
#' @param mode "connectivities" or "distances"
#' @return Modified sparse matrix
#' @keywords internal
select_top_neighbors <- function(conn, n_neighbors, mode = "connectivities") {
  n_cells <- nrow(conn)
  
  # For each row, keep only top n_neighbors
  conn_new <- conn
  
  for (i in seq_len(n_cells)) {
    row_vals <- conn[i, ]
    
    if (length(row_vals@x) > n_neighbors) {
      if (mode == "connectivities") {
        # Keep highest values
        threshold <- sort(row_vals@x, decreasing = TRUE)[n_neighbors]
        row_vals@x[row_vals@x < threshold] <- 0
      } else {
        # Keep smallest non-zero values
        nonzero <- row_vals@x[row_vals@x > 0]
        if (length(nonzero) > n_neighbors) {
          threshold <- sort(nonzero)[n_neighbors]
          row_vals@x[row_vals@x > threshold] <- 0
        }
      }
      conn_new[i, ] <- row_vals
    }
  }
  
  conn_new <- Matrix::drop0(conn_new)
  return(conn_new)
}

#' Get Neighbor Indices
#'
#' @description Get neighbor indices for all or specific cells.
#'
#' @param seurat_obj A Seurat object with computed neighbors
#' @param cells Cell indices or names (NULL for all)
#' @return Matrix of neighbor indices (cells x n_neighbors)
#' @keywords internal
get_neighbor_indices <- function(seurat_obj, cells = NULL) {
  if (is.null(seurat_obj@misc$neighbors)) {
    stop("Neighbors not computed", call. = FALSE)
  }
  
  indices <- seurat_obj@misc$neighbors$indices
  
  if (!is.null(cells)) {
    if (is.character(cells)) {
      cells <- match(cells, colnames(seurat_obj))
    }
    indices <- indices[cells, , drop = FALSE]
  }
  
  return(indices)
}

#' Verify Neighbor Graph
#'
#' @description Check that neighbor graph is valid and complete.
#'
#' @param seurat_obj A Seurat object
#' @return TRUE if valid, throws error otherwise
#' @keywords internal
verify_neighbors <- function(seurat_obj) {
  if (is.null(seurat_obj@misc$neighbors)) {
    stop("Neighbor graph not found. Run compute_neighbors first.", call. = FALSE)
  }

  
  neighbors <- seurat_obj@misc$neighbors
  
  # Check that all cells have neighbors
  n_neighs <- Matrix::rowSums(neighbors$distances_sparse > 0)
  
  if (min(n_neighs) == 0) {
    warning("Some cells have no neighbors. Consider recomputing neighbors.")
  }
  
  # Check for corrupted graph
  if (min(n_neighs) * 2 < max(n_neighs)) {
    warning("Neighbor graph may be corrupted or biased. Consider recomputing.")
  }
  
  TRUE
}

#' Get Iterative Neighbors
#'
#' @description Get neighbors of neighbors (recursive expansion).
#'
#' @param seurat_obj A Seurat object with computed neighbors
#' @param cell_idx Cell index
#' @param n_recurse Number of recursion levels (default: 2)
#' @param max_neighbors Maximum neighbors to return (default: NULL)
#' @return Vector of neighbor indices
#' @keywords internal
get_iterative_neighbors <- function(seurat_obj, 
                                     cell_idx, 
                                     n_recurse = 2,
                                     max_neighbors = NULL) {
  
  indices <- seurat_obj@misc$neighbors$indices
  
  if (is.null(max_neighbors)) {
    max_neighbors <- ncol(indices) * n_recurse
  }
  
  # Convert to 0-based
  indices_0 <- indices - 1L
  
  # Use C++ function
  neighs <- get_iterative_neighbors_cpp(
    indices_0,
    cell_idx - 1L,  # Convert to 0-based
    n_recurse,
    max_neighbors
  )
  
  # Convert back to 1-based
  return(neighs + 1L)
}

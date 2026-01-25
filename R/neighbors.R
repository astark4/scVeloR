#' @title Neighbor Graph Computation
#' @description Compute k-nearest neighbor graphs for velocity analysis.
#' @name neighbors
NULL

#' Compute Neighbors
#'
#' @description Compute k-nearest neighbors using PCA or expression data.
#'
#' @param object Seurat object
#' @param n_neighbors Number of neighbors
#' @param n_pcs Number of principal components to use
#' @param method Method for neighbor computation: "annoy", "hnsw", or "exact"
#' @param use_rep Dimensionality reduction to use (default: "pca")
#' @param metric Distance metric: "euclidean" or "cosine"
#' @param verbose Print progress
#'
#' @return Seurat object with neighbors in misc$scVeloR$neighbors
#' @export
compute_neighbors <- function(object,
                               n_neighbors = 30,
                               n_pcs = 30,
                               method = "exact",
                               use_rep = "pca",
                               metric = "euclidean",
                               verbose = TRUE) {
  
  if (verbose) {
    message(sprintf("Computing %d neighbors...", n_neighbors))
  }
  
  # Get reduced dimension data
  if (use_rep %in% names(object@reductions)) {
    X <- Seurat::Embeddings(object, reduction = use_rep)
    
    if (ncol(X) > n_pcs) {
      X <- X[, 1:n_pcs, drop = FALSE]
    }
  } else {
    # Use expression data directly
    if (!is.null(object@misc$scVeloR$Ms)) {
      X <- object@misc$scVeloR$Ms
    } else {
      X <- get_layer_matrix(object, "spliced")
    }
  }
  
  n_cells <- nrow(X)
  
  # Compute neighbors
  if (method == "annoy" && requireNamespace("RcppAnnoy", quietly = TRUE)) {
    result <- compute_neighbors_annoy(X, n_neighbors, metric)
  } else if (method == "hnsw" && requireNamespace("RcppHNSW", quietly = TRUE)) {
    result <- compute_neighbors_hnsw(X, n_neighbors, metric)
  } else {
    result <- compute_neighbors_exact(X, n_neighbors, metric)
  }
  
  indices <- result$indices
  distances <- result$distances
  
  # Compute connectivities (UMAP-style)
  connectivities <- compute_connectivities(
    indices = indices,
    distances = distances,
    n_cells = n_cells
  )
  
  if (verbose) {
    message(sprintf("  Computed connectivities: %d x %d sparse matrix", 
                    nrow(connectivities), ncol(connectivities)))
  }
  
  # Store results
  if (is.null(object@misc$scVeloR)) {
    object@misc$scVeloR <- list()
  }
  
  object@misc$scVeloR$neighbors <- list(
    indices = indices,
    distances = distances,
    connectivities = connectivities,
    n_neighbors = n_neighbors,
    n_pcs = n_pcs,
    method = method,
    metric = metric
  )
  
  object
}

#' Compute Exact Neighbors
#'
#' @description Compute exact k-nearest neighbors using distance matrix.
#'
#' @param X Data matrix (cells x features)
#' @param n_neighbors Number of neighbors
#' @param metric Distance metric
#'
#' @return List with indices and distances
#' @keywords internal
compute_neighbors_exact <- function(X, n_neighbors, metric = "euclidean") {
  
  n_cells <- nrow(X)
  
  # Compute distance matrix
  if (metric == "cosine") {
    # Normalize rows
    X_norm <- X / sqrt(rowSums(X^2))
    X_norm[!is.finite(X_norm)] <- 0
    
    # Cosine distance = 1 - cosine similarity
    D <- 1 - tcrossprod(X_norm)
  } else {
    D <- as.matrix(dist(X, method = metric))
  }
  
  # Find k-nearest neighbors for each cell
  indices <- matrix(0L, n_cells, n_neighbors)
  distances <- matrix(0, n_cells, n_neighbors)
  
  for (i in seq_len(n_cells)) {
    # Get distances from cell i to all other cells
    d <- D[i, ]
    d[i] <- Inf  # Exclude self
    
    # Find k nearest
    ord <- order(d)[1:n_neighbors]
    indices[i, ] <- ord
    distances[i, ] <- d[ord]
  }
  
  list(indices = indices, distances = distances)
}

#' Compute Neighbors Using Annoy
#'
#' @description Compute approximate neighbors using Annoy.
#'
#' @param X Data matrix
#' @param n_neighbors Number of neighbors
#' @param metric Distance metric
#' @param n_trees Number of trees
#'
#' @return List with indices and distances
#' @keywords internal
compute_neighbors_annoy <- function(X, n_neighbors, metric = "euclidean", n_trees = 50) {
  
  if (!requireNamespace("RcppAnnoy", quietly = TRUE)) {
    warning("RcppAnnoy not available, using exact method")
    return(compute_neighbors_exact(X, n_neighbors, metric))
  }
  
  n_cells <- nrow(X)
  n_dims <- ncol(X)
  
  # Create Annoy index
  if (metric == "cosine") {
    ann <- new(RcppAnnoy::AnnoyAngular, n_dims)
  } else {
    ann <- new(RcppAnnoy::AnnoyEuclidean, n_dims)
  }
  
  # Add items
  for (i in seq_len(n_cells)) {
    ann$addItem(i - 1, X[i, ])
  }
  
  # Build index
  ann$build(n_trees)
  
  # Query neighbors
  indices <- matrix(0L, n_cells, n_neighbors)
  distances <- matrix(0, n_cells, n_neighbors)
  
  for (i in seq_len(n_cells)) {
    result <- ann$getNNsByItemList(i - 1, n_neighbors + 1, -1, TRUE)
    
    # Remove self (first entry)
    idx <- result$item[-1][1:n_neighbors] + 1
    dst <- result$distance[-1][1:n_neighbors]
    
    indices[i, ] <- idx
    distances[i, ] <- dst
  }
  
  list(indices = indices, distances = distances)
}

#' Compute Neighbors Using HNSW
#'
#' @description Compute approximate neighbors using HNSW algorithm.
#'
#' @param X Data matrix
#' @param n_neighbors Number of neighbors
#' @param metric Distance metric
#'
#' @return List with indices and distances
#' @keywords internal
compute_neighbors_hnsw <- function(X, n_neighbors, metric = "euclidean") {
  
  if (!requireNamespace("RcppHNSW", quietly = TRUE)) {
    warning("RcppHNSW not available, using exact method")
    return(compute_neighbors_exact(X, n_neighbors, metric))
  }
  
  n_cells <- nrow(X)
  
  # Create HNSW index
  if (metric == "cosine") {
    ann <- RcppHNSW::hnsw_build(X, distance = "cosine")
  } else {
    ann <- RcppHNSW::hnsw_build(X, distance = "euclidean")
  }
  
  # Search neighbors
  result <- RcppHNSW::hnsw_search(X, ann, k = n_neighbors + 1)
  
  # Remove self (first column)
  indices <- result$idx[, -1, drop = FALSE]
  distances <- result$dist[, -1, drop = FALSE]
  
  list(indices = indices, distances = distances)
}

#' Compute Connectivities
#'
#' @description Compute UMAP-style connectivity weights from neighbor graph.
#'
#' @param indices Neighbor indices
#' @param distances Neighbor distances
#' @param n_cells Number of cells
#' @param local_connectivity Local connectivity parameter
#'
#' @return Sparse connectivity matrix
#' @keywords internal
compute_connectivities <- function(indices,
                                    distances,
                                    n_cells,
                                    local_connectivity = 1) {
  
  n_neighbors <- ncol(indices)
  
  # Compute sigma (bandwidth) for each cell
  sigma <- numeric(n_cells)
  rho <- numeric(n_cells)  # Distance to nearest neighbor
  
  for (i in seq_len(n_cells)) {
    d <- distances[i, ]
    d_sorted <- sort(d[d > 0])
    
    if (length(d_sorted) > 0) {
      rho[i] <- d_sorted[min(local_connectivity, length(d_sorted))]
      
      # Binary search for sigma
      target <- log2(n_neighbors)
      lo <- 1e-10
      hi <- max(d) * 10
      
      for (iter in 1:64) {
        mid <- (lo + hi) / 2
        
        # Compute sum of weights
        val <- sum(exp(-pmax(d - rho[i], 0) / mid))
        
        if (abs(val - target) < 1e-5) break
        
        if (val > target) {
          hi <- mid
        } else {
          lo <- mid
        }
      }
      
      sigma[i] <- mid
    } else {
      sigma[i] <- 1
      rho[i] <- 0
    }
  }
  
  # Build connectivity matrix
  i_idx <- integer(0)
  j_idx <- integer(0)
  x_val <- numeric(0)
  
  for (i in seq_len(n_cells)) {
    neighbors <- indices[i, ]
    d <- distances[i, ]
    
    for (k in seq_along(neighbors)) {
      j <- neighbors[k]
      
      if (j > 0 && j <= n_cells && j != i) {
        # Compute weight
        w <- exp(-max(d[k] - rho[i], 0) / sigma[i])
        
        if (w > 0) {
          i_idx <- c(i_idx, i)
          j_idx <- c(j_idx, j)
          x_val <- c(x_val, w)
        }
      }
    }
  }
  
  # Create sparse matrix
  conn <- sparseMatrix(
    i = i_idx,
    j = j_idx,
    x = x_val,
    dims = c(n_cells, n_cells)
  )
  
  # Symmetrize: C = (C + C') / 2
  conn <- (conn + t(conn)) / 2
  
  # Ensure it's a sparse Matrix
  conn <- drop0(conn)
  
  conn
}

#' Get Neighbor Connectivities
#'
#' @description Extract connectivity matrix from Seurat object.
#'
#' @param object Seurat object
#'
#' @return Sparse connectivity matrix
#' @export
get_connectivities <- function(object) {
  
  if (!is.null(object@misc$scVeloR$neighbors$connectivities)) {
    return(object@misc$scVeloR$neighbors$connectivities)
  }
  
  # Try Seurat's graph
  graphs <- Seurat::Graphs(object)
  if (length(graphs) > 0) {
    graph_name <- graphs[grep("snn", graphs, ignore.case = TRUE)[1]]
    if (!is.na(graph_name)) {
      return(object@graphs[[graph_name]])
    }
  }
  
  NULL
}

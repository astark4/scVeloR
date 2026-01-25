# Helper functions for testing scVeloR

#' Create Simple Test Seurat-like Object
#' 
#' Creates a minimal list structure for testing without Seurat dependency
#' @param n_cells Number of cells
#' @param n_genes Number of genes
#' @param seed Random seed
create_test_object <- function(n_cells = 50, n_genes = 30, seed = 42) {
  set.seed(seed)
  
  # Create correlated spliced/unspliced counts
  spliced <- matrix(rpois(n_cells * n_genes, lambda = 20), 
                   nrow = n_cells, ncol = n_genes)
  
  unspliced <- matrix(0, n_cells, n_genes)
  for (i in seq_len(n_genes)) {
    gamma <- runif(1, 0.3, 2)
    unspliced[, i] <- round(spliced[, i] * gamma + rnorm(n_cells, 0, 2))
  }
  unspliced[unspliced < 0] <- 0
  
  colnames(spliced) <- paste0("Gene", seq_len(n_genes))
  colnames(unspliced) <- paste0("Gene", seq_len(n_genes))
  rownames(spliced) <- paste0("Cell", seq_len(n_cells))
  rownames(unspliced) <- paste0("Cell", seq_len(n_cells))
  
  # Create embedding
  embedding <- matrix(rnorm(n_cells * 2), n_cells, 2)
  colnames(embedding) <- c("UMAP_1", "UMAP_2")
  rownames(embedding) <- rownames(spliced)
  
  list(
    spliced = spliced,
    unspliced = unspliced,
    embedding = embedding,
    n_cells = n_cells,
    n_genes = n_genes
  )
}

#' Create Test Velocity Data
#' 
#' Creates velocity-like data for testing graph construction
#' @param n_cells Number of cells
#' @param n_genes Number of genes  
#' @param seed Random seed
create_test_velocity <- function(n_cells = 50, n_genes = 30, seed = 42) {
  set.seed(seed)
  
  # Create velocity matrix
  velocity <- matrix(rnorm(n_cells * n_genes), n_cells, n_genes)
  
  # Create expression matrix
  expression <- matrix(abs(rnorm(n_cells * n_genes)) + 1, n_cells, n_genes)
  
  colnames(velocity) <- paste0("Gene", seq_len(n_genes))
  colnames(expression) <- paste0("Gene", seq_len(n_genes))
  rownames(velocity) <- paste0("Cell", seq_len(n_cells))
  rownames(expression) <- paste0("Cell", seq_len(n_cells))
  
  list(
    velocity = velocity,
    expression = expression
  )
}

#' Check if result is approximately equal
#' 
#' @param x First value
#' @param y Second value
#' @param tol Tolerance
approx_equal <- function(x, y, tol = 1e-6) {
  all(abs(x - y) < tol, na.rm = TRUE)
}

#' Simulate ODE trajectory
#' 
#' Simulates a trajectory using the analytical solutions
#' @param t_max Maximum time
#' @param n_points Number of time points
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
simulate_trajectory <- function(t_max = 10, n_points = 100,
                                 alpha = 2, beta = 0.5, gamma = 0.3) {
  
  # Time points
  t <- seq(0, t_max, length.out = n_points)
  
  # Induction phase
  u_on <- unspliced(t, u0 = 0, alpha = alpha, beta = beta)
  s_on <- spliced(t, s0 = 0, u0 = 0, alpha = alpha, beta = beta, gamma = gamma)
  
  list(
    t = t,
    u = u_on,
    s = s_on,
    alpha = alpha,
    beta = beta,
    gamma = gamma
  )
}

#' Create simple neighbor indices
#' 
#' @param n_cells Number of cells
#' @param n_neighbors Number of neighbors
#' @param seed Random seed
create_neighbor_indices <- function(n_cells, n_neighbors, seed = 42) {
  set.seed(seed)
  
  indices <- matrix(0L, n_cells, n_neighbors)
  for (i in seq_len(n_cells)) {
    candidates <- setdiff(seq_len(n_cells), i)
    indices[i, ] <- sample(candidates, min(n_neighbors, length(candidates)))
  }
  
  indices
}

#' Check matrix properties
#' 
#' @param mat Matrix to check
check_matrix_valid <- function(mat) {
  list(
    is_finite = all(is.finite(mat)),
    no_na = !any(is.na(mat)),
    dims = dim(mat)
  )
}

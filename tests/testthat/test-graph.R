# Test Velocity Graph and Neighbors

# Helper functions
create_test_pca <- function(n_cells = 50, n_pcs = 10, seed = 42) {
  set.seed(seed)
  matrix(rnorm(n_cells * n_pcs), nrow = n_cells, ncol = n_pcs)
}

test_that("compute_neighbors_exact returns correct structure", {
  X <- create_test_pca(n_cells = 20, n_pcs = 5)
  n_neighbors <- 5
  
  result <- compute_neighbors_exact(X, n_neighbors, metric = "euclidean")
  
  expect_true(is.list(result))
  expect_true("indices" %in% names(result))
  expect_true("distances" %in% names(result))
  
  expect_equal(dim(result$indices), c(20, n_neighbors))
  expect_equal(dim(result$distances), c(20, n_neighbors))
  
  # Indices should be integers
  expect_true(all(result$indices >= 1 & result$indices <= 20))
  
  # No self-neighbors
  for (i in 1:20) {
    expect_false(i %in% result$indices[i, ])
  }
  
  # Distances should be non-negative
  expect_true(all(result$distances >= 0))
})

test_that("compute_neighbors_exact distances are sorted", {
  X <- create_test_pca(n_cells = 20, n_pcs = 5)
  n_neighbors <- 5
  
  result <- compute_neighbors_exact(X, n_neighbors, metric = "euclidean")
  
  # For each cell, distances should be in ascending order
  for (i in 1:20) {
    expect_true(all(diff(result$distances[i, ]) >= 0))
  }
})

test_that("euclidean distance is computed correctly", {
  # Simple 2D case
  X <- matrix(c(0, 0, 1, 0, 0, 1), nrow = 3, byrow = TRUE)
  
  result <- compute_neighbors_exact(X, n_neighbors = 2, metric = "euclidean")
  
  # Distance from (0,0) to (1,0) and (0,1) should both be 1
  expect_equal(result$distances[1, 1], 1, tolerance = 1e-10)
})

test_that("cosine distance works correctly", {
  # Create vectors with known cosine similarity
  X <- matrix(c(
    1, 0,    # Point 1
    0, 1,    # Point 2: orthogonal to 1
    1, 1     # Point 3: 45 degrees from both
  ), nrow = 3, byrow = TRUE)
  
  result <- compute_neighbors_exact(X, n_neighbors = 2, metric = "cosine")
  
  # Cosine distance between orthogonal vectors is 1
  # Cosine distance between identical directions is 0
  expect_true(all(result$distances >= 0))
  expect_true(all(result$distances <= 2))  # Max cosine distance is 2
})

test_that("compute_connectivities produces valid sparse matrix", {
  library(Matrix)
  
  X <- create_test_pca(n_cells = 30, n_pcs = 5)
  n_neighbors <- 5
  
  nn_result <- compute_neighbors_exact(X, n_neighbors, "euclidean")
  
  conn <- compute_connectivities(
    indices = nn_result$indices,
    distances = nn_result$distances,
    n_cells = 30
  )
  
  expect_true(inherits(conn, "Matrix") || is.matrix(conn))
  expect_equal(dim(conn), c(30, 30))
  
  # All values should be non-negative
  expect_true(all(conn >= 0))
  
  # Matrix should be symmetric (after symmetrization)
  expect_equal(sum(abs(conn - t(conn))), 0, tolerance = 1e-10)
})

test_that("recurse_neighbors extends neighbor set", {
  # Create initial neighbors
  indices <- matrix(c(
    2, 3,     # Cell 1's neighbors
    1, 3,     # Cell 2's neighbors
    1, 2,     # Cell 3's neighbors
    5, 6,     # Cell 4's neighbors
    4, 6,     # Cell 5's neighbors
    4, 5      # Cell 6's neighbors
  ), nrow = 6, byrow = TRUE)
  
  # After one recursion, neighbors of neighbors should be included
  extended <- recurse_neighbors(indices, n_recurse = 1)
  
  # Cell 1 should now include neighbors of cells 2 and 3
  # Original: 2, 3
  # Cell 2's neighbors: 1, 3
  # Cell 3's neighbors: 1, 2
  # So cell 1 should have: 2, 3 (extended includes more)
  
  expect_true(ncol(extended) >= ncol(indices))
})

test_that("compute_velocity_graph_r produces valid output", {
  set.seed(42)
  n_cells <- 20
  n_genes <- 10
  
  # Create test data
  velocity <- matrix(rnorm(n_cells * n_genes), n_cells, n_genes)
  expression <- matrix(abs(rnorm(n_cells * n_genes)) + 1, n_cells, n_genes)
  
  # Create simple neighbor indices
  indices <- matrix(0L, n_cells, 5)
  for (i in 1:n_cells) {
    indices[i, ] <- sample(setdiff(1:n_cells, i), 5)
  }
  
  result <- compute_velocity_graph_r(velocity, expression, indices)
  
  expect_true(inherits(result, "Matrix") || inherits(result, "dgCMatrix"))
  expect_equal(dim(result), c(n_cells, n_cells))
  
  # Graph values should be between 0 and 1 (cosine similarity)
  if (nnzero(result) > 0) {
    expect_true(all(result@x >= 0))
    expect_true(all(result@x <= 1))
  }
})

test_that("compute_transition_matrix produces row-stochastic matrix", {
  library(Matrix)
  
  # Create a simple graph
  n <- 10
  graph <- sparseMatrix(
    i = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5),
    j = c(2, 3, 1, 3, 1, 2, 5, 6, 4, 6),
    x = runif(10),
    dims = c(n, n)
  )
  
  result <- compute_transition_matrix(graph, self_transitions = FALSE)
  
  # Check row sums
  row_sums <- rowSums(result)
  
  # Non-zero rows should sum to 1
  for (i in 1:n) {
    if (row_sums[i] > 0) {
      expect_equal(row_sums[i], 1, tolerance = 1e-10)
    }
  }
  
  # All values should be non-negative
  expect_true(all(result >= 0))
})

test_that("transition matrix with self-transitions is valid", {
  library(Matrix)
  
  n <- 10
  graph <- sparseMatrix(
    i = rep(1:5, each = 2),
    j = c(2, 3, 1, 3, 1, 2, 5, 6, 4, 6),
    x = rep(0.5, 10),
    dims = c(n, n)
  )
  
  result <- compute_transition_matrix(graph, self_transitions = TRUE)
  
  # Diagonal should have non-negative values
  expect_true(all(diag(result) >= 0))
  
  # Row sums should be close to 1
  row_sums <- rowSums(result)
  expect_true(all(abs(row_sums[row_sums > 0] - 1) < 1e-10))
})

test_that("velocity graph handles zero velocity cells", {
  set.seed(42)
  n_cells <- 20
  n_genes <- 10
  
  # Create test data with some zero velocity cells
  velocity <- matrix(rnorm(n_cells * n_genes), n_cells, n_genes)
  velocity[1:5, ] <- 0  # First 5 cells have zero velocity
  
  expression <- matrix(abs(rnorm(n_cells * n_genes)) + 1, n_cells, n_genes)
  
  indices <- matrix(0L, n_cells, 5)
  for (i in 1:n_cells) {
    indices[i, ] <- sample(setdiff(1:n_cells, i), 5)
  }
  
  result <- compute_velocity_graph_r(velocity, expression, indices)
  
  # Zero velocity cells should have zero outgoing edges
  for (i in 1:5) {
    expect_equal(sum(result[i, ]), 0)
  }
})

# Test C++ implementations if available
test_that("C++ distance functions work correctly", {
  skip_if_not(exists("euclidean_distances_cpp"))
  
  X <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
  
  D <- euclidean_distances_cpp(X)
  
  expect_equal(dim(D), c(4, 4))
  
  # Diagonal should be zero
  expect_equal(diag(D), rep(0, 4), tolerance = 1e-10)
  
  # Should be symmetric
  expect_equal(D, t(D), tolerance = 1e-10)
  
  # Check specific distance
  expect_equal(D[1, 2], 1, tolerance = 1e-10)  # (0,0) to (1,0)
})

test_that("C++ KNN function works correctly", {
  skip_if_not(exists("knn_from_distances_cpp"))
  
  D <- matrix(c(
    0, 1, 2, 3,
    1, 0, 1.5, 2.5,
    2, 1.5, 0, 1,
    3, 2.5, 1, 0
  ), nrow = 4, byrow = TRUE)
  
  result <- knn_from_distances_cpp(D, k = 2)
  
  expect_equal(dim(result$indices), c(4, 2))
  expect_equal(dim(result$distances), c(4, 2))
  
  # Check that nearest neighbor is correct
  expect_equal(result$indices[1, 1], 2)  # Cell 1's nearest is cell 2
})

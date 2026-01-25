# Test Preprocessing Functions

# Helper function to create test matrices
create_test_matrices <- function(n_cells = 50, n_genes = 30, seed = 42) {
  set.seed(seed)
  
  # Create correlated spliced/unspliced counts
  spliced <- matrix(rpois(n_cells * n_genes, lambda = 20), 
                   nrow = n_cells, ncol = n_genes)
  
  # Unspliced is correlated with spliced (gamma relationship)
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
  
  list(spliced = spliced, unspliced = unspliced)
}

test_that("make_dense converts sparse matrix correctly", {
  # Create sparse matrix
  library(Matrix)
  sparse_mat <- sparseMatrix(
    i = c(1, 2, 3, 1, 2),
    j = c(1, 1, 1, 2, 2),
    x = c(1, 2, 3, 4, 5),
    dims = c(3, 2)
  )
  
  dense_mat <- make_dense(sparse_mat)
  
  expect_true(is.matrix(dense_mat))
  expect_false(inherits(dense_mat, "sparseMatrix"))
  expect_equal(dim(dense_mat), c(3, 2))
  expect_equal(dense_mat[1, 1], 1)
  expect_equal(dense_mat[2, 1], 2)
})

test_that("clip function works correctly", {
  x <- c(-5, 0, 5, 10, 15)
  
  # Clip to [0, 10]
  result <- clip(x, lower = 0, upper = 10)
  expect_equal(result, c(0, 0, 5, 10, 10))
  
  # Only lower bound
  result2 <- clip(x, lower = 0)
  expect_equal(result2, c(0, 0, 5, 10, 15))
  
  # Only upper bound
  result3 <- clip(x, upper = 10)
  expect_equal(result3, c(-5, 0, 5, 10, 10))
})

test_that("scale_minmax scales to [0,1]", {
  x <- c(10, 20, 30, 40, 50)
  
  result <- scale_minmax(x)
  
  expect_equal(min(result), 0)
  expect_equal(max(result), 1)
  expect_equal(result, c(0, 0.25, 0.5, 0.75, 1))
})

test_that("scale_minmax handles constant vector", {
  x <- c(5, 5, 5, 5)
  
  result <- scale_minmax(x)
  
  # Should return 0.5 for all values
  expect_true(all(result == 0.5))
})

test_that("rowVars computes row-wise variance correctly", {
  mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, byrow = TRUE)
  
  result <- rowVars(mat)
  
  # Manual calculation: var of (1,2,3), (4,5,6), (7,8,9)
  expected <- c(var(c(1,2,3)) * 2/3, var(c(4,5,6)) * 2/3, var(c(7,8,9)) * 2/3)  # population variance
  
  expect_length(result, 3)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("colVars computes column-wise variance correctly", {
  mat <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, byrow = TRUE)
  
  result <- colVars(mat)
  
  expect_length(result, 3)
})

test_that("weighted_mean computes correctly", {
  x <- c(1, 2, 3, 4, 5)
  w <- c(1, 1, 1, 1, 1)  # Uniform weights
  
  result <- weighted_mean(x, w)
  expect_equal(result, mean(x))
  
  # Non-uniform weights
  w2 <- c(0, 0, 0, 0, 1)
  result2 <- weighted_mean(x, w2)
  expect_equal(result2, 5)
})

test_that("vcorrcoef computes correlation correctly", {
  set.seed(42)
  x <- rnorm(100)
  Y <- matrix(rnorm(300), nrow = 100, ncol = 3)
  Y[, 1] <- x + rnorm(100, sd = 0.1)  # High correlation with x
  Y[, 2] <- -x + rnorm(100, sd = 0.1)  # Negative correlation
  Y[, 3] <- rnorm(100)  # No correlation
  
  result <- vcorrcoef(x, Y)
  
  expect_length(result, 3)
  expect_true(result[1] > 0.9)  # High positive
  expect_true(result[2] < -0.9)  # High negative
  expect_true(abs(result[3]) < 0.3)  # Near zero
})

test_that("safe_div handles division by zero", {
  result <- safe_div(c(1, 2, 3), c(1, 0, 3))
  
  expect_equal(result[1], 1)
  expect_equal(result[2], 0)  # Default value for 0/0
  expect_equal(result[3], 1)
})

test_that("sparse matrix functions work correctly", {
  library(Matrix)
  
  # Create test sparse matrix
  sparse_mat <- sparseMatrix(
    i = c(1, 1, 2, 2, 3),
    j = c(1, 2, 1, 2, 1),
    x = c(1, 2, 3, 4, 5),
    dims = c(3, 3)
  )
  
  # Test row sums
  rs <- sparse_row_sums(sparse_mat)
  expect_equal(rs, c(3, 7, 5))
  
  # Test column sums
  cs <- sparse_col_sums(sparse_mat)
  expect_equal(cs, c(9, 6, 0))
})

test_that("log1p_offset transforms correctly", {
  x <- c(0, 1, 10, 100)
  
  result <- log1p_offset(x, offset = 1)
  expected <- log(x + 1)
  
  expect_equal(result, expected)
})

test_that("object_size_str returns readable format", {
  small_obj <- 1:10
  large_obj <- matrix(rnorm(1e6), nrow = 1000)
  
  small_str <- object_size_str(small_obj)
  large_str <- object_size_str(large_obj)
  
  expect_true(grepl("B|KB|MB|GB", small_str))
  expect_true(grepl("B|KB|MB|GB", large_str))
})

# Test Matrix operations that are used in preprocessing
test_that("Matrix package operations work as expected", {
  library(Matrix)
  
  # Create connectivity matrix
  n <- 10
  k <- 3
  
  # Simple k-nearest neighbor matrix
  conn <- sparseMatrix(
    i = rep(1:n, each = k),
    j = as.vector(sapply(1:n, function(i) sample(setdiff(1:n, i), k))),
    x = rep(1, n * k),
    dims = c(n, n)
  )
  
  # Row normalize
  row_sums <- Matrix::rowSums(conn)
  conn_norm <- conn / row_sums
  
  # Test that rows sum to 1
  expect_equal(as.vector(Matrix::rowSums(conn_norm)), rep(1, n), tolerance = 1e-10)
  
  # Test matrix multiplication
  X <- matrix(rnorm(n * 5), n, 5)
  result <- as.matrix(conn_norm %*% X)
  
  expect_equal(dim(result), c(n, 5))
})

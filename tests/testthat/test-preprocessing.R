# Test preprocessing functions

test_that("normalize_per_cell works correctly", {
  # Create test sparse matrix
  mat <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2, 3, 3),
    j = c(1, 2, 1, 2, 1, 2),
    x = c(10, 20, 30, 40, 50, 60)
  )
  colnames(mat) <- c("cell1", "cell2")
  rownames(mat) <- c("gene1", "gene2", "gene3")
  
  # Test normalization
  result <- scVeloR:::normalize_per_cell(mat, target_sum = 100, log_transform = FALSE)
  
  # Check column sums
  col_sums <- Matrix::colSums(result)
  expect_true(all(abs(col_sums - 100) < 1e-6))
})

test_that("make_dense works", {
  # Sparse matrix
  sparse_mat <- Matrix::sparseMatrix(
    i = c(1, 2),
    j = c(1, 2),
    x = c(1, 2)
  )
  
  result <- scVeloR:::make_dense(sparse_mat)
  expect_true(is.matrix(result))
  expect_false(inherits(result, "sparseMatrix"))
  
  # Dense matrix should pass through
  dense_mat <- matrix(c(1, 2, 3, 4), nrow = 2)
  result2 <- scVeloR:::make_dense(dense_mat)
  expect_equal(result2, dense_mat)
})

test_that("make_sparse works", {
  # Dense matrix
  dense_mat <- matrix(c(1, 0, 0, 2), nrow = 2)
  
  result <- scVeloR:::make_sparse(dense_mat)
  expect_true(inherits(result, "sparseMatrix"))
})

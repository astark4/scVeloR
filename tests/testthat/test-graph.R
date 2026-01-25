# Test graph functions

test_that("normalize_matrix works for row normalization", {
  mat <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2),
    j = c(1, 2, 1, 2),
    x = c(1, 2, 3, 4)
  )
  
  result <- scVeloR:::normalize_matrix(mat, axis = 1)
  
  # Row sums should be 1
  row_sums <- Matrix::rowSums(result)
  expect_true(all(abs(row_sums - 1) < 1e-6))
})

test_that("normalize_matrix works for column normalization", {
  mat <- Matrix::sparseMatrix(
    i = c(1, 1, 2, 2),
    j = c(1, 2, 1, 2),
    x = c(1, 2, 3, 4)
  )
  
  result <- scVeloR:::normalize_matrix(mat, axis = 2)
  
  # Column sums should be 1
  col_sums <- Matrix::colSums(result)
  expect_true(all(abs(col_sums - 1) < 1e-6))
})

test_that("prod_sum works", {
  a <- matrix(c(1, 2, 3, 4), nrow = 2)
  b <- matrix(c(2, 3, 4, 5), nrow = 2)
  
  # Row sums of element-wise product
  result <- scVeloR:::prod_sum(a, b, axis = 1)
  expected <- c(1*2 + 3*4, 2*3 + 4*5)
  expect_equal(result, expected)
  
  # Column sums
  result2 <- scVeloR:::prod_sum(a, b, axis = 2)
  expected2 <- c(1*2 + 2*3, 3*4 + 4*5)
  expect_equal(result2, expected2)
})

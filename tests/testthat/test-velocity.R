# Test velocity functions

test_that("cosine_similarity works", {
  a <- c(1, 0, 0)
  b <- c(1, 0, 0)
  expect_equal(scVeloR:::cosine_similarity(a, b), 1)
  
  a <- c(1, 0, 0)
  b <- c(0, 1, 0)
  expect_equal(scVeloR:::cosine_similarity(a, b), 0)
  
  a <- c(1, 0, 0)
  b <- c(-1, 0, 0)
  expect_equal(scVeloR:::cosine_similarity(a, b), -1)
})

test_that("scale_01 works", {
  x <- c(0, 5, 10)
  result <- scVeloR:::scale_01(x)
  expect_equal(result, c(0, 0.5, 1))
  
  # Constant vector
  x2 <- c(5, 5, 5)
  result2 <- scVeloR:::scale_01(x2)
  expect_true(all(result2 == 0.5))
})

test_that("R_squared calculation works", {
  residual <- matrix(c(0.1, 0.2, 0.1, 0.2), nrow = 2)
  total <- matrix(c(1, 2, 1, 2), nrow = 2)
  
  r2 <- scVeloR:::R_squared(residual, total)
  expect_true(all(r2 >= 0 & r2 <= 1))
})

test_that("leastsq_generalized returns valid results", {
  set.seed(42)
  s <- runif(100)
  u <- 0.5 * s + rnorm(100, sd = 0.1)
  ss <- s^2 + rnorm(100, sd = 0.01)
  us <- s * u + rnorm(100, sd = 0.01)
  
  result <- scVeloR:::leastsq_generalized(s, u, ss, us, fit_offset = FALSE)
  
  expect_true(!is.na(result$gamma))
  expect_true(result$gamma > 0)
  expect_true(result$r2 >= 0 && result$r2 <= 1)
})

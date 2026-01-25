context("Velocity Computation")

# Helper to create test data
create_test_data <- function(n_cells = 100, n_genes = 50) {
  # Simulate spliced and unspliced counts
  set.seed(42)
  
  # Generate correlated spliced/unspliced counts
  spliced <- matrix(rpois(n_cells * n_genes, lambda = 10), 
                   nrow = n_cells, ncol = n_genes)
  unspliced <- matrix(rpois(n_cells * n_genes, lambda = 5), 
                     nrow = n_cells, ncol = n_genes)
  
  # Add correlation
  for (i in seq_len(n_genes)) {
    gamma <- runif(1, 0.5, 2)
    unspliced[, i] <- round(spliced[, i] * gamma + rnorm(n_cells, 0, 2))
    unspliced[unspliced < 0] <- 0
  }
  
  colnames(spliced) <- paste0("Gene", seq_len(n_genes))
  colnames(unspliced) <- paste0("Gene", seq_len(n_genes))
  rownames(spliced) <- paste0("Cell", seq_len(n_cells))
  rownames(unspliced) <- paste0("Cell", seq_len(n_cells))
  
  list(spliced = spliced, unspliced = unspliced)
}

test_that("dynamics functions work correctly", {
  # Test unspliced dynamics
  tau <- seq(0, 5, length.out = 100)
  u <- unspliced(tau, u0 = 0, alpha = 1, beta = 0.5)
  
  expect_length(u, 100)
  expect_equal(u[1], 0, tolerance = 1e-6)
  expect_true(all(diff(u) >= 0))  # Should be monotonically increasing
  
  # Test spliced dynamics
  s <- spliced(tau, s0 = 0, u0 = 0, alpha = 1, beta = 0.5, gamma = 0.3)
  
  expect_length(s, 100)
  expect_equal(s[1], 0, tolerance = 1e-6)
})

test_that("tau_inv functions work correctly", {
  # Test tau_inv_u
  tau_orig <- 2.0
  u <- unspliced(tau_orig, u0 = 0, alpha = 1, beta = 0.5)
  tau_est <- tau_inv_u(u, u0 = 0, alpha = 1, beta = 0.5)
  
  expect_equal(tau_est, tau_orig, tolerance = 0.01)
})

test_that("steady_state function works", {
  ss <- steady_state(alpha = 2, beta = 1, gamma = 0.5)
  
  expect_equal(ss$u_inf, 2)  # alpha/beta
expect_equal(ss$s_inf, 4)  # alpha/gamma
})

test_that("vectorize_params works correctly", {
  t <- c(0.5, 1.5, 2.5, 3.5)
  t_ <- 2.0
  
  params <- vectorize_params(t, t_, alpha = 1, beta = 0.5, gamma = 0.3)
  
  expect_length(params$tau, 4)
  expect_length(params$alpha, 4)
  
  # First two cells are in "on" phase
  expect_equal(params$alpha[1], 1)
  expect_equal(params$alpha[2], 1)
  
  # Last two cells are in "off" phase
  expect_equal(params$alpha[3], 0)
  expect_equal(params$alpha[4], 0)
})

test_that("mRNA function returns both components", {
  tau <- seq(0, 5, length.out = 50)
  result <- mRNA(tau, u0 = 0, s0 = 0, alpha = 1, beta = 0.5, gamma = 0.3)
  
  expect_true("u" %in% names(result))
  expect_true("s" %in% names(result))
  expect_length(result$u, 50)
  expect_length(result$s, 50)
})

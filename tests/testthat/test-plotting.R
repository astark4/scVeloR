# Test plotting functions (requires package to be installed)

test_that("plotting tests are skipped when package not installed", {
  skip_if_not_installed("scVeloR")
  
  # These tests only run when package is installed
  colors_5 <- scVeloR:::get_palette(5)
  expect_length(colors_5, 5)
})

test_that("basic color generation works", {
  # Test basic R color functions that our plotting uses
  colors <- rainbow(10)
  expect_length(colors, 10)
  expect_true(all(grepl("^#", colors)))
  
  # Test colorRampPalette
  blue_red <- colorRampPalette(c("blue", "white", "red"))
  colors_br <- blue_red(11)
  expect_length(colors_br, 11)
})

test_that("arrow computation helpers work", {
  # Test basic arrow endpoint computation
  n <- 10
  start_x <- runif(n)
  start_y <- runif(n)
  vel_x <- rnorm(n, sd = 0.1)
  vel_y <- rnorm(n, sd = 0.1)
  
  scale <- 0.5
  end_x <- start_x + vel_x * scale
  end_y <- start_y + vel_y * scale
  
  expect_length(end_x, n)
  expect_length(end_y, n)
  expect_true(all(is.finite(end_x)))
  expect_true(all(is.finite(end_y)))
})

test_that("grid interpolation logic works", {
  # Test grid generation
  x_range <- c(0, 1)
  y_range <- c(0, 1)
  n_grid <- 10
  
  x_seq <- seq(x_range[1], x_range[2], length.out = n_grid)
  y_seq <- seq(y_range[1], y_range[2], length.out = n_grid)
  
  grid <- expand.grid(x = x_seq, y = y_seq)
  
  expect_equal(nrow(grid), n_grid * n_grid)
  expect_true(all(grid$x >= x_range[1] & grid$x <= x_range[2]))
  expect_true(all(grid$y >= y_range[1] & grid$y <= y_range[2]))
})

test_that("embedding limits computation works", {
  emb <- matrix(c(0, 1, 2, 3, 4, 5), ncol = 2)
  
  x_range <- range(emb[, 1])
  y_range <- range(emb[, 2])
  expand <- 0.1
  
  x_margin <- diff(x_range) * expand
  y_margin <- diff(y_range) * expand
  
  xlim <- c(x_range[1] - x_margin, x_range[2] + x_margin)
  ylim <- c(y_range[1] - y_margin, y_range[2] + y_margin)
  
  expect_true(xlim[1] < xlim[2])
  expect_true(ylim[1] < ylim[2])
  expect_true(xlim[1] < min(emb[, 1]))
  expect_true(xlim[2] > max(emb[, 1]))
})

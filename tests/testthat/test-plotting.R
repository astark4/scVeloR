# Test plotting functions

test_that("get_palette returns correct number of colors", {
  colors_5 <- scVeloR:::get_palette(5)
  expect_length(colors_5, 5)
  
  colors_30 <- scVeloR:::get_palette(30)
  expect_length(colors_30, 30)
})

test_that("velocity_colors returns diverging palette", {
  colors <- scVeloR:::velocity_colors(11)
  expect_length(colors, 11)
  
  # Check it's a character vector of hex colors
  expect_true(all(grepl("^#", colors)))
})

test_that("create_arrow_grid works with simple data", {
  emb <- matrix(runif(200), ncol = 2)
  rownames(emb) <- paste0("cell", 1:100)
  
  V_emb <- matrix(rnorm(200, sd = 0.1), ncol = 2)
  rownames(V_emb) <- paste0("cell", 1:100)
  
  result <- scVeloR:::create_arrow_grid(emb, V_emb, n_grid = 10, min_mass = 1)
  
  expect_true(is.data.frame(result))
  if (nrow(result) > 0) {
    expect_true(all(c("x", "y", "vx", "vy", "mass") %in% names(result)))
  }
})

test_that("theme_velocity returns valid theme", {
  theme <- scVeloR::theme_velocity()
  expect_s3_class(theme, "theme")
})

test_that("get_embedding_limits returns correct structure", {
  emb <- matrix(c(0, 1, 2, 3, 4, 5), ncol = 2)
  limits <- scVeloR:::get_embedding_limits(emb, expand = 0.1)
  
  expect_true(is.list(limits))
  expect_true(all(c("xlim", "ylim") %in% names(limits)))
  expect_length(limits$xlim, 2)
  expect_length(limits$ylim, 2)
  expect_true(limits$xlim[1] < limits$xlim[2])
  expect_true(limits$ylim[1] < limits$ylim[2])
})

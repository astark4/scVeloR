# Test Dynamics - ODE Solutions

# Test the core mathematical functions for splicing dynamics
# Reference: scvelo equations from Bergen et al., Nature Biotech 2020

test_that("unspliced dynamics at t=0 returns initial value", {
  # At t=0, u(0) = u0
  u0 <- 5
  result <- unspliced(tau = 0, u0 = u0, alpha = 2, beta = 0.5)
  expect_equal(result, u0, tolerance = 1e-10)
})

test_that("unspliced dynamics converges to steady state", {
  # At t->inf, u_inf = alpha/beta
  alpha <- 2
  beta <- 0.5
  u_inf_expected <- alpha / beta  # = 4
  
  # Test at large t
  result <- unspliced(tau = 100, u0 = 0, alpha = alpha, beta = beta)
  expect_equal(result, u_inf_expected, tolerance = 1e-6)
  
  # Test from different starting point
  result2 <- unspliced(tau = 100, u0 = 10, alpha = alpha, beta = beta)
  expect_equal(result2, u_inf_expected, tolerance = 1e-6)
})

test_that("unspliced dynamics follows correct ODE", {
  # du/dt = alpha - beta*u
  # Verify by numerical differentiation
  alpha <- 2
  beta <- 0.5
  u0 <- 1
  
  dt <- 1e-6
  tau <- 2.0
  
  u_t <- unspliced(tau, u0, alpha, beta)
  u_t_dt <- unspliced(tau + dt, u0, alpha, beta)
  
  # Numerical derivative
  du_dt_numerical <- (u_t_dt - u_t) / dt
  
  # Analytical derivative from ODE
  du_dt_analytical <- alpha - beta * u_t
  
  expect_equal(du_dt_numerical, du_dt_analytical, tolerance = 1e-4)
})

test_that("spliced dynamics at t=0 returns initial value", {
  s0 <- 3
  result <- spliced(tau = 0, s0 = s0, u0 = 1, alpha = 2, beta = 0.5, gamma = 0.3)
  expect_equal(result, s0, tolerance = 1e-10)
})

test_that("spliced dynamics converges to steady state", {
  # At t->inf, s_inf = alpha/gamma
  alpha <- 2
  beta <- 0.5
  gamma <- 0.3
  s_inf_expected <- alpha / gamma  # = 6.67
  
  result <- spliced(tau = 200, s0 = 0, u0 = 0, alpha = alpha, beta = beta, gamma = gamma)
  expect_equal(result, s_inf_expected, tolerance = 1e-4)
})

test_that("spliced dynamics follows correct ODE", {
  # ds/dt = beta*u - gamma*s
  alpha <- 2
  beta <- 0.5
  gamma <- 0.3
  u0 <- 1
  s0 <- 0.5
  
  dt <- 1e-6
  tau <- 2.0
  
  u_t <- unspliced(tau, u0, alpha, beta)
  s_t <- spliced(tau, s0, u0, alpha, beta, gamma)
  s_t_dt <- spliced(tau + dt, s0, u0, alpha, beta, gamma)
  
  # Numerical derivative
  ds_dt_numerical <- (s_t_dt - s_t) / dt
  
  # Analytical derivative from ODE
  ds_dt_analytical <- beta * u_t - gamma * s_t
  
  expect_equal(ds_dt_numerical, ds_dt_analytical, tolerance = 1e-4)
})

test_that("mRNA returns both components correctly", {
  tau <- seq(0, 5, length.out = 10)
  result <- mRNA(tau, u0 = 0, s0 = 0, alpha = 2, beta = 0.5, gamma = 0.3)
  
  expect_true(is.list(result))
  expect_true("u" %in% names(result))
  expect_true("s" %in% names(result))
  expect_length(result$u, 10)
  expect_length(result$s, 10)
  
  # Verify values match individual functions
  u_direct <- unspliced(tau, u0 = 0, alpha = 2, beta = 0.5)
  s_direct <- spliced(tau, s0 = 0, u0 = 0, alpha = 2, beta = 0.5, gamma = 0.3)
  
  expect_equal(result$u, u_direct, tolerance = 1e-10)
  expect_equal(result$s, s_direct, tolerance = 1e-10)
})

test_that("tau_inv_u recovers original time", {
  # Forward: compute u from tau
  # Backward: recover tau from u
  alpha <- 2
  beta <- 0.5
  u0 <- 0
  
  tau_original <- c(0.5, 1.0, 2.0, 5.0)
  u_values <- unspliced(tau_original, u0, alpha, beta)
  
  tau_recovered <- tau_inv_u(u_values, u0, alpha, beta)
  
  expect_equal(tau_recovered, tau_original, tolerance = 1e-4)
})

test_that("tau_inv_u handles edge cases", {
  # Test at steady state (u = alpha/beta)
  alpha <- 2
  beta <- 0.5
  u_inf <- alpha / beta
  
  # Should return a large time value
  tau <- tau_inv_u(u_inf - 0.001, u0 = 0, alpha = alpha, beta = beta)
  expect_true(is.finite(tau))
  expect_true(tau > 0)
})

test_that("steady_state returns correct values", {
  alpha <- 3
  beta <- 1.5
  gamma <- 0.5
  
  ss <- steady_state(alpha, beta, gamma)
  
  expect_equal(ss$u_inf, alpha / beta)  # = 2
  expect_equal(ss$s_inf, alpha / gamma)  # = 6
})

test_that("vectorize_params correctly separates on/off phases", {
  t <- c(0.5, 1.0, 1.5, 2.5, 3.0)  # t_ = 2.0
  t_ <- 2.0
  alpha <- 2
  beta <- 0.5
  gamma <- 0.3
  
  params <- vectorize_params(t, t_, alpha, beta, gamma)
  
  # Check structure
  expect_true(is.list(params))
  expect_true(all(c("tau", "alpha", "u0", "s0") %in% names(params)))
  
  # Cells before t_ should be in "on" phase (alpha > 0)
  on_phase <- t < t_
  off_phase <- t >= t_
  
  expect_true(all(params$alpha[on_phase] == alpha))
  expect_true(all(params$alpha[off_phase] == 0))
  
  # Initial states
  expect_true(all(params$u0[on_phase] == 0))
  expect_true(all(params$s0[on_phase] == 0))
})

test_that("compute_velocity_from_params gives correct values", {
  u <- c(1, 2, 3)
  s <- c(2, 4, 6)
  alpha <- 2
  beta <- 0.5
  gamma <- 0.3
  
  vel <- compute_velocity_from_params(u, s, alpha, beta, gamma)
  
  # velocity_u = alpha - beta*u
  expect_equal(vel$velocity_u, alpha - beta * u)
  
  # velocity_s = beta*u - gamma*s
  expect_equal(vel$velocity_s, beta * u - gamma * s)
})

test_that("dynamics are consistent during induction phase", {
  # During induction (alpha > 0), starting from (0,0)
  alpha <- 2
  beta <- 0.5
  gamma <- 0.3
  
  # Use longer time to approach steady state
  tau <- seq(0, 30, by = 0.5)
  
  result <- mRNA(tau, u0 = 0, s0 = 0, alpha, beta, gamma)
  
  # Both should increase initially
  expect_true(all(diff(result$u[1:5]) > 0))
  expect_true(all(diff(result$s[1:5]) > 0))
  
  # Should converge to steady state (with tolerance for numerical convergence)
  expect_equal(result$u[length(tau)], alpha/beta, tolerance = 0.01)
  expect_equal(result$s[length(tau)], alpha/gamma, tolerance = 0.1)
})

test_that("dynamics are consistent during repression phase", {
  # During repression (alpha = 0), starting from steady state
  alpha_on <- 2
  beta <- 0.5
  gamma <- 0.3
  
  # Starting point (steady state from induction)
  u0 <- alpha_on / beta
  s0 <- alpha_on / gamma
  
  # Use longer time for decay
  tau <- seq(0, 40, by = 0.5)
  
  # Repression: alpha = 0
  result <- mRNA(tau, u0 = u0, s0 = s0, alpha = 0, beta, gamma)
  
  # Unspliced should decrease monotonically
  expect_true(all(diff(result$u) <= 0))
  
  # Should converge to near zero (allow small tolerance)
  expect_equal(result$u[length(tau)], 0, tolerance = 0.01)
  expect_equal(result$s[length(tau)], 0, tolerance = 0.05)
})

test_that("safe_divide handles edge cases", {
  # Normal division
  expect_equal(safe_divide(10, 2), 5)
  
  # Division by zero
  expect_equal(safe_divide(10, 0), 0)
  
  # Vector division
  result <- safe_divide(c(1, 2, 3), c(1, 0, 3))
  expect_equal(result, c(1, 0, 1))
})

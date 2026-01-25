# Test Velocity Computation

test_that("dynamics functions are scientifically accurate", {
  # Test case 1: Verify ODE solution matches numerical integration
  # du/dt = alpha - beta*u
  
  alpha <- 2
  beta <- 0.5
  u0 <- 0
  
  # Analytical solution at t=2
  t <- 2
  u_analytical <- unspliced(t, u0, alpha, beta)
  
  # Numerical integration using Euler method
  dt <- 0.0001
  u_numerical <- u0
  for (i in seq_len(t / dt)) {
    du <- alpha - beta * u_numerical
    u_numerical <- u_numerical + du * dt
  }
  
  # Should match within reasonable tolerance
  expect_equal(u_analytical, u_numerical, tolerance = 0.01)
})

test_that("steady state velocity is zero", {
  # At steady state, du/dt = ds/dt = 0
  
  alpha <- 2
  beta <- 0.5
  gamma <- 0.3
  
  ss <- steady_state(alpha, beta, gamma)
  
  vel <- compute_velocity_from_params(
    u = ss$u_inf,
    s = ss$s_inf,
    alpha = alpha,
    beta = beta,
    gamma = gamma
  )
  
  expect_equal(vel$velocity_u, 0, tolerance = 1e-10)
  expect_equal(vel$velocity_s, 0, tolerance = 1e-10)
})

test_that("velocity direction is correct during induction", {
  # During induction (starting from 0), both u and s should increase
  
  alpha <- 2
  beta <- 0.5
  gamma <- 0.3
  
  # At t=0.5, far from steady state
  u <- unspliced(0.5, 0, alpha, beta)
  s <- spliced(0.5, 0, 0, alpha, beta, gamma)
  
  vel <- compute_velocity_from_params(u, s, alpha, beta, gamma)
  
  # Velocities should be positive
  expect_true(vel$velocity_u > 0)
  expect_true(vel$velocity_s > 0)
})

test_that("velocity direction is correct during repression", {
  # During repression (alpha=0), u always decreases
  # s behavior is more complex: initially vel_s ~ 0 at steady state, then becomes negative
  
  alpha_on <- 2
  beta <- 0.5
  gamma <- 0.3
  
  # At steady state: beta*u = gamma*s, so vel_s = 0
  # Test with u decayed more than s to get negative vel_s
  u0 <- 1.0  # Lower than steady state
  s0 <- 4.0  # Higher relative to new u
  
  vel <- compute_velocity_from_params(u0, s0, alpha = 0, beta, gamma)
  
  # velocity_u = -beta*u should be negative
  expect_true(vel$velocity_u < 0)
  # velocity_s = beta*u - gamma*s = 0.5 - 1.2 = -0.7 < 0
  expect_true(vel$velocity_s < 0)
})

test_that("tau_inv recovers time correctly for induction", {
  alpha <- 2
  beta <- 0.5
  gamma <- 0.3
  
  # Generate trajectory
  t_true <- seq(0.1, 5, by = 0.2)
  u <- unspliced(t_true, 0, alpha, beta)
  s <- spliced(t_true, 0, 0, alpha, beta, gamma)
  
  # Recover time
  t_recovered <- tau_inv_u(u, 0, alpha, beta)
  
  expect_equal(t_recovered, t_true, tolerance = 0.01)
})

test_that("assign_timepoints assigns cells to correct phase", {
  alpha <- 2
  beta <- 0.5
  gamma <- 0.3
  t_ <- 3  # Switching time
  
  # Generate cells in induction phase
  t_on <- seq(0.5, 2.5, length.out = 10)
  u_on <- unspliced(t_on, 0, alpha, beta)
  s_on <- spliced(t_on, 0, 0, alpha, beta, gamma)
  
  # Generate cells in repression phase
  u0_ <- unspliced(t_, 0, alpha, beta)
  s0_ <- spliced(t_, 0, 0, alpha, beta, gamma)
  t_off <- seq(0.5, 2.5, length.out = 10)
  u_off <- unspliced(t_off, u0_, 0, beta)
  s_off <- spliced(t_off, s0_, u0_, 0, beta, gamma)
  
  # Combine
  u <- c(u_on, u_off)
  s <- c(s_on, s_off)
  
  result <- assign_timepoints(u, s, alpha, beta, gamma, t_, u0_, s0_, mode = "inverse")
  
  # First 10 cells should be "on" phase (o = 1)
  expect_true(mean(result$o[1:10]) > 0.5)
  
  # Last 10 cells should be "off" phase (o = 0)
  expect_true(mean(result$o[11:20]) < 0.5)
})

test_that("vectorize_params handles phase transition correctly", {
  alpha <- 2
  beta <- 0.5
  gamma <- 0.3
  t_ <- 2
  
  # Cell times spanning both phases
  t <- c(0.5, 1.5, 2.5, 3.5)
  
  params <- vectorize_params(t, t_, alpha, beta, gamma)
  
  # Check tau (time within phase)
  expect_equal(params$tau[1], 0.5)  # On phase
  expect_equal(params$tau[2], 1.5)  # On phase
  expect_equal(params$tau[3], 0.5)  # Off phase (2.5 - 2)
  expect_equal(params$tau[4], 1.5)  # Off phase (3.5 - 2)
  
  # Check alpha values
  expect_equal(params$alpha[1], alpha)  # On phase
  expect_equal(params$alpha[2], alpha)  # On phase
  expect_equal(params$alpha[3], 0)      # Off phase
  expect_equal(params$alpha[4], 0)      # Off phase
})

test_that("mRNA trajectory is continuous at switching point", {
  alpha <- 2
  beta <- 0.5
  gamma <- 0.3
  t_ <- 3
  
  # Just before switching
  result_before <- mRNA(t_ - 0.01, 0, 0, alpha, beta, gamma)
  
  # At switching point
  result_at <- mRNA(t_, 0, 0, alpha, beta, gamma)
  
  # Continuity check
  expect_equal(result_before$u, result_at$u, tolerance = 0.01)
  expect_equal(result_before$s, result_at$s, tolerance = 0.01)
})

test_that("safe_divide prevents numerical issues", {
  # Test division by zero
  expect_equal(safe_divide(1, 0), 0)
  
  # Test normal division
  expect_equal(safe_divide(10, 2), 5)
  
  # Test vector operations
  result <- safe_divide(c(1, 2, 3), c(1, 0, 1))
  expect_equal(result, c(1, 0, 3))
  
  # Test with small denominator
  result2 <- safe_divide(1, 1e-15)
  expect_equal(result2, 0)  # Should return 0 for near-zero
})

test_that("dynamics handle edge case beta = gamma", {
  # When beta = gamma, the formula has a singularity
  # Should handle gracefully
  
  alpha <- 2
  beta <- 0.5
  gamma <- 0.5  # Same as beta
  
  tau <- seq(0, 5, length.out = 50)
  
  # Should not produce NaN or Inf
  s <- spliced(tau, 0, 0, alpha, beta, gamma)
  
  expect_true(all(is.finite(s)))
  expect_true(all(s >= 0))
  
  # Should still converge to steady state
  expect_equal(s[length(s)], alpha/gamma, tolerance = 0.1)
})

test_that("dynamics handle zero alpha (repression)", {
  beta <- 0.5
  gamma <- 0.3
  u0 <- 5
  s0 <- 10
  
  # Use longer time for complete decay
  tau <- seq(0, 40, length.out = 100)
  
  # Repression phase (alpha = 0)
  result <- mRNA(tau, u0, s0, alpha = 0, beta, gamma)
  
  # Should decay toward zero
  expect_true(all(diff(result$u) <= 0))  # Monotonic decrease
  expect_equal(result$u[length(tau)], 0, tolerance = 0.01)
  expect_equal(result$s[length(tau)], 0, tolerance = 0.06)
})

test_that("gamma/beta ratio determines steady state ratio", {
  alpha <- 2
  beta <- 0.5
  gamma <- 0.25  # gamma < beta
  
  ss <- steady_state(alpha, beta, gamma)
  
  # u_inf / s_inf = gamma / beta
  ratio <- ss$u_inf / ss$s_inf
  expected_ratio <- gamma / beta
  
  expect_equal(ratio, expected_ratio, tolerance = 1e-10)
})

test_that("trajectory shape matches expected biology", {
  # During induction:
  # - u rises quickly, then plateaus
  # - s rises more slowly, then plateaus
  # - s lags behind u
  
  alpha <- 2
  beta <- 0.5
  gamma <- 0.3
  
  t <- seq(0, 15, length.out = 200)
  result <- mRNA(t, 0, 0, alpha, beta, gamma)
  
  # Find time to reach 90% of steady state
  u_ss <- alpha / beta
  s_ss <- alpha / gamma
  
  t_u_90 <- min(which(result$u >= 0.9 * u_ss))
  t_s_90 <- min(which(result$s >= 0.9 * s_ss))
  
  # Unspliced should reach steady state faster than spliced
  expect_true(t_u_90 < t_s_90)
})

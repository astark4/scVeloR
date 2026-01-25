#' @title Splicing Dynamics
#' @description Core functions for solving splicing kinetics ODEs.
#' The model describes mRNA dynamics: du/dt = alpha - beta*u, ds/dt = beta*u - gamma*s
#' @name dynamics
NULL

#' Solve Unspliced mRNA Dynamics
#'
#' @description Analytical solution for unspliced mRNA:
#' u(t) = u0 * exp(-beta*t) + alpha/beta * (1 - exp(-beta*t))
#'
#' @param tau Time point(s)
#' @param u0 Initial unspliced abundance
#' @param alpha Transcription rate
#' @param beta Splicing rate
#'
#' @return Unspliced mRNA abundance at time tau
#' @export
unspliced <- function(tau, u0, alpha, beta) {
  # Handle edge case of beta = 0
  if (any(beta == 0)) {
    beta[beta == 0] <- 1e-10
  }
  
  exp_bt <- exp(-beta * tau)
  u0 * exp_bt + alpha / beta * (1 - exp_bt)
}

#' Solve Spliced mRNA Dynamics
#'
#' @description Analytical solution for spliced mRNA:
#' s(t) = s0*exp(-gamma*t) + alpha/gamma*(1-exp(-gamma*t)) + c*(exp(-gamma*t)-exp(-beta*t))
#' where c = (alpha - u0*beta) / (gamma - beta)
#'
#' @param tau Time point(s)
#' @param s0 Initial spliced abundance
#' @param u0 Initial unspliced abundance
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#'
#' @return Spliced mRNA abundance at time tau
#' @export
spliced <- function(tau, s0, u0, alpha, beta, gamma) {
  # Handle edge cases
  if (any(beta == 0)) beta[beta == 0] <- 1e-10
  if (any(gamma == 0)) gamma[gamma == 0] <- 1e-10
  
  exp_bt <- exp(-beta * tau)
  exp_gt <- exp(-gamma * tau)
  
  # Compute c = (alpha - u0*beta) / (gamma - beta)
  # Handle case where gamma == beta
  c <- safe_divide(alpha - u0 * beta, gamma - beta)
  
  s0 * exp_gt + alpha / gamma * (1 - exp_gt) + c * (exp_gt - exp_bt)
}

#' Solve mRNA Dynamics (Both Components)
#'
#' @description Solve both unspliced and spliced dynamics simultaneously.
#'
#' @param tau Time point(s)
#' @param u0 Initial unspliced abundance
#' @param s0 Initial spliced abundance
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#'
#' @return List with u and s vectors
#' @export
mRNA <- function(tau, u0, s0, alpha, beta, gamma) {
  u <- unspliced(tau, u0, alpha, beta)
  s <- spliced(tau, s0, u0, alpha, beta, gamma)
  list(u = u, s = s)
}

#' Inverse Function: Estimate tau from u
#'
#' @description Inverse of unspliced dynamics to estimate time from unspliced abundance:
#' tau = -1/beta * log((u - alpha/beta) / (u0 - alpha/beta))
#'
#' @param u Unspliced abundance
#' @param u0 Initial unspliced abundance
#' @param alpha Transcription rate
#' @param beta Splicing rate
#'
#' @return Estimated time tau
#' @export
tau_inv_u <- function(u, u0, alpha, beta) {
  # Handle edge cases
  if (any(beta == 0)) beta[beta == 0] <- 1e-10
  
  u_inf <- alpha / beta
  
  # Compute ratio with protection against invalid values
  ratio <- (u - u_inf) / (u0 - u_inf)
  
  # Clip ratio to valid range
  ratio <- pmax(ratio, 1e-10)
  ratio <- pmin(ratio, 1 - 1e-10)
  
  tau <- -1 / beta * log(ratio)
  
  # Clip tau to non-negative
  pmax(tau, 0)
}

#' Inverse Function: Estimate tau from u and s
#'
#' @description Estimate time using both unspliced and spliced information.
#' Uses unspliced-based estimate when gamma >= beta, otherwise uses spliced.
#'
#' @param u Unspliced abundance
#' @param s Spliced abundance
#' @param u0 Initial unspliced abundance
#' @param s0 Initial spliced abundance
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#'
#' @return Estimated time tau
#' @export
tau_inv <- function(u, s = NULL, u0 = NULL, s0 = NULL, 
                    alpha = NULL, beta = NULL, gamma = NULL) {
  
  # Determine which formula to use
  use_u_only <- is.null(s) || is.null(gamma) || (gamma >= beta)
  
  if (use_u_only) {
    return(tau_inv_u(u, u0, alpha, beta))
  }
  
  # Use spliced-based estimate when gamma < beta
  # x = s - beta/(gamma-beta) * u - x_inf
  # where x_inf = alpha/gamma - beta/(gamma-beta) * alpha/beta
  # tau = -1/gamma * log(x / x0)
  
  beta_prime <- beta / (gamma - beta)
  x_inf <- alpha / gamma - beta_prime * alpha / beta
  
  x <- s - beta_prime * u - x_inf
  x0 <- s0 - beta_prime * u0 - x_inf
  
  ratio <- x / x0
  ratio <- pmax(ratio, 1e-10)
  ratio <- pmin(ratio, 1 - 1e-10)
  
  tau <- -1 / gamma * log(ratio)
  
  # Combine with unspliced estimate where spliced is unreliable
  tau_u <- tau_inv_u(u, u0, alpha, beta)
  
  # Use unspliced estimate where gamma >= beta
  inv_u_mask <- gamma >= beta
  if (length(inv_u_mask) == 1) {
    if (inv_u_mask) {
      tau <- tau_u
    }
  } else {
    tau[inv_u_mask] <- tau_u[inv_u_mask]
  }
  
  pmax(tau, 0)
}

#' Safe Division
#'
#' @description Division with protection against division by zero.
#'
#' @param a Numerator
#' @param b Denominator
#' @param eps Epsilon for near-zero detection
#'
#' @return a/b with Inf replaced by 0
#' @keywords internal
safe_divide <- function(a, b, eps = 1e-10) {
  result <- a / b
  result[abs(b) < eps] <- 0
  result[!is.finite(result)] <- 0
  result
}

#' Assign Time Points to Cells
#'
#' @description Assign time points to cells by projecting onto dynamics trajectory
#' or using inverse formula.
#'
#' @param u Unspliced values
#' @param s Spliced values
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#' @param t_ Switching time
#' @param u0_ Unspliced at switching
#' @param s0_ Spliced at switching
#' @param mode Assignment mode: "projection" or "inverse"
#'
#' @return List with tau, tau_, and o (on/off indicator)
#' @keywords internal
assign_timepoints <- function(u, s, alpha, beta, gamma, 
                               t_ = NULL, u0_ = NULL, s0_ = NULL,
                               mode = "projection") {
  
  n <- length(u)
  
  # If switching parameters not provided, estimate them
  if (is.null(u0_) || is.null(s0_)) {
    result <- mRNA(t_, 0, 0, alpha, beta, gamma)
    u0_ <- result$u
    s0_ <- result$s
  }
  
  if (mode == "projection") {
    # Project observations onto trajectory curve
    num_points <- min(500, max(200, n %/% 5))
    
    # "On" phase trajectory (induction)
    t_on <- seq(0, t_, length.out = num_points)
    on_traj <- mRNA(t_on, 0, 0, alpha, beta, gamma)
    
    # "Off" phase trajectory (repression)
    t_off_max <- tau_inv_u(min(u[s > 0], na.rm = TRUE), u0_, 0, beta)
    if (!is.finite(t_off_max) || t_off_max <= 0) {
      t_off_max <- t_
    }
    t_off <- seq(0, t_off_max, length.out = num_points)[-1]
    off_traj <- mRNA(t_off, u0_, s0_, 0, beta, gamma)
    
    # Stack observations
    x_obs <- cbind(u, s)
    
    # Stack trajectory points
    x_on <- cbind(on_traj$u, on_traj$s)
    x_off <- cbind(off_traj$u, off_traj$s)
    
    # Find closest trajectory point for each observation
    tau <- numeric(n)
    tau_ <- numeric(n)
    o <- integer(n)
    
    for (i in seq_len(n)) {
      # Distance to "on" curve
      dist_on <- (x_on[, 1] - x_obs[i, 1])^2 + (x_on[, 2] - x_obs[i, 2])^2
      min_on_idx <- which.min(dist_on)
      min_on_dist <- dist_on[min_on_idx]
      
      # Distance to "off" curve
      dist_off <- (x_off[, 1] - x_obs[i, 1])^2 + (x_off[, 2] - x_obs[i, 2])^2
      min_off_idx <- which.min(dist_off)
      min_off_dist <- dist_off[min_off_idx]
      
      if (min_on_dist <= min_off_dist) {
        tau[i] <- t_on[min_on_idx]
        tau_[i] <- 0
        o[i] <- 1L  # "on" phase
      } else {
        tau[i] <- t_
        tau_[i] <- t_off[min_off_idx]
        o[i] <- 0L  # "off" phase
      }
    }
  } else {
    # Use inverse formula
    tau <- tau_inv(u, s, 0, 0, alpha, beta, gamma)
    tau <- pmin(tau, t_)
    
    tau_ <- tau_inv(u, s, u0_, s0_, 0, beta, gamma)
    tau_ <- pmax(tau_, 0)
    
    # Determine on/off state
    result_on <- mRNA(tau, 0, 0, alpha, beta, gamma)
    result_off <- mRNA(tau_, u0_, s0_, 0, beta, gamma)
    
    dist_on <- (u - result_on$u)^2 + (s - result_on$s)^2
    dist_off <- (u - result_off$u)^2 + (s - result_off$s)^2
    
    o <- as.integer(dist_on <= dist_off)
    tau_[o == 1] <- 0
    tau[o == 0] <- t_
  }
  
  list(tau = tau, tau_ = tau_, o = o)
}

#' Vectorize Parameters for Time-Dependent Computation
#'
#' @description Convert scalar parameters to cell-specific vectors based on
#' whether cell is in "on" or "off" phase.
#'
#' @param t Cell times
#' @param t_ Switching time
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#'
#' @return List with tau, alpha_vec, u0, s0
#' @keywords internal
vectorize_params <- function(t, t_, alpha, beta, gamma) {
  # Determine on/off phase
  o <- as.integer(t < t_)
  
  # Compute tau (time within current phase)
  tau <- t * o + (t - t_) * (1 - o)
  
  # Compute state at switching time
  result_switch <- mRNA(t_, 0, 0, alpha, beta, gamma)
  u0_switch <- result_switch$u
  s0_switch <- result_switch$s
  
  # Vectorize initial states and alpha
  u0 <- 0 * o + u0_switch * (1 - o)
  s0 <- 0 * o + s0_switch * (1 - o)
  alpha_vec <- alpha * o  # alpha = 0 in "off" phase
  
  list(tau = tau, alpha = alpha_vec, u0 = u0, s0 = s0)
}

#' Get Steady State Values
#'
#' @description Calculate steady state mRNA levels.
#'
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#'
#' @return List with u_inf and s_inf
#' @export
steady_state <- function(alpha, beta, gamma) {
  list(
    u_inf = alpha / beta,
    s_inf = alpha / gamma
  )
}

#' Compute Velocity from Parameters
#'
#' @description Compute RNA velocity: du/dt = alpha - beta*u, ds/dt = beta*u - gamma*s
#'
#' @param u Unspliced abundance
#' @param s Spliced abundance
#' @param alpha Transcription rate (0 in repression phase)
#' @param beta Splicing rate
#' @param gamma Degradation rate
#'
#' @return List with velocity_u and velocity_s
#' @export
compute_velocity_from_params <- function(u, s, alpha, beta, gamma) {
  velocity_u <- alpha - beta * u
  velocity_s <- beta * u - gamma * s
  list(velocity_u = velocity_u, velocity_s = velocity_s)
}

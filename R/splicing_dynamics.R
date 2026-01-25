#' @title Splicing Dynamics Class
#' @description Complete implementation of SplicingDynamics class
#' equivalent to scvelo.core.SplicingDynamics
#' @name SplicingDynamics
NULL

#' Splicing Dynamics Class (R6)
#' 
#' @description A class that encapsulates the analytical solutions for 
#' splicing kinetics ODEs: du/dt = alpha - beta*u, ds/dt = beta*u - gamma*s
#' 
#' @details This class provides:
#' - Analytical solutions for unspliced and spliced mRNA dynamics
#' - Inverse time estimation from observed states
#' - Trajectory generation for phase portraits
#' 
#' @export
SplicingDynamics <- R6::R6Class(
 "SplicingDynamics",
  
  public = list(
    #' @field alpha Transcription rate
    alpha = NULL,
    #' @field beta Splicing rate  
    beta = NULL,
    #' @field gamma Degradation rate
    gamma = NULL,
    #' @field initial_state Initial state [u0, s0]
    initial_state = NULL,
    #' @field alpha_ Basal transcription rate (for repression phase)
    alpha_ = 0,
    
    #' @description Initialize SplicingDynamics
    #' @param alpha Transcription rate
    #' @param beta Splicing rate
    #' @param gamma Degradation rate
    #' @param initial_state Initial state [u0, s0], default [0, 0]
    #' @param alpha_ Basal transcription rate
    initialize = function(alpha = 1, beta = 1, gamma = 1, 
                          initial_state = c(0, 0), alpha_ = 0) {
      self$alpha <- alpha
      self$beta <- beta
      self$gamma <- gamma
      self$initial_state <- initial_state
      self$alpha_ <- alpha_
    },
    
    #' @description Get unspliced mRNA at time t
    #' @param t Time point(s)
    #' @param u0 Initial unspliced (overrides initial_state if provided)
    #' @param alpha Transcription rate (overrides self$alpha if provided)
    #' @return Unspliced mRNA abundance
    get_unspliced = function(t, u0 = NULL, alpha = NULL) {
      if (is.null(u0)) u0 <- self$initial_state[1]
      if (is.null(alpha)) alpha <- self$alpha
      
      beta <- self$beta
      if (beta == 0) beta <- 1e-10
      
      exp_bt <- exp(-beta * t)
      u0 * exp_bt + alpha / beta * (1 - exp_bt)
    },
    
    #' @description Get spliced mRNA at time t
    #' @param t Time point(s)
    #' @param s0 Initial spliced (overrides initial_state if provided)
    #' @param u0 Initial unspliced (overrides initial_state if provided)
    #' @param alpha Transcription rate (overrides self$alpha if provided)
    #' @return Spliced mRNA abundance
    get_spliced = function(t, s0 = NULL, u0 = NULL, alpha = NULL) {
      if (is.null(s0)) s0 <- self$initial_state[2]
      if (is.null(u0)) u0 <- self$initial_state[1]
      if (is.null(alpha)) alpha <- self$alpha
      
      beta <- self$beta
      gamma <- self$gamma
      
      if (beta == 0) beta <- 1e-10
      if (gamma == 0) gamma <- 1e-10
      
      exp_bt <- exp(-beta * t)
      exp_gt <- exp(-gamma * t)
      
      # c = (alpha - u0*beta) / (gamma - beta)
      # Handle case where gamma == beta
      if (abs(gamma - beta) < 1e-10) {
        # L'Hopital's rule limit
        c_val <- (alpha - u0 * beta) * t * exp_gt
        return(s0 * exp_gt + alpha / gamma * (1 - exp_gt) + c_val)
      }
      
      c_val <- (alpha - u0 * beta) / (gamma - beta)
      s0 * exp_gt + alpha / gamma * (1 - exp_gt) + c_val * (exp_gt - exp_bt)
    },
    
    #' @description Get both unspliced and spliced at time t
    #' @param t Time point(s)
    #' @param stacked Return as matrix (TRUE) or list (FALSE)
    #' @return Matrix or list with u and s values
    get_solution = function(t, stacked = TRUE) {
      u <- self$get_unspliced(t)
      s <- self$get_spliced(t)
      
      if (stacked) {
        cbind(u, s)
      } else {
        list(u = u, s = s)
      }
    },
    
    #' @description Get velocity at time t
    #' @param t Time point(s)
    #' @return List with du/dt (wu) and ds/dt (ws)
    get_velocity = function(t) {
      u <- self$get_unspliced(t)
      s <- self$get_spliced(t)
      
      wu <- self$alpha - self$beta * u  # du/dt
      ws <- self$beta * u - self$gamma * s  # ds/dt
      
      list(wu = wu, ws = ws)
    },
    
    #' @description Get steady state values
    #' @return List with u_inf and s_inf
    get_steady_state = function() {
      list(
        u_inf = self$alpha / self$beta,
        s_inf = self$alpha / self$gamma
      )
    }
  )
)

#' Compute tau from unspliced (inverse function)
#' 
#' @description Estimate time from unspliced abundance using inverse formula:
#' tau = -1/beta * log((u - alpha/beta) / (u0 - alpha/beta))
#' 
#' @param u Unspliced abundance
#' @param u0 Initial unspliced abundance  
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @return Estimated time tau
#' @export
tau_inv_u <- function(u, u0, alpha, beta) {
  # Handle edge cases
  beta <- ifelse(beta == 0, 1e-10, beta)
  
  u_inf <- alpha / beta
  
  # Compute ratio with protection
  ratio <- (u - u_inf) / (u0 - u_inf)
  
  # Clip to valid range for log
  ratio <- pmax(ratio, 1e-10)
  ratio <- pmin(ratio, 1 - 1e-10)
  
  tau <- -1 / beta * log(ratio)
  
  # Return non-negative time
  pmax(tau, 0)
}

#' Compute tau from unspliced and spliced (inverse function)
#'
#' @description Estimate time using both unspliced and spliced information.
#' Uses unspliced-based estimate when gamma >= beta, otherwise uses spliced.
#' This matches exactly the Python scvelo implementation.
#'
#' @param u Unspliced abundance
#' @param s Spliced abundance
#' @param u0 Initial unspliced abundance
#' @param s0 Initial spliced abundance
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#' @return Estimated time tau
#' @export
tau_inv <- function(u, s = NULL, u0 = NULL, s0 = NULL,
                    alpha = NULL, beta = NULL, gamma = NULL) {
  
  # Determine which formula to use
  inv_u <- is.null(gamma) || is.null(s) || (gamma >= beta)
  inv_us <- !inv_u
  
  any_invu <- any(inv_u) || is.null(s)
  any_invus <- any(inv_us) && !is.null(s)
  
  tau <- rep(NA_real_, length(u))
  
  if (any_invus) {
    # tau_inv using spliced: tau = -1/gamma * log((s - beta'*u - x_inf) / (s0 - beta'*u0 - x_inf))
    beta_prime <- beta / (gamma - beta)
    x_inf <- alpha / gamma - beta_prime * (alpha / beta)
    
    ratio <- (s - beta_prime * u - x_inf) / (s0 - beta_prime * u0 - x_inf)
    ratio <- pmax(ratio, 1e-10)
    ratio <- pmin(ratio, 1 - 1e-10)
    
    tau <- -1 / gamma * log(ratio)
  }
  
  if (any_invu) {
    # tau_inv using unspliced
    tau_u <- tau_inv_u(u, u0, alpha, beta)
    
    if (any_invus) {
      # Combine: use unspliced where gamma >= beta
      if (length(inv_u) == 1) {
        if (inv_u) tau <- tau_u
      } else {
        tau[inv_u] <- tau_u[inv_u]
      }
    } else {
      tau <- tau_u
    }
  }
  
  pmax(tau, 0)
}

#' Unspliced mRNA analytical solution (vectorized)
#' 
#' @description u(t) = u0 * exp(-beta*t) + alpha/beta * (1 - exp(-beta*t))
#' @param tau Time point(s)
#' @param u0 Initial unspliced
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @return Unspliced abundance at time tau
#' @export
unspliced <- function(tau, u0, alpha, beta) {
  beta <- ifelse(beta == 0, 1e-10, beta)
  exp_bt <- exp(-beta * tau)
  u0 * exp_bt + alpha / beta * (1 - exp_bt)
}

#' Spliced mRNA analytical solution (vectorized)
#' 
#' @description Full analytical solution for spliced mRNA
#' @param tau Time point(s)
#' @param s0 Initial spliced
#' @param u0 Initial unspliced
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#' @return Spliced abundance at time tau
#' @export
spliced <- function(tau, s0, u0, alpha, beta, gamma) {
  beta <- ifelse(beta == 0, 1e-10, beta)
  gamma <- ifelse(gamma == 0, 1e-10, gamma)
  
  exp_bt <- exp(-beta * tau)
  exp_gt <- exp(-gamma * tau)
  
  # c = (alpha - u0*beta) / (gamma - beta)
  # Use safe division
  denom <- gamma - beta
  denom <- ifelse(abs(denom) < 1e-10, sign(denom) * 1e-10, denom)
  c_val <- (alpha - u0 * beta) / denom
  
  s0 * exp_gt + alpha / gamma * (1 - exp_gt) + c_val * (exp_gt - exp_bt)
}

#' mRNA dynamics (both components)
#' 
#' @description Solve both unspliced and spliced dynamics simultaneously
#' @param tau Time point(s)
#' @param u0 Initial unspliced
#' @param s0 Initial spliced
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#' @return List with u and s vectors
#' @export
mRNA <- function(tau, u0, s0, alpha, beta, gamma) {
  u <- unspliced(tau, u0, alpha, beta)
  s <- spliced(tau, s0, u0, alpha, beta, gamma)
  list(u = u, s = s)
}

#' Vectorize parameters for time-dependent computation
#'
#' @description Convert scalar parameters to cell-specific vectors based on
#' whether cell is in "on" (induction) or "off" (repression) phase.
#' This matches exactly the Python scvelo vectorize() function.
#'
#' @param t Cell times
#' @param t_ Switching time
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate (optional)
#' @param alpha_ Basal transcription (default 0)
#' @param u0 Initial unspliced (default 0)
#' @param s0 Initial spliced (default 0)
#' @param sorted Whether to sort by time
#' @return List with tau, alpha, u0, s0
#' @export
vectorize <- function(t, t_, alpha, beta, gamma = NULL, 
                      alpha_ = 0, u0 = 0, s0 = 0, sorted = FALSE) {
  
  # Determine on/off state
  o <- as.integer(t < t_)
  
  # tau is time within current phase
  tau <- t * o + (t - t_) * (1 - o)
  
  # Compute state at switching time
  u0_ <- unspliced(t_, u0, alpha, beta)
  if (!is.null(gamma)) {
    s0_ <- spliced(t_, s0, u0, alpha, beta, gamma)
  } else {
    s0_ <- spliced(t_, s0, u0, alpha, beta, beta / 2)
  }
  
  # Vectorize u0, s0 and alpha
  u0_vec <- u0 * o + u0_ * (1 - o)
  s0_vec <- s0 * o + s0_ * (1 - o)
  alpha_vec <- alpha * o + alpha_ * (1 - o)
  
  if (sorted) {
    idx <- order(t)
    tau <- tau[idx]
    alpha_vec <- alpha_vec[idx]
    u0_vec <- u0_vec[idx]
    s0_vec <- s0_vec[idx]
  }
  
  list(tau = tau, alpha = alpha_vec, u0 = u0_vec, s0 = s0_vec, o = o)
}

#' Assign time points to cells (projection or inverse method)
#' 
#' @description Assign time points by projecting observations onto dynamics 
#' trajectory or using inverse formula. Matches scvelo's assign_tau().
#' 
#' @param u Unspliced values
#' @param s Spliced values
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#' @param t_ Switching time
#' @param u0_ Unspliced at switching
#' @param s0_ Spliced at switching
#' @param assignment_mode "projection", "partial_projection", "full_projection", or NULL
#' @return List with tau, tau_, t_
#' @export
assign_tau <- function(u, s, alpha, beta, gamma, 
                       t_ = NULL, u0_ = NULL, s0_ = NULL,
                       assignment_mode = NULL) {
  
  n <- length(u)
  
  # Compute switching state if not provided
  if (is.null(u0_) || is.null(s0_)) {
    dynamics <- SplicingDynamics$new(alpha = alpha, beta = beta, gamma = gamma)
    sol <- dynamics$get_solution(t_, stacked = FALSE)
    u0_ <- sol$u
    s0_ <- sol$s
  }
  
  use_projection <- !is.null(assignment_mode) && 
    (assignment_mode %in% c("full_projection", "partial_projection") ||
    (assignment_mode == "projection" && beta < gamma))
  
  if (use_projection) {
    # Project observations onto trajectory curve
    x_obs <- cbind(u, s)
    
    # Estimate maximum time for off phase
    t0 <- tau_inv_u(min(u[s > 0], na.rm = TRUE), u0_, 0, beta)
    if (!is.finite(t0) || t0 <= 0) t0 <- t_
    
    # Generate trajectory points
    num <- as.integer(min(500, max(200, n / 5)))
    tpoints <- seq(0, t_, length.out = num)
    tpoints_ <- seq(0, t0, length.out = num)[-1]
    
    # On phase trajectory
    dynamics_on <- SplicingDynamics$new(alpha = alpha, beta = beta, gamma = gamma)
    xt <- dynamics_on$get_solution(tpoints)
    
    # Off phase trajectory
    dynamics_off <- SplicingDynamics$new(
      alpha = 0, beta = beta, gamma = gamma,
      initial_state = c(u0_, s0_)
    )
    xt_ <- dynamics_off$get_solution(tpoints_)
    
    # Find closest trajectory point for each observation
    tau <- numeric(n)
    tau_ <- numeric(n)
    
    for (i in seq_len(n)) {
      # Distance to on curve
      dist_on <- (xt[, 1] - x_obs[i, 1])^2 + (xt[, 2] - x_obs[i, 2])^2
      min_on_idx <- which.min(dist_on)
      
      # Distance to off curve
      dist_off <- (xt_[, 1] - x_obs[i, 1])^2 + (xt_[, 2] - x_obs[i, 2])^2
      min_off_idx <- which.min(dist_off)
      
      if (dist_on[min_on_idx] <= dist_off[min_off_idx]) {
        tau[i] <- tpoints[min_on_idx]
        tau_[i] <- 0
      } else {
        tau[i] <- t_
        tau_[i] <- tpoints_[min_off_idx]
      }
    }
  } else {
    # Use inverse formula
    tau <- tau_inv(u, s, 0, 0, alpha, beta, gamma)
    tau <- pmin(tau, t_)
    
    tau_ <- tau_inv(u, s, u0_, s0_, 0, beta, gamma)
    tau_ <- pmax(tau_, 0)
    
    # Clip tau_ to max observed
    if (any(s > 0)) {
      tau_ <- pmin(tau_, max(tau_[s > 0], na.rm = TRUE))
    }
  }
  
  list(tau = tau, tau_ = tau_, t_ = t_)
}

#' Clipped logarithm (safe log)
#' 
#' @description Log with clipping to avoid -Inf
#' @param x Input values
#' @param eps Minimum value
#' @return log(max(x, eps))
#' @keywords internal
clipped_log <- function(x, eps = 1e-10) {
  log(pmax(x, eps))
}

#' Invert (safe division)
#' 
#' @description Division that handles near-zero denominators
#' @param x Denominator
#' @param eps Threshold
#' @return 1/x with protection
#' @keywords internal
invert <- function(x, eps = 1e-10) {
  x <- ifelse(abs(x) < eps, sign(x) * eps, x)
  1 / x
}

#' Linear regression (simple)
#' 
#' @description Linear regression through origin: gamma = sum(u*s) / sum(s^2)
#' @param u Unspliced (dependent)
#' @param s Spliced (independent)
#' @return Slope coefficient
#' @keywords internal
linreg <- function(u, s) {
  if (inherits(s, "dgCMatrix") || inherits(s, "sparseMatrix")) {
    ss <- sum(s@x^2)
    us <- sum(s@x * u[which(s@x != 0)])
  } else {
    # Handle NA and non-finite values
    valid <- is.finite(u) & is.finite(s) & (s != 0 | u != 0)
    if (sum(valid) < 2) return(1)
    u <- u[valid]
    s <- s[valid]
    ss <- sum(s^2, na.rm = TRUE)
    us <- sum(s * u, na.rm = TRUE)
  }
  if (ss == 0) return(1)
  us / ss
}

#' Convolve with weights
#' 
#' @description Multiply by weights (sparse-aware)
#' @param x Values
#' @param weights Weights (can be sparse or logical)
#' @return Weighted values
#' @keywords internal
convolve_weights <- function(x, weights = NULL) {
  if (is.null(weights)) {
    return(x)
  }
  if (inherits(weights, "sparseMatrix")) {
    as.vector(weights %*% x)
  } else if (is.logical(weights)) {
    # Filter by logical weights
    x[weights]
  } else {
    weights * x
  }
}

#' Normalize (row or column)
#' 
#' @param X Matrix or vector
#' @param axis 0 for columns, 1 for rows
#' @param min_confidence Minimum value to add to sum
#' @return Normalized values
#' @keywords internal
normalize <- function(X, axis = 0, min_confidence = NULL) {
  if (is.vector(X) || is.null(dim(X))) {
    X_sum <- sum(X, na.rm = TRUE)
  } else if (axis == 0) {
    X_sum <- colSums(X, na.rm = TRUE)
  } else {
    X_sum <- rowSums(X, na.rm = TRUE)
  }
  
  if (!is.null(min_confidence)) {
    X_sum <- X_sum + min_confidence
  }
  X_sum[X_sum == 0] <- 1
  
  if (is.vector(X) || is.null(dim(X))) {
    X / X_sum
  } else if (axis == 0) {
    sweep(X, 2, X_sum, "/")
  } else {
    sweep(X, 1, X_sum, "/")
  }
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

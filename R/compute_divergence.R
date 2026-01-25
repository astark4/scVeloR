#' @title Compute Divergence
#' @description Complete implementation of compute_divergence() function
#' equivalent to scvelo.tools._em_model_utils.compute_divergence
#' This is the CORE function for dynamical model analysis.
#' @name compute_divergence
NULL

#' Compute Divergence Between ODE and Observations
#'
#' @description Estimates the divergence of ODE to observations using various
#' metrics including distance, MSE, likelihood, and log-likelihood.
#' This function is the heart of the dynamical model.
#'
#' @param u Unspliced abundance (scaled)
#' @param s Spliced abundance
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#' @param scaling Scaling factor (default 1)
#' @param t_ Switching time
#' @param u0_ Unspliced at switching
#' @param s0_ Spliced at switching
#' @param tau Time in induction phase
#' @param tau_ Time in repression phase
#' @param std_u Standard deviation of unspliced
#' @param std_s Standard deviation of spliced
#' @param normalized Whether to normalize likelihoods
#' @param mode Output mode (see details)
#' @param assignment_mode Time assignment mode
#' @param var_scale Whether to scale by variance
#' @param kernel_width Kernel width for smoothing
#' @param fit_steady_states Whether to fit steady states
#' @param connectivities Connectivity matrix for smoothing
#' @param constraint_time_increments Whether to constrain time increments
#' @param reg_time Regularization time
#' @param reg_par Regularization parameter
#' @param min_confidence Minimum confidence value
#' @param pval_steady P-value for steady state
#' @param steady_u Steady state unspliced
#' @param steady_s Steady state spliced
#' @param noise_model Noise model ("chi" or "normal")
#' @param time_connectivities Time-based connectivities
#' @param clusters Cluster assignments
#' @param ... Additional arguments
#'
#' @details
#' Available modes:
#' - "distance": Squared Euclidean distance
#' - "mse": Mean squared error
#' - "likelihood": Gaussian likelihood
#' - "nll": Negative log-likelihood
#' - "loglikelihood": Log-likelihood
#' - "confidence": Confidence score
#' - "soft", "soft_eval": Soft state assignment with likelihoods
#' - "hard", "hard_eval": Hard state assignment
#' - "hardsoft", "hardsoft_eval": Mixed assignment
#' - "soft_state": Soft state difference
#' - "hard_state": Hard state indicator
#' - "steady_state": Steady state likelihood
#' - "assign_timepoints", "time": Time assignment
#' - "tau": Return tau values
#' - "dists": Return distance components
#' - "distx": Return total distance
#' - "velocity": Compute velocity
#' - "velocity_residuals": Velocity residuals
#' - "soft_velocity": Soft velocity
#' - "gene_likelihood": Gene-specific likelihood
#' - "unspliced_dists": Unspliced distances only
#' - "outside_of_trajectory": Boolean for outside trajectory
#'
#' @return Result depends on mode (see details)
#' @export
compute_divergence <- function(u, s, alpha, beta, gamma, scaling = 1,
                                t_ = NULL, u0_ = NULL, s0_ = NULL,
                                tau = NULL, tau_ = NULL,
                                std_u = 1, std_s = 1,
                                normalized = FALSE,
                                mode = "distance",
                                assignment_mode = NULL,
                                var_scale = FALSE,
                                kernel_width = NULL,
                                fit_steady_states = TRUE,
                                connectivities = NULL,
                                constraint_time_increments = TRUE,
                                reg_time = NULL, reg_par = NULL,
                                min_confidence = NULL,
                                pval_steady = NULL,
                                steady_u = NULL, steady_s = NULL,
                                noise_model = "chi",
                                time_connectivities = NULL,
                                clusters = NULL,
                                ...) {
  
  # Compute u0_, s0_ at switching time if not provided
  if (is.null(u0_) || is.null(s0_)) {
    dynamics <- SplicingDynamics$new(alpha = alpha, beta = beta, gamma = gamma)
    sol <- dynamics$get_solution(t_, stacked = FALSE)
    u0_ <- sol$u
    s0_ <- sol$s
  }
  
  # Assign time points if not provided
  if (is.null(tau) || is.null(tau_) || is.null(t_)) {
    result <- assign_tau(u, s, alpha, beta, gamma, t_, u0_, s0_, assignment_mode)
    tau <- result$tau
    tau_ <- result$tau_
    t_ <- result$t_
  }
  
  # Scale std_u
  std_u <- std_u / scaling
  
  # Adjust time increments to avoid meaningless jumps
  if (constraint_time_increments) {
    # Compute distances for induction phase
    dynamics_on <- SplicingDynamics$new(alpha = alpha, beta = beta, gamma = gamma)
    sol_on <- dynamics_on$get_solution(tau, stacked = FALSE)
    ut <- sol_on$u
    st <- sol_on$s
    
    # Compute distances for repression phase
    dynamics_off <- SplicingDynamics$new(
      alpha = 0, beta = beta, gamma = gamma,
      initial_state = c(u0_, s0_)
    )
    sol_off <- dynamics_off$get_solution(tau_, stacked = FALSE)
    ut_ <- sol_off$u
    st_ <- sol_off$s
    
    # Compute residuals
    distu <- (u - ut) / std_u
    distu_ <- (u - ut_) / std_u
    dists <- (s - st) / std_s
    dists_ <- (s - st_) / std_s
    
    # Compute distances
    res <- rbind(distu_^2 + dists_^2, distu^2 + dists^2)
    
    # Apply connectivities smoothing
    if (!is.null(connectivities) && !isFALSE(connectivities)) {
      if (is.matrix(res)) {
        res <- rbind(
          as.vector(connectivities %*% res[1, ]),
          as.vector(connectivities %*% res[2, ])
        )
      }
    }
    
    # Determine on/off state
    o <- apply(res, 2, which.min)
    off <- o == 1
    on <- o == 2
    
    # Adjust increments
    if (any(on) && any(off)) {
      adjusted <- adjust_increments(tau[on], tau_[off])
      tau[on] <- adjusted$tau
      tau_[off] <- adjusted$tau_
    } else if (any(on)) {
      tau[on] <- adjust_increments(tau[on])$tau
    } else if (any(off)) {
      tau_[off] <- adjust_increments(tau_[off])$tau
    }
  }
  
  # Compute induction/repression state distances
  dynamics_on <- SplicingDynamics$new(alpha = alpha, beta = beta, gamma = gamma)
  sol_on <- dynamics_on$get_solution(tau, stacked = FALSE)
  ut <- sol_on$u
  st <- sol_on$s
  
  dynamics_off <- SplicingDynamics$new(
    alpha = 0, beta = beta, gamma = gamma,
    initial_state = c(u0_, s0_)
  )
  sol_off <- dynamics_off$get_solution(tau_, stacked = FALSE)
  ut_ <- sol_off$u
  st_ <- sol_off$s
  
  # Flatten if needed
  if (is.matrix(ut) && ncol(ut) == 1) ut <- as.vector(ut)
  if (is.matrix(st) && ncol(st) == 1) st <- as.vector(st)
  if (is.matrix(ut_) && ncol(ut_) == 1) ut_ <- as.vector(ut_)
  if (is.matrix(st_) && ncol(st_) == 1) st_ <- as.vector(st_)
  
  # Compute residuals
  distu <- (u - ut) / std_u
  distu_ <- (u - ut_) / std_u
  dists <- (s - st) / std_s
  dists_ <- (s - st_) / std_s
  
  # Handle special mode: unspliced_dists
  if (mode == "unspliced_dists") {
    return(list(distu = distu, distu_ = distu_))
  }
  
  # Handle special mode: outside_of_trajectory
  if (mode == "outside_of_trajectory") {
    return(sign(distu) * sign(distu_) == 1)
  }
  
  # Compute squared distances
  distx <- distu^2 + dists^2
  distx_ <- distu_^2 + dists_^2
  
  # Default values
  res <- rbind(distx_, distx)
  varx <- 1
  
  # Normal noise model with variance scaling
  if (noise_model == "normal" && var_scale) {
    o <- apply(rbind(distx_, distx), 2, which.min)
    varu <- var(distu * (o == 2) + distu_ * (o == 1), na.rm = TRUE)
    vars <- var(dists * (o == 2) + dists_ * (o == 1), na.rm = TRUE)
    
    distx <- distu^2 / varu + dists^2 / vars
    distx_ <- distu_^2 / varu + dists_^2 / vars
    
    varx <- varu * vars
    std_u <- std_u * sqrt(varu)
    std_s <- std_s * sqrt(vars)
  }
  
  # Compute steady state distances
  if (fit_steady_states) {
    u_inf <- alpha / beta
    s_inf <- alpha / gamma
    distx_steady <- ((u - u_inf) / std_u)^2 + ((s - s_inf) / std_s)^2
    distx_steady_ <- (u / std_u)^2 + (s / std_s)^2
    res <- rbind(distx_, distx, distx_steady_, distx_steady)
  }
  
  # Apply connectivities smoothing
  if (!is.null(connectivities) && !isFALSE(connectivities)) {
    res <- t(apply(res, 1, function(r) as.vector(connectivities %*% r)))
  }
  
  # Compute variances (chi model)
  if (noise_model == "chi") {
    if (var_scale) {
      o <- apply(rbind(distx_, distx), 2, which.min)
      dist <- distx * (o == 2) + distx_ * (o == 1)
      sign_dist <- sign(dists * (o == 2) + dists_ * (o == 1))
      varx <- mean(dist, na.rm = TRUE) - mean(sign_dist * sqrt(dist), na.rm = TRUE)^2
      
      if (!is.null(kernel_width)) {
        varx <- varx * kernel_width^2
      }
      res <- res / varx
    } else if (!is.null(kernel_width)) {
      res <- res / kernel_width^2
    }
  }
  
  # Regularization with time
  if (!is.null(reg_time) && length(reg_time) == length(distu_)) {
    o <- apply(res, 2, which.min)
    t_max <- (t_ + tau_) * (o == 1)
    t_max <- t_max / max(t_max, na.rm = TRUE)
    reg_time <- reg_time / max(reg_time, na.rm = TRUE)
    
    dist_tau <- (tau - reg_time)^2
    dist_tau_ <- (tau_ + t_ - reg_time)^2
    mu_res <- rowMeans(res, na.rm = TRUE)
    
    if (!is.null(reg_par)) {
      mu_res <- mu_res * reg_par
    }
    
    res[1, ] <- res[1, ] + dist_tau_ * mu_res[1]
    res[2, ] <- res[2, ] + dist_tau * mu_res[2]
    if (fit_steady_states) {
      res[3, ] <- res[3, ] + dist_tau * mu_res[2]
      res[4, ] <- res[4, ] + dist_tau_ * mu_res[1]
    }
  }
  
  # Handle different modes
  if (mode == "tau") {
    return(list(tau = tau, tau_ = tau_))
  }
  
  if (mode == "likelihood") {
    res <- 1 / (2 * pi * sqrt(varx)) * exp(-0.5 * res)
    if (normalized) {
      res <- normalize(res, axis = 0, min_confidence = min_confidence)
    }
    return(res)
  }
  
  if (mode == "nll") {
    res <- log(2 * pi * sqrt(varx)) + 0.5 * res
    if (normalized) {
      res <- normalize(res, axis = 0, min_confidence = min_confidence)
    }
    return(res)
  }
  
  if (mode == "confidence") {
    res <- rbind(res[1, ], res[2, ])
    res <- 1 / (2 * pi * sqrt(varx)) * exp(-0.5 * res)
    if (normalized) {
      res <- normalize(res, axis = 0, min_confidence = min_confidence)
    }
    res <- median(apply(res, 2, max) - (colSums(res) - apply(res, 2, max)), na.rm = TRUE)
    return(res)
  }
  
  if (mode %in% c("soft_eval", "soft")) {
    res <- 1 / (2 * pi * sqrt(varx)) * exp(-0.5 * res)
    if (normalized) {
      res <- normalize(res, axis = 0, min_confidence = min_confidence)
    }
    o_ <- res[1, ]
    o <- res[2, ]
    res <- rbind(o_, o, ut * o + ut_ * o_, st * o + st_ * o_)
    return(res)
  }
  
  if (mode %in% c("hardsoft_eval", "hardsoft")) {
    res <- 1 / (2 * pi * sqrt(varx)) * exp(-0.5 * res)
    if (normalized) {
      res <- normalize(res, axis = 0, min_confidence = min_confidence)
    }
    o <- apply(res, 2, which.max)
    o_ <- (o == 1) * res[1, ]
    o <- (o == 2) * res[2, ]
    res <- rbind(o_, o, ut * o + ut_ * o_, st * o + st_ * o_)
    return(res)
  }
  
  if (mode %in% c("hard_eval", "hard")) {
    res <- 1 / (2 * pi * sqrt(varx)) * exp(-0.5 * res)
    if (normalized) {
      res <- normalize(res, axis = 0, min_confidence = min_confidence)
    }
    o <- apply(res, 2, which.max)
    o_ <- o == 1
    o <- o == 2
    res <- rbind(as.numeric(o_), as.numeric(o), ut * o + ut_ * o_, st * o + st_ * o_)
    return(res)
  }
  
  if (mode == "soft_state") {
    res <- 1 / (2 * pi * sqrt(varx)) * exp(-0.5 * res)
    if (normalized) {
      res <- normalize(res, axis = 0, min_confidence = min_confidence)
    }
    return(res[2, ] - res[1, ])
  }
  
  if (mode == "hard_state") {
    return(apply(res, 2, which.min))
  }
  
  if (mode == "steady_state") {
    res <- 1 / (2 * pi * sqrt(varx)) * exp(-0.5 * res)
    if (normalized) {
      res <- normalize(res, axis = 0, min_confidence = min_confidence)
    }
    return(res[3, ] + res[4, ])
  }
  
  if (mode %in% c("assign_timepoints", "time")) {
    o <- apply(res, 2, which.min)
    
    tau_[o != 1] <- 0
    tau[o != 2] <- 0
    
    # Handle steady states
    o[o == 3] <- 2
    o[o == 4] <- 1
    
    t <- tau * (o == 2) + (tau_ + t_) * (o == 1)
    
    if (mode == "assign_timepoints") {
      return(list(t = t, tau = tau, o = as.integer(o == 2)))
    } else {
      return(t)
    }
  }
  
  if (mode == "dists") {
    o <- apply(res, 2, which.min)
    tau_[o != 1] <- 0
    tau[o != 2] <- 0
    o[o == 3] <- 2
    o[o == 4] <- 1
    
    distu_out <- distu * (o == 2) + distu_ * (o == 1)
    dists_out <- dists * (o == 2) + dists_ * (o == 1)
    return(list(distu = distu_out, dists = dists_out))
  }
  
  if (mode == "distx") {
    o <- apply(res, 2, which.min)
    tau_[o != 1] <- 0
    tau[o != 2] <- 0
    o[o == 3] <- 2
    o[o == 4] <- 1
    
    distu_out <- distu * (o == 2) + distu_ * (o == 1)
    dists_out <- dists * (o == 2) + dists_ * (o == 1)
    return(distu_out^2 + dists_out^2)
  }
  
  if (mode == "gene_likelihood") {
    o <- apply(res, 2, which.min)
    tau_[o != 1] <- 0
    tau[o != 2] <- 0
    o[o == 3] <- 2
    o[o == 4] <- 1
    
    distu_out <- distu * (o == 2) + distu_ * (o == 1)
    dists_out <- dists * (o == 2) + dists_ * (o == 1)
    
    # Filter by expression level
    idx <- as.numeric((u > max(u, na.rm = TRUE) / 5) & (s > max(s, na.rm = TRUE) / 5))
    idx[idx == 0] <- NA
    distu_out <- distu_out * idx
    dists_out <- dists_out * idx
    
    distx_out <- distu_out^2 + dists_out^2
    
    # Compute variance
    varx <- mean(distx_out, na.rm = TRUE) - 
            mean(sign(dists_out) * sqrt(distx_out), na.rm = TRUE)^2
    
    if (!is.null(clusters)) {
      cats <- unique(clusters)
      res_list <- lapply(cats, function(cat) {
        idx_cat <- clusters == cat
        distx_cat <- distu_out[idx_cat]^2 + dists_out[idx_cat]^2
        distx_sum <- sum(distx_cat, na.rm = TRUE)
        
        n <- sum(!is.na(distx_cat)) - length(distx_cat) * 0.01
        n <- max(n, 2)
        
        if (n < max(n, na.rm = TRUE) / 5) {
          distx_sum <- max(distx_sum, na.rm = TRUE)
        }
        
        ll <- -1 / 2 / n * distx_sum / varx - 1 / 2 * log(2 * pi * varx)
        ll[distx_sum == 0] <- NA
        exp(ll)
      })
      return(do.call(rbind, res_list))
    } else {
      n <- max(length(distu_out) - length(distu_out) * 0.01, 2)
      ll <- -1 / 2 / n * sum(distx_out, na.rm = TRUE) / varx
      ll <- ll - 1 / 2 * log(2 * pi * varx)
      return(exp(ll))
    }
  }
  
  if (mode == "velocity") {
    res <- 1 / (2 * pi * sqrt(varx)) * exp(-0.5 * res)
    
    if (fit_steady_states) {
      res <- rbind(res[1:4, ], matrix(1e-6, nrow = 1, ncol = ncol(res)))
    } else {
      res <- rbind(res[1:2, ], matrix(ifelse(is.null(min_confidence), 1e-6, min_confidence), 
                                      nrow = 1, ncol = ncol(res)))
    }
    
    if (normalized) {
      res <- normalize(res, axis = 0, min_confidence = min_confidence)
    }
    
    o <- apply(res, 2, which.max)
    o_ <- o == 1
    o <- o == 2
    t <- tau * o + (tau_ + t_) * o_
    
    if (!is.null(time_connectivities) && !isFALSE(time_connectivities)) {
      if (isTRUE(time_connectivities)) time_connectivities <- connectivities
      t <- as.vector(time_connectivities %*% t)
      o <- (o < 3) * (t < t_)
      o_ <- (o < 3) * (t >= t_)
      
      vparams <- vectorize(t, t_, alpha, beta, gamma)
      dynamics_t <- SplicingDynamics$new(alpha = vparams$alpha, beta = beta, gamma = gamma)
      # Note: This would need to handle vector alpha properly
      ut <- unspliced(vparams$tau, vparams$u0, vparams$alpha, beta)
      st <- spliced(vparams$tau, vparams$s0, vparams$u0, vparams$alpha, beta, gamma)
      ut_ <- ut
      st_ <- st
    }
    
    ut <- ut * o + ut_ * o_
    st <- st * o + st_ * o_
    alpha_vec <- alpha * o
    
    vt <- ut * beta - st * gamma  # ds/dt
    wt <- (alpha_vec - beta * ut) * scaling  # du/dt
    
    vt <- pmax(vt, -s)
    wt <- pmax(wt, -u * scaling)
    
    if (is.vector(vt)) {
      vt <- matrix(vt, ncol = 1)
      wt <- matrix(wt, ncol = 1)
    }
    
    return(list(vt = vt, wt = wt))
  }
  
  if (mode == "velocity_residuals") {
    res <- 1 / (2 * pi * sqrt(varx)) * exp(-0.5 * res)
    
    if (fit_steady_states) {
      res <- rbind(res[1:4, ], matrix(1e-6, nrow = 1, ncol = ncol(res)))
    } else {
      res <- rbind(res[1:2, ], matrix(ifelse(is.null(min_confidence), 1e-6, min_confidence),
                                      nrow = 1, ncol = ncol(res)))
    }
    
    if (normalized) {
      res <- normalize(res, axis = 0, min_confidence = min_confidence)
    }
    
    o <- apply(res, 2, which.max)
    o_ <- o == 1
    o <- o == 2
    t <- tau * o + (tau_ + t_) * o_
    
    if (!is.null(time_connectivities) && !isFALSE(time_connectivities)) {
      if (isTRUE(time_connectivities)) time_connectivities <- connectivities
      t <- as.vector(time_connectivities %*% t)
      o <- (o < 3) * (t < t_)
      o_ <- (o < 3) * (t >= t_)
    }
    
    alpha_vec <- alpha * o
    
    vt <- u * beta - s * gamma  # ds/dt using observed values
    wt <- (alpha_vec - beta * u) * scaling  # du/dt
    
    vt <- pmax(vt, -s)
    wt <- pmax(wt, -u * scaling)
    
    return(list(vt = vt, wt = wt))
  }
  
  if (mode == "soft_velocity") {
    res <- 1 / (2 * pi * sqrt(varx)) * exp(-0.5 * res)
    if (normalized) {
      res <- normalize(res, axis = 0, min_confidence = min_confidence)
    }
    
    o_ <- res[1, ]
    o <- res[2, ]
    ut <- ut * o + ut_ * o_
    st <- st * o + st_ * o_
    alpha_vec <- alpha * o
    
    vt <- ut * beta - st * gamma
    wt <- (alpha_vec - beta * ut) * scaling
    
    return(list(vt = vt, wt = wt))
  }
  
  # Default: return distance matrix
  return(res)
}

#' Adjust time increments to avoid meaningless jumps
#' 
#' @description Clips large time jumps based on Poisson statistics
#' @param tau Time values for on phase
#' @param tau_ Time values for off phase (optional)
#' @return Adjusted time values
#' @keywords internal
adjust_increments <- function(tau, tau_ = NULL) {
  tau_new <- tau
  tau_ord <- sort(tau_new)
  dtau <- c(0, diff(tau_ord))
  
  if (is.null(tau_)) {
    ub <- 3 * quantile(dtau, 0.995, na.rm = TRUE)
  } else {
    tau_new_ <- tau_
    tau_ord_ <- sort(tau_new_)
    dtau_ <- c(0, diff(tau_ord_))
    
    ub <- 3 * quantile(c(dtau, dtau_), 0.995, na.rm = TRUE)
    
    idx <- which(dtau_ > ub)
    for (i in idx) {
      ti <- tau_ord_[i]
      dti <- dtau_[i]
      tau_new_[tau_ >= ti] <- tau_new_[tau_ >= ti] - dti
    }
  }
  
  idx <- which(dtau > ub)
  for (i in idx) {
    ti <- tau_ord[i]
    dti <- dtau[i]
    tau_new[tau >= ti] <- tau_new[tau >= ti] - dti
  }
  
  if (is.null(tau_)) {
    list(tau = tau_new)
  } else {
    list(tau = tau_new, tau_ = tau_new_)
  }
}

#' Root time computation
#' 
#' @description Compute rooted time relative to a root cell
#' @param t Time matrix (cells x genes)
#' @param root Root cell index
#' @return List with t_rooted and t_switch
#' @export
root_time <- function(t, root = NULL) {
  # Remove NaN columns
  nans <- apply(t, 2, function(x) any(is.na(x)))
  if (any(nans)) {
    t <- t[, !nans, drop = FALSE]
  }
  
  t_root <- if (is.null(root)) 0 else t[root, ]
  
  # Convert to matrix if vector
  if (is.vector(t_root)) t_root <- matrix(t_root, nrow = 1)
  if (is.vector(t)) t <- matrix(t, nrow = 1)
  
  o <- (t >= t_root[rep(1, nrow(t)), ]) * 1
  t_after <- (t - t_root[rep(1, nrow(t)), ]) * o
  t_origin <- apply(t_after, 2, max, na.rm = TRUE)
  t_before <- sweep((1 - o), 2, t_origin, "*") + t * (1 - o)
  
  t_switch <- apply(t_before, 2, min, na.rm = TRUE)
  t_rooted <- t_after + t_before
  
  list(t_rooted = t_rooted, t_switch = t_switch)
}

#' Compute shared time across genes
#' 
#' @description Compute a gene-shared latent time by coupling gene-specific times
#' @param t Time matrix (cells x genes)
#' @param perc Percentiles to use
#' @param norm Whether to normalize result
#' @return Shared time vector
#' @export
compute_shared_time <- function(t, perc = NULL, norm = TRUE) {
  # Remove NaN columns
  nans <- apply(t, 2, function(x) any(is.na(x)))
  if (any(nans)) {
    t <- t[, !nans, drop = FALSE]
  }
  
  # Shift to start at 0
  t <- t - min(t, na.rm = TRUE)
  
  # Compute percentiles
  if (is.null(perc)) {
    perc <- c(15, 20, 25, 30, 35)
  }
  
  tx_list <- apply(t, 1, function(x) quantile(x, perc / 100, na.rm = TRUE))
  tx_max <- apply(tx_list, 1, max, na.rm = TRUE)
  tx_max[tx_max == 0] <- 1
  tx_list <- tx_list / tx_max
  
  # Find best percentile based on MSE to linear
  mse <- apply(tx_list, 1, function(tx) {
    tx_sorted <- sort(tx)
    linx <- seq(0, 1, length.out = length(tx))
    sum((tx_sorted - linx)^2)
  })
  
  idx_best <- order(mse)[1:2]
  t_shared <- colSums(tx_list[idx_best, , drop = FALSE])
  
  if (norm) {
    t_shared <- t_shared / max(t_shared, na.rm = TRUE)
  }
  
  t_shared
}

#' @title Base Dynamics Class
#' @description Complete implementation of BaseDynamics class
#' equivalent to scvelo.tools._em_model_utils.BaseDynamics
#' @name BaseDynamics
NULL

#' Base Dynamics Class (R6)
#' 
#' @description Base class for dynamics recovery that handles:
#' - Data extraction and weight initialization
#' - Parameter storage and retrieval
#' - Distance/likelihood computation
#' - Time assignment
#' - Differential kinetics testing
#' 
#' @export
BaseDynamics <- R6::R6Class(
  "BaseDynamics",
  
  public = list(
    # Data
    gene = NULL,
    u = NULL,
    s = NULL,
    use_raw = NULL,
    
    # Parameters
    alpha = NULL,
    beta = NULL,
    gamma = NULL,
    scaling = NULL,
    t_ = NULL,
    alpha_ = 0,
    
    # Switching state
    u0_ = NULL,
    s0_ = NULL,
    
    # Time assignments
    t = NULL,
    tau = NULL,
    tau_ = NULL,
    o = NULL,
    
    # Weights
    weights = NULL,
    weights_upper = NULL,
    weights_outer = NULL,
    
    # Statistics
    std_u = NULL,
    std_s = NULL,
    u0 = 0,
    s0 = 0,
    
    # Results
    likelihood = NULL,
    loss = NULL,
    pars = NULL,
    varx = NULL,
    
    # Settings
    max_iter = 10,
    perc = 99,
    recoverable = TRUE,
    refit_time = TRUE,
    
    assignment_mode = NULL,
    steady_state_ratio = NULL,
    steady_state_prior = NULL,
    
    fit_scaling = TRUE,
    fit_steady_states = TRUE,
    fit_connected_states = TRUE,
    connectivities = NULL,
    high_pars_resolution = FALSE,
    init_vals = NULL,
    
    # Differential kinetics
    clusters = NULL,
    cats = NULL,
    orth_beta = NULL,
    diff_kinetics = NULL,
    pval_kinetics = NULL,
    pvals_kinetics = NULL,
    
    # Steady state info
    pval_steady = NULL,
    steady_u = NULL,
    steady_s = NULL,
    
    # Simplex optimizer settings
    simplex_kwargs = NULL,
    
    #' @description Initialize BaseDynamics
    #' @param u Unspliced expression
    #' @param s Spliced expression
    #' @param gene Gene name
    #' @param use_raw Use raw counts
    #' @param perc Percentile for filtering
    #' @param max_iter Maximum iterations
    #' @param fit_time Fit time
    #' @param fit_scaling Fit scaling
    #' @param fit_steady_states Fit steady states
    #' @param fit_connected_states Use connectivity smoothing
    #' @param fit_basal_transcription Fit basal transcription
    #' @param high_pars_resolution High resolution optimization
    #' @param steady_state_prior Prior for steady state
    #' @param init_vals Initial parameter values
    #' @param connectivities Connectivity matrix
    initialize = function(u, s, gene = NULL,
                          use_raw = FALSE,
                          perc = 99,
                          max_iter = 10,
                          fit_time = TRUE,
                          fit_scaling = TRUE,
                          fit_steady_states = TRUE,
                          fit_connected_states = TRUE,
                          fit_basal_transcription = NULL,
                          high_pars_resolution = FALSE,
                          steady_state_prior = NULL,
                          init_vals = NULL,
                          connectivities = NULL) {
      
      self$gene <- gene
      self$u <- as.vector(u)
      self$s <- as.vector(s)
      self$use_raw <- use_raw
      
      # Basal transcription
      if (!is.null(fit_basal_transcription) && fit_basal_transcription) {
        self$u0 <- min(self$u, na.rm = TRUE)
        self$s0 <- min(self$s, na.rm = TRUE)
        self$u <- self$u - self$u0
        self$s <- self$s - self$s0
      }
      
      self$max_iter <- max_iter
      self$simplex_kwargs <- list(
        method = "Nelder-Mead",
        control = list(maxit = as.integer(max_iter / 5))
      )
      
      self$perc <- perc
      self$recoverable <- TRUE
      
      # Initialize weights
      tryCatch({
        self$initialize_weights()
      }, error = function(e) {
        self$recoverable <- FALSE
        warning(paste("Model for", self$gene, "could not be instantiated:", e$message))
      })
      
      self$refit_time <- fit_time
      self$steady_state_prior <- steady_state_prior
      self$fit_scaling <- fit_scaling
      self$fit_steady_states <- fit_steady_states
      self$fit_connected_states <- fit_connected_states
      
      if (isTRUE(fit_connected_states)) {
        self$connectivities <- connectivities
      } else if (!isFALSE(fit_connected_states)) {
        self$connectivities <- fit_connected_states
      }
      
      self$high_pars_resolution <- high_pars_resolution
      self$init_vals <- init_vals
    },
    
    #' @description Initialize weights based on expression values
    #' @param weighted Use weighted filtering
    initialize_weights = function(weighted = TRUE) {
      nonzero_s <- self$s > 0
      nonzero_u <- self$u > 0
      
      weights <- nonzero_s & nonzero_u
      self$recoverable <- sum(weights) > 2
      
      if (self$recoverable) {
        if (weighted) {
          ub_s <- quantile(self$s[weights], self$perc / 100, na.rm = TRUE)
          ub_u <- quantile(self$u[weights], self$perc / 100, na.rm = TRUE)
          
          if (ub_s > 0) weights <- weights & (self$s <= ub_s)
          if (ub_u > 0) weights <- weights & (self$u <= ub_u)
        }
        
        self$weights <- weights
        u_w <- self$u[weights]
        s_w <- self$s[weights]
        
        self$std_u <- sd(u_w)
        self$std_s <- sd(s_w)
        
        # Upper weights (high expression)
        self$weights_upper <- weights
        if (any(weights)) {
          w_upper <- (self$u > max(u_w) / 3) & (self$s > max(s_w) / 3)
          self$weights_upper <- self$weights_upper & w_upper
        }
      }
    },
    
    #' @description Get current weights
    #' @param weighted Weight type: NULL, "outer", "upper"
    #' @param weights_cluster Cluster-specific weights
    get_weights = function(weighted = NULL, weights_cluster = NULL) {
      weights <- if (is.null(weighted)) {
        rep(TRUE, length(self$weights))
      } else if (weighted == "outer") {
        self$weights_outer
      } else if (weighted == "upper") {
        self$weights_upper
      } else {
        self$weights
      }
      
      if (!is.null(weights_cluster) && length(weights) == length(weights_cluster)) {
        weights <- weights & weights_cluster
      }
      
      weights
    },
    
    #' @description Get scaled reads
    #' @param scaling Scaling factor
    #' @param weighted Weight type
    #' @param weights_cluster Cluster weights
    get_reads = function(scaling = NULL, weighted = NULL, weights_cluster = NULL) {
      if (is.null(scaling)) scaling <- self$scaling
      
      u <- self$u / scaling
      s <- self$s
      
      if (!is.null(weighted) || !is.null(weights_cluster)) {
        weights <- self$get_weights(weighted, weights_cluster)
        u <- u[weights]
        s <- s[weights]
      }
      
      list(u = u, s = s)
    },
    
    #' @description Get current parameter values
    #' @param alpha Override alpha
    #' @param beta Override beta
    #' @param gamma Override gamma
    #' @param scaling Override scaling
    #' @param t_ Override t_
    #' @param u0_ Override u0_
    #' @param s0_ Override s0_
    get_vars = function(alpha = NULL, beta = NULL, gamma = NULL,
                        scaling = NULL, t_ = NULL, u0_ = NULL, s0_ = NULL) {
      
      if (is.null(alpha)) alpha <- self$alpha
      if (is.null(beta)) beta <- self$beta
      if (is.null(gamma)) gamma <- self$gamma
      if (is.null(scaling)) scaling <- self$scaling
      
      if (is.null(t_) || t_ == 0) {
        if (is.null(u0_)) {
          t_ <- self$t_
        } else {
          t_ <- tau_inv(u0_, s0_, 0, 0, alpha, beta, gamma)
        }
      }
      
      list(alpha = alpha, beta = beta, gamma = gamma, scaling = scaling, t_ = t_)
    },
    
    #' @description Compute divergence using current parameters
    #' @param alpha Override alpha
    #' @param beta Override beta
    #' @param gamma Override gamma
    #' @param scaling Override scaling
    #' @param t_ Override t_
    #' @param u0_ Override u0_
    #' @param s0_ Override s0_
    #' @param mode Divergence mode
    #' @param ... Additional arguments
    get_divergence = function(alpha = NULL, beta = NULL, gamma = NULL,
                              scaling = NULL, t_ = NULL, u0_ = NULL, s0_ = NULL,
                              mode = NULL, ...) {
      
      vars <- self$get_vars(alpha, beta, gamma, scaling, t_, u0_, s0_)
      
      reads <- self$get_reads(vars$scaling, weighted = FALSE)
      u <- reads$u
      s <- reads$s
      
      compute_divergence(
        u = u, s = s,
        alpha = vars$alpha, beta = vars$beta, gamma = vars$gamma,
        scaling = vars$scaling,
        t_ = vars$t_, u0_ = u0_, s0_ = s0_,
        std_u = self$std_u, std_s = self$std_s,
        mode = mode,
        assignment_mode = self$assignment_mode,
        connectivities = self$connectivities,
        fit_steady_states = self$fit_steady_states,
        ...
      )
    },
    
    #' @description Get time assignment
    #' @param alpha Override alpha
    #' @param beta Override beta
    #' @param gamma Override gamma
    #' @param scaling Override scaling
    #' @param t_ Override t_
    #' @param u0_ Override u0_
    #' @param s0_ Override s0_
    #' @param t Override time
    #' @param refit_time Whether to refit
    #' @param rescale_factor Rescaling factor
    #' @param weighted Weight type
    #' @param weights_cluster Cluster weights
    get_time_assignment = function(alpha = NULL, beta = NULL, gamma = NULL,
                                   scaling = NULL, t_ = NULL, u0_ = NULL, s0_ = NULL,
                                   t = NULL, refit_time = NULL, rescale_factor = NULL,
                                   weighted = NULL, weights_cluster = NULL) {
      
      if (is.null(refit_time)) refit_time <- self$refit_time
      
      if (!is.null(t)) {
        if (is.null(t_)) t_ <- self$t_
        o <- as.integer(t < t_)
        tau <- t * o + (t - t_) * (1 - o)
      } else if (refit_time) {
        if (!is.null(rescale_factor)) {
          vars <- self$get_vars(alpha, beta, gamma, scaling, t_, u0_, s0_)
          
          if (is.null(u0_)) u0_ <- self$u0_
          if (is.null(s0_)) s0_ <- self$s0_
          
          rescale_factor <- rescale_factor * vars$gamma / vars$beta
          vars$scaling <- vars$scaling * rescale_factor
          vars$beta <- vars$beta * rescale_factor
          u0_ <- u0_ / rescale_factor
          
          t_ <- tau_inv(u0_, s0_, 0, 0, vars$alpha, vars$beta, vars$gamma)
        }
        
        result <- self$get_divergence(alpha, beta, gamma, scaling, t_, u0_, s0_,
                                      mode = "assign_timepoints")
        t <- result$t
        tau <- result$tau
        o <- result$o
        
        if (!is.null(rescale_factor)) {
          t <- t * self$t_ / t_
          tau <- tau * self$t_ / t_
        }
      } else {
        t <- self$t
        tau <- self$tau
        o <- self$o
      }
      
      if (!is.null(weighted) || !is.null(weights_cluster)) {
        weights <- self$get_weights(weighted, weights_cluster)
        t <- t[weights]
        tau <- tau[weights]
        o <- o[weights]
      }
      
      list(t = t, tau = tau, o = o)
    },
    
    #' @description Get distances
    #' @param t Override time
    #' @param t_ Override t_
    #' @param alpha Override alpha
    #' @param beta Override beta
    #' @param gamma Override gamma
    #' @param scaling Override scaling
    #' @param u0_ Override u0_
    #' @param s0_ Override s0_
    #' @param refit_time Whether to refit
    #' @param weighted Weight type
    #' @param weights_cluster Cluster weights
    #' @param reg Regularization
    get_dists = function(t = NULL, t_ = NULL, alpha = NULL, beta = NULL, gamma = NULL,
                         scaling = NULL, u0_ = NULL, s0_ = NULL, refit_time = NULL,
                         weighted = TRUE, weights_cluster = NULL, reg = NULL) {
      
      reads <- self$get_reads(scaling, weighted = weighted, weights_cluster = weights_cluster)
      u <- reads$u
      s <- reads$s
      
      vars <- self$get_vars(alpha, beta, gamma, scaling, t_, u0_, s0_)
      time <- self$get_time_assignment(
        vars$alpha, vars$beta, vars$gamma, vars$scaling, vars$t_, u0_, s0_,
        t, refit_time, weighted = weighted, weights_cluster = weights_cluster
      )
      
      vparams <- vectorize(time$t, vars$t_, vars$alpha, vars$beta, vars$gamma)
      
      ut <- unspliced(vparams$tau, vparams$u0, vparams$alpha, vars$beta)
      st <- spliced(vparams$tau, vparams$s0, vparams$u0, vparams$alpha, vars$beta, vars$gamma)
      
      udiff <- (ut - u) / self$std_u * vars$scaling
      sdiff <- (st - s) / self$std_s
      
      if (is.null(reg)) {
        reg <- 0
        if (!is.null(self$steady_state_ratio)) {
          reg <- (vars$gamma / vars$beta - self$steady_state_ratio) * s / self$std_s
        }
      }
      
      list(udiff = udiff, sdiff = sdiff, reg = reg)
    },
    
    #' @description Get total distance
    #' @param noise_model Noise model ("normal" or "laplace")
    #' @param regularize Use regularization
    #' @param ... Additional arguments
    get_distx = function(noise_model = "normal", regularize = TRUE, ...) {
      dists <- self$get_dists(...)
      distx <- dists$udiff^2 + dists$sdiff^2
      
      if (regularize) {
        distx <- distx + dists$reg^2
      }
      
      if (noise_model == "laplace") {
        sqrt(distx)
      } else {
        distx
      }
    },
    
    #' @description Get sum of squared errors
    get_se = function(...) {
      sum(self$get_distx(...), na.rm = TRUE)
    },
    
    #' @description Get mean squared error
    get_mse = function(...) {
      mean(self$get_distx(...), na.rm = TRUE)
    },
    
    #' @description Get loss function value
    #' @param t Override time
    #' @param t_ Override t_
    #' @param alpha Override alpha
    #' @param beta Override beta
    #' @param gamma Override gamma
    #' @param scaling Override scaling
    #' @param u0_ Override u0_
    #' @param s0_ Override s0_
    #' @param refit_time Whether to refit
    get_loss = function(t = NULL, t_ = NULL, alpha = NULL, beta = NULL, gamma = NULL,
                        scaling = NULL, u0_ = NULL, s0_ = NULL, refit_time = NULL) {
      self$get_se(t = t, t_ = t_, alpha = alpha, beta = beta, gamma = gamma,
                  scaling = scaling, u0_ = u0_, s0_ = s0_, refit_time = refit_time)
    },
    
    #' @description Get log-likelihood
    #' @param varx Override variance
    #' @param noise_model Noise model
    #' @param ... Additional arguments
    get_loglikelihood = function(varx = NULL, noise_model = "normal", ...) {
      kwargs <- list(...)
      if (!"weighted" %in% names(kwargs)) {
        kwargs$weighted <- "upper"
      }
      
      dists <- do.call(self$get_dists, kwargs)
      distx <- dists$udiff^2 + dists$sdiff^2 + dists$reg^2
      eucl_distx <- sqrt(distx)
      n <- max(length(distx) - length(self$u) * 0.01, 2)
      
      # Compute variance
      if (is.null(varx)) {
        varx <- mean(distx, na.rm = TRUE) - mean(sign(dists$sdiff) * eucl_distx, na.rm = TRUE)^2
      }
      varx <- varx + (varx == 0)  # Edge case
      
      if (noise_model == "normal") {
        loglik <- -1 / 2 / n * sum(distx, na.rm = TRUE) / varx
        loglik <- loglik - 1 / 2 * log(2 * pi * varx)
      } else if (noise_model == "laplace") {
        loglik <- -1 / sqrt(2) / n * sum(eucl_distx, na.rm = TRUE) / sqrt(varx)
        loglik <- loglik - 1 / 2 * log(2 * varx)
      } else {
        stop("Unsupported noise model")
      }
      
      loglik
    },
    
    #' @description Get likelihood
    get_likelihood = function(...) {
      kwargs <- list(...)
      if (!"weighted" %in% names(kwargs)) {
        kwargs$weighted <- "upper"
      }
      exp(do.call(self$get_loglikelihood, kwargs))
    },
    
    #' @description Get variance
    get_variance = function(...) {
      kwargs <- list(...)
      if (!"weighted" %in% names(kwargs)) {
        kwargs$weighted <- "upper"
      }
      
      dists <- do.call(self$get_dists, kwargs)
      distx <- dists$udiff^2 + dists$sdiff^2
      mean(distx, na.rm = TRUE) - mean(sign(dists$sdiff) * sqrt(distx), na.rm = TRUE)^2
    },
    
    #' @description Get model unspliced at current times
    get_ut = function(...) {
      vars <- self$get_vars(...)
      time <- self$get_time_assignment(...)
      vparams <- vectorize(time$t, vars$t_, vars$alpha, vars$beta, vars$gamma)
      unspliced(vparams$tau, vparams$u0, vparams$alpha, vars$beta)
    },
    
    #' @description Get model spliced at current times
    get_st = function(...) {
      vars <- self$get_vars(...)
      time <- self$get_time_assignment(...)
      vparams <- vectorize(time$t, vars$t_, vars$alpha, vars$beta, vars$gamma)
      spliced(vparams$tau, vparams$s0, vparams$u0, vparams$alpha, vars$beta, vars$gamma)
    },
    
    #' @description Get velocity at current times
    #' @param mode Evaluation mode
    get_vt = function(mode = "soft_eval") {
      vars <- self$get_vars()
      result <- self$get_divergence(mode = mode)
      o_ <- result[1, ]
      o <- result[2, ]
      ut <- result[3, ]
      st <- result[4, ]
      ut * vars$beta - st * vars$gamma
    },
    
    # ========== Differential Kinetics Methods ==========
    
    #' @description Initialize for differential kinetics testing
    #' @param clusters Cluster assignments
    initialize_diff_kinetics = function(clusters) {
      if (is.null(self$varx)) {
        self$varx <- self$get_variance()
      }
      self$initialize_weights(weighted = FALSE)
      self$steady_state_ratio <- NULL
      self$clusters <- clusters
      self$cats <- unique(clusters)
      self$weights_outer <- self$weights & self$get_divergence(mode = "outside_of_trajectory")
    },
    
    #' @description Get orthogonal fit beta
    get_orth_fit = function(...) {
      kwargs <- list(...)
      kwargs$weighted <- TRUE
      
      reads <- do.call(self$get_reads, kwargs)
      u <- reads$u
      s <- reads$s
      
      a <- sum(s * u, na.rm = TRUE)
      b <- sum(u^2 - s^2, na.rm = TRUE)
      
      (b + sqrt(b^2 + 4 * a^2)) / (2 * a)
    },
    
    #' @description Get orthogonal distance
    #' @param orth_beta Orthogonal beta
    #' @param ... Additional arguments
    get_orth_distx = function(orth_beta = NULL, ...) {
      kwargs <- list(...)
      if (!"weighted" %in% names(kwargs)) {
        kwargs$weighted <- "outer"
      }
      
      reads <- do.call(self$get_reads, kwargs)
      u <- reads$u
      s <- reads$s
      
      if (is.null(orth_beta)) {
        orth_beta <- do.call(self$get_orth_fit, kwargs)
      }
      
      s_real <- (s + orth_beta * u) / (1 + orth_beta^2)
      sdiff <- (s_real - s) / self$std_s
      udiff <- (orth_beta * s_real - u) / self$std_u * self$scaling
      
      udiff^2 + sdiff^2
    },
    
    #' @description Get p-value
    #' @param model Model type ("dynamical" or "orthogonal")
    #' @param ... Additional arguments
    get_pval = function(model = "dynamical", ...) {
      kwargs <- list(...)
      if (!"weighted" %in% names(kwargs)) {
        kwargs$weighted <- "outer"
      }
      
      distx <- if (model == "orthogonal") {
        do.call(self$get_orth_distx, kwargs)
      } else {
        do.call(self$get_distx, kwargs) / 2
      }
      
      # Chi-squared test
      pchisq(sum(distx, na.rm = TRUE) / self$varx, df = 2 * length(distx), lower.tail = FALSE)
    },
    
    #' @description Get p-value for differential kinetics
    #' @param orth_beta Orthogonal beta
    #' @param min_cells Minimum cells required
    #' @param ... Additional arguments
    get_pval_diff_kinetics = function(orth_beta = NULL, min_cells = 10, ...) {
      kwargs <- list(...)
      
      if ("weights_cluster" %in% names(kwargs) && sum(kwargs$weights_cluster) < min_cells) {
        return(1)
      }
      
      if (!"weighted" %in% names(kwargs)) {
        kwargs$weighted <- "outer"
      }
      
      distx <- do.call(self$get_distx, kwargs) / 2
      orth_distx <- do.call(self$get_orth_distx, c(list(orth_beta = orth_beta), kwargs))
      
      denom <- self$varx * sqrt(4 * 2 * length(distx))
      
      # Normal approximation for large df
      pnorm((sum(distx, na.rm = TRUE) - sum(orth_distx, na.rm = TRUE)) / denom, lower.tail = FALSE)
    },
    
    #' @description Get cluster-specific MSE
    #' @param clusters Cluster assignments
    #' @param min_cells Minimum cells
    #' @param weighted Weight type
    get_cluster_mse = function(clusters = NULL, min_cells = 10, weighted = "outer") {
      if (is.null(self$clusters) || !is.null(clusters)) {
        self$initialize_diff_kinetics(clusters)
      }
      
      mse <- sapply(self$cats, function(c) {
        self$get_mse(weights_cluster = self$clusters == c, weighted = weighted)
      })
      
      if (!is.null(min_cells)) {
        w <- if (weighted == "outer") {
          self$weights_outer
        } else if (weighted == "upper") {
          self$weights_upper
        } else {
          self$weights
        }
        
        cell_counts <- sapply(self$cats, function(c) sum(w & (self$clusters == c)))
        mse[cell_counts < min_cells] <- 0
      }
      
      mse
    },
    
    #' @description Get cluster-specific p-values
    #' @param clusters Cluster assignments
    #' @param model Model type
    #' @param orth_beta Orthogonal beta
    #' @param ... Additional arguments
    get_cluster_pvals = function(clusters = NULL, model = NULL, orth_beta = NULL, ...) {
      if (is.null(self$clusters) || !is.null(clusters)) {
        self$initialize_diff_kinetics(clusters)
      }
      
      sapply(self$cats, function(c) {
        if (is.null(model)) {
          self$get_pval_diff_kinetics(weights_cluster = self$clusters == c,
                                      orth_beta = orth_beta, ...)
        } else {
          self$get_pval(model = model, weights_cluster = self$clusters == c, ...)
        }
      })
    },
    
    #' @description Run differential kinetic test
    #' @param clusters Cluster assignments
    #' @param as_df Return as data frame
    #' @param min_cells Minimum cells
    #' @param weighted Weight type
    differential_kinetic_test = function(clusters, as_df = FALSE, 
                                         min_cells = 10, weighted = "outer") {
      self$initialize_diff_kinetics(clusters)
      mse <- self$get_cluster_mse(weighted = weighted, min_cells = min_cells)
      
      weights_cluster <- self$clusters == self$cats[which.max(mse)]
      self$orth_beta <- self$get_orth_fit(weights_cluster = weights_cluster, weighted = FALSE)
      
      pval <- self$get_pval_diff_kinetics(weights_cluster = weights_cluster,
                                          orth_beta = self$orth_beta, weighted = weighted)
      
      if (pval > 1e-2) {
        weighted <- "upper"
        mse <- self$get_cluster_mse(weighted = weighted, min_cells = min_cells)
        self$orth_beta <- self$get_orth_fit(weights_cluster = self$clusters == self$cats[which.max(mse)])
      }
      
      self$pvals_kinetics <- self$get_cluster_pvals(orth_beta = self$orth_beta, weighted = weighted)
      
      self$diff_kinetics <- paste(self$cats[self$pvals_kinetics < 1e-2], collapse = ",")
      
      if (any(self$pvals_kinetics < 1e-2)) {
        self$pval_kinetics <- max(self$pvals_kinetics[self$pvals_kinetics < 1e-2])
      }
      
      if (as_df) {
        data.frame(
          cluster = self$cats,
          pval = round(self$pvals_kinetics, 4),
          stringsAsFactors = FALSE
        )
      } else {
        invisible(self)
      }
    }
  )
)

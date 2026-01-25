#' @title Dynamics Recovery Class
#' @description Complete implementation of DynamicsRecovery class
#' equivalent to scvelo.tools._em_model_core.DynamicsRecovery
#' This implements the EM algorithm for parameter inference.
#' @name DynamicsRecovery
NULL

#' Dynamics Recovery Class (R6)
#' 
#' @description Class that recovers transcriptional dynamics using the 
#' Expectation-Maximization (EM) algorithm. Extends BaseDynamics.
#' 
#' @details The EM algorithm iteratively:
#' 1. E-step: Assign time points to cells
#' 2. M-step: Optimize parameters (alpha, beta, gamma, t_, scaling)
#' 
#' @export
DynamicsRecovery <- R6::R6Class(
  "DynamicsRecovery",
  inherit = BaseDynamics,
  
  public = list(
    # Callbacks (for high resolution)
    cb_fit_t_and_alpha = NULL,
    cb_fit_scaling_ = NULL,
    cb_fit_rates = NULL,
    cb_fit_t_ = NULL,
    cb_fit_t_and_rates = NULL,
    cb_fit_rates_all = NULL,
    
    #' @description Initialize DynamicsRecovery and run initialization
    #' @param ... Arguments passed to BaseDynamics
    #' @param load_pars Load parameters from previous fit
    initialize = function(..., load_pars = NULL) {
      super$initialize(...)
      
      if (!is.null(load_pars) && load_pars) {
        # Load parameters logic would go here
        # For now, just initialize
      }
      
      if (self$recoverable) {
        self$init_parameters()
      }
    },
    
    #' @description Initialize parameters for EM
    init_parameters = function() {
      u <- self$u
      s <- self$s
      w <- self$weights
      perc <- 98
      
      u_w <- u[w]
      s_w <- s[w]
      
      # Initialize scaling
      if (self$std_u == 0 || self$std_s == 0) {
        self$std_u <- 1
        self$std_s <- 1
      }
      
      scaling <- if (isTRUE(self$fit_scaling)) {
        self$std_u / self$std_s
      } else if (is.numeric(self$fit_scaling)) {
        self$fit_scaling
      } else {
        1
      }
      
      u_scaled <- u / scaling
      u_w_scaled <- u_w / scaling
      
      # Initialize beta and gamma from extreme quantiles
      weights_s <- s_w >= quantile(s_w, perc / 100, na.rm = TRUE)
      weights_u <- u_w_scaled >= quantile(u_w_scaled, perc / 100, na.rm = TRUE)
      
      weights_g <- if (is.null(self$steady_state_prior)) {
        weights_s
      } else {
        weights_s | self$steady_state_prior[w]
      }
      
      beta <- 1
      gamma <- linreg(
        convolve_weights(u_w_scaled, weights_g),
        convolve_weights(s_w, weights_g)
      ) + 1e-6
      
      # Clip gamma
      if (gamma < 0.05 / scaling) gamma <- gamma * 1.2
      if (gamma > 1.5 / scaling) gamma <- gamma / 1.2
      
      # Initialize alpha
      selected <- weights_u | weights_s
      if (sum(selected, na.rm = TRUE) == 0) {
        selected <- rep(TRUE, length(u_w_scaled))
      }
      u_inf <- mean(u_w_scaled[selected], na.rm = TRUE)
      s_inf <- mean(s_w[weights_s], na.rm = TRUE)
      
      # Handle edge cases
      if (!is.finite(u_inf) || u_inf <= 0) u_inf <- max(u_w_scaled, na.rm = TRUE)
      if (!is.finite(s_inf) || s_inf <= 0) s_inf <- max(s_w, na.rm = TRUE)
      
      u0_ <- u_inf
      s0_ <- s_inf
      alpha <- u_inf * beta
      
      # Test bimodality
      pval_u <- 1
      pval_s <- 1
      means_u <- c(0, 0)
      means_s <- c(0, 0)
      
      tryCatch({
        bimod_u <- test_bimodality_kde(u_w)
        bimod_s <- test_bimodality_kde(s_w)
        pval_u <- bimod_u$pval
        pval_s <- bimod_s$pval
        means_u <- bimod_u$means
        means_s <- bimod_s$means
      }, error = function(e) {
        # Keep default values
      })
      
      self$pval_steady <- max(pval_u, pval_s)
      self$steady_u <- means_u[2]
      self$steady_s <- means_s[2]
      
      if (self$pval_steady < 1e-3) {
        u_inf <- mean(c(u_inf, self$steady_u))
        alpha <- gamma * s_inf
        beta <- alpha / u_inf
        u0_ <- u_inf
        s0_ <- s_inf
      }
      
      # Initialize switching time
      t_ <- tau_inv_u(u0_, 0, alpha, beta)
      if (!is.finite(t_) || t_ <= 0) t_ <- 1
      
      # Apply initial value modifications if provided
      if (!is.null(self$init_vals)) {
        init_vals <- as.numeric(self$init_vals)
        alpha <- alpha * init_vals[1]
        beta <- beta * init_vals[2]
        gamma <- gamma * init_vals[3]
      }
      
      # Store parameters
      self$alpha <- alpha
      self$beta <- beta
      self$gamma <- gamma
      self$scaling <- scaling
      self$alpha_ <- 0
      self$u0_ <- u0_
      self$s0_ <- s0_
      self$t_ <- t_
      
      self$pars <- matrix(c(alpha, beta, gamma, t_, scaling), ncol = 1)
      
      # Get initial time assignment
      time_result <- self$get_time_assignment()
      self$t <- time_result$t
      self$tau <- time_result$tau
      self$o <- time_result$o
      
      self$loss <- list(self$get_loss())
      
      # Initialize scaling search
      if (self$fit_scaling) {
        self$initialize_scaling(sight = 0.5)
        self$initialize_scaling(sight = 0.1)
      }
      
      self$steady_state_ratio <- self$gamma / self$beta
      
      # Set callbacks
      self$set_callbacks()
    },
    
    #' @description Initialize scaling by grid search
    #' @param sight Search range as fraction of current value
    initialize_scaling = function(sight = 0.5) {
      z_vals <- self$scaling + seq(-1, 1, length.out = 4) * self$scaling * sight
      
      for (z in z_vals) {
        self$update(scaling = z, beta = self$beta / self$scaling * z)
      }
    },
    
    #' @description Main fit function - runs EM algorithm
    #' @param assignment_mode Time assignment mode
    fit = function(assignment_mode = NULL) {
      if (self$max_iter > 0) {
        # For comparison with exact time assignment
        if (!is.null(assignment_mode) && assignment_mode == "full_projection") {
          self$assignment_mode <- assignment_mode
        }
        
        # Pre-train with explicit time assignment
        self$fit_t_and_alpha()
        
        if (self$fit_scaling) {
          self$fit_scaling_()
        }
        
        self$fit_rates()
        self$fit_t_()
        
        # Actual EM
        self$fit_t_and_rates()
        
        # Train with optimal time assignment
        self$assignment_mode <- assignment_mode
        self$update(adjust_t_ = FALSE)
        self$fit_t_and_rates(refit_time = FALSE)
      }
      
      # Final update
      self$update()
      
      # Get final tau values
      tau_result <- self$get_divergence(mode = "tau")
      self$tau <- tau_result$tau
      self$tau_ <- tau_result$tau_
      
      # Compute likelihood and variance
      self$likelihood <- self$get_likelihood(refit_time = FALSE)
      self$varx <- self$get_variance()
    },
    
    #' @description Fit t_ and alpha
    #' @param sight Search sight
    #' @param ... Additional arguments
    fit_t_and_alpha = function(sight = 0.5, ...) {
      # Grid search for alpha
      alpha_vals <- self$alpha + seq(-1, 1, length.out = 5) * self$alpha / 10
      for (alpha_val in alpha_vals) {
        self$update(alpha = alpha_val)
      }
      
      # Optimize t_ and alpha
      x0 <- c(self$t_, self$alpha)
      
      result <- optim(
        par = x0,
        fn = function(x) {
          self$get_mse(t_ = x[1], alpha = x[2], ...)
        },
        method = "Nelder-Mead",
        control = list(maxit = self$max_iter %/% 5)
      )
      
      self$update(t_ = result$par[1], alpha = result$par[2])
    },
    
    #' @description Fit alpha only
    #' @param sight Search sight
    #' @param ... Additional arguments
    fit_alpha = function(sight = 0.5, ...) {
      val <- self$alpha
      vals <- val + seq(-1, 1, length.out = 4) * val * sight
      
      for (v in vals) {
        self$update(alpha = v)
      }
      
      result <- optim(
        par = val,
        fn = function(x) self$get_mse(alpha = x[1], ...),
        method = "Nelder-Mead",
        control = list(maxit = self$max_iter %/% 5)
      )
      
      self$update(alpha = result$par[1])
    },
    
    #' @description Fit beta only
    #' @param sight Search sight
    #' @param ... Additional arguments
    fit_beta = function(sight = 0.5, ...) {
      val <- self$beta
      vals <- val + seq(-1, 1, length.out = 4) * val * sight
      
      for (v in vals) {
        self$update(beta = v)
      }
      
      result <- optim(
        par = val,
        fn = function(x) self$get_mse(beta = x[1], ...),
        method = "Nelder-Mead",
        control = list(maxit = self$max_iter %/% 5)
      )
      
      self$update(beta = result$par[1])
    },
    
    #' @description Fit gamma only
    #' @param sight Search sight
    #' @param ... Additional arguments
    fit_gamma = function(sight = 0.5, ...) {
      val <- self$gamma
      vals <- val + seq(-1, 1, length.out = 4) * val * sight
      
      for (v in vals) {
        self$update(gamma = v)
      }
      
      result <- optim(
        par = val,
        fn = function(x) self$get_mse(gamma = x[1], ...),
        method = "Nelder-Mead",
        control = list(maxit = self$max_iter %/% 5)
      )
      
      self$update(gamma = result$par[1])
    },
    
    #' @description Fit rates (alpha and gamma)
    #' @param ... Additional arguments
    fit_rates = function(...) {
      x0 <- c(self$alpha, self$gamma)
      
      result <- optim(
        par = x0,
        fn = function(x) {
          self$get_mse(alpha = x[1], gamma = x[2], ...)
        },
        method = "Nelder-Mead",
        control = list(maxit = self$max_iter %/% 5, reltol = 1e-2)
      )
      
      self$update(alpha = result$par[1], gamma = result$par[2])
    },
    
    #' @description Fit t_ only
    #' @param ... Additional arguments
    fit_t_ = function(...) {
      result <- optim(
        par = self$t_,
        fn = function(x) self$get_mse(t_ = x[1], ...),
        method = "Nelder-Mead",
        control = list(maxit = self$max_iter %/% 5)
      )
      
      self$update(t_ = result$par[1])
    },
    
    #' @description Fit all rates (alpha, beta, gamma)
    #' @param ... Additional arguments
    fit_rates_all = function(...) {
      x0 <- c(self$alpha, self$beta, self$gamma)
      
      result <- optim(
        par = x0,
        fn = function(x) {
          self$get_mse(alpha = x[1], beta = x[2], gamma = x[3], ...)
        },
        method = "Nelder-Mead",
        control = list(maxit = self$max_iter %/% 5, reltol = 1e-2)
      )
      
      self$update(alpha = result$par[1], beta = result$par[2], gamma = result$par[3])
    },
    
    #' @description Fit t_ and all rates
    #' @param ... Additional arguments
    fit_t_and_rates = function(...) {
      x0 <- c(self$t_, self$alpha, self$beta, self$gamma)
      
      result <- optim(
        par = x0,
        fn = function(x) {
          self$get_mse(t_ = x[1], alpha = x[2], beta = x[3], gamma = x[4], ...)
        },
        method = "Nelder-Mead",
        control = list(maxit = self$max_iter %/% 5, reltol = 1e-2)
      )
      
      self$update(t_ = result$par[1], alpha = result$par[2],
                  beta = result$par[3], gamma = result$par[4])
    },
    
    #' @description Fit scaling
    #' @param ... Additional arguments
    fit_scaling_ = function(...) {
      x0 <- c(self$t_, self$beta, self$scaling)
      
      result <- optim(
        par = x0,
        fn = function(x) {
          self$get_mse(t_ = x[1], beta = x[2], scaling = x[3], ...)
        },
        method = "Nelder-Mead",
        control = list(maxit = self$max_iter %/% 5)
      )
      
      self$update(t_ = result$par[1], beta = result$par[2], scaling = result$par[3])
    },
    
    #' @description Set callbacks for high-resolution optimization
    set_callbacks = function() {
      if (!self$high_pars_resolution) {
        self$cb_fit_t_and_alpha <- NULL
        self$cb_fit_scaling_ <- NULL
        self$cb_fit_rates <- NULL
        self$cb_fit_t_ <- NULL
        self$cb_fit_t_and_rates <- NULL
        self$cb_fit_rates_all <- NULL
      }
    },
    
    #' @description Update parameters if improvement achieved
    #' @param t Override time
    #' @param t_ Override switching time
    #' @param alpha Override alpha
    #' @param beta Override beta
    #' @param gamma Override gamma
    #' @param scaling Override scaling
    #' @param u0_ Override u0_
    #' @param s0_ Override s0_
    #' @param adjust_t_ Whether to adjust t_
    update = function(t = NULL, t_ = NULL, alpha = NULL, beta = NULL, gamma = NULL,
                      scaling = NULL, u0_ = NULL, s0_ = NULL, adjust_t_ = TRUE) {
      
      loss_prev <- if (length(self$loss) > 0) self$loss[[length(self$loss)]] else 1e6
      
      # Get current or new values
      vars <- self$get_vars(alpha, beta, gamma, scaling, t_, u0_, s0_)
      alpha <- vars$alpha
      beta <- vars$beta
      gamma <- vars$gamma
      scaling <- vars$scaling
      t_ <- vars$t_
      
      # Get time assignment
      time_result <- self$get_time_assignment(alpha, beta, gamma, scaling, t_, u0_, s0_, t)
      t_new <- time_result$t
      tau <- time_result$tau
      o <- time_result$o
      
      # Compute loss
      loss <- self$get_loss(t_new, t_, alpha, beta, gamma, scaling)
      
      perform_update <- loss < loss_prev
      
      # Adjust t_ if needed
      on <- o == 1
      if (adjust_t_ && any(on)) {
        if (!perform_update) {
          vars <- self$get_vars()
          alpha <- vars$alpha
          beta <- vars$beta
          gamma <- vars$gamma
          scaling <- vars$scaling
          t_ <- vars$t_
          
          time_result <- self$get_time_assignment()
          t_new <- time_result$t
          tau <- time_result$tau
          o <- time_result$o
          
          loss <- self$get_loss()
        }
        
        alt_t_ <- max(t_new[on], na.rm = TRUE)
        
        if (alt_t_ > 0 && alt_t_ < t_) {
          alt_t_ <- alt_t_ + max(t_new, na.rm = TRUE) / length(t_new) * sum(t_new == t_, na.rm = TRUE)
          
          alt_time <- self$get_time_assignment(alpha, beta, gamma, scaling, alt_t_)
          alt_t_new <- alt_time$t
          alt_tau <- alt_time$tau
          alt_o <- alt_time$o
          
          alt_loss <- self$get_loss(alt_t_new, alt_t_, alpha, beta, gamma, scaling)
          
          ut_cur <- unspliced(t_, 0, alpha, beta)
          ut_alt <- unspliced(alt_t_, 0, alpha, beta)
          
          min_loss <- min(loss, loss_prev)
          
          if (alt_loss * 0.99 <= min_loss || ut_cur * 0.99 < ut_alt) {
            t_new <- alt_t_new
            tau <- alt_tau
            o <- alt_o
            t_ <- alt_t_
            loss <- alt_loss
            perform_update <- TRUE
          }
        }
      }
      
      # Perform update if improved
      if (perform_update) {
        if (!is.null(scaling) && scaling != self$scaling) {
          self$steady_u <- self$steady_u * self$scaling / scaling
          self$u0_ <- self$u0_ * self$scaling / scaling
        }
        
        if (!is.null(u0_)) self$u0_ <- u0_
        if (!is.null(s0_)) self$s0_ <- s0_
        
        self$t <- t_new
        self$tau <- tau
        self$o <- o
        self$alpha <- alpha
        self$beta <- beta
        self$gamma <- gamma
        self$scaling <- scaling
        self$t_ <- t_
        
        new_pars <- c(alpha, beta, gamma, t_, scaling)
        self$pars <- cbind(self$pars, new_pars)
        self$loss <- c(self$loss, list(loss))
      }
      
      perform_update
    }
  )
)

#' Test bimodality using KDE
#' 
#' @description Test whether expression is bimodal using kernel density estimation
#' @param x Expression values
#' @param kde Use KDE (TRUE) or histogram (FALSE)
#' @return List with pval, means, and bimodal flag
#' @keywords internal
test_bimodality_kde <- function(x, kde = TRUE) {
  x <- x[x > 0 & is.finite(x)]
  
  if (length(x) < 20) {
    return(list(pval = 1, means = c(0, 0), bimodal = FALSE))
  }
  
  if (kde) {
    # KDE-based test
    dens <- density(x, n = 512)
    
    # Find local maxima
    y <- dens$y
    d <- diff(y)
    local_max_idx <- which(d[-length(d)] > 0 & d[-1] < 0) + 1
    
    if (length(local_max_idx) >= 2) {
      # Get the two highest peaks
      peak_heights <- y[local_max_idx]
      top2 <- order(peak_heights, decreasing = TRUE)[1:2]
      
      means <- sort(dens$x[local_max_idx[top2]])
      
      # Simple bimodality test: is there a significant dip between peaks?
      idx1 <- which.min(abs(dens$x - means[1]))
      idx2 <- which.min(abs(dens$x - means[2]))
      
      if (idx1 > idx2) {
        tmp <- idx1
        idx1 <- idx2
        idx2 <- tmp
      }
      
      min_between <- min(y[idx1:idx2])
      max_peaks <- min(y[local_max_idx[top2]])
      
      # Bimodal if dip is significant
      dip_ratio <- min_between / max_peaks
      
      if (dip_ratio < 0.5) {
        return(list(pval = 1e-3, means = means, bimodal = TRUE))
      }
    }
  }
  
  list(pval = 1, means = c(0, max(x, na.rm = TRUE)), bimodal = FALSE)
}

#' Recover dynamics for multiple genes
#' 
#' @description Main function to recover transcriptional dynamics for multiple genes
#' @param Ms Spliced moment matrix (cells x genes)
#' @param Mu Unspliced moment matrix (cells x genes)
#' @param gene_names Gene names
#' @param var_names Genes to fit
#' @param n_top_genes Number of top genes
#' @param max_iter Maximum EM iterations
#' @param fit_scaling Fit scaling
#' @param fit_time Fit time
#' @param fit_steady_states Fit steady states
#' @param fit_connected_states Use connectivity smoothing
#' @param connectivities Connectivity matrix
#' @param assignment_mode Time assignment mode
#' @param t_max Total time range
#' @param n_cores Number of cores
#' @param verbose Print progress
#' @return List with gene parameters and time matrices
#' @export
recover_dynamics_multi <- function(Ms, Mu, gene_names,
                                   var_names = NULL,
                                   n_top_genes = NULL,
                                   max_iter = 10,
                                   fit_scaling = TRUE,
                                   fit_time = TRUE,
                                   fit_steady_states = TRUE,
                                   fit_connected_states = NULL,
                                   connectivities = NULL,
                                   assignment_mode = "projection",
                                   t_max = 20,
                                   n_cores = 1,
                                   verbose = TRUE) {
  
  # Determine genes to fit
  if (is.null(var_names)) {
    var_names <- gene_names
  } else if (length(var_names) == 1 && var_names == "all") {
    var_names <- gene_names
  }
  
  # Filter to valid genes
  var_names <- intersect(var_names, gene_names)
  
  # Limit genes if requested
  if (!is.null(n_top_genes) && length(var_names) > n_top_genes) {
    gene_sums <- colSums(Ms[, var_names, drop = FALSE], na.rm = TRUE)
    var_names <- var_names[order(gene_sums, decreasing = TRUE)[1:n_top_genes]]
  }
  
  n_genes <- length(var_names)
  n_cells <- nrow(Ms)
  
  if (verbose) {
    message(sprintf("Recovering dynamics for %d genes...", n_genes))
  }
  
  # Use connected states
  if (is.null(fit_connected_states)) {
    fit_connected_states <- !is.null(connectivities)
  }
  
  # Initialize result containers
  gene_indices <- match(var_names, gene_names)
  
  alpha <- rep(NA_real_, n_genes)
  beta <- rep(NA_real_, n_genes)
  gamma <- rep(NA_real_, n_genes)
  t_ <- rep(NA_real_, n_genes)
  scaling <- rep(NA_real_, n_genes)
  std_u <- rep(NA_real_, n_genes)
  std_s <- rep(NA_real_, n_genes)
  likelihood <- rep(NA_real_, n_genes)
  u0_vec <- rep(NA_real_, n_genes)
  s0_vec <- rep(NA_real_, n_genes)
  pval_steady <- rep(NA_real_, n_genes)
  steady_u <- rep(NA_real_, n_genes)
  steady_s <- rep(NA_real_, n_genes)
  variance <- rep(NA_real_, n_genes)
  
  T_matrix <- matrix(NA_real_, n_cells, n_genes)
  Tau_matrix <- matrix(NA_real_, n_cells, n_genes)
  Tau_matrix_ <- matrix(NA_real_, n_cells, n_genes)
  
  # Fit function for single gene
  fit_gene <- function(i) {
    idx <- gene_indices[i]
    gene <- var_names[i]
    
    u <- Mu[, idx]
    s <- Ms[, idx]
    
    tryCatch({
      dm <- DynamicsRecovery$new(
        u = u, s = s,
        gene = gene,
        max_iter = max_iter,
        fit_scaling = fit_scaling,
        fit_time = fit_time,
        fit_steady_states = fit_steady_states,
        fit_connected_states = fit_connected_states,
        connectivities = connectivities
      )
      
      if (dm$recoverable) {
        dm$fit(assignment_mode = assignment_mode)
        
        list(
          success = TRUE,
          alpha = dm$alpha,
          beta = dm$beta / dm$scaling,  # Unscale
          gamma = dm$gamma,
          t_ = dm$t_,
          scaling = dm$scaling,
          std_u = dm$std_u,
          std_s = dm$std_s,
          likelihood = dm$likelihood,
          variance = dm$varx,
          u0 = dm$u0,
          s0 = dm$s0,
          pval_steady = dm$pval_steady,
          steady_u = dm$steady_u * dm$scaling,
          steady_s = dm$steady_s,
          t = dm$t,
          tau = dm$tau,
          tau_ = dm$tau_
        )
      } else {
        list(success = FALSE)
      }
    }, error = function(e) {
      list(success = FALSE, error = e$message)
    })
  }
  
  # Run fitting
  if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    results <- parallel::mclapply(seq_len(n_genes), fit_gene, mc.cores = n_cores)
  } else {
    results <- lapply(seq_len(n_genes), function(i) {
      if (verbose && i %% 100 == 0) {
        message(sprintf("  Fitted %d/%d genes", i, n_genes))
      }
      fit_gene(i)
    })
  }
  
  # Collect results
  n_success <- 0
  for (i in seq_len(n_genes)) {
    res <- results[[i]]
    if (res$success) {
      n_success <- n_success + 1
      alpha[i] <- res$alpha
      beta[i] <- res$beta
      gamma[i] <- res$gamma
      t_[i] <- res$t_
      scaling[i] <- res$scaling
      std_u[i] <- res$std_u
      std_s[i] <- res$std_s
      likelihood[i] <- res$likelihood
      variance[i] <- res$variance
      u0_vec[i] <- res$u0
      s0_vec[i] <- res$s0
      pval_steady[i] <- res$pval_steady
      steady_u[i] <- res$steady_u
      steady_s[i] <- res$steady_s
      T_matrix[, i] <- res$t
      Tau_matrix[, i] <- res$tau
      if (!is.null(res$tau_)) {
        Tau_matrix_[, i] <- res$tau_
      }
    }
  }
  
  if (verbose) {
    message(sprintf("  Successfully fitted %d/%d genes", n_success, n_genes))
  }
  
  # Align dynamics
  if (!is.null(t_max) && t_max != FALSE) {
    align_result <- align_dynamics_multi(
      T_matrix, t_, alpha, beta, gamma,
      t_max = t_max
    )
    
    T_matrix <- align_result$T
    t_ <- align_result$t_
    alpha <- align_result$alpha
    beta <- align_result$beta
    gamma <- align_result$gamma
    Tau_matrix <- align_result$Tau
    Tau_matrix_ <- align_result$Tau_
    alignment_scaling <- align_result$alignment_scaling
  } else {
    alignment_scaling <- rep(1, n_genes)
  }
  
  # Smooth time with connectivities
  if (!is.null(connectivities)) {
    for (i in which(!is.na(t_))) {
      T_matrix[, i] <- as.vector(connectivities %*% T_matrix[, i])
    }
  }
  
  # Create result
  list(
    gene_params = data.frame(
      gene = var_names,
      alpha = alpha,
      beta = beta,
      gamma = gamma,
      t_ = t_,
      scaling = scaling,
      std_u = std_u,
      std_s = std_s,
      likelihood = likelihood,
      u0 = u0_vec,
      s0 = s0_vec,
      pval_steady = pval_steady,
      steady_u = steady_u,
      steady_s = steady_s,
      variance = variance,
      alignment_scaling = alignment_scaling,
      stringsAsFactors = FALSE
    ),
    fit_t = T_matrix,
    fit_tau = Tau_matrix,
    fit_tau_ = Tau_matrix_,
    genes = var_names,
    gene_indices = gene_indices
  )
}

#' Align dynamics to common time scale
#' 
#' @param T_matrix Time matrix
#' @param t_ Switching times
#' @param alpha Transcription rates
#' @param beta Splicing rates
#' @param gamma Degradation rates
#' @param t_max Target max time
#' @return List with aligned values
#' @keywords internal
align_dynamics_multi <- function(T_matrix, t_, alpha, beta, gamma, t_max = 20) {
  n_genes <- ncol(T_matrix)
  n_cells <- nrow(T_matrix)
  
  # Copy for modification
  Tau_matrix <- matrix(NA_real_, n_cells, n_genes)
  Tau_matrix_ <- matrix(NA_real_, n_cells, n_genes)
  
  # Compute T_max for each gene
  T_max <- numeric(n_genes)
  
  for (i in seq_len(n_genes)) {
    if (is.na(t_[i])) next
    
    t_col <- T_matrix[, i]
    valid <- !is.na(t_col)
    
    if (!any(valid)) next
    
    # Time in on phase
    on_mask <- t_col < t_[i]
    T_on <- if (any(on_mask & valid)) max(t_col[on_mask & valid], na.rm = TRUE) else 0
    
    # Time in off phase
    off_mask <- t_col >= t_[i]
    T_off <- if (any(off_mask & valid)) max(t_col[off_mask & valid] - t_[i], na.rm = TRUE) else 0
    
    T_max[i] <- T_on + T_off
    
    # Adjust for steady state cells
    n_steady <- sum(t_col == t_[i] | t_col == 0, na.rm = TRUE)
    n_valid <- sum(valid)
    
    if (n_steady > 0 && n_valid > n_steady && T_max[i] > 0) {
      T_max[i] <- T_max[i] / (1 - n_steady / n_valid)
    }
    
    # Compute tau values
    o <- as.integer(t_col < t_[i])
    Tau_matrix[, i] <- t_col * o + (t_col - t_[i]) * (1 - o)
    Tau_matrix_[, i] <- (t_col - t_[i]) * (1 - o)
  }
  
  T_max[T_max == 0 | !is.finite(T_max)] <- 1
  
  # Compute alignment scaling
  m <- t_max / T_max
  m[!is.finite(m)] <- 1
  
  # Apply scaling
  for (i in seq_len(n_genes)) {
    if (is.na(t_[i]) || m[i] == 1) next
    
    T_matrix[, i] <- T_matrix[, i] * m[i]
    Tau_matrix[, i] <- Tau_matrix[, i] * m[i]
    Tau_matrix_[, i] <- Tau_matrix_[, i] * m[i]
    t_[i] <- t_[i] * m[i]
    alpha[i] <- alpha[i] / m[i]
    beta[i] <- beta[i] / m[i]
    gamma[i] <- gamma[i] / m[i]
  }
  
  alignment_scaling <- m
  alignment_scaling[alignment_scaling == 1] <- NA
  
  list(
    T = T_matrix,
    Tau = Tau_matrix,
    Tau_ = Tau_matrix_,
    t_ = t_,
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    alignment_scaling = alignment_scaling
  )
}

#' @title Dynamical Velocity Model
#' @description Recover full splicing kinetics using the EM algorithm.
#' The dynamical model infers transcription rates (alpha), splicing rates (beta),
#' degradation rates (gamma), as well as cell-specific latent time and
#' transcriptional states.
#' @name velocity_dynamical
NULL

#' Recover Dynamics
#'
#' @description Main function to recover transcriptional dynamics using the 
#' dynamical model with EM optimization.
#'
#' @param object Seurat object
#' @param genes Character vector of gene names, or "velocity_genes"
#' @param n_top_genes Maximum number of top genes to use
#' @param max_iter Maximum EM iterations
#' @param fit_scaling Whether to fit scaling between unspliced and spliced
#' @param fit_time Whether to fit cell times
#' @param t_max Total time range for alignment (default 20)
#' @param min_likelihood Minimum likelihood threshold for genes
#' @param n_cores Number of cores for parallel processing
#' @param verbose Print progress
#' @param spliced_layer Name of spliced layer
#' @param unspliced_layer Name of unspliced layer
#' @param use_moments Whether to use moment-based estimates (Ms, Mu)
#'
#' @return Seurat object with dynamics parameters in misc$scVeloR$dynamics
#' @export
recover_dynamics <- function(object, 
                             genes = "velocity_genes",
                             n_top_genes = NULL,
                             max_iter = 10,
                             fit_scaling = TRUE,
                             fit_time = TRUE,
                             t_max = 20,
                             min_likelihood = 0.1,
                             n_cores = 1,
                             verbose = TRUE,
                             spliced_layer = "spliced",
                             unspliced_layer = "unspliced",
                             use_moments = TRUE) {
  
  # Get data matrices
  data <- get_velocity_data(object, spliced_layer, unspliced_layer, use_moments)
  Ms <- data$Ms
  Mu <- data$Mu
  
  # Determine genes to fit
  gene_names <- colnames(Ms)
  
  if (is.character(genes) && length(genes) == 1) {
    if (genes == "velocity_genes" || genes == "all") {
      # Use velocity genes if available
      if (!is.null(object@misc$scVeloR$velocity$velocity_genes)) {
        velocity_genes <- object@misc$scVeloR$velocity$velocity_genes
        genes <- gene_names[velocity_genes]
      } else {
        genes <- gene_names
      }
    } else if (genes %in% names(object@misc$scVeloR)) {
      genes <- object@misc$scVeloR[[genes]]
    }
  }
  
  # Filter to valid genes
  genes <- intersect(genes, gene_names)
  
  # Limit number of genes
  if (!is.null(n_top_genes) && length(genes) > n_top_genes) {
    gene_totals <- colSums(Ms[, genes, drop = FALSE])
    genes <- genes[order(gene_totals, decreasing = TRUE)[1:n_top_genes]]
  }
  
  n_genes <- length(genes)
  n_cells <- nrow(Ms)
  
  if (verbose) {
    message(sprintf("Recovering dynamics for %d genes using %d cores...", n_genes, n_cores))
  }
  
  # Get connectivities for connected state fitting
  conn <- NULL
  if (!is.null(object@misc$scVeloR$neighbors$connectivities)) {
    conn <- object@misc$scVeloR$neighbors$connectivities
  }
  
  # Initialize result containers
  alpha <- rep(NA_real_, n_genes)
  beta <- rep(NA_real_, n_genes)
  gamma <- rep(NA_real_, n_genes)
  t_ <- rep(NA_real_, n_genes)
  scaling <- rep(NA_real_, n_genes)
  std_u <- rep(NA_real_, n_genes)
  std_s <- rep(NA_real_, n_genes)
  likelihood <- rep(NA_real_, n_genes)
  variance <- rep(NA_real_, n_genes)
  
  T_matrix <- matrix(NA_real_, n_cells, n_genes)
  Tau_matrix <- matrix(NA_real_, n_cells, n_genes)
  O_matrix <- matrix(NA_integer_, n_cells, n_genes)
  
  # Fit dynamics for each gene
  gene_indices <- match(genes, gene_names)
  
  fit_gene <- function(i) {
    idx <- gene_indices[i]
    gene <- genes[i]
    
    u <- Mu[, idx]
    s <- Ms[, idx]
    
    tryCatch({
      dm <- dynamics_recovery(
        u = u, s = s,
        gene = gene,
        max_iter = max_iter,
        fit_scaling = fit_scaling,
        fit_time = fit_time,
        conn = conn
      )
      
      list(
        success = dm$recoverable,
        alpha = dm$alpha,
        beta = dm$beta,
        gamma = dm$gamma,
        t_ = dm$t_,
        scaling = dm$scaling,
        std_u = dm$std_u,
        std_s = dm$std_s,
        likelihood = dm$likelihood,
        variance = dm$variance,
        t = dm$t,
        tau = dm$tau,
        o = dm$o
      )
    }, error = function(e) {
      list(success = FALSE)
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
      beta[i] <- res$beta / res$scaling  # Unscale beta
      gamma[i] <- res$gamma
      t_[i] <- res$t_
      scaling[i] <- res$scaling
      std_u[i] <- res$std_u
      std_s[i] <- res$std_s
      likelihood[i] <- res$likelihood
      variance[i] <- res$variance
      T_matrix[, i] <- res$t
      Tau_matrix[, i] <- res$tau
      O_matrix[, i] <- res$o
    }
  }
  
  if (verbose) {
    message(sprintf("  Successfully fitted %d/%d genes", n_success, n_genes))
  }
  
  # Align dynamics to common time scale
  if (!is.null(t_max) && t_max != FALSE) {
    alignment_result <- align_dynamics(
      T_matrix, t_, alpha, beta, gamma, scaling,
      t_max = t_max
    )
    
    T_matrix <- alignment_result$T
    t_ <- alignment_result$t_
    alpha <- alignment_result$alpha
    beta <- alignment_result$beta
    gamma <- alignment_result$gamma
    alignment_scaling <- alignment_result$alignment_scaling
  } else {
    alignment_scaling <- rep(1, n_genes)
  }
  
  # Smooth time with connectivities
  if (!is.null(conn)) {
    for (i in which(!is.na(t_))) {
      T_matrix[, i] <- as.vector(conn %*% T_matrix[, i])
    }
  }
  
  # Create gene parameter data frame
  gene_params <- data.frame(
    gene = genes,
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    t_ = t_,
    scaling = scaling,
    std_u = std_u,
    std_s = std_s,
    likelihood = likelihood,
    variance = variance,
    alignment_scaling = alignment_scaling,
    stringsAsFactors = FALSE
  )
  rownames(gene_params) <- genes
  
  # Store results
  if (is.null(object@misc$scVeloR)) {
    object@misc$scVeloR <- list()
  }
  
  object@misc$scVeloR$dynamics <- list(
    gene_params = gene_params,
    fit_t = T_matrix,
    fit_tau = Tau_matrix,
    fit_o = O_matrix,
    genes = genes,
    gene_indices = gene_indices,
    t_max = t_max
  )
  
  if (verbose) {
    message("Done. Results stored in object@misc$scVeloR$dynamics")
  }
  
  object
}

#' Single Gene Dynamics Recovery
#'
#' @description Fit dynamics model for a single gene using EM algorithm.
#'
#' @param u Unspliced values
#' @param s Spliced values
#' @param gene Gene name
#' @param max_iter Maximum iterations
#' @param fit_scaling Whether to fit scaling
#' @param fit_time Whether to fit time
#' @param conn Connectivity matrix
#' @param perc Percentile for weight initialization
#'
#' @return List with fitted parameters
#' @keywords internal
dynamics_recovery <- function(u, s, gene = NULL,
                              max_iter = 10,
                              fit_scaling = TRUE,
                              fit_time = TRUE,
                              conn = NULL,
                              perc = 99) {
  
  # Initialize output structure
  result <- list(
    gene = gene,
    recoverable = FALSE,
    alpha = NA, beta = NA, gamma = NA,
    t_ = NA, scaling = NA,
    std_u = NA, std_s = NA,
    likelihood = NA, variance = NA,
    t = NA, tau = NA, o = NA
  )
  
  # Convert to vectors
  u <- as.vector(u)
  s <- as.vector(s)
  n <- length(u)
  
  # Initialize weights (filter extreme values and zeros)
  nonzero <- (s > 0) & (u > 0)
  
  if (sum(nonzero) < 10) {
    return(result)
  }
  
  weights <- nonzero
  
  # Apply percentile filter
  ub_s <- quantile(s[nonzero], perc / 100, na.rm = TRUE)
  ub_u <- quantile(u[nonzero], perc / 100, na.rm = TRUE)
  
  if (ub_s > 0) weights <- weights & (s <= ub_s)
  if (ub_u > 0) weights <- weights & (u <= ub_u)
  
  if (sum(weights) < 10) {
    return(result)
  }
  
  u_w <- u[weights]
  s_w <- s[weights]
  
  # Initialize std
  std_u <- sd(u_w)
  std_s <- sd(s_w)
  
  if (std_u == 0 || std_s == 0 || !is.finite(std_u) || !is.finite(std_s)) {
    return(result)
  }
  
  # Initialize scaling
  if (isTRUE(fit_scaling)) {
    scaling <- std_u / std_s
  } else if (is.numeric(fit_scaling)) {
    scaling <- fit_scaling
  } else {
    scaling <- 1
  }
  
  # Scale u
  u_scaled <- u / scaling
  u_w_scaled <- u_w / scaling
  
  # Initialize beta and gamma from extreme quantiles
  perc_high <- 98
  weights_s <- s_w >= quantile(s_w, perc_high / 100)
  weights_u <- u_w_scaled >= quantile(u_w_scaled, perc_high / 100)
  
  # Initialize gamma from steady-state regression
  beta <- 1
  
  # Linear regression: gamma = u_w / s_w at steady state
  ss_mask <- weights_s
  if (sum(ss_mask) < 3) ss_mask <- rep(TRUE, length(u_w_scaled))
  
  u_ss <- mean(u_w_scaled[ss_mask])
  s_ss <- mean(s_w[ss_mask])
  
  gamma <- u_ss / s_ss
  if (!is.finite(gamma) || gamma <= 0) gamma <- 1
  
  # Clip gamma to reasonable range
  if (gamma < 0.05 / scaling) gamma <- gamma * 1.2
  if (gamma > 1.5 / scaling) gamma <- gamma / 1.2
  
  # Initialize alpha and switching time
  u_inf <- mean(u_w_scaled[weights_u | weights_s])
  s_inf <- mean(s_w[weights_s])
  
  alpha <- u_inf * beta
  
  # Initial estimate of switching time
  t_ <- tau_inv_u(u_inf, 0, alpha, beta)
  if (!is.finite(t_) || t_ <= 0) t_ <- 1
  
  # Compute initial state at switching
  result_switch <- mRNA(t_, 0, 0, alpha, beta, gamma)
  u0_ <- result_switch$u
  s0_ <- result_switch$s
  
  # Get initial time assignment
  assignment <- assign_timepoints(
    u_scaled, s, alpha, beta, gamma, t_, u0_, s0_,
    mode = if (beta < gamma) "projection" else "inverse"
  )
  
  t <- assignment$tau * assignment$o + (assignment$tau_ + t_) * (1 - assignment$o)
  tau <- assignment$tau
  o <- assignment$o
  
  # Calculate initial loss
  loss <- compute_loss(u_scaled, s, t, t_, alpha, beta, gamma, std_u / scaling, std_s)
  
  result$recoverable <- TRUE
  
  if (max_iter > 0) {
    # EM optimization
    for (iter in seq_len(max_iter)) {
      # Store previous values
      alpha_old <- alpha
      beta_old <- beta
      gamma_old <- gamma
      t__old <- t_
      loss_old <- loss
      
      # Optimize t_ and alpha
      opt_result <- optim(
        par = c(t_, alpha),
        fn = function(x) {
          compute_loss(u_scaled, s, t, x[1], x[2], beta, gamma, 
                      std_u / scaling, std_s, refit_time = TRUE)
        },
        method = "Nelder-Mead",
        control = list(maxit = max_iter %/% 5)
      )
      
      if (opt_result$value < loss) {
        t_ <- opt_result$par[1]
        alpha <- opt_result$par[2]
        loss <- opt_result$value
      }
      
      # Optimize scaling if requested
      if (isTRUE(fit_scaling)) {
        opt_result <- optim(
          par = c(t_, beta, scaling),
          fn = function(x) {
            compute_loss(u / x[3], s, t, x[1], alpha, x[2], gamma,
                        std_u / x[3], std_s, refit_time = TRUE)
          },
          method = "Nelder-Mead",
          control = list(maxit = max_iter %/% 5)
        )
        
        if (opt_result$value < loss) {
          t_ <- opt_result$par[1]
          beta <- opt_result$par[2]
          scaling <- opt_result$par[3]
          u_scaled <- u / scaling
          loss <- opt_result$value
        }
      }
      
      # Optimize rates
      opt_result <- optim(
        par = c(alpha, gamma),
        fn = function(x) {
          compute_loss(u_scaled, s, t, t_, x[1], beta, x[2],
                      std_u / scaling, std_s, refit_time = TRUE)
        },
        method = "Nelder-Mead",
        control = list(maxit = max_iter %/% 5)
      )
      
      if (opt_result$value < loss) {
        alpha <- opt_result$par[1]
        gamma <- opt_result$par[2]
        loss <- opt_result$value
      }
      
      # Optimize t_
      opt_result <- optim(
        par = t_,
        fn = function(x) {
          compute_loss(u_scaled, s, t, x[1], alpha, beta, gamma,
                      std_u / scaling, std_s, refit_time = TRUE)
        },
        method = "Nelder-Mead",
        control = list(maxit = max_iter %/% 5)
      )
      
      if (opt_result$value < loss) {
        t_ <- opt_result$par[1]
        loss <- opt_result$value
      }
      
      # Optimize all rates together
      opt_result <- optim(
        par = c(t_, alpha, beta, gamma),
        fn = function(x) {
          compute_loss(u_scaled, s, t, x[1], x[2], x[3], x[4],
                      std_u / scaling, std_s, refit_time = TRUE)
        },
        method = "Nelder-Mead",
        control = list(maxit = max_iter %/% 5)
      )
      
      if (opt_result$value < loss) {
        t_ <- opt_result$par[1]
        alpha <- opt_result$par[2]
        beta <- opt_result$par[3]
        gamma <- opt_result$par[4]
        loss <- opt_result$value
      }
      
      # Update time assignment
      result_switch <- mRNA(t_, 0, 0, alpha, beta, gamma)
      u0_ <- result_switch$u
      s0_ <- result_switch$s
      
      assignment <- assign_timepoints(
        u_scaled, s, alpha, beta, gamma, t_, u0_, s0_,
        mode = if (beta < gamma) "projection" else "inverse"
      )
      
      t <- assignment$tau * assignment$o + (assignment$tau_ + t_) * (1 - assignment$o)
      tau <- assignment$tau
      o <- assignment$o
      
      # Check convergence
      if (abs(loss_old - loss) < 1e-6 * loss_old) {
        break
      }
    }
  }
  
  # Final time assignment
  result_switch <- mRNA(t_, 0, 0, alpha, beta, gamma)
  u0_ <- result_switch$u
  s0_ <- result_switch$s
  
  assignment <- assign_timepoints(
    u_scaled, s, alpha, beta, gamma, t_, u0_, s0_,
    mode = "projection"
  )
  
  t <- assignment$tau * assignment$o + (assignment$tau_ + t_) * (1 - assignment$o)
  tau <- assignment$tau
  o <- assignment$o
  
  # Compute final likelihood
  likelihood <- compute_likelihood(
    u_scaled, s, t, t_, alpha, beta, gamma,
    std_u / scaling, std_s
  )
  
  # Compute variance
  variance <- compute_variance(
    u_scaled, s, t, t_, alpha, beta, gamma,
    std_u / scaling, std_s
  )
  
  # Return results
  result$alpha <- alpha
  result$beta <- beta
  result$gamma <- gamma
  result$t_ <- t_
  result$scaling <- scaling
  result$std_u <- std_u
  result$std_s <- std_s
  result$likelihood <- likelihood
  result$variance <- variance
  result$t <- t
  result$tau <- tau
  result$o <- o
  
  result
}

#' Compute Loss Function
#'
#' @description Compute sum of squared errors between observations and model.
#'
#' @param u Scaled unspliced
#' @param s Spliced
#' @param t Cell times
#' @param t_ Switching time
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#' @param std_u Std of unspliced
#' @param std_s Std of spliced
#' @param refit_time Whether to refit time
#'
#' @return Loss value
#' @keywords internal
compute_loss <- function(u, s, t, t_, alpha, beta, gamma,
                         std_u, std_s, refit_time = FALSE) {
  
  # Ensure positive parameters
  if (alpha <= 0 || beta <= 0 || gamma <= 0 || t_ <= 0) {
    return(Inf)
  }
  
  if (refit_time) {
    # Recompute time assignment
    result_switch <- mRNA(t_, 0, 0, alpha, beta, gamma)
    u0_ <- result_switch$u
    s0_ <- result_switch$s
    
    assignment <- assign_timepoints(
      u, s, alpha, beta, gamma, t_, u0_, s0_,
      mode = "inverse"
    )
    
    t <- assignment$tau * assignment$o + (assignment$tau_ + t_) * (1 - assignment$o)
  }
  
  # Vectorize parameters
  params <- vectorize_params(t, t_, alpha, beta, gamma)
  
  # Compute model predictions
  result <- mRNA(params$tau, params$u0, params$s0, params$alpha, beta, gamma)
  ut <- result$u
  st <- result$s
  
  # Compute normalized residuals
  u_diff <- (u - ut) / std_u
  s_diff <- (s - st) / std_s
  
  # Sum of squared errors
  sum(u_diff^2 + s_diff^2, na.rm = TRUE)
}

#' Compute Likelihood
#'
#' @description Compute model likelihood.
#'
#' @param u Scaled unspliced
#' @param s Spliced
#' @param t Cell times
#' @param t_ Switching time
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#' @param std_u Std of unspliced
#' @param std_s Std of spliced
#'
#' @return Likelihood value
#' @keywords internal
compute_likelihood <- function(u, s, t, t_, alpha, beta, gamma,
                                std_u, std_s) {
  
  # Vectorize parameters
  params <- vectorize_params(t, t_, alpha, beta, gamma)
  
  # Compute model predictions
  result <- mRNA(params$tau, params$u0, params$s0, params$alpha, beta, gamma)
  ut <- result$u
  st <- result$s
  
  # Compute normalized residuals
  u_diff <- (u - ut) / std_u
  s_diff <- (s - st) / std_s
  
  distx <- u_diff^2 + s_diff^2
  
  # Compute variance
  n <- sum(!is.na(distx))
  n <- max(n - n * 0.01, 2)
  
  varx <- mean(distx, na.rm = TRUE) - mean(sign(s_diff) * sqrt(distx), na.rm = TRUE)^2
  varx <- max(varx, 1e-10)
  
  # Log likelihood
  log_lik <- -1 / 2 / n * sum(distx, na.rm = TRUE) / varx - 1 / 2 * log(2 * pi * varx)
  
  exp(log_lik)
}

#' Compute Variance
#'
#' @description Compute model variance for goodness-of-fit.
#'
#' @param u Scaled unspliced
#' @param s Spliced
#' @param t Cell times
#' @param t_ Switching time
#' @param alpha Transcription rate
#' @param beta Splicing rate
#' @param gamma Degradation rate
#' @param std_u Std of unspliced
#' @param std_s Std of spliced
#'
#' @return Variance value
#' @keywords internal
compute_variance <- function(u, s, t, t_, alpha, beta, gamma,
                              std_u, std_s) {
  
  # Vectorize parameters
  params <- vectorize_params(t, t_, alpha, beta, gamma)
  
  # Compute model predictions
  result <- mRNA(params$tau, params$u0, params$s0, params$alpha, beta, gamma)
  ut <- result$u
  st <- result$s
  
  # Compute normalized residuals
  u_diff <- (u - ut) / std_u
  s_diff <- (s - st) / std_s
  
  distx <- u_diff^2 + s_diff^2
  
  mean(distx, na.rm = TRUE) - mean(sign(s_diff) * sqrt(distx), na.rm = TRUE)^2
}

#' Align Dynamics
#'
#' @description Align gene-specific dynamics to a common time scale.
#'
#' @param T_matrix Matrix of cell times (cells x genes)
#' @param t_ Vector of switching times
#' @param alpha Vector of transcription rates
#' @param beta Vector of splicing rates
#' @param gamma Vector of degradation rates
#' @param scaling Vector of scaling factors
#' @param t_max Target maximum time
#'
#' @return List with aligned parameters
#' @keywords internal
align_dynamics <- function(T_matrix, t_, alpha, beta, gamma, scaling,
                           t_max = 20) {
  
  n_genes <- ncol(T_matrix)
  
  # Compute maximum time for each gene
  T_max <- numeric(n_genes)
  
  for (i in seq_len(n_genes)) {
    if (is.na(t_[i])) next
    
    t_col <- T_matrix[, i]
    valid <- !is.na(t_col)
    
    # Time in "on" phase (before switching)
    on_mask <- t_col < t_[i]
    T_on <- if (any(on_mask & valid)) max(t_col[on_mask & valid]) else 0
    
    # Time in "off" phase (after switching)
    off_mask <- t_col >= t_[i]
    T_off <- if (any(off_mask & valid)) max(t_col[off_mask & valid] - t_[i]) else 0
    
    T_max[i] <- T_on + T_off
    
    # Adjust for steady state cells
    n_steady <- sum(t_col == t_[i] | t_col == 0, na.rm = TRUE)
    if (n_steady > 0 && T_max[i] > 0) {
      T_max[i] <- T_max[i] / (1 - n_steady / sum(valid))
    }
  }
  
  T_max[T_max == 0 | !is.finite(T_max)] <- 1
  
  # Compute alignment scaling
  m <- t_max / T_max
  m[!is.finite(m)] <- 1
  
  # Apply scaling
  for (i in seq_len(n_genes)) {
    if (is.na(t_[i]) || m[i] == 1) next
    
    T_matrix[, i] <- T_matrix[, i] * m[i]
    t_[i] <- t_[i] * m[i]
    alpha[i] <- alpha[i] / m[i]
    beta[i] <- beta[i] / m[i]
    gamma[i] <- gamma[i] / m[i]
  }
  
  alignment_scaling <- m
  alignment_scaling[alignment_scaling == 1] <- NA
  
  list(
    T = T_matrix,
    t_ = t_,
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    alignment_scaling = alignment_scaling
  )
}

#' Compute Dynamical Velocity
#'
#' @description Compute velocity using parameters from the dynamical model.
#'
#' @param object Seurat object with dynamics fitted
#' @param mode Velocity mode: "soft", "hard", or "deterministic"
#'
#' @return Seurat object with velocity_s in misc$scVeloR$velocity
#' @export
velocity_from_dynamics <- function(object, mode = "deterministic") {
  
  if (is.null(object@misc$scVeloR$dynamics)) {
    stop("Run recover_dynamics() first")
  }
  
  dynamics <- object@misc$scVeloR$dynamics
  gene_params <- dynamics$gene_params
  fit_t <- dynamics$fit_t
  genes <- dynamics$genes
  gene_indices <- dynamics$gene_indices
  
  # Get data
  data <- get_velocity_data(object)
  Ms <- data$Ms
  Mu <- data$Mu
  
  n_cells <- nrow(Ms)
  n_genes <- length(genes)
  
  # Initialize velocity matrix
  velocity_s <- matrix(0, n_cells, ncol(Ms))
  colnames(velocity_s) <- colnames(Ms)
  
  for (i in seq_len(n_genes)) {
    idx <- gene_indices[i]
    
    alpha <- gene_params$alpha[i]
    beta <- gene_params$beta[i] * gene_params$scaling[i]  # Scale back
    gamma <- gene_params$gamma[i]
    t_ <- gene_params$t_[i]
    scaling <- gene_params$scaling[i]
    
    if (is.na(alpha)) next
    
    t <- fit_t[, i]
    u <- Mu[, idx] / scaling
    s <- Ms[, idx]
    
    # Vectorize based on time
    params <- vectorize_params(t, t_, alpha, beta, gamma)
    
    # Compute model predictions
    result <- mRNA(params$tau, params$u0, params$s0, params$alpha, beta, gamma)
    ut <- result$u
    st <- result$s
    
    # Compute velocity: ds/dt = beta*u - gamma*s
    v_s <- ut * beta - st * gamma
    
    # Clip velocity
    v_s <- pmax(v_s, -s)
    
    velocity_s[, idx] <- v_s
  }
  
  # Store results
  if (is.null(object@misc$scVeloR$velocity)) {
    object@misc$scVeloR$velocity <- list()
  }
  
  object@misc$scVeloR$velocity$velocity_s <- velocity_s
  object@misc$scVeloR$velocity$velocity_genes <- gene_indices
  object@misc$scVeloR$velocity$mode <- "dynamical"
  
  object
}

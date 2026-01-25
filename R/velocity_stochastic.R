#' @title Stochastic Velocity Model
#' @description Functions for estimating RNA velocity using the stochastic model.
#' @name velocity_stochastic
NULL

#' Fit Stochastic Velocity Model
#'
#' @description Estimate RNA velocity using second-order moments.
#' This model accounts for transcriptional variability.
#'
#' @param seurat_obj A Seurat object with computed moments (Ms, Mu)
#' @param fit_offset Whether to fit offset (default: FALSE)
#' @param fit_offset2 Whether to fit second offset (default: FALSE)
#' @param min_r2 Minimum R-squared to consider a gene well-fit (default: 0.01)
#' @param n_jobs Number of parallel jobs (default: 1)
#'
#' @return Modified Seurat object with velocity and fit parameters
#'
#' @export
fit_velocity_stochastic <- function(seurat_obj,
                                     fit_offset = FALSE,
                                     fit_offset2 = FALSE,
                                     min_r2 = 0.01,
                                     n_jobs = 1L) {
  
  .vmessage("Fitting stochastic velocity model")
  
  # Check for moments
  if (!.layer_exists(seurat_obj, "Ms") || !.layer_exists(seurat_obj, "Mu")) {
    stop("Moments not found. Run compute_moments first.", call. = FALSE)
  }
  
  Ms <- .get_layer_data(seurat_obj, "Ms")
  Mu <- .get_layer_data(seurat_obj, "Mu")
  
  Ms <- make_dense(Ms)
  Mu <- make_dense(Mu)
  
  # Compute second-order moments
  .vmessage("Computing second-order moments...")
  second_moments <- compute_second_order_moments(seurat_obj)
  Mss <- second_moments$Mss
  Mus <- second_moments$Mus
  
  # Get genes
  genes <- rownames(Ms)
  n_genes <- length(genes)
  n_cells <- ncol(Ms)
  
  # Get velocity genes mask
  velocity_genes <- get_velocity_genes_mask(seurat_obj, genes)
  
  # Initialize results
  gamma <- numeric(n_genes)
  gamma2 <- numeric(n_genes)
  offset <- numeric(n_genes)
  offset2 <- numeric(n_genes)
  r2 <- numeric(n_genes)
  residual <- numeric(n_genes)
  velocity_genes_result <- rep(FALSE, n_genes)
  
  names(gamma) <- genes
  names(r2) <- genes
  names(velocity_genes_result) <- genes
  
  # Fit each gene
  .vmessage("Fitting ", sum(velocity_genes), " genes...")
  
  fit_gene_stochastic <- function(gene_idx) {
    
    s <- Ms[gene_idx, ]
    u <- Mu[gene_idx, ]
    ss <- Mss[gene_idx, ]
    us <- Mus[gene_idx, ]
    
    # Stochastic model uses second-order moments
    # du = u - gamma * s - offset
    # ds = s^2 - gamma2 * s - offset2
    
    result <- fit_stochastic_single(s, u, ss, us, 
                                     fit_offset, fit_offset2)
    
    return(result)
  }
  
  # Run fitting
  results <- apply_with_progress(
    seq_len(n_genes),
    fit_gene_stochastic,
    n_jobs = n_jobs,
    show_progress = TRUE
  )
  
  # Extract results
  for (i in seq_len(n_genes)) {
    if (!is.null(results[[i]])) {
      gamma[i] <- results[[i]]$gamma
      gamma2[i] <- results[[i]]$gamma2
      offset[i] <- results[[i]]$offset
      offset2[i] <- results[[i]]$offset2
      r2[i] <- results[[i]]$r2
      residual[i] <- results[[i]]$residual
      velocity_genes_result[i] <- velocity_genes[i] && 
                                   !is.na(gamma[i]) && 
                                   gamma[i] > 0 &&
                                   r2[i] >= min_r2
    }
  }
  
  # Calculate velocity
  .vmessage("Computing velocity...")
  
  velocity <- Mu - sweep(Ms, 1, gamma, "*")
  if (fit_offset) {
    velocity <- sweep(velocity, 1, offset, "-")
  }
  
  # Set non-velocity genes to NA
  velocity[!velocity_genes_result, ] <- NA
  
  # Store results
  velocity <- as(velocity, "CsparseMatrix")
  rownames(velocity) <- genes
  colnames(velocity) <- colnames(Ms)
  
  seurat_obj <- .set_layer_data(seurat_obj, "velocity", velocity)
  
  # Store fit parameters
  fit_params <- data.frame(
    row.names = genes,
    fit_gamma = gamma,
    fit_gamma2 = gamma2,
    fit_offset = offset,
    fit_offset2 = offset2,
    fit_r2 = r2,
    fit_residual = residual,
    velocity_genes = velocity_genes_result
  )
  
  seurat_obj <- store_fit_params(seurat_obj, fit_params)
  seurat_obj@misc$velocity_mode <- "stochastic"
  
  n_fit <- sum(velocity_genes_result, na.rm = TRUE)
  .vmessage("Fitted ", n_fit, " velocity genes")
  
  return(seurat_obj)
}

#' Fit Stochastic Model for Single Gene
#'
#' @param s Spliced moments
#' @param u Unspliced moments
#' @param ss Spliced second moments
#' @param us Cross moments
#' @param fit_offset Whether to fit offset
#' @param fit_offset2 Whether to fit second offset
#' @return List with fit results
#' @keywords internal
fit_stochastic_single <- function(s, u, ss, us, fit_offset, fit_offset2) {
  
  n <- length(s)
  
  # Handle zeros and NAs
  valid <- !is.na(s) & !is.na(u) & !is.na(ss) & !is.na(us) &
           (s > 0 | u > 0)
  
  if (sum(valid) < 10) {
    return(list(gamma = NA, gamma2 = NA, offset = NA, offset2 = NA, 
                r2 = 0, residual = Inf))
  }
  
  s <- s[valid]
  u <- u[valid]
  ss <- ss[valid]
  us <- us[valid]
  
  # Compute variances
  var_s <- ss - s^2
  cov_us <- us - u * s
  
  # Weighted regression
  # gamma = cov(u,s) / var(s)
  var_s_sum <- sum(pmax(var_s, 1e-10))
  cov_us_sum <- sum(cov_us)
  
  gamma <- max(0, cov_us_sum / var_s_sum)
  
  # For second moment: ss = gamma2 * s + offset2
  # gamma2 = sum(s * ss) / sum(s^2)
  s2_sum <- sum(s^2)
  if (s2_sum > 0) {
    gamma2 <- sum(s * ss) / s2_sum
  } else {
    gamma2 <- 0
  }
  
  # Offsets
  if (fit_offset) {
    offset <- mean(u) - gamma * mean(s)
  } else {
    offset <- 0
  }
  
  if (fit_offset2) {
    offset2 <- mean(ss) - gamma2 * mean(s)
  } else {
    offset2 <- 0
  }
  
  # Calculate R-squared
  u_pred <- gamma * s + offset
  ss_res <- sum((u - u_pred)^2)
  ss_tot <- sum((u - mean(u))^2)
  r2 <- ifelse(ss_tot > 0, max(0, 1 - ss_res / ss_tot), 0)
  
  residual <- sqrt(ss_res / length(u))
  
  list(
    gamma = gamma,
    gamma2 = gamma2,
    offset = offset,
    offset2 = offset2,
    r2 = r2,
    residual = residual
  )
}

#' Generalized Least Squares for Velocity
#'
#' @description Fit velocity using generalized least squares with variance weighting.
#'
#' @param s Spliced values
#' @param u Unspliced values
#' @param ss Second moment of spliced
#' @param us Cross moment
#' @param fit_offset Whether to fit offset
#' @return List with gamma, offset, and r2
#' @keywords internal
leastsq_generalized <- function(s, u, ss, us, fit_offset = FALSE) {
  
  # Compute variance weights
  var_s <- pmax(ss - s^2, 1e-10)
  w <- 1 / var_s
  w <- w / sum(w)
  
  # Weighted regression
  sw <- s * w
  uw <- u * w
  
  if (fit_offset) {
    # With intercept
    X <- cbind(1, s)
    Xw <- X * sqrt(w)
    yw <- u * sqrt(w)
    
    fit <- tryCatch({
      coef <- solve(t(Xw) %*% Xw) %*% t(Xw) %*% yw
      list(gamma = coef[2], offset = coef[1])
    }, error = function(e) {
      list(gamma = sum(sw * u) / sum(sw * s), offset = 0)
    })
  } else {
    # No intercept
    fit <- list(
      gamma = sum(sw * u) / sum(sw * s),
      offset = 0
    )
  }
  
  # R-squared
  pred <- fit$gamma * s + fit$offset
  ss_res <- sum(w * (u - pred)^2)
  ss_tot <- sum(w * (u - sum(uw))^2)
  r2 <- max(0, 1 - ss_res / ss_tot)
  
  fit$r2 <- r2
  return(fit)
}

#' Maximum Likelihood Velocity Estimation
#'
#' @description Estimate velocity using maximum likelihood with assumed noise model.
#'
#' @param u Unspliced values
#' @param s Spliced values
#' @param fit_offset Whether to fit offset
#' @return List with gamma, offset, and likelihood
#' @keywords internal
maximum_likelihood <- function(u, s, fit_offset = FALSE) {
  
  n <- length(u)
  
  # Initial estimates
  if (fit_offset) {
    init <- stats::lm(u ~ s)
    gamma0 <- stats::coef(init)[2]
    offset0 <- stats::coef(init)[1]
  } else {
    gamma0 <- sum(s * u) / sum(s^2)
    offset0 <- 0
  }
  
  # Negative log-likelihood assuming Gaussian noise
  neg_log_lik <- function(par) {
    if (fit_offset) {
      gamma <- par[1]
      offset <- par[2]
    } else {
      gamma <- par[1]
      offset <- 0
    }
    
    residual <- u - gamma * s - offset
    sigma <- max(1e-10, stats::sd(residual))
    
    -sum(stats::dnorm(residual, mean = 0, sd = sigma, log = TRUE))
  }
  
  # Optimize
  if (fit_offset) {
    opt <- stats::optim(c(gamma0, offset0), neg_log_lik, method = "L-BFGS-B",
                        lower = c(0, -Inf), upper = c(Inf, Inf))
    gamma <- opt$par[1]
    offset <- opt$par[2]
  } else {
    opt <- stats::optim(gamma0, neg_log_lik, method = "L-BFGS-B",
                        lower = 0, upper = Inf)
    gamma <- opt$par[1]
    offset <- 0
  }
  
  likelihood <- -opt$value
  
  list(gamma = gamma, offset = offset, likelihood = likelihood)
}

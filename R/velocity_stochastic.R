#' @title Stochastic Velocity Model
#' @description Estimate RNA velocity using the stochastic model.
#' The stochastic model uses second-order moments (variance and covariance)
#' to estimate velocity, accounting for transcriptional bursting.
#' @name velocity_stochastic
NULL

#' Compute Stochastic Velocity
#'
#' @description Estimate RNA velocity using the stochastic model with second-order moments.
#' This model accounts for noise and bursting dynamics.
#'
#' @param object Seurat object with moments computed
#' @param min_r2 Minimum R-squared for velocity genes
#' @param use_raw Use raw counts
#' @param n_cores Number of cores
#' @param verbose Print progress
#'
#' @return Seurat object with velocity in misc$scVeloR$velocity
#' @export
velocity_stochastic <- function(object,
                                 min_r2 = 0.01,
                                 use_raw = FALSE,
                                 n_cores = 1,
                                 verbose = TRUE) {
  
  if (verbose) {
    message("Computing stochastic velocity...")
  }
  
  # Get first-order moments
  data <- get_velocity_data(object, use_moments = !use_raw)
  Ms <- data$Ms
  Mu <- data$Mu
  
  # Get second-order moments
  if (is.null(object@misc$scVeloR$Mss)) {
    object <- second_order_moments(object, verbose = verbose)
  }
  
  Mss <- object@misc$scVeloR$Mss
  Mus <- object@misc$scVeloR$Mus
  Muu <- object@misc$scVeloR$Muu
  
  n_cells <- nrow(Ms)
  n_genes <- ncol(Ms)
  gene_names <- colnames(Ms)
  
  # Initialize results
  velocity <- matrix(0, n_cells, n_genes)
  gamma <- rep(NA_real_, n_genes)
  r2 <- rep(NA_real_, n_genes)
  velocity_genes <- rep(FALSE, n_genes)
  
  # Fit each gene using stochastic model
  fit_gene_stochastic <- function(i) {
    s <- Ms[, i]
    u <- Mu[, i]
    ss <- Mss[, i]
    us <- Mus[, i]
    uu <- Muu[, i]
    
    # Skip if too few non-zero values
    valid <- (s > 0) & (u > 0)
    if (sum(valid) < 10) {
      return(list(gamma = NA, r2 = NA, v = rep(0, n_cells)))
    }
    
    s_v <- s[valid]
    u_v <- u[valid]
    ss_v <- ss[valid]
    us_v <- us[valid]
    uu_v <- uu[valid]
    
    # Compute variance and covariance
    # Var(s) = <s^2> - <s>^2
    # Cov(u,s) = <us> - <u><s>
    var_s <- ss_v - s_v^2
    cov_us <- us_v - u_v * s_v
    
    # Stochastic gamma estimation:
    # gamma = Cov(u,s) / Var(s)
    # This comes from the fact that at steady state,
    # the covariance captures the correlation structure
    
    # Use cells with positive variance
    pos_var <- var_s > 1e-10
    
    if (sum(pos_var) < 10) {
      # Fall back to steady-state method
      gamma_i <- sum(u_v * s_v) / sum(s_v^2)
    } else {
      gamma_i <- sum(cov_us[pos_var]) / sum(var_s[pos_var])
    }
    
    # Ensure positive gamma
    if (!is.finite(gamma_i) || gamma_i <= 0) {
      return(list(gamma = NA, r2 = NA, v = rep(0, n_cells)))
    }
    
    # Compute velocity
    v <- u - gamma_i * s
    
    # R-squared (using deterministic prediction)
    u_pred <- gamma_i * s_v
    ss_res <- sum((u_v - u_pred)^2)
    ss_tot <- sum((u_v - mean(u_v))^2)
    r2_i <- 1 - ss_res / max(ss_tot, 1e-10)
    
    list(gamma = gamma_i, r2 = r2_i, v = v)
  }
  
  # Run fitting
  if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    results <- parallel::mclapply(seq_len(n_genes), fit_gene_stochastic, mc.cores = n_cores)
  } else {
    results <- lapply(seq_len(n_genes), function(i) {
      if (verbose && i %% 1000 == 0) {
        message(sprintf("  Fitted %d/%d genes", i, n_genes))
      }
      fit_gene_stochastic(i)
    })
  }
  
  # Collect results
  for (i in seq_len(n_genes)) {
    res <- results[[i]]
    gamma[i] <- res$gamma
    r2[i] <- res$r2
    velocity[, i] <- res$v
    
    if (!is.na(res$r2) && res$r2 >= min_r2) {
      velocity_genes[i] <- TRUE
    }
  }
  
  colnames(velocity) <- gene_names
  rownames(velocity) <- rownames(Ms)
  names(gamma) <- gene_names
  names(r2) <- gene_names
  
  n_velocity_genes <- sum(velocity_genes)
  
  if (verbose) {
    message(sprintf("  Found %d velocity genes (R² >= %.2f)", n_velocity_genes, min_r2))
  }
  
  # Store results
  if (is.null(object@misc$scVeloR)) {
    object@misc$scVeloR <- list()
  }
  
  object@misc$scVeloR$velocity <- list(
    velocity_s = velocity,
    gamma = gamma,
    r2 = r2,
    velocity_genes = which(velocity_genes),
    mode = "stochastic",
    min_r2 = min_r2
  )
  
  if (verbose) {
    message("Done. Velocity stored in object@misc$scVeloR$velocity")
  }
  
  object
}

#' Compute Variance Velocity
#'
#' @description Compute velocity based on variance decomposition.
#' This is an alternative formulation of the stochastic model.
#'
#' @param object Seurat object with second-order moments
#' @param verbose Print progress
#'
#' @return Seurat object with variance-based velocity
#' @keywords internal
variance_velocity <- function(object, verbose = TRUE) {
  
  if (is.null(object@misc$scVeloR$Mss)) {
    object <- second_order_moments(object, verbose = verbose)
  }
  
  # Get moments
  Ms <- object@misc$scVeloR$Ms
  Mu <- object@misc$scVeloR$Mu
  Mss <- object@misc$scVeloR$Mss
  Mus <- object@misc$scVeloR$Mus
  Muu <- object@misc$scVeloR$Muu
  
  n_cells <- nrow(Ms)
  n_genes <- ncol(Ms)
  
  # Compute variances and covariances
  var_s <- Mss - Ms^2
  var_u <- Muu - Mu^2
  cov_us <- Mus - Mu * Ms
  
  # Velocity based on variance ratio
  # This captures the directional change in expression variance
  velocity <- matrix(0, n_cells, n_genes)
  
  for (i in seq_len(n_genes)) {
    var_s_i <- var_s[, i]
    cov_us_i <- cov_us[, i]
    
    # Gamma from covariance/variance ratio
    gamma_i <- sum(cov_us_i, na.rm = TRUE) / sum(var_s_i + 1e-10, na.rm = TRUE)
    
    if (is.finite(gamma_i) && gamma_i > 0) {
      velocity[, i] <- Mu[, i] - gamma_i * Ms[, i]
    }
  }
  
  colnames(velocity) <- colnames(Ms)
  rownames(velocity) <- rownames(Ms)
  
  object@misc$scVeloR$variance_velocity <- velocity
  
  object
}

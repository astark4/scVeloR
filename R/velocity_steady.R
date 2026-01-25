#' @title Steady-State Velocity Model
#' @description Estimate RNA velocity using the steady-state model.
#' The steady-state model assumes that cells are near transcriptional equilibrium:
#' ds/dt = beta*u - gamma*s ≈ 0, which gives gamma = beta * u/s.
#' Velocity is then computed as v = beta*u - gamma*s.
#' @name velocity_steady
NULL

#' Compute Steady-State Velocity
#'
#' @description Estimate RNA velocity using linear regression on steady-state assumption.
#' For each gene, fit gamma = u/s ratio using cells at extreme expression values.
#'
#' @param object Seurat object with moments computed
#' @param min_r2 Minimum R-squared for a gene to be considered a velocity gene
#' @param perc_range Percentile range for extreme quantile fitting
#' @param fit_offset Whether to fit intercept (for imputation model)
#' @param use_raw Use raw counts instead of moments
#' @param n_cores Number of cores for parallel computation
#' @param verbose Print progress
#'
#' @return Seurat object with velocity in misc$scVeloR$velocity
#' @export
velocity_steady_state <- function(object,
                                   min_r2 = 0.01,
                                   perc_range = c(5, 95),
                                   fit_offset = FALSE,
                                   use_raw = FALSE,
                                   n_cores = 1,
                                   verbose = TRUE) {
  
  if (verbose) {
    message("Computing steady-state velocity...")
  }
  
  # Get expression data
  data <- get_velocity_data(object, use_moments = !use_raw)
  Ms <- data$Ms
  Mu <- data$Mu
  
  n_cells <- nrow(Ms)
  n_genes <- ncol(Ms)
  gene_names <- colnames(Ms)
  
  # Initialize result matrices
  velocity <- matrix(0, n_cells, n_genes)
  gamma <- rep(NA_real_, n_genes)
  r2 <- rep(NA_real_, n_genes)
  velocity_genes <- rep(FALSE, n_genes)
  
  # Percentile thresholds
  perc_lo <- perc_range[1]
  perc_hi <- perc_range[2]
  
  # Fit each gene
  fit_gene <- function(i) {
    s <- Ms[, i]
    u <- Mu[, i]
    
    # Skip if too few non-zero values
    valid <- (s > 0) & (u > 0)
    if (sum(valid) < 10) {
      return(list(gamma = NA, r2 = NA, v = rep(0, n_cells)))
    }
    
    s_valid <- s[valid]
    u_valid <- u[valid]
    
    # Use extreme quantiles for fitting
    s_lo <- quantile(s_valid, perc_lo / 100)
    s_hi <- quantile(s_valid, perc_hi / 100)
    
    extreme <- (s_valid <= s_lo) | (s_valid >= s_hi)
    
    if (sum(extreme) < 10) {
      extreme <- rep(TRUE, length(s_valid))
    }
    
    s_fit <- s_valid[extreme]
    u_fit <- u_valid[extreme]
    
    # Linear regression: u = gamma * s + offset
    # (with or without offset)
    if (fit_offset) {
      fit <- lm(u_fit ~ s_fit)
      gamma_i <- coef(fit)[2]
      offset_i <- coef(fit)[1]
      
      u_pred <- gamma_i * s_fit + offset_i
    } else {
      # Regression through origin: gamma = sum(u*s) / sum(s^2)
      gamma_i <- sum(u_fit * s_fit) / sum(s_fit^2)
      offset_i <- 0
      
      u_pred <- gamma_i * s_fit
    }
    
    # R-squared
    ss_res <- sum((u_fit - u_pred)^2)
    ss_tot <- sum((u_fit - mean(u_fit))^2)
    r2_i <- 1 - ss_res / max(ss_tot, 1e-10)
    
    # Ensure gamma is positive
    if (!is.finite(gamma_i) || gamma_i <= 0) {
      return(list(gamma = NA, r2 = NA, v = rep(0, n_cells)))
    }
    
    # Compute velocity: v = u - gamma * s
    # (note: this is proportional to ds/dt since ds/dt = beta*u - gamma*s,
    # and we assume beta = 1 or absorb it into the ratio)
    v <- u - gamma_i * s
    
    list(gamma = gamma_i, r2 = r2_i, v = v)
  }
  
  # Run fitting
  if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    results <- parallel::mclapply(seq_len(n_genes), fit_gene, mc.cores = n_cores)
  } else {
    results <- lapply(seq_len(n_genes), function(i) {
      if (verbose && i %% 1000 == 0) {
        message(sprintf("  Fitted %d/%d genes", i, n_genes))
      }
      fit_gene(i)
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
  names(velocity_genes) <- gene_names
  
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
    mode = "steady_state",
    min_r2 = min_r2
  )
  
  if (verbose) {
    message("Done. Velocity stored in object@misc$scVeloR$velocity")
  }
  
  object
}

#' Rank Velocity Genes
#'
#' @description Rank genes by velocity fit quality (R-squared).
#'
#' @param object Seurat object with velocity computed
#' @param n_top Number of top genes to return
#' @param groupby Optional grouping variable
#'
#' @return Data frame with ranked genes
#' @export
rank_velocity_genes <- function(object,
                                n_top = 100,
                                groupby = NULL) {
  
  if (is.null(object@misc$scVeloR$velocity)) {
    stop("Run velocity() first")
  }
  
  r2 <- object@misc$scVeloR$velocity$r2
  
  # Remove NA values
  valid <- !is.na(r2)
  r2_valid <- r2[valid]
  gene_names <- names(r2)[valid]
  
  # Sort by R2
  order_idx <- order(r2_valid, decreasing = TRUE)
  
  n_return <- min(n_top, length(order_idx))
  
  data.frame(
    gene = gene_names[order_idx[1:n_return]],
    r2 = r2_valid[order_idx[1:n_return]],
    rank = 1:n_return,
    stringsAsFactors = FALSE
  )
}

#' Test Velocity Gene Bimodality
#'
#' @description Test whether a gene shows bimodal expression distribution,
#' which may indicate that it has both inducing and repressing populations.
#'
#' @param object Seurat object
#' @param gene Gene name
#'
#' @return List with test results
#' @keywords internal
test_bimodality <- function(object, gene) {
  
  data <- get_velocity_data(object)
  Ms <- data$Ms
  Mu <- data$Mu
  
  gene_idx <- which(colnames(Ms) == gene)
  
  if (length(gene_idx) == 0) {
    return(list(pval = 1, bimodal = FALSE))
  }
  
  s <- Ms[, gene_idx]
  u <- Mu[, gene_idx]
  
  # Use kernel density estimation
  if (sum(s > 0) < 20) {
    return(list(pval = 1, bimodal = FALSE))
  }
  
  s_pos <- s[s > 0]
  
  # Estimate density
  dens <- density(s_pos, n = 512)
  
  # Find local maxima
  d_diff <- diff(dens$y)
  local_max <- which(d_diff[-length(d_diff)] > 0 & d_diff[-1] < 0) + 1
  
  # Bimodal if two or more significant peaks
  bimodal <- length(local_max) >= 2
  
  list(
    pval = if (bimodal) 0.01 else 1,
    bimodal = bimodal,
    n_peaks = length(local_max)
  )
}

#' @title Dynamical Velocity Model
#' @description Functions for recovering RNA dynamics using the EM algorithm.
#' @name velocity_dynamical
NULL

#' Recover Dynamics
#'
#' @description Recover the full splicing dynamics using an expectation-maximization
#' algorithm. This fits transcription, splicing, and degradation rates for each gene.
#'
#' @param seurat_obj A Seurat object with computed moments
#' @param var_names Gene names to fit, or "velocity_genes" to use filtered genes
#' @param n_top_genes Number of top genes to use if var_names is numeric
#' @param max_iter Maximum EM iterations per gene (default: 10)
#' @param fit_time Whether to fit time for each cell (default: TRUE)
#' @param fit_scaling Whether to fit scaling between u and s (default: TRUE)
#' @param fit_steady_states Fit steady states (default: TRUE)
#' @param fit_connected_states Only fit states with neighbors (default: FALSE)
#' @param t_max Maximum time to consider (default: 20)
#' @param eps Convergence threshold (default: 1e-4)
#' @param n_jobs Number of parallel jobs (default: 1)
#'
#' @return Modified Seurat object with fitted dynamics parameters
#'
#' @export
recover_dynamics <- function(seurat_obj,
                               var_names = "velocity_genes",
                               n_top_genes = NULL,
                               max_iter = 10L,
                               fit_time = TRUE,
                               fit_scaling = TRUE,
                               fit_steady_states = TRUE,
                               fit_connected_states = FALSE,
                               t_max = 20,
                               eps = 1e-4,
                               n_jobs = 1L) {
  
  .vmessage("Recovering dynamics")
  
  # Get expression data
  if (!.layer_exists(seurat_obj, "Ms") || !.layer_exists(seurat_obj, "Mu")) {
    stop("Moments not found. Run compute_moments first.", call. = FALSE)
  }
  
  Ms <- make_dense(.get_layer_data(seurat_obj, "Ms"))
  Mu <- make_dense(.get_layer_data(seurat_obj, "Mu"))
  
  # Get genes to fit
  if (is.character(var_names) && var_names == "velocity_genes") {
    genes <- velocity_genes(seurat_obj)
  } else if (is.numeric(var_names)) {
    genes <- velocity_genes(seurat_obj, n_top_genes = var_names)
  } else {
    genes <- var_names
  }
  
  # Limit to top genes if specified
  if (!is.null(n_top_genes) && length(genes) > n_top_genes) {
    genes <- genes[seq_len(n_top_genes)]
  }
  
  genes <- intersect(genes, rownames(Ms))
  n_genes <- length(genes)
  n_cells <- ncol(Ms)
  
  .vmessage("Fitting dynamics for ", n_genes, " genes")
  
  # Initialize parameter storage
  params <- data.frame(
    row.names = rownames(Ms),
    fit_alpha = NA_real_,
    fit_beta = NA_real_,
    fit_gamma = NA_real_,
    fit_t_ = NA_real_,
    fit_scaling = NA_real_,
    fit_alignment = NA_real_,
    fit_likelihood = NA_real_,
    fit_pval_steady = NA_real_,
    fit_pval_dyn = NA_real_
  )
  
  # Storage for latent time
  latent_time <- matrix(NA_real_, nrow = n_genes, ncol = n_cells)
  rownames(latent_time) <- rownames(Ms)
  colnames(latent_time) <- colnames(Ms)
  
  # Fit each gene
  fit_single_gene <- function(gene) {
    u <- Mu[gene, ]
    s <- Ms[gene, ]
    
    # Fit dynamics
    result <- fit_dynamics_gene(u, s, 
                                 max_iter = max_iter,
                                 fit_time = fit_time,
                                 fit_scaling = fit_scaling,
                                 t_max = t_max,
                                 eps = eps)
    
    return(result)
  }
  
  # Run fitting
  results <- apply_with_progress(
    genes,
    fit_single_gene,
    n_jobs = n_jobs,
    show_progress = TRUE
  )
  
  # Extract results
  for (i in seq_along(genes)) {
    gene <- genes[i]
    if (!is.null(results[[i]]) && !is.null(results[[i]]$alpha)) {
      params[gene, "fit_alpha"] <- results[[i]]$alpha
      params[gene, "fit_beta"] <- results[[i]]$beta
      params[gene, "fit_gamma"] <- results[[i]]$gamma
      params[gene, "fit_t_"] <- results[[i]]$t_
      params[gene, "fit_scaling"] <- results[[i]]$scaling
      params[gene, "fit_likelihood"] <- results[[i]]$likelihood
      
      if (!is.null(results[[i]]$tau)) {
        latent_time[gene, ] <- results[[i]]$tau
      }
    }
  }
  
  # Store parameters
  seurat_obj <- store_fit_params(seurat_obj, params)
  seurat_obj@misc$latent_time <- latent_time
  seurat_obj@misc$dynamics_recovered <- TRUE
  
  n_fit <- sum(!is.na(params$fit_alpha))
  .vmessage("Recovered dynamics for ", n_fit, " genes")
  
  return(seurat_obj)
}

#' Fit Dynamics for Single Gene
#'
#' @description Fit splicing dynamics for a single gene using EM algorithm.
#'
#' @param u Unspliced values
#' @param s Spliced values
#' @param max_iter Maximum iterations
#' @param fit_time Whether to fit cell times
#' @param fit_scaling Whether to fit scaling
#' @param t_max Maximum time
#' @param eps Convergence threshold
#' @return List with fitted parameters
#' @keywords internal
fit_dynamics_gene <- function(u, s, max_iter = 10, fit_time = TRUE,
                               fit_scaling = TRUE, t_max = 20, eps = 1e-4) {
  
  n <- length(u)
  
  # Remove NAs and zeros
  valid <- !is.na(u) & !is.na(s) & (u > 0 | s > 0)
  if (sum(valid) < 20) {
    return(NULL)
  }
  
  u <- u[valid]
  s <- s[valid]
  n_valid <- length(u)
  
  # Initialize parameters from steady-state
  gamma_init <- sum(u * s) / sum(s^2)
  if (gamma_init <= 0 || !is.finite(gamma_init)) gamma_init <- 0.1
  
  beta_init <- gamma_init
  alpha_init <- gamma_init * max(s) / 2
  
  # Initial time estimate
  t_init <- t_max / 2
  scaling <- max(s) / max(u)
  if (!is.finite(scaling) || scaling <= 0) scaling <- 1
  
  # Scale u
  u_scaled <- u * scaling
  
  # EM iterations
  alpha <- alpha_init
  beta <- beta_init
  gamma <- gamma_init
  t_ <- t_init
  
  prev_likelihood <- -Inf
  
  for (iter in seq_len(max_iter)) {
    # E-step: Assign cells to time points
    if (fit_time) {
      # Estimate tau for each cell
      tau <- tau_inv_u_cpp(u_scaled, 0, alpha, beta)
      tau <- pmax(0, pmin(tau, t_max))
    } else {
      tau <- rep(t_ / 2, n_valid)
    }
    
    # Get predicted values
    dynamics <- mRNA_dynamics_cpp(tau, 0, 0, alpha, beta, gamma)
    u_pred <- dynamics$u
    s_pred <- dynamics$s
    
    # Compute likelihood (MSE-based)
    mse <- compute_mse_cpp(u_scaled, s, u_pred, s_pred)
    likelihood <- -mse
    
    # Check convergence
    if (abs(likelihood - prev_likelihood) < eps) {
      break
    }
    prev_likelihood <- likelihood
    
    # M-step: Update parameters
    
    # Update gamma: gamma = sum(beta*u*s) / sum(s^2)
    ss_sum <- sum(s^2)
    if (ss_sum > 0) {
      gamma_new <- beta * sum(u_scaled * s) / ss_sum
      if (is.finite(gamma_new) && gamma_new > 0) {
        gamma <- gamma_new
      }
    }
    
    # Update alpha and beta jointly
    # Use simplified update based on steady-state levels
    u_max <- max(u_scaled)
    s_max <- max(s)
    
    if (u_max > 0 && s_max > 0) {
      # At steady state: u_ss = alpha/beta, s_ss = alpha/gamma
      alpha_new <- gamma * s_max * 0.8  # Factor for non-steady cells
      beta_new <- alpha_new / u_max * 1.2
      
      if (is.finite(alpha_new) && alpha_new > 0) alpha <- alpha_new
      if (is.finite(beta_new) && beta_new > 0) beta <- beta_new
    }
    
    # Update t_ (switching time)
    if (fit_time) {
      # Estimate t_ as time when transcription switches off
      t_ <- stats::quantile(tau[tau > 0], 0.9, na.rm = TRUE)
      if (!is.finite(t_)) t_ <- t_max / 2
    }
  }
  
  # Final likelihood
  if (fit_time) {
    tau <- tau_inv_u_cpp(u_scaled, 0, alpha, beta)
    tau <- pmax(0, pmin(tau, t_max))
  } else {
    tau <- rep(t_ / 2, n_valid)
  }
  
  dynamics <- mRNA_dynamics_cpp(tau, 0, 0, alpha, beta, gamma)
  mse <- compute_mse_cpp(u_scaled, s, dynamics$u, dynamics$s)
  likelihood <- exp(-mse)  # Convert to probability-like scale
  
  # Rescale parameters
  if (fit_scaling) {
    alpha <- alpha / scaling
    beta <- beta
    gamma <- gamma
  }
  
  # Expand tau back to full array
  tau_full <- rep(NA_real_, length(valid))
  tau_full[valid] <- tau
  
  list(
    alpha = alpha,
    beta = beta,
    gamma = gamma,
    t_ = t_,
    scaling = scaling,
    likelihood = likelihood,
    tau = tau_full
  )
}

#' Recover Latent Time
#'
#' @description Compute latent time for each cell based on fitted dynamics.
#'
#' @param seurat_obj A Seurat object with recovered dynamics
#' @param vkey Velocity key
#' @param min_likelihood Minimum likelihood threshold
#' @param use_top_genes Number of top genes to use for time computation
#'
#' @return Modified Seurat object with latent_time in metadata
#'
#' @export
recover_latent_time <- function(seurat_obj,
                                  vkey = "velocity",
                                  min_likelihood = 0.1,
                                  use_top_genes = 100) {
  
  .vmessage("Recovering latent time")
  
  if (is.null(seurat_obj@misc$dynamics_recovered)) {
    stop("Dynamics not recovered. Run recover_dynamics first.", call. = FALSE)
  }
  
  # Get stored latent times per gene
  latent_time <- seurat_obj@misc$latent_time
  
  if (is.null(latent_time)) {
    stop("Latent time matrix not found", call. = FALSE)
  }
  
  # Get fit parameters
  fit_params <- get_fit_params(seurat_obj)
  
  # Select genes based on likelihood
  if ("fit_likelihood" %in% colnames(fit_params)) {
    likelihood <- fit_params$fit_likelihood
    good_genes <- !is.na(likelihood) & (likelihood >= min_likelihood)
  } else {
    good_genes <- !is.na(fit_params$fit_alpha)
  }
  
  good_gene_names <- rownames(fit_params)[good_genes]
  good_gene_names <- intersect(good_gene_names, rownames(latent_time))
  
  # Limit to top genes
  if (length(good_gene_names) > use_top_genes) {
    gene_scores <- fit_params[good_gene_names, "fit_likelihood"]
    top_idx <- order(gene_scores, decreasing = TRUE)[seq_len(use_top_genes)]
    good_gene_names <- good_gene_names[top_idx]
  }
  
  .vmessage("Using ", length(good_gene_names), " genes for latent time")
  
  # Compute consensus latent time
  latent_time_subset <- latent_time[good_gene_names, , drop = FALSE]
  
  # Median across genes (robust to outliers)
  cell_time <- apply(latent_time_subset, 2, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    stats::median(x)
  })
  
  # Scale to [0, 1]
  cell_time <- (cell_time - min(cell_time, na.rm = TRUE)) / 
               (max(cell_time, na.rm = TRUE) - min(cell_time, na.rm = TRUE))
  
  # Store
  seurat_obj@meta.data$latent_time <- cell_time
  
  .vmessage("Added latent_time to metadata")
  
  return(seurat_obj)
}

#' Differential Kinetic Test
#'
#' @description Test for differential kinetics between cell groups.
#'
#' @param seurat_obj A Seurat object with recovered dynamics
#' @param groupby Column name for grouping
#' @param groups Groups to compare (default: all pairs)
#' @param n_jobs Number of parallel jobs
#'
#' @return Data frame with differential kinetics results
#'
#' @export
differential_kinetic_test <- function(seurat_obj,
                                        groupby,
                                        groups = NULL,
                                        n_jobs = 1L) {
  
  .vmessage("Testing differential kinetics")
  
  if (!groupby %in% colnames(seurat_obj@meta.data)) {
    stop("Column '", groupby, "' not found in metadata", call. = FALSE)
  }
  
  # Get groups
  cell_groups <- seurat_obj@meta.data[[groupby]]
  unique_groups <- unique(cell_groups)
  
  if (!is.null(groups)) {
    groups <- intersect(groups, unique_groups)
  } else {
    groups <- unique_groups
  }
  
  if (length(groups) < 2) {
    stop("Need at least 2 groups for comparison", call. = FALSE)
  }
  
  # Get parameters
  fit_params <- get_fit_params(seurat_obj)
  genes <- rownames(fit_params)
  
  # For each gene, test if gamma differs between groups
  results <- data.frame(
    gene = genes,
    pvalue = NA_real_,
    log2fc = NA_real_,
    group1_gamma = NA_real_,
    group2_gamma = NA_real_
  )
  
  # Get moments
  Ms <- make_dense(.get_layer_data(seurat_obj, "Ms"))
  Mu <- make_dense(.get_layer_data(seurat_obj, "Mu"))
  
  for (i in seq_along(genes)) {
    gene <- genes[i]
    s <- Ms[gene, ]
    u <- Mu[gene, ]
    
    # Fit gamma per group
    gammas <- numeric(length(groups))
    names(gammas) <- groups
    
    for (g in seq_along(groups)) {
      mask <- cell_groups == groups[g]
      sg <- s[mask]
      ug <- u[mask]
      
      ss <- sum(sg^2)
      if (ss > 0) {
        gammas[g] <- sum(sg * ug) / ss
      }
    }
    
    # Simple t-test (comparing gamma estimates would require bootstrap)
    # Here we compare expression-velocity relationship
    if (length(groups) == 2) {
      mask1 <- cell_groups == groups[1]
      mask2 <- cell_groups == groups[2]
      
      # Velocity correlation as measure
      v1 <- u[mask1] - gammas[1] * s[mask1]
      v2 <- u[mask2] - gammas[2] * s[mask2]
      
      # Test difference in velocity
      if (length(v1) > 5 && length(v2) > 5) {
        test <- tryCatch({
          stats::t.test(v1, v2)
        }, error = function(e) NULL)
        
        if (!is.null(test)) {
          results$pvalue[i] <- test$p.value
        }
      }
      
      results$group1_gamma[i] <- gammas[1]
      results$group2_gamma[i] <- gammas[2]
      
      if (gammas[1] > 0 && gammas[2] > 0) {
        results$log2fc[i] <- log2(gammas[2] / gammas[1])
      }
    }
  }
  
  # Adjust p-values
  results$padj <- stats::p.adjust(results$pvalue, method = "BH")
  
  # Sort by significance
  results <- results[order(results$pvalue), ]
  
  return(results)
}

#' @title Main Velocity Function
#' @description Main interface for RNA velocity estimation.
#' @name velocity
NULL

#' Estimate RNA Velocity
#'
#' @description Main function for estimating RNA velocity. Supports steady-state,
#' stochastic, and dynamical models.
#'
#' @param seurat_obj A Seurat object with computed moments
#' @param mode Velocity mode: "steady_state" (default), "stochastic", or "dynamical"
#' @param fit_offset Whether to fit offset (default: FALSE for deterministic models)
#' @param min_r2 Minimum R-squared to consider a gene well-fit (default: 0.01)
#' @param min_likelihood Minimum likelihood for dynamical model (default: 0.001)
#' @param n_jobs Number of parallel jobs (default: 1)
#' @param copy Return copy instead of modifying in place (default: FALSE)
#'
#' @return Modified Seurat object with velocity estimates
#'
#' @export
#' @examples
#' \dontrun{
#' # Steady-state model
#' seurat_obj <- velocity(seurat_obj, mode = "steady_state")
#'
#' # Stochastic model
#' seurat_obj <- velocity(seurat_obj, mode = "stochastic")
#'
#' # Dynamical model (requires recover_dynamics first)
#' seurat_obj <- recover_dynamics(seurat_obj)
#' seurat_obj <- velocity(seurat_obj, mode = "dynamical")
#' }
velocity <- function(seurat_obj,
                      mode = c("steady_state", "stochastic", "dynamical"),
                      fit_offset = FALSE,
                      min_r2 = 0.01,
                      min_likelihood = 0.001,
                      n_jobs = 1L,
                      copy = FALSE) {
  
  mode <- match.arg(mode)
  
  if (copy) {
    seurat_obj <- seurat_obj  # R creates copy on modification anyway
  }
  
  .vmessage("Estimating velocity using ", mode, " model")
  
  # Check prerequisites
  if (!.layer_exists(seurat_obj, "Ms") || !.layer_exists(seurat_obj, "Mu")) {
    .vmessage("Moments not found, computing...")
    seurat_obj <- compute_moments(seurat_obj)
  }
  
  # Dispatch to appropriate model
  if (mode == "steady_state") {
    seurat_obj <- fit_velocity_steady(seurat_obj,
                                       fit_intercept = fit_offset,
                                       min_r2 = min_r2,
                                       n_jobs = n_jobs)
  } else if (mode == "stochastic") {
    seurat_obj <- fit_velocity_stochastic(seurat_obj,
                                           fit_offset = fit_offset,
                                           min_r2 = min_r2,
                                           n_jobs = n_jobs)
  } else if (mode == "dynamical") {
    # For dynamical model, dynamics must be recovered first
    if (is.null(seurat_obj@misc$dynamics_recovered)) {
      stop("Dynamics not recovered. Run recover_dynamics first.", call. = FALSE)
    }
    seurat_obj <- fit_velocity_dynamical(seurat_obj,
                                          min_likelihood = min_likelihood,
                                          n_jobs = n_jobs)
  }
  
  .vmessage("Velocity estimation complete")
  
  return(seurat_obj)
}

#' Get Velocity Matrix
#'
#' @description Extract velocity matrix from Seurat object.
#'
#' @param seurat_obj A Seurat object with computed velocity
#' @param genes Genes to include (NULL for all)
#' @param remove_na Whether to remove genes with NA velocities
#'
#' @return Velocity matrix (genes x cells)
#'
#' @export
get_velocity <- function(seurat_obj, genes = NULL, remove_na = TRUE) {
  
  if (!.layer_exists(seurat_obj, "velocity")) {
    stop("Velocity not found. Run velocity() first.", call. = FALSE)
  }
  
  velocity <- .get_layer_data(seurat_obj, "velocity")
  
  # Subset genes
  if (!is.null(genes)) {
    genes <- intersect(genes, rownames(velocity))
    velocity <- velocity[genes, , drop = FALSE]
  }
  
  # Remove NA genes
  if (remove_na) {
    na_genes <- Matrix::rowSums(is.na(velocity)) == ncol(velocity)
    velocity <- velocity[!na_genes, , drop = FALSE]
  }
  
  return(velocity)
}

#' Get Velocity Mode
#'
#' @description Get the mode used for velocity estimation.
#'
#' @param seurat_obj A Seurat object
#' @return Character string: "steady_state", "stochastic", or "dynamical"
#'
#' @export
get_velocity_mode <- function(seurat_obj) {
  mode <- seurat_obj@misc$velocity_mode
  if (is.null(mode)) {
    return(NA_character_)
  }
  return(mode)
}

#' Fit Dynamical Velocity
#'
#' @description Calculate velocity from fitted dynamical parameters.
#'
#' @param seurat_obj A Seurat object with recovered dynamics
#' @param min_likelihood Minimum likelihood (default: 0.001)
#' @param n_jobs Number of parallel jobs
#'
#' @return Modified Seurat object with velocity
#'
#' @keywords internal
fit_velocity_dynamical <- function(seurat_obj,
                                    min_likelihood = 0.001,
                                    n_jobs = 1L) {
  
  .vmessage("Computing velocity from dynamical model...")
  
  # Get fit parameters
  fit_params <- get_fit_params(seurat_obj)
  
  if (is.null(fit_params) || !"fit_alpha" %in% colnames(fit_params)) {
    stop("Dynamical parameters not found. Run recover_dynamics first.", call. = FALSE)
  }
  
  # Get moments
  Ms <- make_dense(.get_layer_data(seurat_obj, "Ms"))
  Mu <- make_dense(.get_layer_data(seurat_obj, "Mu"))
  
  genes <- rownames(Ms)
  n_genes <- length(genes)
  n_cells <- ncol(Ms)
  
  # Get parameters
  alpha <- fit_params$fit_alpha
  beta <- fit_params$fit_beta
  gamma <- fit_params$fit_gamma
  likelihood <- fit_params$fit_likelihood
  
  # Determine velocity genes
  velocity_genes <- !is.na(alpha) & !is.na(gamma) &
                     gamma > 0 &
                     (is.na(likelihood) | likelihood >= min_likelihood)
  
  # Calculate velocity: v = alpha/gamma * beta*u - gamma*s
  # Simplified: v = beta * u - gamma * s  (when at steady state)
  velocity <- matrix(NA_real_, nrow = n_genes, ncol = n_cells)
  rownames(velocity) <- genes
  colnames(velocity) <- colnames(Ms)
  
  for (i in which(velocity_genes)) {
    u <- Mu[i, ]
    s <- Ms[i, ]
    
    # Get gene-specific time if available
    if ("fit_t_" %in% colnames(fit_params)) {
      t_ <- fit_params$fit_t_[i]
    } else {
      t_ <- NA
    }
    
    # Compute velocity
    # For on-state: v = alpha - beta*u = beta*(alpha/beta - u)
    # For off-state: v = -beta*u
    # Simplified steady-state approximation: v = beta*u - gamma*s
    
    b <- beta[i]
    g <- gamma[i]
    
    if (!is.na(b) && !is.na(g) && g > 0) {
      velocity[i, ] <- b * u - g * s
    }
  }
  
  # Store velocity
  velocity <- as(velocity, "CsparseMatrix")
  seurat_obj <- .set_layer_data(seurat_obj, "velocity", velocity)
  
  # Update velocity genes
  fit_params$velocity_genes <- velocity_genes
  seurat_obj <- store_fit_params(seurat_obj, fit_params)
  
  seurat_obj@misc$velocity_mode <- "dynamical"
  
  n_fit <- sum(velocity_genes, na.rm = TRUE)
  .vmessage("Computed velocity for ", n_fit, " genes")
  
  return(seurat_obj)
}

#' Rank Velocity Genes
#'
#' @description Rank genes by their velocity fit quality or other metrics.
#'
#' @param seurat_obj A Seurat object with fitted velocity
#' @param groupby Column name to group cells
#' @param groups Groups to compare (default: all)
#' @param n_genes Number of top genes to return
#' @param min_r2 Minimum R-squared
#' @param method Ranking method: "r2" (default), "velocity_var", "likelihood"
#'
#' @return Data frame with ranked genes
#'
#' @export
rank_velocity_genes <- function(seurat_obj,
                                 groupby = NULL,
                                 groups = NULL,
                                 n_genes = 20,
                                 min_r2 = 0.01,
                                 method = "r2") {
  
  # Get fit parameters
  fit_params <- get_fit_params(seurat_obj)
  
  if (is.null(fit_params)) {
    stop("Fit parameters not found. Run velocity() first.", call. = FALSE)
  }
  
  genes <- rownames(fit_params)
  
  # Filter by minimum criteria
  mask <- rep(TRUE, nrow(fit_params))
  
  if ("fit_r2" %in% colnames(fit_params)) {
    mask <- mask & !is.na(fit_params$fit_r2) & (fit_params$fit_r2 >= min_r2)
  }
  
  if ("velocity_genes" %in% colnames(fit_params)) {
    mask <- mask & fit_params$velocity_genes
  }
  
  # Get ranking metric
  if (method == "r2" && "fit_r2" %in% colnames(fit_params)) {
    score <- fit_params$fit_r2
  } else if (method == "likelihood" && "fit_likelihood" %in% colnames(fit_params)) {
    score <- fit_params$fit_likelihood
  } else {
    # Default to r2
    if ("fit_r2" %in% colnames(fit_params)) {
      score <- fit_params$fit_r2
    } else {
      score <- rep(0, nrow(fit_params))
    }
  }
  
  # Apply mask
  score[!mask] <- NA
  
  # Rank
  rank_order <- order(score, decreasing = TRUE, na.last = TRUE)
  
  # Select top genes
  n_genes <- min(n_genes, sum(mask))
  top_indices <- rank_order[seq_len(n_genes)]
  
  # Create result data frame
  result <- data.frame(
    gene = genes[top_indices],
    score = score[top_indices],
    row.names = NULL
  )
  
  # Add other metrics if available
  for (col in c("fit_gamma", "fit_r2", "fit_likelihood")) {
    if (col %in% colnames(fit_params)) {
      result[[col]] <- fit_params[[col]][top_indices]
    }
  }
  
  return(result)
}

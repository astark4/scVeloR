#' @title Latent Time Estimation
#' @description Compute gene-shared latent time based on dynamical model parameters.
#' Gene-specific latent timepoints are coupled to a universal gene-shared latent time,
#' which represents the cell's internal clock based only on transcriptional dynamics.
#' @name latent_time
NULL

#' Compute Latent Time
#'
#' @description Compute gene-shared latent time from dynamical model results.
#' The latent time represents the cell's internal clock and is based only on
#' its transcriptional dynamics.
#'
#' @param object Seurat object with dynamics fitted
#' @param min_likelihood Minimum likelihood for gene inclusion
#' @param min_confidence Confidence threshold for local coherence
#' @param root_key Key in object metadata for root cells
#' @param end_key Key in object metadata for end cells
#' @param t_max Maximum time scale (default NULL, auto-determined)
#' @param verbose Print progress
#'
#' @return Seurat object with latent_time in metadata
#' @export
compute_latent_time <- function(object,
                                 min_likelihood = 0.1,
                                 min_confidence = 0.75,
                                 root_key = NULL,
                                 end_key = NULL,
                                 t_max = NULL,
                                 verbose = TRUE) {
  
  if (is.null(object@misc$scVeloR$dynamics)) {
    stop("Run recover_dynamics() first")
  }
  
  dynamics <- object@misc$scVeloR$dynamics
  gene_params <- dynamics$gene_params
  fit_t <- dynamics$fit_t
  
  if (verbose) {
    message("Computing latent time...")
  }
  
  # Filter genes by likelihood
  valid_genes <- !is.na(gene_params$likelihood)
  if (!is.null(min_likelihood)) {
    valid_genes <- valid_genes & (gene_params$likelihood >= min_likelihood)
  }
  
  if (sum(valid_genes) < 3) {
    stop("Too few valid genes for latent time computation")
  }
  
  # Get time matrix for valid genes
  t <- fit_t[, valid_genes, drop = FALSE]
  
  # Remove any columns with all NAs
  valid_cols <- colSums(!is.na(t)) > 0
  t <- t[, valid_cols, drop = FALSE]
  
  # Get connectivity matrix
  conn <- NULL
  if (!is.null(object@misc$scVeloR$neighbors$connectivities)) {
    conn <- object@misc$scVeloR$neighbors$connectivities
  }
  
  # Determine root cells
  root_cells <- NULL
  
  # Sum of times as initial root indicator (cells with lowest sum are early)
  t_sum <- rowSums(t, na.rm = TRUE)
  
  if (!is.null(root_key)) {
    if (root_key %in% colnames(object@meta.data)) {
      root_indicator <- object@meta.data[[root_key]]
      
      if (is.numeric(root_indicator)) {
        # Use cells with highest values
        root_cells <- which(root_indicator > quantile(root_indicator, 0.95, na.rm = TRUE))
      } else if (is.logical(root_indicator)) {
        root_cells <- which(root_indicator)
      } else {
        # Assume factor/character, use first level
        root_cells <- which(root_indicator == levels(factor(root_indicator))[1])
      }
    }
  }
  
  # If no root key provided or no cells found, use velocity pseudotime
  if (is.null(root_cells) || length(root_cells) == 0) {
    # Use cells with lowest t_sum as roots
    root_cells <- order(t_sum)[1:min(10, ceiling(length(t_sum) * 0.05))]
  }
  
  if (verbose) {
    message(sprintf("  Using %d root cells", length(root_cells)))
  }
  
  # Compute latent time using multiple roots
  n_roots <- min(length(root_cells), 4)
  latent_time_list <- matrix(NA, n_roots, nrow(t))
  
  for (i in seq_len(n_roots)) {
    root <- root_cells[i]
    
    # Root the time matrix
    rooted <- root_time(t, root)
    t_rooted <- rooted$t
    
    # Compute shared time
    latent_time_list[i, ] <- compute_shared_time(t_rooted)
  }
  
  # Average across roots
  latent_time <- colMeans(latent_time_list, na.rm = TRUE)
  
  # Scale to [0, 1]
  latent_time <- scale_to_range(latent_time)
  
  # Incorporate end cell information if provided
  if (!is.null(end_key) && end_key %in% colnames(object@meta.data)) {
    end_indicator <- object@meta.data[[end_key]]
    
    end_cells <- NULL
    if (is.numeric(end_indicator)) {
      end_cells <- which(end_indicator > quantile(end_indicator, 0.95, na.rm = TRUE))
    } else if (is.logical(end_indicator)) {
      end_cells <- which(end_indicator)
    }
    
    if (!is.null(end_cells) && length(end_cells) > 0) {
      # Compute time from end cells
      n_ends <- min(length(end_cells), 4)
      latent_time_end <- matrix(NA, n_ends, nrow(t))
      
      for (i in seq_len(n_ends)) {
        end <- end_cells[i]
        rooted <- root_time(t, end)
        latent_time_end[i, ] <- 1 - compute_shared_time(rooted$t)
      }
      
      latent_time_end_avg <- colMeans(latent_time_end, na.rm = TRUE)
      latent_time_end_avg <- scale_to_range(latent_time_end_avg)
      
      # Combine with weight
      latent_time <- scale_to_range(latent_time + 0.2 * latent_time_end_avg)
    }
  }
  
  # Smooth with connectivity and local coherence filtering
  if (!is.null(conn)) {
    tl <- latent_time
    tc <- as.vector(conn %*% latent_time)
    
    # Local coherence confidence
    z <- sum(tl * tc) / sum(tc * tc)
    tl_conf <- (1 - abs(tl / max(tl, na.rm = TRUE) - tc * z / max(tl, na.rm = TRUE)))^2
    
    idx_low_confidence <- tl_conf < min_confidence
    
    # Smooth low confidence cells
    conn_smooth <- conn
    conn_smooth[, idx_low_confidence] <- 0
    conn_smooth <- drop0(conn_smooth)
    
    latent_time <- as.vector(conn_smooth %*% latent_time)
  }
  
  # Final scaling
  latent_time <- scale_to_range(latent_time)
  
  # Apply t_max if specified
  if (!is.null(t_max)) {
    latent_time <- latent_time * t_max
  }
  
  # Store in metadata
  object@meta.data$latent_time <- latent_time
  
  if (verbose) {
    message("Done. Latent time stored in object@meta.data$latent_time")
  }
  
  object
}

#' Root Time Matrix
#'
#' @description Reroot time matrix based on a specified root cell.
#'
#' @param t Time matrix (cells x genes)
#' @param root Index of root cell
#'
#' @return List with rooted time matrix and switch times
#' @keywords internal
root_time <- function(t, root = NULL) {
  
  # Handle NaN columns
  valid_cols <- colSums(is.na(t)) < nrow(t)
  t_valid <- t[, valid_cols, drop = FALSE]
  
  if (is.null(root)) {
    t_root <- rep(0, ncol(t_valid))
  } else {
    t_root <- t_valid[root, ]
  }
  
  n_cells <- nrow(t_valid)
  n_genes <- ncol(t_valid)
  
  t_rooted <- matrix(NA, n_cells, n_genes)
  t_switch <- numeric(n_genes)
  
  for (j in seq_len(n_genes)) {
    t_col <- t_valid[, j]
    t_r <- t_root[j]
    
    if (is.na(t_r)) {
      t_rooted[, j] <- t_col
      next
    }
    
    # Cells after root (o = 1)
    o <- as.integer(t_col >= t_r)
    
    # Time after root
    t_after <- (t_col - t_r) * o
    
    # Origin time (max time from root)
    t_origin <- max(t_after, na.rm = TRUE)
    if (!is.finite(t_origin)) t_origin <- 0
    
    # Time before root (shifted to positive)
    t_before <- (t_col + t_origin) * (1 - o)
    
    # Switch time
    t_switch[j] <- min(t_before[t_before > 0], na.rm = TRUE)
    if (!is.finite(t_switch[j])) t_switch[j] <- 0
    
    # Rooted time
    t_rooted[, j] <- t_after + t_before
  }
  
  list(t = t_rooted, t_switch = t_switch)
}

#' Compute Shared Time
#'
#' @description Compute gene-shared time by selecting optimal percentiles.
#' The algorithm finds percentiles that produce the most linear distribution
#' of cell times.
#'
#' @param t Rooted time matrix (cells x genes)
#' @param percentiles Vector of percentiles to try
#' @param normalize Whether to normalize to [0, 1]
#'
#' @return Vector of shared times
#' @keywords internal
compute_shared_time <- function(t, 
                                 percentiles = c(15, 20, 25, 30, 35),
                                 normalize = TRUE) {
  
  # Remove NaN columns
  valid_cols <- colSums(is.na(t)) < nrow(t)
  t <- t[, valid_cols, drop = FALSE]
  
  # Shift to start at 0
  t <- t - min(t, na.rm = TRUE)
  
  n_cells <- nrow(t)
  n_perc <- length(percentiles)
  
  # Compute percentile times for each gene
  tx_list <- matrix(NA, n_perc, n_cells)
  
  for (i in seq_len(n_perc)) {
    perc <- percentiles[i]
    
    # For each cell, get the percentile across genes
    for (k in seq_len(n_cells)) {
      tx_list[i, k] <- quantile(t[k, ], perc / 100, na.rm = TRUE)
    }
  }
  
  # Normalize each percentile series
  tx_max <- apply(tx_list, 1, max, na.rm = TRUE)
  tx_max[tx_max == 0] <- 1
  tx_list <- tx_list / tx_max
  
  # Find percentiles with most linear distribution
  mse <- numeric(n_perc)
  
  for (i in seq_len(n_perc)) {
    tx <- tx_list[i, ]
    tx_sorted <- sort(tx, na.last = TRUE)
    linx <- seq(0, 1, length.out = length(tx_sorted))
    mse[i] <- sum((tx_sorted - linx)^2, na.rm = TRUE)
  }
  
  # Use best 2 percentiles
  best_idx <- order(mse)[1:min(2, n_perc)]
  
  # Average the best percentiles
  t_shared <- colMeans(tx_list[best_idx, , drop = FALSE], na.rm = TRUE)
  
  if (normalize) {
    t_shared <- t_shared / max(t_shared, na.rm = TRUE)
  }
  
  t_shared
}

#' Scale to Range
#'
#' @description Scale vector to [0, 1] range.
#'
#' @param x Numeric vector
#'
#' @return Scaled vector
#' @keywords internal
scale_to_range <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  
  if (x_max == x_min) {
    return(rep(0.5, length(x)))
  }
  
  (x - x_min) / (x_max - x_min)
}

#' Terminal States Detection
#'
#' @description Identify terminal states (root and end cells) using velocity graph.
#'
#' @param object Seurat object with velocity graph computed
#' @param n_states Number of terminal states to detect (default 2)
#' @param min_cells Minimum cells for a terminal state
#' @param verbose Print progress
#'
#' @return Seurat object with terminal state annotations
#' @export
terminal_states <- function(object,
                            n_states = 2,
                            min_cells = 10,
                            verbose = TRUE) {
  
  # Get transition matrix
  T <- NULL
  
  if (!is.null(object@misc$scVeloR$velocity$transition_matrix)) {
    T <- object@misc$scVeloR$velocity$transition_matrix
  } else if (!is.null(object@misc$scVeloR$velocity_graph$transition_matrix)) {
    T <- object@misc$scVeloR$velocity_graph$transition_matrix
  }
  
  if (is.null(T)) {
    stop("Compute velocity graph first using velocity_graph()")
  }
  
  if (verbose) {
    message("Detecting terminal states...")
  }
  
  n_cells <- nrow(T)
  
  # Compute stationary distribution for roots (reverse transitions)
  T_rev <- t(T)
  T_rev <- T_rev / rowSums(T_rev)
  T_rev[!is.finite(T_rev)] <- 0
  
  # Power iteration for stationary distribution
  root_probs <- rep(1 / n_cells, n_cells)
  for (i in 1:100) {
    root_probs_new <- as.vector(T_rev %*% root_probs)
    root_probs_new <- root_probs_new / sum(root_probs_new)
    
    if (max(abs(root_probs_new - root_probs)) < 1e-6) break
    root_probs <- root_probs_new
  }
  
  # Compute end probabilities (forward transitions)
  T_fwd <- T / rowSums(T)
  T_fwd[!is.finite(T_fwd)] <- 0
  
  end_probs <- rep(1 / n_cells, n_cells)
  for (i in 1:100) {
    end_probs_new <- as.vector(T_fwd %*% end_probs)
    end_probs_new <- end_probs_new / sum(end_probs_new)
    
    if (max(abs(end_probs_new - end_probs)) < 1e-6) break
    end_probs <- end_probs_new
  }
  
  # Store in metadata
  object@meta.data$root_cells <- root_probs
  object@meta.data$end_cells <- end_probs
  
  # Identify discrete terminal states
  root_threshold <- quantile(root_probs, 1 - min_cells / n_cells)
  end_threshold <- quantile(end_probs, 1 - min_cells / n_cells)
  
  object@meta.data$is_root <- root_probs >= root_threshold
  object@meta.data$is_end <- end_probs >= end_threshold
  
  if (verbose) {
    n_roots <- sum(object@meta.data$is_root)
    n_ends <- sum(object@meta.data$is_end)
    message(sprintf("  Identified %d root cells and %d end cells", n_roots, n_ends))
  }
  
  object
}

#' Velocity Pseudotime
#'
#' @description Compute pseudotime based on velocity-directed diffusion.
#'
#' @param object Seurat object
#' @param root Root cell index or NULL for automatic detection
#' @param n_dcs Number of diffusion components
#' @param verbose Print progress
#'
#' @return Seurat object with velocity_pseudotime in metadata
#' @export
velocity_pseudotime <- function(object,
                                 root = NULL,
                                 n_dcs = 10,
                                 verbose = TRUE) {
  
  # Get transition matrix
  T <- NULL
  
  if (!is.null(object@misc$scVeloR$velocity$transition_matrix)) {
    T <- object@misc$scVeloR$velocity$transition_matrix
  } else if (!is.null(object@misc$scVeloR$velocity_graph$transition_matrix)) {
    T <- object@misc$scVeloR$velocity_graph$transition_matrix
  }
  
  if (is.null(T)) {
    stop("Compute velocity graph first using velocity_graph()")
  }
  
  if (verbose) {
    message("Computing velocity pseudotime...")
  }
  
  n_cells <- nrow(T)
  
  # Determine root if not provided
  if (is.null(root)) {
    if ("root_cells" %in% colnames(object@meta.data)) {
      root <- which.max(object@meta.data$root_cells)
    } else {
      # Use cell with most outgoing transitions
      out_degree <- rowSums(T > 0)
      root <- which.max(out_degree)
    }
  }
  
  if (verbose) {
    message(sprintf("  Using cell %d as root", root))
  }
  
  # Compute diffusion distance from root
  # Symmetric normalized transition matrix
  T_sym <- (T + t(T)) / 2
  D <- rowSums(T_sym)
  D_inv_sqrt <- 1 / sqrt(D)
  D_inv_sqrt[!is.finite(D_inv_sqrt)] <- 0
  
  T_norm <- Diagonal(x = D_inv_sqrt) %*% T_sym %*% Diagonal(x = D_inv_sqrt)
  
  # Compute eigendecomposition
  n_components <- min(n_dcs + 1, n_cells - 1)
  
  eig <- tryCatch({
    RSpectra::eigs_sym(T_norm, k = n_components, which = "LM")
  }, error = function(e) {
    eigen(as.matrix(T_norm), symmetric = TRUE)
  })
  
  # Get diffusion components
  dcs <- eig$vectors[, 2:min(n_components, ncol(eig$vectors)), drop = FALSE]
  
  # Scale by eigenvalues
  eigenvalues <- eig$values[2:min(n_components, length(eig$values))]
  dcs <- dcs %*% diag(eigenvalues)
  
  # Compute pseudotime as distance from root in diffusion space
  root_dc <- dcs[root, ]
  
  pseudotime <- numeric(n_cells)
  for (i in seq_len(n_cells)) {
    pseudotime[i] <- sqrt(sum((dcs[i, ] - root_dc)^2))
  }
  
  # Scale to [0, 1]
  pseudotime <- scale_to_range(pseudotime)
  
  # Store in metadata
  object@meta.data$velocity_pseudotime <- pseudotime
  
  if (verbose) {
    message("Done. Velocity pseudotime stored in object@meta.data$velocity_pseudotime")
  }
  
  object
}

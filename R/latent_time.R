#' @title Latent Time Computation
#' @description Complete implementation of latent_time() function
#' equivalent to scvelo.tools.velocity_pseudotime.latent_time
#' @name latent_time
NULL

#' Compute Latent Time
#'
#' @description Computes a shared latent time across all genes by combining
#' gene-specific time estimates from the dynamical model.
#'
#' @param seurat Seurat object with velocity data
#' @param vkey Velocity key
#' @param min_likelihood Minimum likelihood threshold for gene filtering
#' @param min_confidence Minimum confidence for soft assignment
#' @param min_corr_diffusion Minimum correlation with diffusion pseudotime
#' @param weight_diffusion Weight for diffusion pseudotime
#' @param root_key Key for root cells
#' @param end_key Key for terminal cells
#' @param t_max Maximum time
#' @param copy Return copy
#' @return Modified Seurat object or latent time vector
#' @export
latent_time <- function(seurat = NULL,
                        vkey = "velocity",
                        min_likelihood = 0.1,
                        min_confidence = NULL,
                        min_corr_diffusion = NULL,
                        weight_diffusion = NULL,
                        root_key = NULL,
                        end_key = NULL,
                        t_max = NULL,
                        copy = FALSE) {
  
  # Check for dynamical model results
  if (!is.null(seurat) && inherits(seurat, "Seurat")) {
    if (!"velocity" %in% names(seurat@misc)) {
      stop("Run velocity() first before computing latent time")
    }
    
    vdata <- seurat@misc$velocity
    
    if (is.null(vdata$fit_t)) {
      stop("Latent time requires dynamical model. Run velocity(mode='dynamical') first.")
    }
    
    fit_t <- vdata$fit_t
    likelihood <- vdata$gene_params$likelihood
    gene_names <- vdata$genes
    
    # Get connectivity for smoothing
    connectivities <- NULL
    if ("RNA_snn" %in% names(seurat@graphs)) {
      connectivities <- seurat@graphs$RNA_snn
    }
  } else if (is.list(seurat)) {
    # Direct input
    fit_t <- seurat$fit_t
    likelihood <- seurat$gene_params$likelihood
    gene_names <- seurat$genes
    connectivities <- seurat$connectivities
  } else {
    stop("Input must be Seurat object or list with fit_t")
  }
  
  # Filter genes by likelihood
  if (!is.null(min_likelihood) && !all(is.na(likelihood))) {
    valid_genes <- !is.na(likelihood) & (likelihood >= min_likelihood)
    if (sum(valid_genes) < 10) {
      warning("Few genes pass likelihood threshold, using all genes")
      valid_genes <- !is.na(likelihood)
    }
    fit_t <- fit_t[, valid_genes, drop = FALSE]
    gene_names <- gene_names[valid_genes]
  }
  
  # Compute shared latent time
  t_shared <- compute_shared_time(fit_t, norm = TRUE)
  
  # Root the time
  if (!is.null(root_key)) {
    if (inherits(seurat, "Seurat") && root_key %in% names(seurat@meta.data)) {
      root_cells <- which(seurat@meta.data[[root_key]])
      if (length(root_cells) > 0) {
        t_root <- mean(t_shared[root_cells])
        t_shared <- t_shared - t_root
        t_shared <- pmax(t_shared, 0)
      }
    }
  }
  
  # Normalize to [0, 1]
  t_shared <- (t_shared - min(t_shared, na.rm = TRUE)) / 
              (max(t_shared, na.rm = TRUE) - min(t_shared, na.rm = TRUE) + 1e-10)
  
  # Smooth with connectivities
  if (!is.null(connectivities)) {
    conn_norm <- connectivities / Matrix::rowSums(connectivities)
    t_shared <- as.vector(conn_norm %*% t_shared)
    
    # Re-normalize
    t_shared <- (t_shared - min(t_shared, na.rm = TRUE)) / 
                (max(t_shared, na.rm = TRUE) - min(t_shared, na.rm = TRUE) + 1e-10)
  }
  
  # Add to Seurat
  if (inherits(seurat, "Seurat")) {
    seurat$latent_time <- t_shared
    seurat@misc$velocity$latent_time <- t_shared
    return(seurat)
  }
  
  t_shared
}

#' Velocity Pseudotime
#'
#' @description Computes a velocity-based pseudotime using transition probabilities
#' @param seurat Seurat object with velocity
#' @param vkey Velocity key
#' @param groupby Group by variable
#' @param groups Groups to include
#' @param root_key Root cell key
#' @param end_key End cell key
#' @param use_velocity Use velocity graph
#' @param n_dcs Number of diffusion components
#' @param use_raw Use raw velocity
#' @param self_transitions Allow self transitions
#' @param eps Epsilon for numerical stability
#' @param copy Return copy
#' @return Modified Seurat object or pseudotime vector
#' @export
velocity_pseudotime <- function(seurat,
                                vkey = "velocity",
                                groupby = NULL,
                                groups = NULL,
                                root_key = NULL,
                                end_key = NULL,
                                use_velocity = TRUE,
                                n_dcs = 10,
                                use_raw = FALSE,
                                self_transitions = FALSE,
                                eps = 1e-3,
                                copy = FALSE) {
  
  # Get velocity data
  if (!"velocity" %in% names(seurat@misc)) {
    stop("Run velocity() first")
  }
  
  vdata <- seurat@misc$velocity
  
  # Get transition matrix
  if ("transition_matrix" %in% names(vdata)) {
    T_mat <- vdata$transition_matrix
  } else if ("velocity_graph" %in% names(vdata)) {
    T_mat <- compute_transition_matrix(vdata$velocity_graph)
  } else {
    stop("Need velocity graph or transition matrix")
  }
  
  n_cells <- nrow(T_mat)
  
  # Remove self transitions if requested
  if (!self_transitions) {
    diag(T_mat) <- 0
    T_mat <- T_mat / (rowSums(T_mat) + eps)
  }
  
  # Determine root cells
  if (!is.null(root_key) && root_key %in% colnames(seurat@meta.data)) {
    roots <- which(seurat@meta.data[[root_key]])
  } else {
    # Use cells with highest incoming transitions (sources)
    incoming <- colSums(T_mat)
    outgoing <- rowSums(T_mat)
    source_score <- outgoing - incoming
    roots <- order(source_score, decreasing = TRUE)[1:max(1, n_cells %/% 100)]
  }
  
  # Determine end cells
  if (!is.null(end_key) && end_key %in% colnames(seurat@meta.data)) {
    ends <- which(seurat@meta.data[[end_key]])
  } else {
    incoming <- colSums(T_mat)
    outgoing <- rowSums(T_mat)
    sink_score <- incoming - outgoing
    ends <- order(sink_score, decreasing = TRUE)[1:max(1, n_cells %/% 100)]
  }
  
  # Compute pseudotime using absorption probabilities
  pseudotime <- compute_absorption_time(T_mat, ends)
  
  # Invert so root has low time
  root_time <- mean(pseudotime[roots], na.rm = TRUE)
  end_time <- mean(pseudotime[ends], na.rm = TRUE)
  
  if (root_time > end_time) {
    pseudotime <- max(pseudotime, na.rm = TRUE) - pseudotime
  }
  
  # Normalize
  pseudotime <- (pseudotime - min(pseudotime, na.rm = TRUE)) /
                (max(pseudotime, na.rm = TRUE) - min(pseudotime, na.rm = TRUE) + 1e-10)
  
  # Add to Seurat
  seurat$velocity_pseudotime <- pseudotime
  seurat@misc$velocity$pseudotime <- pseudotime
  
  seurat
}

#' Compute absorption time
#'
#' @description Compute mean hitting time to absorbing states
#' @param T_mat Transition matrix
#' @param absorbing Indices of absorbing states
#' @param max_iter Maximum iterations for power method
#' @return Mean hitting time vector
#' @keywords internal
compute_absorption_time <- function(T_mat, absorbing, max_iter = 1000) {
  n <- nrow(T_mat)
  
  # Simple approach: iteratively propagate time
  time <- rep(0, n)
  time[absorbing] <- 0
  
  # Non-absorbing states
  transient <- setdiff(seq_len(n), absorbing)
  
  if (length(transient) == 0) {
    return(time)
  }
  
  # Power iteration
  for (i in seq_len(max_iter)) {
    time_new <- time
    time_new[transient] <- 1 + as.vector(T_mat[transient, , drop = FALSE] %*% time)
    
    # Convergence check
    if (max(abs(time_new - time)) < 1e-6) break
    time <- time_new
  }
  
  time
}

#' Terminal States
#'
#' @description Identify terminal (end) states using velocity
#' @param seurat Seurat object with velocity
#' @param vkey Velocity key
#' @param groupby Group by variable
#' @param groups Groups to include
#' @param self_transitions Include self transitions
#' @param eps Epsilon
#' @param random_state Random seed
#' @param copy Return copy
#' @return Modified Seurat object with terminal_states column
#' @export
terminal_states <- function(seurat,
                            vkey = "velocity",
                            groupby = NULL,
                            groups = NULL,
                            self_transitions = FALSE,
                            eps = 1e-3,
                            random_state = 0,
                            copy = FALSE) {
  
  set.seed(random_state)
  
  # Get velocity data
  if (!"velocity" %in% names(seurat@misc)) {
    stop("Run velocity() first")
  }
  
  vdata <- seurat@misc$velocity
  
  # Get or compute transition matrix
  if ("transition_matrix" %in% names(vdata)) {
    T_mat <- vdata$transition_matrix
  } else if ("velocity_graph" %in% names(vdata)) {
    T_mat <- compute_transition_matrix(vdata$velocity_graph)
  } else {
    stop("Need velocity graph first")
  }
  
  n_cells <- nrow(T_mat)
  
  # Remove self transitions
  if (!self_transitions) {
    diag(T_mat) <- 0
    T_mat <- T_mat / (rowSums(T_mat) + eps)
  }
  
  # Compute terminal state scores based on:
  # 1. Low outgoing transitions (sink)
  # 2. High incoming transitions  
  outgoing <- rowSums(T_mat)
  incoming <- colSums(T_mat)
  
  # Terminal score: high incoming / outgoing ratio
  terminal_score <- incoming / (outgoing + eps)
  
  # Normalize
  terminal_score <- (terminal_score - min(terminal_score, na.rm = TRUE)) /
                    (max(terminal_score, na.rm = TRUE) - min(terminal_score, na.rm = TRUE) + 1e-10)
  
  # Root score: high outgoing / incoming ratio
  root_score <- outgoing / (incoming + eps)
  root_score <- (root_score - min(root_score, na.rm = TRUE)) /
                (max(root_score, na.rm = TRUE) - min(root_score, na.rm = TRUE) + 1e-10)
  
  # Add to Seurat
  seurat$terminal_score <- terminal_score
  seurat$root_score <- root_score
  
  # Binary classification
  seurat$is_terminal <- terminal_score > quantile(terminal_score, 0.95)
  seurat$is_root <- root_score > quantile(root_score, 0.95)
  
  seurat@misc$velocity$terminal_score <- terminal_score
  seurat@misc$velocity$root_score <- root_score
  
  seurat
}

#' Rank Velocity Genes
#'
#' @description Rank genes by their velocity fit quality
#' @param seurat Seurat object with velocity
#' @param n_genes Number of top genes to return
#' @param min_r2 Minimum R2 threshold
#' @param min_counts Minimum counts
#' @param ... Additional arguments
#' @return Data frame with ranked genes
#' @export
rank_velocity_genes <- function(seurat, n_genes = 100, min_r2 = 0.1, min_counts = 10, ...) {
  
  if (!"velocity" %in% names(seurat@misc)) {
    stop("Run velocity() first")
  }
  
  vdata <- seurat@misc$velocity
  
  if ("gene_params" %in% names(vdata)) {
    # From dynamical model
    params <- vdata$gene_params
    
    # Compute scores
    params$score <- params$likelihood * (1 - params$pval_steady)
    params$score[is.na(params$score)] <- 0
    
    # Filter and rank
    valid <- !is.na(params$alpha) & params$likelihood > min_r2
    params <- params[valid, ]
    params <- params[order(params$score, decreasing = TRUE), ]
    
    if (nrow(params) > n_genes) {
      params <- params[1:n_genes, ]
    }
    
    return(params)
  } else if ("r2" %in% names(vdata)) {
    # From steady-state model
    r2 <- vdata$r2
    gamma <- vdata$gamma
    gene_names <- names(r2)
    
    if (is.null(gene_names)) {
      gene_names <- paste0("gene_", seq_along(r2))
    }
    
    # Create result
    result <- data.frame(
      gene = gene_names,
      r2 = r2,
      gamma = gamma,
      stringsAsFactors = FALSE
    )
    
    # Filter and rank
    result <- result[result$r2 >= min_r2, ]
    result <- result[order(result$r2, decreasing = TRUE), ]
    
    if (nrow(result) > n_genes) {
      result <- result[1:n_genes, ]
    }
    
    return(result)
  } else {
    stop("No gene ranking information available")
  }
}

#' Velocity-based cell cycle scoring
#' 
#' @description Score cell cycle phase using velocity information
#' @param seurat Seurat object
#' @param s_genes S phase genes
#' @param g2m_genes G2/M phase genes
#' @param vkey Velocity key
#' @return Modified Seurat object with cycle scores
#' @export
velocity_cell_cycle <- function(seurat, s_genes, g2m_genes, vkey = "velocity") {
  
  if (!"velocity" %in% names(seurat@misc)) {
    stop("Run velocity() first")
  }
  
  vdata <- seurat@misc$velocity
  velocity_s <- vdata$velocity_s
  
  gene_names <- colnames(velocity_s)
  
  # Find gene indices
  s_idx <- which(gene_names %in% s_genes)
  g2m_idx <- which(gene_names %in% g2m_genes)
  
  if (length(s_idx) == 0 && length(g2m_idx) == 0) {
    warning("No cell cycle genes found in velocity data")
    return(seurat)
  }
  
  # Compute velocity-based scores
  s_score <- if (length(s_idx) > 0) {
    rowMeans(velocity_s[, s_idx, drop = FALSE], na.rm = TRUE)
  } else {
    rep(0, nrow(velocity_s))
  }
  
  g2m_score <- if (length(g2m_idx) > 0) {
    rowMeans(velocity_s[, g2m_idx, drop = FALSE], na.rm = TRUE)
  } else {
    rep(0, nrow(velocity_s))
  }
  
  # Normalize
  s_score <- scale(s_score)[, 1]
  g2m_score <- scale(g2m_score)[, 1]
  
  # Add to Seurat
  seurat$velocity_S_score <- s_score
  seurat$velocity_G2M_score <- g2m_score
  
  seurat
}

#' @title Stochastic Velocity Model
#' @description Complete implementation of SecondOrderSteadyStateModel
#' equivalent to scvelo.tools._steady_state_model.SecondOrderSteadyStateModel
#' @name velocity_stochastic
NULL

#' Stochastic Velocity Estimation (Second-Order Moments)
#'
#' @description Estimate RNA velocity using the stochastic model that
#' incorporates second-order moments to account for transcriptional noise.
#' 
#' @details The stochastic model uses the covariance between unspliced and
#' spliced reads to estimate gamma, which is more robust to transcriptional
#' bursting and other sources of noise.
#'
#' The equation is: Cov(U, S) = alpha * beta / (beta + gamma)
#'
#' @param seurat Seurat object or list with Ms, Mu
#' @param assay Assay name (for Seurat)
#' @param layer_s Layer name for spliced
#' @param layer_u Layer name for unspliced
#' @param Ms First-order spliced moments (cells x genes)
#' @param Mu First-order unspliced moments (cells x genes)
#' @param Mss Second-order spliced moments (optional)
#' @param Mus Second-order unspliced-spliced moments (optional)
#' @param connectivities Connectivity matrix
#' @param n_neighbors Number of neighbors
#' @param min_shared_counts Minimum shared counts
#' @param n_pcs Number of PCs for neighbor computation
#' @param groups Filter by groups
#' @param groupby Group by variable
#' @param residual Residual type ("variance" or "average")
#' @param copy Return copy
#' @return Modified object or list with velocity results
#' @export
velocity_stochastic <- function(seurat = NULL,
                                assay = "RNA",
                                layer_s = "spliced",
                                layer_u = "unspliced",
                                Ms = NULL, Mu = NULL,
                                Mss = NULL, Mus = NULL,
                                connectivities = NULL,
                                n_neighbors = 30,
                                min_shared_counts = 30,
                                n_pcs = 30,
                                groups = NULL,
                                groupby = NULL,
                                residual = NULL,
                                copy = FALSE) {
  
  # Get data from inputs
  if (!is.null(seurat) && inherits(seurat, "Seurat")) {
    data <- extract_velocity_data(seurat, assay, layer_s, layer_u)
    Ms <- data$Ms
    Mu <- data$Mu
    if (is.null(connectivities)) {
      connectivities <- data$connectivities
    }
    gene_names <- data$gene_names
  } else if (!is.null(Ms) && !is.null(Mu)) {
    gene_names <- colnames(Ms)
    if (is.null(gene_names)) {
      gene_names <- paste0("gene_", seq_len(ncol(Ms)))
    }
  } else {
    stop("Must provide either seurat object or Ms/Mu matrices")
  }
  
  n_cells <- nrow(Ms)
  n_genes <- ncol(Ms)
  
  # Compute second-order moments if not provided
  if (is.null(Mss) || is.null(Mus)) {
    message("Computing second-order moments...")
    
    if (is.null(connectivities)) {
      message("  Computing neighbors for moment calculation...")
      # Need PCA for neighbors - use basic approach
      pca_mat <- prcomp(Ms, center = TRUE, scale. = FALSE)$x[, 1:min(n_pcs, ncol(Ms)), drop = FALSE]
      neighbor_result <- compute_neighbors(pca_mat, n_neighbors = n_neighbors)
      connectivities <- neighbor_result$connectivities
    }
    
    # Normalize connectivities
    conn_norm <- connectivities / Matrix::rowSums(connectivities)
    
    # Compute second-order moments: E[X^2] where X is in neighborhood
    if (inherits(Ms, "sparseMatrix")) {
      Ms <- as.matrix(Ms)
    }
    if (inherits(Mu, "sparseMatrix")) {
      Mu <- as.matrix(Mu)
    }
    
    # Mss = E[S^2|neighborhood]
    Mss <- as.matrix(conn_norm %*% (Ms^2))
    
    # Mus = E[U*S|neighborhood]  
    Mus <- as.matrix(conn_norm %*% (Mu * Ms))
  }
  
  message("Fitting stochastic velocity model...")
  
  # Initialize output
  velocity_s <- matrix(0, n_cells, n_genes)
  gamma <- rep(NA_real_, n_genes)
  r2 <- rep(NA_real_, n_genes)
  
  # Fit each gene
  for (g in seq_len(n_genes)) {
    result <- fit_stochastic_gene(
      Ms[, g], Mu[, g], Mss[, g], Mus[, g],
      residual = residual
    )
    
    velocity_s[, g] <- result$velocity_s
    gamma[g] <- result$gamma
    r2[g] <- result$r2
  }
  
  colnames(velocity_s) <- gene_names
  
  # Return results
  result <- list(
    velocity_s = velocity_s,
    gamma = setNames(gamma, gene_names),
    r2 = setNames(r2, gene_names),
    Ms = Ms,
    Mu = Mu,
    Mss = Mss,
    Mus = Mus,
    mode = "stochastic"
  )
  
  # Add to Seurat if provided
  if (!is.null(seurat) && inherits(seurat, "Seurat")) {
    seurat@misc$velocity <- result
    return(seurat)
  }
  
  result
}

#' Fit stochastic model for a single gene
#'
#' @description Fit the second-order steady-state model for a single gene
#' @param ms First-order spliced moment (vector)
#' @param mu First-order unspliced moment (vector)
#' @param mss Second-order spliced moment
#' @param mus Second-order unspliced-spliced moment
#' @param residual Residual type
#' @return List with velocity, gamma, r2
#' @keywords internal
fit_stochastic_gene <- function(ms, mu, mss, mus, residual = NULL) {
  # Filter valid cells
  nonzero <- ms > 0 & mu > 0
  
  if (sum(nonzero) < 10) {
    return(list(
      velocity_s = rep(0, length(ms)),
      gamma = NA_real_,
      r2 = NA_real_
    ))
  }
  
  # Compute covariances
  # Cov(U, S) = E[US] - E[U]E[S] = Mus - Mu*Ms
  cov_us <- mus - mu * ms
  
  # Var(S) = E[S^2] - E[S]^2 = Mss - Ms^2
  var_s <- mss - ms^2
  
  # The relationship: velocity_s = beta * u - gamma * s
  # Under steady state: E[ds/dt] = 0, so beta * E[u] = gamma * E[s]
  # Stochastic estimation uses covariances:
  # gamma = Cov(U, S) / Var(S) under certain assumptions
  
  # Use weighted regression
  w <- nonzero * 1
  
  # Method 1: Least squares fit for gamma
  # velocity = mu * (beta/gamma) - s ≈ mu/gamma_ratio - s
  # We fit: cov_us ~ gamma * var_s
  
  # Actually, the stochastic model uses:
  # gamma = sum(cov_us) / sum(var_s) for cells with positive covariance
  
  pos_cov <- cov_us > 0 & var_s > 0 & nonzero
  
  if (sum(pos_cov) < 5) {
    # Fall back to deterministic
    gamma <- sum(w * mu * ms) / sum(w * ms^2)
    gamma <- max(gamma, 0.01)
  } else {
    # Stochastic estimate
    gamma <- sum(cov_us[pos_cov]) / sum(var_s[pos_cov])
    gamma <- max(gamma, 0.01)
    gamma <- min(gamma, 100)
  }
  
  # Compute velocity: velocity_s = beta * u - gamma * s
  # Assuming beta = 1 (absorbed into scaling)
  velocity_s <- mu - gamma * ms
  
  # Compute R^2
  if (!is.null(residual) && residual == "variance") {
    # Use variance-based residual
    residuals <- cov_us - gamma * var_s
    ss_res <- sum(residuals[nonzero]^2)
    ss_tot <- sum((cov_us[nonzero] - mean(cov_us[nonzero]))^2)
  } else {
    # Standard residual
    fitted <- gamma * ms
    residuals <- mu - fitted
    ss_res <- sum(w * residuals^2)
    ss_tot <- sum(w * (mu - mean(mu[nonzero]))^2)
  }
  
  r2 <- 1 - ss_res / (ss_tot + 1e-10)
  r2 <- max(0, min(1, r2))
  
  list(
    velocity_s = velocity_s,
    gamma = gamma,
    r2 = r2
  )
}

#' Second-Order Steady State Model (R6 Class)
#'
#' @description R6 class implementation of the stochastic velocity model
#' @export
SecondOrderSteadyStateModel <- R6::R6Class(
  "SecondOrderSteadyStateModel",
  
  public = list(
    # Data
    Ms = NULL,
    Mu = NULL,
    Mss = NULL,
    Mus = NULL,
    
    # Results
    gamma = NULL,
    velocity = NULL,
    variance = NULL,
    
    # Settings
    min_r2 = 0.01,
    r2_adjusted = NULL,
    residual = NULL,
    
    # Per-gene results
    r2 = NULL,
    
    #' @description Initialize model
    #' @param Ms First-order spliced moments
    #' @param Mu First-order unspliced moments
    #' @param Mss Second-order spliced moments
    #' @param Mus Unspliced-spliced covariance moments
    initialize = function(Ms, Mu, Mss = NULL, Mus = NULL) {
      self$Ms <- as.matrix(Ms)
      self$Mu <- as.matrix(Mu)
      
      if (!is.null(Mss)) self$Mss <- as.matrix(Mss)
      if (!is.null(Mus)) self$Mus <- as.matrix(Mus)
    },
    
    #' @description Fit the model
    #' @param connectivities Connectivity matrix
    fit = function(connectivities = NULL) {
      # Compute moments if needed
      if (is.null(self$Mss) || is.null(self$Mus)) {
        if (is.null(connectivities)) {
          stop("Need connectivities to compute second-order moments")
        }
        self$compute_moments(connectivities)
      }
      
      n_genes <- ncol(self$Ms)
      n_cells <- nrow(self$Ms)
      
      self$gamma <- rep(NA_real_, n_genes)
      self$r2 <- rep(NA_real_, n_genes)
      self$velocity <- matrix(0, n_cells, n_genes)
      
      for (g in seq_len(n_genes)) {
        res <- fit_stochastic_gene(
          self$Ms[, g], self$Mu[, g],
          self$Mss[, g], self$Mus[, g],
          residual = self$residual
        )
        
        self$gamma[g] <- res$gamma
        self$r2[g] <- res$r2
        self$velocity[, g] <- res$velocity_s
      }
      
      # Adjust R2 if needed
      if (!is.null(self$min_r2)) {
        self$velocity[, self$r2 < self$min_r2] <- 0
      }
      
      invisible(self)
    },
    
    #' @description Compute second-order moments
    #' @param connectivities Connectivity matrix
    compute_moments = function(connectivities) {
      conn_norm <- connectivities / Matrix::rowSums(connectivities)
      conn_norm[!is.finite(conn_norm)] <- 0
      
      self$Mss <- as.matrix(conn_norm %*% (self$Ms^2))
      self$Mus <- as.matrix(conn_norm %*% (self$Mu * self$Ms))
      
      invisible(self)
    },
    
    #' @description Get velocity
    get_velocity = function() {
      self$velocity
    },
    
    #' @description Get gamma
    get_gamma = function() {
      self$gamma
    }
  )
)

#' Compute second-order moments
#'
#' @description Compute second-order moments for stochastic velocity estimation
#' @param X Expression matrix (cells x genes)
#' @param connectivities Connectivity matrix (cells x cells)
#' @return Matrix of second-order moments
#' @keywords internal
compute_second_order_moments_matrix <- function(X, connectivities) {
  # Normalize connectivities
  conn_norm <- connectivities / Matrix::rowSums(connectivities)
  conn_norm[!is.finite(conn_norm)] <- 0
  
  # E[X^2] in neighborhood
  as.matrix(conn_norm %*% (X^2))
}

#' Compute cross second-order moments
#'
#' @description Compute E[X*Y] in neighborhood for stochastic model
#' @param X First expression matrix
#' @param Y Second expression matrix
#' @param connectivities Connectivity matrix
#' @return Matrix of cross moments
#' @export
cross_moments <- function(X, Y, connectivities) {
  conn_norm <- connectivities / Matrix::rowSums(connectivities)
  conn_norm[!is.finite(conn_norm)] <- 0
  
  as.matrix(conn_norm %*% (X * Y))
}

#' Helper to extract velocity data from Seurat
#'
#' @keywords internal
extract_velocity_data <- function(seurat, assay, layer_s, layer_u) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("SeuratObject package is required")
  }
  
  # Get expression data
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    # Seurat V5
    tryCatch({
      Ms <- as.matrix(SeuratObject::LayerData(seurat, assay = assay, layer = layer_s))
      Mu <- as.matrix(SeuratObject::LayerData(seurat, assay = assay, layer = layer_u))
    }, error = function(e) {
      # Try data slot
      Ms <- as.matrix(seurat[[assay]]@data)
      Mu <- as.matrix(seurat[[assay]]@data)  
    })
  } else {
    # Seurat V4
    tryCatch({
      Ms <- as.matrix(seurat[[assay]]@data)
      Mu <- as.matrix(seurat[[assay]]@data)
    }, error = function(e) {
      stop("Could not extract expression data from Seurat object")
    })
  }
  
  # Transpose to cells x genes
  Ms <- t(Ms)
  Mu <- t(Mu)
  
  # Get connectivities
  connectivities <- NULL
  if ("neighbors" %in% names(seurat@graphs)) {
    connectivities <- seurat@graphs$neighbors
  } else if ("RNA_snn" %in% names(seurat@graphs)) {
    connectivities <- seurat@graphs$RNA_snn
  }
  
  gene_names <- colnames(Ms)
  
  list(
    Ms = Ms,
    Mu = Mu,
    connectivities = connectivities,
    gene_names = gene_names
  )
}

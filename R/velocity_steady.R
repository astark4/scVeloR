#' @title Steady-State Velocity Model
#' @description Functions for estimating RNA velocity using the steady-state model.
#' @name velocity_steady
NULL

#' Fit Steady-State Velocity Model
#'
#' @description Estimate RNA velocity assuming cells are at steady-state equilibrium.
#' The model fits: velocity = unspliced - gamma * spliced
#'
#' @param seurat_obj A Seurat object with computed moments (Ms, Mu)
#' @param fit_intercept Whether to fit an intercept term (default: FALSE)
#' @param percentile Percentile of cells to use for fitting (default: c(5, 95))
#' @param min_r2 Minimum R-squared to consider a gene well-fit (default: 0.01)
#' @param n_jobs Number of parallel jobs (default: 1)
#' @param use_raw Use raw counts instead of moments (default: FALSE)
#'
#' @return Modified Seurat object with velocity and fit parameters
#'
#' @export
fit_velocity_steady <- function(seurat_obj,
                                 fit_intercept = FALSE,
                                 percentile = c(5, 95),
                                 min_r2 = 0.01,
                                 n_jobs = 1L,
                                 use_raw = FALSE) {
  
  .vmessage("Fitting steady-state velocity model")
  
  # Get data
  if (use_raw) {
    data <- extract_velocity_data(seurat_obj)
    Ms <- data$spliced
    Mu <- data$unspliced
  } else {
    # Check for moments
    if (!.layer_exists(seurat_obj, "Ms") || !.layer_exists(seurat_obj, "Mu")) {
      stop("Moments not found. Run compute_moments first.", call. = FALSE)
    }
    Ms <- .get_layer_data(seurat_obj, "Ms")
    Mu <- .get_layer_data(seurat_obj, "Mu")
  }
  
  Ms <- make_dense(Ms)
  Mu <- make_dense(Mu)
  
  # Get genes to fit
  genes <- rownames(Ms)
  n_genes <- length(genes)
  n_cells <- ncol(Ms)
  
  # Get velocity genes mask
  velocity_genes <- get_velocity_genes_mask(seurat_obj, genes)
  
  # Initialize results
  gamma <- numeric(n_genes)
  offset <- numeric(n_genes)
  r2 <- numeric(n_genes)
  velocity_genes_result <- rep(FALSE, n_genes)
  
  names(gamma) <- genes
  names(offset) <- genes
  names(r2) <- genes
  names(velocity_genes_result) <- genes
  
  # Fit each gene
  .vmessage("Fitting ", sum(velocity_genes), " genes...")
  
  genes_to_fit <- genes[velocity_genes]
  
  fit_gene_steady <- function(gene_idx) {
    s <- Ms[gene_idx, ]
    u <- Mu[gene_idx, ]
    
    # Filter by percentile
    s_range <- stats::quantile(s, percentile / 100, na.rm = TRUE)
    u_range <- stats::quantile(u, percentile / 100, na.rm = TRUE)
    
    mask <- (s >= s_range[1] | s <= s_range[2]) &
            (u >= u_range[1] | u <= u_range[2])
    
    s_fit <- s[mask]
    u_fit <- u[mask]
    
    if (length(s_fit) < 10) {
      return(list(gamma = NA, offset = NA, r2 = 0))
    }
    
    # Fit linear model
    if (fit_intercept) {
      fit <- tryCatch({
        model <- stats::lm(u_fit ~ s_fit)
        coefs <- stats::coef(model)
        list(
          gamma = unname(coefs[2]),
          offset = unname(coefs[1])
        )
      }, error = function(e) {
        list(gamma = NA, offset = NA)
      })
    } else {
      # No intercept: gamma = sum(s*u) / sum(s^2)
      ss <- sum(s_fit^2)
      if (ss > 0) {
        fit <- list(
          gamma = sum(s_fit * u_fit) / ss,
          offset = 0
        )
      } else {
        fit <- list(gamma = NA, offset = 0)
      }
    }
    
    # Calculate R-squared
    if (!is.na(fit$gamma)) {
      predicted <- fit$gamma * s + fit$offset
      residual <- u - predicted
      ss_res <- sum(residual^2)
      ss_tot <- sum((u - mean(u))^2)
      r2_val <- ifelse(ss_tot > 0, 1 - ss_res / ss_tot, 0)
    } else {
      r2_val <- 0
    }
    
    list(gamma = fit$gamma, offset = fit$offset, r2 = r2_val)
  }
  
  # Run fitting
  results <- apply_with_progress(
    seq_len(n_genes),
    fit_gene_steady,
    n_jobs = n_jobs,
    show_progress = TRUE
  )
  
  # Extract results
  for (i in seq_len(n_genes)) {
    if (!is.null(results[[i]])) {
      gamma[i] <- results[[i]]$gamma
      offset[i] <- results[[i]]$offset
      r2[i] <- results[[i]]$r2
      velocity_genes_result[i] <- velocity_genes[i] && 
                                   !is.na(gamma[i]) && 
                                   gamma[i] > 0 &&
                                   r2[i] >= min_r2
    }
  }
  
  # Calculate velocity: v = u - gamma * s
  .vmessage("Computing velocity...")
  
  velocity <- Mu - sweep(Ms, 1, gamma, "*")
  if (fit_intercept) {
    velocity <- velocity - offset
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
    fit_offset = offset,
    fit_r2 = r2,
    velocity_genes = velocity_genes_result
  )
  
  seurat_obj <- store_fit_params(seurat_obj, fit_params)
  seurat_obj@misc$velocity_mode <- "steady_state"
  
  n_fit <- sum(velocity_genes_result, na.rm = TRUE)
  .vmessage("Fitted ", n_fit, " velocity genes")
  
  return(seurat_obj)
}

#' Get Velocity Genes Mask
#'
#' @description Get boolean mask for velocity genes.
#'
#' @param seurat_obj Seurat object
#' @param genes Gene names
#' @return Logical vector
#' @keywords internal
get_velocity_genes_mask <- function(seurat_obj, genes) {
  
  # Try to get from misc first
  if (!is.null(seurat_obj@misc$velocity_genes)) {
    vg <- seurat_obj@misc$velocity_genes
    if (is.logical(vg)) {
      if (length(vg) == length(genes)) {
        return(vg)
      }
    } else if (!is.null(names(vg))) {
      return(genes %in% names(vg)[vg])
    }
  }
  
  # Try from assay metadata
  assay_name <- SeuratObject::DefaultAssay(seurat_obj)
  tryCatch({
    meta <- seurat_obj[[assay_name]]@meta.features
    if ("velocity_genes" %in% colnames(meta)) {
      return(meta$velocity_genes[match(genes, rownames(meta))])
    }
  }, error = function(e) NULL)
  
  # Default: all genes
  rep(TRUE, length(genes))
}

#' Store Fit Parameters
#'
#' @description Store velocity fit parameters in Seurat object.
#'
#' @param seurat_obj Seurat object
#' @param fit_params Data frame of fit parameters
#' @return Modified Seurat object
#' @keywords internal
store_fit_params <- function(seurat_obj, fit_params) {
  
  assay_name <- SeuratObject::DefaultAssay(seurat_obj)
  
  tryCatch({
    meta <- seurat_obj[[assay_name]]@meta.features
    
    for (col in colnames(fit_params)) {
      idx <- match(rownames(fit_params), rownames(meta))
      meta[[col]] <- NA
      meta[[col]][idx[!is.na(idx)]] <- fit_params[[col]][!is.na(idx)]
    }
    
    seurat_obj[[assay_name]]@meta.features <- meta
  }, error = function(e) {
    # Store in misc as fallback
    seurat_obj@misc$fit_params <- fit_params
  })
  
  return(seurat_obj)
}

#' Identify Velocity Genes
#'
#' @description Identify genes suitable for velocity analysis based on fit quality.
#'
#' @param seurat_obj A Seurat object with fitted velocity
#' @param min_r2 Minimum R-squared (default: 0.01)
#' @param min_corr Minimum correlation (default: 0.3)
#' @param min_ratio Minimum unspliced/spliced ratio (default: 0.1)
#' @param n_top_genes Number of top genes to select (default: NULL for all passing)
#'
#' @return Character vector of velocity gene names
#'
#' @export
velocity_genes <- function(seurat_obj,
                            min_r2 = 0.01,
                            min_corr = 0.3,
                            min_ratio = 0.1,
                            n_top_genes = NULL) {
  
  # Get fit parameters
  fit_params <- get_fit_params(seurat_obj)
  
  if (is.null(fit_params)) {
    stop("Fit parameters not found. Run velocity estimation first.", call. = FALSE)
  }
  
  genes <- rownames(fit_params)
  
  # Apply filters
  mask <- rep(TRUE, nrow(fit_params))
  
  if ("fit_r2" %in% colnames(fit_params)) {
    mask <- mask & !is.na(fit_params$fit_r2) & (fit_params$fit_r2 >= min_r2)
  }
  
  if ("fit_gamma" %in% colnames(fit_params)) {
    mask <- mask & !is.na(fit_params$fit_gamma) & (fit_params$fit_gamma > 0)
  }
  
  if ("velocity_genes" %in% colnames(fit_params)) {
    mask <- mask & fit_params$velocity_genes
  }
  
  velocity_genes <- genes[mask]
  
  # Select top genes if specified
  if (!is.null(n_top_genes) && length(velocity_genes) > n_top_genes) {
    r2_vals <- fit_params$fit_r2[mask]
    order_idx <- order(r2_vals, decreasing = TRUE)
    velocity_genes <- velocity_genes[order_idx[seq_len(n_top_genes)]]
  }
  
  .vmessage("Selected ", length(velocity_genes), " velocity genes")
  
  return(velocity_genes)
}

#' Get Fit Parameters
#'
#' @description Get velocity fit parameters from Seurat object.
#'
#' @param seurat_obj Seurat object
#' @return Data frame of fit parameters or NULL
#' @keywords internal
get_fit_params <- function(seurat_obj) {
  
  # Try misc first
  if (!is.null(seurat_obj@misc$fit_params)) {
    return(seurat_obj@misc$fit_params)
  }
  
  # Try assay meta features
  assay_name <- SeuratObject::DefaultAssay(seurat_obj)
  
  tryCatch({
    meta <- seurat_obj[[assay_name]]@meta.features
    
    fit_cols <- c("fit_gamma", "fit_offset", "fit_r2", "velocity_genes",
                  "fit_alpha", "fit_beta", "fit_t_", "fit_scaling",
                  "fit_likelihood", "fit_alignment")
    
    available_cols <- intersect(fit_cols, colnames(meta))
    
    if (length(available_cols) > 0) {
      return(meta[, available_cols, drop = FALSE])
    }
    
    NULL
  }, error = function(e) NULL)
}

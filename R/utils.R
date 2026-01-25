#' @title Utility Functions
#' @description Internal utility functions for scVeloR
#' @name utils
NULL

#' Check Verbose Option
#' @keywords internal
.is_verbose <- function() {
  getOption("scVeloR.verbose", TRUE)
}

#' Message if Verbose
#' @keywords internal
.vmessage <- function(...) {
  if (.is_verbose()) {
    message(...)
  }
}

#' Safe Log Transform
#'
#' @description Log transform with clipping to avoid -Inf values.
#'
#' @param x Numeric vector or matrix
#' @param lb Lower bound (default: 0)
#' @param ub Upper bound (default: 1)
#' @param eps Epsilon for clipping (default: 1e-6)
#'
#' @return Log-transformed values
#' @keywords internal
clipped_log <- function(x, lb = 0, ub = 1, eps = 1e-6) {
  log(pmax(pmin(x, ub - eps), lb + eps))
}

#' Safe Inversion
#'
#' @description Invert values, setting Inf to 0.
#'
#' @param x Numeric vector or matrix
#' @return Inverted values with Inf replaced by 0
#' @keywords internal
safe_invert <- function(x) {
  result <- 1 / x
  result[!is.finite(result)] <- 0
  result
}

#' L2 Norm
#'
#' @description Calculate L2 norm along rows or columns.
#'
#' @param x Matrix
#' @param axis 1 for row-wise, 2 for column-wise
#' @return Vector of L2 norms
#' @keywords internal
l2_norm <- function(x, axis = 1) {
  if (axis == 1) {
    sqrt(Matrix::rowSums(x^2))
  } else {
    sqrt(Matrix::colSums(x^2))
  }
}

#' Normalize Matrix
#'
#' @description Normalize matrix rows or columns to sum to 1.
#'
#' @param x Matrix (can be sparse)
#' @param axis 1 for row normalization, 2 for column normalization
#' @param min_val Minimum value to add to avoid division by zero
#' @return Normalized matrix
#' @keywords internal
normalize_matrix <- function(x, axis = 1, min_val = 1e-10) {
  if (axis == 1) {
    row_sums <- Matrix::rowSums(x)
    row_sums[row_sums < min_val] <- min_val
    x / row_sums
  } else {
    col_sums <- Matrix::colSums(x)
    col_sums[col_sums < min_val] <- min_val
    t(t(x) / col_sums)
  }
}

#' Sparse Matrix Row Product Sum
#'
#' @description Calculate sum of element-wise products along axis.
#'
#' @param a First matrix
#' @param b Second matrix
#' @param axis Axis for summation (1 = rows, 2 = columns)
#' @return Vector of sums
#' @keywords internal
prod_sum <- function(a, b, axis = 1) {
  if (inherits(a, "sparseMatrix") || inherits(b, "sparseMatrix")) {
    ab <- a * b
    if (axis == 1) {
      return(Matrix::rowSums(ab))
    } else {
      return(Matrix::colSums(ab))
    }
  } else {
    if (axis == 1) {
      return(rowSums(a * b))
    } else {
      return(colSums(a * b))
    }
  }
}

#' Efficient Matrix Sum
#'
#' @description Sum matrix elements along axis, handling sparse matrices.
#'
#' @param x Matrix (can be sparse)
#' @param axis 1 for row sums, 2 for column sums, NULL for total sum
#' @return Vector of sums or scalar
#' @keywords internal
mat_sum <- function(x, axis = NULL) {
  if (is.null(axis)) {
    return(sum(x))
  }
  
  if (inherits(x, "sparseMatrix")) {
    if (axis == 1) {
      return(Matrix::rowSums(x))
    } else {
      return(Matrix::colSums(x))
    }
  } else {
    if (axis == 1) {
      return(rowSums(x))
    } else {
      return(colSums(x))
    }
  }
}

#' Multiply Sparse Matrices Element-wise
#'
#' @description Element-wise multiplication handling sparse matrices.
#'
#' @param a First matrix
#' @param b Second matrix
#' @return Element-wise product
#' @keywords internal
sparse_multiply <- function(a, b) {
  if (inherits(a, "sparseMatrix")) {
    return(a * b)
  } else if (inherits(b, "sparseMatrix")) {
    return(b * a)
  } else {
    return(a * b)
  }
}

#' Scale Data
#'
#' @description Scale data to [0, 1] range.
#'
#' @param x Numeric vector or matrix
#' @return Scaled values
#' @keywords internal
scale_01 <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  if (x_max == x_min) {
    return(rep(0.5, length(x)))
  }
  (x - x_min) / (x_max - x_min)
}

#' Get Percentile Values
#'
#' @description Get values within specified percentile range.
#'
#' @param x Numeric vector or matrix
#' @param perc Percentile(s) - single value or c(lower, upper)
#' @return Logical mask for values within percentile range
#' @keywords internal
get_percentile_mask <- function(x, perc) {
  if (length(perc) == 1) {
    # Upper percentile only
    threshold <- stats::quantile(x, perc / 100, na.rm = TRUE)
    return(x >= threshold)
  } else {
    # Both lower and upper
    lower <- stats::quantile(x, perc[1] / 100, na.rm = TRUE)
    upper <- stats::quantile(x, perc[2] / 100, na.rm = TRUE)
    return(x <= lower | x >= upper)
  }
}

#' R-squared Calculation
#'
#' @description Calculate coefficient of determination.
#'
#' @param residual Residuals
#' @param total Total variation (y - mean(y))
#' @return R-squared value(s)
#' @keywords internal
R_squared <- function(residual, total) {
  ss_res <- Matrix::colSums(residual^2)
  ss_tot <- Matrix::colSums(total^2)
  ss_tot[ss_tot == 0] <- 1  # Avoid division by zero
  1 - ss_res / ss_tot
}

#' Make Dense Matrix
#'
#' @description Convert sparse matrix to dense if needed.
#'
#' @param x Matrix
#' @return Dense matrix
#' @keywords internal
make_dense <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    as.matrix(x)
  } else {
    x
  }
}

#' Make Sparse Matrix
#'
#' @description Convert dense matrix to sparse.
#'
#' @param x Matrix
#' @return Sparse matrix
#' @keywords internal
make_sparse <- function(x) {
  if (!inherits(x, "sparseMatrix")) {
    as(x, "CsparseMatrix")
  } else {
    x
  }
}

#' Groups to Boolean Mask
#'
#' @description Convert group specification to boolean mask.
#'
#' @param seurat_obj Seurat object
#' @param groups Groups to select
#' @param groupby Column name in metadata
#' @return Boolean vector
#' @keywords internal
groups_to_bool <- function(seurat_obj, groups, groupby) {
  if (is.null(groups)) {
    return(NULL)
  }
  
  if (is.null(groupby)) {
    stop("groupby must be specified when groups is provided", call. = FALSE)
  }
  
  meta <- seurat_obj@meta.data
  if (!groupby %in% colnames(meta)) {
    stop("Column '", groupby, "' not found in metadata", call. = FALSE)
  }
  
  meta[[groupby]] %in% groups
}

#' Validate Seurat Object
#'
#' @description Check that a Seurat object has required components.
#'
#' @param seurat_obj A Seurat object
#' @param require_velocity Whether velocity data is required
#' @param require_graph Whether velocity graph is required
#' @param require_embedding Which embedding is required (NULL for none)
#' @return TRUE if valid, otherwise throws error
#' @keywords internal
validate_seurat <- function(seurat_obj, 
                            require_velocity = FALSE,
                            require_graph = FALSE,
                            require_embedding = NULL) {
  
  if (!inherits(seurat_obj, "Seurat")) {
    stop("Input must be a Seurat object", call. = FALSE)
  }
  
  if (require_velocity) {
    if (!.layer_exists(seurat_obj, "velocity")) {
      stop("Velocity data not found. Run velocity() first.", call. = FALSE)
    }
  }
  
  if (require_graph) {
    if (is.null(seurat_obj@misc[["velocity_graph"]])) {
      stop("Velocity graph not found. Run compute_velocity_graph() first.", call. = FALSE)
    }
  }
  
  if (!is.null(require_embedding)) {
    emb_name <- tolower(gsub("^X_", "", require_embedding))
    if (!emb_name %in% names(seurat_obj@reductions)) {
      stop("Embedding '", require_embedding, "' not found.", call. = FALSE)
    }
  }
  
  TRUE
}

#' Get Number of Jobs
#'
#' @description Determine number of parallel jobs to use.
#'
#' @param n_jobs Number of jobs requested (NULL = default, -1 = all cores)
#' @return Integer number of jobs
#' @keywords internal
get_n_jobs <- function(n_jobs = NULL) {
  if (is.null(n_jobs)) {
    n_jobs <- getOption("scVeloR.n_jobs", 1L)
  }
  
  max_cores <- parallel::detectCores()
  
  if (n_jobs < 0) {
    n_jobs <- max(1L, max_cores + n_jobs + 1L)
  } else if (n_jobs > max_cores) {
    n_jobs <- max_cores
  }
  
  as.integer(max(1L, n_jobs))
}

#' Progress Wrapper
#'
#' @description Apply function with optional progress bar.
#'
#' @param X Vector to iterate over
#' @param FUN Function to apply
#' @param ... Additional arguments to FUN
#' @param n_jobs Number of parallel jobs
#' @param show_progress Whether to show progress bar
#' @return List of results
#' @keywords internal
apply_with_progress <- function(X, FUN, ..., n_jobs = 1L, show_progress = TRUE) {
  n_jobs <- get_n_jobs(n_jobs)
  
  if (n_jobs == 1L) {
    # Sequential
    if (show_progress && getOption("scVeloR.progress", TRUE)) {
      progressr::with_progress({
        p <- progressr::progressor(along = X)
        result <- lapply(X, function(x) {
          res <- FUN(x, ...)
          p()
          res
        })
      })
    } else {
      result <- lapply(X, FUN, ...)
    }
  } else {
    # Parallel
    if (show_progress && getOption("scVeloR.progress", TRUE)) {
      progressr::with_progress({
        p <- progressr::progressor(along = X)
        result <- future.apply::future_lapply(X, function(x) {
          res <- FUN(x, ...)
          p()
          res
        }, future.seed = TRUE)
      })
    } else {
      result <- future.apply::future_lapply(X, FUN, ..., future.seed = TRUE)
    }
  }
  
  result
}

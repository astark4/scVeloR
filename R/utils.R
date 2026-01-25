#' @title Utility Functions
#' @description General utility functions for scVeloR.
#' @name utils
NULL

#' Make Dense Matrix
#'
#' @description Convert sparse matrix to dense if needed.
#'
#' @param x Matrix (sparse or dense)
#'
#' @return Dense matrix
#' @keywords internal
make_dense <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    as.matrix(x)
  } else if (is.matrix(x)) {
    x
  } else {
    as.matrix(x)
  }
}

#' Clip Values
#'
#' @description Clip values to specified range.
#'
#' @param x Numeric vector or matrix
#' @param lower Lower bound
#' @param upper Upper bound
#'
#' @return Clipped values
#' @keywords internal
clip <- function(x, lower = -Inf, upper = Inf) {
  pmin(pmax(x, lower), upper)
}

#' Log Transform with Offset
#'
#' @description Log transform with offset to handle zeros: log(x + offset)
#'
#' @param x Numeric vector or matrix
#' @param offset Offset value (default 1)
#'
#' @return Log-transformed values
#' @keywords internal
log1p_offset <- function(x, offset = 1) {
  log(x + offset)
}

#' Scale Vector
#'
#' @description Scale vector to [0, 1] range.
#'
#' @param x Numeric vector
#'
#' @return Scaled vector
#' @keywords internal
scale_minmax <- function(x) {
  x_min <- min(x, na.rm = TRUE)
  x_max <- max(x, na.rm = TRUE)
  
  if (x_max == x_min) {
    return(rep(0.5, length(x)))
  }
  
  (x - x_min) / (x_max - x_min)
}

#' Check If Object Has scVeloR Results
#'
#' @description Check if Seurat object has scVeloR results.
#'
#' @param object Seurat object
#' @param component Which component to check
#'
#' @return Logical
#' @keywords internal
has_scVeloR <- function(object, component = NULL) {
  if (is.null(object@misc$scVeloR)) {
    return(FALSE)
  }
  
  if (is.null(component)) {
    return(TRUE)
  }
  
  !is.null(object@misc$scVeloR[[component]])
}

#' Print Message with Timing
#'
#' @description Print message with timestamp.
#'
#' @param msg Message to print
#' @param verbose Whether to print
#'
#' @keywords internal
msg <- function(msg, verbose = TRUE) {
  if (verbose) {
    message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), msg))
  }
}

#' Row-Wise Variance
#'
#' @description Compute row-wise variance efficiently.
#'
#' @param x Matrix
#'
#' @return Vector of variances
#' @keywords internal
rowVars <- function(x) {
  n <- ncol(x)
  rowMeans(x^2) - rowMeans(x)^2
}

#' Column-Wise Variance
#'
#' @description Compute column-wise variance efficiently.
#'
#' @param x Matrix
#'
#' @return Vector of variances
#' @keywords internal
colVars <- function(x) {
  n <- nrow(x)
  colMeans(x^2) - colMeans(x)^2
}

#' Sparse Matrix Row Sums with Names
#'
#' @description Compute row sums of sparse matrix preserving names.
#'
#' @param x Sparse matrix
#'
#' @return Named vector
#' @keywords internal
sparse_row_sums <- function(x) {
  rs <- Matrix::rowSums(x)
  names(rs) <- rownames(x)
  rs
}

#' Sparse Matrix Column Sums with Names
#'
#' @description Compute column sums of sparse matrix preserving names.
#'
#' @param x Sparse matrix
#'
#' @return Named vector
#' @keywords internal
sparse_col_sums <- function(x) {
  cs <- Matrix::colSums(x)
  names(cs) <- colnames(x)
  cs
}

#' Weighted Mean
#'
#' @description Compute weighted mean.
#'
#' @param x Values
#' @param w Weights
#'
#' @return Weighted mean
#' @keywords internal
weighted_mean <- function(x, w) {
  sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

#' Check Seurat Version
#'
#' @description Check if Seurat V5 or V4.
#'
#' @param object Seurat object
#'
#' @return Integer (4 or 5)
#' @keywords internal
seurat_version <- function(object) {
  if ("Assay5" %in% class(object@assays[[Seurat::DefaultAssay(object)]])) {
    return(5L)
  }
  return(4L)
}

#' Safe Divide
#'
#' @description Division with protection against division by zero.
#'
#' @param a Numerator
#' @param b Denominator
#' @param default Default value when b is zero
#'
#' @return a/b with zeros handled
#' @keywords internal
safe_div <- function(a, b, default = 0) {
  result <- a / b
  result[!is.finite(result)] <- default
  result[b == 0] <- default
  result
}

#' Vector Correlation
#'
#' @description Compute correlation between vector and matrix columns.
#'
#' @param x Vector
#' @param Y Matrix
#'
#' @return Vector of correlations
#' @keywords internal
vcorrcoef <- function(x, Y) {
  n <- length(x)
  
  # Center
  x_c <- x - mean(x)
  Y_c <- Y - colMeans(Y)
  
  # Correlation
  x_sd <- sqrt(sum(x_c^2))
  Y_sd <- sqrt(colSums(Y_c^2))
  
  colSums(x_c * Y_c) / (x_sd * Y_sd)
}

#' Print Object Size
#'
#' @description Print size of R object in human-readable format.
#'
#' @param x R object
#'
#' @return Character string
#' @keywords internal
object_size_str <- function(x) {
  size <- object.size(x)
  
  if (size < 1024) {
    return(sprintf("%d B", size))
  } else if (size < 1024^2) {
    return(sprintf("%.1f KB", size / 1024))
  } else if (size < 1024^3) {
    return(sprintf("%.1f MB", size / 1024^2))
  } else {
    return(sprintf("%.1f GB", size / 1024^3))
  }
}

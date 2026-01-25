#' @title scVeloR: RNA Velocity Analysis for Single-Cell RNA-Seq Data
#'
#' @description
#' scVeloR provides a comprehensive R implementation of RNA velocity analysis
#' for single-cell RNA sequencing data. The package implements steady-state,
#' stochastic, and dynamical velocity models, along with visualization tools.
#'
#' @section Main Functions:
#' \describe{
#'   \item{Preprocessing}{
#'     \code{\link{filter_genes}}, \code{\link{normalize_layers}},
#'     \code{\link{compute_moments}}, \code{\link{filter_and_normalize}}
#'   }
#'   \item{Velocity Estimation}{
#'     \code{\link{velocity}}, \code{\link{fit_velocity_steady}},
#'     \code{\link{fit_velocity_stochastic}}, \code{\link{recover_dynamics}}
#'   }
#'   \item{Velocity Graph}{
#'     \code{\link{compute_velocity_graph}}, \code{\link{project_velocity_embedding}},
#'     \code{\link{transition_matrix}}
#'   }
#'   \item{Visualization}{
#'     \code{\link{velocity_embedding_plot}}, \code{\link{velocity_stream_plot}},
#'     \code{\link{velocity_grid_plot}}, \code{\link{velocity_heatmap}}
#'   }
#' }
#'
#' @section Seurat Compatibility:
#' scVeloR is fully compatible with Seurat V4 and V5 objects. The package
#' automatically detects the Seurat version and uses appropriate methods
#' for data access.
#'
#' @section Parallel Computing:
#' scVeloR supports parallel computation through the \code{future} framework.
#' Use \code{future::plan()} to configure parallel backends.
#'
#' @author Zaoqu Liu \email{liuzaoqu@@163.com}
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/Zaoqu-Liu/scVeloR}
#'   \item Report bugs at \url{https://github.com/Zaoqu-Liu/scVeloR/issues}
#' }
#'
#' @docType package
#' @name scVeloR-package
#' @aliases scVeloR
#' @useDynLib scVeloR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
"_PACKAGE"

# Package-level options
.onLoad <- function(libname, pkgname) {
  # Set default options
  op <- options()
  op.scVeloR <- list(
    scVeloR.verbose = TRUE,
    scVeloR.n_jobs = 1L,
    scVeloR.progress = TRUE
  )
  toset <- !(names(op.scVeloR) %in% names(op))
  if (any(toset)) options(op.scVeloR[toset])
  
  invisible()
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "scVeloR v", utils::packageVersion("scVeloR"), "\n",
    "RNA Velocity Analysis for Single-Cell RNA-Seq Data\n",
    "Author: Zaoqu Liu <liuzaoqu@163.com>\n",
    "GitHub: https://github.com/Zaoqu-Liu/scVeloR"
  )
}

# Global variables to avoid R CMD check notes
utils::globalVariables(c(
  ".", "x", "y", "vx", "vy", "color", "size", "alpha",
  "velocity", "spliced", "unspliced", "Ms", "Mu",
  "gene", "cell", "value", "cluster", "time",
  "start_x", "start_y", "end_x", "end_y",
  "grid_x", "grid_y", "magnitude"
))

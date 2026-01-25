#' @title scVeloR: RNA Velocity Analysis in R
#'
#' @description
#' A comprehensive R implementation of RNA velocity analysis based on the 
#' scvelo Python package. This package enables estimation of RNA velocities 
#' for single-cell RNA-seq data, supporting deterministic (steady-state), 
#' stochastic, and dynamical models.
#'
#' @section Main Functions:
#' \describe{
#'   \item{\code{\link{run_velocity}}}{Full velocity analysis pipeline}
#'   \item{\code{\link{velocity}}}{Compute RNA velocity}
#'   \item{\code{\link{velocity_graph}}}{Build velocity graph}
#'   \item{\code{\link{velocity_embedding}}}{Project velocity onto embedding}
#'   \item{\code{\link{compute_latent_time}}}{Compute gene-shared latent time}
#' }
#'
#' @section Velocity Models:
#' \describe{
#'   \item{Deterministic}{Steady-state model assuming equilibrium}
#'   \item{Stochastic}{Uses second-order moments for noise modeling}
#'   \item{Dynamical}{Full kinetics model with EM algorithm}
#' }
#'
#' @section Preprocessing:
#' \describe{
#'   \item{\code{\link{prepare_velocity}}}{Prepare data for velocity analysis}
#'   \item{\code{\link{filter_genes}}}{Filter genes by expression criteria}
#'   \item{\code{\link{moments}}}{Compute neighbor-averaged moments}
#'   \item{\code{\link{compute_neighbors}}}{Compute k-nearest neighbors}
#' }
#'
#' @section Visualization:
#' \describe{
#'   \item{\code{\link{plot_velocity}}}{Plot velocity arrows}
#'   \item{\code{\link{plot_velocity_stream}}}{Plot velocity streamlines}
#'   \item{\code{\link{plot_velocity_grid}}}{Plot velocity on grid}
#'   \item{\code{\link{plot_phase}}}{Plot phase portraits}
#' }
#'
#' @section Example Usage:
#' \preformatted{
#' library(scVeloR)
#' library(Seurat)
#'
#' # Load Seurat object with spliced and unspliced layers
#' seurat_obj <- readRDS("your_seurat_object.rds")
#'
#' # Run full velocity pipeline
#' seurat_obj <- run_velocity(seurat_obj, mode = "dynamical")
#'
#' # Visualize results
#' plot_velocity(seurat_obj, embedding = "umap")
#' plot_velocity_stream(seurat_obj, embedding = "umap")
#' }
#'
#' @author Zaoqu Liu \email{liuzaoqu@@163.com}
#' @references
#' Bergen, V., Lange, M., Peidli, S. et al. (2020). 
#' Generalizing RNA velocity to transient cell states through dynamical modeling. 
#' Nature Biotechnology. \doi{10.1038/s41587-020-0591-3}
#'
#' La Manno, G., Soldatov, R., Zeisel, A. et al. (2018). 
#' RNA velocity of single cells. 
#' Nature. \doi{10.1038/s41586-018-0414-6}
#'
#' @docType package
#' @name scVeloR-package
#' @aliases scVeloR
#'
#' @useDynLib scVeloR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix sparseMatrix drop0 Diagonal
#' @importFrom stats optim lm coef quantile sd median density dist setNames na.omit
#' @importFrom utils head tail
#' @importFrom grDevices colorRampPalette
#' @importFrom methods new slot slotNames
NULL

#' @importFrom Matrix sparseMatrix drop0 Diagonal t
NULL

#' @title RNA Velocity Estimation
#' @description Main entry point for RNA velocity computation.
#' @name velocity
NULL

#' Compute RNA Velocity
#'
#' @description Main function to compute RNA velocity.
#' Supports three modes: "deterministic" (steady-state), "stochastic", and "dynamical".
#'
#' @param object Seurat object with spliced and unspliced data
#' @param mode Velocity mode: "deterministic", "stochastic", or "dynamical"
#' @param min_r2 Minimum R-squared for velocity genes (for deterministic/stochastic)
#' @param max_iter Maximum iterations (for dynamical mode)
#' @param n_top_genes Number of top genes for dynamical mode
#' @param fit_scaling Fit scaling parameter (for dynamical)
#' @param n_cores Number of cores for parallel computation
#' @param verbose Print progress
#' @param ... Additional arguments passed to specific velocity functions
#'
#' @return Seurat object with velocity results in misc$scVeloR$velocity
#'
#' @details
#' The function supports three velocity estimation modes:
#'
#' \describe{
#'   \item{deterministic}{Steady-state model assuming du/dt ≈ 0 at equilibrium.
#'     Uses linear regression to estimate gamma = u/s ratio.}
#'   \item{stochastic}{Uses second-order moments (variance, covariance) to
#'     account for transcriptional noise and bursting.}
#'   \item{dynamical}{Full kinetics model using EM algorithm to infer
#'     transcription (alpha), splicing (beta), and degradation (gamma) rates.}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic velocity computation
#' seurat_obj <- velocity(seurat_obj, mode = "deterministic")
#'
#' # Dynamical model
#' seurat_obj <- velocity(seurat_obj, mode = "dynamical", max_iter = 10)
#' }
#'
#' @export
velocity <- function(object,
                     mode = "deterministic",
                     min_r2 = 0.01,
                     max_iter = 10,
                     n_top_genes = 2000,
                     fit_scaling = TRUE,
                     n_cores = 1,
                     verbose = TRUE,
                     ...) {
  
  if (verbose) {
    message(sprintf("Computing RNA velocity (mode: %s)...", mode))
  }
  
  # Ensure moments are computed
  if (is.null(object@misc$scVeloR$Ms)) {
    if (verbose) message("  Computing moments first...")
    object <- moments(object, verbose = verbose)
  }
  
  # Compute velocity based on mode
  if (mode == "deterministic" || mode == "steady_state") {
    object <- velocity_steady_state(
      object = object,
      min_r2 = min_r2,
      n_cores = n_cores,
      verbose = verbose,
      ...
    )
    
  } else if (mode == "stochastic") {
    object <- velocity_stochastic(
      object = object,
      min_r2 = min_r2,
      n_cores = n_cores,
      verbose = verbose,
      ...
    )
    
  } else if (mode == "dynamical") {
    # First run steady-state to identify velocity genes
    if (verbose) message("  Identifying velocity genes...")
    object <- velocity_steady_state(
      object = object,
      min_r2 = min_r2,
      n_cores = n_cores,
      verbose = FALSE
    )
    
    # Run dynamical recovery
    object <- recover_dynamics(
      object = object,
      genes = "velocity_genes",
      n_top_genes = n_top_genes,
      max_iter = max_iter,
      fit_scaling = fit_scaling,
      n_cores = n_cores,
      verbose = verbose,
      ...
    )
    
    # Compute velocity from dynamics
    object <- velocity_from_dynamics(object, ...)
    
  } else {
    stop(sprintf("Unknown velocity mode: '%s'. Use 'deterministic', 'stochastic', or 'dynamical'.", mode))
  }
  
  object
}

#' Full Velocity Pipeline
#'
#' @description Run complete velocity analysis pipeline:
#' preprocessing, velocity computation, embedding projection, and visualization.
#'
#' @param object Seurat object with spliced and unspliced data
#' @param mode Velocity mode
#' @param embedding Embedding to project velocity onto
#' @param n_neighbors Number of neighbors
#' @param min_counts Minimum counts for gene filtering
#' @param min_cells Minimum cells for gene filtering
#' @param n_cores Number of cores
#' @param verbose Print progress
#' @param ... Additional arguments
#'
#' @return Seurat object with full velocity analysis results
#'
#' @examples
#' \dontrun{
#' seurat_obj <- run_velocity(seurat_obj, mode = "dynamical")
#' plot_velocity(seurat_obj)
#' }
#'
#' @export
run_velocity <- function(object,
                         mode = "deterministic",
                         embedding = "umap",
                         n_neighbors = 30,
                         min_counts = 20,
                         min_cells = 10,
                         n_cores = 1,
                         verbose = TRUE,
                         ...) {
  
  if (verbose) {
    message("Running full velocity pipeline...")
  }
  
  # Step 1: Preprocessing
  if (verbose) message("\n[Step 1/5] Preprocessing...")
  object <- prepare_velocity(
    object = object,
    n_neighbors = n_neighbors,
    min_counts = min_counts,
    min_cells = min_cells,
    verbose = verbose
  )
  
  # Step 2: Compute velocity
  if (verbose) message("\n[Step 2/5] Computing velocity...")
  object <- velocity(
    object = object,
    mode = mode,
    n_cores = n_cores,
    verbose = verbose,
    ...
  )
  
  # Step 3: Build velocity graph
  if (verbose) message("\n[Step 3/5] Building velocity graph...")
  object <- velocity_graph(
    object = object,
    n_neighbors = n_neighbors,
    n_cores = n_cores,
    verbose = verbose
  )
  
  # Step 4: Project velocity onto embedding
  if (verbose) message("\n[Step 4/5] Projecting velocity onto embedding...")
  
  # Check if embedding exists
  emb <- get_embedding(object, embedding)
  if (is.null(emb)) {
    warning(sprintf("Embedding '%s' not found. Skipping velocity embedding.", embedding))
  } else {
    object <- velocity_embedding(
      object = object,
      embedding_name = embedding,
      verbose = verbose
    )
  }
  
  # Step 5: Compute latent time if using dynamical mode
  if (mode == "dynamical") {
    if (verbose) message("\n[Step 5/5] Computing latent time...")
    object <- compute_latent_time(object, verbose = verbose)
  } else {
    if (verbose) message("\n[Step 5/5] Skipping latent time (not dynamical mode)")
  }
  
  if (verbose) {
    message("\n✓ Velocity pipeline complete!")
    message("  Use plot_velocity() or plot_velocity_stream() to visualize results.")
  }
  
  object
}

#' Get Velocity Summary
#'
#' @description Print summary of velocity analysis results.
#'
#' @param object Seurat object with velocity computed
#'
#' @return Invisible NULL, prints summary
#' @export
velocity_summary <- function(object) {
  
  cat("\n=== scVeloR Velocity Analysis Summary ===\n\n")
  
  if (is.null(object@misc$scVeloR)) {
    cat("No velocity analysis found. Run velocity() first.\n")
    return(invisible(NULL))
  }
  
  scVeloR <- object@misc$scVeloR
  
  # Preprocessing
  cat("Preprocessing:\n")
  if (!is.null(scVeloR$Ms)) {
    cat(sprintf("  - Moments: Ms (%d x %d), Mu (%d x %d)\n",
                nrow(scVeloR$Ms), ncol(scVeloR$Ms),
                nrow(scVeloR$Mu), ncol(scVeloR$Mu)))
  }
  if (!is.null(scVeloR$neighbors)) {
    cat(sprintf("  - Neighbors: %d neighbors computed\n",
                ncol(scVeloR$neighbors$indices)))
  }
  
  # Velocity
  cat("\nVelocity:\n")
  if (!is.null(scVeloR$velocity)) {
    vel <- scVeloR$velocity
    cat(sprintf("  - Mode: %s\n", vel$mode))
    cat(sprintf("  - Velocity genes: %d\n", length(vel$velocity_genes)))
    
    if (!is.null(vel$r2)) {
      r2_valid <- vel$r2[!is.na(vel$r2)]
      cat(sprintf("  - R² range: [%.3f, %.3f], median: %.3f\n",
                  min(r2_valid), max(r2_valid), median(r2_valid)))
    }
  } else {
    cat("  - Not computed\n")
  }
  
  # Dynamics
  if (!is.null(scVeloR$dynamics)) {
    cat("\nDynamical model:\n")
    params <- scVeloR$dynamics$gene_params
    n_fitted <- sum(!is.na(params$alpha))
    cat(sprintf("  - Fitted genes: %d\n", n_fitted))
    
    if (n_fitted > 0) {
      lik_valid <- params$likelihood[!is.na(params$likelihood)]
      cat(sprintf("  - Likelihood range: [%.4f, %.4f]\n",
                  min(lik_valid), max(lik_valid)))
    }
  }
  
  # Velocity graph
  if (!is.null(scVeloR$velocity_graph)) {
    cat("\nVelocity graph:\n")
    cat(sprintf("  - Graph computed: TRUE\n"))
    cat(sprintf("  - Transition matrix: %d x %d\n",
                nrow(scVeloR$velocity_graph$transition_matrix),
                ncol(scVeloR$velocity_graph$transition_matrix)))
  }
  
  # Velocity embedding
  if (!is.null(scVeloR$velocity_embedding)) {
    cat("\nVelocity embedding:\n")
    emb_names <- names(scVeloR$velocity_embedding)
    emb_names <- emb_names[emb_names != "current"]
    cat(sprintf("  - Embeddings: %s\n", paste(emb_names, collapse = ", ")))
  }
  
  # Latent time
  if ("latent_time" %in% colnames(object@meta.data)) {
    cat("\nLatent time:\n")
    lt <- object@meta.data$latent_time
    cat(sprintf("  - Range: [%.3f, %.3f]\n", min(lt, na.rm = TRUE), max(lt, na.rm = TRUE)))
  }
  
  cat("\n")
  invisible(NULL)
}

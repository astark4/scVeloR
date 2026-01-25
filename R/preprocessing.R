#' @title Preprocessing Functions
#' @description Functions for preprocessing single-cell data for velocity analysis.
#' @name preprocessing
NULL

#' Filter Genes for Velocity Analysis
#'
#' @description Filter genes based on minimum counts and expression in cells.
#'
#' @param seurat_obj A Seurat object with spliced and unspliced layers
#' @param min_counts Minimum total counts per gene (default: 10)
#' @param min_cells Minimum number of cells expressing the gene (default: 10)
#' @param min_counts_u Minimum unspliced counts per gene (default: 5)
#' @param min_cells_u Minimum cells with unspliced expression (default: 5)
#' @param spliced_layer Name of spliced layer (default: "spliced")
#' @param unspliced_layer Name of unspliced layer (default: "unspliced")
#' @param retain_genes Genes to always retain regardless of filtering
#'
#' @return Modified Seurat object with velocity_genes in var metadata
#'
#' @export
filter_genes <- function(seurat_obj,
                          min_counts = 10,
                          min_cells = 10,
                          min_counts_u = 5,
                          min_cells_u = 5,
                          spliced_layer = "spliced",
                          unspliced_layer = "unspliced",
                          retain_genes = NULL) {
  
  # Extract data
  data <- extract_velocity_data(seurat_obj, spliced_layer, unspliced_layer)
  spliced <- data$spliced
  unspliced <- data$unspliced
  genes <- data$genes
  
  n_genes_initial <- length(genes)
  
  # Calculate statistics
  spliced_counts <- Matrix::rowSums(spliced)
  unspliced_counts <- Matrix::rowSums(unspliced)
  spliced_cells <- Matrix::rowSums(spliced > 0)
  unspliced_cells <- Matrix::rowSums(unspliced > 0)
  
  # Apply filters
  pass_filter <- (spliced_counts >= min_counts) &
                 (spliced_cells >= min_cells) &
                 (unspliced_counts >= min_counts_u) &
                 (unspliced_cells >= min_cells_u)
  
  # Retain specified genes
  if (!is.null(retain_genes)) {
    retain_idx <- genes %in% retain_genes
    pass_filter <- pass_filter | retain_idx
  }
  
  n_genes_pass <- sum(pass_filter)
  
  # Store in metadata
  version <- .get_seurat_version(seurat_obj)
  assay_name <- SeuratObject::DefaultAssay(seurat_obj)
  
  # Add velocity_genes to var metadata
  if (version == "v5") {
    # V5: Use features metadata
    tryCatch({
      seurat_obj[[assay_name]]@meta.data$velocity_genes <- pass_filter
    }, error = function(e) {
      seurat_obj@misc$velocity_genes <- setNames(pass_filter, genes)
    })
  } else {
    # V4: Add to assay meta.features
    tryCatch({
      seurat_obj[[assay_name]]@meta.features$velocity_genes <- pass_filter
    }, error = function(e) {
      seurat_obj@misc$velocity_genes <- setNames(pass_filter, genes)
    })
  }
  
  .vmessage("Filtered genes: ", n_genes_pass, " / ", n_genes_initial, " passed")
  
  return(seurat_obj)
}

#' Normalize Layers
#'
#' @description Normalize spliced and unspliced counts per cell.
#'
#' @param seurat_obj A Seurat object
#' @param target_sum Target sum for normalization (default: NULL uses median)
#' @param log_transform Whether to log-transform (default: FALSE for velocity)
#' @param spliced_layer Name of spliced layer (default: "spliced")
#' @param unspliced_layer Name of unspliced layer (default: "unspliced")
#'
#' @return Modified Seurat object with normalized layers
#'
#' @export
normalize_layers <- function(seurat_obj,
                              target_sum = NULL,
                              log_transform = FALSE,
                              spliced_layer = "spliced",
                              unspliced_layer = "unspliced") {
  
  # Extract data
  data <- extract_velocity_data(seurat_obj, spliced_layer, unspliced_layer)
  spliced <- data$spliced
  unspliced <- data$unspliced
  
  # Store initial sizes
  initial_size_spliced <- Matrix::colSums(spliced)
  initial_size_unspliced <- Matrix::colSums(unspliced)
  
  seurat_obj@meta.data$initial_size_spliced <- initial_size_spliced
  seurat_obj@meta.data$initial_size_unspliced <- initial_size_unspliced
  seurat_obj@meta.data$initial_size <- initial_size_spliced + initial_size_unspliced
  
  # Determine target sum
  if (is.null(target_sum)) {
    target_sum <- stats::median(initial_size_spliced + initial_size_unspliced)
  }
  
  # Normalize spliced
  spliced_norm <- normalize_per_cell(spliced, target_sum, log_transform)
  
  # Normalize unspliced
  unspliced_norm <- normalize_per_cell(unspliced, target_sum, log_transform)
  
  # Store normalized layers
  seurat_obj <- .set_layer_data(seurat_obj, spliced_layer, spliced_norm)
  seurat_obj <- .set_layer_data(seurat_obj, unspliced_layer, unspliced_norm)
  
  .vmessage("Normalized ", ncol(spliced), " cells to target sum ", round(target_sum, 2))
  
  return(seurat_obj)
}

#' Normalize Per Cell Helper
#'
#' @param x Sparse matrix (genes x cells)
#' @param target_sum Target sum per cell
#' @param log_transform Whether to log-transform
#' @return Normalized matrix
#' @keywords internal
normalize_per_cell <- function(x, target_sum, log_transform = FALSE) {
  cell_sums <- Matrix::colSums(x)
  cell_sums[cell_sums == 0] <- 1  # Avoid division by zero
  
  # Normalize
  x_norm <- t(t(x) / cell_sums * target_sum)
  
  # Log transform if requested
  if (log_transform) {
    x_norm@x <- log1p(x_norm@x)
  }
  
  return(x_norm)
}

#' Compute First-Order Moments
#'
#' @description Compute smoothed expression values (first-order moments) using 
#' the neighbor graph.
#'
#' @param seurat_obj A Seurat object with computed neighbors
#' @param n_neighbors Number of neighbors (default: 30)
#' @param mode Use "connectivities" or "distances" (default: "connectivities")
#' @param spliced_layer Name of spliced layer (default: "spliced")
#' @param unspliced_layer Name of unspliced layer (default: "unspliced")
#'
#' @return Modified Seurat object with Ms and Mu layers
#'
#' @export
compute_moments <- function(seurat_obj,
                             n_neighbors = 30,
                             mode = "connectivities",
                             spliced_layer = "spliced",
                             unspliced_layer = "unspliced") {
  
  .vmessage("Computing moments based on ", mode)
  
  # Check if neighbors exist
  if (is.null(seurat_obj@misc$neighbors)) {
    .vmessage("Computing neighbors first...")
    seurat_obj <- compute_neighbors(seurat_obj, n_neighbors = n_neighbors)
  }
  
  # Get connectivities
  conn <- get_connectivities(seurat_obj, mode = mode, n_neighbors = n_neighbors)
  
  # Ensure row-normalized
  conn <- normalize_matrix(conn, axis = 1)
  
  # Extract data
  data <- extract_velocity_data(seurat_obj, spliced_layer, unspliced_layer)
  spliced <- data$spliced
  unspliced <- data$unspliced
  
  # Compute moments using C++ if sparse
  if (inherits(spliced, "sparseMatrix") && inherits(conn, "sparseMatrix")) {
    # Use Rcpp function
    Ms_dense <- compute_moments_sparse_cpp(
      as(spliced, "dgCMatrix"),
      as(conn, "dgCMatrix")
    )
    Mu_dense <- compute_moments_sparse_cpp(
      as(unspliced, "dgCMatrix"),
      as(conn, "dgCMatrix")
    )
    
    Ms <- as(Ms_dense, "CsparseMatrix")
    Mu <- as(Mu_dense, "CsparseMatrix")
  } else {
    # Dense computation
    spliced_dense <- make_dense(spliced)
    unspliced_dense <- make_dense(unspliced)
    conn_dense <- make_dense(conn)
    
    Ms <- as(spliced_dense %*% t(conn_dense), "CsparseMatrix")
    Mu <- as(unspliced_dense %*% t(conn_dense), "CsparseMatrix")
  }
  
  # Set row and column names
  rownames(Ms) <- rownames(spliced)
  colnames(Ms) <- colnames(spliced)
  rownames(Mu) <- rownames(unspliced)
  colnames(Mu) <- colnames(unspliced)
  
  # Store moments
  seurat_obj <- .set_layer_data(seurat_obj, "Ms", Ms)
  seurat_obj <- .set_layer_data(seurat_obj, "Mu", Mu)
  
  .vmessage("Added Ms and Mu to layers")
  
  return(seurat_obj)
}

#' Compute Second-Order Moments
#'
#' @description Compute second-order moments for stochastic velocity estimation.
#'
#' @param seurat_obj A Seurat object with computed neighbors
#' @param spliced_layer Name of spliced layer (default: "spliced")
#' @param unspliced_layer Name of unspliced layer (default: "unspliced")
#' @param adjusted Whether to compute adjusted moments (default: FALSE)
#'
#' @return List with Mss and Mus matrices
#'
#' @keywords internal
compute_second_order_moments <- function(seurat_obj,
                                          spliced_layer = "spliced",
                                          unspliced_layer = "unspliced",
                                          adjusted = FALSE) {
  
  # Get connectivities
  conn <- get_connectivities(seurat_obj)
  
  # Extract data
  data <- extract_velocity_data(seurat_obj, spliced_layer, unspliced_layer)
  spliced <- data$spliced
  unspliced <- data$unspliced
  
  if (inherits(spliced, "sparseMatrix") && inherits(conn, "sparseMatrix")) {
    # Use C++ functions
    Mss <- compute_second_moments_sparse_cpp(
      as(spliced, "dgCMatrix"),
      as(conn, "dgCMatrix")
    )
    Mus <- compute_cross_moments_sparse_cpp(
      as(unspliced, "dgCMatrix"),
      as(spliced, "dgCMatrix"),
      as(conn, "dgCMatrix")
    )
  } else {
    # Dense computation
    s <- make_dense(spliced)
    u <- make_dense(unspliced)
    conn_dense <- make_dense(conn)
    
    Mss <- (s * s) %*% t(conn_dense)
    Mus <- (s * u) %*% t(conn_dense)
  }
  
  if (adjusted) {
    Ms <- .get_layer_data(seurat_obj, "Ms")
    Mu <- .get_layer_data(seurat_obj, "Mu")
    Mss <- 2 * Mss - make_dense(Ms)
    Mus <- 2 * Mus - make_dense(Mu)
  }
  
  list(Mss = Mss, Mus = Mus)
}

#' Filter and Normalize (One-Step Preprocessing)
#'
#' @description Perform all preprocessing steps in one call.
#'
#' @param seurat_obj A Seurat object with spliced and unspliced layers
#' @param min_counts Minimum total counts per gene (default: 10)
#' @param min_cells Minimum number of cells expressing the gene (default: 10)
#' @param min_counts_u Minimum unspliced counts per gene (default: 5)
#' @param min_cells_u Minimum cells with unspliced expression (default: 5)
#' @param n_neighbors Number of neighbors for moment computation (default: 30)
#' @param target_sum Target sum for normalization (default: NULL uses median)
#' @param spliced_layer Name of spliced layer (default: "spliced")
#' @param unspliced_layer Name of unspliced layer (default: "unspliced")
#'
#' @return Preprocessed Seurat object
#'
#' @export
#' @examples
#' \dontrun{
#' seurat_obj <- filter_and_normalize(seurat_obj)
#' }
filter_and_normalize <- function(seurat_obj,
                                   min_counts = 10,
                                   min_cells = 10,
                                   min_counts_u = 5,
                                   min_cells_u = 5,
                                   n_neighbors = 30,
                                   target_sum = NULL,
                                   spliced_layer = "spliced",
                                   unspliced_layer = "unspliced") {
  
  .vmessage("Preprocessing data...")
  
  # Store proportions
  show_proportions(seurat_obj, spliced_layer, unspliced_layer)
  
  # Filter genes
  seurat_obj <- filter_genes(seurat_obj,
                              min_counts = min_counts,
                              min_cells = min_cells,
                              min_counts_u = min_counts_u,
                              min_cells_u = min_cells_u,
                              spliced_layer = spliced_layer,
                              unspliced_layer = unspliced_layer)
  
  # Normalize
  seurat_obj <- normalize_layers(seurat_obj,
                                  target_sum = target_sum,
                                  spliced_layer = spliced_layer,
                                  unspliced_layer = unspliced_layer)
  
  .vmessage("Preprocessing complete")
  
  return(seurat_obj)
}

#' Check if Data is Already Normalized
#'
#' @param x Sparse matrix
#' @return Logical
#' @keywords internal
is_normalized <- function(x) {
  # Check if values look normalized (most values < 100)
  if (inherits(x, "sparseMatrix")) {
    sample_vals <- x@x[seq(1, min(length(x@x), 10000))]
  } else {
    sample_vals <- as.vector(x[seq(1, min(length(x), 10000))])
  }
  
  # If most values are small floats, probably normalized
  max_val <- max(sample_vals, na.rm = TRUE)
  median_val <- stats::median(sample_vals[sample_vals > 0], na.rm = TRUE)
  
  # Heuristic: raw counts tend to be integers and can be large
  is_float <- any(sample_vals != floor(sample_vals))
  is_small <- median_val < 10
  
  return(is_float && is_small)
}

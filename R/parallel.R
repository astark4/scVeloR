#' @title Parallel Computing Support
#' @description Functions for parallel computation support.
#' @name parallel
NULL

#' Set Parallel Plan
#'
#' @description Configure parallel computing backend.
#'
#' @param n_workers Number of worker processes (default: NULL for auto-detect)
#' @param type Parallel type: "multisession", "multicore", or "sequential"
#'
#' @export
#' @examples
#' \dontrun{
#' set_parallel(n_workers = 4)
#' # ... run computations ...
#' set_parallel(type = "sequential")  # Reset
#' }
set_parallel <- function(n_workers = NULL, type = "multisession") {
  
  if (type == "sequential") {
    future::plan(future::sequential)
    .vmessage("Set to sequential execution")
    return(invisible(NULL))
  }
  
  if (is.null(n_workers)) {
    n_workers <- max(1, parallel::detectCores() - 1)
  }
  
  if (type == "multicore") {
    if (.Platform$OS.type == "windows") {
      warning("multicore not supported on Windows, using multisession")
      type <- "multisession"
    }
  }
  
  if (type == "multisession") {
    future::plan(future::multisession, workers = n_workers)
  } else if (type == "multicore") {
    future::plan(future::multicore, workers = n_workers)
  }
  
  .vmessage("Set to ", type, " with ", n_workers, " workers")
  
  invisible(n_workers)
}

#' Get Current Parallel Status
#'
#' @description Get information about current parallel configuration.
#'
#' @return List with parallel status information
#'
#' @export
parallel_status <- function() {
  plan <- future::plan()
  
  # Get number of workers
  n_workers <- tryCatch({
    future::nbrOfWorkers()
  }, error = function(e) 1L)
  
  # Check if parallel
  is_parallel <- n_workers > 1
  
  list(
    plan = class(plan)[1],
    workers = n_workers,
    is_parallel = is_parallel
  )
}

#' Parallel Apply with Progress
#'
#' @description Apply function in parallel with progress reporting.
#'
#' @param X Vector or list to iterate over
#' @param FUN Function to apply
#' @param ... Additional arguments to FUN
#' @param n_jobs Number of jobs (default: from options)
#' @param chunk_size Chunk size for batching (default: NULL for auto)
#'
#' @return List of results
#'
#' @export
parallel_apply <- function(X, FUN, ..., n_jobs = NULL, chunk_size = NULL) {
  
  n_jobs <- get_n_jobs(n_jobs)
  n_items <- length(X)
  
  if (n_jobs == 1L || n_items < 10) {
    # Sequential with progress
    return(apply_with_progress(X, FUN, ..., n_jobs = 1L, show_progress = TRUE))
  }
  
  # Parallel execution
  if (is.null(chunk_size)) {
    chunk_size <- max(1, ceiling(n_items / (n_jobs * 4)))
  }
  
  # Set up parallel
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  
  if (!inherits(old_plan, "sequential")) {
    # Already have parallel plan
  } else {
    future::plan(future::multisession, workers = n_jobs)
  }
  
  # Run with progress
  result <- apply_with_progress(X, FUN, ..., n_jobs = n_jobs, show_progress = TRUE)
  
  return(result)
}

#' Batch Process Genes
#'
#' @description Process genes in batches for memory efficiency.
#'
#' @param genes Gene names to process
#' @param FUN Function to apply to each gene
#' @param batch_size Batch size (default: 100)
#' @param n_jobs Number of parallel jobs
#' @param ... Additional arguments to FUN
#'
#' @return List of results
#' @keywords internal
batch_process_genes <- function(genes, FUN, batch_size = 100, n_jobs = 1L, ...) {
  
  n_genes <- length(genes)
  n_batches <- ceiling(n_genes / batch_size)
  
  all_results <- list()
  
  for (i in seq_len(n_batches)) {
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, n_genes)
    batch_genes <- genes[start_idx:end_idx]
    
    .vmessage("Processing batch ", i, "/", n_batches, " (", length(batch_genes), " genes)")
    
    batch_results <- parallel_apply(batch_genes, FUN, ..., n_jobs = n_jobs)
    
    all_results <- c(all_results, batch_results)
  }
  
  names(all_results) <- genes
  return(all_results)
}

#' Check Memory Usage
#'
#' @description Check current memory usage.
#'
#' @return Memory usage in MB
#' @keywords internal
check_memory <- function() {
  if (.Platform$OS.type == "unix") {
    # Unix-like systems
    gc_result <- gc(verbose = FALSE)
    mem_used <- sum(gc_result[, 2])  # MB
  } else {
    # Windows
    mem_used <- memory.size()
  }
  
  return(mem_used)
}

#' Estimate Memory Requirements
#'
#' @description Estimate memory needed for computation.
#'
#' @param n_cells Number of cells
#' @param n_genes Number of genes
#' @param operation Type of operation
#'
#' @return Estimated memory in MB
#' @keywords internal
estimate_memory <- function(n_cells, n_genes, operation = "velocity_graph") {
  
  # Base memory per cell-gene (8 bytes for double)
  bytes_per_element <- 8
  
  if (operation == "velocity_graph") {
    # Need: X matrix + V matrix + graph (~n_cells^2 * sparsity)
    matrix_size <- n_cells * n_genes * bytes_per_element * 2
    graph_size <- n_cells * n_cells * 0.1 * bytes_per_element  # 10% sparsity
    total <- matrix_size + graph_size
  } else if (operation == "moments") {
    # Need: expression + connectivity + result
    matrix_size <- n_cells * n_genes * bytes_per_element * 2
    conn_size <- n_cells * 30 * bytes_per_element  # ~30 neighbors
    total <- matrix_size + conn_size
  } else {
    total <- n_cells * n_genes * bytes_per_element
  }
  
  # Convert to MB with safety margin
  mb <- total / (1024^2) * 1.5
  
  return(round(mb, 1))
}

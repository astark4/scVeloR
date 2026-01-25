#' @title Velocity Embedding Projection
#' @description Project velocity vectors onto low-dimensional embeddings (UMAP, tSNE, etc.).
#' @name velocity_embedding
NULL

#' Project Velocity onto Embedding
#'
#' @description Project high-dimensional velocity vectors onto 2D/3D embeddings
#' using the velocity graph transition matrix.
#'
#' @param object Seurat object with velocity computed
#' @param embedding_name Name of embedding to use (e.g., "umap", "tsne")
#' @param sigma Kernel bandwidth for transition smoothing
#' @param scale Scale factor for arrow lengths
#' @param direct Use direct projection instead of transition-based
#' @param autoscale Automatically scale velocities
#' @param verbose Print progress
#'
#' @return Seurat object with velocity_embedding in misc$scVeloR
#' @export
velocity_embedding <- function(object,
                                embedding_name = "umap",
                                sigma = 0.05,
                                scale = 1,
                                direct = FALSE,
                                autoscale = TRUE,
                                verbose = TRUE) {
  
  if (is.null(object@misc$scVeloR$velocity$velocity_s)) {
    stop("Run velocity() first")
  }
  
  if (verbose) {
    message("Projecting velocity onto embedding...")
  }
  
  # Get embedding coordinates
  embedding <- get_embedding(object, embedding_name)
  
  if (is.null(embedding)) {
    stop(sprintf("Embedding '%s' not found", embedding_name))
  }
  
  n_cells <- nrow(embedding)
  n_dims <- ncol(embedding)
  
  if (direct) {
    # Direct projection method
    velocity_embedded <- project_velocity_direct(
      object = object,
      embedding = embedding,
      verbose = verbose
    )
  } else {
    # Transition matrix-based projection
    # Get transition matrix
    T <- NULL
    
    if (!is.null(object@misc$scVeloR$velocity$transition_matrix)) {
      T <- object@misc$scVeloR$velocity$transition_matrix
    } else if (!is.null(object@misc$scVeloR$velocity_graph$transition_matrix)) {
      T <- object@misc$scVeloR$velocity_graph$transition_matrix
    }
    
    if (is.null(T)) {
      # Build velocity graph first
      object <- velocity_graph(object, verbose = FALSE)
      T <- object@misc$scVeloR$velocity_graph$transition_matrix
    }
    
    velocity_embedded <- project_velocity_transition(
      embedding = embedding,
      transition_matrix = T,
      sigma = sigma
    )
  }
  
  # Scale velocities
  if (autoscale) {
    # Scale by median displacement
    embedding_diffs <- diff(as.matrix(embedding[order(object@meta.data$latent_time %||% 1:n_cells), ]))
    median_displacement <- median(sqrt(rowSums(embedding_diffs^2)))
    
    velocity_norms <- sqrt(rowSums(velocity_embedded^2))
    median_velocity <- median(velocity_norms[velocity_norms > 0])
    
    if (median_velocity > 0 && is.finite(median_velocity)) {
      scale_factor <- median_displacement / median_velocity * scale
      velocity_embedded <- velocity_embedded * scale_factor
    }
  } else {
    velocity_embedded <- velocity_embedded * scale
  }
  
  # Store results
  if (is.null(object@misc$scVeloR$velocity_embedding)) {
    object@misc$scVeloR$velocity_embedding <- list()
  }
  
  colnames(velocity_embedded) <- paste0("V", 1:n_dims)
  rownames(velocity_embedded) <- rownames(embedding)
  
  object@misc$scVeloR$velocity_embedding[[embedding_name]] <- velocity_embedded
  object@misc$scVeloR$velocity_embedding$current <- embedding_name
  
  if (verbose) {
    message(sprintf("Done. Velocity embedding stored for '%s'", embedding_name))
  }
  
  object
}

#' Project Velocity Using Transition Matrix
#'
#' @description Project velocity using the transition probability matrix.
#' For each cell, the velocity is the expected displacement based on transition probabilities.
#'
#' @param embedding Embedding coordinates
#' @param transition_matrix Transition probability matrix
#' @param sigma Kernel bandwidth
#'
#' @return Matrix of velocity vectors (n_cells x n_dims)
#' @keywords internal
project_velocity_transition <- function(embedding,
                                         transition_matrix,
                                         sigma = 0.05) {
  
  n_cells <- nrow(embedding)
  n_dims <- ncol(embedding)
  
  velocity_embedded <- matrix(0, n_cells, n_dims)
  
  # Make transition matrix sparse if not already
  if (!inherits(transition_matrix, "sparseMatrix")) {
    transition_matrix <- as(transition_matrix, "sparseMatrix")
  }
  
  # Get current cell positions
  X <- as.matrix(embedding)
  
  # For each cell, compute expected displacement
  for (i in seq_len(n_cells)) {
    # Get transition probabilities to neighbors
    probs <- transition_matrix[i, ]
    
    # Find non-zero transitions
    nonzero <- which(probs > 0)
    
    if (length(nonzero) == 0) next
    
    # Get neighbor positions
    X_neighbors <- X[nonzero, , drop = FALSE]
    p_neighbors <- probs[nonzero]
    
    # Remove self-transition if present
    if (i %in% nonzero) {
      self_idx <- which(nonzero == i)
      X_neighbors <- X_neighbors[-self_idx, , drop = FALSE]
      p_neighbors <- p_neighbors[-self_idx]
    }
    
    if (length(p_neighbors) == 0 || sum(p_neighbors) == 0) next
    
    # Normalize probabilities
    p_neighbors <- p_neighbors / sum(p_neighbors)
    
    # Compute displacements
    dX <- X_neighbors - matrix(X[i, ], nrow = nrow(X_neighbors), ncol = n_dims, byrow = TRUE)
    
    # Weight by probability
    velocity_embedded[i, ] <- colSums(dX * p_neighbors)
  }
  
  velocity_embedded
}

#' Direct Velocity Projection
#'
#' @description Project velocity directly using correlation-based method.
#'
#' @param object Seurat object
#' @param embedding Embedding coordinates
#' @param verbose Print progress
#'
#' @return Matrix of velocity vectors
#' @keywords internal
project_velocity_direct <- function(object,
                                     embedding,
                                     verbose = FALSE) {
  
  # Get expression and velocity data
  data <- get_velocity_data(object)
  Ms <- data$Ms
  
  velocity <- object@misc$scVeloR$velocity$velocity_s
  velocity_genes <- object@misc$scVeloR$velocity$velocity_genes
  
  Ms <- Ms[, velocity_genes, drop = FALSE]
  velocity <- velocity[, velocity_genes, drop = FALSE]
  
  n_cells <- nrow(Ms)
  n_dims <- ncol(embedding)
  
  velocity_embedded <- matrix(0, n_cells, n_dims)
  
  # Get neighbors
  if (!is.null(object@misc$scVeloR$neighbors)) {
    indices <- object@misc$scVeloR$neighbors$indices
    
    for (i in seq_len(n_cells)) {
      v <- velocity[i, ]
      v_norm <- sqrt(sum(v^2))
      
      if (v_norm == 0 || !is.finite(v_norm)) next
      
      neighbors <- indices[i, ]
      neighbors <- neighbors[neighbors > 0 & neighbors <= n_cells & neighbors != i]
      
      if (length(neighbors) == 0) next
      
      # Compute correlations with displacements
      correlations <- numeric(length(neighbors))
      displacements <- matrix(0, length(neighbors), n_dims)
      
      for (k in seq_along(neighbors)) {
        j <- neighbors[k]
        
        # Expression displacement
        dx <- Ms[j, ] - Ms[i, ]
        dx_norm <- sqrt(sum(dx^2))
        
        if (dx_norm > 0 && is.finite(dx_norm)) {
          correlations[k] <- sum(v * dx) / (v_norm * dx_norm)
        }
        
        # Embedding displacement
        displacements[k, ] <- embedding[j, ] - embedding[i, ]
      }
      
      # Weight displacements by positive correlations
      pos_corr <- pmax(correlations, 0)
      
      if (sum(pos_corr) > 0) {
        weights <- pos_corr / sum(pos_corr)
        velocity_embedded[i, ] <- colSums(displacements * weights)
      }
    }
  }
  
  velocity_embedded
}

#' Get Embedding from Seurat Object
#'
#' @description Extract embedding coordinates from Seurat object.
#'
#' @param object Seurat object
#' @param embedding_name Name of embedding
#'
#' @return Embedding matrix (cells x dims)
#' @keywords internal
get_embedding <- function(object, embedding_name = "umap") {
  
  # Try to find the embedding
  reduction_name <- tolower(embedding_name)
  
  if (reduction_name %in% names(object@reductions)) {
    emb <- Seurat::Embeddings(object, reduction = reduction_name)
    return(as.matrix(emb))
  }
  
  # Try common variations
  alternatives <- c(
    embedding_name,
    tolower(embedding_name),
    toupper(embedding_name),
    paste0(embedding_name, "_"),
    paste0("X_", embedding_name)
  )
  
  for (alt in alternatives) {
    if (alt %in% names(object@reductions)) {
      emb <- Seurat::Embeddings(object, reduction = alt)
      return(as.matrix(emb))
    }
  }
  
  # Not found
  NULL
}

#' Compute Grid Velocity Embedding
#'
#' @description Compute velocity vectors on a regular grid for visualization.
#'
#' @param object Seurat object with velocity embedding
#' @param embedding_name Name of embedding
#' @param n_grid Number of grid points per dimension
#' @param min_cells Minimum cells per grid point
#' @param smooth Smoothing factor
#'
#' @return List with grid coordinates and velocities
#' @export
velocity_grid <- function(object,
                          embedding_name = "umap",
                          n_grid = 50,
                          min_cells = 5,
                          smooth = 0.5) {
  
  # Get embedding
  embedding <- get_embedding(object, embedding_name)
  
  if (is.null(embedding)) {
    stop(sprintf("Embedding '%s' not found", embedding_name))
  }
  
  # Get velocity embedding
  if (is.null(object@misc$scVeloR$velocity_embedding[[embedding_name]])) {
    stop("Run velocity_embedding() first")
  }
  
  velocity <- object@misc$scVeloR$velocity_embedding[[embedding_name]]
  
  n_dims <- ncol(embedding)
  
  # Create grid
  x_range <- range(embedding[, 1])
  y_range <- range(embedding[, 2])
  
  x_grid <- seq(x_range[1], x_range[2], length.out = n_grid)
  y_grid <- seq(y_range[1], y_range[2], length.out = n_grid)
  
  grid <- expand.grid(x = x_grid, y = y_grid)
  
  # Compute grid velocities using Gaussian kernel smoothing
  grid_velocity <- matrix(0, nrow(grid), n_dims)
  grid_density <- numeric(nrow(grid))
  
  bandwidth <- smooth * c(diff(x_range), diff(y_range)) / n_grid
  
  for (g in seq_len(nrow(grid))) {
    gx <- grid$x[g]
    gy <- grid$y[g]
    
    # Distances to grid point
    dx <- embedding[, 1] - gx
    dy <- embedding[, 2] - gy
    
    # Gaussian weights
    weights <- exp(-0.5 * ((dx / bandwidth[1])^2 + (dy / bandwidth[2])^2))
    
    if (sum(weights) > 0 && sum(weights > 0.01) >= min_cells) {
      weights <- weights / sum(weights)
      
      grid_velocity[g, ] <- colSums(velocity * weights)
      grid_density[g] <- sum(weights > 0.01)
    }
  }
  
  # Filter out low-density grid points
  valid <- grid_density >= min_cells
  
  list(
    x = grid$x[valid],
    y = grid$y[valid],
    vx = grid_velocity[valid, 1],
    vy = grid_velocity[valid, 2],
    density = grid_density[valid]
  )
}

#' Compute Streamline Data
#'
#' @description Compute streamlines for velocity stream plot.
#'
#' @param object Seurat object with velocity embedding
#' @param embedding_name Name of embedding
#' @param n_grid Grid size for streamline computation
#' @param density Density of streamlines
#' @param min_length Minimum streamline length
#' @param max_length Maximum streamline length
#'
#' @return List with streamline data
#' @export
velocity_streamlines <- function(object,
                                  embedding_name = "umap",
                                  n_grid = 100,
                                  density = 1,
                                  min_length = 0.01,
                                  max_length = 4) {
  
  # Get embedding
  embedding <- get_embedding(object, embedding_name)
  
  if (is.null(embedding)) {
    stop(sprintf("Embedding '%s' not found", embedding_name))
  }
  
  # Get velocity embedding
  if (is.null(object@misc$scVeloR$velocity_embedding[[embedding_name]])) {
    stop("Run velocity_embedding() first")
  }
  
  velocity <- object@misc$scVeloR$velocity_embedding[[embedding_name]]
  
  # Create velocity field on grid
  x_range <- range(embedding[, 1])
  y_range <- range(embedding[, 2])
  
  x_grid <- seq(x_range[1], x_range[2], length.out = n_grid)
  y_grid <- seq(y_range[1], y_range[2], length.out = n_grid)
  
  # Interpolate velocity field
  vx_field <- matrix(0, n_grid, n_grid)
  vy_field <- matrix(0, n_grid, n_grid)
  
  bandwidth <- 0.1 * c(diff(x_range), diff(y_range)) / n_grid
  
  for (i in seq_len(n_grid)) {
    for (j in seq_len(n_grid)) {
      gx <- x_grid[i]
      gy <- y_grid[j]
      
      dx <- embedding[, 1] - gx
      dy <- embedding[, 2] - gy
      
      weights <- exp(-0.5 * ((dx / bandwidth[1])^2 + (dy / bandwidth[2])^2))
      
      if (sum(weights) > 0) {
        weights <- weights / sum(weights)
        vx_field[i, j] <- sum(velocity[, 1] * weights)
        vy_field[i, j] <- sum(velocity[, 2] * weights)
      }
    }
  }
  
  # Generate seed points for streamlines
  n_seeds <- round(n_grid * density / 2)
  seed_x <- runif(n_seeds, x_range[1], x_range[2])
  seed_y <- runif(n_seeds, y_range[1], y_range[2])
  
  # Trace streamlines
  streamlines <- list()
  dt <- min(diff(x_grid)[1], diff(y_grid)[1]) * 0.5
  max_steps <- round(max_length / dt)
  
  for (s in seq_len(n_seeds)) {
    x <- seed_x[s]
    y <- seed_y[s]
    
    line_x <- numeric(max_steps)
    line_y <- numeric(max_steps)
    
    for (step in seq_len(max_steps)) {
      line_x[step] <- x
      line_y[step] <- y
      
      # Interpolate velocity at current position
      i <- which.min(abs(x_grid - x))
      j <- which.min(abs(y_grid - y))
      
      if (i < 1 || i > n_grid || j < 1 || j > n_grid) break
      
      vx <- vx_field[i, j]
      vy <- vy_field[i, j]
      
      v_norm <- sqrt(vx^2 + vy^2)
      if (v_norm < 1e-10) break
      
      # Normalize velocity and take step
      x <- x + vx / v_norm * dt
      y <- y + vy / v_norm * dt
      
      # Check bounds
      if (x < x_range[1] || x > x_range[2] || y < y_range[1] || y > y_range[2]) break
    }
    
    # Store streamline if long enough
    valid_steps <- which(line_x != 0 | line_y != 0)
    if (length(valid_steps) > 0) {
      line_length <- sum(sqrt(diff(line_x[valid_steps])^2 + diff(line_y[valid_steps])^2))
      
      if (line_length >= min_length * max(diff(x_range), diff(y_range))) {
        streamlines[[length(streamlines) + 1]] <- list(
          x = line_x[valid_steps],
          y = line_y[valid_steps]
        )
      }
    }
  }
  
  list(
    streamlines = streamlines,
    x_range = x_range,
    y_range = y_range
  )
}

#' Null coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

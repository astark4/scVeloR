#' @title Stream Plot
#' @description Functions for creating streamline visualizations of velocity.
#' @name plot_stream
NULL

#' Velocity Stream Plot
#'
#' @description Create a streamline plot showing velocity flow on embedding.
#'
#' @param seurat_obj A Seurat object with velocity embedding
#' @param basis Embedding to use (default: "umap")
#' @param color_by Variable for coloring cells
#' @param density Stream density (default: 1)
#' @param smooth Smoothing factor for streamlines (default: 0.5)
#' @param n_grid Number of grid points for interpolation
#' @param arrow_size Size of streamline arrows
#' @param alpha Stream transparency
#' @param linewidth Width of streamlines
#' @param palette Color palette
#' @param title Plot title
#'
#' @return ggplot object
#'
#' @export
velocity_stream_plot <- function(seurat_obj,
                                   basis = "umap",
                                   color_by = NULL,
                                   density = 1,
                                   smooth = 0.5,
                                   n_grid = 50,
                                   arrow_size = 0.3,
                                   alpha = 0.7,
                                   linewidth = 0.5,
                                   palette = NULL,
                                   title = NULL) {
  
  # Get embedding and velocity
  emb <- get_embedding(seurat_obj, basis)
  V_emb <- get_velocity_embedding(seurat_obj)
  
  n_cells <- nrow(emb)
  
  # Create grid
  x_range <- range(emb[, 1])
  y_range <- range(emb[, 2])
  
  # Expand range slightly
  x_margin <- diff(x_range) * 0.05
  y_margin <- diff(y_range) * 0.05
  
  x_seq <- seq(x_range[1] - x_margin, x_range[2] + x_margin, length.out = n_grid)
  y_seq <- seq(y_range[1] - y_margin, y_range[2] + y_margin, length.out = n_grid)
  
  # Interpolate velocity field onto grid
  grid_vx <- matrix(0, n_grid, n_grid)
  grid_vy <- matrix(0, n_grid, n_grid)
  grid_mass <- matrix(0, n_grid, n_grid)
  
  # Gaussian kernel interpolation
  sigma <- smooth * min(diff(x_seq)[1], diff(y_seq)[1]) * 5
  
  for (k in seq_len(n_cells)) {
    xi <- emb[k, 1]
    yi <- emb[k, 2]
    vxi <- V_emb[k, 1]
    vyi <- V_emb[k, 2]
    
    if (is.na(vxi) || is.na(vyi)) next
    
    # Find nearby grid points
    dx <- (x_seq - xi)^2
    dy <- (y_seq - yi)^2
    
    for (i in seq_len(n_grid)) {
      for (j in seq_len(n_grid)) {
        dist_sq <- dx[i] + dy[j]
        weight <- exp(-dist_sq / (2 * sigma^2))
        
        if (weight > 0.01) {
          grid_vx[i, j] <- grid_vx[i, j] + weight * vxi
          grid_vy[i, j] <- grid_vy[i, j] + weight * vyi
          grid_mass[i, j] <- grid_mass[i, j] + weight
        }
      }
    }
  }
  
  # Normalize by mass
  grid_mass[grid_mass == 0] <- 1
  grid_vx <- grid_vx / grid_mass
  grid_vy <- grid_vy / grid_mass
  
  # Generate streamlines
  streamlines <- generate_streamlines(x_seq, y_seq, grid_vx, grid_vy, 
                                       n_lines = round(50 * density),
                                       step_size = 0.3)
  
  # Create point data for background
  point_data <- data.frame(
    x = emb[, 1],
    y = emb[, 2]
  )
  
  if (!is.null(color_by) && color_by %in% colnames(seurat_obj@meta.data)) {
    point_data$color <- seurat_obj@meta.data[[color_by]]
  }
  
  # Create plot
  p <- ggplot2::ggplot() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "right"
    ) +
    ggplot2::labs(
      x = paste0(toupper(basis), "_1"),
      y = paste0(toupper(basis), "_2"),
      title = title
    )
  
  # Add background points
  if (!is.null(color_by) && "color" %in% names(point_data)) {
    p <- p + ggplot2::geom_point(
      data = point_data,
      ggplot2::aes(x = x, y = y, color = color),
      size = 0.5,
      alpha = 0.3
    )
  } else {
    p <- p + ggplot2::geom_point(
      data = point_data,
      ggplot2::aes(x = x, y = y),
      color = "grey80",
      size = 0.5,
      alpha = 0.3
    )
  }
  
  # Add streamlines
  if (nrow(streamlines) > 0) {
    p <- p + ggplot2::geom_path(
      data = streamlines,
      ggplot2::aes(x = x, y = y, group = line_id),
      color = "black",
      alpha = alpha,
      linewidth = linewidth,
      arrow = ggplot2::arrow(
        length = ggplot2::unit(arrow_size, "lines"),
        ends = "last",
        type = "closed"
      )
    )
  }
  
  return(p)
}

#' Generate Streamlines
#'
#' @description Generate streamline paths from velocity field.
#'
#' @param x_seq X grid coordinates
#' @param y_seq Y grid coordinates
#' @param vx X velocity field
#' @param vy Y velocity field
#' @param n_lines Number of streamlines
#' @param step_size Integration step size
#' @param n_steps Maximum steps per line
#'
#' @return Data frame with streamline coordinates
#' @keywords internal
generate_streamlines <- function(x_seq, y_seq, vx, vy, 
                                   n_lines = 50,
                                   step_size = 0.3,
                                   n_steps = 100) {
  
  n_grid <- length(x_seq)
  dx <- x_seq[2] - x_seq[1]
  dy <- y_seq[2] - y_seq[1]
  
  # Random starting points
  x_range <- range(x_seq)
  y_range <- range(y_seq)
  
  starts_x <- stats::runif(n_lines, x_range[1], x_range[2])
  starts_y <- stats::runif(n_lines, y_range[1], y_range[2])
  
  all_lines <- list()
  
  for (line_id in seq_len(n_lines)) {
    x <- starts_x[line_id]
    y <- starts_y[line_id]
    
    line_x <- numeric(n_steps)
    line_y <- numeric(n_steps)
    line_x[1] <- x
    line_y[1] <- y
    
    for (step in 2:n_steps) {
      # Interpolate velocity at current position
      ix <- max(1, min(n_grid - 1, floor((x - x_seq[1]) / dx) + 1))
      iy <- max(1, min(n_grid - 1, floor((y - y_seq[1]) / dy) + 1))
      
      # Bilinear interpolation
      fx <- (x - x_seq[ix]) / dx
      fy <- (y - y_seq[iy]) / dy
      
      vx_interp <- (1-fx)*(1-fy)*vx[ix,iy] + fx*(1-fy)*vx[ix+1,iy] +
                   (1-fx)*fy*vx[ix,iy+1] + fx*fy*vx[ix+1,iy+1]
      vy_interp <- (1-fx)*(1-fy)*vy[ix,iy] + fx*(1-fy)*vy[ix+1,iy] +
                   (1-fx)*fy*vy[ix,iy+1] + fx*fy*vy[ix+1,iy+1]
      
      # Check for valid velocity
      mag <- sqrt(vx_interp^2 + vy_interp^2)
      if (is.na(mag) || mag < 1e-10) break
      
      # Normalize and step
      x <- x + step_size * dx * vx_interp / mag
      y <- y + step_size * dy * vy_interp / mag
      
      # Check bounds
      if (x < x_range[1] || x > x_range[2] || 
          y < y_range[1] || y > y_range[2]) break
      
      line_x[step] <- x
      line_y[step] <- y
    }
    
    # Remove unused steps
    valid <- line_x != 0 | line_y != 0 | seq_along(line_x) == 1
    valid[1] <- TRUE
    
    if (sum(valid) > 5) {
      all_lines[[length(all_lines) + 1]] <- data.frame(
        x = line_x[valid],
        y = line_y[valid],
        line_id = line_id
      )
    }
  }
  
  if (length(all_lines) == 0) {
    return(data.frame(x = numeric(0), y = numeric(0), line_id = integer(0)))
  }
  
  do.call(rbind, all_lines)
}

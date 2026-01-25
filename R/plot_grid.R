#' @title Grid Plot
#' @description Functions for creating grid-based velocity visualizations.
#' @name plot_grid
NULL

#' Velocity Grid Plot
#'
#' @description Create a grid-aggregated velocity arrow plot.
#'
#' @param seurat_obj A Seurat object with velocity embedding
#' @param basis Embedding to use (default: "umap")
#' @param color_by Variable for coloring cells
#' @param n_grid Number of grid cells per dimension (default: 40)
#' @param min_mass Minimum cells per grid cell (default: 1)
#' @param arrow_size Arrow size multiplier (default: 1)
#' @param arrow_color Color for arrows (default: "black")
#' @param smooth Smoothing factor (default: 0.5)
#' @param alpha Point transparency (default: 0.3)
#' @param palette Color palette
#' @param title Plot title
#'
#' @return ggplot object
#'
#' @export
velocity_grid_plot <- function(seurat_obj,
                                 basis = "umap",
                                 color_by = NULL,
                                 n_grid = 40,
                                 min_mass = 1,
                                 arrow_size = 1,
                                 arrow_color = "black",
                                 smooth = 0.5,
                                 alpha = 0.3,
                                 palette = NULL,
                                 title = NULL) {
  
  # Get embedding and velocity
  emb <- get_embedding(seurat_obj, basis)
  V_emb <- get_velocity_embedding(seurat_obj)
  
  # Create grid data
  grid_data <- create_arrow_grid(emb, V_emb, n_grid = n_grid, min_mass = min_mass)
  
  # Apply smoothing to grid velocities if requested
  if (smooth > 0 && nrow(grid_data) > 3) {
    grid_data <- smooth_grid_velocity(grid_data, smooth)
  }
  
  # Scale arrows
  if (nrow(grid_data) > 0) {
    mag <- sqrt(grid_data$vx^2 + grid_data$vy^2)
    median_mag <- stats::median(mag[mag > 0])
    if (median_mag > 0) {
      scale_factor <- 1 / median_mag * (max(emb[,1]) - min(emb[,1])) / n_grid * 0.8
      grid_data$vx <- grid_data$vx * scale_factor * arrow_size
      grid_data$vy <- grid_data$vy * scale_factor * arrow_size
    }
  }
  
  # Create point data
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
      alpha = alpha
    )
    
    if (is.numeric(point_data$color)) {
      p <- p + ggplot2::scale_color_viridis_c()
    }
    p <- p + ggplot2::labs(color = color_by)
  } else {
    p <- p + ggplot2::geom_point(
      data = point_data,
      ggplot2::aes(x = x, y = y),
      color = "grey80",
      size = 0.5,
      alpha = alpha
    )
  }
  
  # Add grid arrows
  if (nrow(grid_data) > 0) {
    p <- p + ggplot2::geom_segment(
      data = grid_data,
      ggplot2::aes(x = x, y = y, xend = x + vx, yend = y + vy),
      arrow = ggplot2::arrow(
        length = ggplot2::unit(0.1 * arrow_size, "inches"),
        type = "closed"
      ),
      color = arrow_color,
      linewidth = 0.5 * arrow_size,
      alpha = 0.8
    )
  }
  
  return(p)
}

#' Smooth Grid Velocity
#'
#' @description Apply spatial smoothing to grid velocity field.
#'
#' @param grid_data Grid data frame with x, y, vx, vy
#' @param smooth Smoothing factor
#' @return Smoothed grid data
#' @keywords internal
smooth_grid_velocity <- function(grid_data, smooth) {
  
  n <- nrow(grid_data)
  if (n < 3) return(grid_data)
  
  # Estimate kernel width from grid spacing
  x_spacing <- min(diff(sort(unique(grid_data$x))))
  y_spacing <- min(diff(sort(unique(grid_data$y))))
  sigma <- smooth * min(x_spacing, y_spacing) * 2
  
  vx_smooth <- numeric(n)
  vy_smooth <- numeric(n)
  
  for (i in seq_len(n)) {
    xi <- grid_data$x[i]
    yi <- grid_data$y[i]
    
    # Compute weights
    dist_sq <- (grid_data$x - xi)^2 + (grid_data$y - yi)^2
    weights <- exp(-dist_sq / (2 * sigma^2))
    weights <- weights / sum(weights)
    
    # Weighted average
    vx_smooth[i] <- sum(weights * grid_data$vx)
    vy_smooth[i] <- sum(weights * grid_data$vy)
  }
  
  grid_data$vx <- vx_smooth
  grid_data$vy <- vy_smooth
  
  return(grid_data)
}

#' Plot Velocity Field
#'
#' @description Plot the full velocity vector field on a grid.
#'
#' @param seurat_obj A Seurat object
#' @param basis Embedding basis
#' @param n_grid Grid resolution
#' @param color_by Color vectors by magnitude or direction
#'
#' @return ggplot object
#'
#' @export
velocity_field_plot <- function(seurat_obj,
                                  basis = "umap",
                                  n_grid = 30,
                                  color_by = "magnitude") {
  
  emb <- get_embedding(seurat_obj, basis)
  V_emb <- get_velocity_embedding(seurat_obj)
  
  # Create dense grid
  x_range <- range(emb[, 1])
  y_range <- range(emb[, 2])
  
  x_seq <- seq(x_range[1], x_range[2], length.out = n_grid)
  y_seq <- seq(y_range[1], y_range[2], length.out = n_grid)
  
  grid_expand <- expand.grid(x = x_seq, y = y_seq)
  
  # Interpolate velocities
  grid_expand$vx <- 0
  grid_expand$vy <- 0
  
  sigma <- diff(x_seq)[1] * 2
  
  for (k in seq_len(nrow(emb))) {
    dist_sq <- (grid_expand$x - emb[k, 1])^2 + (grid_expand$y - emb[k, 2])^2
    weights <- exp(-dist_sq / (2 * sigma^2))
    
    grid_expand$vx <- grid_expand$vx + weights * V_emb[k, 1]
    grid_expand$vy <- grid_expand$vy + weights * V_emb[k, 2]
  }
  
  # Normalize
  grid_expand$magnitude <- sqrt(grid_expand$vx^2 + grid_expand$vy^2)
  max_mag <- max(grid_expand$magnitude)
  if (max_mag > 0) {
    grid_expand$vx <- grid_expand$vx / max_mag * diff(x_range) / n_grid * 0.8
    grid_expand$vy <- grid_expand$vy / max_mag * diff(y_range) / n_grid * 0.8
  }
  
  # Compute direction for coloring
  grid_expand$direction <- atan2(grid_expand$vy, grid_expand$vx)
  
  # Plot
  p <- ggplot2::ggplot(grid_expand, ggplot2::aes(x = x, y = y)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank()
    )
  
  if (color_by == "magnitude") {
    p <- p + ggplot2::geom_segment(
      ggplot2::aes(xend = x + vx, yend = y + vy, color = magnitude),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "inches"), type = "closed")
    ) +
    ggplot2::scale_color_viridis_c()
  } else {
    p <- p + ggplot2::geom_segment(
      ggplot2::aes(xend = x + vx, yend = y + vy, color = direction),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.05, "inches"), type = "closed")
    ) +
    ggplot2::scale_color_gradientn(
      colors = grDevices::rainbow(12),
      limits = c(-pi, pi)
    )
  }
  
  return(p)
}

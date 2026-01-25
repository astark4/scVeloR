#' @title Velocity Embedding Plot
#' @description Functions for visualizing RNA velocity on embeddings.
#' @name plot_velocity
NULL

#' Velocity Plot
#'
#' @description Main plotting function for velocity embedding visualization.
#'
#' @param seurat_obj A Seurat object with velocity embedding
#' @param basis Embedding to use (default: "umap")
#' @param color_by Variable for coloring points (column name in metadata)
#' @param arrows Whether to show velocity arrows (default: TRUE)
#' @param arrow_size Size of arrows (default: 1)
#' @param arrow_length Length multiplier for arrows (default: 1)
#' @param alpha Point transparency (default: 0.5)
#' @param point_size Point size (default: 1)
#' @param density Arrow density (default: 1, use values < 1 for fewer arrows)
#' @param n_sample Number of cells to sample for arrows (alternative to density)
#' @param palette Color palette name or vector
#' @param title Plot title
#' @param ... Additional arguments to ggplot
#'
#' @return ggplot object
#'
#' @export
velocity_plot <- function(seurat_obj,
                           basis = "umap",
                           color_by = NULL,
                           arrows = TRUE,
                           arrow_size = 1,
                           arrow_length = 1,
                           alpha = 0.5,
                           point_size = 1,
                           density = 1,
                           n_sample = NULL,
                           palette = NULL,
                           title = NULL,
                           ...) {
  
  # Get embedding
  emb <- get_embedding(seurat_obj, basis)
  n_cells <- nrow(emb)
  
  # Create data frame
  plot_data <- data.frame(
    x = emb[, 1],
    y = emb[, 2],
    row.names = rownames(emb)
  )
  
  # Add color variable
  if (!is.null(color_by)) {
    if (color_by %in% colnames(seurat_obj@meta.data)) {
      plot_data$color <- seurat_obj@meta.data[[color_by]]
    } else {
      warning("Color variable '", color_by, "' not found")
      plot_data$color <- "cell"
    }
  } else {
    plot_data$color <- "cell"
  }
  
  # Get velocity embedding if arrows requested
  if (arrows) {
    if (is.null(seurat_obj@misc$velocity_embedding)) {
      warning("Velocity embedding not found. Run project_velocity_embedding first.")
      arrows <- FALSE
    } else {
      V_emb <- get_velocity_embedding(seurat_obj)
      plot_data$vx <- V_emb[, 1] * arrow_length
      plot_data$vy <- V_emb[, 2] * arrow_length
    }
  }
  
  # Sample cells for arrows
  if (arrows && !is.null(density) && density < 1) {
    n_sample <- round(n_cells * density)
  }
  
  if (arrows && !is.null(n_sample) && n_sample < n_cells) {
    arrow_idx <- sample(seq_len(n_cells), n_sample)
    arrow_data <- plot_data[arrow_idx, ]
  } else if (arrows) {
    arrow_data <- plot_data
  }
  
  # Create base plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y)) +
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
  
  # Add points
  if (!is.null(color_by) && color_by != "cell") {
    p <- p + ggplot2::geom_point(
      ggplot2::aes(color = color),
      size = point_size,
      alpha = alpha
    )
    
    # Add color scale
    if (is.numeric(plot_data$color)) {
      if (!is.null(palette)) {
        p <- p + ggplot2::scale_color_gradientn(colors = palette)
      } else {
        p <- p + ggplot2::scale_color_viridis_c()
      }
    } else {
      if (!is.null(palette)) {
        p <- p + ggplot2::scale_color_manual(values = palette)
      }
    }
    p <- p + ggplot2::labs(color = color_by)
  } else {
    p <- p + ggplot2::geom_point(
      color = "grey70",
      size = point_size,
      alpha = alpha
    )
  }
  
  # Add arrows
  if (arrows) {
    p <- p + ggplot2::geom_segment(
      data = arrow_data,
      ggplot2::aes(
        xend = x + vx,
        yend = y + vy
      ),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.05 * arrow_size, "inches"),
                              type = "closed"),
      size = 0.3 * arrow_size,
      alpha = 0.7
    )
  }
  
  return(p)
}

#' Velocity Embedding Plot (Alias)
#'
#' @description Plot velocity arrows on embedding. Alias for velocity_plot.
#'
#' @inheritParams velocity_plot
#' @return ggplot object
#'
#' @export
velocity_embedding_plot <- function(seurat_obj,
                                     basis = "umap",
                                     color_by = NULL,
                                     arrow_size = 1,
                                     arrow_length = 1,
                                     alpha = 0.5,
                                     point_size = 1,
                                     density = 1,
                                     palette = NULL,
                                     title = NULL,
                                     ...) {
  
  velocity_plot(
    seurat_obj,
    basis = basis,
    color_by = color_by,
    arrows = TRUE,
    arrow_size = arrow_size,
    arrow_length = arrow_length,
    alpha = alpha,
    point_size = point_size,
    density = density,
    palette = palette,
    title = title,
    ...
  )
}

#' Plot Velocity by Gene
#'
#' @description Plot velocity for a specific gene showing spliced vs unspliced.
#'
#' @param seurat_obj A Seurat object
#' @param gene Gene name
#' @param color_by Variable for coloring
#' @param fit Show fit line (default: TRUE)
#'
#' @return ggplot object
#'
#' @export
plot_velocity_gene <- function(seurat_obj, gene, color_by = NULL, fit = TRUE) {
  
  # Get moments
  Ms <- make_dense(.get_layer_data(seurat_obj, "Ms"))
  Mu <- make_dense(.get_layer_data(seurat_obj, "Mu"))
  
  if (!gene %in% rownames(Ms)) {
    stop("Gene '", gene, "' not found", call. = FALSE)
  }
  
  plot_data <- data.frame(
    spliced = Ms[gene, ],
    unspliced = Mu[gene, ]
  )
  
  if (!is.null(color_by) && color_by %in% colnames(seurat_obj@meta.data)) {
    plot_data$color <- seurat_obj@meta.data[[color_by]]
  }
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = spliced, y = unspliced)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Spliced",
      y = "Unspliced",
      title = gene
    )
  
  if (!is.null(color_by) && "color" %in% names(plot_data)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = color), alpha = 0.5)
  } else {
    p <- p + ggplot2::geom_point(color = "steelblue", alpha = 0.5)
  }
  
  # Add fit line
  if (fit) {
    fit_params <- get_fit_params(seurat_obj)
    if (!is.null(fit_params) && gene %in% rownames(fit_params)) {
      gamma <- fit_params[gene, "fit_gamma"]
      offset <- fit_params[gene, "fit_offset"]
      if (is.na(offset)) offset <- 0
      
      if (!is.na(gamma)) {
        p <- p + ggplot2::geom_abline(
          slope = gamma,
          intercept = offset,
          color = "red",
          linetype = "dashed"
        )
      }
    }
  }
  
  return(p)
}

#' Create Arrow Grid Data
#'
#' @description Create aggregated arrow data on a grid for cleaner visualization.
#'
#' @param emb Embedding coordinates
#' @param V_emb Velocity embedding
#' @param n_grid Number of grid cells per dimension
#' @param min_mass Minimum mass (number of cells) per grid cell
#'
#' @return Data frame with grid arrow coordinates
#' @keywords internal
create_arrow_grid <- function(emb, V_emb, n_grid = 40, min_mass = 1) {
  
  # Create grid
  x_range <- range(emb[, 1])
  y_range <- range(emb[, 2])
  
  x_breaks <- seq(x_range[1], x_range[2], length.out = n_grid + 1)
  y_breaks <- seq(y_range[1], y_range[2], length.out = n_grid + 1)
  
  # Assign cells to grid
  x_bin <- cut(emb[, 1], breaks = x_breaks, labels = FALSE, include.lowest = TRUE)
  y_bin <- cut(emb[, 2], breaks = y_breaks, labels = FALSE, include.lowest = TRUE)
  
  # Aggregate velocities per grid cell
  grid_data <- list()
  
  for (i in seq_len(n_grid)) {
    for (j in seq_len(n_grid)) {
      mask <- x_bin == i & y_bin == j
      n_cells <- sum(mask, na.rm = TRUE)
      
      if (n_cells >= min_mass) {
        grid_data[[length(grid_data) + 1]] <- data.frame(
          x = mean(emb[mask, 1]),
          y = mean(emb[mask, 2]),
          vx = mean(V_emb[mask, 1]),
          vy = mean(V_emb[mask, 2]),
          mass = n_cells
        )
      }
    }
  }
  
  if (length(grid_data) == 0) {
    return(data.frame())
  }
  
  do.call(rbind, grid_data)
}

#' Get Color Palette
#'
#' @description Get a color palette for plotting.
#'
#' @param n Number of colors needed
#' @param name Palette name
#' @return Vector of colors
#' @keywords internal
get_palette <- function(n, name = "default") {
  
  if (name == "default" || name == "scanpy") {
    # Default scanpy-like palette
    cols <- c(
      "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
      "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
      "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
      "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"
    )
  } else if (name == "viridis") {
    cols <- grDevices::colorRampPalette(
      c("#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
        "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725")
    )(n)
  } else {
    cols <- grDevices::rainbow(n)
  }
  
  if (n <= length(cols)) {
    return(cols[seq_len(n)])
  } else {
    return(grDevices::colorRampPalette(cols)(n))
  }
}

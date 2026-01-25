#' @title Plotting Utilities
#' @description Utility functions for plotting.
#' @name plot_utils
NULL

#' Default Velocity Theme
#'
#' @description Get the default ggplot theme for velocity plots.
#'
#' @param base_size Base font size
#' @param axis_text Whether to show axis text
#'
#' @return ggplot theme object
#'
#' @export
theme_velocity <- function(base_size = 11, axis_text = FALSE) {
  theme <- ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "right",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  if (!axis_text) {
    theme <- theme + ggplot2::theme(
      axis.text = ggplot2::element_blank()
    )
  }
  
  return(theme)
}

#' Combine Multiple Velocity Plots
#'
#' @description Combine multiple ggplot objects into one figure.
#'
#' @param ... ggplot objects to combine
#' @param ncol Number of columns
#' @param nrow Number of rows
#' @param widths Relative widths
#' @param heights Relative heights
#'
#' @return Combined plot
#'
#' @export
combine_plots <- function(..., ncol = NULL, nrow = NULL, 
                           widths = NULL, heights = NULL) {
  
  plots <- list(...)
  n_plots <- length(plots)
  
  if (is.null(ncol) && is.null(nrow)) {
    ncol <- ceiling(sqrt(n_plots))
    nrow <- ceiling(n_plots / ncol)
  } else if (is.null(ncol)) {
    ncol <- ceiling(n_plots / nrow)
  } else if (is.null(nrow)) {
    nrow <- ceiling(n_plots / ncol)
  }
  
  # Use patchwork if available, otherwise grid.arrange
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combined <- patchwork::wrap_plots(plots, ncol = ncol, nrow = nrow)
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    combined <- gridExtra::grid.arrange(grobs = plots, ncol = ncol, nrow = nrow)
  } else {
    warning("Install patchwork or gridExtra for plot combination")
    return(plots[[1]])
  }
  
  return(combined)
}

#' Save Velocity Plot
#'
#' @description Save a velocity plot to file.
#'
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Plot width in inches
#' @param height Plot height in inches
#' @param dpi Resolution
#' @param ... Additional arguments to ggsave
#'
#' @export
save_velocity_plot <- function(plot, filename, 
                                 width = 8, height = 6, 
                                 dpi = 300, ...) {
  
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    ...
  )
  
  .vmessage("Saved plot to ", filename)
}

#' Get Embedding Limits
#'
#' @description Get axis limits for embedding plot.
#'
#' @param emb Embedding matrix
#' @param expand Expansion factor
#'
#' @return List with xlim and ylim
#' @keywords internal
get_embedding_limits <- function(emb, expand = 0.05) {
  x_range <- range(emb[, 1], na.rm = TRUE)
  y_range <- range(emb[, 2], na.rm = TRUE)
  
  x_expand <- diff(x_range) * expand
  y_expand <- diff(y_range) * expand
  
  list(
    xlim = c(x_range[1] - x_expand, x_range[2] + x_expand),
    ylim = c(y_range[1] - y_expand, y_range[2] + y_expand)
  )
}

#' Create Density Mask for Arrows
#'
#' @description Create a mask for subsampling arrows based on local density.
#'
#' @param emb Embedding coordinates
#' @param density Target density (fraction of cells)
#' @param method "random" or "grid"
#'
#' @return Logical vector
#' @keywords internal
density_sample <- function(emb, density = 1, method = "random") {
  n <- nrow(emb)
  
  if (density >= 1) {
    return(rep(TRUE, n))
  }
  
  n_sample <- max(1, round(n * density))
  
  if (method == "random") {
    idx <- sample(seq_len(n), n_sample)
    mask <- seq_len(n) %in% idx
  } else {
    # Grid-based sampling
    n_grid <- round(sqrt(n_sample))
    x_breaks <- seq(min(emb[,1]), max(emb[,1]), length.out = n_grid + 1)
    y_breaks <- seq(min(emb[,2]), max(emb[,2]), length.out = n_grid + 1)
    
    x_bin <- cut(emb[,1], x_breaks, labels = FALSE, include.lowest = TRUE)
    y_bin <- cut(emb[,2], y_breaks, labels = FALSE, include.lowest = TRUE)
    
    # Sample one cell per bin
    mask <- rep(FALSE, n)
    for (i in seq_len(n_grid)) {
      for (j in seq_len(n_grid)) {
        in_bin <- which(x_bin == i & y_bin == j)
        if (length(in_bin) > 0) {
          selected <- sample(in_bin, 1)
          mask[selected] <- TRUE
        }
      }
    }
  }
  
  return(mask)
}

#' Color Scale for Velocity
#'
#' @description Create a diverging color scale centered at zero.
#'
#' @param n Number of colors
#' @param low Color for negative values
#' @param mid Color for zero
#' @param high Color for positive values
#'
#' @return Vector of colors
#' @keywords internal
velocity_colors <- function(n = 100, 
                             low = "#313695", 
                             mid = "#FFFFBF", 
                             high = "#A50026") {
  grDevices::colorRampPalette(c(low, mid, high))(n)
}

#' Plot Velocity Metrics
#'
#' @description Create a summary plot of velocity quality metrics.
#'
#' @param seurat_obj A Seurat object with velocity
#'
#' @return ggplot object
#'
#' @export
plot_velocity_metrics <- function(seurat_obj) {
  
  # Collect metrics
  metrics <- data.frame(row.names = colnames(seurat_obj))
  
  # Velocity confidence
  if ("velocity_confidence" %in% colnames(seurat_obj@meta.data)) {
    metrics$confidence <- seurat_obj@meta.data$velocity_confidence
  }
  
  # Velocity length
  if (!is.null(seurat_obj@misc$velocity_embedding)) {
    V_emb <- seurat_obj@misc$velocity_embedding
    metrics$velocity_length <- sqrt(rowSums(V_emb^2))
  }
  
  # Latent time
  if ("latent_time" %in% colnames(seurat_obj@meta.data)) {
    metrics$latent_time <- seurat_obj@meta.data$latent_time
  }
  
  # Convert to long format
  metrics$cell <- rownames(metrics)
  plot_data <- stats::reshape(
    metrics,
    direction = "long",
    varying = setdiff(names(metrics), "cell"),
    v.names = "value",
    timevar = "metric",
    times = setdiff(names(metrics), "cell")
  )
  
  # Create violin plots
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = metric, y = value)) +
    ggplot2::geom_violin(fill = "steelblue", alpha = 0.6) +
    ggplot2::geom_boxplot(width = 0.1, fill = "white") +
    ggplot2::theme_minimal() +
    ggplot2::labs(x = "Metric", y = "Value", title = "Velocity Quality Metrics") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  
  return(p)
}

#' @title Heatmap and Scatter Plots
#' @description Functions for creating heatmaps and scatter plots of velocity data.
#' @name plot_heatmap
NULL

#' Velocity Heatmap
#'
#' @description Create a heatmap of velocity values for selected genes.
#'
#' @param seurat_obj A Seurat object with velocity
#' @param genes Genes to include (default: top velocity genes)
#' @param n_genes Number of genes if not specified (default: 50)
#' @param cluster_rows Cluster rows (genes) (default: TRUE)
#' @param cluster_cols Cluster columns (cells) (default: FALSE)
#' @param order_by Variable to order cells by (default: NULL)
#' @param group_by Variable to group cells (for annotation)
#' @param scale Scale values ("row", "column", "none")
#' @param colors Color palette
#' @param show_rownames Show gene names (default: TRUE)
#' @param show_colnames Show cell names (default: FALSE)
#'
#' @return ggplot object
#'
#' @export
velocity_heatmap <- function(seurat_obj,
                              genes = NULL,
                              n_genes = 50,
                              cluster_rows = TRUE,
                              cluster_cols = FALSE,
                              order_by = NULL,
                              group_by = NULL,
                              scale = "row",
                              colors = NULL,
                              show_rownames = TRUE,
                              show_colnames = FALSE) {
  
  # Get velocity
  velocity <- get_velocity(seurat_obj, remove_na = TRUE)
  
  # Select genes
  if (is.null(genes)) {
    # Use top genes by variance
    gene_var <- apply(velocity, 1, stats::var, na.rm = TRUE)
    genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(n_genes, length(gene_var)))]
  }
  
  genes <- intersect(genes, rownames(velocity))
  velocity <- velocity[genes, , drop = FALSE]
  
  # Convert to dense matrix
  mat <- make_dense(velocity)
  
  # Scale
  if (scale == "row") {
    mat <- t(scale(t(mat)))
  } else if (scale == "column") {
    mat <- scale(mat)
  }
  
  # Handle NAs from scaling
  mat[!is.finite(mat)] <- 0
  
  # Order cells
  if (!is.null(order_by)) {
    if (order_by %in% colnames(seurat_obj@meta.data)) {
      order_var <- seurat_obj@meta.data[[order_by]]
      cell_order <- order(order_var)
      mat <- mat[, cell_order, drop = FALSE]
    } else if (order_by == "pseudotime" && "velocity_pseudotime" %in% colnames(seurat_obj@meta.data)) {
      order_var <- seurat_obj@meta.data$velocity_pseudotime
      cell_order <- order(order_var)
      mat <- mat[, cell_order, drop = FALSE]
    }
  } else if (cluster_cols) {
    hc <- stats::hclust(stats::dist(t(mat)))
    mat <- mat[, hc$order, drop = FALSE]
  }
  
  # Cluster rows
  if (cluster_rows && nrow(mat) > 2) {
    hc <- stats::hclust(stats::dist(mat))
    mat <- mat[hc$order, , drop = FALSE]
  }
  
  # Convert to long format for ggplot
  plot_data <- data.frame(
    gene = rep(rownames(mat), ncol(mat)),
    cell = rep(colnames(mat), each = nrow(mat)),
    value = as.vector(mat)
  )
  
  # Maintain order
  plot_data$gene <- factor(plot_data$gene, levels = rev(rownames(mat)))
  plot_data$cell <- factor(plot_data$cell, levels = colnames(mat))
  
  # Color palette
  if (is.null(colors)) {
    colors <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8",
                "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026")
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = cell, y = gene, fill = value)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(
      colors = colors,
      limits = c(-max(abs(plot_data$value)), max(abs(plot_data$value)))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = if (show_colnames) ggplot2::element_text(angle = 90, hjust = 1, size = 6) 
                    else ggplot2::element_blank(),
      axis.text.y = if (show_rownames) ggplot2::element_text(size = 8) 
                    else ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = "Cells", y = "Genes", fill = "Velocity")
  
  # Add grouping annotation
  if (!is.null(group_by) && group_by %in% colnames(seurat_obj@meta.data)) {
    # This would require additional annotation bar
    # For simplicity, just mention in message
    .vmessage("To add grouping annotations, consider using ComplexHeatmap package")
  }
  
  return(p)
}

#' Velocity Scatter Plot
#'
#' @description Create scatter plots of velocity-related quantities.
#'
#' @param seurat_obj A Seurat object
#' @param x Variable for x-axis
#' @param y Variable for y-axis
#' @param color_by Variable for coloring
#' @param size_by Variable for point sizing
#' @param layer Layer to use for expression values
#' @param alpha Point transparency
#' @param point_size Point size
#'
#' @return ggplot object
#'
#' @export
velocity_scatter <- function(seurat_obj,
                               x,
                               y,
                               color_by = NULL,
                               size_by = NULL,
                               layer = "Ms",
                               alpha = 0.6,
                               point_size = 1) {
  
  # Build plot data
  plot_data <- data.frame(row.names = colnames(seurat_obj))
  
  # Helper to get variable
  get_var <- function(var_name) {
    if (var_name %in% colnames(seurat_obj@meta.data)) {
      return(seurat_obj@meta.data[[var_name]])
    } else if (.layer_exists(seurat_obj, layer)) {
      mat <- .get_layer_data(seurat_obj, layer)
      if (var_name %in% rownames(mat)) {
        return(as.numeric(mat[var_name, ]))
      }
    }
    # Try velocity
    if (.layer_exists(seurat_obj, "velocity")) {
      vel <- .get_layer_data(seurat_obj, "velocity")
      if (var_name %in% rownames(vel)) {
        return(as.numeric(vel[var_name, ]))
      }
    }
    return(NULL)
  }
  
  # Get x and y
  plot_data$x <- get_var(x)
  plot_data$y <- get_var(y)
  
  if (is.null(plot_data$x) || is.null(plot_data$y)) {
    stop("Could not find variables '", x, "' and/or '", y, "'", call. = FALSE)
  }
  
  # Get color
  if (!is.null(color_by)) {
    plot_data$color <- get_var(color_by)
    if (is.null(plot_data$color)) {
      plot_data$color <- seurat_obj@meta.data[[color_by]]
    }
  }
  
  # Get size
  if (!is.null(size_by)) {
    plot_data$size <- get_var(size_by)
  }
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y)) +
    ggplot2::theme_minimal()
  
  # Add points with appropriate aesthetics
  if (!is.null(color_by) && !is.null(size_by)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = color, size = size), alpha = alpha)
  } else if (!is.null(color_by)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = color), size = point_size, alpha = alpha)
  } else if (!is.null(size_by)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(size = size), alpha = alpha)
  } else {
    p <- p + ggplot2::geom_point(color = "steelblue", size = point_size, alpha = alpha)
  }
  
  # Color scale
  if (!is.null(color_by) && "color" %in% names(plot_data)) {
    if (is.numeric(plot_data$color)) {
      p <- p + ggplot2::scale_color_viridis_c()
    }
    p <- p + ggplot2::labs(color = color_by)
  }
  
  p <- p + ggplot2::labs(x = x, y = y)
  
  return(p)
}

#' Plot Phase Portrait
#'
#' @description Plot phase portrait (spliced vs unspliced) for genes.
#'
#' @param seurat_obj A Seurat object
#' @param genes Gene names
#' @param color_by Variable for coloring
#' @param ncol Number of columns in faceted plot
#'
#' @return ggplot object
#'
#' @export
plot_phase_portrait <- function(seurat_obj,
                                  genes,
                                  color_by = NULL,
                                  ncol = 4) {
  
  # Get data
  Ms <- make_dense(.get_layer_data(seurat_obj, "Ms"))
  Mu <- make_dense(.get_layer_data(seurat_obj, "Mu"))
  
  genes <- intersect(genes, rownames(Ms))
  
  if (length(genes) == 0) {
    stop("No valid genes found", call. = FALSE)
  }
  
  # Build plot data
  plot_list <- lapply(genes, function(gene) {
    df <- data.frame(
      spliced = Ms[gene, ],
      unspliced = Mu[gene, ],
      gene = gene
    )
    if (!is.null(color_by) && color_by %in% colnames(seurat_obj@meta.data)) {
      df$color <- seurat_obj@meta.data[[color_by]]
    }
    df
  })
  
  plot_data <- do.call(rbind, plot_list)
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = spliced, y = unspliced)) +
    ggplot2::theme_minimal() +
    ggplot2::facet_wrap(~ gene, scales = "free", ncol = ncol)
  
  if (!is.null(color_by) && "color" %in% names(plot_data)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = color), alpha = 0.5, size = 0.5)
  } else {
    p <- p + ggplot2::geom_point(color = "steelblue", alpha = 0.5, size = 0.5)
  }
  
  # Add fit lines
  fit_params <- get_fit_params(seurat_obj)
  if (!is.null(fit_params)) {
    for (gene in genes) {
      if (gene %in% rownames(fit_params)) {
        gamma <- fit_params[gene, "fit_gamma"]
        offset <- fit_params[gene, "fit_offset"]
        if (is.na(offset)) offset <- 0
        
        if (!is.na(gamma) && gamma > 0) {
          gene_data <- plot_data[plot_data$gene == gene, ]
          x_range <- range(gene_data$spliced, na.rm = TRUE)
          
          line_data <- data.frame(
            spliced = x_range,
            unspliced = gamma * x_range + offset,
            gene = gene
          )
          
          p <- p + ggplot2::geom_line(
            data = line_data,
            color = "red",
            linetype = "dashed"
          )
        }
      }
    }
  }
  
  return(p)
}

#' @title Velocity Visualization
#' @description Functions for visualizing RNA velocity results.
#' @name plot_velocity
NULL

#' Plot Velocity Arrows
#'
#' @description Plot velocity arrows on embedding coordinates.
#'
#' @param object Seurat object with velocity embedding
#' @param embedding Name of embedding to use
#' @param color_by Variable to color cells by (default: cell type if available)
#' @param arrow_scale Scale factor for arrows
#' @param arrow_length Arrow head length
#' @param alpha Transparency
#' @param size Point size
#' @param show_density Show velocity density
#' @param n_arrows Number of arrows to show (subsample if too many cells)
#' @param seed Random seed for subsampling
#' @param title Plot title
#'
#' @return ggplot2 object
#' @export
plot_velocity <- function(object,
                          embedding = "umap",
                          color_by = NULL,
                          arrow_scale = 1,
                          arrow_length = 0.15,
                          alpha = 0.7,
                          size = 0.5,
                          show_density = FALSE,
                          n_arrows = NULL,
                          seed = 42,
                          title = NULL) {
  
  # Get embedding coordinates
  emb <- get_embedding(object, embedding)
  
  if (is.null(emb)) {
    stop(sprintf("Embedding '%s' not found", embedding))
  }
  
  # Get velocity embedding
  vel_emb <- object@misc$scVeloR$velocity_embedding[[embedding]]
  
  if (is.null(vel_emb)) {
    stop("Run velocity_embedding() first")
  }
  
  n_cells <- nrow(emb)
  
  # Create data frame
  df <- data.frame(
    x = emb[, 1],
    y = emb[, 2],
    vx = vel_emb[, 1] * arrow_scale,
    vy = vel_emb[, 2] * arrow_scale
  )
  
  # Add coloring variable
  if (!is.null(color_by)) {
    if (color_by %in% colnames(object@meta.data)) {
      df$color <- object@meta.data[[color_by]]
    } else {
      df$color <- 1  # Default
    }
  } else {
    # Try to find a cluster column
    cluster_cols <- c("seurat_clusters", "clusters", "celltype", "cell_type")
    found <- FALSE
    for (col in cluster_cols) {
      if (col %in% colnames(object@meta.data)) {
        df$color <- object@meta.data[[col]]
        found <- TRUE
        break
      }
    }
    if (!found) {
      df$color <- "cells"
    }
  }
  
  # Subsample if too many cells
  if (!is.null(n_arrows) && n_arrows < n_cells) {
    set.seed(seed)
    arrow_idx <- sample(n_cells, n_arrows)
  } else {
    arrow_idx <- seq_len(n_cells)
  }
  
  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(color = color), size = size, alpha = alpha) +
    ggplot2::geom_segment(
      data = df[arrow_idx, ],
      ggplot2::aes(xend = x + vx, yend = y + vy),
      arrow = ggplot2::arrow(length = ggplot2::unit(arrow_length, "cm"), type = "closed"),
      linewidth = 0.3,
      alpha = 0.6
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = paste0(toupper(embedding), "_1"),
      y = paste0(toupper(embedding), "_2"),
      title = title,
      color = color_by
    ) +
    ggplot2::coord_fixed()
  
  # Color scale
  if (is.numeric(df$color)) {
    p <- p + ggplot2::scale_color_viridis_c()
  } else {
    p <- p + ggplot2::scale_color_manual(values = get_palette(length(unique(df$color))))
  }
  
  p
}

#' Plot Velocity Stream
#'
#' @description Plot velocity streamlines on embedding.
#'
#' @param object Seurat object with velocity embedding
#' @param embedding Name of embedding
#' @param color_by Variable to color cells by
#' @param density Streamline density
#' @param size Point size
#' @param alpha Point transparency
#' @param line_width Streamline width
#' @param line_alpha Streamline transparency
#' @param title Plot title
#'
#' @return ggplot2 object
#' @export
plot_velocity_stream <- function(object,
                                  embedding = "umap",
                                  color_by = NULL,
                                  density = 1,
                                  size = 0.5,
                                  alpha = 0.5,
                                  line_width = 0.5,
                                  line_alpha = 0.8,
                                  title = NULL) {
  
  # Get embedding
  emb <- get_embedding(object, embedding)
  
  if (is.null(emb)) {
    stop(sprintf("Embedding '%s' not found", embedding))
  }
  
  # Get streamlines
  stream_data <- velocity_streamlines(object, embedding, density = density)
  
  # Create cell data frame
  df <- data.frame(
    x = emb[, 1],
    y = emb[, 2]
  )
  
  # Add coloring
  if (!is.null(color_by) && color_by %in% colnames(object@meta.data)) {
    df$color <- object@meta.data[[color_by]]
  } else {
    cluster_cols <- c("seurat_clusters", "clusters", "celltype")
    for (col in cluster_cols) {
      if (col %in% colnames(object@meta.data)) {
        df$color <- object@meta.data[[col]]
        break
      }
    }
    if (!"color" %in% names(df)) df$color <- "cells"
  }
  
  # Create base plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = x, y = y, color = color),
      size = size,
      alpha = alpha
    )
  
  # Add streamlines
  for (i in seq_along(stream_data$streamlines)) {
    sl <- stream_data$streamlines[[i]]
    if (length(sl$x) > 1) {
      stream_df <- data.frame(x = sl$x, y = sl$y, group = i)
      p <- p + ggplot2::geom_path(
        data = stream_df,
        ggplot2::aes(x = x, y = y, group = group),
        linewidth = line_width,
        alpha = line_alpha,
        color = "black",
        arrow = ggplot2::arrow(
          length = ggplot2::unit(0.15, "cm"),
          type = "closed",
          ends = "last"
        )
      )
    }
  }
  
  p <- p +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = paste0(toupper(embedding), "_1"),
      y = paste0(toupper(embedding), "_2"),
      title = title
    ) +
    ggplot2::coord_fixed()
  
  # Color scale
  if (is.numeric(df$color)) {
    p <- p + ggplot2::scale_color_viridis_c()
  } else {
    p <- p + ggplot2::scale_color_manual(values = get_palette(length(unique(df$color))))
  }
  
  p
}

#' Plot Velocity Grid
#'
#' @description Plot velocity arrows on a regular grid.
#'
#' @param object Seurat object with velocity embedding
#' @param embedding Name of embedding
#' @param color_by Variable to color cells by
#' @param n_grid Number of grid points per dimension
#' @param arrow_scale Arrow scale factor
#' @param size Point size
#' @param alpha Point transparency
#' @param title Plot title
#'
#' @return ggplot2 object
#' @export
plot_velocity_grid <- function(object,
                                embedding = "umap",
                                color_by = NULL,
                                n_grid = 40,
                                arrow_scale = 1,
                                size = 0.5,
                                alpha = 0.3,
                                title = NULL) {
  
  # Get embedding
  emb <- get_embedding(object, embedding)
  
  if (is.null(emb)) {
    stop(sprintf("Embedding '%s' not found", embedding))
  }
  
  # Get grid data
  grid_data <- velocity_grid(object, embedding, n_grid = n_grid)
  
  # Create cell data frame
  df <- data.frame(
    x = emb[, 1],
    y = emb[, 2]
  )
  
  # Add coloring
  if (!is.null(color_by) && color_by %in% colnames(object@meta.data)) {
    df$color <- object@meta.data[[color_by]]
  } else {
    cluster_cols <- c("seurat_clusters", "clusters", "celltype")
    for (col in cluster_cols) {
      if (col %in% colnames(object@meta.data)) {
        df$color <- object@meta.data[[col]]
        break
      }
    }
    if (!"color" %in% names(df)) df$color <- "cells"
  }
  
  # Grid arrow data
  grid_df <- data.frame(
    x = grid_data$x,
    y = grid_data$y,
    vx = grid_data$vx * arrow_scale,
    vy = grid_data$vy * arrow_scale
  )
  
  # Create plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_point(
      data = df,
      ggplot2::aes(x = x, y = y, color = color),
      size = size,
      alpha = alpha
    ) +
    ggplot2::geom_segment(
      data = grid_df,
      ggplot2::aes(x = x, y = y, xend = x + vx, yend = y + vy),
      arrow = ggplot2::arrow(length = ggplot2::unit(0.12, "cm"), type = "closed"),
      linewidth = 0.4,
      color = "black",
      alpha = 0.8
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = paste0(toupper(embedding), "_1"),
      y = paste0(toupper(embedding), "_2"),
      title = title
    ) +
    ggplot2::coord_fixed()
  
  # Color scale
  if (is.numeric(df$color)) {
    p <- p + ggplot2::scale_color_viridis_c()
  } else {
    p <- p + ggplot2::scale_color_manual(values = get_palette(length(unique(df$color))))
  }
  
  p
}

#' Plot Phase Portrait
#'
#' @description Plot phase portrait (spliced vs unspliced) for a gene.
#'
#' @param object Seurat object
#' @param gene Gene name
#' @param color_by Variable to color by
#' @param show_fit Show fitted curve from dynamics
#' @param size Point size
#' @param alpha Transparency
#' @param title Plot title
#'
#' @return ggplot2 object
#' @export
plot_phase <- function(object,
                       gene,
                       color_by = NULL,
                       show_fit = TRUE,
                       size = 1,
                       alpha = 0.7,
                       title = NULL) {
  
  # Get expression data
  data <- get_velocity_data(object)
  Ms <- data$Ms
  Mu <- data$Mu
  
  gene_idx <- which(colnames(Ms) == gene)
  
  if (length(gene_idx) == 0) {
    stop(sprintf("Gene '%s' not found", gene))
  }
  
  s <- Ms[, gene_idx]
  u <- Mu[, gene_idx]
  
  df <- data.frame(s = s, u = u)
  
  # Add coloring
  if (!is.null(color_by) && color_by %in% colnames(object@meta.data)) {
    df$color <- object@meta.data[[color_by]]
  } else if ("latent_time" %in% colnames(object@meta.data)) {
    df$color <- object@meta.data$latent_time
    color_by <- "latent_time"
  } else {
    df$color <- 1
  }
  
  # Create plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = s, y = u, color = color)) +
    ggplot2::geom_point(size = size, alpha = alpha) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "Spliced",
      y = "Unspliced",
      title = if (is.null(title)) gene else title,
      color = color_by
    )
  
  # Add fitted curve if dynamics are available
  if (show_fit && !is.null(object@misc$scVeloR$dynamics)) {
    gene_params <- object@misc$scVeloR$dynamics$gene_params
    
    if (gene %in% rownames(gene_params)) {
      params <- gene_params[gene, ]
      
      if (!is.na(params$alpha)) {
        alpha_p <- params$alpha
        beta_p <- params$beta * params$scaling
        gamma_p <- params$gamma
        t_ <- params$t_
        scaling <- params$scaling
        
        # Generate trajectory
        t_points <- seq(0, t_ * 2, length.out = 200)
        
        traj_u <- numeric(length(t_points))
        traj_s <- numeric(length(t_points))
        
        for (i in seq_along(t_points)) {
          t <- t_points[i]
          if (t < t_) {
            result <- mRNA(t, 0, 0, alpha_p, beta_p, gamma_p)
            traj_u[i] <- result$u * scaling
            traj_s[i] <- result$s
          } else {
            result_switch <- mRNA(t_, 0, 0, alpha_p, beta_p, gamma_p)
            u0_ <- result_switch$u
            s0_ <- result_switch$s
            result <- mRNA(t - t_, u0_, s0_, 0, beta_p, gamma_p)
            traj_u[i] <- result$u * scaling
            traj_s[i] <- result$s
          }
        }
        
        traj_df <- data.frame(s = traj_s, u = traj_u)
        
        p <- p + ggplot2::geom_path(
          data = traj_df,
          ggplot2::aes(x = s, y = u),
          color = "red",
          linewidth = 1,
          inherit.aes = FALSE
        )
      }
    }
  }
  
  # Color scale
  if (is.numeric(df$color)) {
    p <- p + ggplot2::scale_color_viridis_c()
  } else {
    p <- p + ggplot2::scale_color_manual(values = get_palette(length(unique(df$color))))
  }
  
  p
}

#' Get Color Palette
#'
#' @description Get a nice color palette for categorical data.
#'
#' @param n Number of colors needed
#'
#' @return Vector of colors
#' @keywords internal
get_palette <- function(n) {
  if (n <= 10) {
    # Use Set3 for small n
    palette <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3",
                 "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD")
  } else if (n <= 20) {
    # Use combined palette
    palette <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
                 "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
                 "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
                 "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5")
  } else {
    # Generate more colors
    palette <- colorRampPalette(
      c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
        "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
    )(n)
  }
  
  palette[1:n]
}

# scVeloR <img src="man/figures/logo.png" align="right" height="139" />
<!-- badges: start -->
[![R-CMD-check](https://github.com/Zaoqu-Liu/scVeloR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Zaoqu-Liu/scVeloR/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/scVeloR)](https://CRAN.R-project.org/package=scVeloR)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

## RNA Velocity Analysis for Single-Cell RNA-Seq Data in R

**scVeloR** is a comprehensive R implementation of RNA velocity analysis, providing a native R alternative to the Python [scvelo](https://github.com/theislab/scvelo) package. It offers three velocity estimation models (steady-state, stochastic, and dynamical), velocity graph construction, embedding projection, and rich visualization capabilities.

### Key Features

- **Multiple Velocity Models**: Steady-state, stochastic, and dynamical (EM) models
- **Seurat Integration**: Full compatibility with Seurat V4 and V5 objects
- **High Performance**: Rcpp-accelerated core computations
- **Cross-Platform Parallel Computing**: via `future` framework
- **Rich Visualizations**: Embedding plots, streamlines, grid plots, heatmaps

## Installation

```r
# Install from GitHub
if (!require("remotes")) install.packages("remotes")
remotes::install_github("Zaoqu-Liu/scVeloR")
```

### Dependencies

scVeloR requires the following packages:

```r
# Core dependencies (installed automatically)
install.packages(c("Seurat", "Matrix", "ggplot2", "future", "future.apply", 
                   "FNN", "uwot", "irlba", "progressr", "scales"))

# For Rcpp compilation
install.packages(c("Rcpp", "RcppArmadillo"))
```

## Quick Start

```r
library(scVeloR)
library(Seurat)

# Load your Seurat object with spliced/unspliced layers
# The object should have 'spliced' and 'unspliced' assays or layers

# Step 1: Preprocessing
seurat_obj <- filter_and_normalize(seurat_obj)

# Step 2: Compute moments (smoothed expression)
seurat_obj <- compute_moments(seurat_obj, n_neighbors = 30)

# Step 3: Estimate velocities (choose one model)
# Option A: Steady-state model (fastest)
seurat_obj <- velocity(seurat_obj, mode = "steady_state")

# Option B: Stochastic model (recommended)
seurat_obj <- velocity(seurat_obj, mode = "stochastic")

# Option C: Dynamical model (most accurate)
seurat_obj <- recover_dynamics(seurat_obj)
seurat_obj <- velocity(seurat_obj, mode = "dynamical")

# Step 4: Compute velocity graph
seurat_obj <- compute_velocity_graph(seurat_obj)

# Step 5: Project to embedding
seurat_obj <- project_velocity_embedding(seurat_obj, basis = "umap")

# Step 6: Visualize
velocity_embedding_plot(seurat_obj, basis = "umap", color_by = "clusters")
velocity_stream_plot(seurat_obj, basis = "umap")
```

## Velocity Models

### Steady-State Model

The simplest model assuming cells are at steady state equilibrium:

```r
seurat_obj <- velocity(seurat_obj, mode = "steady_state")
```

### Stochastic Model

Incorporates second-order moments for better velocity estimation:
```r
seurat_obj <- velocity(seurat_obj, mode = "stochastic")
```

### Dynamical Model

Full splicing kinetics model with EM optimization:

```r
seurat_obj <- recover_dynamics(seurat_obj, 
                                var_names = "velocity_genes",
                                max_iter = 10)
seurat_obj <- velocity(seurat_obj, mode = "dynamical")
seurat_obj <- recover_latent_time(seurat_obj)
```

## Visualization

### Velocity Embedding Plot

```r
# Basic plot
velocity_embedding_plot(seurat_obj)

# With customization
velocity_embedding_plot(seurat_obj, 
                        basis = "umap",
                        color_by = "clusters",
                        arrow_size = 1,
                        density = 0.5)
```

### Stream Plot

```r
velocity_stream_plot(seurat_obj, 
                     basis = "umap",
                     density = 1,
                     smooth = 0.5)
```

### Grid Plot

```r
velocity_grid_plot(seurat_obj, 
                   basis = "umap",
                   n_grid = 40)
```

## Parallel Computing

scVeloR supports parallel computation via the `future` framework:

```r
library(future)

# Use multiple cores
plan(multisession, workers = 4)

# Run velocity computation
seurat_obj <- velocity(seurat_obj, mode = "stochastic")

# Reset to sequential
plan(sequential)
```

## Data Requirements

Your Seurat object should contain:

1. **Spliced counts**: In `seurat_obj[["spliced"]]` or as a layer
2. **Unspliced counts**: In `seurat_obj[["unspliced"]]` or as a layer
3. **Dimensional reduction**: UMAP or tSNE coordinates in `seurat_obj@reductions`

### Creating from velocyto output

```r
library(Seurat)

# Load loom file
loom_data <- ReadVelocity("your_data.loom")

# Create Seurat object
seurat_obj <- CreateSeuratObject(counts = loom_data$spliced)
seurat_obj[["unspliced"]] <- CreateAssayObject(counts = loom_data$unspliced)
seurat_obj[["spliced"]] <- CreateAssayObject(counts = loom_data$spliced)

# Standard Seurat workflow
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
```

## Citation

If you use scVeloR in your research, please cite:

```bibtex
@software{scVeloR,
  author = {Zaoqu Liu},
  title = {scVeloR: RNA Velocity Analysis for Single-Cell RNA-Seq Data in R},
  year = {2026},
  url = {https://github.com/Zaoqu-Liu/scVeloR}
}
```

Also cite the original scvelo paper:

> Bergen et al. (2020). Generalizing RNA velocity to transient cell states through dynamical modeling. Nature Biotechnology.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Author

**Zaoqu Liu** - [GitHub](https://github.com/Zaoqu-Liu) - liuzaoqu@163.com

# scVeloR: RNA Velocity Analysis in R

<img src="man/figures/logo.png" align="right" height="139" />
  
[![R-CMD-check](https://github.com/Zaoqu-Liu/scVeloR/workflows/R-CMD-check/badge.svg)](https://github.com/Zaoqu-Liu/scVeloR/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

**scVeloR** is a comprehensive R implementation of RNA velocity analysis, based on the [scvelo](https://github.com/theislab/scvelo) Python package. It enables estimation of RNA velocities for single-cell RNA-seq data, supporting:

- **Deterministic (Steady-State) Model**: Assumes cells are near transcriptional equilibrium
- **Stochastic Model**: Accounts for transcriptional bursting using second-order moments
- **Dynamical Model**: Full kinetics inference using EM algorithm

The package integrates seamlessly with **Seurat** objects (V4 and V5) and provides visualization functions for velocity embeddings, streamlines, and phase portraits.

## Installation

```r
# Install from GitHub
devtools::install_github("Zaoqu-Liu/scVeloR")
```

### Dependencies

- R (>= 4.0.0)
- Seurat (>= 4.0.0)
- Matrix
- Rcpp, RcppArmadillo
- ggplot2

## Quick Start

```r
library(scVeloR)
library(Seurat)

# Load your Seurat object with spliced and unspliced layers
# seurat_obj should have 'spliced' and 'unspliced' in assay layers

# Run full velocity pipeline
seurat_obj <- run_velocity(
  seurat_obj, 
  mode = "dynamical",  # or "deterministic", "stochastic"
  embedding = "umap"
)

# Visualize results
plot_velocity(seurat_obj, embedding = "umap")
plot_velocity_stream(seurat_obj, embedding = "umap")
```

## Detailed Usage

### Step-by-Step Analysis

```r
# 1. Prepare data
seurat_obj <- prepare_velocity(seurat_obj, n_neighbors = 30)

# 2. Compute velocity
seurat_obj <- velocity(seurat_obj, mode = "dynamical")

# 3. Build velocity graph
seurat_obj <- velocity_graph(seurat_obj)

# 4. Project onto embedding
seurat_obj <- velocity_embedding(seurat_obj, embedding_name = "umap")

# 5. Compute latent time (for dynamical model)
seurat_obj <- compute_latent_time(seurat_obj)
```

### Velocity Models

#### Deterministic (Steady-State)

```r
seurat_obj <- velocity(seurat_obj, mode = "deterministic", min_r2 = 0.01)
```

The steady-state model assumes `du/dt ≈ 0` at equilibrium, giving `γ = u/s`.

#### Stochastic

```r
seurat_obj <- velocity(seurat_obj, mode = "stochastic")
```

Uses second-order moments to account for transcriptional noise.

#### Dynamical

```r
seurat_obj <- velocity(seurat_obj, mode = "dynamical", max_iter = 10)
```

Infers full transcriptional kinetics (α, β, γ) using EM algorithm.

### Visualization

```r
# Velocity arrows
plot_velocity(seurat_obj, color_by = "seurat_clusters", arrow_scale = 1)

# Velocity streamlines
plot_velocity_stream(seurat_obj, density = 1)

# Velocity grid
plot_velocity_grid(seurat_obj, n_grid = 40)

# Phase portrait for specific gene
plot_phase(seurat_obj, gene = "Gene1", show_fit = TRUE)
```

### Gene Ranking

```r
# Rank genes by velocity fit quality
top_genes <- rank_velocity_genes(seurat_obj, n_top = 100)
head(top_genes)
```

## Data Requirements

Your Seurat object should contain:

1. **Spliced counts**: In a layer named "spliced"
2. **Unspliced counts**: In a layer named "unspliced"
3. **Dimensionality reduction**: UMAP or tSNE for visualization

### Creating from loom file

```r
# If you have a loom file from velocyto
library(SeuratDisk)
seurat_obj <- LoadLoom("your_file.loom")
```

### Creating from AnnData

```r
# Import from AnnData h5ad file
seurat_obj <- import_from_anndata("your_file.h5ad")
```

## Output

scVeloR stores results in `object@misc$scVeloR`:

- `velocity`: Velocity matrix and parameters
- `dynamics`: Dynamical model parameters (if used)
- `velocity_graph`: Velocity graph and transition matrix
- `velocity_embedding`: Projected velocity vectors

Cell metadata (`object@meta.data`) includes:

- `latent_time`: Gene-shared latent time
- `velocity_pseudotime`: Diffusion-based pseudotime
- `velocity_confidence`: Confidence scores

## Mathematical Background

### RNA Velocity Model

The model describes mRNA dynamics:

```
du/dt = α - β·u  (unspliced)
ds/dt = β·u - γ·s  (spliced)
```

Where:
- α: Transcription rate
- β: Splicing rate  
- γ: Degradation rate

### Analytical Solutions

**Unspliced:**
```
u(t) = u₀·e^(-βt) + α/β·(1 - e^(-βt))
```

**Spliced:**
```
s(t) = s₀·e^(-γt) + α/γ·(1 - e^(-γt)) + c·(e^(-γt) - e^(-βt))
```

where `c = (α - u₀·β)/(γ - β)`

## Citation

If you use scVeloR in your research, please cite:

```
@software{scVeloR,
  author = {Zaoqu Liu},
  title = {scVeloR: RNA Velocity Analysis in R},
  url = {https://github.com/Zaoqu-Liu/scVeloR},
  year = {2024}
}
```

And the original scvelo paper:

```
@article{bergen2020generalizing,
  title={Generalizing RNA velocity to transient cell states through dynamical modeling},
  author={Bergen, Volker and Lange, Marius and Peidli, Stefan and Wolf, F Alexander and Theis, Fabian J},
  journal={Nature biotechnology},
  year={2020}
}
```

## Contributing

Contributions are welcome! Please feel free to submit issues or pull requests.

## License

MIT License - see [LICENSE](LICENSE) for details.

## Acknowledgments

This package is inspired by and based on the [scvelo](https://github.com/theislab/scvelo) Python package developed by the Theis Lab.

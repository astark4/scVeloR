# Helper functions for testing

# Create a minimal mock Seurat object for testing
# This avoids needing actual Seurat data
create_mock_data <- function(n_cells = 50, n_genes = 100) {
  # Spliced counts
  spliced <- Matrix::rsparsematrix(n_genes, n_cells, density = 0.3)
  spliced@x <- abs(spliced@x) * 10
  rownames(spliced) <- paste0("gene", seq_len(n_genes))
  colnames(spliced) <- paste0("cell", seq_len(n_cells))
  
  # Unspliced counts (correlated with spliced)
  unspliced <- spliced * 0.3 + Matrix::rsparsematrix(n_genes, n_cells, density = 0.2)
  unspliced@x <- abs(unspliced@x) * 5
  rownames(unspliced) <- rownames(spliced)
  colnames(unspliced) <- colnames(spliced)
  
  # Simple PCA embedding
  pca <- matrix(rnorm(n_cells * 30), ncol = 30)
  rownames(pca) <- colnames(spliced)
  colnames(pca) <- paste0("PC_", seq_len(30))
  
  # UMAP embedding
  umap <- matrix(rnorm(n_cells * 2, sd = 5), ncol = 2)
  rownames(umap) <- colnames(spliced)
  colnames(umap) <- c("UMAP_1", "UMAP_2")
  
  list(
    spliced = spliced,
    unspliced = unspliced,
    pca = pca,
    umap = umap,
    n_cells = n_cells,
    n_genes = n_genes
  )
}

# Skip tests if Seurat not available
skip_if_no_seurat <- function() {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    skip("Seurat package not available")
  }
}

# Create test velocity data
create_test_velocity <- function(mock_data) {
  n_genes <- mock_data$n_genes
  n_cells <- mock_data$n_cells
  
  gamma <- runif(n_genes, 0.1, 1)
  velocity <- mock_data$unspliced - sweep(mock_data$spliced, 1, gamma, "*")
  
  list(
    velocity = velocity,
    gamma = gamma
  )
}

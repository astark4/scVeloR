// cosine_similarity.cpp - Fast cosine similarity computation
// Author: Zaoqu Liu
// 
// Computes cosine similarity between velocity vectors and cell-cell differences
// This is a critical hot-spot function in velocity graph computation

#include "scVeloR_types.h"

using namespace Rcpp;

//' Compute L2 norm for each row of a matrix
//'
//' @param X Matrix (n x p)
//' @return Vector of L2 norms (length n)
//' @keywords internal
// [[Rcpp::export]]
NumericVector l2_norm_rows_cpp(const NumericMatrix& X) {
    int n = X.nrow();
    int p = X.ncol();
    NumericVector norms(n);
    
    for (int i = 0; i < n; i++) {
        double sum_sq = 0.0;
        for (int j = 0; j < p; j++) {
            double val = X(i, j);
            if (!NumericVector::is_na(val)) {
                sum_sq += val * val;
            }
        }
        norms[i] = std::sqrt(sum_sq);
    }
    
    return norms;
}

//' Compute cosine correlation between a vector and rows of a matrix
//'
//' @param dX Matrix of differences (n_neighbors x n_genes)
//' @param v Velocity vector (length n_genes)
//' @return Vector of cosine correlations (length n_neighbors)
//' @keywords internal
// [[Rcpp::export]]
NumericVector cosine_correlation_cpp(const NumericMatrix& dX, 
                                      const NumericVector& v) {
    int n = dX.nrow();
    int p = dX.ncol();
    NumericVector result(n);
    
    // Compute norm of v
    double v_norm = 0.0;
    for (int j = 0; j < p; j++) {
        if (!NumericVector::is_na(v[j])) {
            v_norm += v[j] * v[j];
        }
    }
    v_norm = std::sqrt(v_norm);
    
    if (v_norm < 1e-10) {
        // Zero velocity vector
        std::fill(result.begin(), result.end(), 0.0);
        return result;
    }
    
    // Compute cosine correlation for each row
    for (int i = 0; i < n; i++) {
        double dot_prod = 0.0;
        double dx_norm = 0.0;
        
        for (int j = 0; j < p; j++) {
            double dx_val = dX(i, j);
            double v_val = v[j];
            
            if (!NumericVector::is_na(dx_val) && !NumericVector::is_na(v_val)) {
                dot_prod += dx_val * v_val;
                dx_norm += dx_val * dx_val;
            }
        }
        
        dx_norm = std::sqrt(dx_norm);
        
        if (dx_norm < 1e-10) {
            result[i] = 0.0;
        } else {
            result[i] = dot_prod / (dx_norm * v_norm);
        }
    }
    
    return result;
}

//' Compute velocity graph for a single cell
//'
//' @param X Expression matrix (cells x genes)
//' @param V Velocity matrix (cells x genes)
//' @param cell_idx Index of the cell (0-based)
//' @param neighbor_indices Indices of neighbors (0-based)
//' @param sqrt_transform Whether to apply sqrt transform
//' @return Vector of cosine correlations with neighbors
//' @keywords internal
// [[Rcpp::export]]
NumericVector compute_cell_velocity_graph_cpp(
    const NumericMatrix& X,
    const NumericMatrix& V,
    int cell_idx,
    const IntegerVector& neighbor_indices,
    bool sqrt_transform = false) {
    
    int n_neighbors = neighbor_indices.size();
    int n_genes = X.ncol();
    
    // Get velocity vector for this cell
    NumericVector v(n_genes);
    for (int j = 0; j < n_genes; j++) {
        v[j] = V(cell_idx, j);
    }
    
    // Check if velocity is all zeros/NaN
    double v_sum = 0.0;
    for (int j = 0; j < n_genes; j++) {
        if (!NumericVector::is_na(v[j])) {
            v_sum += std::abs(v[j]);
        }
    }
    
    NumericVector result(n_neighbors);
    
    if (v_sum < 1e-10) {
        // Zero velocity - return zeros
        return result;
    }
    
    // Compute differences matrix dX = X[neighbors] - X[cell]
    NumericMatrix dX(n_neighbors, n_genes);
    
    for (int i = 0; i < n_neighbors; i++) {
        int neighbor = neighbor_indices[i];
        for (int j = 0; j < n_genes; j++) {
            double diff = X(neighbor, j) - X(cell_idx, j);
            if (sqrt_transform) {
                // Apply sqrt transform with sign preservation
                diff = std::copysign(std::sqrt(std::abs(diff)), diff);
            }
            dX(i, j) = diff;
        }
    }
    
    // Apply sqrt transform to velocity if needed
    NumericVector v_transformed(n_genes);
    if (sqrt_transform) {
        for (int j = 0; j < n_genes; j++) {
            v_transformed[j] = std::copysign(std::sqrt(std::abs(v[j])), v[j]);
        }
    } else {
        v_transformed = v;
    }
    
    // Center velocity by subtracting mean
    double v_mean = 0.0;
    int v_count = 0;
    for (int j = 0; j < n_genes; j++) {
        if (!NumericVector::is_na(v_transformed[j])) {
            v_mean += v_transformed[j];
            v_count++;
        }
    }
    if (v_count > 0) {
        v_mean /= v_count;
    }
    for (int j = 0; j < n_genes; j++) {
        if (!NumericVector::is_na(v_transformed[j])) {
            v_transformed[j] -= v_mean;
        }
    }
    
    // Compute cosine correlation
    result = cosine_correlation_cpp(dX, v_transformed);
    
    return result;
}

//' Batch compute velocity graph for multiple cells
//'
//' @param X Expression matrix (cells x genes)
//' @param V Velocity matrix (cells x genes)
//' @param indices_list List of neighbor indices for each cell
//' @param sqrt_transform Whether to apply sqrt transform
//' @return List containing:
//'   - vals: correlation values
//'   - rows: row indices
//'   - cols: column indices
//' @keywords internal
// [[Rcpp::export]]
List compute_velocity_graph_batch_cpp(
    const NumericMatrix& X,
    const NumericMatrix& V,
    const List& indices_list,
    bool sqrt_transform = false) {
    
    int n_cells = indices_list.size();
    
    // Estimate total size
    int total_size = 0;
    for (int i = 0; i < n_cells; i++) {
        IntegerVector idx = indices_list[i];
        total_size += idx.size();
    }
    
    // Pre-allocate vectors
    std::vector<double> vals;
    std::vector<int> rows;
    std::vector<int> cols;
    vals.reserve(total_size);
    rows.reserve(total_size);
    cols.reserve(total_size);
    
    // Process each cell
    for (int cell_i = 0; cell_i < n_cells; cell_i++) {
        IntegerVector neighbor_indices = indices_list[cell_i];
        
        if (neighbor_indices.size() == 0) {
            continue;
        }
        
        NumericVector correlations = compute_cell_velocity_graph_cpp(
            X, V, cell_i, neighbor_indices, sqrt_transform
        );
        
        // Store results
        for (int j = 0; j < correlations.size(); j++) {
            if (!NumericVector::is_na(correlations[j])) {
                vals.push_back(correlations[j]);
                rows.push_back(cell_i);
                cols.push_back(neighbor_indices[j]);
            }
        }
    }
    
    return List::create(
        Named("vals") = wrap(vals),
        Named("rows") = wrap(rows),
        Named("cols") = wrap(cols)
    );
}

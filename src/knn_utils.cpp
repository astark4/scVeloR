// knn_utils.cpp - KNN-related utility functions
// Author: Zaoqu Liu
//
// Helper functions for neighbor graph construction

#include "scVeloR_types.h"

using namespace Rcpp;

//' Get iterative neighbor indices (neighbors of neighbors)
//'
//' @param neighbor_matrix Matrix of neighbor indices (cells x n_neighbors), 0-based
//' @param cell_idx Cell index (0-based)
//' @param n_recurse Number of recursions
//' @param max_neighbors Maximum neighbors to return
//' @return Vector of unique neighbor indices
//' @keywords internal
// [[Rcpp::export]]
IntegerVector get_iterative_neighbors_cpp(const IntegerMatrix& neighbor_matrix,
                                           int cell_idx,
                                           int n_recurse = 2,
                                           int max_neighbors = 100) {
    int n_neighbors = neighbor_matrix.ncol();
    
    // Use set to track unique neighbors
    std::set<int> neighbor_set;
    std::vector<int> current_level;
    
    // Start with direct neighbors
    for (int j = 0; j < n_neighbors; j++) {
        int neighbor = neighbor_matrix(cell_idx, j);
        if (neighbor >= 0 && neighbor != cell_idx) {
            neighbor_set.insert(neighbor);
            current_level.push_back(neighbor);
        }
    }
    
    // Recursively add neighbors of neighbors
    for (int r = 1; r < n_recurse; r++) {
        std::vector<int> next_level;
        
        for (int idx : current_level) {
            for (int j = 0; j < n_neighbors; j++) {
                int neighbor = neighbor_matrix(idx, j);
                if (neighbor >= 0 && neighbor != cell_idx) {
                    if (neighbor_set.find(neighbor) == neighbor_set.end()) {
                        neighbor_set.insert(neighbor);
                        next_level.push_back(neighbor);
                    }
                }
            }
        }
        
        current_level = next_level;
        
        // Stop if we have enough neighbors
        if ((int)neighbor_set.size() >= max_neighbors) {
            break;
        }
    }
    
    // Convert to vector
    IntegerVector result(neighbor_set.begin(), neighbor_set.end());
    
    // Randomly sample if too many
    if (result.size() > max_neighbors) {
        // Simple random sampling without replacement
        IntegerVector indices = seq(0, result.size() - 1);
        std::random_shuffle(indices.begin(), indices.end());
        
        IntegerVector sampled(max_neighbors);
        for (int i = 0; i < max_neighbors; i++) {
            sampled[i] = result[indices[i]];
        }
        return sampled;
    }
    
    return result;
}

//' Convert neighbor indices to sparse connectivity matrix
//'
//' @param indices Matrix of neighbor indices (cells x n_neighbors), 0-based
//' @param distances Matrix of distances (cells x n_neighbors)
//' @param n_cells Total number of cells
//' @param mode "distances" or "connectivities"
//' @return Sparse matrix
//' @keywords internal
// [[Rcpp::export]]
arma::sp_mat indices_to_sparse_cpp(const IntegerMatrix& indices,
                                    const NumericMatrix& distances,
                                    int n_cells,
                                    std::string mode = "connectivities") {
    int n_neighbors = indices.ncol();
    
    // Count non-zero entries
    std::vector<arma::uword> row_idx;
    std::vector<arma::uword> col_idx;
    std::vector<double> values;
    
    row_idx.reserve(n_cells * n_neighbors);
    col_idx.reserve(n_cells * n_neighbors);
    values.reserve(n_cells * n_neighbors);
    
    for (int i = 0; i < n_cells; i++) {
        for (int j = 0; j < n_neighbors; j++) {
            int neighbor = indices(i, j);
            if (neighbor >= 0 && neighbor < n_cells) {
                row_idx.push_back(i);
                col_idx.push_back(neighbor);
                
                if (mode == "distances") {
                    values.push_back(distances(i, j));
                } else {
                    // Connectivity: 1 if neighbor, 0 otherwise
                    values.push_back(1.0);
                }
            }
        }
    }
    
    // Create sparse matrix
    arma::umat locations(2, row_idx.size());
    for (size_t i = 0; i < row_idx.size(); i++) {
        locations(0, i) = row_idx[i];
        locations(1, i) = col_idx[i];
    }
    
    arma::vec vals(values);
    
    return arma::sp_mat(locations, vals, n_cells, n_cells);
}

//' Compute adaptive kernel width for connectivity
//'
//' Based on local bandwidth estimation as in UMAP
//'
//' @param distances Matrix of distances (cells x n_neighbors)
//' @param local_connectivity Local connectivity parameter (default: 1)
//' @return Vector of kernel widths (sigma) per cell
//' @keywords internal
// [[Rcpp::export]]
NumericVector compute_sigma_cpp(const NumericMatrix& distances,
                                 double local_connectivity = 1.0) {
    int n_cells = distances.nrow();
    int n_neighbors = distances.ncol();
    
    NumericVector sigma(n_cells);
    
    // Target is to have local_connectivity neighbors "effectively connected"
    double target = std::log2(n_neighbors);
    
    for (int i = 0; i < n_cells; i++) {
        // Get sorted distances for this cell
        std::vector<double> dists(n_neighbors);
        for (int j = 0; j < n_neighbors; j++) {
            dists[j] = distances(i, j);
        }
        std::sort(dists.begin(), dists.end());
        
        // Estimate rho (distance to closest neighbor)
        double rho = dists[0];
        if (rho < 1e-10) {
            rho = dists[1];  // Skip if first is essentially zero
        }
        
        // Binary search for sigma
        double lo = 0.0;
        double hi = 1000.0;  // Upper bound
        double mid = 1.0;
        
        for (int iter = 0; iter < 64; iter++) {
            mid = (lo + hi) / 2.0;
            
            // Compute effective number of neighbors
            double sum = 0.0;
            for (int j = 0; j < n_neighbors; j++) {
                double d = std::max(0.0, dists[j] - rho);
                sum += std::exp(-d / mid);
            }
            
            if (std::abs(sum - target) < 1e-5) {
                break;
            }
            
            if (sum < target) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        
        sigma[i] = mid;
    }
    
    return sigma;
}

//' Compute UMAP-style connectivities from distances
//'
//' @param distances Sparse distance matrix
//' @param set_op_mix_ratio Mix ratio for set operations (default: 1.0)
//' @return Sparse connectivity matrix
//' @keywords internal
// [[Rcpp::export]]
arma::sp_mat compute_connectivities_umap_cpp(const arma::sp_mat& distances,
                                              double set_op_mix_ratio = 1.0) {
    int n = distances.n_rows;
    
    // Find local sigmas
    std::vector<double> sigmas(n, 1.0);
    std::vector<double> rhos(n, 0.0);
    
    // Iterate over rows
    for (int i = 0; i < n; i++) {
        // Get non-zero distances for this row
        std::vector<double> row_dists;
        for (arma::sp_mat::const_row_iterator it = distances.begin_row(i);
             it != distances.end_row(i); ++it) {
            if (*it > 0) {
                row_dists.push_back(*it);
            }
        }
        
        if (row_dists.size() > 0) {
            std::sort(row_dists.begin(), row_dists.end());
            rhos[i] = row_dists[0];
            
            // Estimate sigma using binary search
            double target = std::log2(row_dists.size());
            double lo = 0.0, hi = 1000.0, mid = 1.0;
            
            for (int iter = 0; iter < 64; iter++) {
                mid = (lo + hi) / 2.0;
                double sum = 0.0;
                for (double d : row_dists) {
                    sum += std::exp(-std::max(0.0, d - rhos[i]) / mid);
                }
                if (std::abs(sum - target) < 1e-5) break;
                if (sum < target) lo = mid;
                else hi = mid;
            }
            sigmas[i] = mid;
        }
    }
    
    // Compute connectivities
    arma::sp_mat conn = distances;
    
    for (arma::sp_mat::iterator it = conn.begin(); it != conn.end(); ++it) {
        int i = it.row();
        double d = *it;
        double rho_i = rhos[i];
        double sigma_i = sigmas[i];
        
        double new_val = std::exp(-std::max(0.0, d - rho_i) / sigma_i);
        *it = new_val;
    }
    
    // Symmetrize: fuzzy set union
    arma::sp_mat conn_t = conn.t();
    arma::sp_mat result = conn + conn_t - conn % conn_t;
    
    // Mix with simple average if needed
    if (set_op_mix_ratio < 1.0) {
        arma::sp_mat simple_avg = (conn + conn_t) / 2.0;
        result = set_op_mix_ratio * result + (1.0 - set_op_mix_ratio) * simple_avg;
    }
    
    return result;
}

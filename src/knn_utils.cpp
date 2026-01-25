// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Compute Euclidean Distance Matrix
//' 
//' @param X Data matrix (cells x features)
//' @return Distance matrix
//' @keywords internal
// [[Rcpp::export]]
arma::mat euclidean_distances_cpp(const arma::mat& X) {
    int n = X.n_rows;
    arma::mat D(n, n, arma::fill::zeros);
    
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            double dist = arma::norm(X.row(i) - X.row(j));
            D(i, j) = dist;
            D(j, i) = dist;
        }
    }
    
    return D;
}

//' Find K Nearest Neighbors
//' 
//' @param D Distance matrix
//' @param k Number of neighbors
//' @return List with indices and distances
//' @keywords internal
// [[Rcpp::export]]
List knn_from_distances_cpp(const arma::mat& D, int k) {
    int n = D.n_rows;
    
    arma::imat indices(n, k);
    arma::mat distances(n, k);
    
    for (int i = 0; i < n; i++) {
        arma::vec d = D.row(i).t();
        d(i) = arma::datum::inf;  // Exclude self
        
        arma::uvec sorted_idx = arma::sort_index(d);
        
        for (int j = 0; j < k; j++) {
            indices(i, j) = sorted_idx(j) + 1;  // 1-based indexing
            distances(i, j) = d(sorted_idx(j));
        }
    }
    
    return List::create(
        Named("indices") = indices,
        Named("distances") = distances
    );
}

//' Compute Cosine Distance Matrix
//' 
//' @param X Data matrix (cells x features)
//' @return Cosine distance matrix (1 - cosine similarity)
//' @keywords internal
// [[Rcpp::export]]
arma::mat cosine_distances_cpp(const arma::mat& X) {
    int n = X.n_rows;
    
    // Normalize rows
    arma::mat X_norm = X;
    for (int i = 0; i < n; i++) {
        double norm = arma::norm(X.row(i));
        if (norm > 0) {
            X_norm.row(i) /= norm;
        }
    }
    
    // Cosine similarity
    arma::mat sim = X_norm * X_norm.t();
    
    // Convert to distance
    arma::mat D = 1.0 - sim;
    
    return D;
}

//' Compute Connectivity Weights
//' 
//' @description UMAP-style connectivity computation
//' 
//' @param indices Neighbor indices (1-based)
//' @param distances Neighbor distances
//' @param n_cells Number of cells
//' @param local_connectivity Local connectivity parameter
//' 
//' @return List with i, j, x for sparse matrix construction
//' @keywords internal
// [[Rcpp::export]]
List compute_connectivities_cpp(const arma::imat& indices,
                                 const arma::mat& distances,
                                 int n_cells,
                                 double local_connectivity = 1.0) {
    
    int n_neighbors = indices.n_cols;
    
    std::vector<int> i_vec;
    std::vector<int> j_vec;
    std::vector<double> x_vec;
    
    i_vec.reserve(n_cells * n_neighbors * 2);
    j_vec.reserve(n_cells * n_neighbors * 2);
    x_vec.reserve(n_cells * n_neighbors * 2);
    
    // For each cell, compute bandwidth and weights
    for (int i = 0; i < n_cells; i++) {
        arma::rowvec d = distances.row(i);
        
        // Find rho (distance to nearest neighbor based on local_connectivity)
        arma::vec d_sorted = arma::sort(d.t());
        int rho_idx = std::min((int)local_connectivity - 1, (int)d_sorted.n_elem - 1);
        rho_idx = std::max(rho_idx, 0);
        double rho = d_sorted(rho_idx);
        
        // Binary search for sigma
        double target = std::log2(n_neighbors);
        double lo = 1e-10;
        double hi = arma::max(d) * 10;
        double sigma = 1.0;
        
        for (int iter = 0; iter < 64; iter++) {
            double mid = (lo + hi) / 2.0;
            
            double val = 0.0;
            for (int k = 0; k < n_neighbors; k++) {
                double dk = std::max(d(k) - rho, 0.0);
                val += std::exp(-dk / mid);
            }
            
            if (std::abs(val - target) < 1e-5) {
                sigma = mid;
                break;
            }
            
            if (val > target) {
                hi = mid;
            } else {
                lo = mid;
            }
            sigma = mid;
        }
        
        // Compute weights
        for (int k = 0; k < n_neighbors; k++) {
            int j = indices(i, k) - 1;  // Convert to 0-based
            
            if (j >= 0 && j < n_cells && j != i) {
                double dk = std::max(d(k) - rho, 0.0);
                double w = std::exp(-dk / sigma);
                
                if (w > 0) {
                    i_vec.push_back(i);
                    j_vec.push_back(j);
                    x_vec.push_back(w);
                }
            }
        }
    }
    
    return List::create(
        Named("i") = wrap(i_vec),
        Named("j") = wrap(j_vec),
        Named("x") = wrap(x_vec)
    );
}

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
using namespace Rcpp;

//' Compute Cosine Similarity for Velocity Graph
//' 
//' @description Fast C++ implementation of cosine similarity computation
//' between velocity vectors and cell-cell displacement vectors.
//' 
//' @param velocity Velocity matrix (cells x genes)
//' @param expression Expression matrix (cells x genes)
//' @param indices Neighbor indices matrix (cells x n_neighbors, 1-based)
//' 
//' @return List with i, j, x vectors for sparse matrix construction
//' @keywords internal
// [[Rcpp::export]]
List compute_cosine_similarity_cpp(const arma::mat& velocity,
                                    const arma::mat& expression,
                                    const arma::imat& indices) {
    
    int n_cells = velocity.n_rows;
    int n_neighbors = indices.n_cols;
    
    // Pre-allocate vectors
    std::vector<int> i_vec;
    std::vector<int> j_vec;
    std::vector<double> x_vec;
    
    // Reserve space (estimate)
    i_vec.reserve(n_cells * n_neighbors);
    j_vec.reserve(n_cells * n_neighbors);
    x_vec.reserve(n_cells * n_neighbors);
    
    for (int i = 0; i < n_cells; i++) {
        // Get velocity vector for cell i
        arma::rowvec v = velocity.row(i);
        double v_norm = arma::norm(v);
        
        if (v_norm == 0 || !std::isfinite(v_norm)) {
            continue;
        }
        
        // Iterate over neighbors
        for (int k = 0; k < n_neighbors; k++) {
            int j = indices(i, k) - 1;  // Convert to 0-based indexing
            
            if (j < 0 || j >= n_cells || j == i) {
                continue;
            }
            
            // Compute displacement vector
            arma::rowvec dx = expression.row(j) - expression.row(i);
            double dx_norm = arma::norm(dx);
            
            if (dx_norm == 0 || !std::isfinite(dx_norm)) {
                continue;
            }
            
            // Cosine similarity
            double cos_sim = arma::dot(v, dx) / (v_norm * dx_norm);
            
            // Only store positive correlations
            if (cos_sim > 0 && std::isfinite(cos_sim)) {
                i_vec.push_back(i);  // 0-based
                j_vec.push_back(j);  // 0-based
                x_vec.push_back(cos_sim);
            }
        }
    }
    
    return List::create(
        Named("i") = wrap(i_vec),
        Named("j") = wrap(j_vec),
        Named("x") = wrap(x_vec)
    );
}

//' Compute Correlation-Based Velocity Graph
//' 
//' @description Alternative implementation using Pearson correlation.
//' 
//' @param velocity Velocity matrix (cells x genes)
//' @param expression Expression matrix (cells x genes)
//' @param indices Neighbor indices matrix (cells x n_neighbors, 1-based)
//' 
//' @return List with i, j, x vectors for sparse matrix construction
//' @keywords internal
// [[Rcpp::export]]
List compute_correlation_graph_cpp(const arma::mat& velocity,
                                    const arma::mat& expression,
                                    const arma::imat& indices) {
    
    int n_cells = velocity.n_rows;
    int n_genes = velocity.n_cols;
    int n_neighbors = indices.n_cols;
    
    std::vector<int> i_vec;
    std::vector<int> j_vec;
    std::vector<double> x_vec;
    
    i_vec.reserve(n_cells * n_neighbors);
    j_vec.reserve(n_cells * n_neighbors);
    x_vec.reserve(n_cells * n_neighbors);
    
    for (int i = 0; i < n_cells; i++) {
        arma::rowvec v = velocity.row(i);
        double v_mean = arma::mean(v);
        arma::rowvec v_centered = v - v_mean;
        double v_std = arma::stddev(v);
        
        if (v_std == 0 || !std::isfinite(v_std)) {
            continue;
        }
        
        for (int k = 0; k < n_neighbors; k++) {
            int j = indices(i, k) - 1;
            
            if (j < 0 || j >= n_cells || j == i) {
                continue;
            }
            
            // Displacement vector
            arma::rowvec dx = expression.row(j) - expression.row(i);
            double dx_mean = arma::mean(dx);
            arma::rowvec dx_centered = dx - dx_mean;
            double dx_std = arma::stddev(dx);
            
            if (dx_std == 0 || !std::isfinite(dx_std)) {
                continue;
            }
            
            // Pearson correlation
            double corr = arma::dot(v_centered, dx_centered) / 
                          (n_genes * v_std * dx_std);
            
            if (corr > 0 && std::isfinite(corr)) {
                i_vec.push_back(i);
                j_vec.push_back(j);
                x_vec.push_back(corr);
            }
        }
    }
    
    return List::create(
        Named("i") = wrap(i_vec),
        Named("j") = wrap(j_vec),
        Named("x") = wrap(x_vec)
    );
}

//' Row-Normalize Sparse Matrix
//' 
//' @description Normalize rows of sparse matrix to sum to 1.
//' 
//' @param x Values of sparse matrix
//' @param i Row indices (0-based)
//' @param n_rows Number of rows
//' 
//' @return Normalized values
//' @keywords internal
// [[Rcpp::export]]
NumericVector row_normalize_sparse_cpp(const NumericVector& x,
                                        const IntegerVector& i,
                                        int n_rows) {
    
    int nnz = x.size();
    
    // Compute row sums
    std::vector<double> row_sums(n_rows, 0.0);
    for (int k = 0; k < nnz; k++) {
        row_sums[i[k]] += x[k];
    }
    
    // Normalize
    NumericVector x_norm(nnz);
    for (int k = 0; k < nnz; k++) {
        if (row_sums[i[k]] > 0) {
            x_norm[k] = x[k] / row_sums[i[k]];
        } else {
            x_norm[k] = 0;
        }
    }
    
    return x_norm;
}

// moments.cpp - Fast moment computation
// Author: Zaoqu Liu
//
// Computes first and second order moments for velocity estimation

#include "scVeloR_types.h"

using namespace Rcpp;

//' Compute first-order moments (weighted average) using sparse connectivity matrix
//'
//' @param X Sparse expression matrix (genes x cells)
//' @param conn Sparse connectivity matrix (cells x cells), row-normalized
//' @return Dense matrix of first-order moments (genes x cells)
//' @keywords internal
// [[Rcpp::export]]
arma::mat compute_moments_sparse_cpp(const arma::sp_mat& X, 
                                      const arma::sp_mat& conn) {
    // Result = X * conn^T (genes x cells)
    // conn is cells x cells, we want weighted average across neighbors
    
    // conn should be row-normalized (rows sum to 1)
    arma::mat result = X * conn.t();
    
    return result;
}

//' Compute first-order moments using dense matrices
//'
//' @param X Dense expression matrix (genes x cells)
//' @param conn Dense connectivity matrix (cells x cells), row-normalized
//' @return Dense matrix of first-order moments (genes x cells)
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix compute_moments_dense_cpp(const NumericMatrix& X, 
                                         const NumericMatrix& conn) {
    int n_genes = X.nrow();
    int n_cells = X.ncol();
    
    NumericMatrix result(n_genes, n_cells);
    
    // For each gene
    for (int g = 0; g < n_genes; g++) {
        // For each cell
        for (int i = 0; i < n_cells; i++) {
            double moment = 0.0;
            // Weighted sum across all cells
            for (int j = 0; j < n_cells; j++) {
                double weight = conn(i, j);
                if (weight > 0) {
                    moment += weight * X(g, j);
                }
            }
            result(g, i) = moment;
        }
    }
    
    return result;
}

//' Compute second-order moments (X * X)
//'
//' @param X Sparse expression matrix (genes x cells)
//' @param conn Sparse connectivity matrix (cells x cells), row-normalized
//' @return Dense matrix of second-order moments (genes x cells)
//' @keywords internal
// [[Rcpp::export]]
arma::mat compute_second_moments_sparse_cpp(const arma::sp_mat& X, 
                                             const arma::sp_mat& conn) {
    // X_sq = X .* X (element-wise square)
    arma::sp_mat X_sq = X % X;
    
    // Result = X_sq * conn^T
    arma::mat result = X_sq * conn.t();
    
    return result;
}

//' Compute cross-moments (U * S)
//'
//' @param U Sparse unspliced matrix (genes x cells)
//' @param S Sparse spliced matrix (genes x cells)
//' @param conn Sparse connectivity matrix (cells x cells), row-normalized
//' @return Dense matrix of cross-moments (genes x cells)
//' @keywords internal
// [[Rcpp::export]]
arma::mat compute_cross_moments_sparse_cpp(const arma::sp_mat& U, 
                                            const arma::sp_mat& S,
                                            const arma::sp_mat& conn) {
    // US = U .* S (element-wise product)
    arma::sp_mat US = U % S;
    
    // Result = US * conn^T
    arma::mat result = US * conn.t();
    
    return result;
}

//' Compute moments for a subset of cells using indices
//'
//' @param X Expression matrix (genes x cells)
//' @param neighbor_indices Matrix of neighbor indices (cells x n_neighbors), 0-based
//' @param weights Matrix of neighbor weights (cells x n_neighbors)
//' @return Matrix of moments (genes x cells)
//' @keywords internal
// [[Rcpp::export]]
NumericMatrix compute_moments_indexed_cpp(const NumericMatrix& X,
                                           const IntegerMatrix& neighbor_indices,
                                           const NumericMatrix& weights) {
    int n_genes = X.nrow();
    int n_cells = X.ncol();
    int n_neighbors = neighbor_indices.ncol();
    
    NumericMatrix result(n_genes, n_cells);
    
    // For each cell
    for (int i = 0; i < n_cells; i++) {
        // For each gene
        for (int g = 0; g < n_genes; g++) {
            double moment = 0.0;
            double weight_sum = 0.0;
            
            // Average over neighbors
            for (int k = 0; k < n_neighbors; k++) {
                int neighbor = neighbor_indices(i, k);
                if (neighbor >= 0 && neighbor < n_cells) {
                    double w = weights(i, k);
                    moment += w * X(g, neighbor);
                    weight_sum += w;
                }
            }
            
            // Normalize
            if (weight_sum > 0) {
                result(g, i) = moment / weight_sum;
            } else {
                result(g, i) = X(g, i);
            }
        }
    }
    
    return result;
}

//' Compute connectivity-weighted row sums for normalization
//'
//' @param conn Sparse connectivity matrix (cells x cells)
//' @return Vector of row sums
//' @keywords internal
// [[Rcpp::export]]
arma::vec connectivity_row_sums_cpp(const arma::sp_mat& conn) {
    return arma::vec(arma::sum(conn, 1));
}

//' Normalize connectivity matrix by row sums
//'
//' @param conn Sparse connectivity matrix (cells x cells)
//' @return Row-normalized sparse connectivity matrix
//' @keywords internal
// [[Rcpp::export]]
arma::sp_mat normalize_connectivity_cpp(const arma::sp_mat& conn) {
    arma::vec row_sums = arma::vec(arma::sum(conn, 1));
    
    // Avoid division by zero
    row_sums.transform([](double val) { return (val > 0) ? val : 1.0; });
    
    // Create diagonal matrix for normalization
    arma::sp_mat D_inv(conn.n_rows, conn.n_rows);
    for (arma::uword i = 0; i < conn.n_rows; i++) {
        D_inv(i, i) = 1.0 / row_sums(i);
    }
    
    return D_inv * conn;
}

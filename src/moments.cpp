// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

//' Compute First-Order Moments (Neighbor Averaging)
//' 
//' @description Compute neighbor-averaged expression values.
//' 
//' @param X Expression matrix (cells x genes)
//' @param indices Neighbor indices (1-based)
//' @param weights Neighbor weights (optional, uniform if NULL)
//' 
//' @return Smoothed expression matrix
//' @keywords internal
// [[Rcpp::export]]
arma::mat compute_moments_cpp(const arma::mat& X,
                               const arma::imat& indices,
                               const Rcpp::Nullable<arma::mat>& weights = R_NilValue) {
    
    int n_cells = X.n_rows;
    int n_genes = X.n_cols;
    int n_neighbors = indices.n_cols;
    
    arma::mat M(n_cells, n_genes, arma::fill::zeros);
    
    bool use_weights = weights.isNotNull();
    arma::mat W;
    if (use_weights) {
        W = Rcpp::as<arma::mat>(weights);
    }
    
    for (int i = 0; i < n_cells; i++) {
        double total_weight = 0.0;
        
        for (int k = 0; k < n_neighbors; k++) {
            int j = indices(i, k) - 1;  // Convert to 0-based
            
            if (j >= 0 && j < n_cells) {
                double w = use_weights ? W(i, k) : 1.0;
                M.row(i) += w * X.row(j);
                total_weight += w;
            }
        }
        
        if (total_weight > 0) {
            M.row(i) /= total_weight;
        }
    }
    
    return M;
}

//' Compute Second-Order Moments
//' 
//' @description Compute neighbor-averaged second-order moments:
//' M_ss = <s^2>, M_us = <u*s>, M_uu = <u^2>
//' 
//' @param S Spliced expression matrix (cells x genes)
//' @param U Unspliced expression matrix (cells x genes)
//' @param indices Neighbor indices (1-based)
//' @param weights Neighbor weights (optional)
//' 
//' @return List with Mss, Mus, Muu matrices
//' @keywords internal
// [[Rcpp::export]]
List compute_second_moments_cpp(const arma::mat& S,
                                 const arma::mat& U,
                                 const arma::imat& indices,
                                 const Rcpp::Nullable<arma::mat>& weights = R_NilValue) {
    
    int n_cells = S.n_rows;
    int n_genes = S.n_cols;
    int n_neighbors = indices.n_cols;
    
    arma::mat Mss(n_cells, n_genes, arma::fill::zeros);
    arma::mat Mus(n_cells, n_genes, arma::fill::zeros);
    arma::mat Muu(n_cells, n_genes, arma::fill::zeros);
    
    bool use_weights = weights.isNotNull();
    arma::mat W;
    if (use_weights) {
        W = Rcpp::as<arma::mat>(weights);
    }
    
    for (int i = 0; i < n_cells; i++) {
        double total_weight = 0.0;
        
        for (int k = 0; k < n_neighbors; k++) {
            int j = indices(i, k) - 1;
            
            if (j >= 0 && j < n_cells) {
                double w = use_weights ? W(i, k) : 1.0;
                
                Mss.row(i) += w * (S.row(j) % S.row(j));  // s^2
                Mus.row(i) += w * (U.row(j) % S.row(j));  // u*s
                Muu.row(i) += w * (U.row(j) % U.row(j));  // u^2
                
                total_weight += w;
            }
        }
        
        if (total_weight > 0) {
            Mss.row(i) /= total_weight;
            Mus.row(i) /= total_weight;
            Muu.row(i) /= total_weight;
        }
    }
    
    return List::create(
        Named("Mss") = Mss,
        Named("Mus") = Mus,
        Named("Muu") = Muu
    );
}

//' Row Normalize Matrix
//' 
//' @param X Matrix
//' @return Row-normalized matrix
//' @keywords internal
// [[Rcpp::export]]
arma::mat row_normalize_cpp(const arma::mat& X) {
    arma::mat X_norm = X;
    
    for (int i = 0; i < (int)X.n_rows; i++) {
        double row_sum = arma::accu(X.row(i));
        if (row_sum > 0) {
            X_norm.row(i) /= row_sum;
        }
    }
    
    return X_norm;
}

//' Sparse Matrix-Vector Multiplication
//' 
//' @description Multiply sparse matrix with dense vector efficiently.
//' 
//' @param x Sparse matrix values
//' @param i Row indices (0-based)
//' @param j Column indices (0-based)
//' @param v Dense vector
//' @param n_rows Number of rows
//' 
//' @return Result vector
//' @keywords internal
// [[Rcpp::export]]
arma::vec sparse_matvec_cpp(const arma::vec& x,
                             const arma::ivec& i,
                             const arma::ivec& j,
                             const arma::vec& v,
                             int n_rows) {
    
    arma::vec result(n_rows, arma::fill::zeros);
    
    for (int k = 0; k < (int)x.n_elem; k++) {
        result(i(k)) += x(k) * v(j(k));
    }
    
    return result;
}

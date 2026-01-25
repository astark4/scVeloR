// scVeloR_types.h - Common type definitions
// Author: Zaoqu Liu

#ifndef SCVELOR_TYPES_H
#define SCVELOR_TYPES_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Ensure RcppArmadillo is properly configured
// Disable OpenMP in Armadillo to avoid conflicts
#define ARMA_DONT_USE_OPENMP

// Type aliases for clarity
using arma::mat;
using arma::vec;
using arma::uvec;
using arma::sp_mat;
using arma::rowvec;
using arma::colvec;

#endif // SCVELOR_TYPES_H

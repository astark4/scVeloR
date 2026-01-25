// em_core.cpp - Core C++ functions for EM algorithm and dynamics
// Complete implementation matching scvelo Python

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// =====================================================================
// ODE Analytical Solutions (Splicing Dynamics)
// =====================================================================

//' Unspliced mRNA analytical solution (C++)
//' 
//' @description u(t) = u0 * exp(-beta*t) + alpha/beta * (1 - exp(-beta*t))
//' @param tau Time point(s)
//' @param u0 Initial unspliced
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @return Unspliced abundance at time tau
//' @export
// [[Rcpp::export]]
NumericVector unspliced_cpp(NumericVector tau, double u0, double alpha, double beta) {
  int n = tau.size();
  NumericVector result(n);
  
  // Handle zero beta
  if (beta == 0) beta = 1e-10;
  
  double alpha_beta = alpha / beta;
  
  for (int i = 0; i < n; i++) {
    double exp_bt = std::exp(-beta * tau[i]);
    result[i] = u0 * exp_bt + alpha_beta * (1.0 - exp_bt);
  }
  
  return result;
}

//' Spliced mRNA analytical solution (C++)
//' 
//' @description Full analytical solution for spliced mRNA
//' @param tau Time point(s)
//' @param s0 Initial spliced
//' @param u0 Initial unspliced
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' @return Spliced abundance at time tau
//' @export
// [[Rcpp::export]]
NumericVector spliced_cpp(NumericVector tau, double s0, double u0,
                          double alpha, double beta, double gamma) {
  int n = tau.size();
  NumericVector result(n);
  
  // Handle edge cases
  if (beta == 0) beta = 1e-10;
  if (gamma == 0) gamma = 1e-10;
  
  double alpha_gamma = alpha / gamma;
  double denom = gamma - beta;
  
  // Handle case where gamma ≈ beta (use L'Hopital's rule)
  bool use_limit = std::abs(denom) < 1e-10;
  
  if (use_limit) {
    for (int i = 0; i < n; i++) {
      double t = tau[i];
      double exp_gt = std::exp(-gamma * t);
      double c_val = (alpha - u0 * beta) * t * exp_gt;
      result[i] = s0 * exp_gt + alpha_gamma * (1.0 - exp_gt) + c_val;
    }
  } else {
    double c_val = (alpha - u0 * beta) / denom;
    
    for (int i = 0; i < n; i++) {
      double t = tau[i];
      double exp_bt = std::exp(-beta * t);
      double exp_gt = std::exp(-gamma * t);
      result[i] = s0 * exp_gt + alpha_gamma * (1.0 - exp_gt) + c_val * (exp_gt - exp_bt);
    }
  }
  
  return result;
}

//' mRNA dynamics - both components (C++)
//' 
//' @param tau Time point(s)
//' @param u0 Initial unspliced
//' @param s0 Initial spliced
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' @return Matrix with u and s columns
//' @export
// [[Rcpp::export]]
NumericMatrix mRNA_cpp(NumericVector tau, double u0, double s0,
                       double alpha, double beta, double gamma) {
  int n = tau.size();
  NumericMatrix result(n, 2);
  
  NumericVector u = unspliced_cpp(tau, u0, alpha, beta);
  NumericVector s = spliced_cpp(tau, s0, u0, alpha, beta, gamma);
  
  result(_, 0) = u;
  result(_, 1) = s;
  
  return result;
}

//' Inverse time from unspliced (C++)
//' 
//' @description tau = -1/beta * log((u - alpha/beta) / (u0 - alpha/beta))
//' @param u Unspliced abundance
//' @param u0 Initial unspliced
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @return Estimated time tau
//' @export
// [[Rcpp::export]]
NumericVector tau_inv_u_cpp(NumericVector u, double u0, double alpha, double beta) {
  int n = u.size();
  NumericVector result(n);
  
  if (beta == 0) beta = 1e-10;
  
  double u_inf = alpha / beta;
  double denom = u0 - u_inf;
  
  for (int i = 0; i < n; i++) {
    double ratio = (u[i] - u_inf) / denom;
    
    // Clip to valid range
    if (ratio <= 0) ratio = 1e-10;
    if (ratio >= 1) ratio = 1 - 1e-10;
    
    double tau = -1.0 / beta * std::log(ratio);
    result[i] = std::max(tau, 0.0);
  }
  
  return result;
}

// =====================================================================
// Vectorization Functions
// =====================================================================

//' Vectorize parameters for time-dependent computation (C++)
//' 
//' @param t Cell times
//' @param t_ Switching time
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' @param alpha_ Basal transcription (default 0)
//' @param u0 Initial unspliced (default 0)
//' @param s0 Initial spliced (default 0)
//' @return List with tau, alpha_vec, u0_vec, s0_vec, o
//' @export
// [[Rcpp::export]]
List vectorize_params_cpp(NumericVector t, double t_, double alpha, double beta, double gamma,
                          double alpha_ = 0.0, double u0 = 0.0, double s0 = 0.0) {
  int n = t.size();
  
  NumericVector tau(n);
  NumericVector alpha_vec(n);
  NumericVector u0_vec(n);
  NumericVector s0_vec(n);
  IntegerVector o(n);
  
  // Compute state at switching time
  double u0_ = unspliced_cpp(NumericVector::create(t_), u0, alpha, beta)[0];
  double s0_ = spliced_cpp(NumericVector::create(t_), s0, u0, alpha, beta, gamma)[0];
  
  for (int i = 0; i < n; i++) {
    bool is_on = t[i] < t_;
    o[i] = is_on ? 1 : 0;
    
    if (is_on) {
      tau[i] = t[i];
      alpha_vec[i] = alpha;
      u0_vec[i] = u0;
      s0_vec[i] = s0;
    } else {
      tau[i] = t[i] - t_;
      alpha_vec[i] = alpha_;
      u0_vec[i] = u0_;
      s0_vec[i] = s0_;
    }
  }
  
  return List::create(
    Named("tau") = tau,
    Named("alpha") = alpha_vec,
    Named("u0") = u0_vec,
    Named("s0") = s0_vec,
    Named("o") = o
  );
}

// =====================================================================
// Loss and Distance Computation
// =====================================================================

//' Compute squared distances for EM (C++)
//' 
//' @param u Observed unspliced (scaled)
//' @param s Observed spliced
//' @param ut Model unspliced
//' @param st Model spliced
//' @param std_u Standard deviation of unspliced
//' @param std_s Standard deviation of spliced
//' @return Squared distance vector
//' @export
// [[Rcpp::export]]
NumericVector compute_squared_dist_cpp(NumericVector u, NumericVector s,
                                       NumericVector ut, NumericVector st,
                                       double std_u, double std_s) {
  int n = u.size();
  NumericVector result(n);
  
  for (int i = 0; i < n; i++) {
    double du = (u[i] - ut[i]) / std_u;
    double ds = (s[i] - st[i]) / std_s;
    result[i] = du * du + ds * ds;
  }
  
  return result;
}

//' Compute total loss (sum of squared errors) (C++)
//' 
//' @param u Observed unspliced
//' @param s Observed spliced
//' @param t Cell times
//' @param t_ Switching time
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' @param scaling Scaling factor
//' @param std_u Standard deviation of unspliced
//' @param std_s Standard deviation of spliced
//' @return Total loss
//' @export
// [[Rcpp::export]]
double compute_loss_cpp(NumericVector u, NumericVector s,
                        NumericVector t, double t_,
                        double alpha, double beta, double gamma,
                        double scaling, double std_u, double std_s) {
  int n = u.size();
  double loss = 0.0;
  
  // Scale u
  NumericVector u_scaled(n);
  for (int i = 0; i < n; i++) {
    u_scaled[i] = u[i] / scaling;
  }
  
  // Vectorize parameters
  List params = vectorize_params_cpp(t, t_, alpha, beta, gamma);
  NumericVector tau = params["tau"];
  NumericVector alpha_vec = params["alpha"];
  NumericVector u0_vec = params["u0"];
  NumericVector s0_vec = params["s0"];
  
  // Compute model predictions for each cell
  for (int i = 0; i < n; i++) {
    double ut = unspliced_cpp(NumericVector::create(tau[i]), u0_vec[i], alpha_vec[i], beta)[0];
    double st = spliced_cpp(NumericVector::create(tau[i]), s0_vec[i], u0_vec[i], alpha_vec[i], beta, gamma)[0];
    
    double du = (u_scaled[i] - ut) / std_u;
    double ds = (s[i] - st) / std_s;
    
    loss += du * du + ds * ds;
  }
  
  return loss;
}

// =====================================================================
// Time Point Assignment
// =====================================================================

//' Assign time points to cells (C++)
//' 
//' @description Project cells onto trajectory and assign time points
//' @param u Unspliced values
//' @param s Spliced values
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' @param t_ Switching time
//' @param u0_ Unspliced at switching
//' @param s0_ Spliced at switching
//' @param num Number of trajectory points to generate
//' @return List with tau (on), tau_ (off), and o (state indicator)
//' @export
// [[Rcpp::export]]
List assign_timepoints_cpp(NumericVector u, NumericVector s,
                           double alpha, double beta, double gamma,
                           double t_, double u0_, double s0_,
                           int num = 300) {
  int n = u.size();
  
  NumericVector tau(n);
  NumericVector tau_(n);
  IntegerVector o(n);
  
  // Generate on trajectory
  std::vector<double> t_on(num);
  std::vector<double> u_on(num);
  std::vector<double> s_on(num);
  
  for (int i = 0; i < num; i++) {
    t_on[i] = i * t_ / (num - 1);
    u_on[i] = unspliced_cpp(NumericVector::create(t_on[i]), 0, alpha, beta)[0];
    s_on[i] = spliced_cpp(NumericVector::create(t_on[i]), 0, 0, alpha, beta, gamma)[0];
  }
  
  // Estimate off phase max time
  double t0 = t_;
  if (beta > 0) {
    // tau_inv_u for decay to near zero
    t0 = std::max(t_, -std::log(0.01) / beta);
  }
  
  // Generate off trajectory
  std::vector<double> t_off(num);
  std::vector<double> u_off(num);
  std::vector<double> s_off(num);
  
  for (int i = 0; i < num; i++) {
    t_off[i] = i * t0 / (num - 1);
    u_off[i] = unspliced_cpp(NumericVector::create(t_off[i]), u0_, 0, beta)[0];
    s_off[i] = spliced_cpp(NumericVector::create(t_off[i]), s0_, u0_, 0, beta, gamma)[0];
  }
  
  // Find closest trajectory point for each cell
  for (int i = 0; i < n; i++) {
    double min_dist_on = std::numeric_limits<double>::max();
    int min_idx_on = 0;
    
    for (int j = 0; j < num; j++) {
      double du = u[i] - u_on[j];
      double ds = s[i] - s_on[j];
      double dist = du * du + ds * ds;
      if (dist < min_dist_on) {
        min_dist_on = dist;
        min_idx_on = j;
      }
    }
    
    double min_dist_off = std::numeric_limits<double>::max();
    int min_idx_off = 0;
    
    for (int j = 0; j < num; j++) {
      double du = u[i] - u_off[j];
      double ds = s[i] - s_off[j];
      double dist = du * du + ds * ds;
      if (dist < min_dist_off) {
        min_dist_off = dist;
        min_idx_off = j;
      }
    }
    
    if (min_dist_on <= min_dist_off) {
      tau[i] = t_on[min_idx_on];
      tau_[i] = 0;
      o[i] = 1;  // On phase
    } else {
      tau[i] = t_;
      tau_[i] = t_off[min_idx_off];
      o[i] = 0;  // Off phase
    }
  }
  
  return List::create(
    Named("tau") = tau,
    Named("tau_") = tau_,
    Named("o") = o
  );
}

// =====================================================================
// Likelihood and Divergence
// =====================================================================

//' Compute Gaussian log-likelihood (C++)
//' 
//' @param distx Squared distances
//' @param varx Variance
//' @return Log-likelihood
//' @export
// [[Rcpp::export]]
double gaussian_loglikelihood_cpp(NumericVector distx, double varx) {
  int n = distx.size();
  
  double sum_dist = 0.0;
  for (int i = 0; i < n; i++) {
    sum_dist += distx[i];
  }
  
  double ll = -0.5 / n * sum_dist / varx - 0.5 * std::log(2 * M_PI * varx);
  
  return ll;
}

//' Compute variance from residuals (C++)
//' 
//' @param distu Unspliced residuals
//' @param dists Spliced residuals
//' @return Estimated variance
//' @export
// [[Rcpp::export]]
double compute_variance_cpp(NumericVector distu, NumericVector dists) {
  int n = distu.size();
  
  // distx = distu^2 + dists^2
  // varx = E[distx] - E[sign(dists) * sqrt(distx)]^2
  
  double sum_distx = 0.0;
  double sum_signed_sqrt = 0.0;
  
  for (int i = 0; i < n; i++) {
    double distx = distu[i] * distu[i] + dists[i] * dists[i];
    double signed_sqrt = (dists[i] >= 0 ? 1 : -1) * std::sqrt(distx);
    
    sum_distx += distx;
    sum_signed_sqrt += signed_sqrt;
  }
  
  double mean_distx = sum_distx / n;
  double mean_signed_sqrt = sum_signed_sqrt / n;
  
  return mean_distx - mean_signed_sqrt * mean_signed_sqrt;
}

// =====================================================================
// Velocity Computation
// =====================================================================

//' Compute velocity from parameters (C++)
//' 
//' @description ds/dt = beta * u - gamma * s
//' @param u Unspliced values
//' @param s Spliced values
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' @param o State indicator (1=on, 0=off)
//' @param scaling Scaling factor
//' @return Velocity vector
//' @export
// [[Rcpp::export]]
NumericVector compute_velocity_cpp(NumericVector u, NumericVector s,
                                   double alpha, double beta, double gamma,
                                   IntegerVector o, double scaling = 1.0) {
  int n = u.size();
  NumericVector velocity(n);
  
  for (int i = 0; i < n; i++) {
    double u_scaled = u[i] / scaling;
    // velocity_s = beta * u - gamma * s
    double vel = beta * u_scaled - gamma * s[i];
    // Clip to not exceed current expression
    velocity[i] = std::max(vel, -s[i]);
  }
  
  return velocity;
}

//' Compute velocity unspliced component (C++)
//' 
//' @description du/dt = alpha - beta * u (only in on phase)
//' @param u Unspliced values
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param o State indicator
//' @param scaling Scaling factor
//' @return Velocity unspliced vector
//' @export
// [[Rcpp::export]]
NumericVector compute_velocity_u_cpp(NumericVector u, double alpha, double beta,
                                     IntegerVector o, double scaling = 1.0) {
  int n = u.size();
  NumericVector velocity(n);
  
  for (int i = 0; i < n; i++) {
    double u_scaled = u[i] / scaling;
    double alpha_i = (o[i] == 1) ? alpha : 0.0;
    double vel = (alpha_i - beta * u_scaled) * scaling;
    // Clip
    velocity[i] = std::max(vel, -u[i]);
  }
  
  return velocity;
}

// =====================================================================
// Utility Functions
// =====================================================================

//' Row-wise variance (C++)
//' 
//' @param X Matrix
//' @return Variance vector
//' @export
// [[Rcpp::export]]
NumericVector row_var_cpp(NumericMatrix X) {
  int n = X.nrow();
  int p = X.ncol();
  NumericVector result(n);
  
  for (int i = 0; i < n; i++) {
    double sum = 0.0;
    double sum_sq = 0.0;
    
    for (int j = 0; j < p; j++) {
      double val = X(i, j);
      sum += val;
      sum_sq += val * val;
    }
    
    double mean = sum / p;
    result[i] = sum_sq / p - mean * mean;
  }
  
  return result;
}

//' Column-wise variance (C++)
//' 
//' @param X Matrix
//' @return Variance vector
//' @export
// [[Rcpp::export]]
NumericVector col_var_cpp(NumericMatrix X) {
  int n = X.nrow();
  int p = X.ncol();
  NumericVector result(p);
  
  for (int j = 0; j < p; j++) {
    double sum = 0.0;
    double sum_sq = 0.0;
    
    for (int i = 0; i < n; i++) {
      double val = X(i, j);
      sum += val;
      sum_sq += val * val;
    }
    
    double mean = sum / n;
    result[j] = sum_sq / n - mean * mean;
  }
  
  return result;
}

//' Weighted correlation coefficient (C++)
//' 
//' @param x Vector
//' @param y Vector
//' @param weights Weights
//' @return Correlation coefficient
//' @export
// [[Rcpp::export]]
double weighted_corr_cpp(NumericVector x, NumericVector y, NumericVector weights) {
  int n = x.size();
  
  double sum_w = 0.0;
  double sum_wx = 0.0;
  double sum_wy = 0.0;
  
  for (int i = 0; i < n; i++) {
    sum_w += weights[i];
    sum_wx += weights[i] * x[i];
    sum_wy += weights[i] * y[i];
  }
  
  double mean_x = sum_wx / sum_w;
  double mean_y = sum_wy / sum_w;
  
  double cov_xy = 0.0;
  double var_x = 0.0;
  double var_y = 0.0;
  
  for (int i = 0; i < n; i++) {
    double dx = x[i] - mean_x;
    double dy = y[i] - mean_y;
    cov_xy += weights[i] * dx * dy;
    var_x += weights[i] * dx * dx;
    var_y += weights[i] * dy * dy;
  }
  
  return cov_xy / std::sqrt(var_x * var_y + 1e-10);
}

//' Safe divide (C++)
//' 
//' @param a Numerator
//' @param b Denominator
//' @param fill Fill value for division by zero
//' @return a/b with protection
//' @export
// [[Rcpp::export]]
NumericVector safe_divide_cpp(NumericVector a, NumericVector b, double fill = 0.0) {
  int n = a.size();
  NumericVector result(n);
  
  for (int i = 0; i < n; i++) {
    if (std::abs(b[i]) < 1e-10) {
      result[i] = fill;
    } else {
      result[i] = a[i] / b[i];
    }
  }
  
  return result;
}

//' Clamp values to range (C++)
//' 
//' @param x Values
//' @param min_val Minimum
//' @param max_val Maximum
//' @return Clamped values
//' @export
// [[Rcpp::export]]
NumericVector clamp_cpp(NumericVector x, double min_val, double max_val) {
  int n = x.size();
  NumericVector result(n);
  
  for (int i = 0; i < n; i++) {
    result[i] = std::max(min_val, std::min(max_val, x[i]));
  }
  
  return result;
}

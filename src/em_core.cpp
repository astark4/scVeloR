// em_core.cpp - Core functions for EM dynamical model
// Author: Zaoqu Liu
//
// Implements the splicing dynamics ODE solution and optimization helpers

#include "scVeloR_types.h"

using namespace Rcpp;

// Constants
const double EPS = 1e-10;
const double INF_REPLACEMENT = 0.0;

//' Safe division avoiding Inf
//' @keywords internal
inline double safe_div(double a, double b) {
    if (std::abs(b) < EPS) {
        return INF_REPLACEMENT;
    }
    return a / b;
}

//' Compute unspliced mRNA dynamics
//'
//' u(t) = u0 * exp(-beta*t) + alpha/beta * (1 - exp(-beta*t))
//'
//' @param tau Time points
//' @param u0 Initial unspliced abundance
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @return Unspliced abundance at time tau
//' @keywords internal
// [[Rcpp::export]]
NumericVector unspliced_dynamics_cpp(const NumericVector& tau,
                                      double u0,
                                      double alpha,
                                      double beta) {
    int n = tau.size();
    NumericVector result(n);
    
    double alpha_beta = safe_div(alpha, beta);
    
    for (int i = 0; i < n; i++) {
        double exp_bt = std::exp(-beta * tau[i]);
        result[i] = u0 * exp_bt + alpha_beta * (1.0 - exp_bt);
    }
    
    return result;
}

//' Compute spliced mRNA dynamics
//'
//' s(t) = s0*exp(-gamma*t) + alpha/gamma*(1-exp(-gamma*t)) + c*(exp(-gamma*t)-exp(-beta*t))
//' where c = (alpha - u0*beta) / (gamma - beta)
//'
//' @param tau Time points
//' @param s0 Initial spliced abundance
//' @param u0 Initial unspliced abundance
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' @return Spliced abundance at time tau
//' @keywords internal
// [[Rcpp::export]]
NumericVector spliced_dynamics_cpp(const NumericVector& tau,
                                    double s0,
                                    double u0,
                                    double alpha,
                                    double beta,
                                    double gamma) {
    int n = tau.size();
    NumericVector result(n);
    
    double c = safe_div(alpha - u0 * beta, gamma - beta);
    double alpha_gamma = safe_div(alpha, gamma);
    
    for (int i = 0; i < n; i++) {
        double exp_gt = std::exp(-gamma * tau[i]);
        double exp_bt = std::exp(-beta * tau[i]);
        result[i] = s0 * exp_gt + alpha_gamma * (1.0 - exp_gt) + c * (exp_gt - exp_bt);
    }
    
    return result;
}

//' Compute both unspliced and spliced dynamics
//'
//' @param tau Time points
//' @param u0 Initial unspliced abundance
//' @param s0 Initial spliced abundance
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' @return List with u and s vectors
//' @keywords internal
// [[Rcpp::export]]
List mRNA_dynamics_cpp(const NumericVector& tau,
                       double u0,
                       double s0,
                       double alpha,
                       double beta,
                       double gamma) {
    
    NumericVector u = unspliced_dynamics_cpp(tau, u0, alpha, beta);
    NumericVector s = spliced_dynamics_cpp(tau, s0, u0, alpha, beta, gamma);
    
    return List::create(
        Named("u") = u,
        Named("s") = s
    );
}

//' Inverse function to estimate tau from u
//'
//' tau = -1/beta * log((u - alpha/beta) / (u0 - alpha/beta))
//'
//' @param u Unspliced abundance
//' @param u0 Initial unspliced abundance
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @return Estimated time tau
//' @keywords internal
// [[Rcpp::export]]
NumericVector tau_inv_u_cpp(const NumericVector& u,
                            double u0,
                            double alpha,
                            double beta) {
    int n = u.size();
    NumericVector result(n);
    
    double u_inf = safe_div(alpha, beta);
    double denom = u0 - u_inf;
    
    for (int i = 0; i < n; i++) {
        double numer = u[i] - u_inf;
        double ratio = safe_div(numer, denom);
        
        if (ratio <= 0 || ratio > 1e10) {
            result[i] = 0.0;
        } else {
            result[i] = -safe_div(1.0, beta) * std::log(ratio);
        }
        
        // Clip to non-negative
        if (result[i] < 0) result[i] = 0.0;
    }
    
    return result;
}

//' Inverse function to estimate tau from u and s
//'
//' @param u Unspliced abundance
//' @param s Spliced abundance
//' @param u0 Initial unspliced abundance
//' @param s0 Initial spliced abundance
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' @return Estimated time tau
//' @keywords internal
// [[Rcpp::export]]
NumericVector tau_inv_us_cpp(const NumericVector& u,
                              const NumericVector& s,
                              double u0,
                              double s0,
                              double alpha,
                              double beta,
                              double gamma) {
    int n = u.size();
    NumericVector result(n);
    
    // Check which formula to use based on beta vs gamma
    bool use_u_only = (gamma >= beta);
    
    if (use_u_only) {
        result = tau_inv_u_cpp(u, u0, alpha, beta);
    } else {
        // Use spliced formula
        double beta_prime = safe_div(beta, gamma - beta);
        double x_inf = safe_div(alpha, gamma) - beta_prime * safe_div(alpha, beta);
        double x0 = s0 - beta_prime * u0 - x_inf;
        
        for (int i = 0; i < n; i++) {
            double x = s[i] - beta_prime * u[i] - x_inf;
            double ratio = safe_div(x, x0);
            
            if (ratio <= 0 || ratio > 1e10) {
                result[i] = 0.0;
            } else {
                result[i] = -safe_div(1.0, gamma) * std::log(ratio);
            }
            
            if (result[i] < 0) result[i] = 0.0;
        }
    }
    
    return result;
}

//' Compute mean squared error for dynamics fit
//'
//' @param u_obs Observed unspliced
//' @param s_obs Observed spliced
//' @param u_pred Predicted unspliced
//' @param s_pred Predicted spliced
//' @param weights Optional weights
//' @return MSE value
//' @keywords internal
// [[Rcpp::export]]
double compute_mse_cpp(const NumericVector& u_obs,
                       const NumericVector& s_obs,
                       const NumericVector& u_pred,
                       const NumericVector& s_pred,
                       const NumericVector& weights = NumericVector(0)) {
    int n = u_obs.size();
    
    bool use_weights = (weights.size() == n);
    double total_weight = use_weights ? 0.0 : (double)n;
    double mse = 0.0;
    
    for (int i = 0; i < n; i++) {
        double w = use_weights ? weights[i] : 1.0;
        if (use_weights) total_weight += w;
        
        double u_err = u_obs[i] - u_pred[i];
        double s_err = s_obs[i] - s_pred[i];
        
        mse += w * (u_err * u_err + s_err * s_err);
    }
    
    return mse / total_weight;
}

//' Assign time points by projection onto dynamics curve
//'
//' @param u_obs Observed unspliced (n_obs)
//' @param s_obs Observed spliced (n_obs)
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' @param t_switch Switching time
//' @param n_points Number of points on curve for projection
//' @return List with tau, tau_, and o (on/off indicator)
//' @keywords internal
// [[Rcpp::export]]
List assign_time_by_projection_cpp(const NumericVector& u_obs,
                                    const NumericVector& s_obs,
                                    double alpha,
                                    double beta,
                                    double gamma,
                                    double t_switch,
                                    int n_points = 300) {
    int n_obs = u_obs.size();
    
    // Generate time points for "on" phase
    NumericVector t_on(n_points);
    for (int i = 0; i < n_points; i++) {
        t_on[i] = t_switch * i / (n_points - 1);
    }
    
    // Compute "on" phase dynamics
    List on_dynamics = mRNA_dynamics_cpp(t_on, 0.0, 0.0, alpha, beta, gamma);
    NumericVector u_on = on_dynamics["u"];
    NumericVector s_on = on_dynamics["s"];
    
    // Get state at t_switch
    double u_switch = u_on[n_points - 1];
    double s_switch = s_on[n_points - 1];
    
    // Compute "off" phase - estimate off time from data
    double max_off_time = t_switch;  // Heuristic
    NumericVector t_off(n_points);
    for (int i = 0; i < n_points; i++) {
        t_off[i] = max_off_time * i / (n_points - 1);
    }
    
    // "Off" phase dynamics (alpha = 0)
    List off_dynamics = mRNA_dynamics_cpp(t_off, u_switch, s_switch, 0.0, beta, gamma);
    NumericVector u_off = off_dynamics["u"];
    NumericVector s_off = off_dynamics["s"];
    
    // Results
    NumericVector tau(n_obs);
    NumericVector tau_off(n_obs);
    IntegerVector o(n_obs);
    
    // For each observation, find closest point on curve
    for (int i = 0; i < n_obs; i++) {
        double min_dist_on = R_PosInf;
        int best_j_on = 0;
        
        double min_dist_off = R_PosInf;
        int best_j_off = 0;
        
        // Check "on" phase
        for (int j = 0; j < n_points; j++) {
            double du = u_obs[i] - u_on[j];
            double ds = s_obs[i] - s_on[j];
            double dist = du * du + ds * ds;
            
            if (dist < min_dist_on) {
                min_dist_on = dist;
                best_j_on = j;
            }
        }
        
        // Check "off" phase
        for (int j = 0; j < n_points; j++) {
            double du = u_obs[i] - u_off[j];
            double ds = s_obs[i] - s_off[j];
            double dist = du * du + ds * ds;
            
            if (dist < min_dist_off) {
                min_dist_off = dist;
                best_j_off = j;
            }
        }
        
        // Assign to closer phase
        if (min_dist_on <= min_dist_off) {
            tau[i] = t_on[best_j_on];
            tau_off[i] = 0.0;
            o[i] = 1;  // "on" phase
        } else {
            tau[i] = t_switch;
            tau_off[i] = t_off[best_j_off];
            o[i] = 0;  // "off" phase
        }
    }
    
    return List::create(
        Named("tau") = tau,
        Named("tau_") = tau_off,
        Named("o") = o
    );
}

//' Linear regression for gamma estimation
//'
//' @param u Unspliced values (predictor)
//' @param s Spliced values (response)
//' @param weights Optional weights
//' @return Estimated gamma (slope)
//' @keywords internal
// [[Rcpp::export]]
double linreg_gamma_cpp(const NumericVector& u,
                        const NumericVector& s,
                        const NumericVector& weights = NumericVector(0)) {
    int n = u.size();
    bool use_weights = (weights.size() == n);
    
    double sum_uu = 0.0;
    double sum_us = 0.0;
    
    for (int i = 0; i < n; i++) {
        double w = use_weights ? weights[i] : 1.0;
        sum_uu += w * u[i] * u[i];
        sum_us += w * u[i] * s[i];
    }
    
    return safe_div(sum_us, sum_uu);
}

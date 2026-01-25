// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
using namespace Rcpp;

//' Compute Unspliced mRNA Dynamics
//' 
//' @description Analytical solution for unspliced mRNA:
//' u(t) = u0 * exp(-beta*t) + alpha/beta * (1 - exp(-beta*t))
//' 
//' @param tau Time points
//' @param u0 Initial unspliced
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' 
//' @return Unspliced values
//' @keywords internal
// [[Rcpp::export]]
arma::vec unspliced_cpp(const arma::vec& tau,
                        double u0,
                        double alpha,
                        double beta) {
    
    if (beta < 1e-10) beta = 1e-10;
    
    arma::vec exp_bt = arma::exp(-beta * tau);
    return u0 * exp_bt + alpha / beta * (1.0 - exp_bt);
}

//' Compute Spliced mRNA Dynamics
//' 
//' @description Analytical solution for spliced mRNA.
//' 
//' @param tau Time points
//' @param s0 Initial spliced
//' @param u0 Initial unspliced
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' 
//' @return Spliced values
//' @keywords internal
// [[Rcpp::export]]
arma::vec spliced_cpp(const arma::vec& tau,
                      double s0,
                      double u0,
                      double alpha,
                      double beta,
                      double gamma) {
    
    if (beta < 1e-10) beta = 1e-10;
    if (gamma < 1e-10) gamma = 1e-10;
    
    arma::vec exp_bt = arma::exp(-beta * tau);
    arma::vec exp_gt = arma::exp(-gamma * tau);
    
    // c = (alpha - u0*beta) / (gamma - beta)
    double denom = gamma - beta;
    double c;
    if (std::abs(denom) < 1e-10) {
        // When gamma ≈ beta, use L'Hopital's rule approximation
        c = (alpha - u0 * beta) * tau(0);  // Approximation
    } else {
        c = (alpha - u0 * beta) / denom;
    }
    
    return s0 * exp_gt + alpha / gamma * (1.0 - exp_gt) + c * (exp_gt - exp_bt);
}

//' Compute mRNA Dynamics (Both Components)
//' 
//' @param tau Time points
//' @param u0 Initial unspliced
//' @param s0 Initial spliced
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' 
//' @return Matrix with u and s columns
//' @keywords internal
// [[Rcpp::export]]
arma::mat mRNA_cpp(const arma::vec& tau,
                   double u0,
                   double s0,
                   double alpha,
                   double beta,
                   double gamma) {
    
    arma::vec u = unspliced_cpp(tau, u0, alpha, beta);
    arma::vec s = spliced_cpp(tau, s0, u0, alpha, beta, gamma);
    
    arma::mat result(tau.n_elem, 2);
    result.col(0) = u;
    result.col(1) = s;
    
    return result;
}

//' Inverse Function: Estimate tau from u
//' 
//' @description tau = -1/beta * log((u - alpha/beta) / (u0 - alpha/beta))
//' 
//' @param u Unspliced values
//' @param u0 Initial unspliced
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' 
//' @return Estimated tau
//' @keywords internal
// [[Rcpp::export]]
arma::vec tau_inv_u_cpp(const arma::vec& u,
                        double u0,
                        double alpha,
                        double beta) {
    
    if (beta < 1e-10) beta = 1e-10;
    
    double u_inf = alpha / beta;
    
    arma::vec ratio = (u - u_inf) / (u0 - u_inf);
    
    // Clip ratio
    ratio = arma::clamp(ratio, 1e-10, 1.0 - 1e-10);
    
    arma::vec tau = -1.0 / beta * arma::log(ratio);
    
    // Ensure non-negative
    return arma::clamp(tau, 0.0, arma::datum::inf);
}

//' Compute Loss Function
//' 
//' @description Sum of squared errors between observed and predicted values.
//' 
//' @param u_obs Observed unspliced
//' @param s_obs Observed spliced
//' @param t Cell times
//' @param t_ Switching time
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' @param std_u Standard deviation of unspliced
//' @param std_s Standard deviation of spliced
//' 
//' @return Loss value
//' @keywords internal
// [[Rcpp::export]]
double compute_loss_cpp(const arma::vec& u_obs,
                        const arma::vec& s_obs,
                        const arma::vec& t,
                        double t_,
                        double alpha,
                        double beta,
                        double gamma,
                        double std_u,
                        double std_s) {
    
    if (alpha <= 0 || beta <= 0 || gamma <= 0 || t_ <= 0) {
        return 1e20;
    }
    
    int n = u_obs.n_elem;
    
    // Compute state at switching time
    arma::vec tau_switch(1);
    tau_switch(0) = t_;
    arma::mat switch_state = mRNA_cpp(tau_switch, 0, 0, alpha, beta, gamma);
    double u0_ = switch_state(0, 0);
    double s0_ = switch_state(0, 1);
    
    double total_loss = 0.0;
    
    for (int i = 0; i < n; i++) {
        double ti = t(i);
        double ut, st;
        
        if (ti < t_) {
            // On phase: starts from (0, 0)
            arma::vec tau_i(1);
            tau_i(0) = ti;
            arma::mat state = mRNA_cpp(tau_i, 0, 0, alpha, beta, gamma);
            ut = state(0, 0);
            st = state(0, 1);
        } else {
            // Off phase: starts from (u0_, s0_)
            arma::vec tau_i(1);
            tau_i(0) = ti - t_;
            arma::mat state = mRNA_cpp(tau_i, u0_, s0_, 0, beta, gamma);
            ut = state(0, 0);
            st = state(0, 1);
        }
        
        double u_diff = (u_obs(i) - ut) / std_u;
        double s_diff = (s_obs(i) - st) / std_s;
        
        total_loss += u_diff * u_diff + s_diff * s_diff;
    }
    
    return total_loss;
}

//' Assign Time Points to Cells
//' 
//' @description Project cells onto trajectory or use inverse formula.
//' 
//' @param u Unspliced values
//' @param s Spliced values
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' @param t_ Switching time
//' @param num_points Number of points for projection
//' 
//' @return List with tau, tau_, and o
//' @keywords internal
// [[Rcpp::export]]
List assign_timepoints_cpp(const arma::vec& u,
                           const arma::vec& s,
                           double alpha,
                           double beta,
                           double gamma,
                           double t_,
                           int num_points = 200) {
    
    int n = u.n_elem;
    
    // Compute switch state
    arma::vec tau_switch(1);
    tau_switch(0) = t_;
    arma::mat switch_state = mRNA_cpp(tau_switch, 0, 0, alpha, beta, gamma);
    double u0_ = switch_state(0, 0);
    double s0_ = switch_state(0, 1);
    
    // Generate trajectory points
    // On phase
    arma::vec t_on = arma::linspace(0, t_, num_points);
    arma::mat traj_on(num_points, 2);
    for (int k = 0; k < num_points; k++) {
        arma::vec tau_k(1);
        tau_k(0) = t_on(k);
        arma::mat state = mRNA_cpp(tau_k, 0, 0, alpha, beta, gamma);
        traj_on(k, 0) = state(0, 0);
        traj_on(k, 1) = state(0, 1);
    }
    
    // Off phase - estimate max time
    double t_off_max = t_;  // Default
    arma::vec u_pos = u(arma::find(s > 0));
    if (u_pos.n_elem > 0) {
        double u_min = arma::min(u_pos);
        arma::vec tau_est = tau_inv_u_cpp(arma::vec({u_min}), u0_, 0, beta);
        if (std::isfinite(tau_est(0)) && tau_est(0) > 0) {
            t_off_max = tau_est(0);
        }
    }
    
    arma::vec t_off = arma::linspace(0, t_off_max, num_points);
    arma::mat traj_off(num_points, 2);
    for (int k = 0; k < num_points; k++) {
        arma::vec tau_k(1);
        tau_k(0) = t_off(k);
        arma::mat state = mRNA_cpp(tau_k, u0_, s0_, 0, beta, gamma);
        traj_off(k, 0) = state(0, 0);
        traj_off(k, 1) = state(0, 1);
    }
    
    // Assign each cell
    arma::vec tau(n);
    arma::vec tau_(n);
    arma::ivec o(n);
    
    for (int i = 0; i < n; i++) {
        double ui = u(i);
        double si = s(i);
        
        // Distance to on trajectory
        double min_dist_on = 1e20;
        int min_idx_on = 0;
        
        for (int k = 0; k < num_points; k++) {
            double dist = std::pow(ui - traj_on(k, 0), 2) + 
                          std::pow(si - traj_on(k, 1), 2);
            if (dist < min_dist_on) {
                min_dist_on = dist;
                min_idx_on = k;
            }
        }
        
        // Distance to off trajectory
        double min_dist_off = 1e20;
        int min_idx_off = 0;
        
        for (int k = 0; k < num_points; k++) {
            double dist = std::pow(ui - traj_off(k, 0), 2) + 
                          std::pow(si - traj_off(k, 1), 2);
            if (dist < min_dist_off) {
                min_dist_off = dist;
                min_idx_off = k;
            }
        }
        
        if (min_dist_on <= min_dist_off) {
            tau(i) = t_on(min_idx_on);
            tau_(i) = 0;
            o(i) = 1;  // On phase
        } else {
            tau(i) = t_;
            tau_(i) = t_off(min_idx_off);
            o(i) = 0;  // Off phase
        }
    }
    
    return List::create(
        Named("tau") = tau,
        Named("tau_") = tau_,
        Named("o") = o
    );
}

//' Vectorize Parameters Based on Cell Phase
//' 
//' @param t Cell times
//' @param t_ Switching time
//' @param alpha Transcription rate
//' @param beta Splicing rate
//' @param gamma Degradation rate
//' 
//' @return List with tau, alpha_vec, u0, s0
//' @keywords internal
// [[Rcpp::export]]
List vectorize_params_cpp(const arma::vec& t,
                          double t_,
                          double alpha,
                          double beta,
                          double gamma) {
    
    int n = t.n_elem;
    
    // Compute switch state
    arma::vec tau_switch(1);
    tau_switch(0) = t_;
    arma::mat switch_state = mRNA_cpp(tau_switch, 0, 0, alpha, beta, gamma);
    double u0_ = switch_state(0, 0);
    double s0_ = switch_state(0, 1);
    
    arma::vec tau(n);
    arma::vec alpha_vec(n);
    arma::vec u0_vec(n);
    arma::vec s0_vec(n);
    
    for (int i = 0; i < n; i++) {
        if (t(i) < t_) {
            // On phase
            tau(i) = t(i);
            alpha_vec(i) = alpha;
            u0_vec(i) = 0;
            s0_vec(i) = 0;
        } else {
            // Off phase
            tau(i) = t(i) - t_;
            alpha_vec(i) = 0;
            u0_vec(i) = u0_;
            s0_vec(i) = s0_;
        }
    }
    
    return List::create(
        Named("tau") = tau,
        Named("alpha") = alpha_vec,
        Named("u0") = u0_vec,
        Named("s0") = s0_vec
    );
}

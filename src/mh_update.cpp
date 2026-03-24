// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h> // for matrix/vector types

// NOTE: OpenMP is not configured yet
// The guards are here for future multi-threaded parallel execution
// When OpenMP is available, the parallel phase of mh_update_cpp will
// automatically distribute the N observations across CPU cores
#ifdef _OPENMP
#include <omp.h>
#endif
#include <cmath>

// Helper function that computes normal-approximation-to-Beta proposal parameters
// for a single (observation, dimension) pair.
// Uses R::digamma / R::trigamma (R's internal C functions) which are NOT THREAD-SAFE
// do not call this in a parallel execution region.
static inline void normbeta_params(
    double y_i,      // count for category i
    double y_ref,    // count for reference (last) category
    double mu_i,     // mean for dimension i
    double sigma_ii, // Sigma[i,i] (already scaled by mh_scale in R)
    double &muprop,
    double &sigprop)
{
    double exp_mu = std::exp(mu_i);
    double alpha = (1.0 + exp_mu) / sigma_ii - exp_mu / (1.0 + exp_mu);
    if (alpha < 1e-3) // safety floor preventing numerical blowup
        alpha = 1e-3;
    double beta_val = alpha * std::exp(-mu_i);
    double a_star = y_i + alpha;
    double b_star = y_ref + beta_val;
    muprop = R::digamma(a_star) - R::digamma(b_star);
    sigprop = std::sqrt(R::trigamma(a_star) + R::trigamma(b_star));
}

//  C++ replacement for the three apply() calls in the
//  normbeta Metropolis-Hastings latent-W update shared by FMLN and MMLN.
//
//  @param W     N x d matrix of current latent log-ratio values.
//  @param Y     N x (d+1) matrix of multinomial counts.
//  @param Mu    N x d matrix of linear predictors.
//  @param Sigma d x d proposal covariance.
//
//  @return A named list with three elements:
//      - W_new:        N x d matrix of proposed latent values.
//      - log_q_old:    Length-N vector of proposal log-densities at the current W.
//      - log_q_new:    Length-N vector of proposal log-densities at W_new.
//
// [[Rcpp::export]]
Rcpp::List mh_update_cpp(
    const arma::mat &W,
    const arma::mat &Y,
    const arma::mat &Mu,
    const arma::mat &Sigma)
{
    const int N = static_cast<int>(W.n_rows); // num observations
    const int d = static_cast<int>(W.n_cols); // num outcome categories - 1

    const arma::vec sigma_diag = Sigma.diag(); // only diagonal elements used (no need to work on full matrix)

    // Pre-computation phase (not thread-safe)
    // operational matrices
    arma::mat muprop_mat(N, d);
    arma::mat sigprop_mat(N, d);
    arma::mat Z_mat(N, d);

    for (int i = 0; i < N; i++)
    {
        const double y_ref = Y(i, d); // reference-category count (last column)
        for (int j = 0; j < d; j++)
        {
            // building proposal distribution parameters
            normbeta_params(
                Y(i, j), y_ref, Mu(i, j), sigma_diag(j),
                muprop_mat(i, j), sigprop_mat(i, j));
            Z_mat(i, j) = R::rnorm(0.0, 1.0);
        }
    }

    // Parallel computation phase (pure C++ math, thread-safe)
    // output variables
    arma::mat W_new(N, d);
    arma::vec log_q_old(N);
    arma::vec log_q_new(N);

    const double log_sqrt_2pi = 0.5 * std::log(2.0 * M_PI); // precompute to be reused in loop

// Only runs in parallel if running machine has OpenMP setup,
// ignores and runs single threaded if not
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i = 0; i < N; i++)
    {
        double lq_old = 0.0;
        double lq_new = 0.0;
        for (int j = 0; j < d; j++)
        {
            // pre-computed proposal mean, sd, and random draw for this obs/dim
            const double mp = muprop_mat(i, j);
            const double sp = sigprop_mat(i, j);
            const double z = Z_mat(i, j);

            // proposed value: W_new = mp + sp * z
            // generate proposed value by shifting proposed mean by z
            W_new(i, j) = mp + sp * z;

            // log dnorm(W_old | mp, sp)
            // distance from current value (old) from the proposed center
            const double dev_old = (W(i, j) - mp) / sp;
            lq_old += -log_sqrt_2pi - std::log(sp) - 0.5 * dev_old * dev_old;

            // log dnorm(W_new | mp, sp) — (W_new - mp)/sp == z by construction
            // formula for lq_new:
            // -log_sqrt_2pi - log(sp) - 0.5 * ((W_new - mp) / sp)^2
            //
            // Since W_new = mp + sp * z;
            // (W_new - mp) / sp = (mp + sp * z - mp) / sp = sp*z / sp = z
            //
            // so formula simplifies to:
            // -log_sqrt_2pi - log(sp) - 0.5 * z^2
            lq_new += -log_sqrt_2pi - std::log(sp) - 0.5 * z * z;
        }
        // to be used in MMLN/FMLN for MH accept/reject decision
        log_q_old(i) = lq_old;
        log_q_new(i) = lq_new;
    }

    return Rcpp::List::create(
        Rcpp::Named("W_new") = W_new,
        Rcpp::Named("log_q_old") = log_q_old,
        Rcpp::Named("log_q_new") = log_q_new);
}

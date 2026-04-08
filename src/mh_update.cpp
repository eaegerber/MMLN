// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

// Helper function that computes normal-approximation-to-Beta proposal parameters
// for a single (observation, dimension) pair.
// Uses R::digamma / R::trigamma which are R's internal C functions.
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
Rcpp::List mh_update_normbeta_cpp(
    const arma::mat &W,
    const arma::mat &Y,
    const arma::mat &Mu,
    const arma::mat &Sigma)
{
    const int N = static_cast<int>(W.n_rows); // num observations
    const int d = static_cast<int>(W.n_cols); // num outcome categories - 1

    const arma::vec sigma_diag = Sigma.diag(); // only diagonal elements used

    // Compute proposal parameters and random draws
    arma::mat muprop_mat(N, d);
    arma::mat sigprop_mat(N, d);
    arma::mat Z_mat(N, d);

    for (int i = 0; i < N; i++)
    {
        const double y_ref = Y(i, d); // reference-category count (last column)
        for (int j = 0; j < d; j++)
        {
            normbeta_params(
                Y(i, j), y_ref, Mu(i, j), sigma_diag(j),
                muprop_mat(i, j), sigprop_mat(i, j));
            Z_mat(i, j) = R::rnorm(0.0, 1.0);
        }
    }

    // Compute proposals and log-densities
    arma::mat W_new(N, d);
    arma::vec log_q_old(N);
    arma::vec log_q_new(N);

    const double log_sqrt_2pi = 0.5 * std::log(2.0 * M_PI);

    for (int i = 0; i < N; i++)
    {
        double lq_old = 0.0;
        double lq_new = 0.0;
        for (int j = 0; j < d; j++)
        {
            const double mp = muprop_mat(i, j);
            const double sp = sigprop_mat(i, j);
            const double z = Z_mat(i, j);

            // proposed value: W_new = mp + sp * z
            W_new(i, j) = mp + sp * z;

            // log dnorm(W_old | mp, sp)
            const double dev_old = (W(i, j) - mp) / sp;
            lq_old += -log_sqrt_2pi - std::log(sp) - 0.5 * dev_old * dev_old;

            // log dnorm(W_new | mp, sp)
            // Since W_new = mp + sp * z, we have (W_new - mp)/sp = z, so:
            lq_new += -log_sqrt_2pi - std::log(sp) - 0.5 * z * z;
        }
        log_q_old(i) = lq_old;
        log_q_new(i) = lq_new;
    }

    return Rcpp::List::create(
        Rcpp::Named("W_new") = W_new,
        Rcpp::Named("log_q_old") = log_q_old,
        Rcpp::Named("log_q_new") = log_q_new);
}

// Compute the Beta proposal distribution parameters for a single (obs, dim) pair.
// These are shared by the proposal draw and the log-likelihood evaluation.
//
// The idea: we approximate the posterior for the sigmoid-transformed latent
// variable pstar = sigmoid(w) using a Beta distribution. The shape parameters
// come from matching moments to the prior + data:
//
//   alpha = max( (1 + e^mu) / sigma_jj - e^mu / (1 + e^mu), 1e-3 )
//   beta  = alpha * e^{-mu}
//
// Then we fold in the observed counts to get the posterior-like shapes:
//   a_star = y_j   + alpha     (category j count + prior shape)
//   b_star = y_ref + beta      (reference category count + prior shape)
//
// These a_star, b_star are used both for drawing Beta proposals
// and evaluating their log-densities.
static inline void beta_params(
    double y_j,      // count for category j
    double y_ref,    // count for reference (last) category
    double mu_j,     // mean for dimension j
    double sigma_jj, // Sigma[j,j] (already scaled by mh_scale in R)
    double &a_star,
    double &b_star)
{
    double exp_mu = std::exp(mu_j);
    double alpha = (1.0 + exp_mu) / sigma_jj - exp_mu / (1.0 + exp_mu);
    if (alpha < 1e-3)
        alpha = 1e-3;
    double beta_val = alpha * std::exp(-mu_j);
    a_star = y_j + alpha;
    b_star = y_ref + beta_val;
}

// C++ replacement for the six apply() calls in the "beta" proposal
// Metropolis-Hastings latent-W update shared by FMLN and MMLN.
//
// The beta proposal works in sigmoid space (pstar = sigmoid(w)):
//   1. Map current W to pstar via sigmoid: pstar = e^w / (1 + e^w)
//   2. Draw new pstar from Beta(a_star, b_star)
//   3. Map back to W space via logit:       w = log(pstar / (1 - pstar))
//   4. Evaluate proposal log-densities at both old and new pstar,
//      including the Jacobian correction for the logit transform
//
// The log-density of the proposal at a point pstar is:
//   sum_j [ log Beta(pstar_j | a_star_j, b_star_j) + log(pstar_j) + log(1 - pstar_j) ]
//
// where the log(pstar) + log(1-pstar) terms are the Jacobian from the
// change of variables w -> pstar = sigmoid(w). This ensures the MH ratio
// correctly accounts for the nonlinear transformation.
//
// Expanding the Beta log-density:
//   log Beta(x | a, b) = (a-1)*log(x) + (b-1)*log(1-x) - lbeta(a, b)
//
// So each dimension's contribution to the full log q is:
//   (a_star - 1)*log(pstar) + (b_star - 1)*log(1-pstar) - lbeta(a_star, b_star)
//   + log(pstar) + log(1-pstar)
//
// Combining the log(pstar) and log(1-pstar) terms:
//   = a_star*log(pstar) + b_star*log(1-pstar) - lbeta(a_star, b_star)
//
// This simplified form is what we compute below.
//
// @param W     N x d matrix of current latent log-ratio values.
// @param Y     N x (d+1) matrix of multinomial counts.
// @param Mu    N x d matrix of linear predictors.
// @param Sigma d x d proposal covariance (already scaled by mh_scale).
//
// @return Named list: W_new (N x d), log_q_old (N), log_q_new (N).
//
// [[Rcpp::export]]
Rcpp::List mh_update_beta_cpp(
    const arma::mat &W,
    const arma::mat &Y,
    const arma::mat &Mu,
    const arma::mat &Sigma)
{
    const int N = static_cast<int>(W.n_rows);
    const int d = static_cast<int>(W.n_cols);

    const arma::vec sigma_diag = Sigma.diag();

    // Precompute shape parameters, sigmoid of current W, Beta draws,
    // and lbeta values. R::rbeta and R::lbeta are R C internals.
    arma::mat a_star_mat(N, d);
    arma::mat b_star_mat(N, d);
    arma::mat Pstar_old(N, d);
    arma::mat Pstar_new(N, d);
    arma::mat lbeta_mat(N, d);

    for (int i = 0; i < N; i++)
    {
        const double y_ref = Y(i, d);
        for (int j = 0; j < d; j++)
        {
            beta_params(
                Y(i, j), y_ref, Mu(i, j), sigma_diag(j),
                a_star_mat(i, j), b_star_mat(i, j));

            // sigmoid transform of current W
            double exp_w = std::exp(W(i, j));
            Pstar_old(i, j) = exp_w / (1.0 + exp_w);

            // draw from Beta proposal
            Pstar_new(i, j) = R::rbeta(a_star_mat(i, j), b_star_mat(i, j));

            // precompute lbeta for use in log-density calculation
            lbeta_mat(i, j) = R::lbeta(a_star_mat(i, j), b_star_mat(i, j));
        }
    }

    // Compute W_new (logit of Pstar_new) and log-densities
    arma::mat W_new(N, d);
    arma::vec log_q_old(N);
    arma::vec log_q_new(N);

    for (int i = 0; i < N; i++)
    {
        double lq_old = 0.0;
        double lq_new = 0.0;
        for (int j = 0; j < d; j++)
        {
            const double a = a_star_mat(i, j);
            const double b = b_star_mat(i, j);
            const double lb = lbeta_mat(i, j);

            const double ps_old = Pstar_old(i, j);
            const double ps_new = Pstar_new(i, j);

            // logit transform: map new pstar back to W space
            // w = log(pstar / (1 - pstar))
            W_new(i, j) = std::log(ps_new / (1.0 - ps_new));

            // Log proposal density (see derivation in function header):
            //   log q(pstar) = a_star * log(pstar) + b_star * log(1 - pstar) - lbeta(a_star, b_star)
            lq_old += a * std::log(ps_old) + b * std::log(1.0 - ps_old) - lb;
            lq_new += a * std::log(ps_new) + b * std::log(1.0 - ps_new) - lb;
        }
        log_q_old(i) = lq_old;
        log_q_new(i) = lq_new;
    }

    return Rcpp::List::create(
        Rcpp::Named("W_new") = W_new,
        Rcpp::Named("log_q_old") = log_q_old,
        Rcpp::Named("log_q_new") = log_q_new);
}

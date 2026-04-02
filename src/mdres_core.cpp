// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// Helper: Smithson-Verkuilen zero-count correction applied row-wise.
// Matches R compress_counts(): rows with any zero get
//   y_new = (y * (rowsum - 1) + 1/K) / rowsum
static arma::mat compress_counts_cpp(const arma::mat &Y)
{
    const int N = static_cast<int>(Y.n_rows);
    const int K = static_cast<int>(Y.n_cols);
    arma::mat out = Y;
    for (int i = 0; i < N; ++i)
    {
        bool has_zero = false;
        for (int j = 0; j < K; ++j)
        {
            if (Y(i, j) == 0.0)
            {
                has_zero = true;
                break;
            }
        }
        if (has_zero)
        {
            const double rs = arma::sum(Y.row(i));
            out.row(i) = ((Y.row(i) * (rs - 1.0)) + (1.0 / K)) / rs;
        }
    }
    return out;
}

// Helper: Additive log-ratio transform: N x J -> N x (J-1).
// Last category is reference: W[i,j] = log(P[i,j] / P[i,J-1]).
// Matches R alr().
static arma::mat alr_cpp(const arma::mat &P)
{
    const int N = static_cast<int>(P.n_rows);
    const int J = static_cast<int>(P.n_cols);
    arma::mat W(N, J - 1);
    for (int i = 0; i < N; ++i)
    {
        const double ref = P(i, J - 1);
        for (int j = 0; j < J - 1; ++j)
            W(i, j) = std::log(P(i, j) / ref);
    }
    return W;
}

//  C++ replacement for full MDres pipeline
//
//  Per-observation sample covariance is computed via arma::cov() (divides
//  by P-1, matching R cov()). A diagonal jitter of 1e-8 is always applied
//  before inversion for numerical stability.
//
//  @param Y           N x J matrix of observed multinomial counts.
//  @param Y_pred_list List of P matrices (each N x J) of posterior
//                     predictive counts.
//
//  @return Numeric vector of length N containing z-residuals (class set
//          by the R wrapper MDres()).
//
// [[Rcpp::export]]
Rcpp::NumericVector mdres_core_cpp(
    const arma::mat &Y,
    const Rcpp::List &Y_pred_list)
{
    const int N = static_cast<int>(Y.n_rows);           // num observations
    const int J = static_cast<int>(Y.n_cols);           // num outcome categories
    const int d = J - 1;                                // ALR dimension
    const int P = static_cast<int>(Y_pred_list.size()); // num posterior samples

    // compress and ALR-transform observed counts
    const arma::mat alr_obs = alr_cpp(compress_counts_cpp(Y)); // N x d

    // build pred_cube[obs, dim, sample]
    // pred_cube(i, j, k) = ALR-transformed value for observation i,
    //                       dimension j, posterior sample k
    arma::cube pred_cube(N, d, P);
    for (int k = 0; k < P; ++k)
    {
        const arma::mat Yk = Rcpp::as<arma::mat>(Y_pred_list[k]);
        pred_cube.slice(k) = alr_cpp(compress_counts_cpp(Yk)); // N x d
    }

    // per-observation Mahalanobis distances + ECDF quantile
    int n_singular = 0;
    Rcpp::NumericVector z_resids(N);

    for (int i = 0; i < N; ++i)
    {

        // collect P ALR-transformed predictive samples into P x d matrix
        arma::mat samp(P, d);
        for (int k = 0; k < P; ++k)
        {
            samp.row(k) = pred_cube.slice(k).row(i);
        }

        // per-observation sample covariance via arma::cov (divides by P-1,
        // matching R cov())
        // replaces R: Sigma_all[[i]] <- cov(t(pred_array[i,,]))
        arma::mat S = arma::cov(samp); // d x d

        // singularity check: count obs where min eigenvalue < 1e-8
        const arma::vec ev = arma::eig_sym(S);
        if (ev.min() < 1e-8)
        {
            ++n_singular;
        }

        // diagonal jitter for numerical stability before Cholesky
        S.diag() += 1e-8;

        // map S into Eigen (arma and Eigen both use column-major; zero copy)
        const Eigen::MatrixXd S_eig =
            Eigen::Map<const Eigen::MatrixXd>(S.memptr(), d, d);

        // Cholesky factorization: S = L L^T
        // replaces R: solve(Sigma_all[[i]]) and chol2inv(chol(...))
        Eigen::LLT<Eigen::MatrixXd> llt(S_eig);
        if (llt.info() != Eigen::Success)
        {
            Rcpp::warning(
                "LLT decomposition failed for observation %d; "
                "result set to NA.",
                i + 1);
            z_resids[i] = NA_REAL;
            continue;
        }

        // column means of predictive samples
        // replaces R: mu_all <- apply(pred_array, 1:2, mean)
        const arma::rowvec mu = arma::mean(samp, 0); // 1 x d

        // map samp and mu into Eigen views (zero copy)
        const Eigen::MatrixXd samp_eig =
            Eigen::Map<const Eigen::MatrixXd>(samp.memptr(), P, d);
        const Eigen::RowVectorXd mu_eig =
            Eigen::Map<const Eigen::RowVectorXd>(mu.memptr(), d);

        // build centered (P+1) x d matrix:
        //   row 0     = alr_obs[i,] - mu   (the observed point)
        //   rows 1..P = samp - mu           (the posterior predictive points)
        //
        // the .rowwise() broadcast replaces R sweep(pred_array, 2, mu)
        // no intermediate matrix allocated
        Eigen::MatrixXd centered(P + 1, d);
        for (int j = 0; j < d; ++j)
        {
            centered(0, j) = alr_obs(i, j) - mu(j);
        }
        centered.bottomRows(P).noalias() = samp_eig.rowwise() - mu_eig;

        // row-wise Mahalanobis distances via a single triangular solve:
        //   md_j = || L^{-1} centered_j ||^2
        //
        // solve L * Z = centered^T  ->  Z is d x (P+1)
        // then squared column norms give the Mahalanobis distances
        // Eigen expression template avoids intermediate allocation
        //
        // replaces R: apply(obsi, 1, mahalanobis, center = mu, cov = S)
        const Eigen::MatrixXd Z =
            llt.matrixL().solve(centered.transpose());
        const Eigen::RowVectorXd mds = Z.colwise().squaredNorm(); // 1 x (P+1)

        const double obs_md = mds(0);

        // sort predictive distances for ECDF lookup
        // replaces R: sort(post_vals) and order()
        // std::sort is O(P log P), then binary search is O(log P) per query
        std::vector<double> sorted_post(P);
        for (int k = 0; k < P; ++k)
        {
            sorted_post[k] = mds(k + 1);
        }
        std::sort(sorted_post.begin(), sorted_post.end());

        const double min_val = sorted_post.front();
        const double max_val = sorted_post.back();

        // ecdf(x) = #{sorted_post <= x} / P
        // O(log P) via std::upper_bound instead of R's linear ecdf()
        auto ecdf_val = [&](double x) -> double
        {
            return static_cast<double>(
                       std::upper_bound(sorted_post.begin(), sorted_post.end(), x) - sorted_post.begin()) /
                   P;
        };

        // randomised-quantile CDF step (matches R u_resids logic exactly):
        // find the ECDF interval [minpct, maxpct] around the observed
        // Mahalanobis distance, then draw uniform within that interval
        double minpct, maxpct;
        if (obs_md <= min_val)
        {
            minpct = 0.0;
            maxpct = ecdf_val(min_val);
            if (maxpct == 0.0)
                maxpct = static_cast<double>(P - 1) / P;
        }
        else if (obs_md > max_val)
        {
            minpct = ecdf_val(max_val); // == 1.0
            maxpct = 1.0;
            if (minpct == 1.0)
                minpct = static_cast<double>(P - 1) / P;
        }
        else
        {
            // lower_idx: last 0-based index where sorted_post[k] < obs_md
            const int lower_idx = static_cast<int>(
                                      std::lower_bound(sorted_post.begin(), sorted_post.end(), obs_md) - sorted_post.begin()) -
                                  1;
            minpct = ecdf_val(sorted_post[lower_idx]);
            maxpct = ecdf_val(obs_md);
            if (minpct == 1.0)
                minpct = static_cast<double>(P - 1) / P;
            if (maxpct == 0.0)
                maxpct = static_cast<double>(P - 1) / P;
        }

        // draw uniform CDF quantile, then invert to normal scale
        // to be used in R wrapper for QQ diagnostic plots
        const double u_i = R::runif(minpct, maxpct);
        z_resids[i] = R::qnorm(u_i, 0.0, 1.0, 1, 0);
    }

    if (n_singular > 0)
    {
        Rcpp::warning(
            "Covariance matrix near-singular (min eigenvalue < 1e-8) for "
            "%d observation(s); 1e-8 diagonal jitter was applied.",
            n_singular);
    }

    return z_resids;
}
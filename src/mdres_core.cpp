// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <algorithm>
#include <vector>

using namespace Rcpp;

// ── internal helpers ────────────────────────────────────────────────────────

// Smithson-Verkuilen zero-count correction applied row-wise.
// Matches R compress_counts(): rows with any zero get
//   y_new = (y * (rowsum - 1) + 1/K) / rowsum
static arma::mat compress_counts_cpp(const arma::mat& Y) {
    const int N = static_cast<int>(Y.n_rows);
    const int K = static_cast<int>(Y.n_cols);
    arma::mat out = Y;
    for (int i = 0; i < N; ++i) {
        bool has_zero = false;
        for (int j = 0; j < K; ++j) {
            if (Y(i, j) == 0.0) { has_zero = true; break; }
        }
        if (has_zero) {
            const double rs = arma::sum(Y.row(i));
            out.row(i) = ((Y.row(i) * (rs - 1.0)) + (1.0 / K)) / rs;
        }
    }
    return out;
}

// Additive log-ratio transform: N×J  →  N×(J−1).
// Last category is reference: W[i,j] = log(P[i,j] / P[i,J-1]).
// Matches R alr().
static arma::mat alr_cpp(const arma::mat& P) {
    const int N = static_cast<int>(P.n_rows);
    const int J = static_cast<int>(P.n_cols);
    arma::mat W(N, J - 1);
    for (int i = 0; i < N; ++i) {
        const double ref = P(i, J - 1);
        for (int j = 0; j < J - 1; ++j)
            W(i, j) = std::log(P(i, j) / ref);
    }
    return W;
}

// ── exported function ────────────────────────────────────────────────────────

//' mdres_core_cpp: C++ core for Mahalanobis distance residuals
//'
//' Implements the full MDres pipeline in C++, replacing four R bottlenecks:
//' \itemize{
//'   \item apply over Mahalanobis distances → row-wise Eigen expression template
//'   \item sweep centering → Eigen .rowwise() broadcast
//'   \item order/sort → std::sort + std::lower_bound / std::upper_bound
//'   \item solve / chol2inv(chol()) → Eigen LLT on small-d covariance
//' }
//' arma::cov() computes the per-observation sample covariance over P posterior
//' samples. A diagonal jitter of 1e-8 is always applied before inversion for
//' numerical stability.
//'
//' @param Y           Numeric matrix (N x J) of observed counts.
//' @param Y_pred_list List of P numeric matrices (each N x J) of posterior
//'                    predictive counts.
//' @return Numeric vector of length N containing z-residuals (class set by
//'         the R wrapper \code{MDres}).
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector mdres_core_cpp(
        const arma::mat&  Y,
        const Rcpp::List& Y_pred_list)
{
    const int N = static_cast<int>(Y.n_rows);
    const int J = static_cast<int>(Y.n_cols);
    const int d = J - 1;
    const int P = static_cast<int>(Y_pred_list.size());

    // ── Step 1: compress and ALR-transform observed counts ──────────────────
    const arma::mat alr_obs = alr_cpp(compress_counts_cpp(Y));  // N × d

    // ── Step 2: build pred_cube[obs, dim, sample] ───────────────────────────
    // pred_cube(i, j, k) = ALR-transformed value for observation i,
    //                       dimension j, posterior sample k.
    arma::cube pred_cube(N, d, P);
    for (int k = 0; k < P; ++k) {
        const arma::mat Yk      = Rcpp::as<arma::mat>(Y_pred_list[k]);
        pred_cube.slice(k) = alr_cpp(compress_counts_cpp(Yk));   // N × d
    }

    // ── Step 3: per-observation Mahalanobis distances + ECDF quantile ───────
    int n_singular = 0;
    Rcpp::NumericVector z_resids(N);

    for (int i = 0; i < N; ++i) {

        // Collect P ALR-transformed predictive samples → (P × d) arma matrix
        arma::mat samp(P, d);
        for (int k = 0; k < P; ++k)
            samp.row(k) = pred_cube.slice(k).row(i);

        // Per-observation sample covariance via arma::cov (divides by P−1,
        // matching R cov()).
        arma::mat S = arma::cov(samp);  // d × d

        // Mirror R singularity check: count obs where min eigenvalue < 1e-8
        const arma::vec ev = arma::eig_sym(S);
        if (ev.min() < 1e-8) ++n_singular;

        // Diagonal jitter for numerical stability before LLT
        S.diag() += 1e-8;

        // Map S into Eigen (arma and Eigen both use column-major; zero copy)
        const Eigen::MatrixXd S_eig =
            Eigen::Map<const Eigen::MatrixXd>(S.memptr(), d, d);

        // LLT decomposition: S = L L^T  (replaces solve / chol2inv(chol()))
        Eigen::LLT<Eigen::MatrixXd> llt(S_eig);
        if (llt.info() != Eigen::Success) {
            Rcpp::warning(
                "LLT decomposition failed for observation %d; "
                "result set to NA.", i + 1);
            z_resids[i] = NA_REAL;
            continue;
        }

        // Column means of predictive samples (replaces apply(pred_array,1:2,mean))
        const arma::rowvec mu = arma::mean(samp, 0);   // 1 × d

        // Map samp and mu into Eigen views (zero copy)
        const Eigen::MatrixXd samp_eig =
            Eigen::Map<const Eigen::MatrixXd>(samp.memptr(), P, d);
        const Eigen::RowVectorXd mu_eig =
            Eigen::Map<const Eigen::RowVectorXd>(mu.memptr(), d);

        // Build centered (P+1) × d matrix:
        //   row 0 = alr_obs[i,] − mu
        //   rows 1..P = samp_eig.rowwise() − mu  (Eigen expression template,
        //               avoids intermediate allocation; replaces R sweep())
        Eigen::MatrixXd centered(P + 1, d);
        for (int j = 0; j < d; ++j)
            centered(0, j) = alr_obs(i, j) - mu(j);
        centered.bottomRows(P).noalias() = samp_eig.rowwise() - mu_eig;

        // Row-wise Mahalanobis distances via a single triangular solve:
        //   md_j = ||L^{-1} centered_j||²
        // Solve L * Z = centered^T  →  Z is d × (P+1), no intermediate alloc.
        const Eigen::MatrixXd Z =
            llt.matrixL().solve(centered.transpose());
        const Eigen::RowVectorXd mds = Z.colwise().squaredNorm(); // 1 × (P+1)

        const double obs_md = mds(0);

        // Build sorted predictive distances for ECDF
        // (replaces R order()/sort() — uses std::sort on a local vector)
        std::vector<double> sorted_post(P);
        for (int k = 0; k < P; ++k) sorted_post[k] = mds(k + 1);
        std::sort(sorted_post.begin(), sorted_post.end());

        const double min_val = sorted_post.front();
        const double max_val = sorted_post.back();

        // ecdf(x) = #{sorted_post ≤ x} / P  (O(log P) via upper_bound)
        auto ecdf_val = [&](double x) -> double {
            return static_cast<double>(
                std::upper_bound(sorted_post.begin(), sorted_post.end(), x)
                - sorted_post.begin()) / P;
        };

        // Randomised-quantile CDF step (matches R u_resids logic exactly)
        double minpct, maxpct;
        if (obs_md <= min_val) {
            minpct = 0.0;
            maxpct = ecdf_val(min_val);
            if (maxpct == 0.0) maxpct = static_cast<double>(P - 1) / P;

        } else if (obs_md > max_val) {
            minpct = ecdf_val(max_val);   // == 1.0
            maxpct = 1.0;
            if (minpct == 1.0) minpct = static_cast<double>(P - 1) / P;

        } else {
            // lower_idx: last 0-based index where sorted_post[k] < obs_md
            const int lower_idx = static_cast<int>(
                std::lower_bound(sorted_post.begin(), sorted_post.end(), obs_md)
                - sorted_post.begin()) - 1;
            minpct = ecdf_val(sorted_post[lower_idx]);
            maxpct = ecdf_val(obs_md);
            if (minpct == 1.0) minpct = static_cast<double>(P - 1) / P;
            if (maxpct == 0.0) maxpct = static_cast<double>(P - 1) / P;
        }

        // Draw uniform CDF quantile, then invert to normal scale
        const double u_i = R::runif(minpct, maxpct);
        z_resids[i] = R::qnorm(u_i, 0.0, 1.0, 1, 0);
    }

    if (n_singular > 0) {
        Rcpp::warning(
            "Covariance matrix near-singular (min eigenvalue < 1e-8) for "
            "%d observation(s); 1e-8 diagonal jitter was applied.",
            n_singular);
    }

    return z_resids;
}

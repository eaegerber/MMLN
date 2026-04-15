#' compress_counts: Smithson-Verkuilen Correction
#'
#' This function implements the correction from Smithson and Verkuilen (2006) to handle zero counts. Input is a vector or matrix of multinomial counts to be compressed row-wise
#' @param y A vector or matrix (rows) of multinomial counts, to be compressed in the presence of zero counts
#'
#' @return A vector or matrix (rows) of compressed multinomial counts
#'
#' @examples
#' \dontrun{
#' # Example Usage:
#' # in the presence of zero counts
#' compress_counts(c(1,0,3))
#' # if no zero counts
#' compress_counts(c(1,2,3))
#' }
#'
#' @export
#'
compress_counts <- function(y) {
  if (is.matrix(y)) {
    zero_rows <- apply(y == 0, 1, any)
    rs <- rowSums(y[zero_rows, , drop=FALSE])
    K <- ncol(y)
    y_new <- y
    y_new[zero_rows, ] <- ((y[zero_rows, , drop=FALSE] * (rs - 1)) + 1/K) / rs
    return(y_new)
  } else {
    if (any(y == 0)) {
      N <- sum(y); K <- length(y)
      return((y * (N - 1) + 1/K) / N)
    } else {
      return(y)
    }
  }
}

#' alr: Additive Log Ratio calculation
#'
#' Estimating underlying log-odds given probability vector (last category is treated as baseline)
#'
#' @param P A vector or matrix (rows) of probabilities
#'
#' @return A vector or matrix (rows) of log-odds
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' alr(c(.1,.3,.6))
#' }
#'
#' @export
alr <- function(P) {
  if(!is.matrix(P)) P <- matrix(P, nrow=1)
  log(P[,-ncol(P),drop=FALSE] / P[,ncol(P)])
}

#' alr_inv: Inverse Additive Log Ratio
#'
#' Transform log-odds back to a probability vector or matrix using inverse additive log ratio, treating the last category as baseline.
#'
#' @param W A vector or matrix (rows) of log-odds values.
#'
#' @return A vector or matrix (rows) of probabilities summing to one.
#'
#' @examples
#' \dontrun{
#' alr_inv(c(0.2, -0.5))
#' }
#'
#' @export
alr_inv <- function(W) {
  if(!is.matrix(W)) W <- matrix(W, nrow=1)
  expW <- exp(W)
  den <- rowSums(expW) + 1
  cbind(expW/den, 1/den)
}

#' dmnl_loglik: Multinomial Logistic-Normal Log Likelihood
#'
#' Calculate the log likelihood for the Multinomial Logistic-Normal distribution given latent W and counts Y.
#'
#' @param W A vector or matrix (rows) of the latent log-odds for observations (J-1 columns).
#' @param Y A vector or matrix (rows) of multinomial counts (J columns).
#'
#' @return The log likelihood value (numeric).
#'
#' @export
dmnl_loglik <- function(W, Y) {
  if(!is.matrix(W)) W <- matrix(W, nrow=1)
  if(!is.matrix(Y)) Y <- matrix(Y, nrow=nrow(Y), byrow=TRUE)
  P <- alr_inv(W)
  sum(Y * log(P))
}

# Metropolis-Hastings proposal helpers
#' pstartoy: Logit Transformation
#'
#' Compute log-odds (logit) from probability values.
#'
#' @param pstarvec A numeric vector of probabilities (between 0 and 1).
#'
#' @return A numeric vector of log-odds values.
#'
#' @examples
#' \dontrun{
#' pstartoy(c(0.2, 0.5, 0.8))
#' }
#'
#' @export
pstartoy <- function(pstarvec) {
  log(pstarvec / (1 - pstarvec))
}

#' ytopstar: Inverse Logit Transformation
#'
#' Compute probabilities from log-odds values.
#'
#' @param yvec A numeric vector of log-odds values.
#'
#' @return A numeric vector of probabilities (between 0 and 1).
#'
#' @examples
#' \dontrun{
#' ytopstar(c(-1, 0, 1))
#' }
#'
#' @export
ytopstar <- function(yvec) {
  exp(yvec) / (1 + exp(yvec))
}

#' betapropdist: Beta Proposal Distribution for MCMC
#'
#' Draw proposal samples for a multinomial logistic-normal model using a Beta distribution.
#'
#' @param YMu_vec A numeric vector combining count vector y and latent means mu (length k+1 + k).
#' @param Sigma A covariance matrix of dimension (k+1) x (k+1).
#'
#' @return A numeric vector of length k containing Beta proposal samples.
#'
#' @examples
#' \dontrun{
#' betapropdist(c(5, 3, 0.1, 0.2), diag(3))
#' }
#'
#' @export
betapropdist <- function(YMu_vec, Sigma) {
  k <- ncol(Sigma)
  y_vec <- YMu_vec[1:(k+1)]
  mu_vec <- YMu_vec[(k+2):length(YMu_vec)]
  exp_mu <- exp(mu_vec)
  denom <- diag(Sigma)
  alpha <- pmax((1 + exp_mu) / denom - exp_mu / (1 + exp_mu), 1e-3)
  beta <- alpha * exp(-mu_vec)
  alpha_star <- y_vec[1:k] + alpha
  beta_star <- y_vec[k+1] + beta
  rbeta(k, alpha_star, beta_star)
}

#' betaloglike: Beta Log-Likelihood with Jacobian Correction
#'
#' Compute the log-likelihood of Beta proposal samples including the Jacobian adjustment term.
#'
#' @param YMuPstar_vec A numeric vector combining counts y, means mu, and proposal pstar values.
#' @param Sigma A covariance matrix.
#'
#' @return A numeric scalar of the log-likelihood.
#'
#' @examples
#' \dontrun{
#' betaloglike(c(5, 3, 0.1, 0.2, 0.3), diag(3))
#' }
#'
#' @export
betaloglike <- function(YMuPstar_vec, Sigma) {
  k <- ncol(Sigma)
  y_vec <- YMuPstar_vec[1:(k+1)]
  mu_vec <- YMuPstar_vec[(k+2):(2*k+1)]
  pstar_vec <- YMuPstar_vec[(2*k+2):length(YMuPstar_vec)]
  exp_mu <- exp(mu_vec)
  denom <- diag(Sigma)
  alpha <- pmax((1 + exp_mu) / denom - exp_mu / (1 + exp_mu), 1e-3)
  beta <- alpha * exp(-mu_vec)
  alpha_star <- y_vec[1:k] + alpha
  beta_star <- y_vec[k+1] + beta
  loglike <- dbeta(pstar_vec, alpha_star, beta_star, log = TRUE)
  logjac <- log(pstar_vec) + log(1 - pstar_vec)
  sum(loglike + logjac)
}

#' normbetapropdist: Normal Approximation to Beta Proposal Distribution
#'
#' Draw proposal samples for the Beta proposal distribution using a Normal approximation via digamma and trigamma.
#'
#' @param YMu_vec A numeric vector combining count vector y and latent means mu.
#' @param Sigma A covariance matrix.
#'
#' @return A numeric vector of length k containing normally-approximated proposal samples.
#'
#' @examples
#' \dontrun{
#' normbetapropdist(c(5, 3, 0.1, 0.2), diag(3))
#' }
#'
#' @export
normbetapropdist <- function(YMu_vec, Sigma) {
  k <- ncol(Sigma)
  y_vec <- YMu_vec[1:(k+1)]
  mu_vec <- YMu_vec[(k+2):length(YMu_vec)]
  result <- numeric(length = k)
  for (i in 1:k) {
    alpha <- ((1 + exp(mu_vec[i])) / Sigma[i, i]) - (1 / (1 + exp(-mu_vec[i])))
    beta <- alpha * exp(-mu_vec[i])
    alpha_star <- y_vec[i] + alpha
    beta_star <- y_vec[k + 1] + beta
    muprop <- digamma(alpha_star) - digamma(beta_star)
    sigprop <- sqrt(trigamma(alpha_star) + trigamma(beta_star))
    result[i] <- rnorm(1, muprop, sigprop)
  }
  result
}

#' normbetaloglike: Normal Approximation to Beta Log-Likelihood
#'
#' Compute the log-likelihood under the Normal approximation of the Beta proposal distribution.
#'
#' @param YMuW_vec A numeric vector combining counts y, means mu, and log-ratio proposals w.
#' @param Sigma A covariance matrix.
#'
#' @return A numeric scalar of the log-likelihood.
#'
#' @examples
#' \dontrun{
#' normbetaloglike(c(5, 3, 0.1, 0.2, 0.3), diag(3))
#' }
#'
#' @export
normbetaloglike <- function(YMuW_vec, Sigma) {
  k <- ncol(Sigma)
  y_vec <- YMuW_vec[1:(k+1)]
  mu_vec <- YMuW_vec[(k+2):(2*k+1)]
  w_vec <- YMuW_vec[(2*k+2):length(YMuW_vec)]
  result <- 0
  for (i in 1:k) {
    alpha <- ((1 + exp(mu_vec[i])) / Sigma[i, i]) - (1 / (1 + exp(-mu_vec[i])))
    beta <- alpha * exp(-mu_vec[i])
    alpha_star <- y_vec[i] + alpha
    beta_star <- y_vec[k + 1] + beta
    muprop <- digamma(alpha_star) - digamma(beta_star)
    sigprop <- sqrt(trigamma(alpha_star) + trigamma(beta_star))
    result <- result + dnorm(w_vec[i], muprop, sigprop, log = TRUE)
  }
  result
}


#' simulate_mixed_mln_data: Simulate mixed-effects MLN data with random intercepts
#'
#' Generate a dataset from a multinomial logistic-normal (MLN) mixed-effects model with one random intercept per group.
#'
#' @param m Integer. Number of groups.
#' @param n_i Integer scalar or vector. Number of observations per group; if scalar, recycled across groups.
#' @param p Integer. Number of fixed-effect covariates (including intercept).
#' @param d Integer. Number of non-baseline logits (total categories = J = d + 1).
#' @param beta Numeric matrix (p × d). Fixed-effect coefficient matrix.
#' @param Sigma Numeric matrix (d × d). Covariance matrix for within-group latent errors.
#' @param Phi Numeric matrix (d × d). Covariance matrix for group-level random intercepts.
#' @param n_mean Numeric. Mean of Poisson distribution for total counts per observation (default 200).
#'
#' @return A named list containing:
#' \describe{
#'   \item{Y}{N × J (where J = d+1) matrix of simulated counts.}
#'   \item{X}{N × p design matrix for fixed effects.}
#'   \item{Z}{N × m design matrix for random intercepts.}
#'   \item{id}{Integer vector of group identifiers (length N).}
#'   \item{n_mean}{Integer vector of total counts per observation (length N).}
#'   \item{beta}{Input fixed-effect coefficient matrix.}
#'   \item{Sigma}{Input within-group error covariance matrix.}
#'   \item{Phi}{Input random-intercept covariance matrix.}
#' }
#'
#' @examples
#' \dontrun{
#' sim_data <- simulate_mixed_mln_data(
#'   m       = 5,
#'   n_i     = 10,
#'   p       = 3,
#'   d       = 2,
#'   beta    = matrix(c(0.5, -1, 0.2, 0.3, 0.7, -0.4), nrow = 3, ncol = 2),
#'   Sigma   = diag(2),
#'   Phi     = diag(c(0.5, 0.5)),
#'   n_mean  = 150
#' )
#' str(sim_data)
#' }
#'
#' @export
simulate_mixed_mln_data <- function(
    m,
    n_i,
    p,
    d,
    beta,
    Sigma,
    Phi,
    n_mean = 200
) {
  if (length(n_i) == 1) {
    n_i <- rep(n_i, m)
  }
  id <- rep(seq_len(m), times = n_i)
  N <- length(id)
  X <- cbind(1, matrix(rnorm(N * (p - 1)), nrow = N, ncol = p - 1))
  Z <- model.matrix(~ factor(id) - 1)
  psi_mat <- mvnfast::rmvn(m, mu = rep(0, d), sigma = Phi)
  n <- round(rpois(N, n_mean))
  Y <- matrix(NA, nrow = N, ncol = d + 1)
  for (j in seq_len(N)) {
    eps <- mvnfast::rmvn(1, mu = rep(0, d), sigma = Sigma)
    w_j <- X[j, , drop = FALSE] %*% beta + eps + psi_mat[id[j], ]
    exp_w <- exp(w_j)
    probs <- c(exp_w / (1 + sum(exp_w)), 1 / (1 + sum(exp_w)))
    Y[j, ] <- as.vector(rmultinom(1, size = n[j], prob = probs))
  }
  list(
    Y = Y,
    X = X,
    Z = Z,
    id = id,
    n = n,
    beta = beta,
    Sigma = Sigma,
    Phi = Phi
  )
}

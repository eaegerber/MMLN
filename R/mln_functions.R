#' FMLN: Fixed-effects Multinomial Logistic-Normal Gibbs Sampler
#'
#' Estimate fixed-effects multinomial logistic-normal models via a Gibbs sampler with Metropolis–Hastings updates for latent variables.
#'
#' @param W Numeric matrix (N × (d+1)) of observed multinomial counts.
#' @param X Numeric matrix (N × p) of fixed-effects covariates.
#' @param n_iter Integer. Total number of MCMC iterations (default 1000).
#' @param burn_in Integer. Number of initial iterations to discard (default 0).
#' @param thin Integer. Thinning interval for saving samples (default 1).
#' @param mh_scale Numeric. Scaling factor for Metropolis–Hastings proposal covariance (default 1).
#' @param prior_settings List. Prior settings:
#' \describe{
#'   \item{beta_var}{Prior variance for β coefficients.}
#'   \item{nu_S}{Degrees of freedom for Sigma prior.}
#'   \item{Lambda_S}{Scale matrix for Sigma prior.}
#' }
#' @param verbose Logical. Print progress updates (default TRUE).
#' @param proposal Character. MH proposal type: one of "norm", "beta", or "normbeta" (default "normbeta").
#'
#' @return A list with components:
#' \describe{
#'   \item{beta_chain}{List of saved β matrices (p × d) across MCMC samples.}
#'   \item{sigma_chain}{List of saved Sigma matrices (d × d).}
#'   \item{y_chain}{List of latent Y matrices (N × d).}
#' }
#'
#' @examples
#' \dontrun{
#' # Fit fixed-effects MLN model
#' res <- FMLN(
#'   W            = count_matrix,
#'   X            = design_matrix,
#'   n_iter       = 2000,
#'   burn_in      = 500,
#'   thin         = 2,
#'   mh_scale     = 1,
#'   prior_settings = list(beta_var = 5, nu_S = d+1, Lambda_S = diag(d)),
#'   proposal     = "normbeta",
#'   verbose      = FALSE
#' )
#' str(res)
#' }
#'
#' @export
FMLN <- function(W, X, n_iter = 1000, burn_in = 0, thin = 1, mh_scale = 1, prior_settings = NULL, verbose = TRUE, proposal = "normbeta") {
  match.arg(proposal, c("norm", "beta", "normbeta"))
  N <- nrow(W)
  d <- ncol(W) - 1
  p <- ncol(X)
  PA <- rowSums(W)

  if (is.null(prior_settings)) {
    prior_settings <- list(
      beta_var = 10,
      nu_S = d + 1,
      Lambda_S = diag(d)
    )
  }

  beta <- matrix(0, nrow = p, ncol = d)
  Sigma <- diag(d)

  keep_iters <- seq(burn_in + 1, n_iter, by = thin)
  n_save <- length(keep_iters)

  beta_chain <- vector("list", n_save)
  sigma_chain <- vector("list", n_save)
  y_chain <- vector("list", n_save)

  warned_na_ratio <- FALSE
  if (verbose) pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  start_time <- Sys.time()
  save_idx <- 1

  Y <- alr(compress_counts(W))
  Sigma_inv <- chol2inv(chol(Sigma))
  S_xx_inv <- chol2inv(chol(crossprod(X)))

  for (i in seq_len(n_iter)) {
    Mu <- tcrossprod(X, t(beta))
    wmu <- cbind(W, Mu)

    if (proposal == "norm") {
      Y_new <- Y + rmvn(N, mu = rep(0, d), sigma = mh_scale * Sigma)
      log_q_old <- rep(0, N)
      log_q_prop <- rep(0, N)
    } else if (proposal == "beta") {
      Pstar <- t(apply(Y, 1, ytopstar))
      Pstar_new <- t(apply(wmu, 1, betapropdist, Sigma = mh_scale * Sigma))
      Y_new <- t(apply(Pstar_new, 1, pstartoy))
      wmuPstar <- cbind(wmu, Pstar)
      wmuPstar_new <- cbind(wmu, Pstar_new)
      log_q_old <- apply(wmuPstar, 1, betaloglike, Sigma = mh_scale * Sigma)
      log_q_prop <- apply(wmuPstar_new, 1, betaloglike, Sigma = mh_scale * Sigma)
    } else if (proposal == "normbeta") {
      Y_new <- t(apply(wmu, 1, normbetapropdist, Sigma = mh_scale * Sigma))
      wmuy <- cbind(wmu, Y)
      wmuy_new <- cbind(wmu, Y_new)
      log_q_old <- apply(wmuy, 1, normbetaloglike, Sigma = mh_scale * Sigma)
      log_q_prop <- apply(wmuy_new, 1, normbetaloglike, Sigma = mh_scale * Sigma)
    }

    Y_diff <- Y - Mu
    Y_new_diff <- Y_new - Mu
    newnormpart <- rowSums(tcrossprod(Y_new_diff, t(Sigma_inv)) * Y_new_diff) / 2
    oldnormpart <- rowSums(tcrossprod(Y_diff, t(Sigma_inv)) * Y_diff) / 2
    sexpY <- rowSums(exp(Y))
    sexpYn <- rowSums(exp(Y_new))
    newloglike <- rowSums(W[, 1:d] * Y_new[, 1:d]) - rowSums(W * log1p(sexpYn)) - newnormpart
    oldloglike <- rowSums(W[, 1:d] * Y[, 1:d]) - rowSums(W * log1p(sexpY)) - oldnormpart

    ratio <- newloglike - oldloglike + log_q_old - log_q_prop
    if(!warned_na_ratio && anyNA(ratio)) {
      warning("NA detected in MH acceptance ratio; these proposals will be rejected.")
      warned_na_ratio <- TRUE
    }
    ratio[is.na(ratio)] <- -Inf

    accept <- log(runif(N)) < ratio
    Y[accept, ] <- Y_new[accept, ]

    R <- Y
    beta_hat <- tcrossprod(S_xx_inv, t(crossprod(X, R)))
    post_beta_vec_cov <- kronecker(S_xx_inv, Sigma)
    post_beta_vec <- rmvn(1, as.vector(beta_hat), post_beta_vec_cov)
    beta <- matrix(post_beta_vec, nrow = p)

    Sigma <- chol2inv(chol(rWishart(1, df = prior_settings$nu_S + N,
                                    Sigma = solve(prior_settings$Lambda_S + crossprod(Y - X %*% beta)))[,,1]))
    Sigma_inv <- chol2inv(chol(Sigma))

    if (i %in% keep_iters) {
      beta_chain[[save_idx]] <- beta
      sigma_chain[[save_idx]] <- Sigma
      y_chain[[save_idx]] <- Y
      save_idx <- save_idx + 1
    }

    if (verbose && (i %% max(1, floor(n_iter / 100)) == 0 || i == n_iter)) {
      setTxtProgressBar(pb, i)
      elapsed <- Sys.time() - start_time
      eta <- (as.numeric(elapsed) / i) * (n_iter - i)
      cat(sprintf("\r ETA: %s", format(.POSIXct(eta, tz = "GMT"), "%M:%S")))
      flush.console()
    }
  }
  if (verbose) close(pb)

  list(
    beta_chain = beta_chain,
    sigma_chain = sigma_chain,
    y_chain = y_chain
  )
}



#' MMLN: Mixed-effects Multinomial Logistic-Normal Gibbs Sampler
#'
#' Estimate mixed-effects multinomial logistic-normal models with group-level random intercepts.
#'
#' @param W Numeric matrix (N × (d+1)) of observed counts.
#' @param X Numeric matrix (N × p) of fixed-effects covariates.
#' @param Z Numeric matrix (N × m) of random-effects design (group indicators).
#' @param n_iter Integer. Total MCMC iterations (default 1000).
#' @param burn_in Integer. Burn-in iterations to discard (default 0).
#' @param thin Integer. Thinning interval (default 1).
#' @param mh_scale Numeric. MH proposal scale factor (default 1).
#' @param prior_settings List. Prior settings:
#' \describe{
#'   \item{beta_var}{Prior variance for β.}
#'   \item{nu_S}{Degrees of freedom for Sigma prior.}
#'   \item{Lambda_S}{Scale for Sigma prior.}
#'   \item{nu_P}{Degrees of freedom for Phi prior.}
#'   \item{Lambda_P}{Scale for Phi prior.}
#' }
#' @param verbose Logical. Print progress bar (default TRUE).
#' @param proposal Character. One of "norm", "beta", or "normbeta" (default "normbeta").
#'
#' @return A list with:
#' \describe{
#'   \item{beta_chain}{List of saved fixed-effect β matrices (p × d).}
#'   \item{sigma_chain}{List of saved Sigma matrices (d × d).}
#'   \item{phi_chain}{List of saved Phi matrices (d × d) for random intercepts.}
#'   \item{psi_chain}{List of saved random-intercept matrices (m × d).}
#'   \item{y_chain}{List of latent Y matrices (N × d).}
#' }
#'
#' @examples
#' \dontrun{
#' res_mixed <- MMLN(
#'   W              = count_matrix,
#'   X              = fixed_design,
#'   Z              = random_design,
#'   n_iter         = 1500,
#'   burn_in        = 300,
#'   thin           = 5,
#'   prior_settings = list(beta_var=5, nu_S=d+1, Lambda_S=I, nu_P=d+1, Lambda_P=I),
#'   proposal       = "normbeta",
#'   verbose        = FALSE
#' )
#' str(res_mixed)
#' }
#'
#' @export
MMLN <- function(W, X, Z, n_iter = 1000, burn_in = 0, thin = 1, mh_scale = 1, prior_settings = NULL, verbose = TRUE, proposal = "normbeta") {
  match.arg(proposal, c("norm", "beta", "normbeta"))
  N <- nrow(W); d <- ncol(W) - 1; p <- ncol(X); m <- ncol(Z)
  if(is.null(prior_settings)) prior_settings <- list(
    beta_var = 10,
    nu_S     = d + 1, Lambda_S = diag(d),
    nu_P     = d + 1, Lambda_P = diag(d)
  )
  keep <- seq(burn_in + 1, n_iter, by = thin)
  n_save <- length(keep)
  beta_chain  <- vector("list", n_save)
  sigma_chain <- vector("list", n_save)
  phi_chain   <- vector("list", n_save)
  psi_chain   <- vector("list", n_save)
  y_chain     <- vector("list", n_save)

  # initialize
  Y         <- alr(compress_counts(W))   # from mln_helpers.R
  beta      <- matrix(0, p, d)
  Sigma     <- diag(d)
  psi       <- matrix(0, m, d)
  Phi       <- diag(d)
  Sigma_inv <- chol2inv(chol(Sigma))
  S_xx_inv  <- chol2inv(chol(crossprod(X)))

  warned_na_ratio <- FALSE
  if(verbose) {
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
    start_time <- Sys.time()
  }
  save_i <- 1

  for(it in seq_len(n_iter)) {
    # current mean
    Mu  <- X %*% beta + Z %*% psi   # N x d
    wmu <- cbind(W, Mu)

    # MH update of latent Y
    if(proposal == "norm") {
      Y_prop    <- Y + mvnfast::rmvn(N, mu = rep(0, d), sigma = mh_scale * Sigma)
      log_q_old <- log_q_new <- rep(0, N)
    } else if(proposal == "beta") {
      P_old     <- t(apply(Y, 1, ytopstar))
      P_new     <- t(apply(wmu, 1, betapropdist, Sigma = mh_scale * Sigma))
      Y_prop    <- t(apply(P_new, 1, pstartoy))
      log_q_old <- apply(cbind(wmu, P_old), 1, betaloglike, Sigma = mh_scale * Sigma)
      log_q_new <- apply(cbind(wmu, P_new), 1, betaloglike, Sigma = mh_scale * Sigma)
    } else {
      Y_prop    <- t(apply(wmu, 1, normbetapropdist, Sigma = mh_scale * Sigma))
      log_q_old <- apply(cbind(wmu, Y), 1, normbetaloglike, Sigma = mh_scale * Sigma)
      log_q_new <- apply(cbind(wmu, Y_prop), 1, normbetaloglike, Sigma = mh_scale * Sigma)
    }
    expY   <- rowSums(exp(Y)); expYp <- rowSums(exp(Y_prop))
    ll_old <- rowSums(W[,1:d] * Y[,1:d]) - rowSums(W * log1p(expY)) -
      0.5 * rowSums((Y - Mu) %*% Sigma_inv * (Y - Mu))
    ll_new <- rowSums(W[,1:d] * Y_prop[,1:d]) - rowSums(W * log1p(expYp)) -
      0.5 * rowSums((Y_prop - Mu) %*% Sigma_inv * (Y_prop - Mu))
    ratio  <- ll_new - ll_old + log_q_new - log_q_old
    if(!warned_na_ratio && anyNA(ratio)) {
      warning("NA detected in MH acceptance ratio; these proposals will be rejected.")
      warned_na_ratio <- TRUE
    }
    ratio[is.na(ratio)] <- -Inf

    accepted <- log(runif(N)) < ratio
    Y[accepted, ] <- Y_prop[accepted, ]

    # update random intercepts psi_j
    R_tot <- Y - X %*% beta
    for(j in seq_len(m)) {
      idx <- which(Z[, j] == 1)
      R_j <- R_tot[idx, , drop = FALSE]
      V_j <- solve(chol2inv(chol(Phi)) + length(idx) * Sigma_inv)
      M_j <- V_j %*% (Sigma_inv %*% colSums(R_j))
      psi[j, ] <- mvnfast::rmvn(1, mu = as.vector(M_j), sigma = V_j)
    }

    # update Phi
    S_psi <- t(psi) %*% psi
    Phi   <- solve(rWishart(1,
                            df    = prior_settings$nu_P + m,
                            Sigma = solve(prior_settings$Lambda_P + S_psi))[,,1])

    # update beta
    R     <- Y - Z %*% psi
    beta0 <- S_xx_inv %*% (t(X) %*% R)
    cov_b <- kronecker(S_xx_inv, Sigma)
    beta  <- matrix(mvnfast::rmvn(1,
                                  mu    = as.vector(beta0),
                                  sigma = cov_b),
                    nrow = p)

    # update Sigma
    Eps   <- R - X %*% beta
    S_mat <- t(Eps) %*% Eps
    Sigma <- solve(rWishart(1,
                            df    = prior_settings$nu_S + N,
                            Sigma = solve(prior_settings$Lambda_S + S_mat))[,,1])
    Sigma_inv <- chol2inv(chol(Sigma))

    # save samples
    if(it %in% keep) {
      beta_chain[[save_i]]  <- beta
      sigma_chain[[save_i]] <- Sigma
      phi_chain[[save_i]]   <- Phi
      psi_chain[[save_i]]   <- psi
      y_chain[[save_i]]     <- Y
      save_i <- save_i + 1
    }

    if (verbose && (it %% max(1, floor(n_iter / 100)) == 0 || it == n_iter)) {
      setTxtProgressBar(pb, it)
      elapsed <- Sys.time() - start_time
      eta <- (as.numeric(elapsed) / it) * (n_iter - it)
      cat(sprintf("\r ETA: %s", format(.POSIXct(eta, tz = "GMT"), "%M:%S")))
      flush.console()
    }
  }

  if(verbose) close(pb)

  list(
    beta_chain  = beta_chain,
    sigma_chain = sigma_chain,
    phi_chain   = phi_chain,
    psi_chain   = psi_chain,
    y_chain     = y_chain
  )
}

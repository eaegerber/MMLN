#' FMLN: Fixed-effects Multinomial Logistic-Normal Gibbs Sampler
#'
#' Estimate fixed-effects multinomial logistic-normal models via a Gibbs sampler with Metropolis–Hastings updates for latent variables.
#'
#' @param Y Numeric matrix (N × J) of observed multinomial counts.
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
#'   \item{beta_chain}{List of saved β matrices (p × (J-1)) across MCMC samples.}
#'   \item{sigma_chain}{List of saved Sigma matrices ((J-1) × (J-1)).}
#'   \item{w_chain}{List of latent W matrices (N × (J-1)).}
#'   \item{mhaccept_chain}{List of Metropolis-Hastings acceptance ratios (length = n_iter).}
#' }
#'
#' @examples
#' \dontrun{
#' # Fit fixed-effects MLN model
#' res <- FMLN(
#'   Y            = count_matrix,
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
FMLN <- function(Y, X, n_iter = 1000, burn_in = 0, thin = 1, mh_scale = 1, prior_settings = NULL, verbose = TRUE, proposal = "normbeta") {
  match.arg(proposal, c("norm", "beta", "normbeta"))
  N <- nrow(Y)
  d <- ncol(Y) - 1
  p <- ncol(X)
  PA <- rowSums(Y)

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
  w_chain <- vector("list", n_save)
  mhaccept_chain <- vector("list", n_iter)

  warned_na_ratio <- FALSE
  if (verbose) pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  start_time <- Sys.time()
  save_idx <- 1

  W <- alr(compress_counts(Y))
  Sigma_inv <- chol2inv(chol(Sigma))
  S_xx_inv <- chol2inv(chol(crossprod(X)))

  for (i in seq_len(n_iter)) {
    Mu <- tcrossprod(X, t(beta))
    ymu <- cbind(Y, Mu)

    if (proposal == "norm") {
      W_new <- W + rmvn(N, mu = rep(0, d), sigma = mh_scale * Sigma)
      log_q_old <- rep(0, N)
      log_q_prop <- rep(0, N)
    } else if (proposal == "beta") {
      mh_update_output <- mh_update_beta_cpp(W, Y, Mu, mh_scale * Sigma)
      W_new      <- mh_update_output$W_new
      log_q_old  <- mh_update_output$log_q_old
      log_q_prop <- mh_update_output$log_q_new
    } else if (proposal == "normbeta") {
      mh_update_output    <- mh_update_normbeta_cpp(W, Y, Mu, mh_scale * Sigma)
      W_new     <- mh_update_output$W_new
      log_q_old <- mh_update_output$log_q_old
      log_q_prop <- mh_update_output$log_q_new
    }

    W_diff <- W - Mu
    W_new_diff <- W_new - Mu
    newnormpart <- rowSums(tcrossprod(W_new_diff, t(Sigma_inv)) * W_new_diff) / 2
    oldnormpart <- rowSums(tcrossprod(W_diff, t(Sigma_inv)) * W_diff) / 2
    sexpW <- rowSums(exp(W))
    sexpWn <- rowSums(exp(W_new))
    newloglike <- rowSums(Y[, 1:d] * W_new[, 1:d]) - rowSums(Y * log1p(sexpWn)) - newnormpart
    oldloglike <- rowSums(Y[, 1:d] * W[, 1:d]) - rowSums(Y * log1p(sexpW)) - oldnormpart

    ratio <- newloglike - oldloglike + log_q_old - log_q_prop
    if(!warned_na_ratio && anyNA(ratio)) {
      warning("NA detected in MH acceptance ratio; these proposals will be rejected.")
      warned_na_ratio <- TRUE
    }
    ratio[is.na(ratio)] <- -Inf

    accept <- log(runif(N)) < ratio

    # track acceptance ratio for tuning
    mhaccept_chain[[i]] <- sum(accept) / N

    W[accept, ] <- W_new[accept, ]
    W[is.na(W)] <- 0

    R <- W
    beta_hat <- tcrossprod(S_xx_inv, t(crossprod(X, R)))
    post_beta_vec_cov <- kronecker(S_xx_inv, Sigma)
    post_beta_vec <- rmvn(1, as.vector(beta_hat), post_beta_vec_cov)
    beta <- matrix(post_beta_vec, nrow = p)
    
    # Running into non positive-definite issues, so let's try jittering
    S <- prior_settings$Lambda_S + crossprod(W - X %*% beta)
    S <- (S + t(S)) / 2
    diag(S) <- diag(S) + 1e-8

    Wishscale <- chol2inv(chol(S))
    Wishscale <- (Wishscale + t(Wishscale)) / 2
    diag(Wishscale) <- diag(Wishscale) + 1e-8

    Sigma <- chol2inv(chol(rWishart(1, df = prior_settings$nu_S + N,
                                    Sigma = Wishscale)[,,1]))
    Sigma_inv <- chol2inv(chol(Sigma))

    if (i %in% keep_iters) {
      beta_chain[[save_idx]] <- beta
      sigma_chain[[save_idx]] <- Sigma
      w_chain[[save_idx]] <- W
      save_idx <- save_idx + 1
    }

    if (verbose && (i %% max(1, floor(n_iter / 100)) == 0 || i == n_iter)) {
      setTxtProgressBar(pb, i)
      elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
      eta_sec     <- elapsed_sec / i * (n_iter - i)

      h <- floor(eta_sec / 3600)
      mn <- floor((eta_sec %% 3600) / 60)
      s <- round(eta_sec %% 60)

      eta_str <- sprintf("%02d:%02d:%02d", h, mn, s)
      cat(sprintf("\r ETA: %s", eta_str))
      flush.console()
    }
  }
  if (verbose) close(pb)

  list(
    beta_chain = beta_chain,
    sigma_chain = sigma_chain,
    w_chain = w_chain,
    mhaccept_chain = mhaccept_chain
  )
}


## To Do: Update below function to allow for any number of random effects
## Currently: Z can only be N x m, which means you can only have random intercepts (one for each of m groups)
## To Update: Z needs to be able to be N x (m x q) where each group has q random covariates
## The math for this is established in the Gerber \& Craig (2021) paper, but still needs to be implemented here

#' MMLN: Mixed-effects Multinomial Logistic-Normal Gibbs Sampler
#'
#' Estimate mixed-effects multinomial logistic-normal models with group-level random intercepts.
#'
#' @param Y Numeric matrix (N × J) of observed counts.
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
#'   \item{beta_chain}{List of saved fixed-effect β matrices (p × (J-1)).}
#'   \item{sigma_chain}{List of saved Sigma matrices ((J-1) × (J-1)).}
#'   \item{phi_chain}{List of saved Phi matrices ((J-1) × (J-1)) for random intercepts.}
#'   \item{psi_chain}{List of saved random-intercept matrices (m × (J-1)).}
#'   \item{w_chain}{List of latent W matrices (N × (J-1)).}
#'   \item{mhaccept_chain}{List of Metropolis-Hastings acceptance ratios (length = n_iter).}
#' }
#'
#' @examples
#' \dontrun{
#' res_mixed <- MMLN(
#'   Y              = count_matrix,
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
MMLN <- function(Y, X, Z, n_iter = 1000, burn_in = 0, thin = 1, mh_scale = 1, prior_settings = NULL, verbose = TRUE, proposal = "normbeta") {
  match.arg(proposal, c("norm", "beta", "normbeta"))
  N <- nrow(Y); d <- ncol(Y) - 1; p <- ncol(X); m <- ncol(Z)
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
  w_chain     <- vector("list", n_save)
  mhaccept_chain <- vector("list", n_iter)

  # initialize
  W         <- alr(compress_counts(Y))   # from mln_helpers.R
  beta      <- matrix(0, p, d)
  Sigma     <- diag(d)
  psi       <- matrix(0, m, d)
  Phi       <- diag(d)
  Sigma_inv <- chol2inv(chol(Sigma))
  S_xx_inv  <- chol2inv(chol(crossprod(X)))

  # pre-compute group membership indices once (Z is constant across iterations)
  group_indices <- lapply(seq_len(m), function(j) which(Z[, j] == 1))
  group_sizes   <- lengths(group_indices)
  unique_sizes  <- unique(group_sizes)


  warned_na_ratio <- FALSE
  if(verbose) {
    pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
    start_time <- Sys.time()
  }
  save_i <- 1

  for(it in seq_len(n_iter)) {
    # current mean
    Mu  <- X %*% beta + Z %*% psi   # N x d
    ymu <- cbind(Y, Mu)

    # MH update of latent W
    if(proposal == "norm") {
      W_prop    <- W + mvnfast::rmvn(N, mu = rep(0, d), sigma = mh_scale * Sigma)
      log_q_old <- log_q_new <- rep(0, N)
    } else if(proposal == "beta") {
      mh_update_output <- mh_update_beta_cpp(W, Y, Mu, mh_scale * Sigma)
      W_prop    <- mh_update_output$W_new
      log_q_old <- mh_update_output$log_q_old
      log_q_new <- mh_update_output$log_q_new
    } else {
      mh_update_output    <- mh_update_normbeta_cpp(W, Y, Mu, mh_scale * Sigma)
      W_prop    <- mh_update_output$W_new
      log_q_old <- mh_update_output$log_q_old
      log_q_new <- mh_update_output$log_q_new
    }
    expW   <- rowSums(exp(W)); expWp <- rowSums(exp(W_prop))
    ll_old <- rowSums(Y[,1:d] * W[,1:d]) - rowSums(Y * log1p(expW)) -
      0.5 * rowSums((W - Mu) %*% Sigma_inv * (W - Mu))
    ll_new <- rowSums(Y[,1:d] * W_prop[,1:d]) - rowSums(Y * log1p(expWp)) -
      0.5 * rowSums((W_prop - Mu) %*% Sigma_inv * (W_prop - Mu))
    ratio  <- ll_new - ll_old + log_q_old - log_q_new
    if(!warned_na_ratio && anyNA(ratio)) {
      warning("NA detected in MH acceptance ratio; these proposals will be rejected.")
      warned_na_ratio <- TRUE
    }
    ratio[is.na(ratio)] <- -Inf

    accepted <- log(runif(N)) < ratio

    # track acceptance ratio for tuning
    mhaccept_chain[[it]] <- sum(accepted) / N

    W[accepted, ] <- W_prop[accepted, ]
    W[is.na(W)] <- 0

    # update random intercepts psi_j
    R_tot   <- W - X %*% beta
    Phi_inv <- chol2inv(chol(Phi))

    # V_j depends only on group size; compute once per unique size (not once per group)
    V_by_size <- setNames(
      lapply(unique_sizes, function(n) chol2inv(chol(Phi_inv + n * Sigma_inv))),
      as.character(unique_sizes)
    )

    # group indices pre-computed outside MCMC loop; M_j uses upstream column-vector
    # formula for bit-exact arithmetic; V_j reused from cache above
    for(j in seq_len(m)) {
      R_j      <- R_tot[group_indices[[j]], , drop = FALSE]
      V_j      <- V_by_size[[as.character(group_sizes[j])]]
      M_j      <- V_j %*% (Sigma_inv %*% colSums(R_j))

      # This line is the slowest part, but the RNG chain get's reset for each j
      # so the result would be different, but both valid randomness
      psi[j, ] <- mvnfast::rmvn(1, mu = as.vector(M_j), sigma = V_j)
    }


    # Running into non positive-definite issues, so let's try jittering    
    S1 <- prior_settings$Lambda_P + S_psi
    S1 <- (S1 + t(S1)) / 2
    diag(S1) <- diag(S1) + 1e-8

    Wish1scale <- chol2inv(chol(S1))
    Wish1scale <- (Wish1scale + t(Wish1scale)) / 2
    diag(Wish1scale) <- diag(Wish1scale) + 1e-8

    Phi   <- chol2inv(chol(rWishart(1,
                            df    = prior_settings$nu_P + m,
                            Sigma = Wish1scale)[,,1]))

    # update beta
    R     <- W - Z %*% psi
    beta0 <- S_xx_inv %*% (t(X) %*% R)
    cov_b <- kronecker(S_xx_inv, Sigma)
    beta  <- matrix(mvnfast::rmvn(1,
                                  mu    = as.vector(beta0),
                                  sigma = cov_b),
                    nrow = p)



    # update Sigma
    Eps   <- R - X %*% beta
    S_mat <- t(Eps) %*% Eps

    # Running into non positive-definite issues, so let's try jittering    
    S2 <- prior_settings$Lambda_S + S_mat
    S2 <- (S2 + t(S2)) / 2
    diag(S2) <- diag(S2) + 1e-8

    Wish2scale <- chol2inv(chol(S2))
    Wish2scale <- (Wish2scale + t(Wish2scale)) / 2
    diag(Wish2scale) <- diag(Wish2scale) + 1e-8

    Sigma <- chol2inv(chol(rWishart(1,
                            df    = prior_settings$nu_S + N,
                            Sigma = Wish2scale)[,,1]))
    Sigma_inv <- chol2inv(chol(Sigma))

    # save samples
    if(it %in% keep) {
      beta_chain[[save_i]]  <- beta
      sigma_chain[[save_i]] <- Sigma
      phi_chain[[save_i]]   <- Phi
      psi_chain[[save_i]]   <- psi
      w_chain[[save_i]]     <- W
      save_i <- save_i + 1
    }

    if (verbose && (it %% max(1, floor(n_iter / 100)) == 0 || it == n_iter)) {
      setTxtProgressBar(pb, it)
      elapsed_sec <- as.numeric(difftime(Sys.time(), start_time, units="secs"))
      eta_sec     <- elapsed_sec / it * (n_iter - it)

      h <- floor(eta_sec / 3600)
      mn <- floor((eta_sec %% 3600) / 60)
      s <- round(eta_sec %% 60)

      eta_str <- sprintf("%02d:%02d:%02d", h, mn, s)
      cat(sprintf("\r ETA: %s", eta_str))
      flush.console()
    }
  }

  if(verbose) close(pb)

  list(
    beta_chain  = beta_chain,
    sigma_chain = sigma_chain,
    phi_chain   = phi_chain,
    psi_chain   = psi_chain,
    w_chain     = w_chain,
    mhaccept_chain = mhaccept_chain
  )
}

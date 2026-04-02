#' plot_trace_and_summary: Traceplots and summary tables for MCMC chains
#'
#' This function plots and provides posterior means and credible intervals for the MCMC samples of fitted MLN models.
#'
#' @param chain_array An array (rows by cols by iterations) obtained via simplify2array(fit$chain) from a model fit.
#' @param param_name Character. Name of the parameter to label the plots (default "param").
#' @param max_rows Integer. Maximum number of rows of parameter indices to plot (default 4).
#'
#' @return A 3D array (cols by stats by rows) of summary statistics: mean, sd, 2.5%, 97.5% quantiles.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' sim <- simulate_mixed_mln_data(
#' m       = 30,
#' n_i     = 10,
#' p       = 2,
#' d       = 2,
#' beta    = matrix(c(-2, 0.5, 1, -0.), nrow = 2, ncol = 2),
#' Sigma   = diag(c(1.2, .7)),
#' Phi     = diag(c(0.5, 0.3)),
#' n_mean = 100)
#'
#' res <- MMLN(
#' W       = sim$W,
#' X       = sim$X,
#' Z       = sim$Z,
#' n_iter  = 1000,
#' burn_in = 180,
#' thin    = 2,
#' proposal= "normbeta",
#' verbose = TRUE)
#'
#' beta_chain_array <- simplify2array(res$beta_chain)
#' plot_trace_and_summary(beta_chain_array, "beta")

#'
#' }
#'
#' @export
plot_trace_and_summary <- function(chain_array, param_name = "param", max_rows = 4) {
  dims <- dim(chain_array)
  par(mfrow = c(2, 2))
  n_rows <- min(dims[1], max_rows)
  for (i in 1:n_rows) {
    for (j in 1:dims[2]) {
      trace <- chain_array[i, j, ]
      plot(trace, type = "l", main = paste0("Trace: ", param_name, "[", i, ",", j, "]"), ylab = "Value", xlab = "Iteration")
      abline(h = mean(trace), col = "red")
    }
  }

  summary_array <- array(NA, dim = c(dim(chain_array)[1], dim(chain_array)[2], 4))
  for (i in seq_len(dim(chain_array)[1])) {
    for (j in seq_len(dim(chain_array)[2])) {
      trace <- chain_array[i, j, ]
      summary_array[i, j, ] <- c(
        mean = mean(trace),
        sd = sd(trace),
        `2.5%` = quantile(trace, 0.025),
        `97.5%` = quantile(trace, 0.975)
      )
    }
  }
  dimnames(summary_array) <- list(
    paste0("row", seq_len(dim(chain_array)[1])),
    paste0("col", seq_len(dim(chain_array)[2])),
    c("mean", "sd", "2.5%", "97.5%")
  )

  par(mfrow = c(1, 1))
  return(aperm(summary_array, c(2, 3, 1)))
}




#' compute_dic: Deviance Information Criterion
#'
#' Computes the DIC (Deviance Information Criterion) and effective number of parameters pD.
#'
#' @param loglik_chain Numeric vector of posterior log-likelihood samples across iterations.
#' @param loglik_hat Numeric. The posterior mean of the log-likelihood (evaluated at posterior means).
#'
#' @return A list with elements:
#' \describe{
#'   \item{DIC}{The Deviance Information Criterion value.}
#'   \item{pD}{The effective number of parameters.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' sim <- simulate_mixed_mln_data(
#' m       = 30,
#' n_i     = 10,
#' p       = 2,
#' d       = 2,
#' beta    = matrix(c(-2, 0.5, 1, -0.), nrow = 2, ncol = 2),
#' Sigma   = diag(c(1.2, .7)),
#' Phi     = diag(c(0.5, 0.3)),
#' n_mean = 100)
#'
#' res <- MMLN(
#' Y       = sim$Y,
#' X       = sim$X,
#' Z       = sim$Z,
#' n_iter  = 1000,
#' burn_in = 180,
#' thin    = 2,
#' proposal= "normbeta",
#' verbose = TRUE)
#'
#' ll_chain <- lapply(seq_len(length(res$w_chain)),
#' function(i){
#' dmnl_loglik(Y = res$y_chain[[i]], Y = sim$Y)})
#'
#' W_hat <- alr(compress_counts(sim$Y)/rowSums(sim$Y))
#' ll_hat <- dmnl_loglik(W = W_hat, Y = sim$Y)
#'
#' compute_dic(ll_chain, ll_hat)
#'
#' }
#'
#' @export
compute_dic <- function(loglik_chain, loglik_hat) {
  pD <- 2 * (loglik_hat - mean(loglik_chain))
  DIC <- -2 * loglik_hat + 2 * pD
  list(DIC = DIC, pD = pD)
}







#' sample_posterior_predictive: Simulate posterior predictive counts
#'
#' Generate posterior predictive multinomial counts from a fitted (mixed/fixed effects) MLN model.
#'
#' @param X Numeric matrix (N by p) of fixed-effects design.
#' @param beta Numeric matrix (p by (J-1)) of fixed-effects coefficients.
#' @param Sigma Numeric ((J-1) by (J-1)) covariance matrix of latent variables.
#' @param n Numeric vector or scalar for total counts per observation (length N or 1).
#' @param Z Optional numeric matrix (N by q) of random-effects design (required if mixed = TRUE).
#' @param psi Optional numeric matrix (q by (J-1)) of random-effects coefficients (required if mixed = TRUE).
#' @param mixed Logical; include random effects when TRUE (default TRUE).
#' @param verbose Logical; include progress bar when TRUE (default TRUE).
#'
#' @return Numeric matrix (N by J) of simulated counts for each category.
#'
#' @examples
#' \dontrun{
#' sim_counts <- sample_posterior_predictive(X, beta, Sigma, n,
#'                                           Z = Z, psi = psi, mixed = TRUE)
#' }
#'
#' @export
sample_posterior_predictive <- function(X, beta, Sigma, n,
                                        Z = NULL, psi = NULL,
                                        mixed = TRUE,
                                        verbose = TRUE) {
  N <- nrow(X)
  d <- ncol(Sigma)
  # linear predictor
  mu <- X %*% beta
  if(mixed) {
    if(is.null(Z) || is.null(psi)) stop("For mixed = TRUE, must supply Z and psi")
    mu <- mu + Z %*% psi
  }
  # simulate latent alr-scale outcomes
  W_hat <- mu + mvnfast::rmvn(N, mu = rep(0, d), sigma = Sigma)
  # transform back to simplex
  P_hat <- alr_inv(W_hat)
  # sample counts per observation
  Y <- matrix(NA, nrow = N, ncol = d + 1)
  if (verbose && interactive()) {
    cat("Generating posterior predictive samples:\n")
    start_time <- Sys.time()
    last_pct   <- -1
  }
  for(i in seq_len(N)) {
    pct <- floor(100 * i / N)
    if (verbose && interactive() && (pct != last_pct || i == N)) {
      last_pct <- pct
      elapsed <- Sys.time() - start_time
      eta <- (as.numeric(elapsed) / i) * (N - i)
      cat(sprintf("\r[%3d%%] ETA: %s", pct, format(.POSIXct(eta, tz="GMT"), "%M:%S")))
      flush.console()
    }
    size_i <- if(length(n) > 1) n[i] else n
    Y[i, ] <- rmultinom(1, size = size_i, prob = P_hat[i, ])
  }
  Y
}





#' MDres: Mahalanobis residuals for predictive checks
#'
#' Compute quantile-normalized Mahalanobis residuals comparing observed counts to posterior predictive samples.
#'
#' @param Y Numeric matrix (N by J) of observed counts.
#' @param Y_pred_list A list of P numeric matrices (each N by J) of posterior predictive counts.
#'
#' @return Numeric vector of length N of residuals (class 'mdres').
#'
#' @examples
#' \dontrun{
#' resids <- MDres(Y_obs, Y_pred_list)
#' summary(resids)
#' }
#'
#' @export
MDres <- function(Y, Y_pred_list) {
  z_resids <- mdres_core_cpp(Y, Y_pred_list)
  class(z_resids) <- "mdres"
  return(z_resids)
}





#' summary.mdres: Summary and diagnostics for Mahalanobis residuals
#'
#' Perform a Kolmogorov-Smirnov test and generate a QQ plot for Mahalanobis residuals.
#'
#' @param object Numeric vector of residuals (class 'mdres').
#' @param ... Additional arguments (ignored).
#'
#' @return Invisibly returns the Kolmogorov-Smirnov test object.
#'
#' @importFrom car qqPlot
#'
#' @export
summary.mdres <- function(object, ...) {
  # `object` is just a numeric vector of residuals
  cat("Kolmogorov-Smirnov test for normality of Mahalanobis residuals:\n")
  ks <- suppressWarnings(
    ks.test(object,
            "pnorm",
            mean = 0, sd = 1,
            exact = FALSE)
  )
  print(ks)

  # QQ plot
  qqPlot(object, main = "QQ Plot of Mahalanobis Residuals", ylab = 'Residual Quantiles', xlab = "Theoretical Quantiles", line = "none", col.lines = "red", envelope = list(style = "lines"))
  abline(0, 1, col='red', lwd = 2)

  invisible(ks)
}

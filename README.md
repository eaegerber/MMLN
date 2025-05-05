# MMLN: An R Package for Mixed Effects Multinomial Regression and Model Diagnostics

Tools for fitting and evaluating Mixed Effects Multinomial Logistic Normal Regression Models, with support for both fixed-effects and mixed-effects formulations. Includes utilities for data simulation, trace visualization, model comparison (DIC), posterior predictive checks, calculating squared Mahalanobis distance residuals (`MDres`) and performing Kolmogorov-Smirnov tests.


## Installation

```r
# Install devtools if you don’t have it already
install.packages("devtools")

# Install this package from GitHub
devtools::install_github("eaegerber/MMLN")

# Or install from local source
devtools::install("/path/to/MMLN")
```

## Quick Start

```r
library(MMLN)

# 1. Simulate a small mixed-effects dataset
set.seed(42)

sim <- simulate_mixed_mln_data(
  m       = 10,          # 10 groups
  n_i     = 10,          # 10 observations per group
  p       = 3,           # 3 fixed covariates (incl. intercept)
  d       = 2,           # 3 outcome categories
  beta    = matrix(c(0.5, -1, 0.2, 0.3, 0.7, -0.4), 3, 2),
  Sigma   = diag(2),
  Phi     = diag(2),
  PA_mean = 200
)

# 2. Fit a fixed-effects MLN model
res_f <- FMLN(
  W            = sim$W,
  X            = sim$X,
  n_iter       = 2000,
  burn_in      = 500,
  thin         = 2,
  proposal     = "normbeta",
  verbose      = TRUE
)

# 3. Fit a mixed-effects MLN model
res_m <- MMLN(
  W            = sim$W,
  X            = sim$X,
  Z            = sim$Z,
  n_iter       = 2000,
  burn_in      = 500,
  thin         = 2,
  proposal     = "normbeta",
  verbose      = TRUE
)

# 4. Trace plots & posterior summaries
beta_chain_array <- simplify2array(res_m$beta_chain)
trace_stats      <- plot_trace_and_summary(beta_chain_array, "beta")
trace_stats
sim$beta
par(mfrow=c(1,1))

# 5. Compute model DICs
ll_chain <- sapply(res_m$y_chain,
                   function(Y) dmnl_loglik(Y, sim$W))
Y_hat   <- alr(compress_counts(sim$W) / rowSums(sim$W))
ll_hat  <- dmnl_loglik(Y_hat, sim$W)
dic_res <- compute_dic(ll_chain, ll_hat)

# 6. Posterior predictive simulation and Mahalanobis residuals
W_pred_list <- lapply(seq_along(res_m$y_chain), function(i) {
  sample_posterior_predictive(X = sim$X,
                              beta = res_m$beta_chain[[i]],
                              Sigma = res_m$sigma_chain[[i]],
                              PA = sim$PA,
                              Z = sim$Z,
                              psi = res_m$psi_chain[[i]],
                              mixed = TRUE
  )
})
resids <- MDres(sim$W, W_pred_list)
summary(resids)

# 7. Compare to incorrect model fit (fixed model fit to mixed data; should show some overdispersion)
W_pred_list_ovd <- lapply(seq_along(res_f$y_chain), function(i) {
  sample_posterior_predictive(X = sim$X,
                              beta = res_m$beta_chain[[i]],
                              Sigma = res_m$sigma_chain[[i]],
                              PA = sim$PA,
                              mixed = FALSE
  )
})
resids_ovd <- MDres(sim$W, W_pred_list_ovd)
summary(resids_ovd)

# 8. Should also show that incorrect model fit has higher DIC
ll_chain_ovd <- sapply(res_f$y_chain,
                   function(Y) dmnl_loglik(Y, sim$W))
dic_res_ovd <- compute_dic(ll_chain_ovd, ll_hat)
dic_res_ovd$DIC > dic_res$DIC
```


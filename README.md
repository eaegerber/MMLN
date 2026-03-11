# MMLN: An R Package for Mixed Effects Multinomial Regression and Model Diagnostics

Tools for fitting and evaluating Mixed Effects Multinomial Logistic Normal Regression Models, with support for both fixed-effects and mixed-effects formulations. Includes utilities for data simulation, trace visualization, model comparison (DIC), posterior predictive checks, calculating squared Mahalanobis distance residuals (`MDres`) and performing Kolmogorov-Smirnov tests.

## System Requirements

This package includes compiled C++ code (via Rcpp) for performance-critical operations. The C++ source is compiled automatically when you install the package, but your system must have a C++ compiler available for R to use. Follow the instructions below for your operating system.

### Windows

Install **Rtools** before installing the package. Rtools provides the compiler toolchain that R needs on Windows.

1. Download Rtools from [https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/)
2. Run the installer with default settings
3. Restart R/RStudio after installation

You can verify the installation by running in R:

```r
Sys.which("make")
# Should return a path like "C:\\rtools44\\usr\\bin\\make.exe"
```

### macOS

Install the **Xcode Command Line Tools**, which include the `clang` C++ compiler:

```bash
xcode-select --install
```

Follow the prompts to complete installation. You can verify with:

```bash
clang++ --version
```

### Linux (Debian/Ubuntu)

Install the `build-essential` meta-package and R development headers:

```bash
sudo apt-get install build-essential r-base-dev
```

For Fedora/RHEL-based systems, use:

```bash
sudo dnf install R-devel gcc-c++ make
```

### A Note on Numerical Precision

Because the package uses compiled C++ arithmetic, results may differ from the pure-R implementation at approximately the `1e-8` level (around the 8th decimal place) depending on your operating system and compiler version. This is expected floating-point behavior and is well within acceptable tolerance, it does not indicate a bug or correctness issue.

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
  Phi     = 5*diag(2),
  n_mean = 200
)

# 2. Fit a fixed-effects MLN model
res_f <- FMLN(
  Y            = sim$Y,
  X            = sim$X,
  n_iter       = 1000,
  burn_in      = 300,
  thin         = 2,
  proposal     = "normbeta",
  verbose      = TRUE
)

# 3. Fit a mixed-effects MLN model
res_m <- MMLN(
  Y            = sim$Y,
  X            = sim$X,
  Z            = sim$Z,
  n_iter       = 1000,
  burn_in      = 300,
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
ll_chain <- sapply(res_m$w_chain,
                   function(W) dmnl_loglik(W, sim$Y))
W_hat   <- alr(compress_counts(sim$Y) / rowSums(sim$Y))
ll_hat  <- dmnl_loglik(W_hat, sim$Y)
dic_res <- compute_dic(ll_chain, ll_hat)

# 6. Posterior predictive simulation and Mahalanobis residuals
Y_pred_list <- lapply(seq_along(res_m$w_chain), function(i) {
  sample_posterior_predictive(X = sim$X,
                              beta = res_m$beta_chain[[i]],
                              Sigma = res_m$sigma_chain[[i]],
                              n = sim$n,
                              Z = sim$Z,
                              psi = res_m$psi_chain[[i]],
                              mixed = TRUE,
                              verbose = FALSE
  )
})
resids <- MDres(sim$Y, Y_pred_list)
summary(resids)

# 7. Compare to incorrect model fit (fixed model fit to mixed data; should show some overdispersion)
Y_pred_list_ovd <- lapply(seq_along(res_f$w_chain), function(i) {
  sample_posterior_predictive(X = sim$X,
                              beta = res_f$beta_chain[[i]],
                              Sigma = res_f$sigma_chain[[i]],
                              n = sim$n,
                              mixed = FALSE,
                              verbose = FALSE
  )
})
resids_ovd <- MDres(sim$Y, Y_pred_list_ovd)
summary(resids_ovd)

# 8. Should also show that incorrect model fit has higher DIC
ll_chain_ovd <- sapply(res_f$w_chain,
                       function(W) dmnl_loglik(W, sim$Y))
dic_res_ovd <- compute_dic(ll_chain_ovd, ll_hat)
dic_res_ovd$DIC > dic_res$DIC
```

## Real-Data Example: Pollen Data

This package includes a helper, `run_pollen_models()`, that reproduces the Gerber \& Craig (2024) pollen analysis by fitting three models to the **pollen** counts (from the **MM** package):

1. **Multinomial logit** via `MGLMreg(dist = "MN")`
2. **Dirichlet multinomial** via `MGLMreg(dist = "DM")`
3. **Fixed-effects MLN** via `FMLN()`

It then draws replicates from the fitted model distributions and computes Mahalanobis‐residuals for each fit.

### Usage

```r
# install real-data dependencies (if needed)
install.packages(c("MM", "MGLM"))

# load package
library(MMLN)

# run the pollen data example
pollen_res <- run_pollen_models(
  n_iter   = 1000,    # total MLN iterations
  burn_in  = 400,     # MLN burn-in
  thin     = 2,       # MLN thinning
  proposal = "normbeta",
  P        = 500      # number of posterior predictive replicates
)

# inspect KS-tests and QQ-plots of Mahalanobis residuals
# should show the well established (Mosimann, 1962) result that there is overdispersion
# DM and MLN should fit better than MLR
summary(pollen_res$resids_mlr)
summary(pollen_res$resids_dm)
summary(pollen_res$resids_mln)
```

## Real-Data Example: MLB Data (with FMLN vs. MMLN)

This package includes a helper, `clean_Lahman_data()`, which cleans MLB data from the Lahman package in order to mimic some of the work from Gerber \& Craig (2021). The function is currently very rigid, but in the future will allow for more flexibility (in terms of which outcomes, subset of data, etc.). We will compare two models:

1. **Fixed-effects MLN** via `FMLN()`
2. **Mixed-effects MLN** via `MMLN()`

To determine if there is a benefit in accounting for individual player random effects in modeling MLB batting outcomes.

### Vignette

```r

library(MMLN)

# 0. Load the cleaned MLB data (post-1960, min. PA of 200, four categories: HR, BB, SO, Other)
baseball_example <- clean_Lahman_data()

# 1. The fixed model
mlb_f <- FMLN(
  Y = baseball_example$Y,
  X = cbind(1, baseball_example$X),
  n_iter = 100,
  burn_in = 30,
  proposal = "normbeta",
  mh_scale = .2, # may want to fiddle with this for better acceptance ratios
  verbose = TRUE
)

# 2. Adding in random effect (will take longer)
mlb_m <- MMLN(
  Y = baseball_example$Y,
  X = cbind(1, baseball_example$X),
  Z = baseball_example$Z,
  n_iter = 100,
  burn_in = 30,
  proposal = "norm", # this runs faster than normbeta, and for small n_iter is not very different
  mh_scale = .1, # may want to fiddle with this for better acceptance ratios
  verbose = TRUE
)

# 3. Posterior predictive simulation and Mahalanobis residuals
Y_pred_list_f <- lapply(seq_along(mlb_f$w_chain), function(i) {
  sample_posterior_predictive(X = cbind(1, baseball_example$X),
                              beta = mlb_f$beta_chain[[i]],
                              Sigma = mlb_f$sigma_chain[[i]],
                              n = baseball_example$PA,
                              mixed = FALSE,
                              verbose = FALSE
  )
})
resids_f <- MDres(baseball_example$Y, Y_pred_list_f)
summary(resids_f)


Y_pred_list_m <- lapply(seq_along(mlb_m$w_chain), function(i) {
  sample_posterior_predictive(X = cbind(1, baseball_example$X),
                              beta = mlb_m$beta_chain[[i]],
                              Sigma = mlb_m$sigma_chain[[i]],
                              n = baseball_example$PA,
                              Z = baseball_example$Z,
                              psi = mlb_m$psi_chain[[i]],
                              mixed = TRUE,
                              verbose = FALSE
  )
})
resids_m <- MDres(baseball_example$Y, Y_pred_list_m)
summary(resids_m)

# 4. Evaluate posterior parameter means of the fixed effects model
post_beta <- apply(simplify2array(mlb_f$beta_chain), c(1,2), mean)
row.names(post_beta) <- c("int", colnames(baseball_example$X))
colnames(post_beta) <- colnames(baseball_example$Y)[1:3]
post_beta

# Note: first column represents covariate effect on HR relative to Other
# Example: Right handed hitters tend to hit more HR than Left handed hitters
#          NL batters tend to walk more than AL batters
#          Taller batters tend SO more than shorter batters
#          The Average Batter (avg. weight, height, age, Lefty, AL, C)

avg_batter_pi <- alr_inv(post_beta[1,])
avg_batter_pi

```

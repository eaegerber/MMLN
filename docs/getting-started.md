# Getting Started with MMLN

This guide is for researchers who want to install and use the MMLN package to fit Bayesian mixed-effects multinomial logistic-normal regression models. For contributor and development workflows, see [development.md](development.md).

---

## Prerequisites

### R

R must be installed on your system. Download it from the official R Project page:
[https://www.r-project.org/](https://www.r-project.org/)

### IDE

RStudio is the recommended environment for working with R packages:
[https://posit.co/download/rstudio-desktop/](https://posit.co/download/rstudio-desktop/)

Alternatively, R can be run directly in a terminal by typing `R` after installation.

### C++ Compiler (required)

MMLN includes compiled C++ code for performance-critical operations. Your system needs a C++ compiler available to R before installation.

**Windows:** Install **Rtools** from [https://cran.r-project.org/bin/windows/Rtools/](https://cran.r-project.org/bin/windows/Rtools/). Run the installer with default settings and restart R/RStudio afterward. Without Rtools, the C++ backend will not compile on Windows.

**macOS:** Install the **Xcode Command Line Tools**:

```bash
xcode-select --install
```

**Linux (Debian/Ubuntu):** Install `build-essential` and R development headers:

```bash
sudo apt-get install build-essential r-base-dev
```

For Fedora/RHEL-based systems:

```bash
sudo dnf install R-devel gcc-c++ make
```

### R Package Dependencies

The following packages are required and will be installed automatically with MMLN:

- `car`, `dplyr`, `magrittr`, `MGLM`, `mvnfast`, `Rcpp`

The `devtools` package is needed for installation:

```r
install.packages("devtools")
```

---

## Installation

Install MMLN directly from GitHub:

```r
devtools::install_github("eaegerber/MMLN")
```

Or install from a local source directory:

```r
devtools::install("/path/to/MMLN")
```

---

## Important: Loading the Package for C++ Support

After installing, load the package with `library()` as usual — but you must **also** call `devtools::load_all()` for the C++ backend to be active:

```r
library(MMLN)
devtools::load_all()
```

> **Without `devtools::load_all()`, the C++-optimized functions will not run.** This is a required step every time you start a new R session with MMLN. Skipping it is one of the most common sources of errors when first using the package.

You can verify everything is working correctly with:

```r
verify_installation()
# Should print: "C++ is working w/ Rcpp"
```

---

## Examples

### Quick Start: Simulated Data

This example generates a synthetic mixed-effects dataset, fits both fixed- and mixed-effects models, and evaluates fit via trace plots, DIC, and Mahalanobis distance residuals.

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

---

### Real-Data Example: Pollen Data

`run_pollen_models()` reproduces the Gerber & Craig (2024) pollen analysis by fitting three models to pollen counts (from the **MM** package) and comparing their fit via Mahalanobis distance residuals:

1. **Multinomial logit** via `MGLMreg(dist = "MN")`
2. **Dirichlet multinomial** via `MGLMreg(dist = "DM")`
3. **Fixed-effects MLN** via `FMLN()`

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

---

### Real-Data Example: MLB Data (FMLN vs. MMLN)

`clean_Lahman_data()` prepares post-1960 MLB batting records from the **Lahman** package. This example compares a fixed-effects and mixed-effects MLN model to assess whether accounting for individual player random effects improves fit.

Install the Lahman package first if you haven't already:

```r
install.packages("Lahman")
```

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

---

## Interpreting Output

Both `FMLN()` and `MMLN()` return a named list of MCMC chains representing posterior samples:

| Field | Description |
|---|---|
| `beta_chain` | List of posterior samples for the fixed-effects coefficient matrix |
| `sigma_chain` | List of posterior samples for the latent covariance matrix Σ |
| `w_chain` | List of posterior samples for the latent ALR-transformed variables |
| `mhaccept_chain` | Metropolis-Hastings acceptance rates per iteration |
| `psi_chain` | *(MMLN only)* Posterior samples for group-level random intercepts |

Each element of a `_chain` list corresponds to one saved MCMC iteration (after burn-in and thinning). Use `plot_trace_and_summary()` to visualize convergence, `compute_dic()` for model comparison, and `MDres()` + `summary()` to assess model fit via Mahalanobis distance residuals.

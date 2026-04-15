# MMLN: An R Package for Mixed Effects Multinomial Regression and Model Diagnostics

Tools for fitting and evaluating Mixed Effects Multinomial Logistic Normal Regression Models, with support for both fixed-effects and mixed-effects formulations. Includes utilities for data simulation, trace visualization, model comparison (DIC), posterior predictive checks, calculating squared Mahalanobis distance residuals (`MDres`) and performing Kolmogorov-Smirnov tests.

## Documentation

- [Getting Started](docs/getting-started.md) — installation, setup, and examples **for researchers**
- [Architecture](docs/architecture.md) — system design and package internals
- [Development](docs/development.md) — contributor workflow guide **for developers**

## System Requirements

This package includes compiled C++ code (via Rcpp). A C++ compiler must be available on your system before installation. See [Getting Started](docs/getting-started.md) for platform-specific setup instructions (Rtools on Windows, Xcode CLI on macOS, `build-essential` on Linux).

## Installation

Install devtools if you don't have it already

```r
install.packages("devtools")
```

Install MMLN directly from GitHub:

```r
devtools::install_github("eaegerber/MMLN")
```

**Or** install from a local source directory:

```r
devtools::install("/path/to/MMLN")
```

### Building from a Local Clone

If you cloned this repository and are installing from source, you must first generate the Rcpp bridge files before installing. These files are auto-generated and not tracked in git:

```r
Rcpp::compileAttributes("/path/to/MMLN")
devtools::install("/path/to/MMLN")
```

Without this step, C++ functions like `mh_update_cpp()` will not be available and you’ll see `could not find function` errors.

## Core Workflow

The typical workflow has three stages: prepare your data, fit a model, and evaluate the fit. Every function mentioned below is exported and documented. Run `?function_name` in R for full details.

### 1. Data Preparation

If you're working with your own count data, you need an N × J matrix of multinomial counts `Y` and an N × p covariate matrix `X`. For mixed-effects models you also need a group indicator matrix `Z` (N × m, where m is the number of groups).

`simulate_mixed_mln_data()` generates synthetic datasets with known parameters, which is useful for testing or understanding the model before applying it to real data. `clean_Lahman_data()` and `run_pollen_models()` provide ready-made real-data examples.

Two helper functions handle the transformation between counts and the log-ratio space the model operates in. `compress_counts()` applies the Smithson-Verkuilen correction for zero counts, and `alr()` / `alr_inv()` convert between probability vectors and additive log-ratio coordinates (the last category is always treated as the baseline).

### 2. Model Fitting

`FMLN()` fits a fixed-effects multinomial logistic-normal model and `MMLN()` fits the mixed-effects version with group-level random intercepts. Both are Gibbs samplers with Metropolis-Hastings updates for the latent variables.

Key parameters shared by both functions:

- `n_iter`, `burn_in`, `thin` control the MCMC chain length and which samples are saved
- `proposal` selects the MH proposal distribution: `"norm"`, `"beta"`, or `"normbeta"` (default). `"normbeta"` is generally preferred; `"norm"` is faster per iteration and reasonable for short exploratory runs
- `mh_scale` tunes the proposal step size — adjust this to get reasonable acceptance ratios (printed during fitting when `verbose = TRUE`)

Both functions return a list containing the saved MCMC chains: `beta_chain`, `sigma_chain`, `w_chain`, and `mhaccept_chain`. `MMLN` additionally returns `psi_chain` (the random intercept chain).

### 3. Model Evaluation

`plot_trace_and_summary()` produces trace plots and posterior summary statistics for any parameter chain. Use this to check convergence before interpreting results.

`compute_dic()` computes the Deviance Information Criterion for model comparison. You supply a vector of log-likelihood values across the chain and the log-likelihood evaluated at the posterior mean.

`sample_posterior_predictive()` draws replicate datasets from the fitted model. Pass these to `MDres()` to compute squared Mahalanobis distance residuals, which measure how well the model's predictive distribution matches the observed data. `summary()` on the result runs a Kolmogorov-Smirnov test and produces a QQ plot — if the model fits well, the residuals should be approximately chi-squared distributed.

### Utility Functions

`dmnl_loglik()` computes the multinomial logistic-normal log-likelihood for a given latent W and count matrix Y. `ytopstar()` and `pstartoy()` convert between raw counts and the continuous latent scale. The proposal distribution functions (`normbetapropdist`, `normbetaloglike`, `betapropdist`, `betaloglike`) are exported for transparency but are called internally by the samplers, you generally don't need to use them directly.

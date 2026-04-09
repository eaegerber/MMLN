# System Architecture

## Overview

MMLN is an R package for fitting **Mixed Effects Multinomial Logistic Normal (MLN) Regression Models**. It provides Gibbs samplers for both fixed-effects (`FMLN`) and mixed-effects (`MMLN`) formulations, along with model diagnostics via Mahalanobis distance residuals.

## Repository Structure

```
MMLN/
├── DESCRIPTION                    # Package metadata, dependencies, LinkingTo
├── NAMESPACE                      # Exported functions, useDynLib registration
├── R/
│   ├── mln_functions.R            # FMLN() and MMLN() Gibbs samplers
│   ├── mln_helpers.R              # Transforms, log-likelihood, simulation, R-side proposals
│   ├── multi_res.R                # MDres(), DIC, posterior predictive, traceplots
│   ├── real_data_examples.R       # Pollen and MLB data examples
│   └── RcppExports.R             # Auto-generated .Call() wrappers
├── src/
│   ├── mh_update.cpp              # C++ Metropolis-Hastings proposals
│   ├── mdres_core.cpp             # C++ Mahalanobis distance residuals
│   ├── utils.cpp                  # Shared C++ utilities
│   ├── RcppExports.cpp           # Auto-generated C registration table
│   ├── Makevars                   # Compiler flags, LAPACK/BLAS linking
│   └── Makevars.win               # Windows build flags
├── man/                           # Roxygen-generated documentation (.Rd)
├── optimization-testing/
│   ├── Makefile                   # Test orchestration
│   ├── config.R                   # Scenario definitions and tolerances
│   ├── scripts/                   # generate_data, generate_reference, run_test, compare
│   ├── data/                      # Generated synthetic datasets
│   ├── fixtures/                  # Upstream reference outputs
│   └── results/                   # Fork test outputs
└── profiling/                     # Performance profiling reports (HTML)
```

## C++ Performance Layer

Profiling on the MLB baseball dataset identified **Metropolis-Hastings proposal updates** and **Mahalanobis distance residual computation** as the dominant runtime bottlenecks. These paths are implemented in compiled C++ via Rcpp:

- `mh_update_normbeta_cpp` and `mh_update_beta_cpp` — MH proposals across all `N × d` elements per iteration
- `mdres_core_cpp` — full residual pipeline in a single compiled pass (count compression, ALR transform, covariance estimation, Cholesky solve, Mahalanobis distances, quantile-residual inversion)
- `clamp_vec`, `col_means`, `sym_eigenvalues` — shared utilities in `utils.cpp`

Numerical precision: results match a pure-R reference implementation to within `~1e-8` per operation, with cumulative MCMC chain drift bounded at `~1e-6` over 1000 iterations.

## R Layer

### Core Samplers (`R/mln_functions.R`)

- **`FMLN()`** — Fixed-effects Gibbs sampler. Each iteration: compute `Mu = Xβ`, MH-update W (via C++), accept/reject, sample β (conjugate normal), sample Σ (inverse-Wishart).
- **`MMLN()`** — Mixed-effects Gibbs sampler. Adds random intercepts ψ per group: `Mu = Xβ + Zψ`. Additional steps sample ψ (conditional normal per group) and Φ (inverse-Wishart on ψ covariance).

### Helpers (`R/mln_helpers.R`)

- **Transforms:** `alr()` / `alr_inv()` (additive log-ratio ↔ simplex), `compress_counts()` (zero correction)
- **Log-likelihood:** `dmnl_loglik(W, Y)`
- **R-side proposal functions:** `betapropdist()`, `normbetapropdist()` and their log-likelihoods (used for the `"norm"` proposal path, which remains pure-R)
- **Data simulation:** `simulate_mixed_mln_data()` generates synthetic mixed-effects MLN datasets

### Diagnostics (`R/multi_res.R`)

- **`MDres()`** — wraps `mdres_core_cpp()`, returns S3 class `"mdres"`
- **`summary.mdres()`** — KS test for normality + QQ plot of quantile residuals
- **`compute_dic()`** — Deviance Information Criterion from log-likelihood chain
- **`sample_posterior_predictive()`** — simulate replicate counts under the fitted model
- **`plot_trace_and_summary()`** — MCMC traceplots with posterior mean/CI

### Real-Data Examples (`R/real_data_examples.R`)

- **`run_pollen_models()`** — fits MN, DM, and MLN to pollen counts (Mosimann 1962; Gerber & Craig 2024)
- **`clean_Lahman_data()`** — prepares MLB batting data (HR, BB, SO, Other) with player random effects

## C++ Modules

All C++ source lives in `src/` and is compiled into a shared library (`MMLN.so`) at package install time.

### `mh_update.cpp` — Metropolis-Hastings Proposals

Two proposal strategies for updating the latent log-ratio matrix **W**:

| Function                 | Strategy              | How it works                                                                                                                        |
| ------------------------ | --------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `mh_update_normbeta_cpp` | Normal approx to Beta | Computes Beta shape parameters from prior + data, uses digamma/trigamma to find the mode and variance, draws from a Normal proposal |
| `mh_update_beta_cpp`     | Beta in sigmoid space | Maps W through sigmoid, draws new probabilities from Beta(a*, b*), maps back via logit with Jacobian correction                     |

Both iterate over all `N × d` elements, computing per-element proposal densities. The accept/reject decision itself remains in R (it depends on the full log-likelihood and prior).

### `mdres_core.cpp` — Mahalanobis Distance Residuals

Per-observation pipeline:

1. Apply Smithson-Verkuilen zero-count correction
2. ALR-transform observed and posterior predictive counts
3. Estimate per-observation covariance from `P` posterior samples
4. Check for singular covariance (min eigenvalue < 1e-8)
5. Cholesky factorization via Armadillo with 1e-8 diagonal jitter
6. Compute Mahalanobis distance: `‖L⁻¹(x - μ)‖²`
7. Binary-search ECDF to find quantile rank among predictive distances
8. Invert to standard normal via randomized quantile residual

### `utils.cpp` — Shared Utilities

- `clamp_vec(x, lo, hi)` — element-wise clamping (Armadillo)
- `col_means(X)` — column means (Armadillo)
- `sym_eigenvalues(X)` — eigenvalues of symmetric matrix (Eigen)

## How C++ Fits into R with Rcpp

User-facing R functions like `FMLN()`, `MMLN()`, and `MDres()` delegate to C++ by calling thin R wrappers in `R/RcppExports.R`, which in turn invoke the compiled C++ via `.Call()`. These wrappers and their corresponding C registration table (`src/RcppExports.cpp`) are auto-generated by `Rcpp::compileAttributes()`, which scans `src/*.cpp` for `// [[Rcpp::export]]` annotations. At package install time, R's build system compiles all C++ source into a shared library (`MMLN.so`) using the flags in `src/Makevars`, and `NAMESPACE` declares `useDynLib(MMLN, .registration = TRUE)` to load it at runtime. Matrix arguments are passed between R and C++ with zero-copy mappings through Armadillo and Eigen types.

## Rcpp Sub-Libraries

| Library           | Role                                                                             | Where used                                     |
| ----------------- | -------------------------------------------------------------------------------- | ---------------------------------------------- |
| **Rcpp**          | R ↔ C++ type marshalling (`NumericMatrix`, `List`, etc.)                         | All `.cpp` files                               |
| **RcppArmadillo** | Armadillo linear algebra (`arma::mat`, `arma::cov`, `arma::chol`, `arma::clamp`) | `mh_update.cpp`, `mdres_core.cpp`, `utils.cpp` |
| **RcppEigen**     | Eigen linear algebra (`SelfAdjointEigenSolver` for eigenvalue checks)            | `mdres_core.cpp`, `utils.cpp`                  |

Both are declared in `LinkingTo:` in DESCRIPTION. `src/Makevars` links LAPACK/BLAS and suppresses Eigen SSE/SIMD warnings.

## Optimization Testing Infrastructure

The `optimization-testing/` directory contains an end-to-end comparison suite that verifies the C++ implementation produces numerically equivalent results to a pure-R reference implementation.

**Pipeline** (orchestrated by `Makefile`):

1. **`generate_data.R`** — creates 10 synthetic datasets from `config.R` scenario definitions
2. **`generate_reference.R`** — installs the upstream package, runs FMLN/MMLN/MDres, saves fixture outputs
3. **`run_test.R`** — runs the same models on the same data using the C++ implementation
4. **`compare.R`** — statistically compares C++ vs reference outputs against defined tolerances

**Scenarios** stress different numerical edge cases: zero cells, tiny/large counts, high dimensionality (d=8), near-singular covariance (ρ=0.92), large coefficients (overflow risk), unbalanced groups (n=1..80), many groups (m=500), and dominant categories (~90% probability mass).

**Tolerance tiers:**

| Constant     | Value | Scope                                           |
| ------------ | ----- | ----------------------------------------------- |
| `TOL_EXACT`  | 0     | Pure-R proposal paths must be bitwise identical |
| `TOL_FLOAT`  | 1e-8  | Single C++ operation vs R                       |
| `TOL_CHAIN`  | 1e-6  | Cumulative drift over 1000 MCMC iterations      |
| `TOL_ACCEPT` | 1e-10 | Acceptance ratio precision                      |
| `TOL_MDRES`  | 1e-6  | MDres downstream from chain differences         |

Reference fixtures and test results must be generated on the same machine due to cross-platform LAPACK non-determinism (documented in `CROSS_PLATFORM_DIVERGENCE.md`).

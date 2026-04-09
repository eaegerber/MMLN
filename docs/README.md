# MMLN Local Development Guide

This guide is for developing `MMLN` locally, with emphasis on:

- day-to-day R package iteration
- adding and validating C++ improvements
- running and profiling from RStudio

## Repository Map

Key directories and files:

- `R/`: user-facing R functions and wrappers
- `src/`: C++ implementation (`mdres_core.cpp`, `mh_update.cpp`, `utils.cpp`, etc.)
- `man/`: generated `.Rd` docs
- `testing/`: local testing and experimentation scripts
- `profiling/`: profiling scripts and saved profile outputs
- `DESCRIPTION`, `NAMESPACE`: package metadata and exports
- `src/Makevars`, `src/Makevars.win`: compiler flags

## Local Setup

## 1) Prerequisites

Install:

- a recent R version
- RStudio (recommended for interactive development)
- build toolchain for your OS (Rtools on Windows, Xcode CLI on macOS, build-essential on Linux)

## 2) Install development dependencies

In R:

```r
install.packages(c(
  "devtools",
  "roxygen2",
  "testthat",
  "profvis"
))
```

Optional dependencies used by examples:

```r
install.packages(c("Lahman", "MM", "MGLM"))
```

## 3) Install package in dev mode

From repo root in R:

```r
devtools::load_all()
```

Use `load_all()` while iterating, then do a clean install check before merge:

```r
devtools::check()
```

## Running in RStudio

Recommended RStudio workflow:

1. Open the project folder as an RStudio Project.
2. Run `devtools::load_all()` after code changes.
3. Run a quick smoke test (for example from `testing/` scripts).
4. If you changed exported functions or docs, run:

```r
devtools::document()
```

5. Re-run `devtools::check()` before opening a PR.

Useful shortcuts:

- **Build > Load All**
- **Build > Check Package**
- **Build > Test Package**

## How to Add More C++

This package uses `Rcpp` + `RcppArmadillo` + `RcppEigen`.

## Typical pattern

1. Add or update implementation in `src/*.cpp`.
2. Expose function(s) to R via `Rcpp::export` where appropriate.
3. Add or update an R wrapper in `R/` if needed.
4. Run:

```r
Rcpp::compileAttributes()
devtools::document()
devtools::load_all()
```

5. Validate with targeted tests and one end-to-end script.

## C++ coding notes for this repo

- Keep matrix/vector dimensions explicit and validated early.
- Prefer stable numeric transforms (`log-sum-exp`, bounded probabilities, etc.) to reduce under/overflow.
- Preserve reproducibility behavior expected by existing MCMC utilities.
- Keep interface contracts unchanged unless intentionally making an API change.

## Build flags

`src/Makevars` and `src/Makevars.win` currently include:

- `PKG_CXXFLAGS = -Wno-ignored-attributes`

This suppresses an Eigen-specific warning class on newer compilers while leaving other warnings active.

## Validating C++ Improvements

When changing C++ paths, do all of the following:

1. **Correctness checks**
   - compare outputs to prior behavior on fixed seeds
   - verify expected tolerances for floating-point differences

2. **Performance checks**
   - profile before/after on same input and seed
   - record wall time and where hotspots moved

3. **Stability checks**
   - run multiple chains / representative datasets
   - confirm acceptance rates and diagnostics still look reasonable

4. **Package checks**
   - run `devtools::check()`
   - confirm no new warnings from your code path

## Profiling with `profvis`

Use `profvis` for interactive hotspot analysis in R or RStudio.

Minimal example:

```r
library(profvis)
library(MMLN)

prof <- profvis({
  # Replace with representative workload
  set.seed(1)
  sim <- simulate_mixed_mln_data(
    m = 8, n_i = 8, p = 3, d = 2,
    beta = matrix(c(0.5, -1, 0.2, 0.3, 0.7, -0.4), 3, 2),
    Sigma = diag(2),
    Phi = 5 * diag(2),
    n_mean = 100
  )

  fit <- FMLN(
    Y = sim$Y,
    X = sim$X,
    n_iter = 300,
    burn_in = 100,
    thin = 2,
    proposal = "normbeta",
    verbose = FALSE
  )
})

prof
```

Tips:

- warm up once before profiling to avoid first-run effects
- profile realistic workloads (too small can hide true bottlenecks)
- store profile artifacts in `profiling/` with clear names

## Other Cool Tools

In addition to `profvis`, these are useful:

- `bench`: micro/mid-size benchmarking in R
- `microbenchmark`: fast comparative timing for small kernels
- `lintr`: style and static checks for R code quality
- `goodpractice`: package-level quality checks
- `covr`: test coverage reports
- `valgrind` (Linux/macOS): native memory diagnostics for C++ paths
- `gdb`/`lldb`: step-level native debugging when crashes occur

Example benchmark pattern:

```r
# install.packages("bench")
library(bench)

mark(
  baseline = {
    # old path / behavior
  },
  candidate = {
    # new path / behavior
  },
  iterations = 20,
  check = FALSE
)
```

## Suggested Pre-PR Checklist

- code compiles cleanly on your platform
- docs regenerated (`devtools::document()`) when APIs changed
- key scripts in `testing/` still run
- C++ changes are validated for both correctness and speed
- `devtools::check()` passes

## Notes on Numerical Differences

Small cross-platform numeric differences are expected with compiled C++ due to compiler and BLAS/LAPACK variation. Treat tiny differences as normal unless they materially affect diagnostics, rankings, or substantive inference.

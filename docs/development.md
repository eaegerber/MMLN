# MMLN Development and Contribution Guide

This guide covers how to set up a development environment, add new R or C++ code, and validate changes before opening a PR. For a full picture of the package's structure and design decisions, see [architecture.md](architecture.md) first.

---

## Setup

### Prerequisites

- A recent R version
- RStudio (recommended for interactive development)
- Build toolchain for your OS: Rtools on Windows, Xcode CLI on macOS, `build-essential` on Linux

### Development dependencies

```r
install.packages(c(
  "devtools",
  "roxygen2",
  "testthat",
  "profvis"
))
```

Optional — needed to run real-data examples:

```r
install.packages(c("Lahman", "MM", "MGLM"))
```

### Load the package for development

```r
devtools::load_all()
```

Run this after any code change to reload the package in your current session. Before opening a PR, do a clean check:

```r
devtools::check()
```

### RStudio shortcuts

- **Build > Load All** — equivalent to `devtools::load_all()`
- **Build > Check Package** — equivalent to `devtools::check()`
- **Build > Test Package** — equivalent to `devtools::test()`

---

## Adding New R Code

The R layer is organized into four files — see [architecture.md](architecture.md) for what belongs where. The typical workflow:

1. Add or modify a function in the appropriate `R/*.R` file.
2. Document it with a roxygen2 block (`#' @export`, `@param`, `@return`, etc.).
3. Regenerate docs and reload:

```r
devtools::document()
devtools::load_all()
```

4. Add a test in `tests/testthat/` or a script in `testing/` and verify it runs.
5. Run `devtools::check()` before submitting.

Key conventions (also in [architecture.md](architecture.md)):

- `Y` is always a count matrix; pass it through `compress_counts()` before ALR-transforming.
- Latent variables live in ALR space (`d = ncol(Y) - 1`).
- The `proposal` argument controls MH behavior: `"norm"`, `"beta"`, or `"normbeta"`.

---

## Adding New C++ Code

The package uses `Rcpp` + `RcppArmadillo` + `RcppEigen`. See [architecture.md](architecture.md) for how the C++ modules are structured and how they connect to R via `.Call()`.

### Typical pattern

1. Add or update implementation in `src/*.cpp`. New utility functions belong in `utils.cpp`; new performance-critical paths should get their own file.
2. Annotate functions to expose to R:

```cpp
// [[Rcpp::export]]
arma::vec my_function(arma::mat X) { ... }
```

3. Regenerate the Rcpp wrappers, rebuild, and reload:

```r
Rcpp::compileAttributes()
devtools::document()
devtools::load_all()
```

4. Validate — see [Validating C++ Changes](#validating-c-changes) below.

### C++ coding notes

- Keep matrix/vector dimensions explicit and validated early.
- Prefer numerically stable transforms (`log-sum-exp`, bounded probabilities) to reduce under/overflow risk.
- Don't break MCMC reproducibility — existing utilities depend on deterministic behavior given a seed.
- Keep R-facing interface contracts stable unless you're intentionally making an API change.

### Build flags

`src/Makevars` and `src/Makevars.win` set:

```
PKG_CXXFLAGS = -Wno-ignored-attributes
```

This suppresses an Eigen-specific warning on newer compilers while leaving all other warnings active. Add new flags here if needed.

---

## Validating C++ Changes

When modifying any C++ path, work through all four of these before opening a PR:

1. **Correctness** — compare outputs to prior behavior on fixed seeds; verify floating-point differences stay within expected tolerances (`~1e-8` per operation, `~1e-6` over 1000 iterations).
2. **Performance** — profile before and after on the same input and seed; record where hotspots moved.
3. **Stability** — run multiple chains on representative datasets; confirm acceptance rates and diagnostics still look reasonable.
4. **Package check** — `devtools::check()` with no new warnings from your code path.

---

## Profiling

Use `profvis` for interactive hotspot analysis:

```r
library(profvis)
library(MMLN)

prof <- profvis({
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
- Warm up once before profiling to avoid first-run overhead.
- Profile realistic workloads — small inputs can hide true bottlenecks.
- Save profile artifacts to `profiling/` with descriptive names.

---

## Helpful Tools

| Tool | Purpose |
|---|---|
| `profvis` | Interactive flame graph profiling in RStudio |
| `bench` | Mid-size benchmarking with memory tracking |
| `microbenchmark` | Fast comparative timing for small kernels |
| `lintr` | R style and static analysis |
| `goodpractice` | Package-level quality checks |
| `covr` | Test coverage reports |
| `valgrind` (Linux/macOS) | Native memory diagnostics for C++ |
| `gdb` / `lldb` | Step-level native debugging for C++ crashes |

Benchmark pattern with `bench`:

```r
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

---

## Pre-PR Checklist

- [ ] Code compiles cleanly on your platform
- [ ] `devtools::document()` run if any exported functions or docs changed
- [ ] Scripts in `testing/` still run end-to-end
- [ ] C++ changes validated for correctness, performance, and stability
- [ ] `devtools::check()` passes with no new warnings

---

## Numerical Differences

Small cross-platform numeric differences are expected with compiled C++ due to compiler and BLAS/LAPACK variation. Treat differences below `~1e-6` as normal unless they materially affect diagnostics, rankings, or substantive inference.

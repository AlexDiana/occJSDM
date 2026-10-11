# Spatial Recoverability Pilot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and run the oracle-only positive-control study that maps, for two kinds of spatial truth, the site counts and ranges at which held-out spatial structure is recoverable from single-season binary occupancy data.

**Architecture:** A standalone R study under `dev/simstudy/spatial-recoverability-pilot/`, sharing no code with the occJSDM package. Six small R files (kernel, generator, oracle, predict, score, run) plus one new C++ joint elliptical-slice sampler; the truth 2 oracle reuses the sweep's `ellipse.cpp` unmodified. Fits run as forked single-threaded workers under the pinned R 4.5.0 runtime, with thinned draws stored on the research volume and compact results committed to the repo.

**Tech Stack:** R 4.5.0 (pinned private home from the first comparison), Rcpp 1.1.2, RcppArmadillo 15.4.0.1, posterior 1.7.0, testthat 3.3.2, base `parallel::mclapply`.

**Spec:** `dev/superpowers/specs/2026-10-10-spatial-recoverability-pilot-design.md`

## Global Constraints

- Markdown files in this repo: one line per paragraph, no pipe tables, no em-dashes, escape `-\>` and `\~` in prose, never inside code.
- Interpreter: always the pinned launcher `dev/simstudy/spatial-recoverability-pilot/run-r`, never system `Rscript` (system R is 4.6.1). Single-threaded BLAS and OpenMP via the launcher's environment.
- Nothing in `R/`, `src/`, priors, defaults or lessons changes.
- Squared-exponential kernel `exp(-d^2 / (2 l^2))`; half-correlation distance `l * sqrt(2 * log(2))`.
- Spatial SD 1 for every species under both truths; truth 1 split 90% shared, 10% independent.
- Truth 2 trait `d_s = log(range)`, values `log(0.03)`, `log(0.08)`, `log(0.15)`.
- Site counts 100, 300, 900; truth 1 ranges 0.03, 0.08, 0.15; 10 communities per cell; two arrangements per community.
- Clustered centre separation `0.2 * sqrt(10 / k)` for `k = n / 10` clusters of 10 sites, radius 0.02.
- Lattice 40 by 40, cell centres, 1,600 points.
- MCMC: 4 chains, 3,000 burn-in, 4,000 retained; one extension to 8,000 burn-in, 16,000 retained; screen Rhat at most 1.01, bulk and tail ESS at least 400 on field average, field-truth projection and probability average per species.
- Numerical conditioning: Cholesky with jitter 1e-8 when it succeeds and the condition number is at most 1e12, otherwise eigen route with relative truncation 1e-10; lattice conditional always by eigen truncation; reconstruction tolerance 1e-5 absolute; every route and condition number logged.
- Thresholds (pattern median at least 0.6 and at least 8 of 10 at least 0.4; slope median in 0.5 to 1.5; headroom mean at least 1 percentage point with 95% interval excluding zero) are written to the study config in `prepare` before any community is generated, and never edited.
- No community, seed or input from earlier studies is reused; no community is discarded after generation.
- Long runs launched detached with `nohup`; raw job files on `/Volumes/occJSDM Research/occJSDM/spatial-recoverability-20261010/jobs`; run index at `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010/` (gitignored).
- Commit after each task on branch `claude/spatial-recoverability-pilot`; commit messages end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`.

## Review Focus

- A community whose kernel at all 3,400 points cannot be Cholesky-factorised: the generator must fall through to the eigen route and record it, not stop. Test in Task 2.
- A joint-sampler state with the independent scale set to zero and one shared field: must reproduce the existing sampler's draws exactly under the same seed, so a mismatch is a bug and not RNG drift. Test in Task 3.
- An oracle fit that fails the screen after extension: the job must save with status `unresolved`, the scorer must carry that status, and the gate must count the cell as unresolved, never as fail. Test in Task 5 and Task 6.
- A lattice conditional with a negative eigenvalue: truncation must drop it and the recorded reconstruction error must stay under 1e-5, otherwise the job stops with a clear error. Test in Task 4.
- A rerun of `fit` after an interruption: existing job files must be reused unchanged and only missing jobs run, so no community is ever refitted with a different state. Test in Task 6.

---

### Task 1: Study scaffold, launcher and seeds

**Files:**
- Create: `dev/simstudy/spatial-recoverability-pilot/run-r`
- Create: `dev/simstudy/spatial-recoverability-pilot/seeds.R`
- Create: `dev/simstudy/spatial-recoverability-pilot/test-seeds.R`

**Interfaces:**
- Produces: `pilot_seed(label)` returning a distinct positive integer seed for a character label; `pilot_label(...)` joining parts with `-`.

- [ ] **Step 1: Write the launcher**

```sh
#!/bin/sh
# Pinned R 4.5.0 for the spatial recoverability pilot. Reuses the private
# R home built for the first spatial comparison; never follows the system
# framework symlink (system R is 4.6.1). Single-threaded numerics.
set -eu
export R_HOME="${PILOT_R_HOME:-/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-package-comparison-20261008/runtime/R-home}"
export R_SHARE_DIR="$R_HOME/share" R_INCLUDE_DIR="$R_HOME/include" R_DOC_DIR="$R_HOME/doc"
export R_LIBS_USER="/Users/douglasyu/Library/R/arm64/4.5/library"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
exec "$R_HOME/bin/exec/R" --vanilla --slave "$@"
```

- [ ] **Step 2: Make it executable and confirm the version**

Run: `chmod +x dev/simstudy/spatial-recoverability-pilot/run-r && dev/simstudy/spatial-recoverability-pilot/run-r -e 'cat(R.version.string, "\n")'`
Expected: `R version 4.5.0 (2025-04-11)`

- [ ] **Step 3: Write the failing seed test**

```r
# Run: dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-seeds.R
library(testthat)
source("dev/simstudy/spatial-recoverability-pilot/seeds.R")
source("dev/simstudy/spatial-design-sweep/generator.R")   # for sweep_seed()

test_that("seeds are positive integers and distinct across labels", {
  labels <- c("sites-t1-r0.03-n100-rep1", "sites-t1-r0.03-n100-rep2", "data-t2-n900-rep10",
              "oracle-t2-n300-rep4-clustered-species07")
  s <- vapply(labels, pilot_seed, integer(1))
  expect_true(all(s > 0L)); expect_false(anyDuplicated(s) > 0)
})

test_that("pilot seeds do not collide with the sweep's seed namespace", {
  sweep <- unlist(lapply(c("spatial-design-sites", "spatial-design-data", "spatial-design-detection",
                           "spatial-design-survey-spread", "spatial-design-survey-clustered"),
                         function(l) vapply(1:3, function(r) sweep_seed(l, r), integer(1))))
  pilot <- vapply(sprintf("sites-t1-r0.03-n100-rep%d", 1:10), pilot_seed, integer(1))
  expect_length(intersect(sweep, pilot), 0L)
})

test_that("pilot_label joins parts", {
  expect_identical(pilot_label("sites", "t2", "n100", "rep3"), "sites-t2-n100-rep3")
})
cat("Seed tests passed.\n")
```

- [ ] **Step 4: Run to verify it fails**

Run: `dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-seeds.R`
Expected: error, `seeds.R` not found.

- [ ] **Step 5: Implement seeds.R**

```r
# Seed namespace for the spatial recoverability pilot. Labels are hashed by a
# position-weighted byte sum so that every random stream is reproducible and
# distinct from the sweep's label*1000+replicate scheme.
pilot_label <- function(...) paste(c(...), collapse = "-")

pilot_seed <- function(label) {
  stopifnot(is.character(label), length(label) == 1L, nzchar(label))
  b <- as.numeric(utf8ToInt(paste0("recoverability-pilot:", label)))
  as.integer((sum(b * seq_along(b)) * 7919 + sum(b) + 1) %% 2147483646 + 1)
}
```

- [ ] **Step 6: Run to verify it passes**

Run: same command. Expected: `Seed tests passed.`

- [ ] **Step 7: Commit**

```bash
git add dev/simstudy/spatial-recoverability-pilot dev/superpowers/specs/2026-10-10-spatial-recoverability-pilot-design.md dev/superpowers/plans/2026-10-10-spatial-recoverability-pilot.md
git commit -m "Add spatial recoverability pilot spec, plan, launcher and seeds"
```

---

### Task 2: Kernel utilities and stable factorisation

**Files:**
- Create: `dev/simstudy/spatial-recoverability-pilot/kernel.R`
- Create: `dev/simstudy/spatial-recoverability-pilot/test-kernel.R`

**Interfaces:**
- Produces: `sq_exp(a, b, range)`; `half_correlation(range)`; `stable_factor(K, jitter = 1e-8, rel_tol = 1e-10, cond_ceiling = 1e12)` returning `list(lower, route, condition, rank, error)` with `lower` always n by n (zero-padded on the eigen route); `lattice_conditional(xy_sites, xy_lattice, range, jitter, rel_tol)` returning `list(A, lower, rank, min_eigen, error, route)`; `gauss_hermite(n)` returning `list(nodes, weights)`.

- [ ] **Step 1: Write the failing tests**

```r
# Run: dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-kernel.R
library(testthat)
source("dev/simstudy/spatial-recoverability-pilot/kernel.R")

test_that("half-correlation distance matches the spec", {
  expect_equal(half_correlation(0.03), 0.03 * sqrt(2 * log(2)))
  expect_equal(round(half_correlation(c(0.03, 0.08, 0.15)), 3), c(0.035, 0.094, 0.177))
})

test_that("stable_factor takes the Cholesky route on a well-conditioned kernel", {
  set.seed(1); xy <- matrix(runif(200), 100, 2)
  f <- stable_factor(sq_exp(xy, xy, 0.03))
  expect_identical(f$route, "cholesky"); expect_equal(dim(f$lower), c(100L, 100L))
  expect_equal(f$lower %*% t(f$lower), sq_exp(xy, xy, 0.03) + diag(1e-8, 100), tolerance = 1e-10)
})

test_that("stable_factor falls through to the eigen route and pads to square", {
  set.seed(2); xy <- rbind(matrix(runif(600), 300, 2), matrix(runif(600), 300, 2))
  K <- sq_exp(xy, xy, 0.5)   # extremely smooth: Cholesky fails or condition exceeds ceiling
  f <- stable_factor(K)
  expect_identical(f$route, "eigen"); expect_equal(dim(f$lower), c(600L, 600L))
  expect_lt(f$error, 1e-6); expect_lt(f$rank, 600L)
})

test_that("lattice_conditional reproduces the conditional covariance within tolerance", {
  set.seed(7); xy <- matrix(runif(1800), 900, 2)
  g <- (1:40 - 0.5) / 40; lat <- as.matrix(expand.grid(x = g, y = g))
  cf <- lattice_conditional(xy, lat, 0.15)
  expect_lt(cf$error, 1e-5); expect_lt(cf$rank, 1600L)
  expect_equal(dim(cf$A), c(1600L, 900L)); expect_equal(nrow(cf$lower), 1600L)
})

test_that("Gauss-Hermite integrates a logistic-normal mean to 1e-6", {
  gh <- gauss_hermite(40L)
  gh_value <- function(eta, sd) sum(gh$weights * plogis(eta + sqrt(2) * sd * gh$nodes)) / sqrt(pi)
  for (eta in c(-3, -0.5, 0, 1.2)) {
    exact <- integrate(function(f) plogis(eta + f) * dnorm(f, 0, 1), -Inf, Inf)$value
    expect_equal(gh_value(eta, 1), exact, tolerance = 1e-6)
  }
})
cat("Kernel tests passed.\n")
```

- [ ] **Step 2: Run to verify it fails**

Run: `dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-kernel.R`
Expected: error, `kernel.R` not found.

- [ ] **Step 3: Implement kernel.R**

```r
# Kernel, factorisation and quadrature utilities for the recoverability pilot.
sq_exp <- function(a, b, range) {
  exp(-(outer(a[, 1], b[, 1], "-")^2 + outer(a[, 2], b[, 2], "-")^2) / (2 * range^2))
}

half_correlation <- function(range) range * sqrt(2 * log(2))

stable_factor <- function(K, jitter = 1e-8, rel_tol = 1e-10, cond_ceiling = 1e12) {
  K <- (K + t(K)) / 2; n <- nrow(K)
  ev <- eigen(K, symmetric = TRUE)
  condition <- if (min(ev$values) > 0) max(ev$values) / min(ev$values) else Inf
  ch <- tryCatch(chol(K + diag(jitter, n)), error = function(e) NULL)
  if (!is.null(ch) && condition <= cond_ceiling)
    return(list(lower = t(ch), route = "cholesky", condition = condition, rank = n, error = 0))
  keep <- ev$values > rel_tol * max(ev$values)
  L <- ev$vectors[, keep, drop = FALSE] %*% diag(sqrt(ev$values[keep]), sum(keep))
  list(lower = cbind(L, matrix(0, n, n - sum(keep))), route = "eigen", condition = condition,
       rank = sum(keep), error = max(abs(K - L %*% t(L))))
}

lattice_conditional <- function(xy_sites, xy_lattice, range, jitter = 1e-8, rel_tol = 1e-10) {
  Kss <- sq_exp(xy_sites, xy_sites, range) + diag(jitter, nrow(xy_sites))
  Kls <- sq_exp(xy_lattice, xy_sites, range)
  ch <- tryCatch(chol(Kss), error = function(e) NULL)
  if (!is.null(ch)) {
    A <- t(backsolve(ch, forwardsolve(t(ch), t(Kls)))); route <- "cholesky"
  } else {
    ev <- eigen(Kss, symmetric = TRUE); keep <- ev$values > rel_tol * max(ev$values)
    V <- ev$vectors[, keep, drop = FALSE]
    A <- Kls %*% V %*% diag(1 / ev$values[keep], sum(keep)) %*% t(V); route <- "eigen"
  }
  C <- sq_exp(xy_lattice, xy_lattice, range) - A %*% t(Kls); C <- (C + t(C)) / 2
  ev <- eigen(C, symmetric = TRUE); keep <- ev$values > rel_tol * max(ev$values)
  L <- ev$vectors[, keep, drop = FALSE] %*% diag(sqrt(ev$values[keep]), sum(keep))
  list(A = A, lower = L, rank = sum(keep), min_eigen = min(ev$values),
       error = max(abs(C - L %*% t(L))), route = route)
}

gauss_hermite <- function(n) {
  # Golub-Welsch on the Jacobi matrix of the physicists' Hermite polynomials.
  J <- matrix(0, n, n); off <- sqrt(seq_len(n - 1) / 2)
  J[cbind(1:(n - 1), 2:n)] <- off; J[cbind(2:n, 1:(n - 1))] <- off
  e <- eigen(J, symmetric = TRUE)
  list(nodes = e$values, weights = sqrt(pi) * e$vectors[1, ]^2)
}
```

- [ ] **Step 4: Run to verify it passes**

Run: same command. Expected: `Kernel tests passed.`

- [ ] **Step 5: Commit**

```bash
git add dev/simstudy/spatial-recoverability-pilot/kernel.R dev/simstudy/spatial-recoverability-pilot/test-kernel.R
git commit -m "Add kernel, stable factorisation and quadrature utilities for the recoverability pilot"
```

---

### Task 3: Generator

**Files:**
- Create: `dev/simstudy/spatial-recoverability-pilot/generator.R`
- Create: `dev/simstudy/spatial-recoverability-pilot/test-generator.R`

**Interfaces:**
- Consumes: `kernel.R`, `seeds.R`.
- Produces: `pilot_species()` data frame with columns `species, occupancy, class, log_range, slope`; `pilot_lattice()`; `cluster_separation(k)`; `make_arrangements(n, seed)` returning `list(uniform, clustered, centres, separation, rejections)`; `make_community(truth, n, replicate, range = NA)` returning the community list described in Step 3; `pilot_cells()` data frame of all 12 cells with columns `truth, range, n, cell`.

- [ ] **Step 1: Write the failing tests**

```r
# Run: dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-generator.R
library(testthat)
source("dev/simstudy/spatial-recoverability-pilot/kernel.R")
source("dev/simstudy/spatial-recoverability-pilot/seeds.R")
source("dev/simstudy/spatial-recoverability-pilot/generator.R")

test_that("species table crosses three classes with four occupancies", {
  sp <- pilot_species()
  expect_equal(nrow(sp), 12L)
  expect_equal(unname(table(sp$class, sp$occupancy)), matrix(1L, 3, 4))
  expect_equal(sort(unique(exp(sp$log_range))), c(0.03, 0.08, 0.15))
  expect_equal(sp$slope, rep(c(0.4, -0.4), 6))
})

test_that("arrangements have the declared counts, separations and axis balance", {
  for (n in c(100L, 300L, 900L)) {
    a <- make_arrangements(n, seed = 11L)
    expect_equal(nrow(a$uniform), n); expect_equal(nrow(a$clustered), n)
    expect_true(all(a$uniform >= 0 & a$uniform <= 1)); expect_true(all(a$clustered >= 0 & a$clustered <= 1))
    expect_equal(nrow(a$centres), n / 10L)
    d <- as.matrix(dist(a$centres)); diag(d) <- Inf
    expect_true(min(d) >= a$separation - 1e-12)
    expect_equal(a$separation, 0.2 * sqrt(10 / (n / 10)))
    for (m in list(a$uniform, a$clustered)) expect_lt(abs(sd(m[, 1]) / sd(m[, 2]) - 1), 0.1 + 1e-12)
  }
})

test_that("truth 2 community hits its occupancy targets and spatial SD on the lattice", {
  cm <- make_community("t2", n = 100L, replicate = 1L)
  lat <- cm$index$lattice
  expect_equal(colMeans(cm$psi[lat, ]), cm$species$occupancy, tolerance = 1e-6, ignore_attr = TRUE)
  expect_equal(nrow(cm$points), 100L + 100L + 1600L); expect_false(anyDuplicated(cm$points) > 0)
  expect_true(all(cm$z %in% 0:1))
  expect_equal(cm$spatial_sd, rep(1, 12))
  expect_true(all(cm$conditioning$route %in% c("cholesky", "eigen")))
})

test_that("truth 1 community has the 90/10 split in expectation and unit spatial variance", {
  vars <- replicate(20, { cm <- make_community("t1", n = 100L, replicate = sample.int(1e6, 1), range = 0.08)
    c(total = mean(apply(cm$spatial[cm$index$lattice, ], 2, var)),
      shared = mean(apply(cm$shared_part[cm$index$lattice, ], 2, var))) })
  expect_equal(mean(vars["total", ]), 1, tolerance = 0.2)
  expect_equal(mean(vars["shared", ]) / mean(vars["total", ]), 0.9, tolerance = 0.1)
})

test_that("communities reproduce from their replicate and differ across replicates", {
  a <- make_community("t2", 100L, 3L); b <- make_community("t2", 100L, 3L); c3 <- make_community("t2", 100L, 4L)
  expect_identical(a$z, b$z); expect_false(identical(a$z, c3$z))
})

test_that("the cell table has 9 truth 1 cells and 3 truth 2 cells", {
  cells <- pilot_cells(); expect_equal(nrow(cells), 12L)
  expect_equal(sum(cells$truth == "t1"), 9L); expect_true(all(is.na(cells$range[cells$truth == "t2"])))
})
cat("Generator tests passed.\n")
```

- [ ] **Step 2: Run to verify it fails**

Run: `dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-generator.R`
Expected: error, `generator.R` not found.

- [ ] **Step 3: Implement generator.R**

```r
# Generator for the spatial recoverability pilot. Standalone; no occJSDM code.
pilot_species <- function() {
  data.frame(species = sprintf("species%02d", 1:12),
             occupancy = rep(c(0.05, 0.25, 0.50, 0.75), each = 3),
             class = rep(c("short", "medium", "long"), 4),
             log_range = rep(log(c(0.03, 0.08, 0.15)), 4),
             slope = rep(c(0.4, -0.4), 6), stringsAsFactors = FALSE)
}

pilot_lattice <- function(m = 40L) {
  g <- (seq_len(m) - 0.5) / m
  xy <- as.matrix(expand.grid(x = g, y = g)); dimnames(xy) <- list(NULL, c("x", "y")); xy
}

pilot_cells <- function() {
  t1 <- expand.grid(range = c(0.03, 0.08, 0.15), n = c(100L, 300L, 900L))
  cells <- rbind(data.frame(truth = "t1", range = t1$range, n = t1$n),
                 data.frame(truth = "t2", range = NA_real_, n = c(100L, 300L, 900L)))
  cells$cell <- ifelse(cells$truth == "t1", sprintf("t1-r%.2f-n%d", cells$range, cells$n), sprintf("t2-n%d", cells$n))
  cells
}

cluster_separation <- function(k) 0.2 * sqrt(10 / k)

axis_ratio_ok <- function(xy, max_ratio = 0.1) abs(sd(xy[, 1]) / sd(xy[, 2]) - 1) <= max_ratio

draw_separated_centres <- function(k, min_sep, max_tries = 200000L) {
  centres <- matrix(numeric(0), 0, 2)
  for (i in seq_len(max_tries)) {
    cand <- runif(2)
    ok <- nrow(centres) == 0L || min(sqrt(colSums((t(centres) - cand)^2))) >= min_sep
    if (ok) centres <- rbind(centres, cand)
    if (nrow(centres) == k) return(unname(centres))
  }
  stop("Could not place ", k, " cluster centres at separation ", min_sep)
}

make_arrangements <- function(n, seed, cluster_radius = 0.02, max_ratio = 0.1) {
  stopifnot(n %% 10L == 0L)
  set.seed(seed); rejections <- 0L
  repeat { uniform <- matrix(runif(2 * n), n, 2); if (axis_ratio_ok(uniform, max_ratio)) break; rejections <- rejections + 1L }
  k <- n %/% 10L; separation <- cluster_separation(k)
  repeat {
    centres <- draw_separated_centres(k, separation)
    clustered <- do.call(rbind, lapply(seq_len(k), function(j) {
      r <- cluster_radius * sqrt(runif(10)); a <- runif(10, 0, 2 * pi)
      pmin(pmax(cbind(centres[j, 1] + r * cos(a), centres[j, 2] + r * sin(a)), 0), 1)
    }))
    if (axis_ratio_ok(clustered, max_ratio)) break
    rejections <- rejections + 1L
  }
  name <- function(m) { m <- unname(m); dimnames(m) <- list(NULL, c("x", "y")); m }
  list(uniform = name(uniform), clustered = name(clustered), centres = centres,
       separation = separation, rejections = rejections)
}

solve_intercepts <- function(offset_lattice, targets) {
  vapply(seq_along(targets), function(s) uniroot(function(b) mean(plogis(b + offset_lattice[, s])) - targets[s],
                                                 c(-30, 30), tol = 1e-12)$root, numeric(1))
}

make_community <- function(truth, n, replicate, range = NA_real_, env_range = 0.5, env_noise_sd = 0.3,
                           shared_fraction = 0.9) {
  stopifnot(truth %in% c("t1", "t2"), truth == "t2" || is.finite(range))
  cell <- if (truth == "t1") sprintf("t1-r%.2f-n%d", range, n) else sprintf("t2-n%d", n)
  sp <- pilot_species(); S <- nrow(sp)
  arr <- make_arrangements(n, pilot_seed(pilot_label("sites", cell, paste0("rep", replicate))))
  lattice <- pilot_lattice()
  points <- rbind(arr$uniform, arr$clustered, lattice)
  index <- list(uniform = seq_len(n), clustered = n + seq_len(n), lattice = 2L * n + seq_len(nrow(lattice)))
  stopifnot(!anyDuplicated(points))
  N <- nrow(points)
  set.seed(pilot_seed(pilot_label("data", cell, paste0("rep", replicate))))
  conditioning <- list()
  factor_for <- function(r, label) {
    f <- stable_factor(sq_exp(points, points, r))
    conditioning[[length(conditioning) + 1L]] <<- data.frame(what = label, range = r, route = f$route,
                                                             condition = f$condition, rank = f$rank, error = f$error)
    f$lower
  }
  environment <- as.vector(factor_for(env_range, "environment") %*% rnorm(N)) + rnorm(N, 0, env_noise_sd)
  if (truth == "t1") {
    Lf <- factor_for(range, "shared-and-independent")
    shared_fields <- Lf %*% matrix(rnorm(N * 2), N, 2)
    independent <- Lf %*% matrix(rnorm(N * S), N, S)
    loadings <- matrix(rnorm(2 * S), S, 2)
    loadings <- loadings / sqrt(rowSums(loadings^2))
    shared_part <- sqrt(shared_fraction) * shared_fields %*% t(loadings)
    spatial <- shared_part + sqrt(1 - shared_fraction) * independent
    ranges <- rep(range, S)
  } else {
    ranges <- exp(sp$log_range)
    factors <- lapply(sort(unique(ranges)), function(r) factor_for(r, paste0("class-", r)))
    names(factors) <- as.character(sort(unique(ranges)))
    spatial <- vapply(seq_len(S), function(s) as.vector(factors[[as.character(ranges[s])]] %*% rnorm(N)), numeric(N))
    shared_fields <- NULL; independent <- NULL; loadings <- NULL; shared_part <- NULL
  }
  offset <- environment %o% sp$slope + spatial
  B0 <- solve_intercepts(offset[index$lattice, , drop = FALSE], sp$occupancy)
  psi <- plogis(sweep(offset, 2, B0, "+"))
  z <- matrix(as.integer(runif(N * S) < psi), N, S)
  colnames(spatial) <- colnames(psi) <- colnames(z) <- sp$species
  list(truth = truth, cell = cell, n = n, replicate = replicate, range = range, species = sp,
       points = points, index = index, arrangements = arr, environment = environment,
       spatial = spatial, shared_fields = shared_fields, independent = independent, loadings = loadings,
       shared_part = shared_part, shared_fraction = shared_fraction, ranges = ranges, spatial_sd = rep(1, S),
       psi = psi, z = z, B0 = B0, B = sp$slope, env_range = env_range, env_noise_sd = env_noise_sd,
       conditioning = do.call(rbind, conditioning),
       seeds = c(sites = pilot_seed(pilot_label("sites", cell, paste0("rep", replicate))),
                 data = pilot_seed(pilot_label("data", cell, paste0("rep", replicate)))))
}
```

- [ ] **Step 4: Run to verify it passes**

Run: same command. Expected: `Generator tests passed.` (the truth 1 variance test takes about a minute because it builds 20 communities).

- [ ] **Step 5: Commit**

```bash
git add dev/simstudy/spatial-recoverability-pilot/generator.R dev/simstudy/spatial-recoverability-pilot/test-generator.R
git commit -m "Add the recoverability pilot generator for both spatial truths"
```

---

### Task 4: Joint elliptical-slice sampler and oracle wrappers

**Files:**
- Create: `dev/simstudy/spatial-recoverability-pilot/joint-ellipse.cpp`
- Create: `dev/simstudy/spatial-recoverability-pilot/oracle.R`
- Create: `dev/simstudy/spatial-recoverability-pilot/test-oracle.R`

**Interfaces:**
- Consumes: `kernel.R`, `generator.R`, the sweep's `dev/simstudy/spatial-amplitude-prior/diagnosis/ellipse.cpp` (MD5 `8154a611d1b7ab926eb32bf441e9af3e`, sourced unmodified).
- Produces: `compile_pilot_samplers(repo, cache)`; `joint_ellipse_draws_cpp(lower, Y, offset, weights, n_shared, indep_scale, truth_centred, initial, nburn, niter, thin)` returning `list(draws, traces, mean_proposals)` where `draws` is a cube `d x nkeep x nchain`, `traces` a cube `niter x nchain x 3S` ordered field average, field-truth projection, probability average per species; `run_species_oracle(cm, arrangement, species, nburn, niter, thin, seed)`; `run_joint_oracle(cm, arrangement, nburn, niter, thin, seed)`; both return `list(site_draws, traces, offset, factor_record, mean_proposals)` where `site_draws` is a list of block matrices each `n x (nkeep * nchain)`; `screen_traces(traces)` returning a data frame with `quantity, species, rhat, ess_bulk, ess_tail, pass`.

- [ ] **Step 1: Write the failing tests**

```r
# Run: dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-oracle.R
library(testthat)
for (f in c("kernel.R", "seeds.R", "generator.R", "oracle.R")) source(file.path("dev/simstudy/spatial-recoverability-pilot", f))
compile_pilot_samplers(".")

test_that("the sweep sampler source is unmodified", {
  expect_identical(unname(tools::md5sum("dev/simstudy/spatial-amplitude-prior/diagnosis/ellipse.cpp")),
                   "8154a611d1b7ab926eb32bf441e9af3e")
})

test_that("the joint sampler reduces exactly to the single-field sampler", {
  set.seed(5); n <- 30L; xy <- matrix(runif(2 * n), n, 2)
  L <- stable_factor(sq_exp(xy, xy, 0.1))$lower
  y <- as.numeric(runif(n) < 0.5); offset <- rnorm(n)
  initial <- cbind(rep(0, n), L %*% matrix(rnorm(n * 3), n, 3))
  set.seed(99); a <- ellipse_draws_cpp(L, y, 1L, offset, FALSE, initial, 20L, 50L)
  set.seed(99); b <- joint_ellipse_draws_cpp(L, matrix(y, n, 1), matrix(offset, n, 1), matrix(1, 1, 1),
                                             1L, 0, matrix(0, n, 1), initial, 20L, 50L, 1L)
  expect_equal(a$draws, b$draws, tolerance = 1e-12)
  expect_equal(a$mean_proposals, b$mean_proposals)
})

test_that("the joint sampler recovers its prior when there is no data", {
  set.seed(6); n <- 20L; S <- 3L; xy <- matrix(runif(2 * n), n, 2)
  L <- stable_factor(sq_exp(xy, xy, 0.1))$lower
  w <- matrix(rnorm(2 * S), S, 2); w <- sqrt(0.9) * w / sqrt(rowSums(w^2))
  Y <- matrix(0, n, S); offset <- matrix(-50, n, S)   # likelihood flat in the state
  d <- (2L + S) * n
  initial <- matrix(rnorm(d * 4), d, 4)
  fit <- joint_ellipse_draws_cpp(L, Y, offset, w, 2L, sqrt(0.1), matrix(0, n, S), initial, 200L, 4000L, 1L)
  draws <- matrix(fit$draws, d, 4000 * 4)
  block1 <- draws[1:n, ]
  expect_equal(cov(t(block1)), L %*% t(L), tolerance = 0.15)
})

test_that("the joint sampler matches quadrature on a one-site logistic-normal case", {
  L <- matrix(1.2, 1, 1); y <- matrix(1, 1, 1); offset <- matrix(-0.5, 1, 1)
  set.seed(1); initial <- cbind(0, 1.2 * matrix(rnorm(3), 1, 3))
  fit <- joint_ellipse_draws_cpp(L, y, offset, matrix(1, 1, 1), 1L, 0, matrix(0, 1, 1), initial, 2000L, 20000L, 1L)
  num <- integrate(function(f) f * plogis(-0.5 + f) * dnorm(f, 0, 1.2), -Inf, Inf)$value
  den <- integrate(function(f) plogis(-0.5 + f) * dnorm(f, 0, 1.2), -Inf, Inf)$value
  expect_equal(mean(fit$draws), num / den, tolerance = 0.02)
})

test_that("joint draws reproduce under a fixed seed and thinning keeps every thin-th draw", {
  cm <- make_community("t1", 100L, 1L, range = 0.08)
  a <- run_joint_oracle(cm, "uniform", nburn = 10L, niter = 40L, thin = 4L, seed = 3L)
  b <- run_joint_oracle(cm, "uniform", nburn = 10L, niter = 40L, thin = 4L, seed = 3L)
  expect_identical(a$site_draws, b$site_draws)
  expect_length(a$site_draws, 14L); expect_equal(dim(a$site_draws[[1]]), c(100L, 10L * 4L))
  expect_equal(dim(a$traces), c(40L, 4L, 36L))
})

test_that("species oracle returns the agreed shapes and a screen table", {
  cm <- make_community("t2", 100L, 1L)
  r <- run_species_oracle(cm, "clustered", species = 5L, nburn = 10L, niter = 40L, thin = 4L, seed = 2L)
  expect_length(r$site_draws, 1L); expect_equal(dim(r$site_draws[[1]]), c(100L, 40L))
  expect_equal(dim(r$traces), c(40L, 4L, 3L))
  sc <- screen_traces(r$traces)
  expect_equal(nrow(sc), 3L); expect_true(all(c("rhat", "ess_bulk", "ess_tail", "pass") %in% names(sc)))
})
cat("Oracle tests passed.\n")
```

- [ ] **Step 2: Run to verify it fails**

Run: `dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-oracle.R`
Expected: error, `oracle.R` not found.

- [ ] **Step 3: Write joint-ellipse.cpp**

```cpp
// Joint elliptical-slice sampler for the recoverability pilot's truth 1 oracle.
// State = [f_1 .. f_K (shared fields), g_1 .. g_S (independent parts)], each block n long.
// eta_s = offset_s + sum_k weights(s,k) f_k + indep_scale * g_s. Independent blocks are
// present only when indep_scale > 0. With K = 1, S = 1, indep_scale = 0 the random-number
// sequence and arithmetic match ellipse.cpp exactly.
// Murray, Adams and MacKay (2010), Algorithm 1, full 2*pi bracket.
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

static arma::mat species_effects(const arma::vec& state, const arma::mat& weights, int n_shared,
                                 double indep_scale, arma::uword n) {
  const arma::uword S = weights.n_rows;
  arma::mat h(n, S, arma::fill::zeros);
  for (arma::uword s = 0; s < S; ++s) {
    for (int k = 0; k < n_shared; ++k) h.col(s) += weights(s, k) * state.subvec(k * n, (k + 1) * n - 1);
    if (indep_scale > 0) h.col(s) += indep_scale * state.subvec((n_shared + s) * n, (n_shared + s + 1) * n - 1);
  }
  return h;
}

static double joint_loglik(const arma::vec& state, const arma::mat& Y, const arma::mat& offset,
                           const arma::mat& weights, int n_shared, double indep_scale) {
  const arma::uword n = Y.n_rows, S = Y.n_cols;
  arma::mat h = species_effects(state, weights, n_shared, indep_scale, n);
  double value = 0;
  for (arma::uword s = 0; s < S; ++s)
    for (arma::uword i = 0; i < n; ++i) {
      double eta = offset(i, s) + h(i, s);
      double softplus = std::max(eta, 0.0) + std::log1p(std::exp(-std::abs(eta)));
      value += Y(i, s) * eta - softplus;
    }
  return value;
}

// [[Rcpp::export]]
Rcpp::List joint_ellipse_draws_cpp(arma::mat lower, arma::mat Y, arma::mat offset, arma::mat weights,
                                   int n_shared, double indep_scale, arma::mat truth_centred,
                                   arma::mat initial, int nburn, int niter, int thin) {
  const arma::uword n = Y.n_rows, S = Y.n_cols;
  const int n_indep = indep_scale > 0 ? (int) S : 0;
  const arma::uword blocks = n_shared + n_indep, d = blocks * n;
  if (lower.n_rows != n || lower.n_cols != n || offset.n_rows != n || offset.n_cols != S ||
      weights.n_rows != S || weights.n_cols != (arma::uword) n_shared || truth_centred.n_rows != n ||
      truth_centred.n_cols != S || initial.n_rows != d || initial.n_cols < 2 || nburn < 0 || niter < 4 ||
      thin < 1 || n_shared < 1 || !lower.is_finite() || !initial.is_finite() || !offset.is_finite() ||
      arma::any(arma::vectorise(Y) < 0) || arma::any(arma::vectorise(Y) > 1))
    Rcpp::stop("Invalid joint sampler dimensions or values");
  const int nkeep = (niter + thin - 1) / thin;
  arma::cube draws(d, nkeep, initial.n_cols);
  arma::cube traces(niter, initial.n_cols, 3 * S);
  arma::vec proposals(initial.n_cols, arma::fill::zeros);
  arma::vec truth_ss(S);
  for (arma::uword s = 0; s < S; ++s) truth_ss[s] = arma::dot(truth_centred.col(s), truth_centred.col(s));
  for (arma::uword chain = 0; chain < initial.n_cols; ++chain) {
    arma::vec state = initial.col(chain);
    double loglik = joint_loglik(state, Y, offset, weights, n_shared, indep_scale);
    for (int iteration = 0; iteration < nburn + niter; ++iteration) {
      if (iteration % 100 == 0) Rcpp::checkUserInterrupt();
      arma::vec white(d);
      for (arma::uword j = 0; j < d; ++j) white[j] = R::rnorm(0, 1);
      arma::vec direction(d);
      for (arma::uword b = 0; b < blocks; ++b)
        direction.subvec(b * n, (b + 1) * n - 1) = lower * white.subvec(b * n, (b + 1) * n - 1);
      double threshold = loglik + std::log(R::runif(0, 1));
      double theta = R::runif(0, 2 * M_PI);
      double left = theta - 2 * M_PI, right = theta;
      int attempts = 0;
      while (true) {
        if (++attempts > 10000) Rcpp::stop("Elliptical bracket failed to contract");
        arma::vec candidate = state * std::cos(theta) + direction * std::sin(theta);
        double candidate_loglik = joint_loglik(candidate, Y, offset, weights, n_shared, indep_scale);
        if (candidate_loglik >= threshold) { state = candidate; loglik = candidate_loglik; break; }
        if (theta < 0) left = theta; else right = theta;
        theta = R::runif(left, right);
      }
      proposals[chain] += attempts;
      if (iteration >= nburn) {
        const int kept = iteration - nburn;
        arma::mat h = species_effects(state, weights, n_shared, indep_scale, n);
        for (arma::uword s = 0; s < S; ++s) {
          traces(kept, chain, 3 * s) = arma::mean(h.col(s));
          traces(kept, chain, 3 * s + 1) = truth_ss[s] > 0 ? arma::dot(h.col(s), truth_centred.col(s)) / truth_ss[s] : 0.0;
          double p = 0; for (arma::uword i = 0; i < n; ++i) p += 1.0 / (1.0 + std::exp(-(offset(i, s) + h(i, s))));
          traces(kept, chain, 3 * s + 2) = p / n;
        }
        if (kept % thin == 0) draws.slice(chain).col(kept / thin) = state;
      }
    }
  }
  return Rcpp::List::create(Rcpp::Named("draws") = draws, Rcpp::Named("traces") = traces,
                            Rcpp::Named("mean_proposals") = proposals / (nburn + niter));
}
```

- [ ] **Step 4: Write oracle.R**

```r
# Oracle wrappers for the recoverability pilot. Truth 2 uses the sweep's ellipse.cpp
# unmodified; truth 1 uses joint-ellipse.cpp. Both know every parameter but the fields.
compile_pilot_samplers <- function(repo = ".", cache = NULL) {
  if (is.null(cache)) cache <- file.path(tempdir(), "recoverability-pilot-cpp")
  dir.create(cache, showWarnings = FALSE, recursive = TRUE)
  Rcpp::sourceCpp(file.path(repo, "dev/simstudy/spatial-amplitude-prior/diagnosis/ellipse.cpp"), cacheDir = cache, env = globalenv())
  Rcpp::sourceCpp(file.path(repo, "dev/simstudy/spatial-recoverability-pilot/joint-ellipse.cpp"), cacheDir = cache, env = globalenv())
  invisible(tools::md5sum(c(file.path(repo, "dev/simstudy/spatial-amplitude-prior/diagnosis/ellipse.cpp"),
                            file.path(repo, "dev/simstudy/spatial-recoverability-pilot/joint-ellipse.cpp"))))
}

training_offset <- function(cm, idx) sweep(cm$environment[idx] %o% cm$B, 2, cm$B0, "+")

factor_record <- function(f, range, what) data.frame(what = what, range = range, route = f$route,
                                                     condition = f$condition, rank = f$rank, error = f$error)

trace_matrix <- function(draws_full, offset, truth_centred) {
  # draws_full: n x niter x nchain for one species. Returns niter x nchain x 3.
  n <- dim(draws_full)[1]; ni <- dim(draws_full)[2]; nc <- dim(draws_full)[3]
  out <- array(NA_real_, c(ni, nc, 3))
  tss <- sum(truth_centred^2)
  for (ch in seq_len(nc)) {
    m <- draws_full[, , ch]
    out[, ch, 1] <- colMeans(m)
    out[, ch, 2] <- if (tss > 0) colSums(m * truth_centred) / tss else 0
    out[, ch, 3] <- colMeans(plogis(offset + m))
  }
  out
}

run_species_oracle <- function(cm, arrangement, species, nburn = 3000L, niter = 4000L, thin = 10L, seed, nchain = 4L) {
  idx <- cm$index[[arrangement]]; xy <- cm$points[idx, , drop = FALSE]; n <- length(idx)
  f <- stable_factor(sq_exp(xy, xy, cm$ranges[species]))
  offset <- training_offset(cm, idx)[, species]
  set.seed(seed)
  initial <- cbind(rep(0, n), f$lower %*% matrix(rnorm(n * (nchain - 1L)), n, nchain - 1L))
  fit <- ellipse_draws_cpp(f$lower, cm$z[idx, species], 1L, offset, FALSE, initial, nburn, niter)
  truth <- cm$spatial[idx, species]; tc <- truth - mean(truth)
  traces <- trace_matrix(fit$draws, offset, tc)
  keep <- seq(1L, niter, by = thin)
  site_draws <- list(matrix(fit$draws[, keep, ], n, length(keep) * nchain))
  list(site_draws = site_draws, traces = traces, offset = offset,
       factor_record = factor_record(f, cm$ranges[species], paste0("prior-", arrangement, "-", cm$species$species[species])),
       mean_proposals = fit$mean_proposals, nburn = nburn, niter = niter, thin = thin, seed = seed)
}

joint_weights <- function(cm) sqrt(cm$shared_fraction) * cm$loadings

run_joint_oracle <- function(cm, arrangement, nburn = 3000L, niter = 4000L, thin = 10L, seed, nchain = 4L) {
  stopifnot(cm$truth == "t1")
  idx <- cm$index[[arrangement]]; xy <- cm$points[idx, , drop = FALSE]; n <- length(idx); S <- nrow(cm$species)
  f <- stable_factor(sq_exp(xy, xy, cm$range))
  offset <- training_offset(cm, idx)
  w <- joint_weights(cm); c_indep <- sqrt(1 - cm$shared_fraction)
  truth <- cm$spatial[idx, , drop = FALSE]; tc <- sweep(truth, 2, colMeans(truth), "-")
  blocks <- 2L + S; d <- blocks * n
  set.seed(seed)
  initial <- cbind(rep(0, d), do.call(rbind, replicate(blocks, f$lower %*% matrix(rnorm(n * (nchain - 1L)), n, nchain - 1L), simplify = FALSE)))
  fit <- joint_ellipse_draws_cpp(f$lower, cm$z[idx, , drop = FALSE], offset, w, 2L, c_indep, tc, initial, nburn, niter, thin)
  nkeep <- dim(fit$draws)[2]
  site_draws <- lapply(seq_len(blocks), function(b) matrix(fit$draws[(b - 1L) * n + seq_len(n), , ], n, nkeep * nchain))
  list(site_draws = site_draws, traces = fit$traces, offset = offset,
       factor_record = factor_record(f, cm$range, paste0("prior-", arrangement, "-joint")),
       mean_proposals = fit$mean_proposals, nburn = nburn, niter = niter, thin = thin, seed = seed)
}

screen_traces <- function(traces, rhat_max = 1.01, ess_min = 400) {
  q <- dim(traces)[3]; S <- q %/% 3L
  labels <- c("field_average", "field_truth_projection", "probability_average")
  rows <- lapply(seq_len(q), function(j) {
    m <- traces[, , j]
    constant <- any(apply(m, 2, function(z) all(z == z[1])))
    rhat <- if (constant) NA_real_ else suppressWarnings(posterior::rhat(m))
    bulk <- if (constant) NA_real_ else suppressWarnings(posterior::ess_bulk(m))
    tail <- if (constant) NA_real_ else suppressWarnings(posterior::ess_tail(m))
    data.frame(quantity = labels[(j - 1L) %% 3L + 1L], species = (j - 1L) %/% 3L + 1L,
               rhat = rhat, ess_bulk = bulk, ess_tail = tail,
               pass = is.finite(rhat) && is.finite(bulk) && is.finite(tail) && rhat <= rhat_max && bulk >= ess_min && tail >= ess_min)
  })
  do.call(rbind, rows)
}
```

- [ ] **Step 5: Run to verify it passes**

Run: same command. Expected: `Oracle tests passed.` The reduction test must pass at tolerance 1e-12; if it does not, the random-number order in the C++ differs from `ellipse.cpp` and must be fixed, not the tolerance.

- [ ] **Step 6: Commit**

```bash
git add dev/simstudy/spatial-recoverability-pilot/joint-ellipse.cpp dev/simstudy/spatial-recoverability-pilot/oracle.R dev/simstudy/spatial-recoverability-pilot/test-oracle.R
git commit -m "Add the joint elliptical-slice oracle and wrappers for the recoverability pilot"
```

---

### Task 5: Lattice prediction, uninformed reference and per-species scores

**Files:**
- Create: `dev/simstudy/spatial-recoverability-pilot/predict.R`
- Create: `dev/simstudy/spatial-recoverability-pilot/score.R`
- Create: `dev/simstudy/spatial-recoverability-pilot/test-score.R`

**Interfaces:**
- Consumes: `kernel.R`, `generator.R`, `oracle.R` outputs (`site_draws`, `offset`).
- Produces: `uninformed_reference(eta, sd = 1, gh = gauss_hermite(40L))`; `predict_lattice(cm, arrangement, oracle, seed)` returning `list(probability_mean (1600 x S), field_mean (1600 x S), site_field_mean (n x S), site_probability_mean (n x S), conditional_records (data frame))`; `score_species(cm, arrangement, prediction, screen)` returning one row per species with the columns listed in Step 3; `aggregate_scores(scores)`; `gate_cells(aggregated, thresholds)`; `read_map(gate)`.

- [ ] **Step 1: Write the failing tests**

```r
# Run: dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-score.R
library(testthat)
for (f in c("kernel.R", "seeds.R", "generator.R", "oracle.R", "predict.R", "score.R")) source(file.path("dev/simstudy/spatial-recoverability-pilot", f))
compile_pilot_samplers(".")

test_that("uninformed reference agrees with Monte Carlo integration", {
  set.seed(1); eta <- c(-2, 0, 1.5)
  mc <- vapply(eta, function(e) mean(plogis(e + rnorm(2e6))), numeric(1))
  expect_equal(uninformed_reference(eta), mc, tolerance = 2e-3)
  expect_true(all(abs(uninformed_reference(eta) - plogis(eta)) > 1e-3))   # not the biased plug-in
})

test_that("score_species reproduces known correlation, slope and RMSE on a synthetic case", {
  cm <- make_community("t2", 100L, 2L); lat <- cm$index$lattice; S <- 12L
  truth <- cm$spatial[lat, ]
  # Fake prediction: field = 0.5 * truth, probability = true psi exactly.
  prediction <- list(field_mean = 0.5 * truth, probability_mean = cm$psi[lat, ],
                     site_field_mean = 0.5 * cm$spatial[cm$index$uniform, ], site_probability_mean = cm$psi[cm$index$uniform, ],
                     conditional_records = data.frame())
  screen <- data.frame(species = 1:S, status = "pass")
  sc <- score_species(cm, "uniform", prediction, screen)
  expect_equal(nrow(sc), S)
  expect_equal(sc$field_slope, rep(0.5, S), tolerance = 1e-10)
  expect_equal(sc$field_correlation, rep(1, S), tolerance = 1e-10)
  expect_equal(sc$probability_rmse_oracle, rep(0, S), tolerance = 1e-12)
  expect_equal(sc$contrast_correlation, rep(1, S), tolerance = 1e-10)
  expect_equal(sc$headroom_fraction, rep(1, S), tolerance = 1e-10)
  expect_true(all(sc$status == "pass"))
})

test_that("predict_lattice returns probabilities in (0,1) and records the conditional", {
  cm <- make_community("t2", 100L, 1L)
  o <- run_species_oracle(cm, "uniform", 6L, nburn = 20L, niter = 80L, thin = 4L, seed = 4L)
  p <- predict_lattice(cm, "uniform", list(species = 6L, oracles = list(o)), seed = 9L)
  expect_equal(dim(p$probability_mean), c(1600L, 1L)); expect_true(all(p$probability_mean > 0 & p$probability_mean < 1))
  expect_lt(max(p$conditional_records$error), 1e-5)
})

test_that("gate logic counts unresolved cells as unresolved and applies every threshold", {
  th <- default_thresholds()
  make <- function(pattern, slope, headroom, status = "pass") data.frame(
    truth = "t2", cell = "t2-n300", class = "medium", n = 300L, range = NA, replicate = 1:10, arrangement = "uniform",
    contrast_correlation = pattern, field_slope = slope, headroom_pp = headroom, status = status)
  pass <- gate_cells(make(rep(0.7, 10), rep(0.8, 10), rnorm(10, 3, 0.2)), th)
  expect_identical(pass$decision, "pass")
  weak <- gate_cells(make(c(rep(0.7, 7), rep(0.3, 3)), rep(0.8, 10), rnorm(10, 3, 0.2)), th)
  expect_identical(weak$decision, "fail")      # only 7 of 10 at least 0.4
  flat <- gate_cells(make(rep(0.7, 10), rep(0.3, 10), rnorm(10, 3, 0.2)), th)
  expect_identical(flat$decision, "fail")      # slope below 0.5
  noroom <- gate_cells(make(rep(0.7, 10), rep(0.8, 10), rnorm(10, 0.2, 0.5)), th)
  expect_identical(noroom$decision, "fail")    # interval includes zero
  unres <- gate_cells(make(rep(0.7, 10), rep(0.8, 10), rnorm(10, 3, 0.2), status = c(rep("pass", 9), "unresolved")), th)
  expect_identical(unres$decision, "unresolved")
})

test_that("read_map covers all four branches", {
  g <- function(t1_pass_n, t2_pass_n) {
    rows <- rbind(
      expand.grid(truth = "t1", range = c(0.03, 0.08, 0.15), class = NA, n = c(100L, 300L, 900L), stringsAsFactors = FALSE),
      expand.grid(truth = "t2", range = NA, class = c("short", "medium", "long"), n = c(100L, 300L, 900L), stringsAsFactors = FALSE))
    rows$decision <- "fail"
    rows$decision[rows$truth == "t1" & rows$range == 0.08 & rows$n %in% t1_pass_n] <- "pass"
    rows$decision[rows$truth == "t2" & rows$class == "medium" & rows$n %in% t2_pass_n] <- "pass"
    rows
  }
  expect_identical(read_map(g(c(300L, 900L), 900L))$branch, "both")
  expect_identical(read_map(g(c(300L, 900L), 900L))$recommended_n, 900L)
  expect_identical(read_map(g(c(300L, 900L), integer(0)))$branch, "truth1_only")
  expect_identical(read_map(g(c(300L, 900L), integer(0)))$recommended_n, 300L)
  expect_identical(read_map(g(integer(0), 900L))$branch, "truth2_only")
  expect_identical(read_map(g(integer(0), integer(0)))$branch, "stop")
})
cat("Score tests passed.\n")
```

- [ ] **Step 2: Run to verify it fails**

Run: `dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-score.R`
Expected: error, `predict.R` not found.

- [ ] **Step 3: Write predict.R**

```r
# Lattice prediction by Gaussian-process conditioning and the spatially uninformed reference.
uninformed_reference <- function(eta, sd = 1, gh = gauss_hermite(40L)) {
  vapply(eta, function(e) sum(gh$weights * plogis(e + sqrt(2) * sd * gh$nodes)) / sqrt(pi), numeric(1))
}

# `oracle` is list(species = integer vector, oracles = list of run_*_oracle results) for truth 2
# (one oracle per species) or list(joint = run_joint_oracle result) for truth 1.
predict_lattice <- function(cm, arrangement, oracle, seed) {
  idx <- cm$index[[arrangement]]; xy <- cm$points[idx, , drop = FALSE]; lat <- cm$points[cm$index$lattice, , drop = FALSE]
  offset_lat <- training_offset(cm, cm$index$lattice); n <- length(idx); S <- nrow(cm$species)
  set.seed(seed)
  records <- list()
  conditional_for <- function(range, what) {
    cf <- lattice_conditional(xy, lat, range)
    if (cf$error > 1e-5) stop("Lattice conditional reconstruction error ", cf$error, " exceeds 1e-5 for ", what)
    records[[length(records) + 1L]] <<- data.frame(what = what, range = range, route = cf$route, rank = cf$rank,
                                                   min_eigen = cf$min_eigen, error = cf$error)
    cf
  }
  draw_lattice <- function(cf, site_block) {
    # site_block: n x ndraw. Returns 1600 x ndraw conditional draws.
    nd <- ncol(site_block)
    cf$A %*% site_block + cf$lower %*% matrix(rnorm(ncol(cf$lower) * nd), ncol(cf$lower), nd)
  }
  if (cm$truth == "t1") {
    o <- oracle$joint; cf <- conditional_for(cm$range, paste0("conditional-", arrangement, "-joint"))
    blocks <- lapply(o$site_draws, function(b) draw_lattice(cf, b))
    w <- joint_weights(cm); c_indep <- sqrt(1 - cm$shared_fraction)
    species_lattice <- function(s) w[s, 1] * blocks[[1]] + w[s, 2] * blocks[[2]] + c_indep * blocks[[2L + s]]
    species_site <- function(s) w[s, 1] * o$site_draws[[1]] + w[s, 2] * o$site_draws[[2]] + c_indep * o$site_draws[[2L + s]]
    species_index <- seq_len(S); offsets_site <- o$offset
  } else {
    species_index <- oracle$species
    species_lattice <- function(j) { o <- oracle$oracles[[j]]; r <- cm$ranges[species_index[j]]
      key <- paste0("conditional-", arrangement, "-r", r)
      cf <- if (is.null(cache[[key]])) { cache[[key]] <<- conditional_for(r, key); cache[[key]] } else cache[[key]]
      draw_lattice(cf, o$site_draws[[1]]) }
    cache <- list()
    species_site <- function(j) oracle$oracles[[j]]$site_draws[[1]]
    offsets_site <- do.call(cbind, lapply(oracle$oracles, function(o) o$offset))
  }
  m <- length(species_index)
  probability_mean <- matrix(NA_real_, nrow(lat), m); field_mean <- probability_mean
  site_field_mean <- matrix(NA_real_, n, m); site_probability_mean <- site_field_mean
  for (j in seq_len(m)) {
    s <- species_index[j]
    H <- species_lattice(j)
    probability_mean[, j] <- rowMeans(plogis(offset_lat[, s] + H)); field_mean[, j] <- rowMeans(H)
    Hs <- species_site(j)
    site_field_mean[, j] <- rowMeans(Hs); site_probability_mean[, j] <- rowMeans(plogis(offsets_site[, j] + Hs))
  }
  colnames(probability_mean) <- colnames(field_mean) <- colnames(site_field_mean) <- colnames(site_probability_mean) <- cm$species$species[species_index]
  list(species = species_index, probability_mean = probability_mean, field_mean = field_mean,
       site_field_mean = site_field_mean, site_probability_mean = site_probability_mean,
       conditional_records = do.call(rbind, records))
}
```

- [ ] **Step 4: Write score.R**

```r
# Per-species scores, aggregation, gate and map reading for the recoverability pilot.
default_thresholds <- function() list(pattern_median = 0.6, pattern_floor = 0.4, pattern_count = 8L,
                                      slope_low = 0.5, slope_high = 1.5, headroom_pp = 1.0, communities = 10L)

score_one <- function(truth_field, pred_field, truth_psi, pred_psi, uninformed) {
  tc <- truth_field - mean(truth_field); fc <- pred_field - mean(pred_field)
  true_contrast <- truth_psi - uninformed; pred_contrast <- pred_psi - uninformed
  rmse <- function(a, b) sqrt(mean((a - b)^2))
  rmse_unif <- rmse(uninformed, truth_psi); rmse_or <- rmse(pred_psi, truth_psi)
  c(contrast_correlation = cor(pred_contrast, true_contrast), contrast_rms = sqrt(mean(pred_contrast^2)),
    true_contrast_rms = sqrt(mean(true_contrast^2)),
    field_correlation = cor(fc, tc), field_slope = sum(fc * tc) / sum(tc^2), field_rmse = rmse(fc, tc),
    zero_field_rmse = sqrt(mean(tc^2)),
    probability_rmse_oracle = rmse_or, probability_rmse_uninformed = rmse_unif,
    probability_mae_oracle = mean(abs(pred_psi - truth_psi)), probability_mae_uninformed = mean(abs(uninformed - truth_psi)),
    headroom_pp = 100 * (rmse_unif - rmse_or), headroom_fraction = if (rmse_unif > 0) (rmse_unif - rmse_or) / rmse_unif else NA_real_)
}

score_species <- function(cm, arrangement, prediction, screen) {
  lat <- cm$index$lattice; idx <- cm$index[[arrangement]]
  eta_lat <- training_offset(cm, lat); eta_site <- training_offset(cm, idx)
  sp_idx <- if (!is.null(prediction$species)) prediction$species else seq_len(nrow(cm$species))
  rows <- lapply(seq_along(sp_idx), function(j) {
    s <- sp_idx[j]
    held <- score_one(cm$spatial[lat, s], prediction$field_mean[, j], cm$psi[lat, s], prediction$probability_mean[, j],
                      uninformed_reference(eta_lat[, s], cm$spatial_sd[s]))
    train <- score_one(cm$spatial[idx, s], prediction$site_field_mean[, j], cm$psi[idx, s], prediction$site_probability_mean[, j],
                       uninformed_reference(eta_site[, s], cm$spatial_sd[s]))
    names(train) <- paste0("train_", names(train))
    data.frame(truth = cm$truth, cell = cm$cell, n = cm$n, range = cm$range, replicate = cm$replicate,
               arrangement = arrangement, species = cm$species$species[s], occupancy = cm$species$occupancy[s],
               class = if (cm$truth == "t2") cm$species$class[s] else NA_character_,
               status = screen$status[screen$species == s][1], as.list(held), as.list(train), stringsAsFactors = FALSE)
  })
  do.call(rbind, rows)
}

# Average species within occupancy level, then pool the 25, 50 and 75% levels, per community.
aggregate_scores <- function(scores) {
  gate_levels <- c(0.25, 0.50, 0.75)
  s <- scores[scores$occupancy %in% gate_levels, ]
  s$group <- if (all(s$truth == "t2")) s$class else "all"
  key <- c("truth", "cell", "n", "range", "replicate", "arrangement", "group")
  metrics <- c("contrast_correlation", "field_slope", "headroom_pp", "headroom_fraction", "field_correlation")
  by_level <- aggregate(s[metrics], s[c(key, "occupancy")], mean)
  out <- aggregate(by_level[metrics], by_level[key], mean)
  status <- aggregate(list(status = s$status), s[key], function(x) if (any(x == "unresolved")) "unresolved" else "pass")
  merge(out, status, by = key)
}

gate_cells <- function(agg, thresholds) {
  # agg: one community per row for a single cell (and class under truth 2), uniform arrangement.
  stopifnot(nrow(agg) >= 2L)
  if (any(agg$status == "unresolved")) return(data.frame(decision = "unresolved", pattern_median = NA, pattern_count = NA,
                                                           slope_median = NA, headroom_mean = NA, headroom_low = NA, headroom_high = NA))
  pm <- median(agg$contrast_correlation); pc <- sum(agg$contrast_correlation >= thresholds$pattern_floor)
  sm <- median(agg$field_slope)
  hm <- mean(agg$headroom_pp); se <- sd(agg$headroom_pp) / sqrt(nrow(agg)); tq <- qt(0.975, nrow(agg) - 1L)
  lo <- hm - tq * se; hi <- hm + tq * se
  ok <- pm >= thresholds$pattern_median && pc >= thresholds$pattern_count &&
    sm >= thresholds$slope_low && sm <= thresholds$slope_high && hm >= thresholds$headroom_pp && lo > 0
  data.frame(decision = if (ok) "pass" else "fail", pattern_median = pm, pattern_count = pc, slope_median = sm,
             headroom_mean = hm, headroom_low = lo, headroom_high = hi)
}

read_map <- function(map) {
  t1 <- map[map$truth == "t1" & map$range == 0.08 & map$decision == "pass", ]
  t2 <- map[map$truth == "t2" & map$class == "medium" & map$decision == "pass", ]
  t1_n <- if (nrow(t1)) min(t1$n) else NA_integer_; t2_n <- if (nrow(t2)) min(t2$n) else NA_integer_
  if (!is.na(t1_n) && !is.na(t2_n)) list(branch = "both", recommended_n = max(t1_n, t2_n), truth1_n = t1_n, truth2_n = t2_n)
  else if (!is.na(t1_n)) list(branch = "truth1_only", recommended_n = t1_n, truth1_n = t1_n, truth2_n = NA_integer_)
  else if (!is.na(t2_n)) list(branch = "truth2_only", recommended_n = t2_n, truth1_n = NA_integer_, truth2_n = t2_n)
  else list(branch = "stop", recommended_n = NA_integer_, truth1_n = NA_integer_, truth2_n = NA_integer_)
}
```

- [ ] **Step 5: Run to verify it passes**

Run: same command. Expected: `Score tests passed.`

- [ ] **Step 6: Commit**

```bash
git add dev/simstudy/spatial-recoverability-pilot/predict.R dev/simstudy/spatial-recoverability-pilot/score.R dev/simstudy/spatial-recoverability-pilot/test-score.R
git commit -m "Add lattice prediction, uninformed reference, scores and gate for the recoverability pilot"
```

---

### Task 6: Driver with prepare, timing, fit and score modes

**Files:**
- Create: `dev/simstudy/spatial-recoverability-pilot/run.R`
- Create: `dev/simstudy/spatial-recoverability-pilot/test-run.R`

**Interfaces:**
- Consumes: every earlier file.
- Produces: command-line modes `prepare`, `timing`, `fit`, `score`; a study directory layout `STUDY/config.rds`, `STUDY/inputs/<cell>-rep<k>.rds`, `STUDY/jobs/<cell>-rep<k>-<arrangement>.rds` (symlink to the volume), `STUDY/status/<job>.json`, `STUDY/scores/species.csv`, `STUDY/scores/aggregate.csv`, `STUDY/scores/gate.csv`, `STUDY/scores/map.json`, `STUDY/timing.json`; `job_table(cells, communities)`; `run_job(study, job, budget)`.

- [ ] **Step 1: Write the failing smoke test**

```r
# Run: dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-run.R
library(testthat)
source("dev/simstudy/spatial-recoverability-pilot/run.R")
study <- file.path(tempdir(), paste0("pilot-smoke-", as.integer(Sys.time()))); dir.create(study)

test_that("prepare writes config before inputs and inputs have hashes", {
  pilot_prepare(study, communities = 2L, cells = pilot_cells()[pilot_cells()$n == 100L, ], jobs_dir = file.path(study, "jobs-local"))
  cfg <- readRDS(file.path(study, "config.rds"))
  expect_true(file.mtime(file.path(study, "config.rds")) <= min(file.mtime(list.files(file.path(study, "inputs"), full.names = TRUE))))
  expect_equal(cfg$thresholds, default_thresholds())
  expect_equal(length(list.files(file.path(study, "inputs"))), 4L * 2L)
  man <- read.csv(file.path(study, "input-manifest.csv"))
  expect_true(all(nchar(man$md5) == 32L))
})

test_that("fit runs a tiny budget, reuses existing jobs and scores", {
  jobs <- job_table(study)
  expect_equal(nrow(jobs), 4L * 2L * 2L)
  pilot_fit(study, truth = "t2", workers = 2L, budget = list(nburn = 10L, niter = 40L, thin = 4L, ext_nburn = 20L, ext_niter = 80L))
  done <- list.files(file.path(study, "jobs-local"), pattern = "^t2-")
  expect_equal(length(done), 2L * 2L)
  mt <- file.mtime(file.path(study, "jobs-local", done))
  pilot_fit(study, truth = "t2", workers = 2L, budget = list(nburn = 10L, niter = 40L, thin = 4L, ext_nburn = 20L, ext_niter = 80L))
  expect_identical(file.mtime(file.path(study, "jobs-local", done)), mt)   # reused, not refitted
  pilot_fit(study, truth = "t1", workers = 2L, budget = list(nburn = 10L, niter = 40L, thin = 4L, ext_nburn = 20L, ext_niter = 80L))
  pilot_score(study)
  sp <- read.csv(file.path(study, "scores", "species.csv"))
  expect_equal(nrow(sp), 16L * 12L)
  expect_true(all(sp$status %in% c("pass", "unresolved")))
  expect_true(file.exists(file.path(study, "scores", "map.json")))
})
cat("Run tests passed.\n")
```

- [ ] **Step 2: Run to verify it fails**

Run: `dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-run.R`
Expected: error, `run.R` not found.

- [ ] **Step 3: Write run.R**

```r
# Driver for the spatial recoverability pilot.
# Usage (always through run-r):
#   run-r -f run.R --args --mode=prepare --study=DIR --jobs=VOLUME_DIR
#   run-r -f run.R --args --mode=timing  --study=DIR
#   run-r -f run.R --args --mode=fit     --study=DIR --truth=t2 --workers=8
#   run-r -f run.R --args --mode=score   --study=DIR
PILOT_DIR <- "dev/simstudy/spatial-recoverability-pilot"
for (f in c("kernel.R", "seeds.R", "generator.R", "oracle.R", "predict.R", "score.R")) source(file.path(PILOT_DIR, f))

default_budget <- function() list(nburn = 3000L, niter = 4000L, thin = 10L, ext_nburn = 8000L, ext_niter = 16000L)

pilot_prepare <- function(study, communities = 10L, cells = pilot_cells(), jobs_dir) {
  dir.create(study, showWarnings = FALSE, recursive = TRUE)
  stopifnot(!file.exists(file.path(study, "config.rds")))
  dir.create(jobs_dir, showWarnings = FALSE, recursive = TRUE)
  if (!file.exists(file.path(study, "jobs"))) file.symlink(normalizePath(jobs_dir), file.path(study, "jobs"))
  for (d in c("inputs", "status", "scores", "logs")) dir.create(file.path(study, d), showWarnings = FALSE)
  config <- list(created = Sys.time(), communities = communities, cells = cells, thresholds = default_thresholds(),
                 budget = default_budget(), r_version = R.version.string,
                 packages = sapply(c("Rcpp", "RcppArmadillo", "posterior"), function(p) as.character(packageVersion(p))),
                 sources = tools::md5sum(c(list.files(PILOT_DIR, pattern = "\\.(R|cpp)$", full.names = TRUE),
                                           "dev/simstudy/spatial-amplitude-prior/diagnosis/ellipse.cpp")),
                 git = system("git rev-parse HEAD", intern = TRUE))
  saveRDS(config, file.path(study, "config.rds"))
  manifest <- list()
  for (i in seq_len(nrow(cells))) for (r in seq_len(communities)) {
    cm <- make_community(cells$truth[i], cells$n[i], r, range = cells$range[i])
    path <- file.path(study, "inputs", sprintf("%s-rep%d.rds", cells$cell[i], r))
    saveRDS(cm, path)
    manifest[[length(manifest) + 1L]] <- data.frame(cell = cells$cell[i], replicate = r, file = basename(path),
                                                    md5 = unname(tools::md5sum(path)), seed_sites = cm$seeds["sites"], seed_data = cm$seeds["data"])
    write.csv(cbind(cell = cells$cell[i], replicate = r, cm$conditioning), file.path(study, "logs", sprintf("generator-conditioning-%s-rep%d.csv", cells$cell[i], r)), row.names = FALSE)
  }
  write.csv(do.call(rbind, manifest), file.path(study, "input-manifest.csv"), row.names = FALSE)
  invisible(config)
}

job_table <- function(study) {
  man <- read.csv(file.path(study, "input-manifest.csv"), stringsAsFactors = FALSE)
  jobs <- merge(man, data.frame(arrangement = c("uniform", "clustered")), by = NULL)
  jobs$truth <- substr(jobs$cell, 1, 2)
  jobs$job <- sprintf("%s-rep%d-%s", jobs$cell, jobs$replicate, jobs$arrangement)
  jobs[order(jobs$job), ]
}

fit_with_extension <- function(cm, arrangement, budget, fit_fun, seed_label) {
  seed <- pilot_seed(seed_label)
  first <- fit_fun(nburn = budget$nburn, niter = budget$niter, thin = budget$thin, seed = seed)
  sc <- screen_traces(first$traces); phase <- "initial"
  if (!all(sc$pass)) {
    first <- fit_fun(nburn = budget$ext_nburn, niter = budget$ext_niter, thin = budget$thin * 4L, seed = pilot_seed(paste0(seed_label, "-extension")))
    sc <- screen_traces(first$traces); phase <- "extension"
  }
  list(fit = first, screen = sc, phase = phase, status = if (all(sc$pass)) "pass" else "unresolved")
}

run_job <- function(study, job, budget) {
  out <- file.path(study, "jobs", paste0(job$job, ".rds"))
  if (file.exists(out)) return(invisible("existing"))
  cm <- readRDS(file.path(study, "inputs", job$file)); arrangement <- job$arrangement
  t0 <- Sys.time(); S <- nrow(cm$species)
  if (cm$truth == "t1") {
    r <- fit_with_extension(cm, arrangement, budget, function(nburn, niter, thin, seed) run_joint_oracle(cm, arrangement, nburn, niter, thin, seed),
                            pilot_label("oracle", job$job))
    screen <- data.frame(species = seq_len(S), status = vapply(seq_len(S), function(s) if (all(r$screen$pass[r$screen$species == s])) "pass" else "unresolved", character(1)))
    oracle <- list(joint = r$fit); phases <- r$phase; factor_records <- r$fit$factor_record
    draws <- r$fit$site_draws; proposals <- r$fit$mean_proposals
  } else {
    fits <- lapply(seq_len(S), function(s) fit_with_extension(cm, arrangement, budget,
      function(nburn, niter, thin, seed) run_species_oracle(cm, arrangement, s, nburn, niter, thin, seed), pilot_label("oracle", job$job, cm$species$species[s])))
    screen <- data.frame(species = seq_len(S), status = vapply(fits, `[[`, character(1), "status"))
    oracle <- list(species = seq_len(S), oracles = lapply(fits, `[[`, "fit")); phases <- vapply(fits, `[[`, character(1), "phase")
    factor_records <- do.call(rbind, lapply(fits, function(f) f$fit$factor_record))
    draws <- lapply(fits, function(f) f$fit$site_draws[[1]]); proposals <- vapply(fits, function(f) mean(f$fit$mean_proposals), numeric(1))
  }
  prediction <- predict_lattice(cm, arrangement, oracle, seed = pilot_seed(pilot_label("lattice", job$job)))
  scores <- score_species(cm, arrangement, prediction, screen)
  screen_detail <- if (cm$truth == "t1") r$screen else do.call(rbind, lapply(seq_len(S), function(s) cbind(species = s, fits[[s]]$screen[, c("quantity", "rhat", "ess_bulk", "ess_tail", "pass")])))
  result <- list(job = job$job, cell = cm$cell, replicate = cm$replicate, arrangement = arrangement, truth = cm$truth,
                 site_draws = draws, prediction = prediction[c("probability_mean", "field_mean", "site_field_mean", "site_probability_mean")],
                 scores = scores, screen = screen, screen_detail = screen_detail, phases = phases, mean_proposals = proposals,
                 factor_records = factor_records, conditional_records = prediction$conditional_records,
                 budget = budget, elapsed_seconds = as.numeric(Sys.time() - t0, units = "secs"), finished = Sys.time())
  tmp <- paste0(out, ".tmp"); saveRDS(result, tmp); file.rename(tmp, out)
  writeLines(sprintf('{"job":"%s","status":"%s","elapsed_seconds":%.1f,"phases":"%s"}', job$job,
                     if (all(screen$status == "pass")) "pass" else "unresolved", result$elapsed_seconds, paste(phases, collapse = ",")),
             file.path(study, "status", paste0(job$job, ".json")))
  invisible("fitted")
}

pilot_fit <- function(study, truth, workers = 8L, budget = readRDS(file.path(study, "config.rds"))$budget, jobs = NULL) {
  compile_pilot_samplers(".")
  all_jobs <- job_table(study); all_jobs <- all_jobs[all_jobs$truth == truth, ]
  if (!is.null(jobs)) all_jobs <- all_jobs[all_jobs$job %in% jobs, ]
  res <- parallel::mclapply(seq_len(nrow(all_jobs)), function(i) tryCatch(run_job(study, all_jobs[i, ], budget),
    error = function(e) { writeLines(conditionMessage(e), file.path(study, "status", paste0(all_jobs$job[i], ".error"))); "error" }),
    mc.cores = workers, mc.preschedule = FALSE)
  table(unlist(res))
}

pilot_timing <- function(study, workers = 1L) {
  compile_pilot_samplers(".")
  jobs <- job_table(study); budget <- readRDS(file.path(study, "config.rds"))$budget
  pick <- c(jobs$job[jobs$truth == "t2" & grepl("n900", jobs$cell) & jobs$arrangement == "uniform"][1],
            jobs$job[jobs$truth == "t1" & grepl("r0.15-n900", jobs$cell) & jobs$arrangement == "uniform"][1])
  times <- lapply(pick, function(j) { t0 <- Sys.time(); run_job(study, jobs[jobs$job == j, ], budget); as.numeric(Sys.time() - t0, units = "secs") })
  names(times) <- pick
  n_jobs <- table(jobs$truth)
  projection <- list(timing = times, jobs = as.list(n_jobs),
                     projected_cpu_hours_upper = (times[[1]] * n_jobs[["t2"]] + times[[2]] * n_jobs[["t1"]]) / 3600,
                     note = "Upper bound: every job priced at the 900-site cost; extensions add up to 3x for the jobs that need them.")
  writeLines(jsonlite_free(projection), file.path(study, "timing.json"))
  projection
}

jsonlite_free <- function(x) paste(capture.output(dput(x)), collapse = "\n")

pilot_score <- function(study) {
  files <- list.files(file.path(study, "jobs"), pattern = "\\.rds$", full.names = TRUE)
  results <- lapply(files, readRDS)
  species <- do.call(rbind, lapply(results, `[[`, "scores"))
  write.csv(species, file.path(study, "scores", "species.csv"), row.names = FALSE)
  agg <- do.call(rbind, lapply(split(species, species$truth), aggregate_scores))
  write.csv(agg, file.path(study, "scores", "aggregate.csv"), row.names = FALSE)
  th <- readRDS(file.path(study, "config.rds"))$thresholds
  uni <- agg[agg$arrangement == "uniform", ]
  keys <- unique(uni[c("truth", "cell", "n", "range", "group")])
  gate <- do.call(rbind, lapply(seq_len(nrow(keys)), function(i) {
    rows <- merge(uni, keys[i, ]); cbind(keys[i, ], gate_cells(rows, th)) }))
  gate$class <- ifelse(gate$truth == "t2", gate$group, NA_character_)
  write.csv(gate, file.path(study, "scores", "gate.csv"), row.names = FALSE)
  writeLines(jsonlite_free(read_map(gate)), file.path(study, "scores", "map.json"))
  cond <- do.call(rbind, lapply(results, function(r) cbind(job = r$job, rbind(r$factor_records, within(r$conditional_records, condition <- NA)[names(r$factor_records)]))))
  write.csv(cond, file.path(study, "scores", "conditioning.csv"), row.names = FALSE)
  status <- do.call(rbind, lapply(results, function(r) data.frame(job = r$job, truth = r$truth, status = if (all(r$screen$status == "pass")) "pass" else "unresolved",
                                                                  phases = paste(r$phases, collapse = ","), elapsed_seconds = r$elapsed_seconds, mean_proposals = mean(r$mean_proposals))))
  write.csv(status, file.path(study, "scores", "jobs.csv"), row.names = FALSE)
  invisible(gate)
}

if (sys.nframe() == 0L && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  get <- function(key, default = NULL) { v <- sub(paste0("^--", key, "="), "", grep(paste0("^--", key, "="), args, value = TRUE)); if (length(v)) v else default }
  mode <- get("mode"); study <- get("study")
  switch(mode,
    prepare = pilot_prepare(study, jobs_dir = get("jobs")),
    timing = print(pilot_timing(study)),
    fit = print(pilot_fit(study, truth = get("truth"), workers = as.integer(get("workers", "8")))),
    score = print(pilot_score(study)),
    stop("Unknown mode ", mode))
  cat("Completed mode", mode, "at", format(Sys.time()), "\n")
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-run.R`
Expected: `Run tests passed.` Note that `sys.nframe() == 0L` is false when the file is sourced from the test, so the command-line block does not run there.

- [ ] **Step 5: Commit**

```bash
git add dev/simstudy/spatial-recoverability-pilot/run.R dev/simstudy/spatial-recoverability-pilot/test-run.R
git commit -m "Add the recoverability pilot driver with prepare, timing, fit and score modes"
```

---

### Task 7: Prepare the study, run the timing probe and launch the grid

**Files:**
- Create (ignored): `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010/`
- Create (volume): `/Volumes/occJSDM Research/occJSDM/spatial-recoverability-20261010/jobs/`

**Interfaces:**
- Consumes: `run.R` modes.
- Produces: a frozen config, 120 inputs with manifest, `timing.json`, then 240 job files.

- [ ] **Step 1: Run every test file once more under the launcher**

Run, from the worktree root:
```bash
for t in seeds kernel generator oracle score run; do dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/test-$t.R || exit 1; done
```
Expected: six "passed" lines.

- [ ] **Step 2: Prepare**

```bash
STUDY=/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010
JOBS="/Volumes/occJSDM Research/occJSDM/spatial-recoverability-20261010/jobs"
dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/run.R --args --mode=prepare --study=$STUDY --jobs="$JOBS"
ls $STUDY/inputs | wc -l
```
Expected: `120`, `config.rds` older than every input, `input-manifest.csv` with 120 rows.

- [ ] **Step 3: Timing probe**

```bash
dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/run.R --args --mode=timing --study=$STUDY
cat $STUDY/timing.json
```
Expected: two elapsed times and a projected upper bound in CPU hours. Decision rule, fixed here: proceed if the upper bound divided by 8 workers is under 12 wall-clock hours. Otherwise stop and report to Doug with the numbers; do not reduce the grid unilaterally.

- [ ] **Step 4: Launch truth 2 then truth 1, detached**

```bash
nohup dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/run.R --args --mode=fit --study=$STUDY --truth=t2 --workers=8 > $STUDY/logs/fit-t2.log 2>&1 &
```
When `fit-t2.log` ends with `Completed mode fit`, launch truth 1 the same way with `--truth=t1` into `fit-t1.log`. Monitor with `ls $STUDY/status | wc -l` (target 60 then 240) and `ls $STUDY/status/*.error 2>/dev/null`.

- [ ] **Step 5: Score**

```bash
dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/run.R --args --mode=score --study=$STUDY
cat $STUDY/scores/map.json
```
Expected: `species.csv` with 240 × 12 = 2,880 rows, `gate.csv` with 9 truth 1 rows and 9 truth 2 rows, `map.json` naming a branch.

- [ ] **Step 6: Record the run in the spec's revision record and commit**

Add one line to the spec's "Revision record" naming the study directory, the timing numbers and the launch times. Commit:
```bash
git add dev/superpowers/specs/2026-10-10-spatial-recoverability-pilot-design.md
git commit -m "Record the recoverability pilot timing probe and launch"
```

---

### Task 8: Independent verification, compact results, report and README

**Files:**
- Create: `dev/simstudy/spatial-recoverability-pilot/verify.R`
- Create: `dev/simstudy/spatial-recoverability-pilot/plot.R`
- Create: `dev/simstudy/spatial-recoverability-pilot/results/` (compact CSVs copied from the study, plus `recoverability-map.png`)
- Create: `dev/simstudy/spatial-recoverability-pilot/REPORT.md`
- Create: `dev/simstudy/spatial-recoverability-pilot/README.md`
- Modify: `TODO.md` (one bullet in the spatial section pointing at the report)

**Interfaces:**
- Consumes: `STUDY/scores/*.csv`, `STUDY/config.rds`, `STUDY/input-manifest.csv`, `STUDY/logs/generator-conditioning-*.csv`.
- Produces: `verify.R` exits non-zero on any failed check and writes `results/verification.csv`.

- [ ] **Step 1: Write verify.R**

```r
# Independent verification of the recoverability pilot. Recomputes every gate decision from
# species.csv with its own code, checks seeds, inputs and conditioning logs. Run after score.
# run-r -f verify.R --args --study=DIR
args <- commandArgs(trailingOnly = TRUE)
study <- sub("^--study=", "", grep("^--study=", args, value = TRUE))
source("dev/simstudy/spatial-recoverability-pilot/seeds.R")
checks <- list(); add <- function(name, ok, detail = "") checks[[length(checks) + 1L]] <<- data.frame(check = name, pass = isTRUE(ok), detail = detail)

cfg <- readRDS(file.path(study, "config.rds")); th <- cfg$thresholds
species <- read.csv(file.path(study, "scores", "species.csv"), stringsAsFactors = FALSE)
gate <- read.csv(file.path(study, "scores", "gate.csv"), stringsAsFactors = FALSE)
man <- read.csv(file.path(study, "input-manifest.csv"), stringsAsFactors = FALSE)

add("species rows", nrow(species) == nrow(man) * 2L * 12L, nrow(species))
add("input hashes unchanged", all(unname(tools::md5sum(file.path(study, "inputs", man$file))) == man$md5))
add("seeds distinct", !anyDuplicated(c(man$seed_sites, man$seed_data)) > 0)
add("config predates inputs", file.mtime(file.path(study, "config.rds")) <= min(file.mtime(file.path(study, "inputs", man$file))))

# Independent gate recomputation, written without reference to score.R.
u <- species[species$arrangement == "uniform" & species$occupancy %in% c(0.25, 0.5, 0.75), ]
u$group <- ifelse(u$truth == "t2", u$class, "all")
per_level <- aggregate(u[c("contrast_correlation", "field_slope", "headroom_pp")], u[c("cell", "group", "replicate", "occupancy")], mean)
per_comm <- aggregate(per_level[c("contrast_correlation", "field_slope", "headroom_pp")], per_level[c("cell", "group", "replicate")], mean)
unres <- aggregate(list(unresolved = u$status == "unresolved"), u[c("cell", "group", "replicate")], any)
per_comm <- merge(per_comm, unres)
decide <- function(d) {
  if (any(d$unresolved)) return("unresolved")
  p <- d$contrast_correlation; s <- d$field_slope; h <- d$headroom_pp
  ci <- mean(h) + c(-1, 1) * qt(0.975, length(h) - 1) * sd(h) / sqrt(length(h))
  ok <- median(p) >= th$pattern_median && sum(p >= th$pattern_floor) >= th$pattern_count &&
    median(s) >= th$slope_low && median(s) <= th$slope_high && mean(h) >= th$headroom_pp && ci[1] > 0
  if (ok) "pass" else "fail"
}
recomputed <- do.call(rbind, lapply(split(per_comm, per_comm[c("cell", "group")], drop = TRUE), function(d) data.frame(cell = d$cell[1], group = d$group[1], decision = decide(d))))
cmp <- merge(gate[c("cell", "group", "decision")], recomputed, by = c("cell", "group"), suffixes = c("_scored", "_recomputed"))
add("gate decisions recomputed identically", nrow(cmp) == nrow(gate) && all(cmp$decision_scored == cmp$decision_recomputed), paste(sum(cmp$decision_scored != cmp$decision_recomputed), "mismatches"))

cond <- read.csv(file.path(study, "scores", "conditioning.csv"), stringsAsFactors = FALSE)
add("conditional reconstruction within 1e-5", all(cond$error[grepl("^conditional", cond$what)] <= 1e-5))
gen <- do.call(rbind, lapply(list.files(file.path(study, "logs"), pattern = "^generator-conditioning", full.names = TRUE), read.csv))
add("generator factor errors within 1e-5", all(gen$error <= 1e-5))
add("every job has a status and no errors", length(list.files(file.path(study, "status"), pattern = "\\.json$")) == nrow(man) * 2L &&
      length(list.files(file.path(study, "status"), pattern = "\\.error$")) == 0L)

out <- do.call(rbind, checks)
dir.create("dev/simstudy/spatial-recoverability-pilot/results", showWarnings = FALSE)
write.csv(out, "dev/simstudy/spatial-recoverability-pilot/results/verification.csv", row.names = FALSE)
print(out)
if (!all(out$pass)) stop("Verification failed")
cat("Verification passed.\n")
```

- [ ] **Step 2: Write plot.R**

```r
# Copies compact results into the study folder and draws the recoverability map.
# run-r -f plot.R --args --study=DIR
args <- commandArgs(trailingOnly = TRUE)
study <- sub("^--study=", "", grep("^--study=", args, value = TRUE))
out <- "dev/simstudy/spatial-recoverability-pilot/results"; dir.create(out, showWarnings = FALSE)
for (f in c("species.csv", "aggregate.csv", "gate.csv", "map.json", "conditioning.csv", "jobs.csv")) file.copy(file.path(study, "scores", f), file.path(out, f), overwrite = TRUE)
file.copy(file.path(study, "input-manifest.csv"), file.path(out, "input-manifest.csv"), overwrite = TRUE)
file.copy(file.path(study, "timing.json"), file.path(out, "timing.json"), overwrite = TRUE)
agg <- read.csv(file.path(study, "scores", "aggregate.csv")); gate <- read.csv(file.path(study, "scores", "gate.csv"))
png(file.path(out, "recoverability-map.png"), width = 1800, height = 900, res = 150)
par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
cols <- c(pass = "#2a9d8f", fail = "#e76f51", unresolved = "#999999")
for (tr in c("t1", "t2")) {
  g <- gate[gate$truth == tr, ]; a <- agg[agg$truth == tr, ]
  rows <- if (tr == "t1") sort(unique(g$range)) else c("short", "medium", "long"); key <- if (tr == "t1") g$range else g$group
  plot(NA, xlim = c(0.5, 3.5), ylim = c(0.5, length(rows) + 0.5), xaxt = "n", yaxt = "n",
       xlab = "training sites", ylab = if (tr == "t1") "range" else "dispersal class",
       main = if (tr == "t1") "Truth 1: common gradients" else "Truth 2: dispersal fields")
  axis(1, 1:3, c(100, 300, 900)); axis(2, seq_along(rows), rows, las = 1)
  for (i in seq_len(nrow(g))) {
    x <- match(g$n[i], c(100, 300, 900)); y <- match(key[i], rows)
    rect(x - 0.45, y - 0.45, x + 0.45, y + 0.45, col = cols[g$decision[i]], border = NA)
    text(x, y + 0.2, sprintf("r = %.2f, slope %.2f", g$pattern_median[i], g$slope_median[i]), cex = 0.7)
    text(x, y - 0.2, sprintf("%.1f pp [%.1f, %.1f]", g$headroom_mean[i], g$headroom_low[i], g$headroom_high[i]), cex = 0.7)
    ca <- a[a$arrangement == "clustered" & a$n == g$n[i] & (if (tr == "t1") a$range == g$range[i] else a$group == g$group[i]), ]
    text(x, y - 0.38, sprintf("clustered r = %.2f", median(ca$contrast_correlation)), cex = 0.6, col = "grey20")
  }
}
dev.off()
cat("Wrote", file.path(out, "recoverability-map.png"), "\n")
```

- [ ] **Step 3: Run verify and plot**

```bash
STUDY=/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010
dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/verify.R --args --study=$STUDY
dev/simstudy/spatial-recoverability-pilot/run-r -f dev/simstudy/spatial-recoverability-pilot/plot.R --args --study=$STUDY
```
Expected: `Verification passed.` and the PNG written. Inspect the PNG with the Read tool for clipped labels.

- [ ] **Step 4: Write REPORT.md**

Structure, each paragraph on one line, numbers read from `results/gate.csv`, `results/aggregate.csv` and `results/jobs.csv`, never typed from memory:

- Title and one-paragraph summary: the branch from `map.json`, the recommended setting, and the sentence that truth 1's advantage follows from pooling twelve species onto two fields.
- "Design as run": cells, communities, arrangements, scaled separations, MCMC budget, pinned runtime, study directory, timing probe numbers.
- "Recoverability map": bullet list, one bullet per gate row, giving decision, pattern median and count, slope median, headroom mean and interval, and the clustered-arrangement pattern median alongside.
- "Reading the map": apply the spec's branch rule and state the stage 2 scope it implies.
- "Diagnostics and conditioning": how many jobs needed extension, how many are unresolved, factorisation routes used by cell, the largest reconstruction error.
- "Mechanism attribution": the spec's paragraph, verbatim.
- "Limits": oracle-only, no package fits, no detection, three ranges, one landscape shape, 10 communities.
- "Files": what is in `results/` and where the raw jobs live.

- [ ] **Step 5: Write README.md**

Short: purpose, one-paragraph results pointer to REPORT.md, reproduction commands (the exact `run-r` commands from Task 7 in a fenced block), file list.

- [ ] **Step 6: Add the TODO.md bullet and commit**

In `TODO.md`, under the spatial section, add one line: `- **Spatial recoverability pilot complete, <date>.** See dev/simstudy/spatial-recoverability-pilot/REPORT.md. Branch <branch>; stage 2 specification <is / is not> proposed.` (single line, no table).

```bash
git add dev/simstudy/spatial-recoverability-pilot TODO.md
git commit -m "Add the recoverability pilot verification, compact results and report"
```

- [ ] **Step 7: Request a whole-branch review**

Use superpowers:requesting-code-review on the branch diff against `main`, then address findings before handing Doug the report.

---

## Deviations recorded during execution

- `stable_factor` computes eigenvalues only (not vectors) for the condition number when the Cholesky route is taken, and `prepare` runs communities in parallel with `--workers`, because generating a 3,400-point community takes minutes of dense eigendecomposition.
- `predict.R` gained `layout_geometry()`, which the spec's section 4 requires and the first draft of Task 5 omitted; each job stores it and `score` writes `geometry.csv`.
- `score` writes `factors.csv` and `conditionals.csv` separately instead of one combined conditioning file, plus `screen.csv` with the per-quantity diagnostics.
- Status files and `map.json` are written with `jsonlite`, which is in the pinned library.
- Tests were written together with their implementations and run once per task rather than first run against the missing file.
- Two spec section 6 verification items have no committed test: the empirical variogram check of each class's range, and the agreement of the Cholesky route with the eigendecomposition route on prior draws with a fixed seed. The whole-branch review on 11 October ran a one-off variogram check on the stored inputs and reported empirical lag correlations matching the nominal kernels (scratch script, not committed); the route-agreement check has no test and was not run.
- After the plan was written, the truth 1 sampler was swapped for `block-ellipse.cpp`, the resolution extension (`resolve.R`, modes `prescribe-resolution` and `resolve`) was added, and `earlier-readings.R` was added to commit the earlier readings and the agreement summary. All three are recorded in the spec's Revision record.

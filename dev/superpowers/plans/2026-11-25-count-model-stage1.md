# Count Model Stage 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and test, in isolation, the likelihood pieces a count model needs: a Polya-Gamma sampler for any positive shape, a negative-binomial log density with effort offsets, a per-species update of the size, and one name for the count model. No fit changes and `runOccJSDM()` still refuses count data.

**Architecture:** Deterministic formulas (Polya-Gamma moments, negative-binomial log density) go in a new header, `src/count.h`. The sampler goes in `src/jsdm.cpp`, next to the unit sampler it reuses and on the same random stream. The size update is plain R in a new file, `R/count_pieces.R`, written like `R/b0_shift.R`: a Metropolis step on R's random stream with a step size adapted during burn-in. Nothing in this stage is called by the sampler loop; stage 2 does the wiring.

**Tech Stack:** R 4.6.1 (`/usr/local/bin/Rscript`), Rcpp and RcppArmadillo, testthat 3e, devtools for `load_all()` and `test()`.

**Spec:** `dev/superpowers/specs/2026-11-12-count-model-design.md`, sections *How it is fitted*, *Staged plan* (stage 1) and questions 5 and 9. Until PR 42 is merged, the spec and this plan are on the branch `claude/occjsdm-count-model-plan-992fd8`.

## Global Constraints

- No change to any existing fit. Every existing test passes unchanged, and no existing function's draws change. `runOccJSDM()` still stops with "Counts model not supported yet".
- The count model's name is "counts" everywhere in R and C++. The string "count" as a model name must not remain.
- Every random draw in C++ uses `runif()` and `rnorm()` from `src/rng.h`, on the main thread. No R random numbers in C++, no draws inside a TBB worker.
- One size per species (Alex, question 5). The prior on the size is proposed here and not yet confirmed by Alex: a half-normal with standard deviation 1 on the inverse square root of the size, independent between species.
- Above a shape of 170 the Polya-Gamma draw is a normal draw with the Polya-Gamma mean and variance (Alex, question 9). No saddlepoint sampler in this stage; Task 5 measures whether one is needed.
- No Poisson family.
- Never run `R CMD build`, `R CMD check`, `devtools::check()` or `devtools::build()` on a checkout. Test with `devtools::load_all()` and `devtools::test()`.
- After changing any `// [[Rcpp::export]]` signature, run `Rscript -e 'Rcpp::compileAttributes()'`. Never edit `R/RcppExports.R` or `src/RcppExports.cpp` by hand.
- No em-dashes anywhere, including commit messages and code comments. Markdown follows the project `CLAUDE.md`: one line per paragraph, no pipe tables.
- Commit messages end with `Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>`. Do not push or open a pull request until Doug says so.
- Work on a new branch, `claude/count-model-stage1`, cut from `origin/main` in its own worktree. Never switch branches in the main checkout.

## Decisions made in this plan

These go beyond the spec. Each was checked numerically on 25 November 2026 before the plan was written, and each is for Doug and Alex to confirm when they review the plan.

- **The fractional part of the shape** is drawn as the first 20 terms of the Polya-Gamma series plus the exact mean of the terms left out. The variance left out is under 0.04% of the total wherever this route is used.
- **An extreme predictor, beyond 12 in absolute value, uses an inverse Gaussian draw** at every shape, with mean `b / (2 * abs(z))` and shape parameter `b^2 / 4`. This is the large-predictor limit of the Polya-Gamma distribution. It matched a 1,000-term series reference in mean, variance and distribution at shapes from 0.3 to 200. A normal draw there would be negative 8% of the time at a shape of 0.3.
- **The size is kept between `1e-6` and `1e6`.** A proposal outside is rejected.
- **The size of a rare species is barely identified.** In the prototype, a species with 74 non-zero counts in 400 and a true size of 5 had a posterior median of 24 and an upper limit above 30,000. Task 5 reports this by rarity. **Corrected 26 November 2026:** this bullet first said that the prior then pulls a rare species towards the Poisson, as intended. That is wrong, and it generalised from the one species above. The prior's median size is about 2.2 and about a third of its mass is on sizes below 1. Stage 1 measured it: with fewer than 10 non-zero counts the posterior median size was 1 to 2 whatever the truth, close to the prior and not towards the Poisson. See `dev/count-model/stage1-report.md` on the stage 1 branch (pull request 47) and question 5 of the design.

## Out of scope: a defect found in the existing sampler

`samplepg_fast()`, in both `src/jsdm.cpp` and `src/functions.cpp`, switches to a normal draw when the predictor is beyond 12 in absolute value, with variance `1 / (4 z^3)`. The true Polya-Gamma variance there is `1 / (2 z^3)`, twice as large. Checked against the defining series on 25 November 2026: the ratio is 0.500 at 13 and at 20. It affects every binary, occupancy and two-stage fit, but only in cells whose occupancy predictor is beyond 12, where the weight is small. Fixing it changes existing fits draw for draw, so it is not part of this stage. This stage must not copy it: the count sampler never calls `samplepg_fast()` beyond 12. Leave both copies untouched.

## Review Focus

- **A shape below 1** (a zero count with a size below 1): the draw has no whole-number part and must still be positive with the right mean. Pinned in Task 3.
- **An extreme predictor** (a tiny effort or a very rare species gives a predictor of -40 or beyond): moments and draws must stay finite and positive, with no overflow at 700 or 5,000. Pinned in Tasks 2 and 3.
- **A predictor of exactly zero or very near it:** the moment formulas divide by the predictor and must switch to their series. Pinned in Task 2.
- **Missing counts:** an `NA` count gives an `NA` log density and is ignored by the size update. Pinned in Tasks 2 and 4.
- **A species never counted:** all its counts are zero, so the data say almost nothing about its size. The update must return a finite size and never fail. Pinned in Task 4.

---

### Task 1: One name for the count model

**Files:**
- Modify: `R/jsdmfun.R` (the three tests of `model == "count"`, near lines 775, 1315 and 1729)
- Create: `tests/testthat/test-count-pieces.R`

**Interfaces:**
- Produces: `SD(eta, "counts")` works and "count" is no longer a model name.

- [ ] **Step 1: Record the baseline.** Run the whole suite on the untouched branch in the background and keep the last lines for the pull request description.

Run: `Rscript -e 'devtools::test()' 2>&1 | tail -n 8`
Expected: the pass, fail, warn and skip counts of `main`. Write them down.

- [ ] **Step 2: Write the failing test** in `tests/testthat/test-count-pieces.R`

```r
# Stage 1 of the count model: likelihood pieces tested alone. Nothing here is
# called by the sampler loop yet. Design:
# dev/superpowers/specs/2026-11-12-count-model-design.md.

test_that("the count model has one name, counts", {
  eta <- matrix(c(0, log(2), log(3), 0, 0, 0), 3, 2)
  expect_equal(occJSDM:::SD(eta, "counts"), c(stats::sd(c(1, 2, 3)), 0))
  expect_error(occJSDM:::SD(eta, "count"), "Using unavailable model")
})
```

- [ ] **Step 3: Run it and see it fail**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-count-pieces.R")'`
Expected: FAIL, "Using unavailable model" raised for "counts".

- [ ] **Step 4: Rename.** In `R/jsdmfun.R` change the three `"count"` model tests to `"counts"`. Change nothing else.

- [ ] **Step 5: Verify**

Run: `grep -n '"count"' R/*.R src/*.cpp`
Expected: no output.

Run the test file again. Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add R/jsdmfun.R tests/testthat/test-count-pieces.R
git commit -m "Use one name, counts, for the count model"
```

### Task 2: Polya-Gamma moments and the negative-binomial log density

**Files:**
- Create: `src/count.h`
- Modify: `src/jsdm.cpp` (include the header; two exported wrappers in a new section headed "Count model" at the end of the file)
- Regenerate: `R/RcppExports.R`, `src/RcppExports.cpp`
- Test: `tests/testthat/test-count-pieces.R`

**Interfaces:**
- Produces, in `src/count.h`, all `inline`:
  - `double pg_mean(double b, double z)`: mean of a Polya-Gamma variable with shape `b` and predictor `z`.
  - `double pg_var(double b, double z)`: its variance.
  - `double nb_logpmf(double y, double log_mean, double r)`: log probability of count `y` under a negative binomial with mean `exp(log_mean)` and size `r`.
- Produces, exported to R (internal, no roxygen):
  - `pgMoments_cpp(b, z)`: numeric vectors of equal length; returns a two-column matrix, mean then variance.
  - `nbLogPmf_cpp(y, log_mean, r)`: `y` and `log_mean` are matrices with one row per sample and one column per species, `r` has one size per species; returns a matrix of the same shape.

- [ ] **Step 1: Write the failing tests**

```r
test_that("Polya-Gamma moments match the defining series", {
  d <- function(k, z) 1 / (2 * pi^2 * ((k - 0.5)^2 + z^2 / (4 * pi^2)))
  for (z in c(0, 1e-4, 0.005, 0.5, 2, 12, 30)) {
    m <- pgMoments_cpp(c(1, 2.5), c(z, -z))
    expect_equal(m[, 1], c(1, 2.5) * sum(d(1:2e6, z)), tolerance = 1e-5)
    expect_equal(m[, 2], c(1, 2.5) * sum(d(1:2e6, z)^2), tolerance = 1e-5)
  }
  big <- pgMoments_cpp(c(1, 1), c(700, -5000))
  expect_true(all(is.finite(big)) && all(big > 0))
  expect_equal(big[1, ], c(1 / 1400, 1 / (2 * 700^3)), tolerance = 1e-8)
})

test_that("the negative-binomial log density matches dnbinom", {
  y <- matrix(c(0, 1, 7, 250, 1e5, 0, 3, NA), 4, 2)
  log_mean <- matrix(c(-30, -2, 1.5, 40, 11, -700, 0, 1), 4, 2)
  r <- c(0.3, 50)
  got <- nbLogPmf_cpp(y, log_mean, r)
  want <- stats::dnbinom(y, size = rep(r, each = 4), mu = exp(log_mean), log = TRUE)
  expect_equal(got, matrix(want, 4, 2), tolerance = 1e-10)
  expect_true(is.na(got[4, 2]))
  expect_identical(nbLogPmf_cpp(matrix(c(-1, 2.5), 2, 1), matrix(0, 2, 1), 1),
                   matrix(-Inf, 2, 1))
  expect_error(nbLogPmf_cpp(y, log_mean, c(0.3, 0)), "size must be positive")
  expect_error(nbLogPmf_cpp(y, log_mean, 0.3), "one size per species")
  expect_error(nbLogPmf_cpp(y, log_mean[1:2, ], r), "same dimensions")
})
```

- [ ] **Step 2: Run and see both fail** with "could not find function".

- [ ] **Step 3: Write `src/count.h`.** The formulas are fixed; the thresholds are decisions.

```
t = tanh(|z| / 2)

pg_mean(b, z) = b * (1/4 - z^2/48)                        when |z| < 1e-3
              = b * t / (2 |z|)                           otherwise

pg_var(b, z)  = b * (1/24 - z^2/120)                      when |z| < 1e-2
              = b * (2 t - |z| (1 - t^2)) / (4 |z|^3)     otherwise

nb_logpmf(y, log_mean, r) =
    lgamma(y + r) - lgamma(r) - lgamma(y + 1) + r log(r) + y * log_mean
    - (y + r) * logsumexp(log_mean, log(r))
```

`nb_logpmf` returns `NA_REAL` when `y` is `NA`, and `R_NegInf` when `y` is negative or not a whole number. Write `logsumexp(a, b)` as `max(a, b) + log1p(exp(-abs(a - b)))`. Use the `tanh` form of the variance, not `sinh` and `cosh`, which overflow.

- [ ] **Step 4: Write the two wrappers in `src/jsdm.cpp`** and run `Rscript -e 'Rcpp::compileAttributes()'`. `pgMoments_cpp` stops with "b and z must have the same length". `nbLogPmf_cpp` stops with "y and log_mean must have the same dimensions", "one size per species" when `r` is not one per column, and "size must be positive" when any size is not finite and positive.

- [ ] **Step 5: Run the test file.** Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/count.h src/jsdm.cpp src/RcppExports.cpp R/RcppExports.R tests/testthat/test-count-pieces.R
git commit -m "Add Polya-Gamma moments and a negative-binomial log density"
```

### Task 3: The Polya-Gamma sampler for any positive shape

**Files:**
- Modify: `src/jsdm.cpp` (in the "Count model" section)
- Regenerate: `R/RcppExports.R`, `src/RcppExports.cpp`
- Test: `tests/testthat/test-count-pieces.R`

**Interfaces:**
- Consumes: `pg_mean()`, `pg_var()` from `src/count.h`; `samplepg_fast()`, `randinvg()` already in `src/jsdm.cpp`; `runif()`, `rnorm()` from `src/rng.h`; `pgMoments_cpp()` in the tests.
- Produces: `samplePGcount_cpp(shape, phi)`, exported to R. Two matrices of the same dimensions; returns a matrix of Polya-Gamma draws with shape `shape[i, j]` at predictor `phi[i, j]`, drawn in column-major order. Stage 2 calls it with `shape = y + r` and `phi` the log mean minus the log size.

- [ ] **Step 1: Write the failing tests**

```r
test_that("at whole-number shapes the count sampler is the existing sampler, summed", {
  z <- c(-11, -2, 0, 0.3, 5, 12)
  setOccJSDMSeed(4242L)
  unit <- samplePGvariables_parallel(matrix(z, 3, 6, byrow = TRUE))
  setOccJSDMSeed(4242L)
  got <- samplePGcount_cpp(matrix(3, 1, 6), matrix(z, 1, 6))
  expect_equal(as.vector(got), colSums(unit), tolerance = 1e-12)
})

test_that("count sampler draws have the Polya-Gamma mean and variance at every kind of shape", {
  grid <- expand.grid(b = c(0.1, 0.5, 0.9, 1, 2.3, 13.5, 169.5, 171, 5000.25),
                      z = c(0, 1, -5, 12, -20, 40))
  setOccJSDMSeed(99L)
  for (i in seq_len(nrow(grid))) {
    b <- grid$b[i]; z <- grid$z[i]
    n <- if (b > 3) 2e4 else 2e5
    x <- as.vector(samplePGcount_cpp(matrix(b, n, 1), matrix(z, n, 1)))
    m <- pgMoments_cpp(b, z)
    expect_true(all(is.finite(x)) && all(x > 0), info = paste(b, z))
    expect_lt(abs(mean(x) - m[1]) / sqrt(m[2] / n), 4.5)
    expect_equal(stats::var(x), m[2], tolerance = 0.08, info = paste(b, z))
  }
})

test_that("the fractional and extreme-predictor routes have the right distribution, not only the right moments", {
  series <- function(n, b, z, K = 1000) {
    d <- 1 / (2 * pi^2 * ((seq_len(K) - 0.5)^2 + z^2 / (4 * pi^2)))
    as.vector(matrix(stats::rgamma(n * K, b), n, K) %*% d)
  }
  set.seed(20261125)
  setOccJSDMSeed(20261125L)
  for (case in list(c(0.5, 1), c(0.1, -8), c(0.3, 13), c(1, -20), c(5, 12.5))) {
    b <- case[1]; z <- case[2]
    x <- as.vector(samplePGcount_cpp(matrix(b, 5000, 1), matrix(z, 5000, 1)))
    p <- suppressWarnings(stats::ks.test(x, series(5000, b, z))$p.value)
    expect_gt(p, 0.001)
  }
})

test_that("count sampler draws are reproducible and validate their inputs", {
  shape <- matrix(c(0.3, 2, 7.5, 200), 2, 2)
  phi <- matrix(c(-1, 0, 3, -14), 2, 2)
  setOccJSDMSeed(7L); a <- samplePGcount_cpp(shape, phi)
  setOccJSDMSeed(7L); b <- samplePGcount_cpp(shape, phi)
  expect_identical(a, b)
  expect_error(samplePGcount_cpp(shape, phi[, 1, drop = FALSE]), "same dimensions")
  expect_error(samplePGcount_cpp(replace(shape, 1, 0), phi), "shape must be positive")
  expect_error(samplePGcount_cpp(shape, replace(phi, 1, NA)), "phi must be finite")
})
```

- [ ] **Step 2: Run and see them fail** with "could not find function".

- [ ] **Step 3: Write `static double rpg_shape(double b, double z)` and `samplePGcount_cpp()`.** The route for each draw is fixed by the plan:

```
constants: PG_EXTREME_Z = 12.0, PG_NORMAL_SHAPE = 170.0, PG_FRAC_TERMS = 20

if |z| > PG_EXTREME_Z:                       inverse Gaussian, any shape
    lambda = b^2 / 4
    return lambda * randinvg( (b / (2 |z|)) / lambda )

if b > PG_NORMAL_SHAPE:                      normal, exact moments
    x = pg_mean(b, z) + sqrt(pg_var(b, z)) * rnorm()
    return x > 1e-12 ? x : pg_mean(b, z)

n = floor(b), f = b - n                      exact route
x = sum of n draws of samplepg_fast(z)
if f > 1e-8 or n == 0:
    d_k = 1 / (2 pi^2 ((k - 1/2)^2 + z^2 / (4 pi^2))),   k = 1 .. PG_FRAC_TERMS
    x += sum_k G_k d_k + ( pg_mean(f, z) - f * sum_k d_k ),   G_k ~ Gamma(f, 1)
return x
```

Draw `Gamma(f, 1)` for `f` below 1 with Marsaglia and Tsang's method at shape `f + 1`, using `rnorm()` and `runif()` from `src/rng.h`, then multiply by `runif()` raised to `1 / f`. Do not use `std::gamma_distribution` or `R::rgamma`: the first keeps hidden state between calls and the second is R's stream. The `n == 0` clause keeps a draw at a tiny shape positive. Loop over cells on the calling thread in column-major order. `samplePGcount_cpp` stops with "shape and phi must have the same dimensions", "shape must be positive" and "phi must be finite". Do not edit `samplepg_fast()`.

- [ ] **Step 4: Run `Rcpp::compileAttributes()` and the test file.** Expected: PASS. The moments test takes under a minute.

- [ ] **Step 5: Commit**

```bash
git add src/jsdm.cpp src/RcppExports.cpp R/RcppExports.R tests/testthat/test-count-pieces.R
git commit -m "Add a Polya-Gamma sampler for any positive shape"
```

### Task 4: Effort offsets and the size update

**Files:**
- Create: `R/count_pieces.R`
- Test: `tests/testthat/test-count-pieces.R`

**Interfaces:**
- Consumes: `nbLogPmf_cpp(y, log_mean, r)` from Task 2.
- Produces, all internal:
  - `count_log_mean(eta, idx_site, log_effort)`: `eta` is sites by species, `idx_site` maps each sample to its site, `log_effort` has one value per sample. Returns the samples-by-species matrix `eta[idx_site, ] + log_effort`.
  - `count_size_logprior(log_r, prior_sd)`: log prior density of the log size, up to a constant, when the inverse square root of the size is half-normal: `-1 / (2 * prior_sd^2 * exp(log_r)) - log_r / 2`.
  - `count_size_step(r, y, log_mean, include, prior_sd, step, adapt, t)`: one Metropolis proposal per species. `include` is a logical matrix shaped like `y`, `TRUE` where the cell counts towards the likelihood (every cell in the plain count JSDM, occupied cells in the shared-predictor model). Returns `list(r, accepted, step, t)`, with `accepted` logical and `step` numeric, one per species.
  - Constants: `count_size_target_acceptance <- 0.44`, `count_size_initial_step <- 0.5`, `count_size_bounds <- c(1e-6, 1e6)`.

- [ ] **Step 1: Write the failing tests**

```r
test_that("count_log_mean adds each sample's log effort to its site's predictor", {
  eta <- matrix(1:6, 3, 2)
  expect_equal(occJSDM:::count_log_mean(eta, c(1L, 1L, 3L), log(c(1, 2, 4))),
               eta[c(1, 1, 3), ] + log(c(1, 2, 4)))
  expect_error(occJSDM:::count_log_mean(eta, c(1L, 4L), c(0, 0)), "idx_site")
  expect_error(occJSDM:::count_log_mean(eta, c(1L, 2L), c(0, NA)), "log_effort")
})

test_that("the size prior is a half-normal on the inverse square root of the size", {
  log_r <- seq(-3, 5, length.out = 9)
  r <- exp(log_r)
  direct <- log(2 * stats::dnorm(r^(-1/2), 0, 0.7) * 0.5 * r^(-1/2))
  gap <- direct - occJSDM:::count_size_logprior(log_r, 0.7)
  expect_equal(gap, rep(gap[1], 9), tolerance = 1e-10)
})

size_run <- function(y, log_mean, include = matrix(TRUE, nrow(y), ncol(y)),
                     iterations = 3000, burn = 1000) {
  S <- ncol(y)
  r <- rep(1, S); step <- rep(occJSDM:::count_size_initial_step, S); t <- 0L
  keep <- matrix(NA_real_, iterations - burn, S)
  accepted <- keep
  for (it in seq_len(iterations)) {
    o <- occJSDM:::count_size_step(r, y, log_mean, include, 1, step, it <= burn, t)
    r <- o$r; step <- o$step; t <- o$t
    if (it > burn) { keep[it - burn, ] <- r; accepted[it - burn, ] <- o$accepted }
  }
  list(draws = keep, acceptance = colMeans(accepted), step = step, t = t)
}

test_that("the size update recovers the size from counts with a known mean", {
  set.seed(20261125)
  N <- 400
  centre <- c(2, 2, 0, 0, -2, -2)
  truth <- c(0.5, 5, 0.5, 5, 0.5, 5)
  log_mean <- sapply(centre, function(m) stats::rnorm(N, m, 1))
  y <- sapply(1:6, function(s) stats::rnbinom(N, size = truth[s], mu = exp(log_mean[, s])))
  fit <- size_run(y, log_mean)
  q <- apply(fit$draws, 2, stats::quantile, c(0.005, 0.5, 0.995))
  expect_true(all(truth >= q[1, ] & truth <= q[3, ]))
  expect_true(all(abs(log(q[2, 1:4] / truth[1:4])) < log(1.5)))
  expect_true(all(fit$acceptance > 0.25 & fit$acceptance < 0.65))
  expect_identical(fit$t, 1000L)
})

test_that("the size update ignores missing and excluded cells and survives a species never counted", {
  set.seed(1)
  N <- 60
  log_mean <- matrix(stats::rnorm(N * 3), N, 3)
  y <- matrix(stats::rnbinom(N * 3, size = 2, mu = exp(log_mean)), N, 3)
  y[, 3] <- 0
  y[1:10, 1] <- NA
  include <- matrix(TRUE, N, 3)
  include[11:20, 2] <- FALSE
  changed <- y
  changed[11:20, 2] <- 999
  set.seed(5); a <- occJSDM:::count_size_step(rep(1, 3), y, log_mean, include, 1, rep(0.5, 3), TRUE, 0L)
  set.seed(5); b <- occJSDM:::count_size_step(rep(1, 3), changed, log_mean, include, 1, rep(0.5, 3), TRUE, 0L)
  expect_identical(a, b)
  expect_true(all(is.finite(a$r)) && all(a$r > 0))
  expect_identical(a$t, 1L)
  fixed <- occJSDM:::count_size_step(rep(1, 3), y, log_mean, include, 1, rep(0.5, 3), FALSE, 7L)
  expect_identical(fixed$step, rep(0.5, 3))
  expect_identical(fixed$t, 7L)
  expect_error(occJSDM:::count_size_step(rep(1, 2), y, log_mean, include, 1, rep(0.5, 3), TRUE, 0L),
               "one size per species")
})
```

- [ ] **Step 2: Run and see them fail** with "object not found".

- [ ] **Step 3: Write `R/count_pieces.R`.** Open the file with a comment block in the style of `R/b0_shift.R` saying what the update is and that stage 2 calls it. Fixed by the plan:

  - Each call draws `stats::rnorm(S)` and then `stats::runif(S)`, in that order and nothing else, so that a run is reproducible from `set.seed()`.
  - The proposal for species `s` is `log(r[s]) + step[s]` times its normal draw. A proposal outside `count_size_bounds` is rejected.
  - The log acceptance ratio is the difference in the species' summed log density, over cells that are included and not `NA`, plus the difference in `count_size_logprior()`. Accept when the log of the uniform draw is below it and it is finite.
  - `t` is an integer count of adapting calls. With `adapt = TRUE`, add `1L` to `t`, then move each species' log step by its acceptance probability minus the target, divided by the square root of `t`, and keep the step between `1e-3` and `1e2`. With `adapt = FALSE` leave `step` and `t` alone. This is the recursion in `b0_shift_step()`.
  - Stop with "one size per species" unless `r` and `step` have one value per column of `y`, and with "y, log_mean and include must have the same dimensions" otherwise.
  - `count_log_mean()` stops with a message naming `idx_site` when an index is outside the rows of `eta`, and naming `log_effort` when it is not finite or not one per sample.

- [ ] **Step 4: Run the test file.** Expected: PASS, in a few seconds.

- [ ] **Step 5: Commit**

```bash
git add R/count_pieces.R tests/testthat/test-count-pieces.R
git commit -m "Add effort offsets and the per-species size update for counts"
```

### Task 5: The stage 1 measurements

**Files:**
- Create: `dev/count-model/stage1-checks.R`
- Create: `dev/count-model/stage1-report.md` (written by hand from the script's output)

**Interfaces:**
- Consumes: `samplePGcount_cpp()`, `count_size_step()`, `count_log_mean()`.
- Produces: the two measurements the spec asks stage 1 for, and a recommendation on the saddlepoint sampler.

- [ ] **Step 1: Write the script.** One process, no parallel workers. It starts with `devtools::load_all(quiet = TRUE)`, `set.seed(20261125)` and `setOccJSDMSeed(20261125L)`, and prints:

  - The power state: the output of `pmset -g batt | head -n 1` and of `pmset -g | grep lowpowermode`, and `R.version.string`.
  - **Time per draw by shape.** For shapes 0.5, 1.5, 5.5, 13.5, 50.5, 120.5, 169.5 and 171, at predictors 0 and -3: seconds for 200,000 draws.
  - **Where the time goes on a realistic matrix.** Simulate 250 sites by 300 species, one sample each. Species' centres are normal with mean -1 and standard deviation 2, each site's log mean is normal around its species' centre with standard deviation 1.5, and sizes are log-normal with log mean 0 and log standard deviation 1. With `shape = y + r` and `phi` the log mean minus the log size, report for each route (extreme predictor; shape up to 13; shape from 13 to 170; shape above 170) the share of cells and the share of the total draw time. Time each route by drawing its cells alone, 20 repeats.
  - **Size recovery by rarity.** 60 species by 250 sites with a known log mean and true sizes 0.5 and 5 alternating, centres evenly spaced from -4 to 3. Run 5,000 iterations of `count_size_step()`, the first 2,000 adapting. For species grouped by their number of non-zero counts (under 10, 10 to 49, 50 to 149, 150 and over) report the median error of the posterior median on the log scale, the share of 95% intervals that cover the truth, and the median ratio of the interval's upper to lower limit.

- [ ] **Step 2: Run it once, in the foreground, and do not edit the script while it runs.**

Run: `Rscript dev/count-model/stage1-checks.R 2>&1 | tee /tmp/stage1-checks.log`
Expected: all three blocks print, in a few minutes.

- [ ] **Step 3: Write `dev/count-model/stage1-report.md`** from the log, as bullets with no tables. It states the power state, the three sets of numbers, and one recommendation, using this rule: if the route for shapes from 13 to 170 takes under a quarter of the count sampler's draw time on the realistic matrix, no saddlepoint sampler is needed; if it takes more, the first thing to try is the 20-term series for the whole shape, which is already built, and a saddlepoint sampler only after that. It also states plainly how well the size is recovered for species with fewer than 10 non-zero counts.

- [ ] **Step 4: Commit**

```bash
git add dev/count-model/stage1-checks.R dev/count-model/stage1-report.md
git commit -m "Measure the count sampler's cost by shape and size recovery by rarity"
```

### Task 6: Confirm nothing else changed

**Files:** none modified.

- [ ] **Step 1: Run the whole suite in the background.**

Run: `Rscript -e 'devtools::test()' 2>&1 | tail -n 8`
Expected: the baseline's fail, warn and skip counts from Task 1, and its pass count plus the new expectations in `test-count-pieces.R`. Any other difference is a regression: stop and investigate.

- [ ] **Step 2: Confirm the scope of the change.**

Run: `git diff --stat origin/main`
Expected: only `R/jsdmfun.R` (three lines), `R/count_pieces.R`, `R/RcppExports.R`, `src/count.h`, `src/jsdm.cpp`, `src/RcppExports.cpp`, `tests/testthat/test-count-pieces.R` and the two files under `dev/count-model/`. `src/functions.cpp` and `samplepg_fast()` are untouched.

Run: `git diff origin/main -- src/jsdm.cpp | grep -c '^-[^-]'`
Expected: `0`. The change to `src/jsdm.cpp` only adds lines.

- [ ] **Step 3: Hand over to Doug.** Report the suite counts against the baseline, the report's recommendation, and the four decisions under *Decisions made in this plan*. Do not push. The pull request to `main` is opened only when Doug says so; its description opens with a line saying Claude wrote it and lists every question that is still open in full.

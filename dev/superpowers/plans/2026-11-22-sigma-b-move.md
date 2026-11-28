# Scale move for `sigma_b`: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a scale-expansion move for `sigma_b`, the measured group's scale under the grouped prior, so that a chain whose `sigma_b` falls near zero comes back within tens to a few hundred sweeps.

**Architecture:** A new function `sigma_b_move()` beside `sigma_emb_move()` in `R/coef-prior.R`, built the same way, with two proposals per sweep (an adapted step and a fixed large one). It runs first of the three scale moves inside the block of `update_jSDMcoef()` that already runs the other two, and its state and acceptance are carried in the chain's existing `axis_move` record. A simulation on fully simulated data, run before and after, says whether the low state is reproduced and whether the move shortens it.

**Tech Stack:** Base R, roxygen2 8.1.0, testthat edition 3, coda (already a dependency) for R-hat in the study script.

**Spec:** `dev/superpowers/specs/2026-11-22-sigma-b-move-design.md`. Read it first, in particular "Mathematics", "Where it goes" and "Evidence before it is kept". This plan does not repeat the formulas.

## Global Constraints

- This plan changes the sampler. Do not start Task 1 until Doug has said to build. His approval of the spec covered writing this plan only.
- Branch `claude/sigma-b-move`, worktree `~/src/occJSDM/.worktrees/sigma-b-move`, cut from `claude/embedding-gradients` at a2a64a1. Do not touch any other branch or worktree. Do not push, and do not open a PR, until Doug says so.
- The target density is unchanged. No prior, default or other sampler step is changed. A fit under the fixed prior, and a grouped fit of a model type other than "binary", must draw exactly the random numbers it drew at a2a64a1.
- Existing tests are not weakened. Three edits to existing test files are planned: the helper `axis_move_log_target()` gains one term (Task 2), `sigma_move_fit()` moves unchanged into the helper file (Task 4), and one test is appended to `test-sigma-emb-hold.R` (Task 4). Any other edit to an existing test is a ruling to record.
- Never `devtools::load_all()` and never `R CMD build` in a checkout. R is 4.6. Install and test with the loop below, with `SP` set to a scratch directory outside the repository.
- At most 5 fitting workers on the laptop in total, counting the test runs. While a simulation of Task 1 or Task 7 is running on 3 workers, tests use at most 2. Long runs are detached with `nohup`. Check `pmset -g` first and record power source and Low Power Mode with any timing. Never edit an R script while `Rscript` jobs started from it are running.
- No dataset content anywhere. Every fit in this plan is on data from `simulateOccJSDMData()`.
- Notation in comments and help: the embedding block is `X_emb` and its coefficients `B_emb`.
- Markdown (`NEWS.md`, this plan, the spec, the study's record): one paragraph per line, no pipe tables, no em-dashes, and in prose outside code write `-\>`, `\~`, `\<`, `\>`, `\|`. Re-read a markdown file before and after editing it.
- Help text is plain language in the style of the `grouped-prior` page.
- Commit messages end with `Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>`.

The loop, referred to below as "install", "document" and "test FILTER". `SP` is exported, since Task 5 reads it from R:

```bash
# install
mkdir -p $SP/build $SP/lib
git ls-files -co --exclude-standard -z | rsync -a --from0 --files-from=- ./ $SP/build/
R_PROFILE_USER=/dev/null R CMD INSTALL --no-docs --library=$SP/lib $SP/build > $SP/install.log 2>&1 || tail -30 $SP/install.log

# document (after install; then install again so the new help is in the library)
R_PROFILE_USER=/dev/null R_LIBS=$SP/lib Rscript -e 'roxygen2::roxygenise(".", load_code = roxygen2::load_installed)'

# test FILTER
(cd tests && R_PROFILE_USER=/dev/null R_LIBS=$SP/lib NOT_CRAN=true Rscript -e 'testthat::test_dir("testthat", package = "occJSDM", load_package = "installed", filter = "FILTER", stop_on_failure = TRUE)')
```

The baseline library, built once and used by Tasks 1 and 5:

```bash
mkdir -p $SP/base-build $SP/lib-base
git archive a2a64a1 | tar -x -C $SP/base-build
R_PROFILE_USER=/dev/null R CMD INSTALL --no-docs --library=$SP/lib-base $SP/base-build > $SP/install-base.log 2>&1 || tail -30 $SP/install-base.log
```

Before Task 1, commit the spec and this plan:

```bash
git add dev/superpowers/specs/2026-11-22-sigma-b-move-design.md dev/superpowers/plans/2026-11-22-sigma-b-move.md
git commit -m "Add the design and plan for a scale move on sigma_b"
```

## Review Focus

Inputs and conditions the spec implies that are most likely to bite. Each has a test in the task named.

- `sigma_b` and the measured leftovers so small that the likelihood terms are exactly zero. The move must return a finite ratio and still be able to accept upward. Task 2.
- A grouped binary fit with no measured column. The move must change nothing and draw nothing. Task 2.
- The staged burn-in, while `sigma_emb` is held. The `sigma_b` move must still run, and the `sigma_emb` move must still be skipped. Task 4.
- Chains run on separate cores (`ncores`). The record must have one step and one pair of acceptance rates per chain. Task 4.
- A first proposal that is accepted, followed by a second. The second must start from the state the first left, including the leftovers handed on to `sigma_emb_move()`. Tasks 2 and 4.

---

### Task 1: Start the baseline simulation

It runs for hours and needs only the baseline library, so it is launched first and read in Task 7.

**Files:**
- Create: `dev/simstudy/sigma-b-move-20261122/repro.R`
- Create: `dev/simstudy/sigma-b-move-20261122/summarise.R`
- Output (ignored by git): `dev/simstudy/results/sigma-b-move-20261122/`

**Interfaces:**
- Produces: `repro.R --lib PATH --label LABEL --datasets 1:5 --workers N`, which writes `LABEL-d<k>.rds` per data set, a list with `sigma_b` (kept sweeps by chains), `rhat` (see step 1), `axis_move` (the fit's `infos$axis_move`), `minutes` and `power`. `summarise.R --label LABEL` prints one line per chain and one per fit.

- [ ] **Step 1: Write `repro.R`**

For data set k: `set.seed(1000L + k)`, then `simulateOccJSDMData()` with `model = "binary"` and

```r
datasettings <- list(n = 256L, S = 48L, g = 2L, M = rep(1L, 256L), P = 1L,
                     K = rep(1L, 256L), ncov_psi = 5L, ncov_theta = 1L)
jsdmParams <- list(gt = 4L, d = 2L, ds = 0, sigma_b = 0.1, sigma_bs = 0.5, sigma_ts = 0.5,
                   sigma_h = 1, sigma_s = 0.5, l_s = 0.3, tau = rep(1, 48L),
                   useSpatField = FALSE, embedding = list(n_emb = 64L))
```

Then `set.seed(2000L + k)`, `occJSDM:::setOccJSDMSeed(2000L + k)`, and `runOccJSDM()` with `listParams = list(n_factors = 2L, n_lattrait = 4L, embeddingVars = dat$embeddingVars)`, the measured columns then the embedding columns as `occCovariates` (as `grouped_occ_covariates()` in `tests/testthat/helper-grouped-prior.R` orders them), and `MCMCparams = list(nchain = 4L, nburn = 5000L, niter = 20000L, nthin = 1L)`. Threads are set to 1 (`RCPP_PARALLEL_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, `OMP_NUM_THREADS`). Data sets run through `parallel::mclapply()` with `mc.cores` equal to the workers. The package is loaded from `--lib` only, and the script stops if `packageDescription("occJSDM")` is not from that path.

`rhat`: for the first 5,000 kept sweeps and for all 20,000, over all 4 chains, the largest coda R-hat (unsplit, as `returnGroupedConvergence()` uses) among the measured rows of `B_output`, and the number over 1.1. The fit object is not saved.

- [ ] **Step 2: Write `summarise.R`**

A sweep is "low" when `sigma_b` is under 0.02. Per chain: the share of kept sweeps low, the longest run of low sweeps, and the number of runs of 500 or more. Per fit: the two `rhat` pairs and the minutes. Last line: the longest low run over all chains of the label.

- [ ] **Step 3: Check the script on a short run**

Run: `R_PROFILE_USER=/dev/null Rscript --vanilla dev/simstudy/sigma-b-move-20261122/repro.R --lib $SP/lib-base --label smoke --datasets 1:1 --workers 1 --nburn 50 --niter 100`
Expected: one file `smoke-d1.rds`, and `summarise.R --label smoke` prints 4 chain lines and 1 fit line with finite numbers. (`--nburn` and `--niter` exist for this check only.) Delete the smoke file afterwards.

- [ ] **Step 4: Launch, detached**

Run: `pmset -g batt; pmset -g | grep lowpowermode`, then
`nohup env R_PROFILE_USER=/dev/null Rscript --vanilla dev/simstudy/sigma-b-move-20261122/repro.R --lib $SP/lib-base --label base --datasets 1:5 --workers 3 > $SP/repro-base.log 2>&1 &`
Expected: 3 `Rscript` children fitting. A 4-chain fit of this size took 72 to 82 minutes and about 6.5 GB on the AlphaEco design, so about 3 hours for the five. Three workers, not four, for the memory.

- [ ] **Step 5: Commit**

```bash
git add dev/simstudy/sigma-b-move-20261122
git commit -m "Add the simulation that looks for chains with sigma_b near zero"
```

### Task 2: `sigma_b_move()`

**Files:**
- Modify: `R/coef-prior.R` (new function after `sigma_emb_move()`, new constant)
- Modify: `tests/testthat/helper-scale-moves.R` (`axis_move_log_target()`)
- Test: `tests/testthat/test-sigma-b-move.R` (new)

**Interfaces:**
- Produces: `sigma_b_move(B, R, sigma_b, X, kappa, Omega, eta, coef_prior, step)`. `step` is a numeric vector, one proposal per element, made in order. All `delta` are drawn first (`rnorm(length(step), 0, step)`), then all uniforms. Returns `list(B, R, sigma_b, eta, accepted, log_ratio)`, where `R` is the leftovers of the returned `B`, and `accepted` (logical) and `log_ratio` (numeric) have one element per proposal. With no measured column it returns its inputs, `accepted = rep(FALSE, length(step))`, `log_ratio = rep(NA_real_, length(step))`, and draws no random number.
- Produces: `SIGMA_B_MOVE_LARGE_STEP <- 1`.

- [ ] **Step 1: Add the prior of `sigma_b` to the test helper's target**

In `axis_move_log_target()` add `dcauchy(x$sigma_b, 0, x$cp$sigma_b_scale, log = TRUE) + log(x$sigma_b)` to the sum, and correct the comment above it, which says the prior of `sigma_b` is left out. The existing tests take differences at an unchanged `sigma_b`, so they are not affected.

Run: test `scale-move|sigma-emb`
Expected: PASS, as before the edit.

- [ ] **Step 2: Write the failing tests**

```r
sigma_b_leftovers <- function(x, B = x$B) t(B) - x$Tr_b %*% x$G - x$A %*% x$C

run_sigma_b_move <- function(x, step = 0.5) with(x, occJSDM:::sigma_b_move(
  B, sigma_b_leftovers(x), sigma_b, X, kappa, Omega, eta, cp, step))

# The state after scaling the measured leftovers and sigma_b by r.
sigma_b_scaled <- function(x, r) {
  meas <- which(!x$cp$is_emb)
  R <- sigma_b_leftovers(x); R[, meas] <- R[, meas] * r
  x$B <- t(x$Tr_b %*% x$G + x$A %*% x$C + R); x$sigma_b <- x$sigma_b * r
  x$eta <- x$rest + x$X %*% x$B
  x
}

test_that("the sigma_b move accepts on the exact ratio of the target density", {
  n_accept <- 0L; n_reject <- 0L
  for (seed in 1:40) {
    x <- axis_move_state(seed)
    m <- ncol(x$B) * sum(!x$cp$is_emb)
    step <- c(0.5, 0.05)[1L + seed %% 2L]
    seed_both(300L + seed); got <- run_sigma_b_move(x, step)
    seed_both(300L + seed); delta <- rnorm(1, 0, step); logu <- log(runif(1))
    prop <- sigma_b_scaled(x, exp(delta))
    log_ratio <- axis_move_log_target(prop) - axis_move_log_target(x) + m * delta
    expect_equal(got$log_ratio, log_ratio, tolerance = 1e-8, info = paste("seed", seed))
    expect_identical(got$accepted, logu < log_ratio, info = paste("seed", seed))
    if (logu < log_ratio) {
      expect_equal(got[c("B", "sigma_b", "eta")], prop[c("B", "sigma_b", "eta")])
      expect_equal(got$R, sigma_b_leftovers(prop))
      n_accept <- n_accept + 1L
    } else {
      expect_identical(got[c("B", "sigma_b", "eta")], x[c("B", "sigma_b", "eta")])
      n_reject <- n_reject + 1L
    }
  }
  expect_gt(n_accept, 5L); expect_gt(n_reject, 5L)
})

test_that("two proposals are made in order, the second from the state the first leaves", {
  n_first <- 0L
  for (seed in 1:30) {
    x <- axis_move_state(seed); m <- ncol(x$B) * sum(!x$cp$is_emb)
    step <- c(0.05, 1)
    seed_both(500L + seed); got <- run_sigma_b_move(x, step)
    seed_both(500L + seed); delta <- rnorm(2, 0, step); logu <- log(runif(2))
    cur <- x; want_ratio <- numeric(2); want_acc <- logical(2)
    for (j in 1:2) {
      prop <- sigma_b_scaled(cur, exp(delta[j]))
      want_ratio[j] <- axis_move_log_target(prop) - axis_move_log_target(cur) + m * delta[j]
      want_acc[j] <- logu[j] < want_ratio[j]
      if (want_acc[j]) cur <- prop
    }
    expect_equal(got$log_ratio, want_ratio, tolerance = 1e-8, info = paste("seed", seed))
    expect_identical(got$accepted, want_acc, info = paste("seed", seed))
    expect_equal(got[c("B", "sigma_b", "eta")], cur[c("B", "sigma_b", "eta")], tolerance = 1e-10)
    expect_equal(got$R, sigma_b_leftovers(cur), tolerance = 1e-10)
    n_first <- n_first + want_acc[1]
  }
  expect_gt(n_first, 5L)                               # the second proposal did start from a moved state
})

test_that("an accepted sigma_b move scales every measured leftover and nothing else", {
  x <- axis_move_state(3L); x$kappa[] <- 0; x$Omega[] <- 0
  got <- NULL
  for (s in 1:20) { seed_both(s); got <- run_sigma_b_move(x); if (got$accepted) break }
  expect_true(got$accepted)
  r <- got$sigma_b / x$sigma_b
  R0 <- sigma_b_leftovers(x); R1 <- sigma_b_leftovers(x, got$B)
  expect_equal(R1[, 1:2], r * R0[, 1:2], tolerance = 1e-10)
  expect_identical(got$B[3:5, ], x$B[3:5, ])           # embedding rows bit for bit
  expect_equal(got$eta, x$rest + x$X %*% got$B, tolerance = 1e-12)
})

test_that("the sigma_b move does nothing and draws nothing without measured columns", {
  x <- axis_move_state(4L); x$cp$is_emb[] <- TRUE
  seed_both(5L); before <- .Random.seed
  got <- run_sigma_b_move(x, c(0.5, 1))
  expect_identical(.Random.seed, before)
  expect_identical(got[c("B", "sigma_b", "eta")], x[c("B", "sigma_b", "eta")])
  expect_identical(got$accepted, c(FALSE, FALSE))
  expect_identical(got$log_ratio, c(NA_real_, NA_real_))
})

test_that("with sigma_b and its leftovers near zero the ratio is finite and upward moves are accepted", {
  x <- sigma_b_scaled(axis_move_state(6L), 1e-12)     # far below any value a chain was seen at
  up <- 0L; total <- 0L
  for (s in 1:60) {
    seed_both(700L + s); got <- run_sigma_b_move(x, 1)
    seed_both(700L + s); delta <- rnorm(1, 0, 1)
    expect_true(is.finite(got$log_ratio))
    expect_equal(got$log_ratio, delta, tolerance = 1e-6)   # the likelihood and prior terms vanish
    if (delta > 0) { total <- total + 1L; up <- up + got$accepted }
  }
  expect_identical(up, total)                          # every upward proposal accepted
})
```

- [ ] **Step 3: Run them to see them fail**

Run: install, then test `sigma-b-move`
Expected: FAIL, `sigma_b_move` is not an object of the package.

- [ ] **Step 4: Implement `sigma_b_move()` and the constant in `R/coef-prior.R`**

As `sigma_emb_move()`, with the measured columns (`which(!coef_prior$is_emb)`) in place of the embedding columns and `coef_prior$sigma_b_scale` in place of `sigma_emb_scale`, looping over the elements of `step`. After an accepted proposal, `D` and the measured columns of `R` are scaled by `r` so that the next proposal starts from the new state. A header comment in the style of the one on `sigma_emb_move()` says what it is for, cites the spec, and gives the ratio.

- [ ] **Step 5: Run the tests**

Run: install, then test `sigma-b-move|sigma-emb|axis-scale`
Expected: PASS. If the first test does not exercise both outcomes more than 5 times each, change its two step sizes and record the ruling.

- [ ] **Step 6: Commit**

```bash
git add R/coef-prior.R tests/testthat/helper-scale-moves.R tests/testthat/test-sigma-b-move.R
git commit -m "Add sigma_b_move(), a scale-expansion move on the measured group's scale"
```

### Task 3: The three moves together leave the priors unchanged

**Files:**
- Test: `tests/testthat/test-sigma-b-move.R`

**Interfaces:**
- Consumes: `sigma_b_move()` from Task 2.

- [ ] **Step 1: Write the test**

`test_that("the scale updates with all three scale moves leave the half-Cauchy priors unchanged", ...)`. It is the loop of the test of the same kind in `tests/testthat/test-sigma-emb-move.R` (lines 79 to 132), with one addition after `sample_coef_scales()` and before `sigma_emb_move()`:

```r
sb <- occJSDM:::sigma_b_move(B, t(B) - centre, sigma_b, X, zero, zero, zero, cp,
                             c(0.5, occJSDM:::SIGMA_B_MOVE_LARGE_STEP))
B <- sb$B; sigma_b <- sb$sigma_b
```

and `sigma_emb_move()` then takes `t(B) - centre` of the new `B`. Same seed, sizes, 50,000 kept sweeps, quartiles and tolerances. It also expects that more than half of the `sigma_b` proposals at each step were accepted.

- [ ] **Step 2: Run it**

Run: test `sigma-b-move`
Expected: PASS. This test is the check that the ratio is right, so it is not made to pass by loosening a tolerance. If a share misses, stop: the move is wrong, or the comment at the end of `test-sigma-emb-move.R` on how the tolerances were measured needs repeating for this loop, with seeds 1 to 8, and the result recorded as a ruling.

- [ ] **Step 3: Commit**

```bash
git add tests/testthat/test-sigma-b-move.R
git commit -m "Test that the three scale moves leave the half-Cauchy priors unchanged"
```

### Task 4: Run the move in the sweep and record it

**Files:**
- Modify: `R/coef-prior.R` (`axis_move_start()`, `axis_move_record()`, their comments)
- Modify: `R/jsdmfun.R:1535-1585` (the block that runs the moves)
- Modify: `R/chains.R:633-645`, `R/chains.R:733-739`, `R/chains.R:888-895`
- Test: `tests/testthat/test-sigma-b-move.R`, `tests/testthat/test-sigma-emb-hold.R` (one test appended)

**Interfaces:**
- Consumes: `sigma_b_move()`, `SIGMA_B_MOVE_LARGE_STEP`.
- Produces, in the chain's move state: `sigma_b_step` (starts at `AXIS_MOVE_STEP`), `sigma_b_accepted` (logical, length 2), `sigma_b_t` (integer).
- Produces, in `infos$axis_move`: `sigma_b_step` (one per chain), `sigma_b_large_step` (the constant), and `sigma_b_acceptance` and `sigma_b_acceptance_burnin`, each a matrix with rows `"adapted"` and `"large"` and one column per chain.

- [ ] **Step 1: Write the failing tests**

```r
test_that("a grouped binary fit runs the sigma_b move and records its acceptance", {
  for (mcmc in list(FIXTURE_MCMC, c(FIXTURE_MCMC, list(ncores = 2L)))) {
    fit <- sigma_move_fit(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L), mcmc)
    am <- fit$infos$axis_move
    expect_length(am$sigma_b_step, 2L)
    expect_true(all(am$sigma_b_step >= 0.01 & am$sigma_b_step <= 5))
    expect_false(all(am$sigma_b_step == 0.5))
    expect_identical(am$sigma_b_large_step, 1)
    for (nm in c("sigma_b_acceptance", "sigma_b_acceptance_burnin")) {
      expect_identical(dim(am[[nm]]), c(2L, 2L))
      expect_identical(rownames(am[[nm]]), c("adapted", "large"))
      expect_true(all(am[[nm]] >= 0 & am[[nm]] <= 1))
    }
    expect_true(all(is.finite(fit$results_output$jsdm_output$sigmab_output)))
  }
})
```

`sigma_move_fit()` is defined in `test-sigma-emb-move.R`; move it to `helper-scale-moves.R` unchanged so both files can use it.

A second test, `"the coefficient update runs the sigma_b move first and hands its state on"`, spies on `sigma_b_move()` and `sigma_emb_move()` inside `update_jSDMcoef()` exactly as the last test of `test-sigma-emb-move.R` spies on its two moves, in a fit of `nchain = 1, nburn = 3, niter = 3`. It expects, for each of the 6 sweeps: the `R` given to `sigma_b_move()` is the leftovers of the `B` given; its `eta` is `B0 + X B + U L` of that `B`; its `step` is `c(adapted, 1)`; the `B`, `R` and `eta` given to `sigma_emb_move()` are identical to those `sigma_b_move()` returned; and the sweep's returned `sigma_b` is identical to the move's. It also expects the adapted step to start at 0.5, to follow `axis_step_adapt(step, log_ratio[1], i)` through the three burn-in sweeps, and to be frozen afterwards.

A third, appended to `test-sigma-emb-hold.R`, uses that file's own way of making a chain with a hold, and expects that in every sweep where `sigma_emb` is held `sigma_b_move()` was called and `sigma_emb_move()` was not.

- [ ] **Step 2: Run them to see them fail**

Run: test `sigma-b-move|sigma-emb-hold`
Expected: FAIL, `sigma_b_step` is NULL and the spy on `sigma_b_move` is never called.

- [ ] **Step 3: Implement**

- `axis_move_start()`: the three new fields.
- `update_jSDMcoef()`: after `eta_drawn` is built, call `sigma_b_move(B, leftovers, sigma_b, X, k, Omega, eta_drawn, coef_prior, c(axis_move$sigma_b_step, SIGMA_B_MOVE_LARGE_STEP))`, take its `B`, `sigma_b` and `R` (as `leftovers`), and give its `eta` to `sigma_emb_move()`, or to the stand-in list when `sigma_emb` is held. In the adapting block, when there is a measured column, add 1 to `sigma_b_t` and adapt `sigma_b_step` from `log_ratio[1]` with `axis_step_adapt()`. The large step is never adapted.
- `R/chains.R`: a counter pair for the two proposals in each phase, as `sigma_accepted` is kept, and the new entries in the chain's returned `axis_move` list.
- `axis_move_record()`: the four new entries.
- Comments that say "two scale moves" now say three.

- [ ] **Step 4: Run the tests**

Run: install, then test `sigma-b-move|sigma-emb|axis-scale|grouped`
Expected: PASS, including the existing spy test in `test-sigma-emb-move.R`, which checks that `sigma_emb_move()` is still given the leftovers of the `B` it is given.

- [ ] **Step 5: Commit**

```bash
git add R/coef-prior.R R/jsdmfun.R R/chains.R tests/testthat
git commit -m "Run the sigma_b move in grouped binary fits and record its acceptance"
```

### Task 5: Nothing else changed

**Files:**
- Create: `dev/simstudy/sigma-b-move-20261122/identity.R`

- [ ] **Step 1: Write `identity.R`**

`identity.R --lib PATH --out FILE` loads the package from `PATH`, builds three fixtures with the simulator calls of `tests/testthat/helper-grouped-prior.R` and `helper-fixtures.R` copied into the script, and fits each with `set.seed(11L)`, `occJSDM:::setOccJSDMSeed(11L)` and `MCMCparams = list(nchain = 2, nburn = 20, niter = 20, nthin = 1)`: (1) a binary fit under the fixed prior, (2) a grouped fit of model type "occupancy", (3) a grouped binary fit. It saves the three `B_output` arrays in a list.

- [ ] **Step 2: Run it against both libraries and compare**

Run:
```bash
for L in lib-base lib; do R_PROFILE_USER=/dev/null Rscript --vanilla dev/simstudy/sigma-b-move-20261122/identity.R --lib $SP/$L --out $SP/identity-$L.rds; done
Rscript -e 'a <- readRDS(file.path(Sys.getenv("SP"), "identity-lib-base.rds")); b <- readRDS(file.path(Sys.getenv("SP"), "identity-lib.rds")); print(mapply(identical, a, b))'
```
Expected: `TRUE TRUE FALSE`. The first two say that fits without the move are unchanged bit for bit. The third says the move is running.

- [ ] **Step 3: Run the whole suite**

Run: test `.` (the empty filter)
Expected: 0 failures. The 6 skips that predate the branch remain.

- [ ] **Step 4: Commit**

```bash
git add dev/simstudy/sigma-b-move-20261122/identity.R
git commit -m "Add the check that fits without the sigma_b move are unchanged"
```

### Task 6: Help and NEWS

**Files:**
- Modify: `R/coef-prior.R` (roxygen, section "Sampling and chain length", near line 621)
- Modify: `NEWS.md` (development version, "New features")
- Modify: `man/` (by "document")

- [ ] **Step 1: Help**

In "Sampling and chain length", "three extra steps" becomes "four", and one bullet is added after the bullet on the embedding leftover scale:

```
#' * The scale of the measured covariates' species leftovers (`sigma_b`) gets
#'   the same kind of step, together with those leftovers, twice in each
#'   iteration: once with a step tuned during burn-in and once with a large
#'   fixed step. Without it a chain could stay for thousands of iterations
#'   with this scale near zero, and disagree with the other chains on the
#'   measured coefficients meanwhile.
```

- [ ] **Step 2: NEWS**

One bullet at the end of "New features" of the development version, in the file's own list style:

```
-   Under the grouped prior, a `binary` fit with an embedding takes one more sampling step, for `sigma_b`, the scale of the measured covariates' species leftovers. Without it a chain could stay for thousands of iterations with that scale near zero. The model is unchanged. A seeded fit of this kind gives different numbers from the previous development version; every other fit is unchanged.
```

- [ ] **Step 3: Document, install, test**

Run: document, install, then test `grouped|sigma-b-move`
Expected: PASS, and `git status` shows only the help file of the grouped prior changed under `man/`.

- [ ] **Step 4: Commit**

```bash
git add R/coef-prior.R NEWS.md man
git commit -m "Describe the sigma_b move in the help and NEWS"
```

### Task 7: Read the baseline, and run the same data with the move

**Files:**
- Create: `dev/simstudy/sigma-b-move-20261122/RECORD.md`

**Interfaces:**
- Consumes: the `base-d<k>.rds` files of Task 1, and `repro.R` and `summarise.R` unchanged.

- [ ] **Step 1: Read the baseline**

Run: `Rscript --vanilla dev/simstudy/sigma-b-move-20261122/summarise.R --label base`, once Task 1's run has ended (check that 5 files exist; do not rely on the exit code).
Expected, by the spec's rule: "reproduced" if any chain of any data set has a run of 500 or more low sweeps. If it is not reproduced, write that in `RECORD.md`, skip steps 2 and 3, and report to Doug. Do not change the simulated truth to look for the low state without his say.

- [ ] **Step 2: Run the same data and seeds with the move**

Run: `pmset` as in Task 1, then
`nohup env R_PROFILE_USER=/dev/null Rscript --vanilla dev/simstudy/sigma-b-move-20261122/repro.R --lib $SP/lib --label move --datasets 1:5 --workers 3 > $SP/repro-move.log 2>&1 &`
Expected: 5 files `move-d<k>.rds`, about 3 hours.

- [ ] **Step 3: Read it**

Run: `summarise.R --label move`
Expected, by the spec's rule: "stays are short" if the longest low run in any chain is under a fifth of the baseline's longest and none exceeds 500 sweeps. Otherwise the result is "stays are not short", and it is reported as that.

- [ ] **Step 4: Write `RECORD.md`**

The commit of each library, the power state of each run, the two summaries in full, the two readings in the spec's words, the share of kept sweeps low in each (expected not to fall), the two `rhat` pairs per fit for both labels, and the acceptance rates of the adapted and large proposals. Nothing is claimed beyond the two readings.

- [ ] **Step 5: Commit**

```bash
git add dev/simstudy/sigma-b-move-20261122/RECORD.md
git commit -m "Record the simulation of the sigma_b move"
```

---

## After the last task

- A whole-branch review by a fresh reviewer, as the executing skill describes, with the Review Focus list above.
- Report to Doug: the two readings, the rulings, and the deferred minors. Do not push.
- The third piece of evidence in the spec, the replay by the AlphaEco project with fresh seeds, is not a task of this plan. It needs a request that Doug reads first, and a library built from this branch's head.
- If the move is kept, the pinned commit of version 2 of the simulation contract (`claude/validation-contracts`) must move before the contract is issued. That is Doug's to say.

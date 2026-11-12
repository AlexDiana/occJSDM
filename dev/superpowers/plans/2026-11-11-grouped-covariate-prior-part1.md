# Grouped Covariate Prior, Part 1 (Covariate Side) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in prior, switched on by naming embedding columns in `listParams$embeddingVars` (or by `listPriors$b_prior = "grouped"` when there is no embedding), in which measured covariates and embedding columns get separate learned scales, each embedding column gets its own learned scale, and every covariate gets a community-average response.

**Architecture:** The scalar `sigma_b` that sets the prior standard deviation of every occupancy coefficient becomes a vector with one entry per covariate, threaded through the joint coefficient draw, both spatial range functions and the two trait-regression updates. Three new conjugate updates (two group scales and the axis scales) replace `sample_sigmab()` when the option is on. With the option off the vector is constant and every existing fit is reproduced bit for bit.

**Tech Stack:** R 4.6.1, Rcpp and RcppArmadillo (`src/jsdm.cpp`), testthat edition 3, devtools.

**Spec:** `dev/superpowers/specs/2026-11-11-grouped-covariate-prior-design.md`. Read it before starting. This plan covers the sections *The model, part 1*, *Interface*, *Outputs*, *Simulator* and the unit tests under *Staging and tests*. It does not cover taxonomy (part 2), the real embedding extract or the validation study; each of those gets its own plan after this one is merged into the branch.

## Global Constraints

- Work only in the worktree `.worktrees/grouped-covariate-prior` on branch `claude/grouped-covariate-prior`. Never switch branches in the main checkout.
- `b_prior` defaults to `"fixed"`. No default changes anywhere.
- With `b_prior = "fixed"`, fits are bit-identical to commit `4f2ee6c` for `binary`, `continuous`, `occupancy` and `two_stage`. Task 1's test enforces this and must pass at the end of every task.
- The feature works for `binary`, `occupancy` and `two_stage` fits, spatial or not, serial or parallel chains.
- No change to the spatial coefficients (`Bs`, `Gs`, `As`, `Cs`, `sigma_bs`), the factor model, the detection stages or the intercept prior.
- Naming an embedding column in `splineVars` is an error.
- Default prior scales: `sigma_b_scale = 1`, `sigma_emb_scale = 1`, `axis_scale = 1 / sqrt(number of embedding columns)`. On embedding columns the base prior standard deviation of the entries of `C` is `1 / sqrt(gt)`; on measured columns it stays 1.
- Run tests with `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/<file>")'` from the worktree root. `Rscript` must be R 4.6.1.
- Markdown follows `CLAUDE.md`: one line per paragraph, no em-dashes. Commit messages use no em-dashes.
- Do not run `R CMD build` in a checkout that has an untracked `dev/` tree; use `devtools::test()` and `devtools::document()`.

## Review Focus

Inputs the spec implies but does not spell out, most likely first. Each has a test in the task named.

1. More embedding columns than sites (for example 128 columns on 100 sites). The fit must run and return finite draws. Test in Task 8.
2. `n_lattrait = 0` with embedding columns. There are no blends, so the `1 / sqrt(gt)` base must not divide by zero and the fit must run. Tests in Tasks 3 and 8.
3. An embedding column that is constant across sites, or `embeddingVars` naming a factor, a spline variable or a name not in `occCovariates`. Each must stop with a message naming the offending columns, not fail later inside the sampler. Tests in Tasks 3 and 6.
4. A single species with the grouped prior, whether reached through `embeddingVars` or `b_prior = "grouped"`. The scales cannot be learned; it must stop with a clear message, as `b0_prior = "community"` does. Test in Task 3.
5. Downstream functions on a grouped fit: `predictNewSites()` with new embedding values, `thinOutput()`, `returnTraitsCoeff()` (whose `G_output` now has one more row) and `computeDiagnostics()`. Each must work without an index error. Tests in Task 7.

---

## File Structure

- `R/coef-prior.R` (new): `read_coef_prior()`, `coef_prior_sds()`, `sample_scale_half_cauchy()`, `sample_coef_scales()`. Everything that is specific to the grouped prior and is not a change to an existing function.
- `R/jsdmfun.R` (modify): `computePsiCoef()`, `sample_BBsL()`, `sample_GC()`, `sample_A()`, both range functions, `update_jSDMcoef()`, `simulateData()`.
- `src/jsdm.cpp` (modify): `sample_BBsL_cpp()`, `sample_BBsL_parallel()`.
- `R/runOccJSDM.R`, `R/chains.R` (modify): reading the option, chain state, storage, combination.
- `R/output.R`, `R/diagnostics.R` (modify): accessors and block selection.
- `R/simulateData.R` (modify): the embedding block.
- `tests/testthat/helper-grouped-prior.R`, `tests/testthat/test-grouped-prior.R`, `tests/testthat/test-grouped-prior-fits.R`, `tests/testthat/fixtures/` (new).

---

### Task 1: Reference fits for the fixed prior

**Files:**
- Create: `tests/testthat/helper-grouped-prior.R`
- Create: `tests/testthat/fixtures/make-fixed-prior-reference.R`
- Create: `tests/testthat/fixtures/fixed-prior-reference.rds`
- Test: `tests/testthat/test-grouped-prior.R`

**Interfaces:**
- Produces: `fixed_prior_reference_fit(model)`, returning a named list of the arrays `B0_output`, `B_output`, `G_output`, `A_output`, `C_output`, `L_output`, `sigmab_output` from `jsdm_output` of a small fit. `seed_both(s)`, which seeds R's generator and the C++ stream.

- [ ] **Step 1: Confirm the base.** Run `git log -1 --format=%h -- R src` and confirm no commit after `4f2ee6c` touches `R/` or `src/`. If one does, stop and report.

- [ ] **Step 2: Write the helper.**

```r
# tests/testthat/helper-grouped-prior.R
seed_both <- function(s) { set.seed(s); occJSDM:::setOccJSDMSeed(as.integer(s)) }

fixed_prior_reference_fit <- function(model) {
  i <- match(model, c("binary", "continuous", "occupancy", "two_stage"))
  spatial <- model %in% c("binary", "two_stage")   # covers both range functions
  sim <- simulate_fixture(model = model, useSpatField = spatial, seed = 100L + i)
  set.seed(200L + i)
  fit <- fit_fixture(sim, spatial = spatial)
  fit$results_output$jsdm_output[c("B0_output", "B_output", "G_output", "A_output",
                                   "C_output", "L_output", "sigmab_output")]
}
```

- [ ] **Step 3: Write the generator** `tests/testthat/fixtures/make-fixed-prior-reference.R`. It calls `devtools::load_all()`, sources the test helpers with `testthat::source_test_helpers("tests/testthat")`, builds `fits` as a named list over the four models with `fixed_prior_reference_fit()`, and saves `list(fits = fits, la = La_library(), r = R.version.string, commit = "4f2ee6c")` to `tests/testthat/fixtures/fixed-prior-reference.rds` with `saveRDS(..., compress = "xz")`.

- [ ] **Step 4: Generate the reference.** Run `Rscript tests/testthat/fixtures/make-fixed-prior-reference.R`. Expected: the file exists and is under 2 MB. If it is larger, reduce what is kept to `B_output`, `B0_output` and `sigmab_output`.

- [ ] **Step 5: Write the test.**

```r
# tests/testthat/test-grouped-prior.R
test_that("b_prior = 'fixed' reproduces fits made before the option", {
  ref <- readRDS(test_path("fixtures", "fixed-prior-reference.rds"))
  skip_if_not(identical(La_library(), ref$la), "BLAS or LAPACK differs from the reference")
  for (model in names(ref$fits)) {
    expect_identical(fixed_prior_reference_fit(model), ref$fits[[model]], info = model)
  }
})
```

- [ ] **Step 6: Run it.** Run the test file. Expected: PASS, 4 expectations, not skipped. A skip here means the reference was built under a different BLAS; rebuild it.

- [ ] **Step 7: Commit.**

```bash
git add tests/testthat/helper-grouped-prior.R tests/testthat/fixtures tests/testthat/test-grouped-prior.R
git commit -m "Add reference fits that pin the fixed coefficient prior"
```

---

### Task 2: Embedding block in the simulator

**Files:**
- Modify: `R/simulateData.R` (`simulateOccJSDMData()`, roxygen)
- Modify: `R/jsdmfun.R:620-830` (`simulateData()`)
- Modify: `tests/testthat/helper-grouped-prior.R`
- Test: `tests/testthat/test-grouped-prior.R`

**Interfaces:**
- Produces: `list_jsdmParams$embedding`, `NULL` by default (no change to any existing call) or a list with these fields:
  - `n_emb`: number of embedding columns, added after the `ncov_psi` measured columns.
  - `X`: optional `n` by `n_emb` matrix used as the embedding columns as given. When absent, columns are simulated as below.
  - `B`: optional `n_emb` by `S` matrix of true coefficients, used as given. When absent, coefficients are built from the next four fields.
  - `n_blends` (default 3), `n_active` (default `min(10, n_emb)`), `blend_sd` (default 0.5), `sigma_emb` (default 0.1).
- Produces: in the returned `data_list`, embedding columns named `X_emb.1`, `X_emb.2` and so on beside the measured columns, and `data_list$embeddingVars` holding those names. In `true_params$jsdmParams_true`, `B` covers all columns in design order and `covariate_group` is a character vector of `"measured"` and `"embedding"`. Measured columns keep their existing names, which start with `X_psi`.
- Produces: `simulate_grouped_fixture(model = "binary", n = 60L, S = 12L, p_meas = 2L, p_emb = 8L, gt = 2L, seed = 1L, embedding = list())` in the helper, which calls `simulateOccJSDMData()` with no spatial field, `d = 1`, and `embedding = modifyList(list(n_emb = p_emb), embedding)`. It must work with `gt = 0`.
- Produces: `grouped_occ_covariates(dat)` in the helper, returning the measured covariate names followed by `dat$embeddingVars`: `c(setdiff(grep("^X_psi", names(dat$info), value = TRUE), dat$embeddingVars), dat$embeddingVars)`.

- [ ] **Step 1: Write the failing tests.**

```r
test_that("the simulator adds an embedding block without touching existing calls", {
  a <- simulate_fixture(model = "binary", useSpatField = FALSE, seed = 3L)
  b <- simulate_fixture(model = "binary", useSpatField = FALSE, seed = 3L)
  expect_identical(a, b)
  expect_null(a$data_list$embeddingVars)

  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  expect_identical(sim$data_list$embeddingVars, paste0("X_emb.", 1:8))
  expect_true(all(sim$data_list$embeddingVars %in% names(sim$data_list$info)))
  expect_identical(sim$true_params$jsdmParams_true$covariate_group,
                   rep(c("measured", "embedding"), c(2L, 8L)))
  expect_no_error(simulate_grouped_fixture(gt = 0L))   # no latent traits must still simulate
  expect_identical(dim(sim$true_params$jsdmParams_true$B), c(10L, 12L))
})

test_that("the embedding truth is low rank on the active axes and zero elsewhere", {
  sim <- simulate_grouped_fixture(p_emb = 8L, S = 12L,
    embedding = list(n_blends = 2L, n_active = 3L, sigma_emb = 0))
  Bemb <- sim$true_params$jsdmParams_true$B[3:10, ]
  active <- rowSums(abs(Bemb)) > 0
  expect_identical(sum(active), 3L)
  expect_lte(qr(Bemb)$rank, 2L)
})

test_that("a supplied embedding matrix and coefficient matrix are used as given", {
  X <- matrix(rnorm(60 * 8), 60, 8); B <- matrix(0, 8, 12); B[1, ] <- 1
  sim <- simulate_grouped_fixture(n = 60L, p_emb = 8L, S = 12L, embedding = list(X = X, B = B))
  expect_equal(unname(as.matrix(unique(sim$data_list$info[, paste0("X_emb.", 1:8)]))), X)
  expect_identical(unname(sim$true_params$jsdmParams_true$B[3:10, ]), B)
  expect_error(simulate_grouped_fixture(n = 60L, p_emb = 8L, embedding = list(X = X[1:10, ])),
               "embedding\\$X must have one row per site and n_emb columns")
})
```

- [ ] **Step 2: Run them.** Expected: FAIL, `simulate_grouped_fixture` not found.

- [ ] **Step 3: Implement.** In `simulateData()`, generate the measured block exactly as now, so the random stream of existing calls is untouched; draw everything for the embedding block after all existing draws. Simulated columns use a factor model with a decaying spectrum:

```r
Fx <- min(8L, n_emb)
W <- matrix(rnorm(n_emb * Fx), n_emb, Fx) %*% diag(1 / seq_len(Fx), Fx)
X_emb <- scale(matrix(rnorm(n * Fx), n, Fx) %*% t(W) + 0.3 * matrix(rnorm(n * n_emb), n, n_emb))
```

Built coefficients: choose `n_active` columns at random; `A_emb` is `S` by `n_blends` standard Normal; `C_emb` is `n_blends` by `n_emb`, Normal with standard deviation `blend_sd` on the active columns and zero elsewhere; the leftover is Normal with standard deviation `sigma_emb` on the active columns and zero elsewhere; the block is `t(A_emb %*% C_emb + leftover)`. The embedding block adds `X_emb %*% B_emb` to the linear predictor.

- [ ] **Step 4: Run the tests and Task 1's test.** Expected: PASS.

- [ ] **Step 5: Document** the `embedding` field in `?simulateOccJSDMData`, run `devtools::document()`, and commit.

```bash
git add R/simulateData.R R/jsdmfun.R man tests/testthat
git commit -m "Let the simulator add an embedding block of occupancy covariates"
```

---

### Task 3: Reading the option

**Files:**
- Create: `R/coef-prior.R`
- Test: `tests/testthat/test-grouped-prior.R`

**Interfaces:**
- The switch, decided by Doug on 12 November 2026: naming `embeddingVars` is enough to turn the grouped prior on. The type is resolved as follows. With `b_prior` absent: `"grouped"` if `embeddingVars` is non-empty, otherwise `"fixed"`. With `b_prior = "grouped"`: grouped, with or without `embeddingVars`. With `b_prior = "fixed"` and a non-empty `embeddingVars`: an error. Embedding columns listed in `occCovariates` but not in `embeddingVars` are ordinary measured covariates, which is how the study's *raw* arm is fitted.
- Produces: `read_coef_prior(priors, params, S, gt, columns)`. `columns` is a data frame with one row per column of the occupancy design matrix and the fields `column` (its name), `variable` (the `occCovariates` entry it came from) and `kind` (`"numeric"`, `"factor"` or `"spline"`). Returns `list(type = "fixed")`, or `list(type = "grouped", is_emb, sigma_b_scale, sigma_emb_scale, axis_scale, c_base)`, where `is_emb` is a logical vector over design columns and `c_base` is a numeric vector over design columns: 1 on measured columns, and on embedding columns `1 / sqrt(gt)`, or 1 when `gt` is 0. `sigma_emb_scale` and `axis_scale` are `NA_real_` when there are no embedding columns.
- Produces: `coef_prior_sds(coef_prior, sigma_b, sigma_emb, lambda, g_b, p)`, returning `list(b, g, c)`: `b` and `c` are length `p`, `g` is `g_b` by `p`. Fixed: `b = rep(sigma_b, p)`, `g` all 1, `c` all 1. Grouped: `b = ifelse(is_emb, sigma_emb, sigma_b) * lambda`, `g = matrix(lambda, g_b, p, byrow = TRUE)`, `c = c_base * lambda`.

- [ ] **Step 1: Write the failing tests.**

```r
cols <- function(p_meas = 2L, p_emb = 3L, kind = "numeric") data.frame(
  column = c(paste0("m", seq_len(p_meas)), paste0("e", seq_len(p_emb))),
  variable = c(paste0("m", seq_len(p_meas)), paste0("e", seq_len(p_emb))),
  kind = kind)
rcp <- function(priors = list(), params = list(), S = 10L, gt = 2L, columns = cols())
  occJSDM:::read_coef_prior(priors, params, S, gt, columns)

test_that("read_coef_prior resolves defaults", {
  expect_identical(rcp(), list(type = "fixed"))
  expect_identical(rcp(NULL, NULL), list(type = "fixed"))
  got <- rcp(list(), list(embeddingVars = paste0("e", 1:3)))   # embeddingVars alone switches it on
  expect_identical(got$type, "grouped")
  expect_identical(rcp(list(b_prior = "grouped"), list(embeddingVars = paste0("e", 1:3))), got)
  expect_identical(rcp(list(sigma_b_scale = 2), list(embeddingVars = "e1"))$sigma_b_scale, 2)
  expect_identical(got$is_emb, c(FALSE, FALSE, TRUE, TRUE, TRUE))
  expect_identical(got$sigma_b_scale, 1)
  expect_identical(got$sigma_emb_scale, 1)
  expect_equal(got$axis_scale, 1 / sqrt(3))
  expect_equal(got$c_base, c(1, 1, rep(1 / sqrt(2), 3)))
  expect_equal(rcp(list(b_prior = "grouped"), list(embeddingVars = "e1"), gt = 0L)$c_base,
               rep(1, 5))
  none <- rcp(list(b_prior = "grouped"))
  expect_false(any(none$is_emb)); expect_true(is.na(none$axis_scale))
})

test_that("read_coef_prior rejects settings that do not apply", {
  g <- list(b_prior = "grouped"); e <- list(embeddingVars = "e1")
  for (bad in list("Grouped", NA_character_, c("fixed", "grouped"), 1, TRUE))
    expect_error(rcp(list(b_prior = bad)), "b_prior must be 'fixed' or 'grouped'")
  expect_error(rcp(list(b_prior = "fixed"), e), "embeddingVars cannot be used with b_prior = 'fixed'")
  expect_error(rcp(g, e, S = 1L), "the grouped prior needs more than one species")
  expect_error(rcp(list(), e, S = 1L), "the grouped prior needs more than one species")
  expect_error(rcp(g, list(embeddingVars = c("e1", "zz"))),
               "embeddingVars must name occupancy covariates: zz")
  expect_error(rcp(g, e, columns = cols(kind = c("numeric", "numeric", "spline", "numeric", "numeric"))),
               "embeddingVars must be numeric covariates without splines: e1")
  expect_error(rcp(g, e, columns = cols(kind = c("numeric", "numeric", "factor", "numeric", "numeric"))),
               "embeddingVars must be numeric covariates without splines: e1")
  expect_error(rcp(list(sigma_b_scale = 2)), "sigma_b_scale requires the grouped prior")
  expect_error(rcp(c(g, axis_scale = 1)), "axis_scale requires embeddingVars")
  expect_error(rcp(c(g, sigma_emb_scale = 1)), "sigma_emb_scale requires embeddingVars")
  for (nm in c("sigma_b_scale", "sigma_emb_scale", "axis_scale")) for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "a"))
    expect_error(rcp(c(g, stats::setNames(list(bad), nm)), e),
                 paste(nm, "must be a finite positive number"))
})

test_that("coef_prior_sds builds the per-covariate prior scales", {
  f <- occJSDM:::coef_prior_sds(list(type = "fixed"), 0.3, NA, NULL, g_b = 2L, p = 4L)
  expect_identical(f, list(b = rep(0.3, 4), g = matrix(1, 2, 4), c = rep(1, 4)))
  cp <- rcp(list(b_prior = "grouped"), list(embeddingVars = paste0("e", 1:3)))
  lam <- c(1, 1, 0.5, 2, 0.1)
  s <- occJSDM:::coef_prior_sds(cp, 0.3, 0.2, lam, g_b = 2L, p = 5L)
  expect_equal(s$b, c(0.3, 0.3, 0.2 * 0.5, 0.2 * 2, 0.2 * 0.1))
  expect_equal(s$g, matrix(lam, 2, 5, byrow = TRUE))
  expect_equal(s$c, cp$c_base * lam)
})
```

- [ ] **Step 2: Run them.** Expected: FAIL, `read_coef_prior` not found.

- [ ] **Step 3: Implement both functions in `R/coef-prior.R`.** Follow `read_intercept_prior()` in `R/runOccJSDM.R` for style: `get_param()` for defaults, `unname()` on returned numbers, the error messages exactly as in the tests. Resolve the type first, by the rule under *Interfaces*: validate an explicit `b_prior`, then apply the `embeddingVars` rule. Then check species count, then names, then kinds, then the scale settings.

- [ ] **Step 4: Run the tests and Task 1's test.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add R/coef-prior.R tests/testthat/test-grouped-prior.R
git commit -m "Read and validate the grouped coefficient prior settings"
```

---

### Task 4: A prior scale per covariate in the joint draw and the range step

**Files:**
- Modify: `src/jsdm.cpp:841-960` (`sample_BBsL_cpp()`) and `src/jsdm.cpp:1063-1147` (`sample_BBsL_parallel()`)
- Modify: `R/jsdmfun.R:839` (`computePsiCoef()`), `:950` (`sample_BBsL()`), `:1120` and `:1179` (range functions)
- Test: `tests/testthat/test-grouped-prior.R`

**Interfaces:**
- Produces: in `sample_BBsL_cpp()` and `sample_BBsL_parallel()`, the eighth argument becomes `arma::vec sd_b` in place of `double sigma_b`: one prior standard deviation per covariate, with a length-1 value recycled to all covariates. A new last argument `Rcpp::Nullable<Rcpp::NumericMatrix> Tr_b = R_NilValue` is the trait matrix used for the covariate side (the prior mean of `B` and the returned `Bt`); when `NULL`, `Tr` is used. `Tr` alone still serves the spatial side.
- Produces: the R reference `sample_BBsL()` takes the same two changes, with `Tr_b = Tr` as its default.
- Produces: `computePsiCoef(..., U, L, Tr_b = Tr)`.
- Produces: `spatial_range_logweights()` and `spatial_range_logweights_grouped()` take `sd_b` in the position of `sigma_b`, recycled with `rep_len(sd_b, p)`.

- [ ] **Step 1: Write the failing tests.**

```r
bbsl_inputs <- function(n = 40L, S = 5L, p = 4L, d = 1L, g = 1L, gt = 1L) {
  set.seed(11)
  list(k = matrix(rbinom(n * S, 1, .4) - .5, n, S), X = matrix(rnorm(n * p), n, p),
       Tr = matrix(rnorm(S * g), S, g), U = matrix(rnorm(n * d), n, d),
       G = matrix(rnorm(g * p), g, p), A = matrix(rnorm(S * gt), S, gt),
       C = matrix(rnorm(gt * p), gt, p), Omega = matrix(rexp(n * S) / 4, n, S),
       empty = list(Gs = matrix(0, g, 0), As = matrix(0, S, gt), Cs = matrix(0, gt, 0),
                    Ks = matrix(0, n, 0), Xs_centers = matrix(0, n, 0)))
}
draw_bbsl <- function(f, x, sd_b, Tr_b = NULL, seed = 5L) {
  seed_both(seed)
  f(x$k, x$X, x$Tr, x$U, x$G, x$A, x$C, sd_b, x$empty$Gs, x$empty$As, x$empty$Cs, 1,
    x$empty$Ks, x$empty$Xs_centers, x$Omega, "binary", 1, 0, Tr_b)
}

test_that("a recycled scalar and a constant vector give the same joint draw", {
  x <- bbsl_inputs()
  for (f in list(occJSDM:::sample_BBsL_cpp, occJSDM:::sample_BBsL_parallel))
    expect_identical(draw_bbsl(f, x, 0.3), draw_bbsl(f, x, rep(0.3, 4)))
})

test_that("each covariate's prior scale is applied to that covariate", {
  x <- bbsl_inputs()
  M <- t(x$Tr %*% x$G + x$A %*% x$C)
  for (f in list(occJSDM:::sample_BBsL_cpp, occJSDM:::sample_BBsL_parallel)) {
    for (tight in 1:4) {
      sd_b <- rep(1, 4); sd_b[tight] <- 1e-6
      B <- draw_bbsl(f, x, sd_b)$B
      expect_equal(B[tight, ], M[tight, ], tolerance = 1e-4)
      expect_gt(max(abs(B[-tight, ] - M[-tight, ])), 1e-2)
    }
  }
  expect_identical(draw_bbsl(occJSDM:::sample_BBsL_cpp, x, c(.2, .5, 1, 2)),
                   draw_bbsl(occJSDM:::sample_BBsL_parallel, x, c(.2, .5, 1, 2)))
})

test_that("the covariate side can use its own trait matrix", {
  x <- bbsl_inputs()
  Tr_b <- cbind(1, x$Tr); G_b <- rbind(rnorm(4), x$G); x2 <- x; x2$G <- G_b
  out <- draw_bbsl(occJSDM:::sample_BBsL_cpp, x2, rep(1e-6, 4), Tr_b = Tr_b)
  expect_equal(out$B, t(Tr_b %*% G_b + x$A %*% x$C), tolerance = 1e-4)
  expect_equal(out$Bt, t(out$B) - Tr_b %*% G_b - x$A %*% x$C)
})

test_that("the range log-weights use one prior scale per covariate", {
  # Dense reference: for each species, 0.5 * h' Q^-1 h - 0.5 * log det Q, with
  # Q = Z' diag(omega) Z + P and h = Z' kappa + P m, summed over species, plus
  # the Gamma(a_l_s, b_l_s) log density of the range at that grid point.
  sim <- simulate_fixture(model = "binary", useSpatField = TRUE, seed = 9L)
  setup <- range_reference_setup(sim)          # helper, see Step 3
  sd_b <- c(0.2, 1.5)
  got <- with(setup, occJSDM:::spatial_range_logweights(X, U, M_B, M_Bs, sd_b, sigma_bs,
    kappa, Omega, Xs_centers, list_SoRSummaries, 1, 1))
  expect_equal(got, range_reference(setup, sd_b), tolerance = 1e-8)
  expect_identical(
    with(setup, occJSDM:::spatial_range_logweights(X, U, M_B, M_Bs, 0.3, sigma_bs,
      kappa, Omega, Xs_centers, list_SoRSummaries, 1, 1)),
    with(setup, occJSDM:::spatial_range_logweights(X, U, M_B, M_Bs, c(0.3, 0.3), sigma_bs,
      kappa, Omega, Xs_centers, list_SoRSummaries, 1, 1)))
})
```

- [ ] **Step 2: Run them.** Expected: FAIL (wrong number of arguments, or a length error on `sd_b`).

- [ ] **Step 3: Add the two range helpers to `tests/testthat/helper-grouped-prior.R`.** `range_reference_setup(sim)` builds the inputs the way `tests/testthat/test-spatial-range.R` builds them for its own dense reference; reuse that file's construction. `range_reference(setup, sd_b)` evaluates the formula in the test's comment with `P = diag(c(1, 1 / sd_b^2, rep(1, d), rep(1 / sigma_bs^2, ps)))` and `m = rbind(0, M_B, 0, M_Bs)`. Add a matching test for `spatial_range_logweights_grouped()` on a fixture with repeated sites, following `tests/testthat/test-spatial-grouped.R`.

- [ ] **Step 4: Implement.** In C++, build the prior diagonal entry by entry with `std::pow(sd_b(i), 2)`, the same operation the scalar code used, so constant vectors reproduce old fits exactly. Stop with `Rcpp::stop("sd_b must have length 1 or one entry per covariate")` otherwise. In R, `rep_len(sd_b, p)` then `1 / sd_b^2` in place of `rep(1 / sigma_b^2, p)`. Thread `Tr_b` only into the covariate-side prior mean and `Bt`. Update every caller in `R/jsdmfun.R` and `R/chains.R` to pass what it passed before.

- [ ] **Step 5: Run this file, then the existing suites that call these functions.**

Run: the test file, then `Rscript -e 'devtools::test(filter = "spatial|intercept|b0|symmetric|rng")'`
Expected: PASS, including Task 1's test.

- [ ] **Step 6: Commit.**

```bash
git add src/jsdm.cpp R/RcppExports.R src/RcppExports.cpp R/jsdmfun.R R/chains.R tests/testthat
git commit -m "Give the coefficient draw and range step a prior scale per covariate"
```

---

### Task 5: The new conditional updates

**Files:**
- Modify: `R/coef-prior.R`
- Modify: `R/jsdmfun.R:1016` (`sample_GC()`), `:1064` (`sample_A()`)
- Test: `tests/testthat/test-grouped-prior.R`

**Interfaces:**
- Consumes: `coef_prior_sds()` and the `coef_prior` list from Task 3.
- Produces: `sample_scale_half_cauchy(ss, n, sigma, scale)`, returning a new standard deviation: `aux <- rinvgamma_cpp(1, 1 / sigma^2 + 1 / scale^2)` then `sqrt(rinvgamma_cpp((n + 1) / 2, ss / 2 + 1 / aux))`, in that order.
- Produces: `sample_coef_scales(R, G, C, coef_prior, sigma_b, sigma_emb, lambda)`, where `R` is the species by covariates matrix of leftovers. Returns `list(sigma_b, sigma_emb, lambda)`. Draw order: `sigma_b`, then `sigma_emb`, then the embedding entries of `lambda` in column order. With no embedding columns, `sigma_emb` and `lambda` are returned unchanged and only `sigma_b` is drawn.
  - `sigma_b`: `ss` is the sum of squared leftovers over measured columns, `n` is `S` times the number of measured columns, `scale` is `sigma_b_scale`.
  - `sigma_emb`: `ss` is the sum over embedding columns of squared leftovers divided by `lambda^2`, `n` is `S` times the number of embedding columns, `scale` is `sigma_emb_scale`.
  - `lambda[k]`: `ss` is `sum(R[, k]^2) / sigma_emb^2 + sum(G[, k]^2) + sum(C[, k]^2) / c_base[k]^2` using the new `sigma_emb`, `n` is `S + nrow(G) + nrow(C)`, `scale` is `axis_scale`.
- Produces: `sample_GC(B, Tr_b, A, sd_b, sd_g, sd_c)`: for covariate `k`, the regression of `B[k, ]` on `cbind(Tr_b, A)` with noise `sd_b[k]` and prior variance `diag(c(sd_g[, k]^2, rep(sd_c[k]^2, gt)))`. Returns `list(G, C, Bt)` as now.
- Produces: `sample_A(B, C, Tr_b, G, sd_b)`. When every entry of `sd_b` is equal it calls `sampleBuniv(t(C), I, 0, Btilde[, s], sd_b[1])` exactly as now. Otherwise it calls `sampleBuniv(t(C) / sd_b, I, 0, Btilde[, s] / sd_b, 1)`.

- [ ] **Step 1: Write the failing tests.**

```r
scale_inputs <- function() {
  set.seed(1); S <- 7L; gb <- 2L; gt <- 2L
  list(S = S, gb = gb, gt = gt, R = matrix(rnorm(S * 5), S, 5),
       G = matrix(rnorm(gb * 5), gb, 5), C = matrix(rnorm(gt * 5), gt, 5),
       cp = list(type = "grouped", is_emb = c(FALSE, FALSE, TRUE, TRUE, TRUE),
                 sigma_b_scale = 1, sigma_emb_scale = 1, axis_scale = 1 / sqrt(3),
                 c_base = c(1, 1, rep(1 / sqrt(gt), 3))),
       sigma_b = 0.7, sigma_emb = 0.4, lambda = c(1, 1, 0.5, 0.2, 0.9))
}
hc <- function(ss, n, sigma, scale) {
  aux <- occJSDM:::rinvgamma_cpp(1, 1 / sigma^2 + 1 / scale^2)
  sqrt(occJSDM:::rinvgamma_cpp((n + 1) / 2, ss / 2 + 1 / aux))
}

test_that("coefficient scales are drawn from their half-Cauchy conditionals", {
  x <- scale_inputs()
  seed_both(42)
  want_b <- hc(sum(x$R[, 1:2]^2), x$S * 2, x$sigma_b, 1)
  want_e <- hc(sum(sweep(x$R[, 3:5], 2, x$lambda[3:5], "/")^2), x$S * 3, x$sigma_emb, 1)
  want_l <- vapply(3:5, function(k) hc(
    sum(x$R[, k]^2) / want_e^2 + sum(x$G[, k]^2) + sum(x$C[, k]^2) / x$cp$c_base[k]^2,
    x$S + x$gb + x$gt, x$lambda[k], x$cp$axis_scale), numeric(1))
  seed_both(42)
  got <- with(x, occJSDM:::sample_coef_scales(R, G, C, cp, sigma_b, sigma_emb, lambda))
  expect_equal(got, list(sigma_b = want_b, sigma_emb = want_e, lambda = c(1, 1, want_l)))
})

test_that("with no embedding columns only the measured scale moves", {
  x <- scale_inputs(); x$cp$is_emb[] <- FALSE; x$cp$c_base[] <- 1
  seed_both(42); want <- hc(sum(x$R^2), x$S * 5, x$sigma_b, 1)
  seed_both(42)
  got <- with(x, occJSDM:::sample_coef_scales(R, G, C, cp, sigma_b, NA_real_, rep(1, 5)))
  expect_equal(got, list(sigma_b = want, sigma_emb = NA_real_, lambda = rep(1, 5)))
})

test_that("an axis scale tracks the size of its column", {
  x <- scale_inputs(); x$R[, 4] <- 0; x$G[, 4] <- 0; x$C[, 4] <- 0
  x$R[, 5] <- 30 * x$R[, 5]
  draws <- replicate(400, with(x, occJSDM:::sample_coef_scales(R, G, C, cp, sigma_b, sigma_emb, lambda))$lambda)
  expect_lt(median(draws[4, ]), 0.2)
  expect_gt(median(draws[5, ]), 5 * median(draws[4, ]))
})

test_that("sample_GC uses a prior and a noise scale per covariate", {
  set.seed(2); S <- 9L; p <- 3L; gb <- 2L; gt <- 2L
  B <- matrix(rnorm(p * S), p, S); Tr_b <- cbind(1, rnorm(S)); A <- matrix(rnorm(S * gt), S, gt)
  sd_b <- c(.3, 1, 2); sd_g <- matrix(c(.5, 1, 2), gb, p, byrow = TRUE); sd_c <- c(.2, .7, 1.5)
  seed_both(7); got <- occJSDM:::sample_GC(B, Tr_b, A, sd_b, sd_g, sd_c)
  seed_both(7)
  for (k in 1:p) {
    want <- occJSDM:::sampleBuniv(cbind(Tr_b, A), diag(c(sd_g[, k]^2, rep(sd_c[k]^2, gt))),
                                  rep(0, gb + gt), B[k, ], sd_b[k])
    expect_equal(c(got$G[, k], got$C[, k]), as.vector(want))
  }
  expect_equal(got$Bt, t(B) - Tr_b %*% got$G - A %*% got$C)
})

test_that("sample_A weights covariates by their noise scale", {
  set.seed(3); S <- 4L; p <- 6L; gt <- 2L
  B <- matrix(rnorm(p * S), p, S); C <- matrix(rnorm(gt * p), gt, p)
  Tr_b <- matrix(1, S, 1); G <- matrix(rnorm(p), 1, p); sd_b <- c(.2, .2, 1, 1, 3, 3)
  Btilde <- B - t(Tr_b %*% G)
  seed_both(8); got <- occJSDM:::sample_A(B, C, Tr_b, G, sd_b)
  seed_both(8)
  for (s in 1:S) expect_equal(got$A[s, ], as.vector(occJSDM:::sampleBuniv(
    t(C) / sd_b, diag(gt), rep(0, gt), Btilde[, s] / sd_b, 1)))
  seed_both(8); same <- occJSDM:::sample_A(B, C, Tr_b, G, rep(.4, p))
  seed_both(8)
  for (s in 1:S) expect_identical(same$A[s, ], as.vector(occJSDM:::sampleBuniv(
    t(C), diag(gt), rep(0, gt), Btilde[, s], .4)))
})
```

- [ ] **Step 2: Run them.** Expected: FAIL, `sample_coef_scales` not found.

- [ ] **Step 3: Implement.** Update the two call sites of `sample_GC()` and `sample_A()` for the spatial coefficients in `update_jSDMcoef()` to pass `rep(sigma_bs, ps)`, a matrix of ones and a vector of ones, so the spatial side is unchanged.

- [ ] **Step 4: Run the tests and Task 1's test.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add R/coef-prior.R R/jsdmfun.R tests/testthat/test-grouped-prior.R
git commit -m "Add the conditional updates for group and axis scales"
```

---

### Task 6: Wire the prior into the sampler, the chains and the returned object

**Files:**
- Modify: `R/runOccJSDM.R` (covariate processing near line 905, priors near 1062, allocation near 1131, storage near 1217, `infos` near 1365)
- Modify: `R/chains.R:86` (`chain_setup_fields()`), `:470-495` (per-chain storage), `:566-611` (initial state), `:776-791`, `:843-853`
- Modify: `R/jsdmfun.R:1243-1495` (`update_jSDMcoef()`)
- Test: `tests/testthat/test-grouped-prior-fits.R`

**Interfaces:**
- Consumes: everything from Tasks 3 to 5.
- Produces, in `runOccJSDM()`: a data frame of design columns (`column`, `variable`, `kind`) built from the result of `process_covariates()`; `coef_prior <- read_coef_prior(listPriors, listParams, S, gt, columns)`; `Tr_b`, which is `cbind(1, Tr)` under `"grouped"` and `Tr` otherwise; `g_b <- ncol(Tr_b)`. `coef_prior`, `Tr_b` and `g_b` are added to `chain_setup_fields()`, and `coef_prior` to `list_priors`, `Tr_b` to `list_data`.
- Produces, in the chain state `list_jSDMparams`: `sigma_emb` and `lambda`. Initial values under `"grouped"`: `sigma_b = 1`, `sigma_emb = sigma_emb_scale`, `lambda` equal to `axis_scale` on embedding columns and 1 elsewhere. Under `"fixed"`: `sigma_emb = NA_real_`, `lambda = rep(1, p)`. `G` is initialised with `g_b` rows.
- Produces, in `update_jSDMcoef()`: at the top, `sds <- coef_prior_sds(...)`. `sds$b` goes to the range function and the joint draw, with `Tr_b`. After the joint draw: under `"fixed"`, `sample_sigmab()` exactly as now; under `"grouped"`, `sample_coef_scales()` on the leftovers `t(B) - Tr_b %*% G - A %*% C`, then `sds` is recomputed. Then `sample_GC(B, Tr_b, A, sds$b, sds$g, sds$c)` and `sample_A(B, C, Tr_b, G, sds$b)`. Every `computePsiCoef()` call passes `Tr_b = Tr_b`.
- Produces, in the returned object: `G_output` has `g_b` rows. When the fit has embedding columns, `jsdm_output` gains `sigmaemb_output` (`niter` by `nchain`) and `lambda_output` (embedding columns by `niter` by `nchain`, with the design column names as its first dimnames). Neither exists otherwise, so other fits keep their shape, as `mu_b0_output` does. `infos` gains `coef_prior` and `covariate_group`, a character vector of `"measured"` and `"embedding"` named by design column.
- Produces: `runOccJSDM()` stops with `"embeddingVars columns must vary across sites: <names>"` when an embedding column is constant.

- [ ] **Step 1: Write the failing tests.**

```r
# tests/testthat/test-grouped-prior-fits.R
# With an embedding, no b_prior is passed: naming embeddingVars must be enough.
fit_grouped <- function(sim, priors = if (embedding) list() else list(b_prior = "grouped"), embedding = TRUE,
                        params = list(), MCMCparams = FIXTURE_MCMC, spatCovariates = NULL) {
  dat <- sim$data_list
  occ <- grouped_occ_covariates(dat)
  if (embedding) params$embeddingVars <- dat$embeddingVars
  suppressMessages(suppressWarnings(runOccJSDM(dat, listParams = c(list(n_factors = 1L), params),
    listPriors = priors, occCovariates = occ, spatCovariates = spatCovariates, MCMCparams = MCMCparams)))
}

test_that("a grouped binary fit returns the new scales with the right shapes", {
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  fit <- fit_grouped(sim)
  jo <- fit$results_output$jsdm_output; mc <- fit$infos$MCMCparams
  expect_identical(unname(fit$infos$covariate_group), rep(c("measured", "embedding"), c(2L, 8L)))
  expect_identical(fit$infos$coef_prior$type, "grouped")
  expect_identical(dim(jo$lambda_output), c(8L, as.integer(mc$niter), as.integer(mc$nchain)))
  expect_identical(dimnames(jo$lambda_output)[[1]], names(fit$infos$covariate_group)[3:10])
  expect_identical(dim(jo$sigmaemb_output), dim(jo$sigmab_output))
  expect_identical(dim(jo$G_output)[1], ncol(sim$data_list$traits) + 1L)
  for (nm in c("lambda_output", "sigmaemb_output", "sigmab_output", "B_output"))
    expect_true(all(is.finite(jo[[nm]])), info = nm)
  expect_true(all(jo$lambda_output > 0) && all(jo$sigmaemb_output > 0))
  expect_gt(length(unique(as.vector(jo$sigmab_output))), 10L)   # no longer pinned
})

test_that("a grouped fit with no embedding columns learns only the measured scale", {
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  fit <- fit_grouped(sim, embedding = FALSE)
  jo <- fit$results_output$jsdm_output
  expect_null(jo$lambda_output); expect_null(jo$sigmaemb_output)
  expect_true(all(fit$infos$covariate_group == "measured"))
})

test_that("a fixed fit has none of the new elements", {
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  fit <- fit_grouped(sim, priors = list(), embedding = FALSE)
  expect_identical(fit$infos$coef_prior, list(type = "fixed"))
  expect_null(fit$results_output$jsdm_output$lambda_output)
  expect_identical(dim(fit$results_output$jsdm_output$G_output)[1], ncol(sim$data_list$traits))
})

test_that("bad embedding inputs stop before sampling", {
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  sim$data_list$info$X_emb.3 <- 1
  expect_error(fit_grouped(sim), "embeddingVars columns must vary across sites: X_emb.3")
  sim2 <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  expect_error(fit_grouped(sim2, params = list(splineVars = "X_emb.1")),
               "embeddingVars must be numeric covariates without splines: X_emb.1")
})
```

- [ ] **Step 2: Run them.** Expected: FAIL, the grouped fit either errors on an unused argument or returns no `lambda_output`.

- [ ] **Step 3: Implement the wiring** as the Interfaces block describes. Allocate and copy the new per-chain arrays beside `mu_b0_output`, using the same present-or-`NULL` pattern.

- [ ] **Step 4: Run the tests, Task 1's test and the chain tests.**

Run: both grouped test files, then `Rscript -e 'devtools::test(filter = "parallel-chains|smoke|api-contracts|rng")'`
Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add R tests/testthat
git commit -m "Wire the grouped coefficient prior into the sampler and its output"
```

---

### Task 7: Accessors, thinning and diagnostics

**Files:**
- Modify: `R/output.R` (`returnTraitsCoeff()` near line 169, new accessors, roxygen)
- Modify: `R/diagnostics.R:366`
- Modify: `NAMESPACE` through `devtools::document()`
- Test: `tests/testthat/test-grouped-prior-fits.R`

**Interfaces:**
- Consumes: the fit structure from Task 6.
- Produces: `returnCovariateScales(fitModel, confidence = 0.95)`, exported. A data frame with columns `scale`, `covariate`, `mean`, `lower`, `median`, `upper`: one row with `scale = "measured group"`, and for fits with embedding columns one row with `scale = "embedding leftover"` and one row per embedding column with `scale = "axis"`. `covariate` is `NA` for the group rows. Stops with `"returnCovariateScales() needs a fit made with the grouped prior"` on a fixed fit.
- Produces: `returnCommunityResponse(fitModel, confidence = 0.95)`, exported. One row per design column with `covariate`, `group`, `mean`, `lower`, `median`, `upper`, from the first row of `G_output`. Same error on a fixed fit.
- Produces: `returnEmbeddingBlends(fitModel, n_blends = NULL)`, exported. For each draw, form the species by embedding-columns matrix `A %*% C[, embedding columns]`, take its singular value decomposition, keep the first `n_blends` right singular vectors (default: `gt`), and flip each so that its entry of largest absolute value is positive. Returns `list(weights, site_scores, singular_values)`: `weights` is embedding columns by blends, the posterior mean of the flipped vectors; `site_scores` is sites by blends, the embedding columns of `fitModel$X_psi` times `weights`; `singular_values` is a data frame of posterior summaries. The documentation states that blends are defined only up to rotation and must not be interpreted one by one. Stops with `"returnEmbeddingBlends() needs a fit with embeddingVars and at least one latent trait"` otherwise.
- Produces: `returnEmbeddingShare(fitModel, confidence = 0.95)`, exported. Per species, the share of the variance of the covariates' contribution to the linear predictor that comes from the embedding columns: for each draw, `v_e / (v_e + v_m)` with `v_e` the variance across sites of the embedding columns times their coefficients and `v_m` the same for the measured columns. Computed from `B_output` and `fitModel$X_psi`; nothing new is stored.
- Produces: under `"grouped"`, `returnTraitsCoeff()` and `plotTraitsCoefficients()` drop the community-average row and return the trait rows as before. `computeDiagnostics()` skips `A_output` and `C_output` when `lambda_output` is present.
- Two deliberate departures from the spec's wording, both within what it leaves to the builder. The split of the covariates' share between measured and embedding is its own accessor, `returnEmbeddingShare()`, and `returnVariationPartitioning()` is left as it is, because its three-way partition and ternary plot cannot hold a fourth component. `computeDiagnostics()` prints R-hat on the scales as stored; the log-scale R-hat the spec asks for is computed by the validation study, not by the package.

- [ ] **Step 1: Write the failing tests.**

```r
grouped_fit <- function() cached_fixture("grouped_binary", function()
  fit_grouped(simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)))

test_that("the scale and community accessors return one row per quantity", {
  fit <- grouped_fit()
  sc <- returnCovariateScales(fit)
  expect_identical(names(sc), c("scale", "covariate", "mean", "lower", "median", "upper"))
  expect_identical(sc$scale, c("measured group", "embedding leftover", rep("axis", 8)))
  expect_identical(sc$covariate[3:10], names(fit$infos$covariate_group)[3:10])
  expect_true(all(sc$lower <= sc$median & sc$median <= sc$upper))
  cr <- returnCommunityResponse(fit)
  expect_identical(cr$covariate, names(fit$infos$covariate_group))
  expect_identical(cr$group, unname(fit$infos$covariate_group))
})

test_that("blends are sign-fixed and reproduce the low-rank part", {
  fit <- grouped_fit(); bl <- returnEmbeddingBlends(fit, n_blends = 2L)
  expect_identical(dim(bl$weights), c(8L, 2L))
  expect_identical(nrow(bl$site_scores), nrow(fit$X_psi))
  expect_true(all(apply(bl$weights, 2, function(v) v[which.max(abs(v))]) > 0))
  expect_error(returnEmbeddingBlends(fit, n_blends = 99L), "n_blends cannot exceed the number of latent traits")
})

test_that("the embedding share lies between 0 and 1", {
  sh <- returnEmbeddingShare(grouped_fit())
  expect_identical(nrow(sh), 12L)
  expect_true(all(sh$mean >= 0 & sh$mean <= 1))
})

test_that("existing downstream functions work on a grouped fit", {
  fit <- grouped_fit(); sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  expect_identical(nrow(returnTraitsCoeff(fit)) > 0, TRUE)
  thin <- thinOutput(fit, thin = 2)
  expect_identical(dim(thin$results_output$jsdm_output$lambda_output)[2],
                   length(seq(1, fit$infos$MCMCparams$niter, by = 2)))
  expect_no_error(suppressWarnings(utils::capture.output(computeDiagnostics(fit$results_output))))
  newX <- unique(sim$data_list$info[, grouped_occ_covariates(sim$data_list)])[1:5, ]
  pred <- predictNewSites(fit, X_psi = newX, useSpatial = FALSE, useBiotic = FALSE, verbose = FALSE)
  expect_false(anyNA(unlist(pred)))
})

test_that("the accessors refuse fits they do not apply to", {
  sim <- simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)
  fixed <- fit_grouped(sim, priors = list(), embedding = FALSE)
  expect_error(returnCovariateScales(fixed), "needs a fit made with the grouped prior")
  expect_error(returnCommunityResponse(fixed), "needs a fit made with the grouped prior")
  expect_error(returnEmbeddingBlends(fixed), "needs a fit with embeddingVars and at least one latent trait")
})
```

- [ ] **Step 2: Run them.** Expected: FAIL, `returnCovariateScales` not found.

- [ ] **Step 3: Implement**, then run `devtools::document()`.

- [ ] **Step 4: Run both grouped test files and `devtools::test(filter = "diagnostics|posterior-draws|covariate-response|api-contracts")`.** Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add R NAMESPACE man tests/testthat
git commit -m "Add accessors for the grouped prior and keep downstream functions working"
```

---

### Task 8: Prior check, recovery and the remaining model types

**Files:**
- Create: `dev/simstudy/grouped-prior-checks-20261111/prior_check.R`
- Create: `dev/simstudy/grouped-prior-checks-20261111/README.md`
- Test: `tests/testthat/test-grouped-prior-fits.R`

**Interfaces:**
- Consumes: the finished feature.
- Produces: `prior_check.R`, which draws from the prior only and prints, for 64 embedding columns with `gt = 8` and for 128 with `gt = 8`, the median and the 5% and 95% points of the standard deviation across sites of the embedding's contribution to the linear predictor.

- [ ] **Step 1: Write the prior check.** For each of 2,000 draws: `lambda` from a half-Cauchy with scale `axis_scale`, `sigma_emb` from a half-Cauchy with scale 1, a community row and `gt` rows of `C` with standard deviations `lambda` and `lambda / sqrt(gt)`, `A` standard Normal for one species, a leftover with standard deviation `sigma_emb * lambda`; form that species' coefficients, multiply by a standardised 250 by K matrix built with the factor model of Task 2, and record the standard deviation across sites. Run it with `Rscript`.

- [ ] **Step 2: Judge it against the spec.** The spec asks for a prior median between 1 and 2. If the median is inside that range, record the numbers in the README and go on. If it is outside, do not change the defaults: record the numbers and report to Doug with the scale that would bring the median to 1.5.

- [ ] **Step 3: Write the fit tests.**

```r
test_that("axes with no effect are shrunk more than axes with an effect", {
  skip_on_cran()
  sim <- simulate_grouped_fixture(n = 150L, S = 30L, p_meas = 2L, p_emb = 10L, gt = 2L, seed = 4L,
    embedding = list(n_blends = 2L, n_active = 3L, blend_sd = 1.2, sigma_emb = 0.1))
  fit <- fit_grouped(sim, MCMCparams = list(nchain = 2, nburn = 1500, niter = 1500, nthin = 1))
  lam <- apply(fit$results_output$jsdm_output$lambda_output, 1, median)
  active <- rowSums(abs(sim$true_params$jsdmParams_true$B[3:12, ])) > 0
  expect_gt(min(lam[active]), max(lam[!active]))
})

test_that("grouped fits run for every supported model type and layout", {
  skip_on_cran()
  for (model in c("occupancy", "two_stage")) {
    sim <- simulate_grouped_fixture(model = model, p_emb = 6L, S = 8L)
    jo <- fit_grouped(sim)$results_output$jsdm_output
    expect_true(all(is.finite(jo$lambda_output)), info = model)
  }
  wide <- simulate_grouped_fixture(n = 30L, p_emb = 40L, S = 8L)         # more columns than sites
  expect_true(all(is.finite(fit_grouped(wide)$results_output$jsdm_output$B_output)))
  noblend <- simulate_grouped_fixture(p_emb = 6L, S = 8L, gt = 0L)
  fit0 <- fit_grouped(noblend, params = list(n_lattrait = 0L))
  expect_true(all(is.finite(fit0$results_output$jsdm_output$lambda_output)))
  par <- fit_grouped(simulate_grouped_fixture(p_emb = 6L, S = 8L),
                     MCMCparams = c(FIXTURE_MCMC, list(ncores = 2)))
  ser <- fit_grouped(simulate_grouped_fixture(p_emb = 6L, S = 8L))
  expect_identical(dim(par$results_output$jsdm_output$lambda_output),
                   dim(ser$results_output$jsdm_output$lambda_output))
})

test_that("a grouped fit runs with a spatial field", {
  skip_on_cran()
  sim <- simulate_fixture(model = "binary", useSpatField = TRUE, seed = 6L)
  fit <- suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
    listParams = list(n_factors = 1L, n_supportpoints = FIXTURE_KNOTS, embeddingVars = fixture_occ_covariates()[2]),
    listPriors = list(b_prior = "grouped"), occCovariates = fixture_occ_covariates(),
    spatCovariates = fixture_spat_covariates(), MCMCparams = FIXTURE_MCMC)))
  expect_true(all(is.finite(fit$results_output$jsdm_output$lambda_output)))
})
```

- [ ] **Step 4: Run them.** Expected: PASS. If the shrinkage test fails, do not loosen it: rerun with three other seeds, record all four outcomes in the README, and report to Doug.

- [ ] **Step 5: Commit.**

```bash
git add dev/simstudy/grouped-prior-checks-20261111 tests/testthat/test-grouped-prior-fits.R
git commit -m "Check the grouped prior's defaults and its fits across model types"
```

---

### Task 9: Documentation and the full suite

**Files:**
- Modify: `R/coef-prior.R` (a roxygen block for the new help page), `R/runOccJSDM.R` (two pointer sentences only), `R/output.R` (a `@seealso` on each accessor), `man/`, `_pkgdown.yml`
- Modify: `dev/superpowers/specs/2026-11-11-grouped-covariate-prior-design.md` (status line only)

**Interfaces:**
- Consumes: everything above.
- Produces: the help page `grouped-prior`, with the aliases `groupedPrior`, `embeddingVars` and `b_prior`, so that `?groupedPrior` and `?embeddingVars` both open it. Doug decided on 12 November 2026 that the feature has its own help page, so that `?runOccJSDM` stays as standard users know it.

- [ ] **Step 1: Write the help page.** Add a roxygen block ending in `NULL` to `R/coef-prior.R`, with `@name grouped-prior` and the three aliases, following the block for `variance-partitioning-deprecated` at `R/output.R` line 2198 but without its `@keywords internal`, which would hide the page from the index. Title: "Grouped prior for occupancy coefficients (experimental)". Sections, in this order:
  - *What it is for:* embedding columns as occupancy covariates beside measured ones.
  - *How to switch it on:* naming `listParams$embeddingVars` is enough; `listPriors$b_prior = "grouped"` is needed only without an embedding. Include one example call with measured covariates and an embedding.
  - *What else changes:* the prior on the measured covariates becomes the learned measured-group scale, and every covariate gets a community-average response.
  - *Settings:* `sigma_b_scale`, `sigma_emb_scale` and `axis_scale`, with their defaults, and `n_lattrait` as the number of blends.
  - *What the fit holds:* `lambda_output`, `sigmaemb_output`, `infos$coef_prior`, `infos$covariate_group`, and the leading community-average row of `G_output`; then the four accessors.
  - *Reading a management covariate:* two or three sentences saying that its coefficient, fitted alongside the embedding, is a direct effect (the effect not passing through habitat the embedding describes), not the total effect, and that the total effect comes from a second fit without the embedding and without other covariates that management can change. The wording follows the spec's section *What the management coefficient estimates*.
  - *Status:* experimental, and not yet through its validation study.

- [ ] **Step 2: Leave `?runOccJSDM` almost untouched.** Add one sentence to `@details`, beside the `b0_prior` text: embeddings as occupancy covariates use an experimental prior that is off unless `embeddingVars` or `b_prior` is set, documented in the `grouped-prior` page. Add one sentence to `@return`: fits made with that prior hold further elements, listed on the same page. Describe none of the settings or outputs here. Add `@seealso` pointing to the page on each of the four accessors.

- [ ] **Step 3: List the new topics on the site.** In `_pkgdown.yml`, add a reference section titled "Grouped prior for embeddings (experimental)" after *Covariate effects*, containing `grouped-prior`, `returnCovariateScales`, `returnCommunityResponse`, `returnEmbeddingBlends` and `returnEmbeddingShare`. The site build (`.github/workflows/pkgdown.yaml`) fails on an exported topic that is missing from the index.

- [ ] **Step 4: Regenerate and run everything.**

Run: `Rscript -e 'devtools::document()'` then `Rscript -e 'devtools::test()' 2>&1 | tail -25`
Expected: 0 failures. Record the counts of passes, skips and warnings.

Run: `test -f man/grouped-prior.Rd && grep -c 'grouped-prior' man/runOccJSDM.Rd`
Expected: the file exists and the count is 2.

Run, if `pkgdown` is installed: `Rscript -e 'pkgdown::check_pkgdown()'`
Expected: no problems found. If it is not installed, say so in the report; do not install it.

- [ ] **Step 5: Check the fixed prior once more, unskipped.** Run Task 1's test file alone and confirm 4 passes and 0 skips.

- [ ] **Step 6: Update the spec's status line** to say part 1's package code is implemented and unit-tested as of the date, and that the extract, the validation study and part 2 are not started.

- [ ] **Step 7: Commit.**

```bash
git add R man _pkgdown.yml dev/superpowers/specs/2026-11-11-grouped-covariate-prior-design.md
git commit -m "Document the grouped coefficient prior on its own help page"
```

---

## After this plan

In order, each with its own plan: the real embedding extract; the validation study on truths 1 to 4, which gates everything else; then part 2, taxonomy. Nothing here reaches `main` until Doug decides how, once part 1 is validated.

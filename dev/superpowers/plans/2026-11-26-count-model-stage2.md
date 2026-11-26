# Count Model Stage 2 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Fit the plain count JSDM through `runOccJSDM()`: negative-binomial counts of individuals with a log link, an effort offset, optional replicate samples per site and one size per species, asked for with `countModel = list(occupancy = FALSE)`. Fits of the other four models stay identical draw for draw.

**Architecture:** The count likelihood enters the sampler only through the two inputs of the existing Gaussian block. For each iteration the count branch of `update_jSDMcoef()` draws Polya-Gamma weights with the stage 1 sampler, sums them over a site's samples into `Omega` and `k`, and hands them to the unchanged coefficient, range and factor updates under the name the C++ already uses for that form ("binary"). The size is then updated by the stage 1 Metropolis step. Everything else is plumbing: reading `countModel`, grouping sample rows to sites, storing the size, and teaching the post-fit functions about a log link.

**Tech Stack:** R 4.6.1 (`/usr/local/bin/Rscript`), Rcpp and RcppArmadillo, testthat 3e, devtools, roxygen2 8.1.0, MASS for one comparison test.

**Spec:** `dev/superpowers/specs/2026-11-12-count-model-design.md`: *The model*, *How it is fitted*, *Proposed interface*, *Staged plan* (stage 2), and questions 4, 8, 10, 11, 12 and 20. Stage 1's plan and report: `dev/superpowers/plans/2026-11-25-count-model-stage1.md` and `dev/count-model/stage1-report.md`.

## Global Constraints

- With `countModel = NULL`, every fit of the binary, continuous, occupancy and two-stage models is identical draw for draw to the same fit before this stage, and its stored output has the same elements. Task 1 saves reference fits before any edit and later tasks check against them.
- A count model is fitted only when the user supplies `countModel`. It is never chosen from the data. Data of 0 and 1 without `countModel` are still fitted as "binary".
- `countModel` is a list. `occupancy` has no default and must be `TRUE` or `FALSE`. `family` defaults to "negbin", the only family. `offset` names an effort column of `data$info` on the natural scale; without it the effort is 1.
- `occupancy = TRUE` (the shared-predictor model) stops with a clear error in this stage. It is stage 3.
- One size per species. The prior on the size is the stage 1 one, a half-normal on its inverse square root with standard deviation `listPriors$size_prior_sd`, default 1. Alex has not chosen between this, a community prior and a tighter half-normal (question 5); see *Decisions made in this plan*.
- The intercepts of a count fit get the community prior by default (question 8).
- No Poisson family. No embedding block. No taxonomic structure. `collCovariates` are not used by a count model.
- Count code in `update_jSDMcoef()` only forms `k` and `Omega`. Do not edit `sample_BBsL_cpp()`, `sample_U_cpp()`, `spatial_range_logweights()` or `samplepg_fast()`. Keep every edit to `R/jsdmfun.R`, `R/chains.R`, `R/runOccJSDM.R` and `src/jsdm.cpp` small and additive, because the grouped-prior pull requests (41, 44, 45, 46) rewrite the same functions and will be rebased onto this.
- Never run `R CMD build`, `R CMD check`, `devtools::check()` or `devtools::build()` on a checkout. Test with `devtools::load_all()` and `devtools::test()`.
- After changing any `// [[Rcpp::export]]` signature run `Rscript -e 'Rcpp::compileAttributes()'`. After changing roxygen run `Rscript -e 'roxygen2::roxygenise()'`. Never edit the generated files by hand.
- Do not edit an R script while an `Rscript` process is running it. Runs longer than two minutes go in the background; the stage 2 measurements run detached with `nohup`.
- No em-dashes anywhere, including commit messages and code comments. Markdown follows the project `CLAUDE.md`: one line per paragraph, no pipe tables.
- Commit messages end with `Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>`. Do not push or open a pull request until Doug says so.
- **Branch.** If pull request 47 has merged, cut `claude/count-model-stage2` from `origin/main`. If not, cut it from `origin/claude/count-model-stage1`, merge `origin/main` into it at once so the design and plans are present, and target the stage 1 branch with the eventual pull request. Work in its own worktree and never switch branches in the main checkout.

## Decisions made in this plan

Each of these goes beyond the design and is for Doug and Alex to confirm.

- **The size prior is the stage 1 independent half-normal,** with the standard deviation exposed as `listPriors$size_prior_sd`. Stage 1 showed that it sets a rare species' size at about 2. If Alex chooses the community prior, it is one further task after this plan: a normal distribution on the log sizes with an estimated mean and spread. Nothing here blocks it.
- **The intercept prior's spread stays at the package default,** `sigma_b0 = 1`. A count fit's intercepts are log expected counts, and species' abundances can differ by far more than that. Task 11 measures what the default costs when the true spread is 2.
- **The stored parameter is called `size`:** `size_output` in the fit, "size" in `returnPosteriorDraws()`.
- **A new exported function, `returnZeroCheck()`,** gives the posterior predictive check on zeros. Doug to confirm the name.
- **The independent comparison is with `MASS::glm.nb()`,** on a model with intercepts only. The design asked for Hmsc or gllvm. Neither is installed under R 4.6, and Hmsc has no negative binomial. A comparison with gllvm on a model with latent factors needs Doug's go-ahead to install it and is not in this plan.
- **The complete-data WAIC of a count fit has one term per sample row and species,** not per site and species.
- **`predictNewSites()` returns the expected count per unit of effort** for a count fit.
- **The wording of the existing "Counts model not supported yet" error changes,** and Lesson 2 quotes it, so one line of the lesson changes too.
- **Missing counts are not allowed,** as for every model except the two-stage one. Delete the row.
- **`collCovariates` with a count model is an error,** not silently ignored as it is for a binary fit.

## Out of scope: a second defect found while planning

The stored WAIC of a continuous fit is computed with a residual standard deviation of 1 for every species. `run_one_chain()` passes the local `tau` to `computeModelLoglikJSDM_cpp()` (`R/chains.R` line 797), and that local is set once (line 589) and never updated; the sampled value is in `list_jSDMparams$tau`. No draw is affected, only `results_output$WAIC` of continuous fits. Fixing it changes that stored number, so it is not part of this stage. The count code must not repeat it: read the size from `list_jSDMparams`.

## Review Focus

- **Rows that are not grouped by site:** count data with replicate rows in any order, or with no `Site` column at all, must reach the sampler with each row tied to its own site's covariates. Pinned in Task 4.
- **Effort that is missing or unusable:** no `offset` means an effort of 1; an effort of zero, a negative effort or `NA` is an error naming the column. Pinned in Task 2.
- **A species never counted, and a species counted in the hundreds:** the first has no information about its size, the second sends the sampler down its normal route. Both must fit with finite draws. Pinned in Task 6.
- **Data of 0 and 1:** without `countModel` they are still a binary fit; with it they are a count fit. Pinned in Task 4.
- **One species:** the community prior needs two, so a single-species count fit must fall back to the fixed prior and still run. Pinned in Task 6.

---

### Task 1: Reference fits and the baseline

**Files:**
- Create: `dev/count-model/stage2-unchanged.R`

**Interfaces:**
- Produces: `Rscript dev/count-model/stage2-unchanged.R save <dir>` and `... check <dir>`. `<dir>` is outside the repository. `check` prints one line per model and exits with status 1 if any differs.

- [ ] **Step 1: Write the script.** It runs `devtools::load_all(quiet = TRUE)`, sources `tests/testthat/helper-fixtures.R`, and for each of "binary", "continuous", "occupancy" and "two_stage" does `sim <- simulate_fixture(model, seed = 11)`, `set.seed(77)`, `fit <- fit_fixture(sim)`. In `save` mode it writes `list(data = sim$data_list, results = fit$results_output)` to `<dir>/<model>.rds`. In `check` mode it recomputes and compares both elements with `identical()`.

- [ ] **Step 2: Save the references before any other edit.**

Run: `Rscript dev/count-model/stage2-unchanged.R save /tmp/occjsdm-stage2-reference`
Expected: four `.rds` files.

Run: `Rscript dev/count-model/stage2-unchanged.R check /tmp/occjsdm-stage2-reference`
Expected: four lines ending "identical", exit status 0. If this fails on an untouched tree the fits are not reproducible on this machine: stop and tell Doug.

- [ ] **Step 3: Record the baseline suite** in the background and keep the counts for the pull request.

Run: `Rscript -e 'devtools::test()' 2>&1 | tail -n 8`

- [ ] **Step 4: Commit**

```bash
git add dev/count-model/stage2-unchanged.R
git commit -m "Add a before-and-after check that existing fits are unchanged"
```

### Task 2: Read and validate `countModel`

**Files:**
- Create: `R/count_model.R`
- Create: `tests/testthat/test-count-model.R`

**Interfaces:**
- Produces, internal:
  - `read_count_model(countModel, data_info)`: `NULL` when `countModel` is `NULL`; otherwise `list(family = "negbin", occupancy = FALSE, offset = <name or NULL>, log_effort = <one value per row of data_info>)`.
  - `check_count_response(OTU)`: stops unless every value is a whole number of at least 0 with none missing.
  - `read_size_prior(listPriors, counts)`: the prior standard deviation, default 1.

- [ ] **Step 1: Write the failing tests**

```r
# Stage 2 of the count model: the plain count JSDM fitted through
# runOccJSDM(). Design: dev/superpowers/specs/2026-11-12-count-model-design.md.

info <- data.frame(Site = c("a", "a", "b"), days = c(1, 2, 4), bad = c(1, 0, 2))

test_that("countModel is read with family defaulting and occupancy required", {
  expect_null(occJSDM:::read_count_model(NULL, info))
  got <- occJSDM:::read_count_model(list(occupancy = FALSE, offset = "days"), info)
  expect_identical(got$family, "negbin")
  expect_false(got$occupancy)
  expect_equal(got$log_effort, log(c(1, 2, 4)))
  expect_equal(occJSDM:::read_count_model(list(occupancy = FALSE), info)$log_effort, c(0, 0, 0))
})

test_that("countModel fails with a message that says what to change", {
  f <- function(x) occJSDM:::read_count_model(x, info)
  expect_error(f("negbin"), "must be NULL or a list")
  expect_error(f(list(family = "negbin")), "occupancy must be TRUE or FALSE")
  expect_error(f(list(occupancy = NA)), "occupancy must be TRUE or FALSE")
  expect_error(f(list(occupancy = TRUE)), "not available yet")
  expect_error(f(list(occupancy = FALSE, family = "poisson")), "family must be \"negbin\"")
  expect_error(f(list(occupancy = FALSE, effort = "days")), "unknown entries: effort")
  expect_error(f(list(occupancy = FALSE, offset = "hours")), "must name one column of data\\$info")
  expect_error(f(list(occupancy = FALSE, offset = "bad")), "data\\$info\\$bad must be finite and positive")
})

test_that("the count response must be whole numbers of at least zero", {
  expect_silent(occJSDM:::check_count_response(matrix(c(0, 3, 250, 1), 2)))
  for (bad in list(c(0, 1.5), c(0, -1), c(0, NA)))
    expect_error(occJSDM:::check_count_response(matrix(bad, 2)), "whole numbers of at least 0")
})

test_that("the size prior is a count-only setting", {
  expect_identical(occJSDM:::read_size_prior(list(), TRUE), 1)
  expect_identical(occJSDM:::read_size_prior(list(size_prior_sd = 0.5), TRUE), 0.5)
  expect_error(occJSDM:::read_size_prior(list(size_prior_sd = 0), TRUE), "finite positive number")
  expect_error(occJSDM:::read_size_prior(list(size_prior_sd = 1), FALSE), "requires countModel")
  expect_null(occJSDM:::read_size_prior(list(), FALSE))
})
```

- [ ] **Step 2: Run and see them fail.**

Run: `Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-count-model.R")'`

- [ ] **Step 3: Write `R/count_model.R`.** Messages, in full, each raised with `call. = FALSE`:

  - "countModel must be NULL or a list; see ?countModel"
  - "countModel has unknown entries: " followed by the names, comma separated. The known entries are `family`, `occupancy` and `offset`.
  - "countModel$family must be \"negbin\""
  - "countModel$occupancy must be TRUE or FALSE: TRUE fits the shared-predictor model, in which the covariates act on occupancy; FALSE fits the plain count JSDM, in which they act on the expected count"
  - "countModel$occupancy = TRUE (the shared-predictor model) is not available yet; use occupancy = FALSE"
  - "countModel$offset must name one column of data$info"
  - "The effort in data$info$" followed by the column name and " must be finite and positive"
  - "With countModel set, data$OTU must hold whole numbers of at least 0, with no missing values"
  - "size_prior_sd must be a finite positive number" and "size_prior_sd requires countModel"

- [ ] **Step 4: Run the test file.** Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add R/count_model.R tests/testthat/test-count-model.R
git commit -m "Read and validate the countModel argument"
```

### Task 3: Simulate counts

**Files:**
- Modify: `R/simulateData.R` (`simulateOccJSDMData()`: the model check near line 220, the model mapping near line 260, the `info` block near line 379, the returned truths near line 406, and the roxygen for `model` and `rnb`)
- Modify: `tests/testthat/helper-fixtures.R` (`simulate_fixture()`)
- Test: `tests/testthat/test-count-model.R`

**Interfaces:**
- Produces: `simulateOccJSDMData(list_datasettings, list_params, list_jsdmParams, model = "counts")`.
  - Reads `list_jsdmParams$rnb`, one positive size per species, and `list_datasettings$M`, samples per site. An optional `list_datasettings$effort` has one positive value per sample row, default 1.
  - `data_list$info` has columns `Site` (1 to n, repeated `M[i]` times), `Effort`, then the `X_psi` and `Xs` columns of the row's site. `data_list$OTU` has one row per sample.
  - `true_params` gains `size_true` and `effort`.
- Produces: `simulate_fixture(model = "counts", ...)`, which uses `M` samples per site as the replicated designs do, sets `rnb = rep(2, S)`, and draws an effort from `stats::runif(sum(M), 0.5, 3)` after its other random draws.

- [ ] **Step 1: Write the failing test**

```r
test_that("the simulator draws counts per sample row with an effort column", {
  sim <- simulate_fixture(model = "counts", n = 30, S = 3, M = 2, seed = 4)
  info <- sim$data_list$info
  y <- sim$data_list$OTU
  expect_identical(nrow(y), 60L)
  expect_identical(info$Site, rep(1:30, each = 2))
  expect_true(all(info$Effort > 0))
  expect_true(all(y >= 0 & y == round(y)))
  expect_equal(sim$true_params$size_true, rep(2, 3))
  expect_identical(simulate_fixture(model = "counts", n = 30, S = 3, M = 2, seed = 4), sim)
  eta <- sim$true_params$jsdmParams_true$eta
  expect_gt(stats::cor(as.vector(log1p(y)), as.vector(eta[info$Site, ] + log(info$Effort))), 0.5)
})

test_that("a count simulation needs one positive size per species", {
  settings <- list(n = 5, S = 2, g = 1, M = rep(1, 5), P = 1, K = 1, ncov_psi = 1, ncov_theta = 1)
  jsdm <- list(gt = 1, d = 0, ds = 0, sigma_b = 1, sigma_bs = 1, sigma_ts = 1, sigma_h = 1,
               sigma_s = 1, l_s = 1, tau = rep(1, 2), rnb = c(1, -1), useSpatField = FALSE)
  expect_error(simulateOccJSDMData(settings, list(), jsdm, "counts"),
               "rnb must give one positive size per species")
})
```

- [ ] **Step 2: Run and see it fail** with "The model used is not supported".

- [ ] **Step 3: Implement.** Add "counts" to the allowed models and map it to `simulateData(..., model = "counts")`. Ignore the `z` that `simulateData()` returns for counts and draw the counts afterwards, as the last random draws of the function: for sample row `m` at site `i` and species `j`, a negative binomial with `size = rnb[j]` and mean `effort[m] * exp(eta[i, j])`. Drawing last keeps every other model's simulated data unchanged. Rewrite the roxygen for `rnb` (it now sets the size for `model = "counts"` and is ignored otherwise) and for `model`.

- [ ] **Step 4: Run the test file and the unchanged check.**

Run: `Rscript dev/count-model/stage2-unchanged.R check /tmp/occjsdm-stage2-reference`
Expected: four lines ending "identical".

- [ ] **Step 5: Commit**

```bash
git add R/simulateData.R man tests/testthat/helper-fixtures.R tests/testthat/test-count-model.R
git commit -m "Simulate count data with effort and replicate samples"
```

### Task 4: Take count data into `runOccJSDM()`

**Files:**
- Modify: `R/runOccJSDM.R` (signature near line 636; `inferDataModel()` lines 153 to 254; the model tests near lines 720, 727, 768, 828, 927, 954 and 1475; `list_data` near line 1111)
- Modify: `vignettes/occJSDM-lesson-2.Rmd` and `vignettes/occJSDM-lesson-2.md` (the one quoted error message)
- Test: `tests/testthat/test-count-model.R`

**Interfaces:**
- Consumes: `read_count_model()`, `check_count_response()`, `read_size_prior()`.
- Produces: `runOccJSDM(..., listPriors = list(), countModel = NULL)`, with the new argument last. `inferDataModel(data, countModel = NULL)` returns "counts" whenever `countModel` is given.
- Produces, for the chain: `list_data$count`, a list with `y` (sample rows by species, rows sorted by site), `site_index` (an integer vector, the site of each row, 1 to `n`), `log_effort`, and `size_prior_sd`. For a count fit `z` is the `n` by `S` matrix of site totals and is used only for its dimensions.
- Produces: an internal test hook. With `options(occJSDM.return_count_setup = TRUE)`, a count fit returns `list(count = list_data$count, X_psi = X_psi, siteNames = siteNames, model = model, intercept_prior = intercept_prior)` just before the chains start. The hook stays in the package after this stage, undocumented, so that the data handling can be tested without a fit.
- Until Task 6, a count fit without the hook stops just before the chains start with `stop("count sampler not wired yet")`.

- [ ] **Step 1: Write the failing tests.** `count_setup()` uses the test hook to get the prepared data back without fitting.

```r
count_setup <- function(data, ...) {
  old <- options(occJSDM.return_count_setup = TRUE)
  on.exit(options(old))
  suppressMessages(runOccJSDM(data, ...))
}

test_that("count rows are tied to their own site whatever their order", {
  info <- data.frame(Site = c("b", "a", "b", "c", "a"), elev = c(20, 10, 20, 30, 10),
                     days = c(1, 2, 3, 4, 5))
  y <- matrix(c(5, 0, 7, 2, 1, 0, 0, 3, 9, 4), 5, 2, dimnames = list(NULL, c("s1", "s2")))
  got <- count_setup(list(info = info, OTU = y), occCovariates = "elev",
                     countModel = list(occupancy = FALSE, offset = "days"))
  expect_identical(got$model, "counts")
  expect_identical(as.character(got$siteNames), c("a", "b", "c"))
  expect_identical(got$count$site_index, c(1L, 1L, 2L, 2L, 3L))
  expect_equal(unname(got$count$y[, "s1"]), c(0, 1, 5, 7, 2))
  expect_equal(got$count$log_effort, log(c(2, 5, 1, 3, 4)))
  expect_identical(nrow(got$X_psi), 3L)
  expect_identical(got$intercept_prior$type, "community")
})

test_that("count data without a Site column are one sample per site", {
  y <- matrix(c(3, 0, 12, 1), 4, 1, dimnames = list(NULL, "s1"))
  got <- count_setup(list(info = data.frame(elev = 1:4), OTU = y), occCovariates = "elev",
                     countModel = list(occupancy = FALSE))
  expect_identical(got$count$site_index, 1:4)
  expect_equal(got$count$log_effort, rep(0, 4))
  expect_null(got$intercept_prior$type)
})

test_that("a count model is fitted only when asked for", {
  info <- data.frame(Site = 1:6, elev = 1:6)
  counts <- matrix(c(3, 0, 12, 1, 0, 7), 6, 1, dimnames = list(NULL, "s1"))
  expect_error(suppressMessages(runOccJSDM(list(info = info, OTU = counts))),
               "fits counts only when asked")
  zero_one <- matrix(c(1, 0, 1, 1, 0, 0), 6, 1, dimnames = list(NULL, "s1"))
  expect_identical(suppressMessages(occJSDM:::inferDataModel(list(info = info, OTU = zero_one))), "binary")
  expect_identical(count_setup(list(info = info, OTU = zero_one),
                               countModel = list(occupancy = FALSE))$model, "counts")
})

test_that("count fits refuse what they cannot use", {
  info <- data.frame(Site = 1:4, elev = 1:4)
  y <- matrix(c(3, 0, 12, 1), 4, 1, dimnames = list(NULL, "s1"))
  expect_error(count_setup(list(info = info, OTU = y), collCovariates = "elev",
                           countModel = list(occupancy = FALSE)), "collCovariates are not used")
  expect_error(count_setup(list(info = info, OTU = y + 0.5), countModel = list(occupancy = FALSE)),
               "whole numbers of at least 0")
})
```

- [ ] **Step 2: Run and see them fail** with "unused argument (countModel".

- [ ] **Step 3: Implement.**

  - `inferDataModel()`: with `countModel` given, run `check_count_response()` and return "counts" without looking at the sampling design. With `countModel = NULL` it behaves as now, except that the stop at line 250 becomes "Count data detected, but occJSDM fits counts only when asked. Set countModel, for example countModel = list(occupancy = FALSE); see ?countModel".
  - Count data follow the occupancy path for grouping rows to sites: add "counts" where the code adds a `Sample` column (line 720), sorts by site (line 727), computes `M`, `n`, `N` and `siteNames` (line 828), calls `createDataIdx()` with `twostage = FALSE` for `idx_z_w` (line 927), and groups `X_psi` by site (line 954). Do not add it where the threshold is applied (line 750) or where `X_theta` is built (line 977).
  - `read_count_model()` is called on the sorted `data_info`, so that `log_effort` is in the row order of `y`.
  - Stop with "collCovariates are not used by a count model; differences between samples enter through countModel$offset" when `collCovariates` is not empty.
  - `default_b0_prior()` gives "community" for "counts" with more than one species.
  - Build `list_data$count` and the site totals `z` only when the model is "counts". Nothing is added to `list_data` otherwise.
  - Change the one quoted error message in the two Lesson 2 files to the new wording and nothing else in them.

- [ ] **Step 4: Run the test file and the unchanged check.** Expected: PASS and four lines ending "identical".

- [ ] **Step 5: Commit**

```bash
git add R/runOccJSDM.R vignettes/occJSDM-lesson-2.Rmd vignettes/occJSDM-lesson-2.md tests/testthat/test-count-model.R
git commit -m "Take count data into runOccJSDM when countModel is set"
```

### Task 5: The count branch of the coefficient update

**Files:**
- Modify: `R/count_pieces.R` (new function)
- Modify: `R/jsdmfun.R` (`update_jSDMcoef()`, lines 1243 to 1495)
- Test: `tests/testthat/test-count-model.R`

**Interfaces:**
- Consumes: `samplePGcount_cpp(shape, phi)` from stage 1; `list_data$count` from Task 4.
- Produces: `count_working_response(y, eta, size, site_index, log_effort)`, internal, returning `list(k, Omega)`, two unnamed `n` by `S` matrices.
- Produces: `update_jSDMcoef(..., model = "counts")` works. It reads `list_params$size` and returns it unchanged as element `size`, returns `psi = NULL`, and returns `eta` as the log of the expected count per unit of effort.

- [ ] **Step 1: Write the failing test**

```r
test_that("count_working_response sums the Polya-Gamma terms of a site's samples", {
  set.seed(3)
  site_index <- c(1L, 1L, 2L, 3L, 3L, 3L)
  y <- matrix(stats::rnbinom(12, size = 2, mu = 4), 6, 2)
  eta <- matrix(stats::rnorm(6), 3, 2)
  size <- c(0.7, 3)
  log_effort <- log(c(1, 2, 1, 0.5, 1, 4))
  setOccJSDMSeed(5L)
  got <- occJSDM:::count_working_response(y, eta, size, site_index, log_effort)
  o <- outer(log_effort, log(size), "-")
  setOccJSDMSeed(5L)
  omega <- samplePGcount_cpp(sweep(y, 2, size, "+"), eta[site_index, ] + o)
  expect_equal(got$Omega, unname(rowsum(omega, site_index)))
  expect_equal(got$k, unname(rowsum(sweep(y, 2, size, "-") / 2 - omega * o, site_index)))
  expect_identical(dim(got$k), c(3L, 2L))
  expect_true(all(got$Omega > 0))
})
```

- [ ] **Step 2: Run and see it fail.**

- [ ] **Step 3: Write `count_working_response()`.** For sample row `m` at site `i` and species `j`: the offset is `o = log_effort[m] - log(size[j])`, the predictor is `phi = eta[i, j] + o`, the weight is one Polya-Gamma draw of shape `y + size[j]` at `phi`, and the row contributes the weight to `Omega[i, j]` and `(y - size[j]) / 2 - weight * o` to `k[i, j]`.

- [ ] **Step 4: Wire it into `update_jSDMcoef()`.** Five changes, and nothing else in the function:

  - Where `k` is defined (line 1311), the empty stub for counts goes; `k` and `Omega` for counts are set together after `psiCoef` is computed, from `count_working_response(list_data$count$y, psiCoef, list_params$size, list_data$count$site_index, list_data$count$log_effort)`.
  - Define once, near the top, `gaussian_model <- if (model == "counts") "binary" else model`, and pass `gaussian_model` in place of `model` to `sample_BBsL_cpp()` and `sample_U_cpp()`. The existing `kappa <- if (model == "continuous") k*Omega else k` is already right for counts.
  - `computeVariancePartitioning_stdevs()` keeps receiving `model`, so counts are partitioned on the scale of `exp(eta)` (question 12).
  - `psi` is `NULL` for counts.
  - The returned list gains `"size" = list_params$size`.

- [ ] **Step 5: Run the test file and the unchanged check.** Expected: PASS and four lines ending "identical".

- [ ] **Step 6: Commit**

```bash
git add R/count_pieces.R R/jsdmfun.R tests/testthat/test-count-model.R
git commit -m "Form the working response of a count fit from Polya-Gamma draws"
```

### Task 6: The size update, storage and the first count fit

**Files:**
- Modify: `R/chains.R` (`run_one_chain()`: state near line 589, the loop after line 700, storage near lines 490 and 787, the WAIC call at line 797, the return near line 837)
- Modify: `R/runOccJSDM.R` (remove Task 4's stop; output assembly near lines 1183, 1286, 1365 and 1409)
- Modify: `src/functions.cpp` (`computeModelLoglikJSDM_cpp()`: delete the Poisson "counts" branch, lines 1187 to 1193)
- Modify: `tests/testthat/helper-fixtures.R` (`fit_fixture()` gains `countModel = NULL`, passed on)
- Test: `tests/testthat/test-count-model.R`

**Interfaces:**
- Consumes: `count_log_mean()`, `count_size_step()`, `count_size_initial_step`, `nbLogPmf_cpp()` from stage 1.
- Produces, in a count fit only:
  - `results_output$jsdm_output$size_output`: species by iterations by chains.
  - `results_output$jsdm_output$zeros_output`: species by iterations by chains, the number of zero counts in one replicate dataset drawn at each kept iteration.
  - `infos$count`: `list(family, occupancy, offset, y, site_index, log_effort, size_prior_sd, observed_zeros, size_step, size_acceptance)`. `size_step` and `size_acceptance` are species by chains; the acceptance rate is over kept iterations.
- A fit of any other model has none of these elements.

- [ ] **Step 1: Write the failing tests**

```r
count_fit <- function(sim, ...) fit_fixture(sim, spatial = FALSE,
  countModel = list(occupancy = FALSE, offset = "Effort"), ...)

test_that("a count fit runs, stores the size and reproduces from its seed", {
  sim <- simulate_fixture(model = "counts", n = 40, S = 4, M = 2, useSpatField = FALSE, seed = 8)
  sim$data_list$OTU[, 1] <- 0
  sim$data_list$OTU[, 2] <- sim$data_list$OTU[, 2] + 400
  set.seed(21); fit <- count_fit(sim)
  jsdm <- fit$results_output$jsdm_output
  expect_identical(fit$infos$model, "counts")
  expect_identical(dim(jsdm$size_output), c(4L, FIXTURE_MCMC$niter, FIXTURE_MCMC$nchain))
  expect_true(all(is.finite(jsdm$size_output)) && all(jsdm$size_output > 0))
  expect_true(all(is.finite(jsdm$B0_output)))
  expect_true(is.finite(fit$results_output$WAIC))
  expect_null(fit$results_output$psi_output)
  expect_identical(fit$infos$count$observed_zeros[1], 80L)
  expect_identical(dim(fit$infos$count$size_acceptance), c(4L, FIXTURE_MCMC$nchain))
  set.seed(21); again <- count_fit(sim)
  expect_identical(again$results_output, fit$results_output)
})

test_that("the effort offset changes a count fit and other models gain nothing", {
  sim <- simulate_fixture(model = "counts", n = 40, S = 3, M = 1, useSpatField = FALSE, seed = 9)
  set.seed(5); with_effort <- count_fit(sim)
  set.seed(5); without <- fit_fixture(sim, spatial = FALSE, countModel = list(occupancy = FALSE))
  expect_false(identical(with_effort$results_output$jsdm_output$B0_output,
                         without$results_output$jsdm_output$B0_output))
  binary <- fit_fixture(simulate_fixture(model = "binary", seed = 2))
  expect_null(binary$results_output$jsdm_output$size_output)
  expect_null(binary$results_output$jsdm_output$zeros_output)
  expect_null(binary$infos$count)
})

test_that("a count fit with one species uses the fixed intercept prior", {
  set.seed(3)
  info <- data.frame(Site = 1:30, elev = stats::rnorm(30))
  y <- matrix(stats::rnbinom(30, size = 2, mu = 5), 30, 1, dimnames = list(NULL, "s1"))
  fit <- suppressMessages(runOccJSDM(list(info = info, OTU = y), occCovariates = "elev",
                                     countModel = list(occupancy = FALSE), MCMCparams = FIXTURE_MCMC))
  expect_null(fit$infos$intercept_prior$type)
  expect_true(all(is.finite(fit$results_output$jsdm_output$size_output)))
})
```

- [ ] **Step 2: Run and see them fail** with "count sampler not wired yet".

- [ ] **Step 3: Implement the chain.** For a count fit only:

  - Before the loop: `list_jSDMparams$size <- rep(1, S)`; a step of `count_size_initial_step` per species; an integer adaptation count of 0; acceptance counts per species for burn-in and for kept iterations.
  - In the loop, directly after `update_jSDMcoef()` returns: compute `log_mean <- count_log_mean(eta, site_index, log_effort)`, call `count_size_step()` with every cell included, `adapt = iter <= nburn`, and write the new size back to `list_jSDMparams$size`. One proposal per species per iteration. The Polya-Gamma weights are drawn afresh from the new size at the top of the next `update_jSDMcoef()`, which is what makes updating the size with the weights summed out valid: no step may use the old weights after this one.
  - On a kept iteration: store the size; draw the replicate zeros as `colSums` of `stats::rbinom()` with success probability `exp(nbLogPmf_cpp(0 * y, log_mean, size))`; and for WAIC use `as.vector(nbLogPmf_cpp(y, log_mean, list_jSDMparams$size))`, with accumulators of length `N * S`. Read the size from `list_jSDMparams`, not from a local set before the loop.
  - Return the size and zero draws and the adaptation record with the chain's other results.

- [ ] **Step 4: Implement the assembly in `runOccJSDM()`.** Allocate, fill and store `size_output`, `zeros_output` and `infos$count` only when the model is "counts". `observed_zeros` is an integer vector. Remove Task 4's stop. Keep the test hook.

- [ ] **Step 5: Delete the Poisson branch** of `computeModelLoglikJSDM_cpp()` so that "counts" reaches "Unknown model type" there, and run `Rcpp::compileAttributes()`.

- [ ] **Step 6: Run the test file, the unchanged check and the whole suite** (the suite in the background). Expected: PASS; four lines ending "identical"; the baseline's fail and skip counts.

- [ ] **Step 7: Commit**

```bash
git add R/chains.R R/runOccJSDM.R src/functions.cpp src/RcppExports.cpp R/RcppExports.R tests/testthat/helper-fixtures.R tests/testthat/test-count-model.R
git commit -m "Fit the plain count JSDM: size update, storage and WAIC"
```

### Task 7: Agreement with an independent fit, and recovery

**Files:**
- Modify: `DESCRIPTION` (add `MASS` to `Suggests`)
- Modify: `tests/testthat/helper-simstudy.R` (`draw_truth()`, `simstudy_fit()`, `simstudy_param_blocks()`)
- Create: `tests/testthat/test-count-recovery.R`

**Interfaces:**
- Produces, in the helper: a scenario with `model = "counts"` simulates with `rnb` and an effort, fits with `countModel = list(occupancy = FALSE, offset = "Effort")`, and `simstudy_param_blocks()` gains a block named "size" that compares `size_output` with `true_params$size_true`.

- [ ] **Step 1: Write the comparison test.** No covariates, no factors, a wide fixed intercept prior, so that the posterior should sit on the maximum-likelihood fit.

```r
test_that("an intercept-only count fit agrees with MASS::glm.nb", {
  skip_on_cran()
  skip_if_not_installed("MASS")
  set.seed(31)
  n <- 200
  site <- rep(seq_len(n), each = 2)
  effort <- stats::runif(2 * n, 0.5, 3)
  means <- c(8, 5, 3); sizes <- c(0.8, 1.5, 3)
  y <- sapply(1:3, function(s) stats::rnbinom(2 * n, size = sizes[s], mu = effort * means[s]))
  colnames(y) <- paste0("sp", 1:3)
  fit <- suppressMessages(runOccJSDM(
    list(info = data.frame(Site = site, days = effort), OTU = y),
    countModel = list(occupancy = FALSE, offset = "days"),
    listPriors = list(b0_prior = "fixed", sigma_b0 = 10),
    MCMCparams = list(nchain = 2, nburn = 1500, niter = 1500, nthin = 1)))
  jsdm <- fit$results_output$jsdm_output
  for (s in 1:3) {
    ml <- MASS::glm.nb(y[, s] ~ 1 + offset(log(effort)))
    b0 <- as.vector(jsdm$B0_output[s, , ])
    log_size <- log(as.vector(jsdm$size_output[s, , ]))
    expect_lt(abs(mean(b0) - stats::coef(ml)[[1]]), 3 * stats::sd(b0) / sqrt(50) + 0.03)
    expect_lt(abs(mean(log_size) - log(ml$theta)), 3 * stats::sd(log_size) / sqrt(50) + 0.08)
  }
})
```

- [ ] **Step 2: Write the recovery test,** on the pattern of the canary in `tests/testthat/test-recovery.R` (lines 31 to 73).

```r
.count_recovery_cache <- new.env(parent = emptyenv())
count_recovery_rows <- function() {
  if (is.null(.count_recovery_cache$bundle)) {
    scenario <- simstudy_scenarios()[[1]]
    scenario$model <- "counts"
    scenario$n <- 80L
    scenario$S <- 6L
    scenario$M <- 2L
    scenario$useSpatField <- FALSE
    scenario$rnb <- rep(c(0.8, 3), 3)
    .count_recovery_cache$bundle <- simstudy_scenario(
      scenario, R = 4L, verbose = FALSE,
      MCMCparams = list(nchain = 2, nburn = 600, niter = 600, nthin = 1))
  }
  .count_recovery_cache$bundle$params
}

test_that("a count fit recovers its coefficients and sizes (canary)", {
  skip_on_cran()
  rows <- count_recovery_rows()
  expect_true("size" %in% rows$block)
  expect_gt(mean(rows$covered), 0.70)
  expect_gt(min(simstudy_summarise(rows)$coverage), 0.40)
  expect_gt(stats::cor(rows$post_mean, rows$truth), 0.30)
})
```

- [ ] **Step 3: Extend the helper** until both tests run. `model`, `n`, `S`, `M` and `useSpatField` are existing fields of `SIMSTUDY_BASE`; `rnb` is new, and a count scenario draws its effort from `stats::runif(sum(M), 0.5, 3)`. If the name of the block column in `simstudy_rows()` is not `block`, use the helper's own name for it in the test. The floors are the canary's: they catch a broken model, not a slightly miscalibrated one, and Task 11 reports the real numbers.

- [ ] **Step 4: Run both tests.** Expected: PASS. If the comparison test fails, the likelihood wiring is wrong: do not loosen its tolerances. If the recovery floors fail, report the numbers to Doug before changing anything.

- [ ] **Step 5: Commit**

```bash
git add DESCRIPTION tests/testthat/helper-simstudy.R tests/testthat/test-count-recovery.R
git commit -m "Check count fits against glm.nb and for recovery"
```

### Task 8: Post-fit functions

**Files:**
- Modify: `R/count_model.R` (`stop_if_count_fit()`, `returnZeroCheck()`)
- Modify: `R/output.R` (`returnCovariateEffect()` line 573, `plotCovariateEffect()` line 624, and six functions that must refuse)
- Modify: `R/jsdmfun.R` (`returnCovariateEffect_base()` line 372, `plot_covariate_response()` line 403)
- Modify: `R/posterior-draws.R` (`returnPosteriorDraws()`, lines 57, 75 and 86)
- Modify: `R/diagnostics.R` (`returnConvergenceDiagnostics()` near line 325; the exclusion list in `computeDiagnostics()` at line 492)
- Modify: `src/jsdm.cpp` (`computeNewOutputs()`, lines 489 to 493)
- Create: `tests/testthat/test-count-output.R`

**Interfaces:**
- Produces: `returnZeroCheck(fitModel)`, exported. A data frame with one row per species: `species`, `observed`, `predicted_mean`, `predicted_lower`, `predicted_upper` (the 2.5% and 97.5% quantiles of the replicate zeros) and `p_more_zeros`, the share of replicates with at least as many zeros as were observed.
- Produces: `stop_if_count_fit(fitModel, what)`, internal.

- [ ] **Step 1: Write the failing tests.** Build one count fit (`count_fit_obj`) and one binary fit (`binary_fit_obj`) at the top of the file with `fit_fixture()`. Where a test calls an existing function, its arguments follow that function's existing tests in `tests/testthat/test-api-contracts.R`; what this plan fixes is the assertions.

```r
test_that("probability-scale functions refuse a count fit and say why", {
  for (f in list(returnOccupancyRates, plotOccupancyRates, computePredictiveOccupancyProbs,
                 computeConditionalOccupancyProbs))
    expect_error(f(count_fit_obj), "not defined for a count fit")
})

test_that("the size is a named posterior parameter", {
  draws <- returnPosteriorDraws(count_fit_obj, "size")
  values <- if (is.data.frame(draws)) unlist(draws[vapply(draws, is.numeric, logical(1))]) else draws
  expect_true(all(values > 0))
  expect_error(returnPosteriorDraws(binary_fit_obj, "size"), "size")
  expect_true(any(grepl("size", unlist(returnConvergenceDiagnostics(count_fit_obj)))))
})

test_that("covariate effects of a count fit are expected counts", {
  effect <- returnCovariateEffect(count_fit_obj, covName = "X_psi.EnvCov.1", idx_species = 1:2)
  numeric_cols <- vapply(effect, is.numeric, logical(1))
  expect_true(all(unlist(effect[numeric_cols]) >= 0))
  expect_true(max(unlist(effect[numeric_cols])) > 1)
})

test_that("returnZeroCheck compares observed zeros with replicates", {
  check <- returnZeroCheck(count_fit_obj)
  expect_identical(nrow(check), count_fit_obj$infos$S)
  expect_identical(check$observed, count_fit_obj$infos$count$observed_zeros)
  expect_true(all(check$predicted_lower <= check$predicted_upper))
  expect_true(all(check$p_more_zeros >= 0 & check$p_more_zeros <= 1))
  expect_error(returnZeroCheck(binary_fit_obj), "count fit")
})
```

Add a test of `predictNewSites()` on the count fit by copying "predictNewSites() handles a fit with no spatial field" (`tests/testthat/test-api-contracts.R` near line 189) and asserting on its numeric output: every predicted value is finite and at least 0, and at least one exceeds 1. Add tests that `returnOccupancyGradient()` and `plotOccupancyGradient()` refuse a count fit, called with the arguments their existing tests use.

- [ ] **Step 2: Run and see them fail.**

- [ ] **Step 3: Implement.**

  - `stop_if_count_fit()` raises, with `call. = FALSE`: the function's name followed by " is not defined for a count fit: its predictor is the log of the expected count, not an occupancy probability". Call it first in `returnOccupancyRates()`, `plotOccupancyRates()`, `returnOccupancyGradient()`, `plotOccupancyGradient()`, `computePredictiveOccupancyProbs()` and `computeConditionalOccupancyProbs()`.
  - The covariate-effect functions choose the link from `infos$jsdmModel`: "identity" for "continuous", "log" for "counts", "logit" otherwise. `returnCovariateEffect_base()` accepts "log" and returns `exp(eta)`. `plot_covariate_response()` labels a log link "Expected count per unit of effort".
  - `computeNewOutputs()` gains a "counts" branch that stores `arma::exp(linpred)`, and a final `else` that stops with "Unknown model". Today an unknown model returns uninitialised memory.
  - `returnPosteriorDraws()` accepts "size" as a species-only parameter read from `jsdm_output$size_output`. On a fit without it the existing "not available" handling applies; if there is none, stop with "size is stored only in count fits".
  - `returnConvergenceDiagnostics()` gains a block for the size, on the pattern of the `theta0` block.
  - `computeDiagnostics()` skips `zeros_output`, which is not a parameter.
  - `returnZeroCheck()` stops on anything else with "returnZeroCheck needs a count fit".

- [ ] **Step 4: Run `Rcpp::compileAttributes()`, the test file and the unchanged check.**

- [ ] **Step 5: Commit**

```bash
git add R/count_model.R R/output.R R/jsdmfun.R R/posterior-draws.R R/diagnostics.R src/jsdm.cpp src/RcppExports.cpp R/RcppExports.R tests/testthat/test-count-output.R
git commit -m "Teach the post-fit functions about count fits"
```

### Task 9: Site-level WAIC for count fits

**Files:**
- Modify: `src/site_waic.cpp` (new exported function; include `count.h`)
- Modify: `R/site-waic.R` (`computeSiteWAIC()`: the model check at line 121, the response checks near line 143, the per-draw evaluation near line 232)
- Test: `tests/testthat/test-site-waic.R`

**Interfaces:**
- Consumes: `nb_logpmf()` from `src/count.h`; `infos$count` and `size_output` from Task 6.
- Produces: `site_loglik_count_cpp(eta, y, row_site, log_effort, size, loadings, nodes, log_weights)`. `eta` is sites by species; `y` is sample rows by species; `row_site` gives each row's site as a 1-based row of `eta`; the rest are as in `site_loglik_cpp()`. It returns one log-likelihood per site: the log of the sum over nodes of the node's weight times the product, over species and over the site's rows, of the negative-binomial probability at log mean `eta + shift + log_effort`.
- Produces: `computeSiteWAIC()` accepts a non-spatial count fit.

- [ ] **Step 1: Write the failing tests** in `tests/testthat/test-site-waic.R`. `waic_count_fit()` builds a fake count fit by hand on the pattern of `waic_test_fit()` (lines 3 to 17): `B0 = 0`, loadings 0.8, `sigmah = 1`, sizes 2, three iterations, one chain, empty `X_psi`, and `infos$count` holding `y`, `site_index` and `log_effort`.

```r
test_that("count site likelihood integrates the factor over all of a site's rows", {
  y <- matrix(c(3, 0, 5, 1, 2, 2), 3, 2)
  fit <- waic_count_fit(y, site_index = c(1L, 1L, 2L), effort = c(1, 2, 0.5), factors = 1)
  score <- computeSiteWAIC(fit)
  site1 <- function(u) vapply(u, function(v) prod(stats::dnbinom(
    y[1:2, ], size = 2, mu = c(1, 2) * exp(0.8 * v))) * stats::dnorm(v), numeric(1))
  expect_equal(score$log_lik[1, 1],
               log(stats::integrate(site1, -Inf, Inf, rel.tol = 1e-10)$value), tolerance = 1e-6)
})

test_that("with no factors the count site likelihood is a sum of log densities", {
  y <- matrix(c(3, 0, 5, 1, 2, 2), 3, 2)
  fit <- waic_count_fit(y, site_index = c(1L, 1L, 2L), effort = c(1, 2, 0.5), factors = 0)
  score <- computeSiteWAIC(fit)
  want <- c(sum(stats::dnbinom(y[1:2, ], size = 2, mu = c(1, 2), log = TRUE)),
            sum(stats::dnbinom(y[3, ], size = 2, mu = 0.5, log = TRUE)))
  expect_equal(unname(score$log_lik[1, ]), want, tolerance = 1e-10)
})

test_that("a real count fit can be scored and a spatial one is refused as before", {
  sim <- simulate_fixture(model = "counts", n = 30, S = 3, M = 2, useSpatField = FALSE, seed = 6)
  fit <- fit_fixture(sim, spatial = FALSE, countModel = list(occupancy = FALSE, offset = "Effort"))
  score <- computeSiteWAIC(fit)
  expect_true(is.finite(score$WAIC))
  expect_identical(nrow(score$pointwise), 30L)
})
```

- [ ] **Step 2: Run and see them fail** with "Site WAIC supports binary, occupancy and two_stage models only."

- [ ] **Step 3: Implement.** In `computeSiteWAIC()`, a count fit takes its response, row-to-site map and efforts from `infos$count`, pulls `size_output` with the function's own array helper, skips the threshold and the 0-or-1 check, and calls the new integrator in place of `site_loglik_cpp()`. The quadrature ladder, the convergence test, the summaries and the returned object are unchanged. The `observations` element gains `log_effort`, so that `compareSiteWAIC()` refuses fits scored at different efforts. Change the message of the model check to name counts as well.

- [ ] **Step 4: Run `Rcpp::compileAttributes()` and the whole of `test-site-waic.R`.** Expected: PASS, including every existing test.

- [ ] **Step 5: Commit**

```bash
git add src/site_waic.cpp src/RcppExports.cpp R/RcppExports.R R/site-waic.R tests/testthat/test-site-waic.R
git commit -m "Score count fits with the site-level WAIC"
```

### Task 10: Help pages

**Files:**
- Modify: `R/runOccJSDM.R` (roxygen: `@param countModel`, the `size_prior_sd` prior, one line pointing to `?countModel`)
- Modify: `R/count_model.R` (a documentation-only topic named `countModel`; roxygen for `returnZeroCheck()`)
- Modify: `R/output.R` and `R/posterior-draws.R` (the roxygen sentences that say "occupancy probability" or list the parameters, where a count fit now differs)
- Modify: `_pkgdown.yml` (add `countModel` and `returnZeroCheck` to the reference index)
- Modify: `NEWS.md`
- Regenerate: `man/`, `NAMESPACE`

- [ ] **Step 1: Write the `countModel` help topic.** It must say, in this order: what the list holds and that `occupancy` has no default; that with `occupancy = FALSE` the `occCovariates` act on the log of the expected count per unit of effort, despite the argument's name; that the model is for counts of individuals and must not be used for sequence read counts, which are not counts of organisms; that without `countModel`, counts with repeated rows per site are fitted as an occupancy model and cut to detections at `threshold`, with no warning; that `collCovariates` and missing counts are not accepted; that the shared-predictor model, `occupancy = TRUE`, is not available yet; and what a count fit stores (`size_output`, `infos$count`) and which functions read it (`returnPosteriorDraws()`, `returnZeroCheck()`, `computeSiteWAIC()`).

- [ ] **Step 2: Regenerate and check.**

Run: `Rscript -e 'roxygen2::roxygenise()'`
Expected: `man/countModel.Rd` and `man/returnZeroCheck.Rd` are new; `NAMESPACE` gains `export(returnZeroCheck)` and nothing else changes in it.

Run: `git status --short man NAMESPACE`
Expected: only help pages of functions touched in this plan.

- [ ] **Step 3: Add a `NEWS.md` entry** under the development version: counts of individuals can be fitted with `countModel = list(occupancy = FALSE)`; fits of the other models are unchanged.

- [ ] **Step 4: Commit**

```bash
git add R man NAMESPACE _pkgdown.yml NEWS.md
git commit -m "Document countModel and the count outputs"
```

### Task 11: The stage 2 measurements

**Files:**
- Create: `dev/count-model/stage2-checks.R`
- Create: `dev/count-model/stage2-report.md`

- [ ] **Step 1: Write the script.** One process, no parallel chains (`ncores` 1). It starts with `devtools::load_all(quiet = TRUE)` and fixed seeds, prints the power state as the stage 1 script does, and runs:

  - **Recovery by rarity in a full fit.** Simulate 250 sites by 40 species, one sample per site, two covariates, two latent factors, no spatial field, sizes alternating 0.8 and 3, with species' intercepts spread with a standard deviation of 2. Fit with 2 chains of 3,000 burn-in and 3,000 kept. By species' number of non-zero counts (under 10, 10 to 49, 50 to 149, 150 and over): the error of the posterior mean of the intercept, the coverage of its 95% interval, and the same for the log size.
  - **The intercept prior's spread.** The same data fitted again with `listPriors = list(sigma_b0 = 2)`. Report the intercept error by rarity under both.
  - **Mixing.** The effective sample size and R-hat of the size, the intercepts, the coefficients and `mu_b0`, from `returnConvergenceDiagnostics()`; the acceptance rates in `infos$count$size_acceptance`.
  - **The zero check.** The output of `returnZeroCheck()`: how many species have `p_more_zeros` below 0.025 or above 0.975.
  - **Cost.** Seconds per 1,000 iterations of the count fit, and of a binary fit of the same data cut to presence and absence with the same settings.

- [ ] **Step 2: Run it detached and do not edit it while it runs.**

Run: `nohup Rscript dev/count-model/stage2-checks.R > /tmp/stage2-checks.log 2>&1 &`
Expected: the log ends with the cost block. Check that the output is there, not only that the process exited.

- [ ] **Step 3: Write `dev/count-model/stage2-report.md`** from the log, as bullets with no tables: the power state and what else was running, each block's numbers, what the run does not show, and plain statements of whether the sizes of rare species are recovered, whether `sigma_b0 = 1` costs anything when the true spread is 2, and whether any block failed to converge. Say "fails" where something fails.

- [ ] **Step 4: Commit**

```bash
git add dev/count-model/stage2-checks.R dev/count-model/stage2-report.md
git commit -m "Measure recovery, mixing and cost of the plain count JSDM"
```

### Task 12: Confirm nothing else changed

**Files:** none modified.

- [ ] **Step 1: The unchanged check, one last time.**

Run: `Rscript dev/count-model/stage2-unchanged.R check /tmp/occjsdm-stage2-reference`
Expected: four lines ending "identical", exit status 0.

- [ ] **Step 2: The whole suite, in the background.**

Run: `Rscript -e 'devtools::test()' 2>&1 | tail -n 8`
Expected: the baseline's fail and skip counts from Task 1, and its pass count plus the new expectations.

- [ ] **Step 3: The untouched functions.**

Run: `git diff "$BASE" -- src/jsdm.cpp | grep -c '^-[^-]'`, where `BASE` is the commit this branch was cut from: `origin/main` if stage 1 had merged, otherwise the tip of the stage 1 branch.
Expected: at most the two lines of `computeNewOutputs()` that gained an `else`. `sample_BBsL_cpp()`, `sample_U_cpp()` and `samplepg_fast()` have no removed or changed line.

- [ ] **Step 4: Hand over to Doug.** Report the suite counts against the baseline, the report's findings, and every item under *Decisions made in this plan*. Do not push. The pull request is opened only when Doug says so; its description opens with a line saying Claude wrote it and writes out every question still open in full.

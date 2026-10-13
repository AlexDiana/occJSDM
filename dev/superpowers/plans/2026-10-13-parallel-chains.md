# Parallel MCMC Chains Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run `runOccJSDM()`'s MCMC chains on separate cores through base `parallel` PSOCK workers, with per-chain seeds, an explicit core budget, exact combination of draws, latent-state means and WAIC accumulators, and a fit that is identical whatever the worker count.

**Architecture:** The body of the chain loop moves into one internal function, `run_one_chain()`, which takes a prepared `setup` list and a chain's two seeds. `runOccJSDM()` builds `setup`, draws the seeds, then runs chains one at a time in process or in batches on a PSOCK cluster through `clusterApplyLB()`, copying each returned chain into the full arrays and merging per-chain sums and WAIC partials in chain order. Every chain's randomness is fixed by its own seeds, so which process ran it cannot matter.

**Tech Stack:** R 4.6.1 (system `Rscript`, `/usr/local/bin/Rscript`), base `parallel` and `tools`, testthat 3e, roxygen2, pkgload for the development-mode worker branch only. No C++ changes.

**Spec:** `dev/superpowers/specs/2026-10-13-parallel-chains-design.md` (revision `bc3694b` and later). The plan argues from the spec; read both.

## Global Constraints

- Work only in the worktree `~/src/occJSDM/.worktrees/parallel-chains`, branch `claude/parallel-chains`. Never switch branches in the main checkout, which is open in RStudio and holds Doug's uncommitted edits.
- Never run `R CMD build`, `R CMD check`, `devtools::check()` or `devtools::build()` on a checkout: they copy the whole tree first and nearly filled the disk once. Build and install from a `git archive` export (Task 6 shows the exact commands).
- Base `parallel` only. `pkgload` is used only inside a worker when the parent itself was loaded by pkgload. No new dependency beyond adding `parallel` and `tools` to `Imports`.
- R CMD check allows at most two cores: every always-on test uses at most two workers. The always-on tier-1 tests must stay within the 30-second budget in `tests/testthat/helper-fixtures.R`; move cases under `skip_on_cran()` if they do not.
- Markdown in this repository: one line per paragraph, no hard wrapping, no pipe tables, no em-dashes anywhere (use `--`), escape `->` as `-\>` and `~` as `\~` in prose, never inside code spans or fences.
- Commit messages: imperative first line, no em-dashes, end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`.
- Long jobs (benchmark, refit) run detached with `nohup` plus `setsid` from a non-interactive shell, mains power, lid open; they die with the tool session otherwise.
- The shell is zsh: quote `ls` as `'ls'` (it is aliased), and remember unquoted variables are not word-split.
- The hashed set for the refit provenance is `git ls-files R src DESCRIPTION NAMESPACE`. Nothing in it changes after Task 6 without redoing Task 8.
- Sampling order within an iteration, every sampler, every prior and every default other than `ncores` stay unchanged.

## Review Focus

- A fit with `nchain = 1` and `ncores > 1` stays in process: no cluster, no progress file, the session's RcppParallel count set to `ncores` for the fit and restored afterwards, and the execution record saying so. Test added to Task 5.
- `summarisedLatentPresences = FALSE` with workers: the full `z` and `psi` arrays cross the socket and land in the right chain slot; the fit is identical to the in-process fit. Test added to Task 5.
- A session in a non-default RNG kind (for example `L'Ecuyer-CMRG` after `mclapply()`, with Box-Muller normals) gets the same fit at every `ncores` and has its kinds back afterwards. Test added to Task 5.
- `ncores` given as a double such as `2`, or larger than `nchain`, is accepted and the budget remainder goes to threads, never to extra workers. Tests added to Tasks 1 and 5.
- After any exit, normal, error or interrupt, nothing leaks: the setup environment is empty, the BLAS thread variables and the RcppParallel variable are back to their previous values, and no worker process survives. Tests added to Tasks 4 and 5.

---

### Task 1: Scaffolding, `ncores` validation and the core budget

**Files:**
- Create: `R/chains.R`
- Create: `tests/testthat/test-parallel-chains.R`
- Modify: `DESCRIPTION` (the `Imports:` field)
- Modify: `R/occJSDM-package.R` (the `@importFrom` lines)
- Modify: `R/globals.R` (the `globalVariables()` vector)

**Interfaces:**
- Consumes: `get_param(params, key, default)` from `R/runOccJSDM.R:1`.
- Produces: `chain_env` (package-level environment); `resolve_ncores(MCMCparams)` returning an integer; `plan_chain_budget(nchain, ncores)` returning `list(ncores, workers, threads_per_worker, thread_setting)` with integer fields and `thread_setting` either `"ambient"` or `"budget"`; `ambient_rcpp_threads()` returning an integer.

- [ ] **Step 1: Write the failing tests**

Create `tests/testthat/test-parallel-chains.R` with this content:

```r
# Parallel MCMC chains. Design: dev/superpowers/specs/2026-10-13-parallel-chains-design.md.
# Rules of this file: at most two workers (R CMD check's limit), structural
# assertions in the always-on tier, fit identity only where the spec allows it.

alive <- function(pids) vapply(pids, function(p) isTRUE(tools::pskill(p, signal = 0L)), logical(1))

test_that("MCMCparams$ncores is validated", {
  for (bad in list(0, 1.5, "2", c(1, 2), NA_real_, Inf)) {
    expect_error(occJSDM:::resolve_ncores(list(ncores = bad)), "single whole number")
  }
  expect_identical(occJSDM:::resolve_ncores(list()), 1L)
  expect_identical(occJSDM:::resolve_ncores(list(ncores = 2)), 2L)      # a double is accepted
  expect_identical(occJSDM:::resolve_ncores(list(ncores = 3L)), 3L)
  expect_warning(occJSDM:::resolve_ncores(list(ncores = parallel::detectCores() + 1L)), "exceeds")
})

test_that("the core budget splits into workers and threads per worker", {
  b <- occJSDM:::plan_chain_budget
  pick <- function(x) x[c("workers", "threads_per_worker", "thread_setting")]
  expect_identical(pick(b(4L, 1L)), list(workers = 1L, threads_per_worker = 1L, thread_setting = "ambient"))
  expect_identical(pick(b(4L, 4L)), list(workers = 4L, threads_per_worker = 1L, thread_setting = "budget"))
  expect_identical(pick(b(4L, 8L)), list(workers = 4L, threads_per_worker = 2L, thread_setting = "budget"))
  expect_identical(pick(b(2L, 3L)), list(workers = 2L, threads_per_worker = 1L, thread_setting = "budget"))
  expect_identical(pick(b(1L, 4L)), list(workers = 1L, threads_per_worker = 4L, thread_setting = "budget"))
  expect_identical(b(4L, 8L)$ncores, 8L)
})

test_that("the ambient RcppParallel count follows the environment variable", {
  old <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS") else Sys.setenv(RCPP_PARALLEL_NUM_THREADS = old))
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "3")
  expect_identical(occJSDM:::ambient_rcpp_threads(), 3L)
  Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS")
  expect_identical(occJSDM:::ambient_rcpp_threads(), as.integer(RcppParallel::defaultNumThreads()))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

Run from the worktree root:

```bash
Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-parallel-chains.R")'
```

Expected: three failures, each "object 'resolve_ncores' not found" or similar.

- [ ] **Step 3: Create `R/chains.R` with the budget helpers**

```r
# Parallel MCMC chains: core budget, per-chain seeds, the single-chain
# function, the PSOCK driver and the combination helpers.
# Design: dev/superpowers/specs/2026-10-13-parallel-chains-design.md.

# Holds the prepared setup list while chains run, in the parent and in every
# worker. Package level, so that a worker can reach it by name after the task
# function arrives as a reference to the namespace.
chain_env <- new.env(parent = emptyenv())

resolve_ncores <- function(MCMCparams) {
  ncores <- get_param(MCMCparams, "ncores", 1)
  if (!is.numeric(ncores) || length(ncores) != 1L || !is.finite(ncores) ||
      ncores < 1 || ncores != floor(ncores)) {
    stop("MCMCparams$ncores must be a single whole number of at least 1", call. = FALSE)
  }
  ncores <- as.integer(ncores)
  available <- parallel::detectCores()
  if (!is.na(available) && ncores > available) {
    warning(sprintf("MCMCparams$ncores = %d exceeds the %d cores detected on this machine",
                    ncores, available), call. = FALSE)
  }
  ncores
}

# ncores is the total core budget: workers = min(nchain, ncores), each with
# floor(ncores / workers) RcppParallel threads. ncores = 1 is today's fit and
# leaves the session's thread setting alone ("ambient").
plan_chain_budget <- function(nchain, ncores) {
  workers <- as.integer(min(nchain, ncores))
  list(ncores = as.integer(ncores),
       workers = workers,
       threads_per_worker = max(1L, as.integer(floor(ncores / workers))),
       thread_setting = if (ncores == 1L) "ambient" else "budget")
}

ambient_rcpp_threads <- function() {
  env <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = "")
  if (nzchar(env)) as.integer(env) else as.integer(RcppParallel::defaultNumThreads())
}
```

- [ ] **Step 4: Declare the imports**

In `DESCRIPTION`, add `parallel` and `tools` to `Imports:` (alphabetical placement is not required; append after `utils`):

```
    utils,
    parallel,
    tools
```

In `R/occJSDM-package.R`, add two lines beside the existing `@importFrom RcppParallel RcppParallelLibs`:

```r
#' @importFrom parallel makeCluster clusterCall clusterApplyLB stopCluster detectCores
#' @importFrom tools md5sum pskill
```

Regenerate `NAMESPACE` only (no roxygen comment on a documented object changes in this task):

```bash
Rscript -e 'roxygen2::roxygenise(roclets = "namespace")'
```

Then confirm `git diff NAMESPACE` shows only the two new `importFrom()` lines and that `DESCRIPTION`'s `RoxygenNote` did not change (if roxygen2 8.1.0 bumped it, revert that one line with `git checkout -- DESCRIPTION` and re-apply the `Imports:` edit; Task 6 settles the roxygen2 version).

- [ ] **Step 5: Run the tests to verify they pass**

```bash
Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-parallel-chains.R")'
```

Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add R/chains.R tests/testthat/test-parallel-chains.R DESCRIPTION NAMESPACE R/occJSDM-package.R
git commit -m "Add the core budget for parallel chains and validate MCMCparams\$ncores

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 2: WAIC and latent-sum combination helpers

**Files:**
- Modify: `R/chains.R`
- Modify: `tests/testthat/test-parallel-chains.R`

**Interfaces:**
- Consumes: `create_waic_quantities(n_obs)` (`R/runOccJSDM.R:236`), `update_waic_summary(logliks, list_waic, iter)` and `compute_waic(list_waic, numIters)` (`R/mcmcfun.R:35` and `:56`). A WAIC part is `list(M2, mean_log, mean_lik, n)`, the first three numeric vectors of one length, `n` the number of draws accumulated.
- Produces: `merge_waic_pair(a, b)`, `merge_waic_parts(parts)` (fold in list order, `NULL` parts ignored, all-`NULL` gives `NULL`), `add_latent_sums(total, part)` where each is `NULL` or a named list of matrices or `NULL`s.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-parallel-chains.R`:

```r
test_that("per-chain WAIC accumulators merge to the single-pass values", {
  set.seed(11)
  loglik <- matrix(rnorm(7 * 12, mean = -1), nrow = 7, ncol = 12)   # 7 observations, 12 draws
  single <- occJSDM:::create_waic_quantities(7)
  for (t in 1:12) single <- occJSDM:::update_waic_summary(loglik[, t], single, t)
  accumulate <- function(cols) {
    part <- occJSDM:::create_waic_quantities(7)
    for (i in seq_along(cols)) part <- occJSDM:::update_waic_summary(loglik[, cols[i]], part, i)
    c(part, list(n = length(cols)))
  }
  merged <- occJSDM:::merge_waic_parts(list(accumulate(1:5), accumulate(6:12)))   # unequal parts on purpose
  expect_equal(merged$n, 12L)
  expect_equal(merged$mean_log, single$mean_log, tolerance = 1e-12)
  expect_equal(merged$mean_lik, single$mean_lik, tolerance = 1e-12)
  expect_equal(merged$M2, single$M2, tolerance = 1e-12)
  expect_equal(occJSDM:::compute_waic(merged, 12), occJSDM:::compute_waic(single, 12), tolerance = 1e-12)
  three <- occJSDM:::merge_waic_parts(list(accumulate(1:4), accumulate(5:8), accumulate(9:12)))
  expect_equal(three$M2, single$M2, tolerance = 1e-12)
  expect_null(occJSDM:::merge_waic_parts(list(NULL, NULL)))
  expect_equal(occJSDM:::merge_waic_parts(list(NULL, accumulate(1:12)))$M2, single$M2, tolerance = 1e-12)
})

test_that("latent-state sums add elementwise and tolerate absent parts", {
  a <- list(z = matrix(1, 2, 2), psi = matrix(2, 2, 2), theta = NULL, w = NULL)
  b <- list(z = matrix(3, 2, 2), psi = matrix(4, 2, 2), theta = NULL, w = NULL)
  s <- occJSDM:::add_latent_sums(NULL, a)
  expect_identical(s, a)
  s <- occJSDM:::add_latent_sums(s, b)
  expect_equal(s$z, matrix(4, 2, 2))
  expect_equal(s$psi, matrix(6, 2, 2))
  expect_null(s$theta)
  expect_null(occJSDM:::add_latent_sums(NULL, NULL))
})
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-parallel-chains.R")'
```

Expected: the two new tests fail with "object 'merge_waic_parts' not found" and "object 'add_latent_sums' not found".

- [ ] **Step 3: Add the helpers to `R/chains.R`**

```r
# Merge two Welford accumulators (running mean, M2 and mean likelihood per
# observation, each over n draws) exactly, with the standard pairwise formula.
# Chains are merged in chain order by merge_waic_parts(), so the result cannot
# depend on which worker produced which chain.
merge_waic_pair <- function(a, b) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  n <- a$n + b$n
  delta <- b$mean_log - a$mean_log
  list(M2 = a$M2 + b$M2 + delta^2 * a$n * b$n / n,
       mean_log = a$mean_log + delta * b$n / n,
       mean_lik = (a$n * a$mean_lik + b$n * b$mean_lik) / n,
       n = n)
}

merge_waic_parts <- function(parts) {
  Reduce(merge_waic_pair, parts, accumulate = FALSE)
}

# Per-chain sums of the retained latent states, added in chain order.
add_latent_sums <- function(total, part) {
  if (is.null(part)) return(total)
  if (is.null(total)) return(part)
  out <- total
  for (nm in names(part)) {
    out[[nm]] <- if (is.null(total[[nm]]) || is.null(part[[nm]])) NULL else total[[nm]] + part[[nm]]
  }
  out
}
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-parallel-chains.R")'
```

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add R/chains.R tests/testthat/test-parallel-chains.R
git commit -m "Add exact merging of per-chain WAIC accumulators and latent-state sums

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 3: Extract `run_one_chain()`, per-chain seeds, and the in-process streaming loop

This is the large, mechanical move. After it, `runOccJSDM()` still runs every chain in process, but through the new architecture, with per-chain seeds, per-chain sums and WAIC partials, the RNG cleanup handler, and the three new `infos` fields. Every existing test must still pass.

**Files:**
- Modify: `R/chains.R`
- Modify: `R/runOccJSDM.R:1040-1424` (from `# run MCMC` through the end of the `# compute WAIC` block) and `R/runOccJSDM.R:1488-1517` (the `infos` list)
- Modify: `R/globals.R`
- Modify: `tests/testthat/test-parallel-chains.R`

**Interfaces:**
- Consumes: Tasks 1 and 2.
- Produces: `chain_setup_fields()` (character vector); `set_chain_setup(setup)`, `clear_chain_setup()`; `draw_chain_seeds(nchain)` (integer matrix 2 by nchain, rownames `r` and `cpp`); `save_session_rng()`, `restore_session_rng(saved)`; `seed_chain_rng(seed_r, seed_cpp)`; `run_one_chain(setup, seed_r, seed_cpp, chain)` returning `list(chain, seed_r, seed_cpp, elapsed, draws, latent_sums, waic)`; `run_chain_task(k, seeds)`; `chain_failure(k, message)`, `stop_if_chain_failed(res)`; `run_chain_batch(handle, batch, seeds, task = run_chain_task)` (serial when `handle` is `NULL`; Task 4 adds the cluster branch); `drop_execution(fit)`.
- `draws` is a named list with elements `beta_theta, theta0, z, psi, p, q, B0, B, G, A, C, Bs, Gs, As, Cs, U, L, sigmab, sigmabs, sigmah, idx_ls, varPart, tau`, each the per-chain array exactly as today's `*_output_chain` object, or `NULL` where the model has none. `latent_sums` is `list(z, psi, theta, w)` of matrices or `NULL`. `waic` is `list(jsdm, w, y)` of WAIC parts or `NULL`.

- [ ] **Step 1: Confirm the line anchors and the absence of other `on.exit()` calls**

```bash
grep -n "# run MCMC\|setOccJSDMSeed(sample\|for (chain in 1:nchain)\|# compute WAIC\|# reparametrise factor\|infos <- list\|on.exit" R/runOccJSDM.R
```

Expected: `# run MCMC` at 1040, the seed call at 1051, the loop at 1053, `# compute WAIC` at 1406, `# reparametrise factor` at 1426, `infos <- list(` at 1488, and no `on.exit` line. If an `on.exit()` exists, every `on.exit()` this task adds must use `add = TRUE` (they do) and the existing one must be changed to `add = TRUE` as well.

- [ ] **Step 2: Write the failing tests**

Append to `tests/testthat/test-parallel-chains.R`:

```r
test_that("chain seeds are a 2 by nchain integer matrix drawn from R's RNG", {
  set.seed(3); a <- occJSDM:::draw_chain_seeds(3L)
  set.seed(3); b <- occJSDM:::draw_chain_seeds(3L)
  expect_identical(a, b)
  expect_identical(dim(a), c(2L, 3L))
  expect_identical(rownames(a), c("r", "cpp"))
  expect_type(a, "integer")
  expect_false(anyDuplicated(as.vector(a)) > 0)
  set.seed(4); expect_false(identical(occJSDM:::draw_chain_seeds(3L), a))
})

test_that("a chain's seeds do not depend on how many chains follow it", {
  set.seed(8); one <- occJSDM:::draw_chain_seeds(1L)
  set.seed(8); two <- occJSDM:::draw_chain_seeds(2L)
  expect_identical(one[, 1], two[, 1])
})

test_that("restoring the session RNG puts the value and the kinds back", {
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)))
  RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
  set.seed(9); invisible(runif(2))
  saved <- occJSDM:::save_session_rng()
  expected <- runif(3)
  occJSDM:::seed_chain_rng(123L, 456L)       # pins Mersenne-Twister etc., as a chain does
  invisible(rnorm(5))
  occJSDM:::restore_session_rng(saved)
  expect_identical(RNGkind()[1:3], c("L'Ecuyer-CMRG", "Inversion", "Rejection"))
  expect_identical(runif(3), expected)
})

test_that("chain 1 of a two-chain fit equals the one-chain fit under the same seed", {
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  set.seed(515); one <- fit_fixture(sim, MCMCparams = list(nchain = 1, nburn = 20, niter = 20, nthin = 1))
  set.seed(515); two <- fit_fixture(sim, MCMCparams = list(nchain = 2, nburn = 20, niter = 20, nthin = 1))
  expect_identical(one$results_output$jsdm_output$B0_output[, , 1], two$results_output$jsdm_output$B0_output[, , 1])
  expect_identical(one$results_output$p_output[, , , 1], two$results_output$p_output[, , , 1])
  expect_identical(one$infos$chain_seeds[, 1], two$infos$chain_seeds[, 1])
  expect_false(identical(two$results_output$jsdm_output$B0_output[, , 1], two$results_output$jsdm_output$B0_output[, , 2]))
})

test_that("a fit leaves the session RNG just after its seed draws and clears the setup", {
  # The state after a fit is the state after the chain seeds were drawn, so a
  # third chain's seeds, had one been requested, are exactly the next draws.
  # Comparing two fits rather than replaying the draws by hand keeps the test
  # true even if the data preparation before the loop uses R's RNG.
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  set.seed(21); two <- fit_fixture(sim, MCMCparams = list(nchain = 2, nburn = 20, niter = 20, nthin = 1))
  next_pair <- sample.int(.Machine$integer.max, 2L)
  set.seed(21); three <- fit_fixture(sim, MCMCparams = list(nchain = 3, nburn = 20, niter = 20, nthin = 1))
  expect_identical(unname(three$infos$chain_seeds[, 3]), next_pair)
  expect_identical(two$infos$chain_seeds, three$infos$chain_seeds[, 1:2])
  expect_false(exists("setup", envir = occJSDM:::chain_env, inherits = FALSE))
})

test_that("the new infos fields record the sampling settings, seeds and execution", {
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  set.seed(6); fit <- fit_fixture(sim, MCMCparams = list(nchain = 2, nburn = 20, niter = 20, nthin = 1))
  expect_identical(fit$infos$MCMCparams, list(nchain = 2L, nburn = 20L, niter = 20L, nthin = 1L))
  expect_identical(dim(fit$infos$chain_seeds), c(2L, 2L))
  ex <- fit$infos$execution
  expect_identical(ex$ncores, 1L); expect_identical(ex$workers, 1L)
  expect_identical(ex$thread_setting, "ambient")
  expect_identical(ex$threads_per_worker, occJSDM:::ambient_rcpp_threads())
  expect_null(ex$progress_file); expect_null(ex$blas_threads_requested)
  expect_identical(ex$blas_declared, nzchar(Sys.getenv("OCCJSDM_REFERENCE_BLAS")))
  expect_length(ex$elapsed, 2L)
  expect_null(occJSDM:::drop_execution(fit)$infos$execution)
})

test_that("summarised latent means equal the mean over all retained draws", {
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  mc <- list(nchain = 2, nburn = 20, niter = 20, nthin = 1)
  set.seed(31); summarised <- suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
    listParams = list(n_factors = 2, n_supportpoints = FIXTURE_KNOTS),
    occCovariates = fixture_occ_covariates(), spatCovariates = fixture_spat_covariates(),
    MCMCparams = mc, summarisedLatentPresences = TRUE)))
  set.seed(31); full <- suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
    listParams = list(n_factors = 2, n_supportpoints = FIXTURE_KNOTS),
    occCovariates = fixture_occ_covariates(), spatCovariates = fixture_spat_covariates(),
    MCMCparams = mc, summarisedLatentPresences = FALSE)))
  z_mean <- apply(full$results_output$z_output, c(1, 2), mean)
  expect_equal(summarised$results_output$z_output, z_mean, tolerance = 1e-12)
  expect_identical(summarised$results_output$WAIC, full$results_output$WAIC)
  expect_identical(summarised$results_output$jsdm_output$B0_output, full$results_output$jsdm_output$B0_output)
})
```

- [ ] **Step 3: Run the tests to verify they fail**

```bash
Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-parallel-chains.R")'
```

Expected: the new tests fail ("object 'draw_chain_seeds' not found", `infos$MCMCparams` is `NULL`, and so on). The existing suite is untouched so far.

- [ ] **Step 4: Add the seed, RNG, setup and task helpers to `R/chains.R`**

```r
# The fields the chain body reads from runOccJSDM()'s frame. Computed on
# 13 October 2026 with codetools::findGlobals() on the loop body plus the
# three names it both reads and assigns (list_data) or reads inside a call
# (Tr, N3). A PSOCK worker enforces completeness: a missing field is an error.
chain_setup_fields <- function() c(
  "model", "jsdmModel", "threshold",
  "y", "y_NA", "M", "K", "P", "sumM", "sumP", "sumK", "maxP", "primerId_p",
  "idx_w_k", "idx_p_k", "idx_z_w", "list_idx",
  "n", "N", "N3", "S", "d", "g", "gt", "ps", "ncov_theta", "ncov_psi",
  "X_theta", "X_psi", "Tr",
  "list_Xs", "list_SoRSummaries",
  "list_data", "list_priors", "spatial_sd_prior",
  "a_theta0", "b_theta0", "a_p", "b_p", "a_q", "b_q", "b_betatheta", "B_betatheta",
  "nburn", "niter", "nthin", "summarisedLatentPresences")

set_chain_setup <- function(setup) {
  assign("setup", setup, envir = chain_env)
  invisible(NULL)
}

clear_chain_setup <- function() {
  if (exists("setup", envir = chain_env, inherits = FALSE)) rm("setup", envir = chain_env)
  invisible(NULL)
}

# Two seeds per chain, drawn in chain order so that chain k's seeds do not
# depend on how many chains follow it. Distinct by construction except for an
# astronomically unlikely collision, which is redrawn.
draw_chain_seeds <- function(nchain) {
  draw <- function() {
    seeds <- matrix(NA_integer_, nrow = 2L, ncol = nchain, dimnames = list(c("r", "cpp"), NULL))
    for (k in seq_len(nchain)) seeds[, k] <- sample.int(.Machine$integer.max, 2L)
    seeds
  }
  seeds <- draw()
  while (anyDuplicated(as.vector(seeds)) > 0L) seeds <- draw()
  seeds
}

save_session_rng <- function() {
  if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) invisible(sample.int(1L))
  get(".Random.seed", envir = globalenv(), inherits = FALSE)
}

# Reset first, then restore: set.seed() with the pinned kinds clears R's
# cached Box-Muller deviate the same way in every path, and the saved state
# carries the user's kinds back with it. Registered with on.exit() in
# runOccJSDM() so it runs on every exit. Spec: "After the chains".
restore_session_rng <- function(saved) {
  set.seed(0L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  assign(".Random.seed", saved, envir = globalenv())
  invisible(NULL)
}

seed_chain_rng <- function(seed_r, seed_cpp) {
  set.seed(seed_r, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  setOccJSDMSeed(seed_cpp)
  invisible(NULL)
}

chain_failure <- function(k, message) {
  structure(list(chain = k, message = message), class = "occJSDM_chain_failure")
}

stop_if_chain_failed <- function(res) {
  if (inherits(res, "occJSDM_chain_failure")) {
    stop(sprintf("Chain %d failed: %s", res$chain, res$message), call. = FALSE)
  }
  res
}

# The task a worker (or the in-process loop) runs for one chain. Package
# level, so it crosses the socket as a namespace reference; setup is read
# from chain_env, which the bootstrap filled.
run_chain_task <- function(k, seeds) {
  tryCatch(run_one_chain(get("setup", envir = chain_env, inherits = FALSE),
                         seeds["r", k], seeds["cpp", k], k),
           error = function(e) chain_failure(k, conditionMessage(e)))
}

# Task 4 adds the cluster branch; with handle = NULL this is the in-process path.
run_chain_batch <- function(handle, batch, seeds, task = run_chain_task) {
  if (is.null(handle)) return(lapply(batch, task, seeds = seeds))
  parallel::clusterApplyLB(handle$cl, batch, task, seeds = seeds)
}

drop_execution <- function(fit) {
  fit$infos$execution <- NULL
  fit
}
```

- [ ] **Step 5: Write `run_one_chain()` by moving the loop body**

Add to `R/chains.R` the function below. The part marked "moved" is `R/runOccJSDM.R` lines 1055 to 1365 inclusive (from the `# chain output` comment to the closing brace of the `for (iter ...)` loop), pasted verbatim, then edited exactly as listed in Step 6. Do not retype it; cut it from `R/runOccJSDM.R` and paste it here, so that Alex's diff shows a move.

```r
# One chain, from its starting values to its last retained draw. Everything
# the body reads arrives in setup and is unpacked into this frame, so the
# moved code runs unchanged. The chain's randomness is fixed by its two seeds
# alone: it does not continue any other chain's stream, in process or not.
run_one_chain <- function(setup, seed_r, seed_cpp, chain) {
  started <- proc.time()[["elapsed"]]
  list2env(setup, envir = environment())
  seed_chain_rng(seed_r, seed_cpp)

  # ---- moved from runOccJSDM(): chain output, starting values, iteration loop ----
  # (lines 1055 to 1365 of R/runOccJSDM.R before this change, edited per the plan)
  # ---- end of moved code ----

  if (!model %in% c("occupancy", "two_stage")) {
    beta_theta_output_chain <- NULL; theta0_output_chain <- NULL
    z_output_chain <- NULL; psi_output_chain <- NULL
  }
  if (model != "two_stage") { p_output_chain <- NULL; q_output_chain <- NULL }
  if (model %in% c("occupancy", "two_stage") && summarisedLatentPresences) {
    z_output_chain <- NULL; psi_output_chain <- NULL
  }

  list(
    chain = chain, seed_r = seed_r, seed_cpp = seed_cpp,
    elapsed = proc.time()[["elapsed"]] - started,
    draws = list(
      beta_theta = beta_theta_output_chain, theta0 = theta0_output_chain,
      z = z_output_chain, psi = psi_output_chain,
      p = p_output_chain, q = q_output_chain,
      B0 = B0_output_chain, B = B_output_chain, G = G_output_chain, A = A_output_chain,
      C = C_output_chain, Bs = Bs_output_chain, Gs = Gs_output_chain, As = As_output_chain,
      Cs = Cs_output_chain, U = U_output_chain, L = L_output_chain,
      sigmab = sigmab_output_chain, sigmabs = sigmabs_output_chain, sigmah = sigmah_output_chain,
      idx_ls = idx_ls_output_chain, varPart = varPart_output_chain, tau = tau_output_chain),
    latent_sums = list(z = z_sum, psi = psi_sum, theta = theta_sum, w = w_sum),
    waic = list(jsdm = c(list_waic_jsdm, list(n = niter)),
                w = if (is.null(list_waic_w)) NULL else c(list_waic_w, list(n = niter)),
                y = if (is.null(list_waic_y)) NULL else c(list_waic_y, list(n = niter))))
}
```

- [ ] **Step 6: Apply the five edits inside the moved body**

Edit 1, the dead allocation. In the `# chain output` block, delete the single line

```r
          theta_output_chain <- array(NA, dim = c(N, S, niter))
```

Edit 2, per-chain sums and WAIC partials. Immediately after the `# chain output` block's closing brace (before `# starting values`), insert:

```r
    # Per-chain sums of the retained latent states and per-chain WAIC
    # accumulators. The parent combines them in chain order.
    z_sum <- psi_sum <- theta_sum <- w_sum <- NULL
    if (model %in% c("occupancy", "two_stage")) {
      if (summarisedLatentPresences) {
        z_sum <- matrix(0, n, S)
        psi_sum <- matrix(0, n, S)
      }
      theta_sum <- matrix(0, N, S)
    }
    if (model == "two_stage") w_sum <- matrix(0, N, S)
    list_waic_jsdm <- create_waic_quantities(n * S)
    list_waic_w <- if (model %in% c("occupancy", "two_stage")) create_waic_quantities(N * S) else NULL
    list_waic_y <- if (model == "two_stage") create_waic_quantities(N3 * S) else NULL
    currentWAICiter <- 1
```

Edit 3, the latent-state accumulation. In the retained-iteration block, replace

```r
            if(summarisedLatentPresences){
              z_output_mean <- z_output_mean +
                (1 / (niter * nchain)) * z
              psi_output_mean <- psi_output_mean +
                (1 / (niter * nchain)) * psi
            } else {
              z_output_chain[,,currentIter] <- z
              psi_output_chain[,,currentIter] <- psi
            }
            theta_output_mean <- theta_output_mean +
              (1 / (niter * nchain)) * theta
```

with

```r
            if(summarisedLatentPresences){
              z_sum <- z_sum + z
              psi_sum <- psi_sum + psi
            } else {
              z_output_chain[,,currentIter] <- z
              psi_output_chain[,,currentIter] <- psi
            }
            theta_sum <- theta_sum + theta
```

Edit 4, the collection-state accumulation. Replace

```r
            w_output_mean <- w_output_mean +
              (1 / (niter * nchain)) * w
```

with

```r
            w_sum <- w_sum + w
```

Edit 5, nothing else. The WAIC update block (`# update WIAC`) stays exactly as it is: it now updates the per-chain accumulators created in Edit 2 with `currentWAICiter` running from 1 within the chain. The lines that copied `*_output_chain` into the full arrays (old lines 1367 to 1402) are not moved into the function; they are rewritten in the parent in Step 7.

Check the move with a diff of the old file against the new function: every line of the old body is accounted for, and no sampler call changed.

- [ ] **Step 7: Rewrite the parent side of `runOccJSDM()`**

Replace everything from the `# run MCMC` comment (old line 1040) through the closing brace of the `# compute WAIC` block (old line 1424) with:

```r
  # run MCMC
  message("Running MCMC")

  # Every chain gets its own R seed and C++ seed, drawn here from R's RNG, so
  # set.seed() before the fit reproduces it and consecutive fits differ. A
  # chain's draws then depend on its two seeds alone, not on which chains ran
  # before it or in which process, which is what makes a fit identical whatever
  # MCMCparams$ncores says. The session's RNG is reset and restored on every
  # exit. Spec: dev/superpowers/specs/2026-10-13-parallel-chains-design.md.
  seeds <- draw_chain_seeds(nchain)
  rng_saved <- save_session_rng()
  on.exit(restore_session_rng(rng_saved), add = TRUE)

  setup <- mget(chain_setup_fields(), envir = environment(), ifnotfound = list(NULL))
  set_chain_setup(setup)
  on.exit(clear_chain_setup(), add = TRUE)

  budget <- plan_chain_budget(nchain, ncores)
  handle <- NULL
  progress_file <- NULL

  latent_sums <- NULL
  waic_parts <- vector("list", nchain)
  elapsed <- numeric(nchain)
  batches <- split(seq_len(nchain), ceiling(seq_len(nchain) / budget$workers))
  for (batch in batches) {
    results <- run_chain_batch(handle, batch, seeds)
    for (res in results) {
      res <- stop_if_chain_failed(res)
      chain <- res$chain
      dr <- res$draws

      if(model %in% c("occupancy","two_stage")){
        beta_theta_output[,,,chain] <- dr$beta_theta
        theta0_output[,,chain] <- dr$theta0
        if(!summarisedLatentPresences){
          z_output[,,,chain] <- dr$z
          psi_output[,,,chain] <- dr$psi
        }
      }

      if(model == "two_stage"){
        p_output[,,,chain] <- dr$p
        q_output[,,,chain] <- dr$q
      }

      # save jsdm params
      {
        B0_output[,,chain] <- dr$B0
        B_output[,,,chain] <- dr$B
        G_output[,,,chain] <- dr$G
        A_output[,,,chain] <- dr$A
        C_output[,,,chain] <- dr$C
        Bs_output[,,,chain] <- dr$Bs
        Gs_output[,,,chain] <- dr$Gs
        As_output[,,,chain] <- dr$As
        Cs_output[,,,chain] <- dr$Cs
        U_output[,,,chain] <- dr$U
        L_output[,,,chain] <- dr$L
        sigmab_output[,chain] <- dr$sigmab
        sigmabs_output[,chain] <- dr$sigmabs
        sigmah_output[,chain] <- dr$sigmah
        idx_ls_output[,chain] <- dr$idx_ls
        varPart_output[,,,chain] <- dr$varPart
        tau_output[,,chain] <- dr$tau
      }

      latent_sums <- add_latent_sums(latent_sums, res$latent_sums)
      waic_parts[[chain]] <- res$waic
      elapsed[chain] <- res$elapsed
      rm(dr)
    }
    rm(results, res)
  }
  for (chain in seq_len(nchain)) message(sprintf("Chain %d finished in %.1f s", chain, elapsed[chain]))

  # latent-state means, from the per-chain sums
  if (model %in% c("occupancy","two_stage")) {
    if (summarisedLatentPresences) {
      z_output_mean <- latent_sums$z / (niter * nchain)
      psi_output_mean <- latent_sums$psi / (niter * nchain)
    }
    theta_output_mean <- latent_sums$theta / (niter * nchain)
  }
  if (model == "two_stage") w_output_mean <- latent_sums$w / (niter * nchain)

  # compute WAIC
  {
    numIters <- niter * nchain
    for (part in waic_parts) if (part$jsdm$n != niter) stop("Current wAIC iter wrong")
    WAIC_jsdm <- compute_waic(merge_waic_parts(lapply(waic_parts, `[[`, "jsdm")), numIters)

    if(model %in% c("occupancy","two_stage")){
      WAIC_w <- compute_waic(merge_waic_parts(lapply(waic_parts, `[[`, "w")), numIters)
    } else {
      WAIC_w <- 0
    }
    if(model %in% "two_stage"){
      WAIC_y <- compute_waic(merge_waic_parts(lapply(waic_parts, `[[`, "y")), numIters)
    } else {
      WAIC_y <- 0
    }

    WAIC <- WAIC_jsdm + WAIC_w + WAIC_y
  }
```

Then, in the `# chain parameters` block (old lines 901 to 907), add one line after `nthin <- ...`:

```r
    ncores <- resolve_ncores(MCMCparams)
```

Then, after the `infos$intercept_prior <- intercept_prior` line near the end, add:

```r
  infos$MCMCparams <- list(nchain = as.integer(nchain), nburn = as.integer(nburn),
                           niter = as.integer(niter), nthin = as.integer(nthin))
  infos$chain_seeds <- seeds
  infos$execution <- list(
    ncores = budget$ncores,
    workers = budget$workers,
    threads_per_worker = if (budget$thread_setting == "ambient") ambient_rcpp_threads() else budget$threads_per_worker,
    thread_setting = budget$thread_setting,
    blas_threads_requested = if (is.null(handle)) NULL else handle$blas_threads_requested,
    blas_declared = nzchar(Sys.getenv("OCCJSDM_REFERENCE_BLAS")),
    blas = unname(extSoftVersion()[["BLAS"]]),
    lapack = La_library(),
    elapsed = elapsed,
    progress_file = progress_file)
```

- [ ] **Step 8: Declare the unpacked names for R CMD check**

`run_one_chain()` reads the setup fields from its own frame after `list2env()`, which `R CMD check`'s code analysis cannot see. Append the fields to the vector in `R/globals.R`, with a comment, keeping the existing entries:

```r
# Fields of the setup list that run_one_chain() in R/chains.R unpacks into its
# frame with list2env(); the moved chain body reads them as plain names.
# Keep in step with chain_setup_fields().
utils::globalVariables(c(
  "model", "jsdmModel", "threshold",
  "y", "y_NA", "M", "K", "P", "sumM", "sumP", "sumK", "maxP", "primerId_p",
  "idx_w_k", "idx_p_k", "idx_z_w", "list_idx",
  "n", "N", "N3", "S", "d", "g", "gt", "ps", "ncov_theta", "ncov_psi",
  "X_theta", "X_psi", "Tr",
  "list_Xs", "list_SoRSummaries",
  "list_data", "list_priors", "spatial_sd_prior",
  "a_theta0", "b_theta0", "a_p", "b_p", "a_q", "b_q", "b_betatheta", "B_betatheta",
  "nburn", "niter", "nthin", "summarisedLatentPresences"
))
```

(`y` is already declared; a duplicate is harmless.)

- [ ] **Step 9: Run the new tests, then the whole suite**

```bash
Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-parallel-chains.R")'
Rscript -e 'devtools::test()' 2>&1 | tail -25
```

Expected: everything passes. The seed tests in `test-regression-bugs.R` and `test-rng-safety.R` still pass because they assert relations, not values. If a test fails with "object 'X' not found" inside `run_one_chain()`, `X` is a missing setup field: add it to `chain_setup_fields()` and `R/globals.R`, and re-run.

- [ ] **Step 10: Commit**

```bash
git add R/chains.R R/runOccJSDM.R R/globals.R tests/testthat/test-parallel-chains.R
git commit -m "Move the chain loop body into run_one_chain() with per-chain seeds and partial accumulators

Every chain now sets its own R and C++ seeds, returns its draws, latent-state sums and WAIC partials, and the parent combines them in chain order. Serial fits only, so far; the session RNG is reset and restored on every exit.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 4: PSOCK workers: bootstrap, build identity, dispatch, shutdown

**Files:**
- Modify: `R/chains.R`
- Modify: `tests/testthat/test-parallel-chains.R`

**Interfaces:**
- Consumes: Task 3's `run_chain_batch()`, `chain_failure()`, `stop_if_chain_failed()`, `set_chain_setup()`, `save_session_rng()`, `restore_session_rng()`.
- Produces: `blas_thread_variables()` (named character); `with_worker_environment(expr)`; `is_dev_loaded()`; `bootstrap_chain_worker(libpaths, threads, dev_path)` (runs in the worker, environment set to `globalenv()` before sending); `start_chain_cluster(workers, threads, progress_file)` returning a handle `list(cl, pids, worker_threads, progress_file, blas_threads_requested, threads)`; `send_chain_setup(handle, setup)`; `package_fingerprint()` returning `list(path, version, dll_md5, code_md5)`; `compare_fingerprints(parent, worker, index)`; `check_worker_builds(handle)`; `stop_chain_cluster(handle)`; `pin_session_threads(threads)` returning a restore function; `run_chain_batch()` gains its cluster branch (already written in Task 3).

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-parallel-chains.R`:

```r
test_that("the R-code fingerprint matches across fresh workers and changes when a body changes", {
  skip_on_cran()
  handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile("occJSDM-progress-", fileext = ".log"))
  on.exit(occJSDM:::stop_chain_cluster(handle), add = TRUE)
  parent <- occJSDM:::package_fingerprint()
  workers <- parallel::clusterCall(handle$cl, occJSDM:::package_fingerprint)
  for (w in workers) expect_identical(w, parent)
  expect_true(occJSDM:::check_worker_builds(handle))
  expect_identical(handle$worker_threads, c("1", "1"))
  expect_identical(length(handle$pids), 2L)
  expect_true(all(alive(handle$pids)))
  # negative case on a toy namespace-like environment
  fp <- function(env) {
    nms <- sort(ls(env, all.names = TRUE))
    text <- unlist(lapply(nms, function(n) c(n, deparse(get(n, envir = env)))))
    tmp <- tempfile(); on.exit(unlink(tmp)); writeLines(text, tmp); unname(tools::md5sum(tmp))
  }
  e1 <- new.env(); e1$f <- function(x) x + 1
  e2 <- new.env(); e2$f <- function(x) x + 2
  expect_false(identical(fp(e1), fp(e2)))
})

test_that("a build mismatch is reported with both paths and the differing field", {
  parent <- list(path = "/a/occJSDM", version = "0.1.1", dll_md5 = "x", code_md5 = "y")
  expect_true(occJSDM:::compare_fingerprints(parent, parent, 1L))
  worker <- parent; worker$code_md5 <- "z"; worker$path <- "/b/occJSDM"
  expect_error(occJSDM:::compare_fingerprints(parent, worker, 2L), "Worker 2 .*path.*/a/occJSDM.*/b/occJSDM")
  worker <- parent; worker$dll_md5 <- "q"
  expect_error(occJSDM:::compare_fingerprints(parent, worker, 1L), "dll_md5")
})

test_that("a chain that fails in a worker is reported by number and no worker survives", {
  handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile("occJSDM-progress-", fileext = ".log"))
  pids <- handle$pids
  occJSDM:::clear_chain_setup()                      # the real task will fail: no setup in the workers
  res <- occJSDM:::run_chain_batch(handle, 1:2, matrix(1:4, 2, dimnames = list(c("r", "cpp"), NULL)))
  expect_s3_class(res[[2]], "occJSDM_chain_failure")
  expect_error(occJSDM:::stop_if_chain_failed(res[[2]]), "^Chain 2 failed: ")
  occJSDM:::stop_chain_cluster(handle)
  Sys.sleep(0.5)
  expect_false(any(alive(pids)))
})

test_that("cleanup kills every worker even when stopCluster() itself errors", {
  skip_on_cran()
  handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile("occJSDM-progress-", fileext = ".log"))
  pids <- handle$pids
  close(handle$cl[[1]]$con)                          # a broken connection makes stopCluster() throw
  expect_true(alive(pids[2]))
  expect_silent(occJSDM:::stop_chain_cluster(handle))
  Sys.sleep(0.5)
  expect_false(any(alive(pids)))
})

test_that("abandoning running workers leaves no worker process alive", {
  skip_on_cran(); skip_on_os("windows")
  handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile("occJSDM-progress-", fileext = ".log"))
  on.exit(occJSDM:::stop_chain_cluster(handle), add = TRUE)
  pids <- handle$pids
  system(sprintf("(sleep 1; kill -INT %d) > /dev/null 2>&1 &", Sys.getpid()))   # timer starts after the cluster is up
  t0 <- Sys.time()
  outcome <- tryCatch({
    occJSDM:::run_chain_batch(handle, 1:2, matrix(1:4, 2), task = function(k, seeds) { Sys.sleep(10); k })
    "completed"
  }, interrupt = function(c) "interrupted")
  waited <- as.numeric(Sys.time() - t0, units = "secs")
  expect_identical(outcome, "interrupted")
  expect_lt(waited, 8)                               # the workers were still busy when cleanup begins
  expect_true(all(alive(pids)))
  occJSDM:::stop_chain_cluster(handle)
  Sys.sleep(0.5)
  expect_false(any(alive(pids)))
})

test_that("the BLAS thread variables are set for the workers only", {
  skip_on_cran()
  vars <- names(occJSDM:::blas_thread_variables())
  before <- Sys.getenv(vars, unset = NA_character_, names = TRUE)
  handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile("occJSDM-progress-", fileext = ".log"))
  on.exit(occJSDM:::stop_chain_cluster(handle), add = TRUE)
  expect_identical(Sys.getenv(vars, unset = NA_character_, names = TRUE), before)
  seen <- parallel::clusterCall(handle$cl, function(v) Sys.getenv(v, names = TRUE), vars)
  for (s in seen) expect_identical(unname(s), unname(occJSDM:::blas_thread_variables()))
  expect_identical(handle$blas_threads_requested, occJSDM:::blas_thread_variables())
})

test_that("the session RNG is the same after an error exit or an interrupt in either driver", {
  skip_on_cran(); skip_on_os("windows")
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  run_guarded <- function(driver) {
    saved <- occJSDM:::save_session_rng()
    on.exit(occJSDM:::restore_session_rng(saved))
    driver()
  }
  seeds <- matrix(1:4, 2, dimnames = list(c("r", "cpp"), NULL))
  failing <- function(k, seeds) occJSDM:::chain_failure(k, "boom")
  serial_error <- function() for (r in occJSDM:::run_chain_batch(NULL, 1:2, seeds, task = failing)) occJSDM:::stop_if_chain_failed(r)
  parallel_error <- function() {
    handle <- occJSDM:::start_chain_cluster(2L, 1L, tempfile(fileext = ".log"))
    on.exit(occJSDM:::stop_chain_cluster(handle))
    for (r in occJSDM:::run_chain_batch(handle, 1:2, seeds, task = failing)) occJSDM:::stop_if_chain_failed(r)
  }
  next_draws <- function() list(rnorm(3), runif(3), sample(100, 3))
  for (normal_kind in c("Inversion", "Box-Muller")) {
    RNGkind("Mersenne-Twister", normal_kind, "Rejection")
    set.seed(77); invisible(rnorm(1))                 # under Box-Muller this leaves a cached deviate pending
    expect_error(run_guarded(serial_error), "Chain 1 failed: boom"); after_serial <- next_draws()
    set.seed(77); invisible(rnorm(1))
    expect_error(run_guarded(parallel_error), "Chain 1 failed: boom"); after_parallel <- next_draws()
    expect_identical(after_serial, after_parallel)
  }
})
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-parallel-chains.R")'
```

Expected: the new tests fail with "object 'start_chain_cluster' not found" and similar.

- [ ] **Step 3: Add the cluster helpers to `R/chains.R`**

```r
# Workers are asked to run their BLAS single-threaded whatever ncores is; a
# BLAS thread count that varied with ncores would make results vary with it.
# These are requests the package cannot verify (spec, "Interface").
blas_thread_variables <- function() c(
  OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
  MKL_DOMAIN_NUM_THREADS = "MKL_DOMAIN_ALL=1", VECLIB_MAXIMUM_THREADS = "1")

# Set the thread variables in this process only while expr runs (the workers
# read them at process start), then put the previous values back.
with_worker_environment <- function(expr) {
  vars <- blas_thread_variables()
  old <- Sys.getenv(names(vars), unset = NA_character_, names = TRUE)
  on.exit({
    for (nm in names(old)) {
      if (is.na(old[[nm]])) Sys.unsetenv(nm) else do.call(Sys.setenv, as.list(stats::setNames(old[[nm]], nm)))
    }
  }, add = TRUE)
  do.call(Sys.setenv, as.list(vars))
  force(expr)
}

is_dev_loaded <- function() {
  exists(".__DEVTOOLS__", envir = asNamespace("occJSDM"), inherits = FALSE)
}

# Runs inside a fresh worker before occJSDM exists there, so it is sent with
# globalenv() as its environment and uses only base, RcppParallel and pkgload.
bootstrap_chain_worker <- function(libpaths, threads, dev_path) {
  .libPaths(libpaths)
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = as.character(threads))
  RcppParallel::setThreadOptions(numThreads = threads)
  if (is.null(dev_path)) {
    loadNamespace("occJSDM")
  } else {
    # compile = FALSE: never compile in a worker; load the shared object the
    # parent built. Other arguments stay at pkgload's defaults so that the
    # worker's namespace holds the same functions as a parent loaded with
    # devtools::load_all(); the code fingerprint compares exactly that.
    # (pkgload 1.5.3 sources test helpers into the attached package
    # environment, not the namespace, so the helpers argument cannot
    # change the fingerprint; checked 13 October 2026.)
    pkgload::load_all(dev_path, compile = FALSE, quiet = TRUE, attach = FALSE)
  }
  list(pid = Sys.getpid(), threads = Sys.getenv("RCPP_PARALLEL_NUM_THREADS"))
}

start_chain_cluster <- function(workers, threads, progress_file) {
  cl <- with_worker_environment(parallel::makeCluster(workers, outfile = progress_file))
  handle <- list(cl = cl, pids = integer(0), worker_threads = character(0),
                 progress_file = progress_file, blas_threads_requested = blas_thread_variables(),
                 threads = as.integer(threads))
  boot <- bootstrap_chain_worker
  environment(boot) <- globalenv()
  dev_path <- if (is_dev_loaded()) getNamespaceInfo("occJSDM", "path") else NULL
  info <- tryCatch(parallel::clusterCall(cl, boot, .libPaths(), threads, dev_path),
                   error = function(e) { stop_chain_cluster(handle); stop(e) })
  handle$pids <- vapply(info, function(x) as.integer(x$pid), integer(1))
  handle$worker_threads <- vapply(info, function(x) x$threads, character(1))
  handle
}

send_chain_setup <- function(handle, setup) {
  parallel::clusterCall(handle$cl, set_chain_setup, setup)   # a namespace function, so no frame is serialised
  invisible(NULL)
}

# Four things identify a build: namespace path, version, the loaded shared
# object's md5, and a digest of every function's deparsed code. Functions
# only: native-routine bindings carry process-specific pointers (spec,
# "Build identity check").
package_fingerprint <- function() {
  ns <- asNamespace("occJSDM")
  nms <- sort(ls(ns, all.names = TRUE))
  funs <- nms[vapply(nms, function(n) is.function(get(n, envir = ns)), logical(1))]
  text <- unlist(lapply(funs, function(n) c(n, deparse(get(n, envir = ns)))))
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines(text, tmp)
  dll <- getLoadedDLLs()[["occJSDM"]]
  list(path = normalizePath(getNamespaceInfo(ns, "path")),
       version = as.character(getNamespaceVersion(ns)),
       dll_md5 = unname(tools::md5sum(dll[["path"]])),
       code_md5 = unname(tools::md5sum(tmp)))
}

compare_fingerprints <- function(parent, worker, index) {
  for (field in c("path", "version", "dll_md5", "code_md5")) {
    if (!identical(parent[[field]], worker[[field]])) {
      stop(sprintf(paste0("Worker %d loaded a different occJSDM build (%s differs): parent %s, worker %s. ",
                          "Reload or reinstall so that the parent and the workers run the same code."),
                   index, field, parent$path, worker$path), call. = FALSE)
    }
  }
  invisible(TRUE)
}

check_worker_builds <- function(handle) {
  parent <- package_fingerprint()
  workers <- parallel::clusterCall(handle$cl, package_fingerprint)
  for (i in seq_along(workers)) compare_fingerprints(parent, workers[[i]], i)
  invisible(TRUE)
}

# stopCluster() throws on a broken connection; the try() keeps the PID step
# reachable. pskill() ends a worker that is mid-chain, which would otherwise
# only read the stop message when its chain finished.
stop_chain_cluster <- function(handle) {
  if (is.null(handle)) return(invisible(NULL))
  try(parallel::stopCluster(handle$cl), silent = TRUE)
  if (length(handle$pids)) tools::pskill(handle$pids)
  invisible(NULL)
}

# ncores > 1 with one chain runs in process with the session's RcppParallel
# count set to the budget for the fit; the returned function puts it back.
pin_session_threads <- function(threads) {
  old <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
  Sys.setenv(RCPP_PARALLEL_NUM_THREADS = as.character(threads))
  RcppParallel::setThreadOptions(numThreads = threads)
  function() {
    if (is.na(old)) {
      Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS")
      RcppParallel::setThreadOptions(numThreads = "auto")
    } else {
      Sys.setenv(RCPP_PARALLEL_NUM_THREADS = old)
      RcppParallel::setThreadOptions(numThreads = as.integer(old))
    }
    invisible(NULL)
  }
}
```

- [ ] **Step 4: Run the tests to verify they pass**

```bash
Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-parallel-chains.R")'
```

Expected: all pass. This run exercises the pkgload branch of the bootstrap (the parent is dev-loaded). If the fingerprint test fails on `path`, check that `getNamespaceInfo("occJSDM", "path")` in the parent is the worktree path and that the worker's `load_all()` used the same path; `normalizePath()` is applied on both sides. If it fails on `code_md5`, the parent's and the worker's namespaces hold different function sets: print `setdiff()` of `ls(asNamespace("occJSDM"))` between the two (through `clusterCall()`), and align the worker's `load_all()` arguments with how the parent was loaded (for example the `helpers` argument, if test helpers turn out to be sourced into the namespace).

- [ ] **Step 5: Commit**

```bash
git add R/chains.R tests/testthat/test-parallel-chains.R
git commit -m "Add PSOCK worker bootstrap, build identity check, batch dispatch and shutdown for chains

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 5: Wire the parallel path into `runOccJSDM()`

**Files:**
- Modify: `R/runOccJSDM.R` (the parent block written in Task 3)
- Modify: `tests/testthat/test-parallel-chains.R`

**Interfaces:**
- Consumes: Task 4's helpers.
- Produces: the user-facing behaviour: `MCMCparams$ncores`, the start line with the `tail -f` command, the progress file, the execution record's cluster fields, the in-process thread pin for `nchain = 1`.

- [ ] **Step 1: Write the failing tests**

Append to `tests/testthat/test-parallel-chains.R`:

```r
skip_unless_declared_blas <- function() {
  if (!nzchar(Sys.getenv("OCCJSDM_REFERENCE_BLAS"))) {
    skip(paste0("fit identity across worker counts needs a provenance-checked single-threaded BLAS ",
                "(set OCCJSDM_REFERENCE_BLAS after checking); BLAS in use, as an indication only: ",
                extSoftVersion()[["BLAS"]]))
  }
}
same_fit <- function(a, b) expect_identical(occJSDM:::drop_execution(a), occJSDM:::drop_execution(b))
next_draws <- function() list(rnorm(3), runif(3), sample(100, 3))
fit_at <- function(sim, ncores, nchain = 2, seed = 2024, spatial = TRUE, ...) {
  set.seed(seed)
  fit <- fit_fixture(sim, spatial = spatial,
                     MCMCparams = list(nchain = nchain, nburn = 20, niter = 20, nthin = 1, ncores = ncores), ...)
  list(fit = fit, seed_state = get(".Random.seed", envir = globalenv()), draws = next_draws())
}

test_that("a two-stage spatial fit is identical at ncores = 1 and ncores = 2", {
  skip_unless_declared_blas()
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  serial <- fit_at(sim, 1L)
  parallel <- fit_at(sim, 2L)
  same_fit(serial$fit, parallel$fit)
  expect_identical(serial$seed_state, parallel$seed_state)
  expect_identical(serial$draws, parallel$draws)
  ex <- parallel$fit$infos$execution
  expect_identical(ex[c("ncores", "workers", "threads_per_worker", "thread_setting")],
                   list(ncores = 2L, workers = 2L, threads_per_worker = 1L, thread_setting = "budget"))
  expect_true(file.exists(ex$progress_file))
  expect_identical(ex$blas_threads_requested, occJSDM:::blas_thread_variables())
  expect_true(ex$blas_declared)
  expect_identical(serial$fit$infos$execution$thread_setting, "ambient")
  expect_null(serial$fit$infos$execution$progress_file)
})

test_that("a binary fit is identical at ncores = 1 and ncores = 2", {
  skip_unless_declared_blas()
  sim <- simulate_fixture(model = "binary", M = 1L, P = 1L, K = 1L)
  same_fit(fit_at(sim, 1L, spatial = FALSE)$fit, fit_at(sim, 2L, spatial = FALSE)$fit)
})

test_that("the same seed reproduces a parallel fit and a different seed changes it", {
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  a <- fit_at(sim, 2L, seed = 4242); b <- fit_at(sim, 2L, seed = 4242); d <- fit_at(sim, 2L, seed = 9999)
  same_fit(a$fit, b$fit)
  expect_false(identical(a$fit$results_output$jsdm_output$B0_output, d$fit$results_output$jsdm_output$B0_output))
  expect_false(exists("setup", envir = occJSDM:::chain_env, inherits = FALSE))
})

test_that("a parallel fit announces the progress file and the workers write to it", {
  sim <- simulate_fixture(model = "two_stage", P = 1L)
  set.seed(5)
  msgs <- capture_messages(fit <- suppressWarnings(runOccJSDM(sim$data_list,
    listParams = list(n_factors = 2), occCovariates = fixture_occ_covariates(),
    MCMCparams = list(nchain = 2, nburn = 100, niter = 100, nthin = 1, ncores = 2))))
  expect_true(any(grepl("Running 2 chains on 2 workers, 1 RcppParallel thread each. Follow progress with: tail -f ", msgs, fixed = TRUE)))
  expect_true(any(grepl("^Chain 2 finished in [0-9.]+ s", msgs)))
  pf <- fit$infos$execution$progress_file
  expect_true(file.exists(pf))
  lines <- readLines(pf)
  expect_true(any(grepl("Chain 1 - Burn in Iteration 100", lines, fixed = TRUE)))
  expect_true(any(grepl("Chain 2 - Iteration 100", lines, fixed = TRUE)))
})

test_that("one chain with ncores > 1 stays in process with the thread budget applied and restored", {
  sim <- simulate_fixture(model = "binary", M = 1L, P = 1L, K = 1L)
  before <- Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_)
  set.seed(12)
  msgs <- capture_messages(fit <- fit_fixture(sim, spatial = FALSE,
    MCMCparams = list(nchain = 1, nburn = 20, niter = 20, nthin = 1, ncores = 2)))
  expect_false(any(grepl("Follow progress", msgs)))
  ex <- fit$infos$execution
  expect_identical(ex[c("ncores", "workers", "threads_per_worker", "thread_setting")],
                   list(ncores = 2L, workers = 1L, threads_per_worker = 2L, thread_setting = "budget"))
  expect_null(ex$progress_file); expect_null(ex$blas_threads_requested)
  expect_identical(Sys.getenv("RCPP_PARALLEL_NUM_THREADS", unset = NA_character_), before)
})

test_that("non-default RNG kinds give the same fit at every ncores and come back afterwards", {
  skip_unless_declared_blas(); skip_on_cran()
  old_kind <- RNGkind()
  on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
  RNGkind("L'Ecuyer-CMRG", "Box-Muller", "Rejection")
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  set.seed(99); invisible(rnorm(1))                   # a cached Box-Muller deviate is pending
  serial <- fit_at(sim, 1L, seed = 99)
  expect_identical(RNGkind()[1:3], c("L'Ecuyer-CMRG", "Box-Muller", "Rejection"))
  parallel <- fit_at(sim, 2L, seed = 99)
  same_fit(serial$fit, parallel$fit)
  expect_identical(serial$draws, parallel$draws)
})

test_that("three chains on two workers and unsummarised latent states are identical to the in-process fit", {
  skip_unless_declared_blas(); skip_on_cran()
  sim <- simulate_fixture(model = "two_stage", P = 2L)
  same_fit(fit_at(sim, 1L, nchain = 3)$fit, fit_at(sim, 2L, nchain = 3)$fit)
  run <- function(ncores) {
    set.seed(8)
    suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
      listParams = list(n_factors = 2, n_supportpoints = FIXTURE_KNOTS),
      occCovariates = fixture_occ_covariates(), spatCovariates = fixture_spat_covariates(),
      MCMCparams = list(nchain = 2, nburn = 20, niter = 20, nthin = 1, ncores = ncores),
      summarisedLatentPresences = FALSE)))
  }
  a <- run(1L); b <- run(2L)
  same_fit(a, b)
  expect_identical(dim(b$results_output$z_output), c(FIXTURE_N, FIXTURE_S, 20L, 2L))
  for (model in c("continuous", "occupancy")) {
    s <- simulate_fixture(model = model, M = if (model == "occupancy") 2L else 1L, P = 1L, K = 1L)
    same_fit(fit_at(s, 1L, spatial = FALSE, seed = 3)$fit, fit_at(s, 2L, spatial = FALSE, seed = 3)$fit)
  }
})
```

- [ ] **Step 2: Run the tests to verify they fail**

```bash
OCCJSDM_REFERENCE_BLAS=1 Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-parallel-chains.R")'
```

Expected: the identity tests fail because `ncores = 2` still runs in process (`workers` is 2 but `handle` is `NULL`, so the execution record and the progress-file assertions fail). Setting `OCCJSDM_REFERENCE_BLAS=1` here is legitimate on this machine: its BLAS provenance was checked on 13 October 2026 (spec, "Tests"); Task 7's script re-checks it.

- [ ] **Step 3: Add the cluster branch to the parent block in `R/runOccJSDM.R`**

Replace the two lines

```r
  handle <- NULL
  progress_file <- NULL
```

with

```r
  handle <- NULL
  progress_file <- NULL
  if (budget$workers > 1L) {
    progress_file <- tempfile("occJSDM-progress-", fileext = ".log")
    handle <- start_chain_cluster(budget$workers, budget$threads_per_worker, progress_file)
    on.exit(stop_chain_cluster(handle), add = TRUE)
    check_worker_builds(handle)
    send_chain_setup(handle, setup)
    message(sprintf("Running %d chains on %d workers, %d RcppParallel thread%s each. Follow progress with: tail -f %s",
                    nchain, budget$workers, budget$threads_per_worker,
                    if (budget$threads_per_worker == 1L) "" else "s", progress_file))
  } else if (budget$ncores > 1L) {
    restore_threads <- pin_session_threads(budget$ncores)
    on.exit(restore_threads(), add = TRUE)
  }
```

- [ ] **Step 4: Run the tests, then the whole suite, then time the always-on tier**

```bash
OCCJSDM_REFERENCE_BLAS=1 Rscript -e 'devtools::load_all(quiet = TRUE); testthat::test_file("tests/testthat/test-parallel-chains.R")'
OCCJSDM_REFERENCE_BLAS=1 Rscript -e 'devtools::test()' 2>&1 | tail -25
NOT_CRAN=false Rscript -e 'devtools::load_all(quiet = TRUE); t0 <- Sys.time(); testthat::test_file("tests/testthat/test-parallel-chains.R"); print(Sys.time() - t0)'
```

Expected: all pass; the third command runs only the always-on tests (`skip_on_cran()` cases skip because `NOT_CRAN` is false) and must finish in well under 30 seconds. If it does not, move the slowest always-on case (the progress-file fit, then the binary identity) under `skip_on_cran()` and record the measured time in the commit message.

- [ ] **Step 5: Commit**

```bash
git add R/runOccJSDM.R tests/testthat/test-parallel-chains.R
git commit -m "Run chains on PSOCK workers when MCMCparams\$ncores allows it

Workers get the budget's RcppParallel threads and a single-threaded BLAS request, the parent announces the progress file with its tail command, and a fit is identical at every ncores under a declared single-threaded BLAS. Always-on tier measured at N s.

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

Replace `N` with the measured seconds.

---

### Task 6: Documentation inside the hashed set, the full suite, and R CMD check from an export

After this task nothing under `R/`, `src/`, `DESCRIPTION` or `NAMESPACE` changes again before the refit.

**Files:**
- Modify: `R/runOccJSDM.R:316-320` (`@param MCMCparams`), `R/runOccJSDM.R:424-446` (`@return`), `R/runOccJSDM.R:449-455` (`@details`)
- Modify: `R/data.R:48-55` (the `sampleresults` description)
- Modify: `man/runOccJSDM.Rd`, `man/sampleresults.Rd` (regenerated)

- [ ] **Step 1: Replace the `@param MCMCparams` text**

```r
#' @param MCMCparams (Optional) list of MCMC settings: \code{nchain} (number
#' of chains), \code{nburn} (number of burn-in iterations to discard),
#' \code{niter} (number of post-burn-in iterations to keep), \code{nthin}
#' (thinning interval, default \code{1}) and \code{ncores} (the total core
#' budget, default \code{1}). With \code{ncores = 1} the chains run one after
#' another in this process. With \code{ncores > 1} they run on
#' \code{min(nchain, ncores)} worker processes, each with
#' \code{floor(ncores / workers)} RcppParallel threads; for example four chains
#' with \code{ncores = 4} run on four single-threaded workers and with
#' \code{ncores = 8} on four workers of two threads. Defaults to
#' \code{list(nchain = 2, nburn = 5000, niter = 5000, nthin = 1, ncores = 1)}.
```

- [ ] **Step 2: Replace the `@details` paragraph**

```r
#' @details
#' Every chain is seeded on its own: \code{runOccJSDM()} draws two seeds per
#' chain from R's random number generator, one for the R-level draws and one
#' for the package's C++ sampler, so \code{set.seed()} before fitting
#' reproduces a fit on the same platform and build, and the fit is the same
#' whatever \code{ncores} says: a chain's draws depend on its seeds alone, not
#' on which chains ran before it or in which process. Consecutive fits without
#' resetting the seed consume new draws. Because of the per-chain seeding, a
#' seeded fit's numbers differ from those of version 0.1.1, in which chains
#' continued one stream.
#'
#' Random sampling runs serially on R's main thread; deterministic probability
#' calculations may use several RcppParallel threads, and the thread count
#' does not change the sampling stream. The calling session's RcppParallel
#' thread setting is preserved. Identity across \code{ncores} also requires
#' that the BLAS does not change its arithmetic with the thread count: R's
#' reference BLAS is single-threaded, worker processes are asked to run any
#' other BLAS single-threaded, but the package cannot verify which BLAS it is
#' running on, and \code{extSoftVersion()} only indicates it. With an
#' unpinned OpenBLAS, MKL or Accelerate an in-process fit can differ from
#' worker fits.
#'
#' After a fit the session's random number state is restored, so later draws
#' are as if the fit had not happened, with one exception: under the
#' Box-Muller normal kind R keeps a cached deviate outside
#' \code{.Random.seed}, which the fit discards. A parallel fit prints one line
#' naming its progress file, which the workers' iteration messages go to; use
#' the \code{tail -f} command it shows to watch them.
```

- [ ] **Step 3: Extend the `@return` entry for `infos`**

Replace the `\item{infos}{...}` entry with:

```r
#'   \item{infos}{Metadata about the fitted model (e.g. \code{S}, \code{n},
#'   \code{M}, \code{K}, \code{speciesNames}, \code{primerNames},
#'   \code{siteNames}, \code{ncov_psi}, \code{ncov_theta}, the \code{OTU}
#'   matrix used for fitting, and standardisation info for \code{X_psi}/
#'   \code{Xs}), used internally by the plotting/output helpers (e.g.
#'   \code{\link{plotOccupancyRates}}, \code{\link{plotDetectionRates}}). It
#'   also records \code{MCMCparams} (the resolved \code{nchain}, \code{nburn},
#'   \code{niter} and \code{nthin}), \code{chain_seeds} (a 2 by \code{nchain}
#'   matrix of the R and C++ seeds of each chain) and \code{execution} (how
#'   the fit was run: \code{ncores}, \code{workers}, \code{threads_per_worker},
#'   the BLAS thread settings requested of the workers, the BLAS and LAPACK in
#'   use, the elapsed seconds per chain and the progress file). Only
#'   \code{execution} may differ between two fits of the same seed at
#'   different \code{ncores}.}
```

- [ ] **Step 4: Update the `sampleresults` help in `R/data.R`**

Replace lines 50 to 55 (the first paragraph of the `sampleresults` block) with:

```r
#' Output of \code{\link{runOccJSDM}} fit to \code{\link{sampledata}} (10
#' species, 100 sites, up to 3 primers per sample), with
#' \code{MCMCparams = list(nchain = 2, nburn = 5000, niter = 5000, nthin = 1,
#' ncores = 2)}, seed \code{20261007} and R 4.6.1, as in the quickstart's
#' fit-live chunk. Shipped so that \code{vignettes/occJSDM-00-quickstart.Rmd}
#' can load a precomputed fit rather than re-running \code{runOccJSDM()}
#' (which takes several minutes) at vignette build time.
```

- [ ] **Step 5: Regenerate `man/` and `NAMESPACE`**

The repository records roxygen2 8.0.0; 8.1.0 is installed. Try the recorded version first, in a scratch library, so that `man/` does not churn:

```bash
SCRATCH=/private/tmp/claude-501/-Users-douglasyu-src-occJSDM/20c0db62-035a-40e6-b4ff-8dc2f5ddd47d/scratchpad
mkdir -p "$SCRATCH/roxygen800"
Rscript -e 'remotes::install_version("roxygen2", "8.0.0", lib = "'"$SCRATCH"'/roxygen800", repos = "https://cloud.r-project.org", upgrade = "never")'
R_LIBS="$SCRATCH/roxygen800" Rscript -e 'cat(as.character(packageVersion("roxygen2")), "\n"); roxygen2::roxygenise()'
git status --short man NAMESPACE DESCRIPTION
```

Expected: only `man/runOccJSDM.Rd` and `man/sampleresults.Rd` change, and `DESCRIPTION` keeps `RoxygenNote: 8.0.0`. If 8.0.0 cannot be installed, run `Rscript -e 'roxygen2::roxygenise()'` with 8.1.0, accept the `RoxygenNote` bump, inspect `git diff --stat man/` to confirm no unrelated file changed, and say so in the pull request.

- [ ] **Step 6: Run the whole suite under pkgload with the BLAS declaration**

```bash
OCCJSDM_REFERENCE_BLAS=1 Rscript -e 'devtools::test()' 2>&1 | tail -30
```

Expected: 0 failures. Record the counts in the commit message.

- [ ] **Step 7: Commit the hashed-set changes**

```bash
git add R/runOccJSDM.R R/data.R man/runOccJSDM.Rd man/sampleresults.Rd NAMESPACE DESCRIPTION
git commit -m "Document MCMCparams\$ncores, the per-chain seeding contract and the new infos fields

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

- [ ] **Step 8: Build and check from a git archive export, never from the checkout**

```bash
SCRATCH=/private/tmp/claude-501/-Users-douglasyu-src-occJSDM/20c0db62-035a-40e6-b4ff-8dc2f5ddd47d/scratchpad
EXPORT="$SCRATCH/check-export"; rm -rf "$EXPORT"; mkdir -p "$EXPORT/src" "$EXPORT/lib" "$EXPORT/build"
git archive HEAD | tar -x -C "$EXPORT/src"
cd "$EXPORT/build" && R CMD build "$EXPORT/src" > build.log 2>&1 && 'ls' *.tar.gz
cd "$EXPORT/build" && OCCJSDM_REFERENCE_BLAS=1 _R_CHECK_LIMIT_CORES_=TRUE NOT_CRAN=false R CMD check --as-cran occJSDM_*.tar.gz > check.log 2>&1; tail -40 check.log
```

Expected: the recorded baseline, 0 errors, with the same warnings and notes as `AGENTS.md` documents (the `R_ext/Boolean.h` pragma warning and the LICENSE and globals notes); nothing new. This run exercises the `loadNamespace()` branch of the worker bootstrap, because the checked package is installed, not dev-loaded. Keep `check.log` and `build.log` in the scratchpad; Task 9 repeats this from a fresh export for the final evidence.

---

### Task 7: Freeze the build, check BLAS provenance, and run the benchmark

**Files:**
- Create: `dev/simstudy/parallel-chains/README.md`
- Create: `dev/simstudy/parallel-chains/blas-provenance.R`
- Create: `dev/simstudy/parallel-chains/benchmark.R`
- Create: `dev/simstudy/parallel-chains/sample-rss.sh`
- Create: `dev/simstudy/parallel-chains/launch.sh`
- Create: `dev/simstudy/parallel-chains/summarise.R`
- Create: `dev/simstudy/parallel-chains/REPORT.md` (written by hand from `summarise.R`'s output)
- Results go under the ignored `dev/simstudy/results/parallel-chains-<date>/`.

**Interfaces:**
- Consumes: the installed package from the frozen export; `occJSDM:::drop_execution()`.
- Produces: the frozen export directory and library that Task 8 reuses; `results.csv`, `rss.csv`, `provenance.json`; `REPORT.md`.

- [ ] **Step 1: Write the provenance check**

`dev/simstudy/parallel-chains/blas-provenance.R`:

```r
# Establish, by inspection rather than by name, that this R uses its reference
# BLAS, which is single-threaded. Writes the facts as JSON and exits non-zero
# otherwise. Spec: "Tests", the OCCJSDM_REFERENCE_BLAS declaration.
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
out <- args[[1]]
blas <- unname(extSoftVersion()[["BLAS"]])
lapack <- La_library()
info <- Sys.readlink(blas)
is_link <- nzchar(info)
linkage <- if (Sys.info()[["sysname"]] == "Darwin") system2("otool", c("-L", shQuote(blas)), stdout = TRUE) else system2("ldd", shQuote(blas), stdout = TRUE)
foreign <- any(grepl("accelerate|veclib|openblas|mkl|blis|atlas", linkage, ignore.case = TRUE))
defines_dgemm <- if (Sys.info()[["sysname"]] == "Darwin") {
  any(grepl(" T _dgemm_$", system2("nm", c("-gU", shQuote(blas)), stdout = TRUE)))
} else {
  any(grepl(" T dgemm_$", system2("nm", c("-D", shQuote(blas)), stdout = TRUE)))
}
facts <- list(blas = blas, lapack = lapack, blas_is_symlink = is_link, symlink_target = info,
              linkage = linkage, foreign_blas_linked = foreign, defines_dgemm = defines_dgemm,
              blas_md5 = unname(tools::md5sum(blas)), r_version = R.version.string,
              checked_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
reference <- grepl("libRblas", basename(blas)) && !is_link && !foreign && defines_dgemm
facts$verdict <- if (reference) "reference" else "unverified"
writeLines(jsonlite::toJSON(facts, auto_unbox = TRUE, pretty = TRUE), out)
cat("BLAS provenance:", facts$verdict, "\n")
if (!reference) quit(status = 1L)
```

If `jsonlite` is not installed, replace the `writeLines(...)` line with `saveRDS(facts, sub("[.]json$", ".rds", out))` and adjust the readers below.

- [ ] **Step 2: Write the memory sampler**

`dev/simstudy/parallel-chains/sample-rss.sh`:

```sh
#!/bin/sh
# Sample resident memory (KB) of a parent R process and its children once a
# second, until OUT.csv.stop appears or the parent exits.
# Usage: sample-rss.sh PID OUT.csv
PID="$1"; OUT="$2"
echo "time,pid,role,rss_kb" > "$OUT"
while [ ! -f "$OUT.stop" ] && kill -0 "$PID" 2> /dev/null; do
  now=$(date +%s)
  ps -o rss= -p "$PID" | awk -v t="$now" -v p="$PID" '{print t","p",parent,"$1}' >> "$OUT"
  for c in $(pgrep -P "$PID"); do
    ps -o rss= -p "$c" | awk -v t="$now" -v p="$c" '{print t","p",worker,"$1}' >> "$OUT"
  done
  sleep 1
done
```

- [ ] **Step 3: Write the benchmark runner**

`dev/simstudy/parallel-chains/benchmark.R`:

```r
# Time and measure the quickstart fit at several core budgets, twice each, and
# check identity against the in-process fit. Run detached through launch.sh.
# Rscript benchmark.R OUTDIR
args <- commandArgs(TRUE)
stopifnot(length(args) == 1L)
outdir <- normalizePath(args[[1]], mustWork = TRUE)
lib <- file.path(outdir, "library")
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package("occJSDM")) == file.path(lib, "occJSDM"),
          nzchar(Sys.getenv("OCCJSDM_REFERENCE_BLAS")))
sampler <- file.path(outdir, "scripts", "sample-rss.sh")
configs <- list(
  list(label = "serial-ambient", ncores = 1L, pin_one_thread = FALSE),
  list(label = "serial-1thread", ncores = 1L, pin_one_thread = TRUE),
  list(label = "ncores-2", ncores = 2L, pin_one_thread = FALSE),
  list(label = "ncores-4", ncores = 4L, pin_one_thread = FALSE),
  list(label = "ncores-8", ncores = 8L, pin_one_thread = FALSE))
fit_once <- function(ncores) {
  set.seed(20261007)
  fit <- runOccJSDM(data = sampledata, listParams = list(n_factors = 2), threshold = 1,
                    occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
                    collCovariates = c("X_theta.1", "X_theta.2"),
                    spatCovariates = c("Xs.1", "Xs.2"),
                    MCMCparams = list(nchain = 4, nburn = 5000, niter = 5000, nthin = 1, ncores = ncores),
                    summarisedLatentPresences = TRUE)
  list(fit = fit, seed_state = get(".Random.seed", envir = globalenv()), draws = list(rnorm(3), runif(3), sample(100, 3)))
}
reference <- NULL
rows <- list()
for (cfg in configs) for (rep in 1:2) {
  tag <- sprintf("%s-rep%d", cfg$label, rep)
  if (cfg$pin_one_thread) { Sys.setenv(RCPP_PARALLEL_NUM_THREADS = "1"); RcppParallel::setThreadOptions(numThreads = 1L) }
  rss <- file.path(outdir, "rss", paste0(tag, ".csv"))
  system2("sh", c(sampler, Sys.getpid(), shQuote(rss)), wait = FALSE)
  t0 <- Sys.time()
  res <- suppressMessages(suppressWarnings(fit_once(cfg$ncores)))
  wall <- as.numeric(Sys.time() - t0, units = "secs")
  file.create(paste0(rss, ".stop"))               # ends this configuration's sampler
  Sys.sleep(2)
  if (cfg$pin_one_thread) { Sys.unsetenv("RCPP_PARALLEL_NUM_THREADS"); RcppParallel::setThreadOptions(numThreads = "auto") }
  if (is.null(reference)) reference <- res
  identical_fit <- identical(occJSDM:::drop_execution(res$fit), occJSDM:::drop_execution(reference$fit))
  identical_rng <- identical(res$seed_state, reference$seed_state) && identical(res$draws, reference$draws)
  ex <- res$fit$infos$execution
  rows[[tag]] <- data.frame(tag = tag, label = cfg$label, rep = rep, ncores = cfg$ncores, workers = ex$workers,
                            threads_per_worker = ex$threads_per_worker, wall_s = wall,
                            chain_s = paste(round(ex$elapsed, 1), collapse = ";"),
                            identical_fit = identical_fit, identical_rng = identical_rng,
                            blas = ex$blas, stringsAsFactors = FALSE)
  saveRDS(res$fit$infos$execution, file.path(outdir, "execution", paste0(tag, ".rds")))
  write.csv(do.call(rbind, rows), file.path(outdir, "results.csv"), row.names = FALSE)
  cat("DONE", tag, round(wall, 1), "s identical_fit", identical_fit, "\n")
  rm(res); invisible(gc())
}
writeLines(capture.output(sessionInfo()), file.path(outdir, "session-info.txt"))
cat("ALL DONE\n")
```

- [ ] **Step 4: Write the launcher**

`dev/simstudy/parallel-chains/launch.sh`:

```sh
#!/bin/sh
# Freeze the build from a git archive export, check BLAS provenance, run the
# benchmark detached. Run from the worktree root, from a non-interactive shell:
#   nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' \
#     sh dev/simstudy/parallel-chains/launch.sh OUTDIR \
#     > OUTDIR/launch.log 2>&1 < /dev/null &
set -eu
OUT="$1"
REPO="$(pwd -P)"
RSCRIPT=/usr/local/bin/Rscript
mkdir -p "$OUT/src-export" "$OUT/library" "$OUT/scripts" "$OUT/rss" "$OUT/execution"
cp "$REPO"/dev/simstudy/parallel-chains/*.R "$REPO"/dev/simstudy/parallel-chains/*.sh "$OUT/scripts/"
git -C "$REPO" rev-parse HEAD > "$OUT/revision.txt"
git -C "$REPO" archive HEAD | tar -x -C "$OUT/src-export"
R CMD INSTALL --library="$OUT/library" "$OUT/src-export" > "$OUT/install.log" 2>&1
"$RSCRIPT" "$OUT/scripts/blas-provenance.R" "$OUT/provenance.json"
export OCCJSDM_REFERENCE_BLAS=1
caffeinate -i -w $$ &
echo $$ > "$OUT/launcher.pid"
"$RSCRIPT" "$OUT/scripts/benchmark.R" "$OUT" > "$OUT/benchmark.log" 2>&1
echo "ALL DONE $(date '+%F %T')"
```

The `export OCCJSDM_REFERENCE_BLAS=1` line runs only if `blas-provenance.R` exited zero, because of `set -e`.

- [ ] **Step 5: Write the summariser**

`dev/simstudy/parallel-chains/summarise.R`:

```r
# Rscript summarise.R OUTDIR  -> prints the bullet lists for REPORT.md
args <- commandArgs(TRUE)
outdir <- normalizePath(args[[1]], mustWork = TRUE)
r <- read.csv(file.path(outdir, "results.csv"), stringsAsFactors = FALSE)
peak <- function(tag, role) {
  f <- file.path(outdir, "rss", paste0(tag, ".csv"))
  x <- read.csv(f); x <- x[x$role == role, ]
  if (!nrow(x)) return(NA_real_)
  by_t <- tapply(x$rss_kb, x$time, sum)           # all workers alive at that second, summed
  max(by_t) / 1024
}
r$parent_peak_mb <- vapply(r$tag, peak, numeric(1), role = "parent")
r$workers_peak_mb <- vapply(r$tag, peak, numeric(1), role = "worker")
base <- mean(r$wall_s[r$label == "serial-ambient"])
for (lab in unique(r$label)) {
  s <- r[r$label == lab, ]
  cat(sprintf("- **%s** (ncores %d, %d worker%s, %d thread%s each): wall %.0f s and %.0f s (speed-up %.2fx over serial-ambient), parent peak %.0f MB, workers peak %.0f MB, identical fit %s, identical RNG %s.\n",
              lab, s$ncores[1], s$workers[1], if (s$workers[1] == 1) "" else "s", s$threads_per_worker[1],
              if (s$threads_per_worker[1] == 1) "" else "s", s$wall_s[1], s$wall_s[2], base / mean(s$wall_s),
              max(s$parent_peak_mb), max(s$workers_peak_mb, na.rm = TRUE),
              paste(s$identical_fit, collapse = "/"), paste(s$identical_rng, collapse = "/")))
}
```

- [ ] **Step 6: Launch, from a non-interactive shell, on mains power with the lid open**

```bash
cd /Users/douglasyu/src/occJSDM/.worktrees/parallel-chains
OUT=/Users/douglasyu/src/occJSDM/dev/simstudy/results/parallel-chains-$(date +%Y%m%d)
mkdir -p "$OUT"
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh dev/simstudy/parallel-chains/launch.sh "$OUT" > "$OUT/launch.log" 2>&1 < /dev/null &
sleep 5; cat "$OUT/launch.log"; cat "$OUT/provenance.json" | head -20
```

Expected: `BLAS provenance: reference`, then the benchmark starts. Wait with a `Monitor` or by checking `tail -3 "$OUT/benchmark.log"` every 20 minutes; it ends with `ALL DONE`. Roughly 90 minutes.

- [ ] **Step 7: Summarise and write the report**

```bash
Rscript dev/simstudy/parallel-chains/summarise.R "$OUT"
```

Write `dev/simstudy/parallel-chains/REPORT.md` by hand from the output, one line per paragraph, bullet lists only, no tables, no em-dashes, with these sections: what was run (the quickstart call, four chains, the five configurations twice, this laptop, R 4.6.1, the frozen revision from `revision.txt`, the install from the export); BLAS provenance (every fact from `provenance.json`, stated as facts, the verdict, and that identity claims hold on this BLAS); results (the bullet per configuration from `summarise.R`); the chains-outnumber-workers case (the `ncores-2` run's worker and parent peaks on their own); limits (Windows not run, no threaded BLAS available, two repeats only). Copy `results.csv` and `provenance.json` next to the report. Write `README.md` in the same folder saying how to launch and where results land.

- [ ] **Step 8: Commit**

```bash
git add dev/simstudy/parallel-chains/
git commit -m "Benchmark parallel chains on the quickstart fit and record BLAS provenance

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 8: Refit `sampleresults` from the frozen build

Nothing under `R/`, `src/`, `DESCRIPTION` or `NAMESPACE` has changed since Task 6's commit; verify with `git log --oneline -1 -- R src DESCRIPTION NAMESPACE` before starting. The frozen library from Task 7 (`$OUT/library`) is reused.

**Files:**
- Create: `dev/release/2026-10-13-parallel-chains-refit/refit-sampleresults.R`
- Create: `dev/release/2026-10-13-parallel-chains-refit/verify-release.R`
- Create: `dev/release/2026-10-13-parallel-chains-refit/refit-README.md`
- Create: `dev/release/2026-10-13-parallel-chains-refit/REPORT.md`
- Modify: `data/sampleresults.rda` (replaced by the refit)
- Evidence files written by the scripts: `provenance.rds`, `session-info.txt`, `diagnostics.csv`, `diagnostics-flagged.csv`, `blas-provenance.json`.

- [ ] **Step 1: Write the refit script**

`dev/release/2026-10-13-parallel-chains-refit/refit-sampleresults.R` (modelled on the 7 October script; the differences are the two fits, the identity gate, and the new provenance fields):

```r
# Refit the quickstart's shipped example under per-chain seeding, twice (one
# core and two cores), require identity apart from the execution record, and
# ship the two-core object. Rscript refit-sampleresults.R REPO LIB EVIDENCE PROVENANCE_JSON
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
repo <- normalizePath(args[[1]], mustWork = TRUE)
lib <- normalizePath(args[[2]], mustWork = TRUE)
evidence <- normalizePath(args[[3]], mustWork = TRUE)
blas_json <- normalizePath(args[[4]], mustWork = TRUE)
stopifnot(!file.exists(file.path(evidence, "provenance.rds")), nzchar(Sys.getenv("OCCJSDM_REFERENCE_BLAS")))
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(occJSDM))
stopifnot(normalizePath(find.package("occJSDM")) == file.path(lib, "occJSDM"))
RNGkind("Mersenne-Twister", "Inversion", "Rejection")
seed <- 20261007L
settings <- list(
  listParams = list(n_factors = 2), threshold = 1,
  occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
  collCovariates = c("X_theta.1", "X_theta.2"),
  spatCovariates = c("Xs.1", "Xs.2"),
  summarisedLatentPresences = TRUE, listPriors = list())
mcmc <- list(nchain = 2, nburn = 5000, niter = 5000, nthin = 1)
source_files <- system2("git", c("-C", shQuote(repo), "ls-files", "R", "src", "DESCRIPTION", "NAMESPACE"), stdout = TRUE)
source_hashes <- tools::md5sum(file.path(repo, source_files)); names(source_hashes) <- source_files
old_fit_hash <- unname(tools::md5sum(file.path(repo, "data/sampleresults.rda")))
input <- new.env(parent = emptyenv()); load(file.path(repo, "data/sampledata.rda"), envir = input)
stopifnot(is.list(input$sampledata), ncol(input$sampledata$OTU) == 10L)
fit_at <- function(ncores) {
  warnings_seen <- character(); started <- Sys.time()
  set.seed(seed)
  fit <- withCallingHandlers(
    do.call(runOccJSDM, c(list(data = input$sampledata), settings, list(MCMCparams = c(mcmc, list(ncores = ncores))))),
    warning = function(w) warnings_seen <<- c(warnings_seen, conditionMessage(w)))
  list(fit = fit, warnings = warnings_seen, started = started, finished = Sys.time(),
       seed_state = get(".Random.seed", envir = globalenv()))
}
one <- fit_at(1L)
two <- fit_at(2L)
stopifnot(identical(occJSDM:::drop_execution(one$fit), occJSDM:::drop_execution(two$fit)),
          identical(one$seed_state, two$seed_state),
          two$fit$infos$execution$workers == 2L, two$fit$infos$execution$threads_per_worker == 1L,
          one$fit$infos$execution$workers == 1L)
sampleresults <- two$fit
stopifnot(
  identical(sampleresults$infos$model, "two_stage"), sampleresults$infos$n == 100L,
  sampleresults$infos$S == 10L, sampleresults$infos$n_factors == 2L, sampleresults$infos$threshold == 1,
  identical(dim(sampleresults$results_output$jsdm_output$B_output), c(2L, 10L, 5000L, 2L)),
  all(is.finite(sampleresults$results_output$psi_output)),
  identical(source_hashes, stats::setNames(tools::md5sum(file.path(repo, source_files)), source_files)))
diagnostics <- returnConvergenceDiagnostics(sampleresults)
write.csv(diagnostics, file.path(evidence, "diagnostics.csv"), row.names = FALSE)
candidate <- tempfile(".sampleresults-", tmpdir = file.path(repo, "data"), fileext = ".rda")
save(sampleresults, file = candidate, compress = "xz")
stopifnot(file.rename(candidate, file.path(repo, "data/sampleresults.rda")))
file.copy(blas_json, file.path(evidence, "blas-provenance.json"), overwrite = TRUE)
provenance <- list(
  revision = system2("git", c("-C", shQuote(repo), "rev-parse", "HEAD"), stdout = TRUE),
  source_hashes = source_hashes,
  input_md5 = unname(tools::md5sum(file.path(repo, "data/sampledata.rda"))),
  previous_fit_md5 = old_fit_hash,
  fit_md5 = unname(tools::md5sum(file.path(repo, "data/sampleresults.rda"))),
  seed = seed, rng_kind = RNGkind(), settings = settings, MCMCparams = c(mcmc, list(ncores = 2L)),
  chain_seeds = sampleresults$infos$chain_seeds,
  execution_one_core = one$fit$infos$execution, execution_two_cores = two$fit$infos$execution,
  identical_apart_from_execution = TRUE,
  blas_provenance_file = "blas-provenance.json",
  elapsed_seconds = c(one_core = as.numeric(difftime(one$finished, one$started, units = "secs")),
                      two_cores = as.numeric(difftime(two$finished, two$started, units = "secs"))),
  warnings = list(one_core = one$warnings, two_cores = two$warnings),
  installed_library = lib, session = sessionInfo())
saveRDS(provenance, file.path(evidence, "provenance.rds"), compress = "xz")
writeLines(capture.output(sessionInfo()), file.path(evidence, "session-info.txt"))
flagged <- diagnostics[!is.finite(diagnostics$rhat) | diagnostics$rhat > 1.01 | !is.finite(diagnostics$ess) | diagnostics$ess < 400, ]
write.csv(flagged, file.path(evidence, "diagnostics-flagged.csv"), row.names = FALSE)
cat("REFIT COMPLETE: ", nrow(diagnostics), " coefficient diagnostics; ", nrow(flagged), " flagged.\n", sep = "")
```

- [ ] **Step 2: Write the verify script**

`dev/release/2026-10-13-parallel-chains-refit/verify-release.R`:

```r
# Check the refitted data against its provenance, optionally against an installation.
# Rscript verify-release.R REPO [LIB]
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% c(1L, 2L))
repo <- normalizePath(args[[1]], mustWork = TRUE)
evidence <- file.path(repo, "dev/release/2026-10-13-parallel-chains-refit")
p <- readRDS(file.path(evidence, "provenance.rds"))
stopifnot(
  p$seed == 20261007L,
  identical(p$MCMCparams, list(nchain = 2, nburn = 5000, niter = 5000, nthin = 1, ncores = 2L)),
  isTRUE(p$identical_apart_from_execution),
  p$execution_two_cores$workers == 2L, p$execution_one_core$workers == 1L,
  identical(p$execution_two_cores$blas_threads_requested,
            c(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1",
              MKL_DOMAIN_NUM_THREADS = "MKL_DOMAIN_ALL=1", VECLIB_MAXIMUM_THREADS = "1")),
  p$fit_md5 == unname(tools::md5sum(file.path(repo, "data/sampleresults.rda"))),
  p$input_md5 == unname(tools::md5sum(file.path(repo, "data/sampledata.rda"))),
  identical(p$source_hashes, stats::setNames(tools::md5sum(file.path(repo, names(p$source_hashes))), names(p$source_hashes))))
e <- new.env(parent = emptyenv()); load(file.path(repo, "data/sampleresults.rda"), envir = e)
fit <- e$sampleresults
stopifnot(length(fit) == 6L, fit$infos$model == "two_stage", fit$infos$n == 100L,
          identical(fit$infos$chain_seeds, p$chain_seeds),
          identical(dim(fit$results_output$jsdm_output$B_output), c(2L, 10L, 5000L, 2L)))
if (length(args) == 2L) {
  lib <- normalizePath(args[[2]], mustWork = TRUE)
  .libPaths(c(lib, .libPaths()))
  suppressPackageStartupMessages(library(occJSDM))
  stopifnot(normalizePath(find.package("occJSDM")) == file.path(lib, "occJSDM"))
  stopifnot(identical(occJSDM::sampleresults, fit))
}
cat("Refit provenance, source hashes, seed, settings and data identity verified.\n")
```

- [ ] **Step 3: Write the README and run the refit detached**

`refit-README.md`: the launch recipe below, one line per paragraph, and the statement that the script replaces `data/sampleresults.rda` in the checkout it is given, so it must be run in the worktree. Then launch, reusing Task 7's frozen library and provenance file (`$OUT` from Task 7):

```bash
cd /Users/douglasyu/src/occJSDM/.worktrees/parallel-chains
EV=/Users/douglasyu/src/occJSDM/dev/simstudy/results/parallel-chains-refit-$(date +%Y%m%d); mkdir -p "$EV"
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' env OCCJSDM_REFERENCE_BLAS=1 /usr/local/bin/Rscript --vanilla dev/release/2026-10-13-parallel-chains-refit/refit-sampleresults.R "$(pwd -P)" "$OUT/library" "$EV" "$OUT/provenance.json" > "$EV/refit.log" 2>&1 < /dev/null &
```

Wait for `REFIT COMPLETE` in `$EV/refit.log` (two fits of the quickstart size, so roughly 20 to 30 minutes). Then copy `provenance.rds`, `session-info.txt`, `diagnostics.csv`, `diagnostics-flagged.csv` and `blas-provenance.json` from `$EV` into `dev/release/2026-10-13-parallel-chains-refit/`, and run the verify script against the frozen library:

```bash
cp "$EV"/provenance.rds "$EV"/session-info.txt "$EV"/diagnostics*.csv "$EV"/blas-provenance.json dev/release/2026-10-13-parallel-chains-refit/
Rscript dev/release/2026-10-13-parallel-chains-refit/verify-release.R "$(pwd -P)"
```

Expected: "Refit provenance, source hashes, seed, settings and data identity verified." (The installed-library form of the check runs in Task 9 against the fresh export, because the frozen library holds the pre-refit data.)

- [ ] **Step 4: Write `REPORT.md`** for the refit folder, in the 7 October report's shape: what was refitted and why (per-chain seeding changes every seeded fit), the settings, the two fits and their identity apart from the execution record, the chain seeds, the BLAS provenance facts, the elapsed times, the warnings seen, the diagnostics summary (how many flagged at Rhat above 1.01 or ESS below 400, read from `diagnostics-flagged.csv`), the md5s before and after, and the statement that the 7 October verify script now fails its source-hash check by design.

- [ ] **Step 5: Commit**

```bash
git add data/sampleresults.rda dev/release/2026-10-13-parallel-chains-refit/
git commit -m "Refit the shipped sampleresults under per-chain seeding with recorded provenance

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 9: Documentation outside the hashed set, and the final validation from a fresh export

**Files:**
- Modify: `NEWS.md`
- Modify: `vignettes/occJSDM-00-quickstart.Rmd`, `vignettes/occJSDM-00-quickstart.md`, `vignettes/occJSDM-00-quickstart_files/` (rebuilt)
- Modify: `AGENTS.md`

- [ ] **Step 1: NEWS.md**

Under `# occJSDM (development version)` add:

```markdown
## New features

-   `runOccJSDM()` can run its chains on separate cores: set `ncores` in `MCMCparams`. `ncores` is the total core budget, so `ncores = 4` runs four chains on four single-threaded workers and `ncores = 8` on four workers of two RcppParallel threads each. A parallel fit prints the `tail -f` command for its progress file. The returned `infos` records the sampling settings, the seeds of every chain and how the fit was run.

## Behaviour changes

-   Every chain is now seeded on its own, so a fit is identical whatever `ncores` says and no chain continues another's random stream. Because of this, a fit under a given `set.seed()` value gives different numbers from version 0.1.1; the shipped `sampleresults` was refitted under the new seeding with the same seed and settings. Identity across `ncores` is stated for a single-threaded BLAS; see `?runOccJSDM`.
```

- [ ] **Step 2: The quickstart**

In `vignettes/occJSDM-00-quickstart.Rmd`: at line 54, change `R 4.5.0` to `R 4.6.1` in the "same seed" sentence. At line 78, replace the `MCMCparams` bullet with:

```markdown
- `MCMCparams` sets the number of chains, the burn-in and kept iterations, and the thinning. Adding `ncores = 2` to the list runs the two chains on two cores at once; `ncores` is the total number of cores the fit may use.
```

Then re-read every sentence that describes the shipped fit (line 134's warnings paragraph and the diagnostics discussion around it) against the refitted object: load `data/sampleresults.rda`, run `returnConvergenceDiagnostics(sampleresults)` and read the fit-time warnings from `provenance.rds` (`p$warnings$two_cores`), and correct any block name or number the text states. The quickstart loads `sampleresults` from the installed package, so install the branch head, which now carries the refitted data, into a scratch library from an export first, then render with that library in front:

```bash
SCRATCH=/private/tmp/claude-501/-Users-douglasyu-src-occJSDM/20c0db62-035a-40e6-b4ff-8dc2f5ddd47d/scratchpad
DOC="$SCRATCH/doc-export"; rm -rf "$DOC"; mkdir -p "$DOC/src" "$DOC/lib"
git archive HEAD | tar -x -C "$DOC/src"
R CMD INSTALL --library="$DOC/lib" "$DOC/src" > "$DOC/install.log" 2>&1
R_LIBS="$DOC/lib" Rscript -e 'rmarkdown::render("vignettes/occJSDM-00-quickstart.Rmd", output_format = rmarkdown::github_document(toc = TRUE, html_preview = FALSE, pandoc_args = "--wrap=none"), output_file = "occJSDM-00-quickstart.md", output_dir = "vignettes")'
git diff --stat vignettes/
```

Expected: the `.md` changes only where the text or the fit's numbers changed, plus the regenerated figures. If the diff shows format churn throughout the file, the committed `.md` was produced with different options: find them with `git log -p --follow -- vignettes/occJSDM-00-quickstart.md | head -80`, match them, and re-render.

- [ ] **Step 3: AGENTS.md**

After the paragraph that begins `**What it does not fix.** Bit-reproducibility above one thread.` (line 168), insert one paragraph:

```markdown
**Superseded, 13 October 2026.** The paragraph above, and the Fixed bugs 41 advice in `TODO.md` to set `RCPP_PARALLEL_NUM_THREADS=1` for reproducible runs, describe the August 2026 state. Fixed bugs 50 (14 September 2026) moved every random draw onto the main thread, and `tests/testthat/test-rng-safety.R` asserts that fits are bit-identical at 1 and 4 requested threads for all four model types; it was re-run on the installed 0.1.1 build on 13 October 2026 and passes. Thread count no longer affects results. Parallel chains (`MCMCparams$ncores`, the parallel-chains spec under `dev/superpowers/specs/`) seed every chain on its own; what still affects reproducibility is the BLAS, as `?runOccJSDM` explains.
```

- [ ] **Step 4: Final validation from a fresh export**

Commit the documentation first (Step 5 below can be reordered: commit, then validate, then amend nothing). Then:

```bash
git add NEWS.md vignettes/occJSDM-00-quickstart.Rmd vignettes/occJSDM-00-quickstart.md vignettes/occJSDM-00-quickstart_files AGENTS.md
git commit -m "Document parallel chains in NEWS and the quickstart, and mark the stale thread note superseded

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
SCRATCH=/private/tmp/claude-501/-Users-douglasyu-src-occJSDM/20c0db62-035a-40e6-b4ff-8dc2f5ddd47d/scratchpad
FINAL="$SCRATCH/final-export"; rm -rf "$FINAL"; mkdir -p "$FINAL/src" "$FINAL/lib" "$FINAL/build"
git archive HEAD | tar -x -C "$FINAL/src"
cd "$FINAL/build" && R CMD build "$FINAL/src" > build.log 2>&1
cd "$FINAL/build" && OCCJSDM_REFERENCE_BLAS=1 _R_CHECK_LIMIT_CORES_=TRUE NOT_CRAN=false R CMD check --as-cran occJSDM_*.tar.gz > check.log 2>&1; tail -40 check.log
R CMD INSTALL --library="$FINAL/lib" "$FINAL"/build/occJSDM_*.tar.gz > "$FINAL/install.log" 2>&1
cd /Users/douglasyu/src/occJSDM/.worktrees/parallel-chains
Rscript dev/release/2026-10-13-parallel-chains-refit/verify-release.R "$(pwd -P)" "$FINAL/lib"
OCCJSDM_REFERENCE_BLAS=1 Rscript -e 'devtools::test()' 2>&1 | tail -20
```

Expected: the check matches the recorded baseline (0 errors; the known warnings and notes, nothing new); the verify script prints its success line twice (provenance, then installed data); the pkgload suite passes. Save `check.log` into `dev/release/2026-10-13-parallel-chains-refit/` as `final-check.log` and commit it with a one-line addition to that folder's `REPORT.md` naming the checked revision.

---

### Task 10: The pull request

- [ ] **Step 1: Push and open the PR**

```bash
git push -u origin claude/parallel-chains
gh pr create --repo AlexDiana/occJSDM --base main --head claude/parallel-chains --title "Run MCMC chains in parallel with per-chain seeds (MCMCparams\$ncores)" --body-file /private/tmp/claude-501/-Users-douglasyu-src-occJSDM/20c0db62-035a-40e6-b4ff-8dc2f5ddd47d/scratchpad/pr-body.md
```

Write `pr-body.md` first, in the PR #14 style, with these sections and nothing invented: **What this adds** (`ncores`, the budget rule, the progress file, the `infos` fields); **The seeding contract, for Alex** (two seeds per chain drawn from R's RNG, pinned kinds, why a chain's draws depend on its seeds alone, why every seeded fit's numbers change, and that `sampleresults` was refitted); **Exact combination** (draws by chain index; latent sums; the pairwise WAIC merge formula written out, with its unit test against the single pass); **Workers** (PSOCK, bootstrap, the four-part build check, one chain per task, shutdown with `pskill`); **Reproducibility limits** (BLAS: the promise is for a single-threaded BLAS, the package cannot verify which BLAS it runs on, the test gate and the provenance check); **Benchmark** (the headline speed-up and memory overhead copied from `dev/simstudy/parallel-chains/REPORT.md`, with the BLAS it was measured on); **How to review the moved loop** (`run_one_chain()` is lines 1055 to 1365 of the old file moved verbatim plus five listed edits; point at the commit); **Checks** (the final `R CMD check` result, the two suite runs, the verify script); the roxygen2 version used; and the closing line `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.

- [ ] **Step 2: Bind the PR in the app and read CI**

Use the ccd_pr tools: `get_status`, and `bind_pr` if the PR is not reported. The repository has no CI on pull requests (`pkgdown.yaml` is manual), so there is nothing to wait for. Do not enable auto-merge.

- [ ] **Step 3: Record the state**

Update the memory file `parallel-chains-spec-status.md` to say the PR is open and awaiting Alex, and that `TODO.md` is edited only after the merge, in its own commit on `main`, as the spec's "Documentation changes" describes.

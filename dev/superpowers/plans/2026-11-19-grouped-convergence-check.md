# Grouped Prior Convergence Check Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give fits made with `embeddingVars` a convergence check that says whether the chains found gradients and found the same ones, as a function and in the report printed at the end of a fit.

**Architecture:** Pure functions in a new file compute four checks from pieces of a fit (draws of the scales, draws of the measured coefficients, and each chain's posterior mean of the embedding coefficients with the design matrix). One exported function assembles them from a fitted model. `computeDiagnostics()` gains one optional argument, and `runOccJSDM()` passes it and prints the new summary.

**Tech Stack:** R, testthat edition 3, roxygen2 (Rd files regenerated in a scratch copy of the package, never with `devtools::load_all()` in the worktree).

**Spec:** `dev/superpowers/specs/2026-11-19-grouped-convergence-check-design.md`

## Global Constraints

- Branch `claude/grouped-mixing`, worktree `.worktrees/grouped-mixing`. Commit locally in small commits; never push.
- A fit made without `embeddingVars` must print and warn exactly as before. `computeDiagnostics(results_output)` with no second argument is unchanged.
- No random numbers are drawn by any new code.
- Limits, as fixed: state 1; R-hat 1.1; gradients 0.90; disagreement 0.15; small effect 0.5; share of the effect that sets the number of gradients 0.90. The two provisional limits of the spec (0.90 and 0.15) were measured against the 18 saved fits before this plan was written and stand as stated (Task 5 records the figures).
- The word for a flagged fit is "review".
- No em-dashes anywhere. Markdown files: one line per paragraph, no pipe tables.
- Tests are written first and watched to fail.
- Run tests against an installed copy: `<scratchpad>/devtest.sh <filter>` (syncs `R/` and `tests/` to a scratch build, installs, runs `testthat::test_dir(..., load_package = "installed")`).

## Review Focus

- One chain only: R-hat and the two effect checks cannot be computed. Expected: the state is reported, the rest are not assessed, the verdict rests on the state.
- Every occupancy covariate is an embedding column: there are no measured coefficients and `sigma_b` is a prior draw. Expected: the R-hat check covers the embedding scales only and says nothing about measured coefficients.
- The embedding's effect is exactly zero in every chain. Expected: checks 3 and 4 are not assessed, with no division by zero and no `NaN` in the output.
- Thousands of species. Expected: each chain's mean of B is taken without copying the array of draws per coefficient (`rowMeans(x, dims = 2)` on one chain's slice).
- A fit of another model type made with `embeddingVars` (occupancy, two-stage, continuous). Expected: the function works, since nothing in it depends on the model type.

---

### Task 1: The two effect checks

**Files:**
- Create: `R/grouped-convergence.R`
- Test: `tests/testthat/test-grouped-convergence.R`

**Interfaces:**
- Produces: `GROUPED_LIMITS`, a list with `state = 1`, `rhat = 1.1`, `gradients = 0.90`, `disagreement = 0.15`, `small_effect = 0.5`, `hold = 0.90`.
- Produces: `grouped_effect_checks(E, limits = GROUPED_LIMITS)`, where `E` is a sites by species by chains array of each chain's posterior mean of the embedding's combined effect. Returns a list: `assessed` (logical), `why_not` (character, `NA` when assessed), `size` (median over species of the sd across sites of the pooled effect), `n_gradients` (integer), `agreement` (the smallest share over ordered pairs of chains), `agreement_pairs` (chains by chains matrix, `NA` on the diagonal; entry a, b is the share of chain a's effect inside chain b's gradients), `disagreement` (mean over chains), `disagreement_by_chain`.

- [ ] **Step 1: Write the failing tests** in `tests/testthat/test-grouped-convergence.R`

```r
effect_from <- function(scores, loadings) scores %*% t(loadings)
set.seed(1); n <- 40L; S <- 12L
U <- qr.Q(qr(matrix(rnorm(n * 3), n, 3))); V <- matrix(rnorm(S * 3), S, 3)

test_that("identical chains agree fully and do not disagree", {
  E1 <- effect_from(U[, 1:2] %*% diag(c(9, 5)), V[, 1:2])
  got <- occJSDM:::grouped_effect_checks(array(c(E1, E1, E1), c(n, S, 3L)))
  expect_true(got$assessed)
  expect_identical(got$n_gradients, 2L)
  expect_equal(got$agreement, 1); expect_equal(got$disagreement, 0)
  expect_true(all(is.na(diag(got$agreement_pairs))))
})

test_that("chains with the same gradients mixed differently agree fully", {
  turn <- matrix(c(cos(0.7), sin(0.7), -sin(0.7), cos(0.7)), 2, 2)
  scores <- U[, 1:2] %*% diag(c(9, 5))
  E1 <- effect_from(scores, V[, 1:2])
  E2 <- effect_from(scores %*% turn, V[, 1:2] %*% turn)      # the same effect, other blends
  got <- occJSDM:::grouped_effect_checks(array(c(E1, E2), c(n, S, 2L)))
  expect_equal(got$agreement, 1); expect_equal(got$disagreement, 0)
})

test_that("a chain with other gradients is seen by both checks", {
  E1 <- effect_from(U[, 1:2] %*% diag(c(9, 5)), V[, 1:2])
  E3 <- effect_from(U[, c(1, 3)] %*% diag(c(9, 5)), V[, 1:2])   # second gradient replaced
  got <- occJSDM:::grouped_effect_checks(array(c(E1, E1, E3), c(n, S, 3L)))
  expect_lt(got$agreement, 0.90); expect_gt(got$disagreement, 0.15)
  expect_equal(got$agreement_pairs[1, 2], 1)
  expect_lt(got$agreement_pairs[1, 3], 0.90)
})

test_that("a small or absent effect, or a single chain, is not assessed", {
  E1 <- effect_from(U[, 1:2] %*% diag(c(9, 5)), V[, 1:2])
  small <- occJSDM:::grouped_effect_checks(array(c(E1, E1) * 1e-3, c(n, S, 2L)))
  expect_false(small$assessed); expect_match(small$why_not, "small")
  zero <- occJSDM:::grouped_effect_checks(array(0, c(n, S, 2L)))
  expect_false(zero$assessed); expect_identical(zero$size, 0)
  expect_false(anyNA(c(zero$size, zero$n_gradients)) && FALSE)
  expect_true(is.na(zero$agreement) && is.na(zero$disagreement))
  one <- occJSDM:::grouped_effect_checks(array(E1, c(n, S, 1L)))
  expect_false(one$assessed); expect_match(one$why_not, "chain")
})
```

- [ ] **Step 2: Run** `devtest.sh grouped-convergence`. Expected: errors, `grouped_effect_checks` not found.

- [ ] **Step 3: Implement** `GROUPED_LIMITS` and `grouped_effect_checks()` in `R/grouped-convergence.R`.

The algorithm, which the tests do not fully determine:

```r
pooled <- apply(E, c(1, 2), mean)
size <- stats::median(apply(pooled, 2, stats::sd))
# not assessed: fewer than 2 chains ("needs two chains or more"), or size < limits$small_effect
# ("the embedding's effect is small")
d <- svd(pooled, nu = 0, nv = 0)$d
n_gradients <- which(cumsum(d^2) / sum(d^2) >= limits$hold)[1]
U_c <- svd(E[, , c], nu = n_gradients, nv = 0)$u            # chain c's site gradients
agreement_pairs[a, b] <- sum((t(U_b) %*% E[, , a])^2) / sum(E[, , a]^2)
disagreement_by_chain[c] <- sqrt(mean((E[, , c] - pooled)^2)) / sqrt(mean(pooled^2))
```

- [ ] **Step 4: Run** `devtest.sh grouped-convergence`. Expected: all pass.

- [ ] **Step 5: Commit** `R/grouped-convergence.R` and the test file.

### Task 2: The state and R-hat checks, the verdict, and the exported function

**Files:**
- Modify: `R/grouped-convergence.R`
- Test: `tests/testthat/test-grouped-convergence.R`
- Create by roxygen: `man/returnGroupedConvergence.Rd`; modify `NAMESPACE`

**Interfaces:**
- Consumes: `grouped_effect_checks()`, `GROUPED_LIMITS`; `blockRhatESS(arr, n_params, niter, nchain)` from `R/diagnostics.R`; `is_grouped_fit()` from `R/output.R`.
- Produces: `grouped_convergence(B_output, sigmab_output, sigmaemb_output, lambda_output, X, is_emb, limits = GROUPED_LIMITS)`. `B_output` is covariates by species by iterations by chains; the three scale arguments are as stored in `results_output$jsdm_output`; `X` is the design matrix; `is_emb` a logical per design column. Returns a list: `verdict` ("ok" or "review"), `checks` (data frame with columns `check`, `value`, `limit`, `flagged`, `meaning`, and rows `state`, `rhat_measured`, `rhat_scales`, `gradients`, `disagreement`, in that order), `chains` (data frame: `chain`, `sigma_emb_median`, `state`, with `state` either `"shared structure"` or `"leftovers dominate"`), `effect` (what `grouped_effect_checks()` returns).
- Produces: exported `returnGroupedConvergence(fitModel)`, which stops with "returnGroupedConvergence() needs a fit with embeddingVars" for any other fit and otherwise returns `grouped_convergence()` of the fit's pieces.

Rules the tests pin:
- `state`: `value` is the number of chains whose median `sigma_emb` is at least 1. Flagged when it is neither 0 nor the number of chains ("the chains are in different states"), and when it equals the number of chains ("in every chain the leftovers are at least as large as the shared part; right if species respond independently to the embedding, wrong if the chains missed shared gradients, which is the likely reading with fewer than 5,000 burn-in iterations or a model type other than binary").
- `rhat_measured`: the largest R-hat over the measured rows of B. `NA` and not flagged when there are no measured rows or one chain.
- `rhat_scales`: the largest R-hat over `log(sigma_emb)`, the logs of the axis scales and, when there are measured rows, `log(sigma_b)`.
- `gradients` and `disagreement`: from the effect checks; `value` `NA`, not flagged, and `meaning` giving `why_not` when not assessed.
- `verdict` is "review" when any row is flagged.

- [ ] **Step 1: Write the failing tests.** Build the inputs by hand: `p = 5` design columns of which the last three are embedding, `S = 6`, 200 iterations. A helper `made_fit(nchain, sigma_emb_by_chain, shift_measured = 0, gradients = "same")` returns the argument list for `grouped_convergence()`: `B_output` whose embedding rows are a rank-2 pattern plus small noise per draw (the same pattern in every chain, or with the second direction replaced in the last chain), measured rows that are noise plus `shift_measured` in the last chain, `sigmaemb_output` lognormal draws around the given medians, `lambda_output` and `sigmab_output` lognormal draws common to the chains, and `X` a 30 by 5 standardised random matrix. Tests:
  - all chains at medians 0.3: `verdict` "ok", no row flagged, `chains$state` all "shared structure", `checks$check` equal to the five names in order.
  - medians `c(0.3, 6)`: `state` flagged, `value` 1, `meaning` matches "different states", `verdict` "review".
  - medians `c(6, 7)`: `state` flagged, `value` 2, `meaning` matches "respond independently".
  - `shift_measured = 5`: `rhat_measured` flagged with `value` over 1.1.
  - `gradients = "other"` with three chains: `gradients` flagged.
  - one chain: `rhat_measured`, `rhat_scales`, `gradients`, `disagreement` all `NA` and not flagged; `state` still reported.
  - `is_emb` all `TRUE`: `rhat_measured` `NA`, `rhat_scales` finite.
  - `returnGroupedConvergence()` on a grouped binary fit (`simulate_grouped_fixture(p_meas = 2L, p_emb = 8L, S = 12L)`, `FIXTURE_MCMC`) returns the four names and a five-row `checks`; on a fit without `embeddingVars` it stops with the message above; on a grouped occupancy fit it returns without error.

- [ ] **Step 2: Run** `devtest.sh grouped-convergence`. Expected: the new tests fail with `grouped_convergence` not found.

- [ ] **Step 3: Implement** `grouped_convergence()` and `returnGroupedConvergence()`, with roxygen for the exported function (`@export`, `@seealso` the grouped-prior page; the Value section describes the four elements and names the limits). Each chain's mean of the embedding rows is `rowMeans(B_output[is_emb, , , c, drop = FALSE], dims = 2)`.

- [ ] **Step 4: Regenerate** `NAMESPACE` and `man/returnGroupedConvergence.Rd` in the scratch build with `roxygen2::roxygenise(".", load_code = roxygen2::load_installed)`, check that only those two files and no other Rd changed, and copy them to the worktree.

- [ ] **Step 5: Run** `devtest.sh grouped-convergence`. Expected: all pass.

- [ ] **Step 6: Commit.**

### Task 3: The automatic report

**Files:**
- Modify: `R/diagnostics.R` (`computeDiagnostics()`), `R/grouped-convergence.R`, `R/runOccJSDM.R:1447`
- Test: `tests/testthat/test-grouped-convergence.R`

**Interfaces:**
- Consumes: `grouped_convergence()`.
- Produces: `computeDiagnostics(results_output, embedding_rows = NULL)`. With `embedding_rows` (a logical per row of `B_output`) and a `lambda_output` present: `B_output` is reported as two sets, `B (measured)` with its warnings as before and `B (embedding)` with its R-hat and ESS lines, no warning, and the line `  [Not a test of convergence under the grouped prior: see returnGroupedConvergence()]`; `sigmab_output`, `sigmaemb_output` and `lambda_output` are reported on the log scale under the names `log sigmab`, `log sigmaemb`, `log lambda`. Otherwise unchanged.
- Produces: `print_grouped_convergence(x)`, which prints with `message()` a block headed `Grouped prior: convergence checks`, one line per check, the verdict, and for each flagged check a line starting `  ! ` with its meaning.
- `runOccJSDM()`: passes `embedding_rows` when the fit is grouped with embedding columns, then calls `grouped_convergence()` and `print_grouped_convergence()`, and gives one warning, `"The grouped prior's convergence checks flag this fit for review: see returnGroupedConvergence()."`, when the verdict is "review".

- [ ] **Step 1: Write the failing tests.**
  - `computeDiagnostics()` on a hand-built `results_output` (two chains; `B_output` with 2 measured and 3 embedding rows; one embedding coefficient shifted by 5 in chain 2; `lambda_output`, `sigmaemb_output`, `sigmab_output`, `B0_output` well mixed): without the argument it warns with "Convergence issue in 'B_output'"; with `embedding_rows = c(FALSE, FALSE, TRUE, TRUE, TRUE)` it gives no warning, and its messages contain "B (measured)", "B (embedding)", "Not a test of convergence" and "log lambda".
  - The same with a measured coefficient shifted: with the argument it still warns, and the warning names the measured set.
  - `print_grouped_convergence()` on a flagged result prints a line starting "  ! " and "Verdict"; on an unflagged one, no "  ! " line.
  - A grouped binary fit's messages contain "Grouped prior: convergence checks" once; a fit without `embeddingVars` contains none.

- [ ] **Step 2: Run** `devtest.sh grouped-convergence`. Expected: the new tests fail.

- [ ] **Step 3: Implement.** In `computeDiagnostics()` move the body of the loop into a local function taking the array, the name to print, the name used in the warnings, and whether to warn, and call it once per block as now, or twice for `B_output` and on the log for the three scales when `embedding_rows` is given. Messages and warnings of the unchanged path must be character-for-character those of today.

- [ ] **Step 4: Run** `devtest.sh grouped-convergence`, then `devtest.sh diagnostics` (the existing tests of `computeDiagnostics()`, which compare it with the legacy implementation). Expected: all pass.

- [ ] **Step 5: Commit.**

### Task 4: Documentation

**Files:**
- Modify: `R/coef-prior.R` (the grouped-prior help page), `R/diagnostics.R` (the `computeDiagnostics()` help), `NEWS.md`
- Regenerate: `man/grouped-prior.Rd`, `man/computeDiagnostics.Rd`

- [ ] **Step 1:** In "Sampling and chain length", replace the sentences from "Before relying on a fit, compare the chains" to "runs high on such draws." with: what `returnGroupedConvergence()` checks and that `runOccJSDM()` prints it; that it reports "review", not failure, and why; that its checks of the gradients compare chain means and do not say whether each chain explored its own uncertainty; that single embedding coefficients are reported without a warning. In "What the fit holds", add the function to the list of summaries.
- [ ] **Step 2:** Document `embedding_rows` in `computeDiagnostics()`.
- [ ] **Step 3:** Add one entry to the development section of `NEWS.md`.
- [ ] **Step 4:** Regenerate the two Rd files in the scratch build, check that no other Rd changed, copy them over, and render the grouped-prior page to text to read it.
- [ ] **Step 5:** Run the whole suite: `devtest.sh`. Expected: 0 failures, no new skips.
- [ ] **Step 6: Commit.**

### Task 5: The check on the saved fits

**Files:**
- Create: `dev/simstudy/grouped-mixing-20261118/checks/convergence_check_on_saved_fits.R`
- Modify: `dev/simstudy/grouped-mixing-20261118/README.md`
- Create: `dev/simstudy/grouped-mixing-20261118/results/convergence_check_on_saved_fits.csv`

- [ ] **Step 1:** The script takes a library as its first argument, reads each `_summary.rds` of `results/rotation` and `results/orientation`, and calls `occJSDM:::grouped_effect_checks()` on its `contribution_mean`. It writes one row per fit: the number of gradients, the agreement, the disagreement, the size, and which of the two would be flagged.
- [ ] **Step 2:** Run it with the installed prototype. Expected, from the measurement made before this plan: agreement 0.93 to 0.99 in 16 fits and 0.87 and 0.80 in the two fits known to have chains apart or in the leftover state; disagreement over 0.15 in those two and in one fit of seed 3 (0.156).
- [ ] **Step 3:** Record in the README, under a new dated heading: what was built, the table's figures, that the two limits were set against these same fits, and what the function does not check.
- [ ] **Step 4: Commit.**

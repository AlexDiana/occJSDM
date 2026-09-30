# Occupancy-intercept prior widening: implementation and comparison plan

> **For agentic workers:** REQUIRED SUB-SKILL: use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Test whether replacing the hard-coded Normal(0, SD 1) prior on the species occupancy intercept `B0` with a wider prior reduces occupancy-probability bias (rare-species overestimation, and low/high probabilities pulled toward the middle) without harming interval coverage or convergence, then decide whether to change the default.

**Architecture:** Add an opt-in `listPriors$sigma_b0` (prior SD of `B0`, default 1, so defaults are bit-for-bit unchanged) threaded through every place that encodes the `B0` prior precision: the two C++ samplers, the R reference sampler and the two spatial-range log-weight functions. Then run paired full-model fits on the saved datasets, reusing the existing default-prior fits as the control arm, and score with the existing study scorers.

**Tech Stack:** R, Rcpp/Armadillo, testthat, the existing `dev/simstudy` scoring and audit scripts.

**Spec:** [TODO.md](../../../TODO.md), "Required bias recheck and release preparation", item 2, and the release criterion in the same section (retain all features, fix materially biased point estimates, coverage may wait). Background: [CLAUDE-HANDOVER.md](../../../CLAUDE-HANDOVER.md), "Completed diagnosis" (conditional experiment C, which shows estimating the intercept under the current prior gives +3.62 points signed error for 1%-occupancy species, and explicitly does not test a wider prior).

## Global Constraints

- Default priors and default output must be unchanged. `sigma_b0 = 1` must reproduce the current fits within 1e-12 (cross-installation) and identically within one installation.
- Widening is a candidate remedy, not an established fix. The default changes only if the prespecified gate below passes, and then only in a separate reviewed PR.
- Retain informative detection priors (`a_p/b_p`, `a_q/b_q`, `a_theta0/b_theta0` unchanged). Do not touch `b_betatheta_slope_var`, the spatial-amplitude prior or `tau_prior`.
- Never edit frozen archives, inputs or existing fits. Write new fits, hashes and logs to a separate ignored archive.
- At most eight concurrent single-threaded fits (Doug's standing authorisation on this 10-core, 32 GB machine). Launch long runs detached with `nohup` and a waiter script, per the memory note on detached simulation jobs.
- Markdown rules for this repo: one line per paragraph, no pipe tables, escape `-\>` and `\~` in prose, no em-dashes.
- Work on a new branch `codex/occupancy-intercept-prior` from main `7e6568a`. Do not commit to main.

## Design decisions to confirm with Doug before Task 3 (protocol freeze)

These are choices, not facts. The recommendations are what the plan below assumes.

- **Candidate prior SDs.** Recommend 2, 3 and 5 against the control SD 1. Rationale: logit(0.01) is about -4.6 and logit(0.99) about +4.6, so SD 1 puts almost no mass on rare or ubiquitous species, SD 2 covers them at about 2.3 SDs, SD 3 comfortably, and SD 5 is close to flat on the probability scale's useful range. The current mean is 0 and stays 0; a hierarchical prior on the intercept scale is a separate, larger change and is out of scope.
- **Confounding with the collection intercept.** In occupancy and two-stage fits `B0` and the `beta_theta` intercept are only separated by the data at low replication (see PLAN.md section 16 in `dev/simstudy`). Widening `B0` alone may shift bias into `theta`. Recommend reporting `beta_theta` intercept bias and the `B0` versus `theta` posterior correlation as secondary outcomes rather than widening both.
- **Datasets.** Recommend phase A on binary fits (isolates the prior from detection error) and phase B on two-stage fits (the real use case). Both reuse saved inputs and RNG states.

## Data and arms

Reuse, without modification, the inputs and current-main control fits from these archives (paths in `dev/simstudy/current-main-recheck/README.md` and `dev/simstudy/spatial-amplitude-prior/README.md`):

- **Phase A1, non-spatial binary JSDM:** ten communities at each of 100 and 300 sites (the 1,000-site cells are dropped: the prior matters least there, and they are the most expensive).
- **Phase A2, spatial binary:** the nine saved communities at 100 supports (three ranges, three replicates), inverse-gamma spatial amplitude, the longer four-chain schedule where the control used it.
- **Phase B, two-stage:** the ten communities at 300 sites for each of the two contamination levels of the K6 sampling-design study, run only for the arms that pass the phase A gate.

Arms: control `sigma_b0 = 1` (existing fits, verified reproducible in Task 2), then `sigma_b0` in 2, 3 and 5. Initial schedule two chains, 3,000 burn-in and 5,000 retained iterations. Flagged fits repeat once on four chains, 6,000 burn-in and 12,000 retained, with selection rules taken verbatim from the current-main recheck protocol (Rhat above 1.05, ESS below 100, non-finite diagnostics, or a fitting warning), and applied identically to all arms before any error is compared.

## Outcomes and gate (fix before any fit is scored)

Primary outcome: mean signed occupancy-probability error in percentage points, computed against generating truth, in three truth bands (below 20%, 20-80%, above 80%) plus a separate rare-species group (species with about 1% occupancy), with communities weighted equally. This is the same estimand as the agreed provisional spatial target (within five points per band).

Secondary outcomes: occupancy MAE, 95% interval coverage of occupancy probabilities, `B0` and `beta_theta` intercept bias and coverage, convergence flags per arm, and for phase B the `B0`/`theta` posterior correlation.

Gate for recommending a new default, all required:

- Mean absolute signed band error is lower than the control in the below-20% band and the rare-species group, in at least two thirds of communities.
- No band's mean absolute signed error worsens by more than 1 point, and overall occupancy MAE does not worsen by more than 0.5 points.
- Coverage of occupancy probabilities does not fall by more than 0.03 in any band (coverage may otherwise remain below nominal, per the release criterion).
- Convergence flags per arm are no more numerous than control after the single longer repeat.

If several SDs pass, pick the smallest passing SD. If none pass, keep the default, record the negative result, and document the residual bias as a beta limitation (TODO release item 5).

## File structure

- Modify `src/jsdm.cpp` (`sample_BBsL_cpp` around line 1063 and `sample_BBsL_parallel` around line 1276): accept `double sigma_b0` and use it for the intercept element of `diag_B`.
- Modify `R/RcppExports.R` and `src/RcppExports.cpp`: regenerated by `Rcpp::compileAttributes()`, never hand-edited.
- Modify `R/jsdmfun.R`: the R reference `sample_BBsL` (line 1060), both `spatial_range_logweights*` functions (lines 1229 and 1280), the caller near line 1435, and the parameter list at line 698.
- Modify `R/runOccJSDM.R`: add `read_intercept_prior()`, store the value in `list_jSDMparams`, record it in `infos$intercept_prior`, and document `sigma_b0` in the roxygen for `listPriors`. Regenerate `man/runOccJSDM.Rd` with `devtools::document()`.
- Create `tests/testthat/test-intercept-prior.R`: unit and integration tests.
- Create `dev/simstudy/occupancy-intercept-prior/` containing this `PLAN.md`, `run.R`, `analysis.R`, `summarise.R`, `README.md`. Reuse the scorers in `dev/simstudy/spatial-amplitude-prior/metrics.R` and `dev/simstudy/current-main-recheck/summarise.R` by sourcing them, not copying.

## Task 1: Thread `sigma_b0` through every place the prior is encoded

**Files:**
- Modify: `src/jsdm.cpp`, `R/jsdmfun.R`, `R/runOccJSDM.R`
- Test: `tests/testthat/test-intercept-prior.R`

**Interfaces:**
- Produces: `listPriors$sigma_b0` (finite positive scalar, default 1); `list_jSDMparams$sigma_b0`; `sample_BBsL_cpp(..., sigma_b0, ...)` and `sample_BBsL_parallel(..., sigma_b0, ...)` with the new argument placed directly after `sigma_b`; `infos$intercept_prior` equal to `list(mean = 0, sd = sigma_b0)`.

The intercept prior precision is currently the literal `1` at four sites: `diag_B` ones-fill in both C++ samplers, `diag(B_current)` initial value in the R sampler, and `c(1, rep(1/sigma_b^2,p), ...)` in both range log-weight functions. Missing the last two would make the spatial range update marginalise `B0` under a different prior from the one used to sample it, a silent inconsistency, so the tests below target it explicitly.

- [ ] **Step 1: Write failing tests** in `tests/testthat/test-intercept-prior.R`:

```r
test_that("read_intercept_prior validates and defaults", {
  expect_equal(occJSDM:::read_intercept_prior(list()), list(mean = 0, sd = 1))
  expect_equal(occJSDM:::read_intercept_prior(list(sigma_b0 = 3))$sd, 3)
  for (bad in list(0, -1, NA_real_, Inf, c(1, 2), "a")) {
    expect_error(occJSDM:::read_intercept_prior(list(sigma_b0 = bad)),
                 "sigma_b0 must be a finite positive number")
  }
})

test_that("B0 conditional draws match the analytic posterior under a wide prior", {
  # One species, intercept only, Gaussian likelihood surrogate via Polya-Gamma
  # weights fixed at 1: posterior precision = sum(omega) + 1/sigma_b0^2,
  # mean = sum(k) / precision. Compare the sampler mean and variance.
  n <- 40L
  k <- rep(0.5, n)
  omega <- matrix(1, n, 1)
  X <- matrix(0, n, 0); Tr <- matrix(0, 1, 0); U <- matrix(0, n, 0)
  zero <- function(r, c) matrix(0, r, c)
  for (s0 in c(1, 3)) {
    setOccJSDMSeed(11)
    draws <- replicate(4000, occJSDM:::sample_BBsL_cpp(
      matrix(k, n, 1), X, Tr, U, zero(0, 0), zero(1, 0), zero(0, 0), 1,
      zero(0, 0), zero(1, 0), zero(0, 0), 1, zero(n, 0), zero(0, 0),
      omega, "binary")$B0)
    prec <- n + 1 / s0^2
    expect_equal(mean(draws), sum(k) / prec, tolerance = 0.05)
    expect_equal(var(draws), 1 / prec, tolerance = 0.1)
  }
})

test_that("spatial range log-weights use the same B0 prior as the sampler", {
  # Build a tiny fixture, call spatial_range_logweights with sigma_b0 = 1 and 3,
  # and check the difference equals the analytic change in the marginal
  # likelihood of an intercept-only Gaussian block. Write the expected value
  # from the fixture once the signature is confirmed in Step 3.
  skip("fill from fixture in Step 3")
})

test_that("default fit is identical with and without an explicit sigma_b0 = 1", {
  sim <- simulate_fixture(model = "binary")
  run <- function(priors) {
    set.seed(41)
    suppressMessages(suppressWarnings(runOccJSDM(sim$data_list,
      occCovariates = fixture_occ_covariates(),
      listParams = list(n_factors = 0L, n_lattrait = 0L),
      listPriors = priors,
      MCMCparams = list(nchain = 2L, nburn = 12L, niter = 15L, nthin = 1L))))
  }
  expect_identical(run(list(sigma_b0 = 1))$results_output, run(list())$results_output)
  expect_identical(run(list())$infos$intercept_prior, list(mean = 0, sd = 1))
  expect_false(identical(run(list(sigma_b0 = 4))$results_output,
                         run(list())$results_output))
})
```

The third test's `skip` is the only deferred piece and is resolved in Step 3, not left in the committed file.

- [ ] **Step 2: Run them to confirm they fail**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-intercept-prior.R")'`

Expected: FAIL (`read_intercept_prior` not found, unused argument `sigma_b0`).

- [ ] **Step 3: Implement.** In `src/jsdm.cpp`, add `double sigma_b0,` after `double sigma_b,` in both signatures and replace the fill in both bodies:

```cpp
arma::vec diag_B(total_dim, arma::fill::ones);
diag_B(0) = std::pow(sigma_b0, 2);
```

In `R/jsdmfun.R` add the same argument to `sample_BBsL` (setting `diag(B_current)[1] <- sigma_b0^2`), to the two range-log-weight functions (replacing the leading `1` in `prior_precision` and `precision` with `1/sigma_b0^2`), and pass `list_params$sigma_b0` at the caller near line 1435 and in the parameter plumbing near line 698. Write the real expected value for the range log-weight test now, using the fixture the existing `tests/testthat/test-spatial-range.R` builds, and delete the `skip`. In `R/runOccJSDM.R` add:

```r
# The occupancy intercept prior is Normal(0, sd^2); it is not the collection intercept.
read_intercept_prior <- function(priors) {
  sd <- get_param(priors, "sigma_b0", 1)
  if (!is.numeric(sd) || length(sd) != 1L || !is.finite(sd) || sd <= 0)
    stop("sigma_b0 must be a finite positive number")
  list(mean = 0, sd = unname(sd))
}
```

Call it beside `read_spatial_sd_prior`, put `sigma_b0 = intercept_prior$sd` in the initial `list_jSDMparams`, set `infos$intercept_prior`, and add the roxygen paragraph. Then run `Rscript -e 'Rcpp::compileAttributes(); devtools::document()'`.

- [ ] **Step 4: Run the new tests, then the full suite**

Run: `Rscript -e 'devtools::load_all(); testthat::test_file("tests/testthat/test-intercept-prior.R")'` then `Rscript -e 'devtools::test()'`

Expected: new tests PASS; full suite unchanged from the last recorded count (858 passing, three expected skips) plus the new expectations.

- [ ] **Step 5: Independent review, then commit.** Have a fresh reviewer check that every site holding a literal `1` for the intercept precision was changed (`grep -n "c(1, rep\|fill::ones\|B_current" R/jsdmfun.R src/jsdm.cpp`).

```bash
git add src R man tests
git commit -m "Add opt-in occupancy-intercept prior SD sigma_b0, default unchanged"
```

## Task 2: Cross-installation default equivalence and pilot

**Files:**
- Create: `dev/simstudy/occupancy-intercept-prior/run.R`, `equivalence.R`

**Interfaces:**
- Consumes: `runOccJSDM(..., listPriors = list(sigma_b0 = <sd>))` from Task 1.
- Produces: `run.R --repo=REPO --study=STUDY --reference=REFERENCE --phase=A1|A2|B --sd=1|2|3|5 --workers=N`, writing one complete fit RDS plus a result CSV per community into `STUDY/fits/<phase>/sd<sd>/`.

- [ ] **Step 1:** Install the frozen control revision (`7e6568a`) and the new branch into two separate libraries. Copy `dev/simstudy/spatial-amplitude-prior/equivalence.R` and adapt it: short fit on one saved community with the saved RNG state, default priors, compare posterior values, warnings and final RNG state.
- [ ] **Step 2:** Run it. Expected: maximum absolute posterior difference at most 1e-12, identical warnings, identical final RNG state. Any larger difference blocks the study.
- [ ] **Step 3:** Copy `spatial-amplitude-prior/run.R` to `run.R`, replacing the prior switch with `listPriors$sigma_b0`, keeping input hashing, per-fit fingerprints, one sampler thread and the detached-launch pattern. Run one pilot fit at `--sd=3` on one community per phase; verify `infos$intercept_prior$sd == 3`, that inputs are unchanged by hash, and that the control result CSV for that community reproduces from the reused fit.
- [ ] **Step 4:** Commit the scripts and the equivalence result.

## Task 3: Freeze the protocol and run phase A

**Files:**
- Create: `dev/simstudy/occupancy-intercept-prior/README.md` (prespecified protocol, the Outcomes and gate section above, copied verbatim)

- [ ] **Step 1:** Confirm the three design decisions above with Doug, then commit the frozen protocol and record the commit hash. No fitting begins before this hash exists, and nothing above is edited afterwards except by a dated amendment file.
- [ ] **Step 2:** Launch phase A1 and A2 for `--sd=2,3,5` detached, at most eight workers, with a waiter script that writes a completion marker (do not poll from the session).
- [ ] **Step 3:** When complete, apply the selection rule identically to all arms, launch the single longer repeat for flagged fits, and record the selection list before opening any error table.

## Task 4: Score phase A, audit, decide

**Files:**
- Create: `dev/simstudy/occupancy-intercept-prior/analysis.R`, `summarise.R`, `test-analysis.R`

**Interfaces:**
- Consumes: fits from Task 3 and the control results from the archived studies.
- Produces: `STUDY/summary/band-error.csv` (arm by band by community, signed error in points), `coverage.csv`, `convergence.csv`, `gate.csv` (one row per SD, one column per gate criterion, plus overall PASS or FAIL).

- [ ] **Step 1:** Write `test-analysis.R` first, with a hand-computed three-species example for each metric (signed error by band, rare group, coverage), and watch it fail.
- [ ] **Step 2:** Implement `analysis.R` by sourcing the existing scorers; add only the gate logic. Run the tests to green.
- [ ] **Step 3:** Independently audit every summary from the saved draws, as the spatial-amplitude study did (`verify.R`), and require exact reproduction of stored diagnostics.
- [ ] **Step 4:** Read `gate.csv`. If no SD passes, stop: record the negative result and skip phase B. If an SD passes, proceed to Task 5 for that SD (and any smaller passing SD) only.
- [ ] **Step 5:** Commit the results and a plain-language report.

## Task 5: Phase B (two-stage), conditional on the gate

- [ ] **Step 1:** Run the two-stage fits for the surviving SDs on the ten communities per contamination level, same schedule and selection rules.
- [ ] **Step 2:** Score with the same gate plus the `B0`/`theta` secondary outcomes. A pass in phase A that fails here is reported as a binary-only improvement and does not change the default.
- [ ] **Step 3:** Commit results.

## Task 6: Record and hand off

- [ ] **Step 1:** Update TODO.md item 2 (mark done, state the outcome and link the report), and AGENTS.md minimally, observing the markdown rules. Re-read each file immediately before and after editing.
- [ ] **Step 2:** If the gate passed in both phases, open a separate small PR that changes only the default of `sigma_b0`, regenerates `sampleresults`, and rebuilds affected vignettes (this joins TODO release item 4). Otherwise open the PR containing only the opt-in option and the negative result. Do not merge or change defaults without Doug's approval.

## Self-review

- Spec coverage: compares the current prior with wider alternatives (Tasks 1 to 3), in paired full-model fits on saved datasets (Data and arms), assesses bias, coverage and convergence before choosing a default (Outcomes and gate, Task 4), retains informative detection priors (Global Constraints), and treats widening as a candidate rather than a fix (gate, Task 4 Step 4).
- No placeholders remain except the range log-weight expected value, which Task 1 Step 3 requires the implementer to compute from the existing spatial-range test fixture.
- Names are consistent: `sigma_b0`, `read_intercept_prior`, `infos$intercept_prior`.

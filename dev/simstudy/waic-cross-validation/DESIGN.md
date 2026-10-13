# Site cross-validation of the Lesson 4 factor-count comparison: step 1

Written 13 October 2026, before any study fit. The decision rule below is fixed before the results are seen.

## Question

Lesson 4 compares a one-factor and a two-factor model of the same survey. The observed-data site WAIC (`extractWAIC()`, PR #14) gives a difference of 1.31 on a score of about 23,000 (paired SE 0.76), and flags 99 of 100 sites as unreliable. WAIC is a shortcut for leave-one-site-out cross-validation. This study computes cross-validation directly, by refitting, to answer three things: (1) whether the WAIC difference agrees with it, and so whether users can be told to rank models by WAIC; (2) whether one or two factors predict new sites better on Lesson 4's survey, where "no meaningful difference" is a valid answer; (3) refitting and scoring code that could become a user-facing site-level cross-validation workflow. It does not validate WAIC in general, find a "true" factor count, or cover spatial models.

## Data and models

- Data: Lesson 4's survey, `input$sim$data_list` in `vignettes/teaching-data/nonspatial-lesson.rds` (identical to the archived `input.rds`): 100 sites, 10 species, three field samples per site, two primers, six PCRs per primer and sample, read threshold 1.
- Models: the archived Lesson 4 arguments (`helpers.R` and `prediction-build.R`), differing only in `n_factors` (1 or 2): one latent trait, occupancy covariates `X_psi.EnvCov.1` and `X_psi.EnvCov.2`, collection covariate `X_theta`, no spatial term, default priors, summarised latent presences.
- MCMC: four chains, 3,000 burn-in and 6,000 retained iterations (24,000 draws), as in the archived fits.
- Software: R 4.6.1 and occJSDM at the worktree's `HEAD` (`c102745`), installed into a dedicated study library from a `git archive` export. The package source has changed since the archived fits (`a56f554`), so both full-data models are refitted too, and WAIC and cross-validation come from identical software.
- Threads: every fit runs with `RCPP_PARALLEL_NUM_THREADS=1`, the only setting at which fits are bit-reproducible (AGENTS.md), and several fits run side by side.

## Design

- Folds: the 100 sites are shuffled with seed 20261013 and split into ten folds of ten. Both models use the same folds, so the comparison is paired by site.
- Fits: 22 in total. Two full-data fits, with the archived seeds (one factor 20260926, two factors 20260921), give WAIC. Twenty fold fits, each on the other 90 sites, use seed `2026101300 + 10 * fold + m`, where m is 1 for one factor and 2 for two.
- Held-out score: for each withheld site, the log posterior predictive density of its complete observed survey under the fold's fit, with occupancy and collection states summed out and the two (or one) shared site factors integrated out. This is exactly `computeSiteWAIC()`'s per-draw likelihood, applied to the withheld sites' data prepared as `runOccJSDM()` prepares data, with covariates standardized by the training fit's constants. All 24,000 draws are used. The scoring route is tested to reproduce the fitted design and `computeSiteWAIC()` exactly on training sites, and to match a brute-force likelihood from the raw PCR rows on withheld sites (`test-cv-lib.R`).
- WAIC: `computeSiteWAIC()` on each full-data refit, all 24,000 draws.
- Numerical checks: each score is recomputed on 1,000 evenly spaced draws with stricter quadrature (orders 31, 51, 81, 121; tolerance 1e-6); site log likelihoods must agree within 1e-4.
- Convergence: the project's strict screen on every refit (nonfinite Rhat, Rhat above 1.05, or bulk or tail ESS below 400), on the occupancy, collection and detection parameters and on the fitted sites' occupancy probabilities, which are identified although the factors are not. Reported for every fit. The screen does not by itself decide the outcome; the target-level check below does.

## Quantities

All in elpd units (log predictive density summed over sites, higher is better), with deviance units (`-2 *` elpd) alongside to connect with the WAIC numbers Lesson 4 quotes. Differences are one factor minus two factors, so a positive difference favours one factor.

- Cross-validation: elpd_cv per model, the paired difference, and its across-site SE, `sqrt(100 * var(d_i))`.
- Monte Carlo SE of the cross-validation difference: per fold, a delta-method MCSE of the fold's summed held-out lppd from the four chains; combined across the 20 independent fold fits.
- WAIC: elpd_waic per model from the full refits, the paired difference and its SE.
- Per site: correlation between elpd_cv and elpd_waic, their mean difference (WAIC's optimism), and the correlation of the paired per-site differences.

## Decision rule

- Each method's verdict: "clear preference" for the model with higher elpd when the absolute difference is at least 4 and at least twice its paired SE; otherwise "no meaningful difference". The threshold of 4 follows the cross-validation FAQ of Vehtari and colleagues (Sivula, Magnusson and Vehtari 2022): smaller elpd differences are small regardless of their SE.
- Agreement: WAIC agrees with cross-validation when the two verdicts are the same, including the direction of any preference.
- Step 1 is unclear if the cross-validation verdict changes when its difference is moved by two Monte Carlo SEs in either direction, or when it is recomputed from any single chain's draws alone (four separate estimates).
- Step 2 (full leave-one-site-out, 200 fits) is run only if step 1 is unclear or WAIC disagrees with it. Otherwise step 1 answers questions (1) and (2) for this survey.

## What a result can and cannot support

Agreement on one survey shows that WAIC gave the same answer as direct cross-validation here, not that it is reliable in general; the 99 of 100 flagged sites remain a reason for caution. Ten-fold cross-validation trains on 90 sites rather than 99, which biases both models' scores slightly downwards; the paired difference is much less affected. A "no meaningful difference" verdict is not evidence for either factor count, only that new-site predictions do not separate them on this survey.

## Files

- `cv-lib.R`: folds, subsetting, fitting, the held-out scoring copy, summaries and the convergence screen.
- `test-cv-lib.R`: tests, run before launch.
- `prepare.R`, `fit-job.R`, `score-job.R`, `summarise.R`, `launch.sh`: the run. Raw fits, scores and logs go to `dev/simstudy/results/waic-cv-20261013/` (ignored); the compact results come back here.

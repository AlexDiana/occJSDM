# Stage 2: spatial representation comparison on the recoverability pilot's communities

Date: 11 October 2026. Status: design agreed in conversation with Doug on 11 October; written for his review before an implementation plan is made. Nothing has been implemented or run. The design follows from the spatial recoverability pilot (`dev/superpowers/specs/2026-10-10-spatial-recoverability-pilot-design.md`, report `dev/simstudy/spatial-recoverability-pilot/REPORT.md`) and reuses infrastructure from the first spatial package comparison (branch `codex/spatial-comparison-plan`, design `dev/superpowers/specs/2026-10-08-spatial-package-comparison-design.md`).

## Purpose

The pilot established, with an oracle that knew every parameter but the fields, where spatial structure is recoverable from single-season binary occupancy: common gradients (truth 1) from 300 sites at range 0.08 and from 100 at range 0.15; per-species dispersal fields (truth 2) from 300 sites for the long class and 900 for the medium class. Its map-reading rule recommended 900 sites as the positive-control setting, with 300 as the harder secondary.

Stage 2 fits models that must estimate what the oracle was told, and asks which existing spatial representation recovers the most of the oracle's held-out headroom under each truth, and how much each loses when the truth is the other kind. Its output is a development recommendation between two representations, shared spatial factors with species loadings and per-species fields, and a measured target for a later trait-linked prototype. Doug's stated preference for the eventual model is a continuous trait controlling each species' range; that prototype is deferred to a stage 3 so that stage 2 can be run with code that already exists.

Decisions Doug made on 11 October: compare existing representations only; fit four pilot cells; five arms; the smallest meaningful difference between arms is 0.10 of the oracle's headroom.

## Non-goals

- No new occJSDM sampler, prior, default or API change. The trait-linked per-species range prototype is stage 3.
- No new communities, seeds, ranges, arrangements or site counts. Only the pilot's inputs are used.
- No imperfect detection and no eDNA replicate structure.
- No mixed truth.
- No sjSDM arm and no fixed-basis arm (trend, Moran eigenvectors). The first comparison found the coordinate network unfit under its acceptance screen and the fixed bases worse than environment-only in every package.
- No tuning of any setting against the generating range or the lattice truth.
- No change to the pilot's thresholds, scores or results.

## Hypotheses

1. Under truth 1, the factor arms capture more of the oracle headroom than occJSDM's per-species GP, because the truth is two shared fields and the GP must estimate twelve fields with one range and one amplitude.
2. Under truth 2, occJSDM's GP captures more than the two-factor arms, because twelve independent fields cannot be represented by two shared ones, and the six-factor arm narrows but does not close that gap.
3. Every arm falls well short of the oracle in every cell, because the oracle knows intercepts, slopes, range, amplitude and, under truth 1, loadings. The size of that shortfall is the development target.

The pairwise comparisons in "Decision rule" are fixed so that each hypothesis has a reading whichever way it falls.

## Datasets

- The pilot's communities, unchanged, from `dev/simstudy/results/spatial-recoverability-20261010/inputs/` with the MD5s in its `input-manifest.csv`: cells `t1-r0.08-n300`, `t1-r0.08-n900`, `t2-n300` and `t2-n900`, replicates 1 to 10, random-uniform arrangement only. 40 datasets.
- Each dataset supplies the training occupancy states, the measured environment and the raw coordinates of its training sites, and the same 40 by 40 held-out lattice of 1,600 points with its true fields and probabilities. Training data carry no traits.
- The oracle reference for every dataset is the pilot's final job file: the resolution refit in `jobs-resolution/` where one exists, else the re-run file in `jobs/`. No oracle is fitted in stage 2. The oracle's screen status is carried into the scores.
- No community is discarded or replaced after any fit or score, and no favourable community is selected.

## Arms

Five arms are fitted on every dataset. Every arm receives the same training occurrences, species identifiers and measured environment; coordinates go only to spatial arms. Configurations follow the first comparison's adapters where they exist; the two settings marked "frozen by the timing pilot" are chosen from cost alone, before any comparison fit, and recorded in the run config.

- **env**: occJSDM binary model, `occCovariates = "environment"`, no `spatCovariates`, `n_factors = 0`, `n_lattrait = 0`, default priors. The spatially uninformed baseline that every spatial arm is measured against.
- **gp**: occJSDM binary model with `spatCovariates = c("longitude", "latitude")`, `n_factors = 0`, `n_lattrait = 0`, default inverse-gamma amplitude prior and default range grid. Per-species fields with one shared range and one shared amplitude. Support points: one per site at n = 300 (300); at n = 900 the count is frozen by the timing pilot, 900 if affordable, otherwise the largest of 600, 450 and 300 that is, with the realised count recorded. Prediction at new sites through the package's spatial predictor, as certified in the first comparison.
- **hmsc2**: Hmsc probit, one spatial random level at site, `nfMin = nfMax = 2`, raw coordinates. The truth 1 generating dimension.
- **hmsc6**: as hmsc2 with `nfMin = nfMax = 6`. Room for a factor model to approximate independent fields.
- **gllvm2**: gllvm binomial logit, variational approximation, two spatial latent variables with exponential correlation, `num.lv` otherwise zero, no additional row effects. A second implementation of the two-factor representation, so that a factor-arm finding can be read as representation rather than fitting method.

Hmsc's spatial method is frozen by the timing pilot: `sMethod = "Full"` if affordable, otherwise the predictive-process approximation with a knot count recorded in the config. Hmsc's realised default range grid and gllvm's range parameterisation are recorded as in the first comparison. Hmsc's two or six exponentially correlated factors differ from the squared-exponential generating kernel; this is part of what the representation comparison measures and is not corrected for.

## Prediction targets and scores

The common target is the marginal latent occurrence probability at every lattice point, as in the first comparison: means are averages of draw-level probabilities, never an inverse link of a mean linear predictor. For each selected fit, two versions: the full prediction, and the spatially constant prediction with the draw-level training mean spatial effect moved into the baseline and the spatial variation removed. The arm's contrast is full minus spatially constant. Hmsc uses its certified conditional latent-field prediction; gllvm uses the first comparison's certified conditional-moment adapter; occJSDM uses its own spatial predictor. Where a certified new-site prediction is unavailable for a fit, the cell is marked unsupported and stays in every denominator.

Per species, per dataset and arm, on the held-out lattice, against the true probability and the true field:

- Probability RMSE and MAE of the full prediction.
- The headroom share, the primary metric: `H = (RMSE_env - RMSE_arm) / (RMSE_env - RMSE_oracle)`, where every RMSE is the held-out probability RMSE against the true probability for that species and dataset, `RMSE_env` is the fitted env arm's, and `RMSE_oracle` is the pilot oracle's. H is 1 when the arm matches the oracle, 0 when it adds nothing over the fitted environment, and negative when it is worse. Because the oracle also knows intercepts and slopes, H above 1 is possible only through noise.
- Pattern: the correlation between the arm's contrast and the true contrast, where the true contrast is the true probability minus the pilot's spatially uninformed reference. This is the pilot's pattern score applied to a fitted arm.
- Attenuation slope and correlation of the arm's centred species logit field on the true centred field, for logit-scale arms only (gp, gllvm2). Hmsc's probit fields are not scored against logit truth.
- Training-site versions of the same quantities, reported and never used in decisions.

Aggregation follows the pilot: average species within occupancy level, pool the 25, 50 and 75 percent levels, then one value per dataset and arm. The 5 percent species are fitted and reported and are not part of any decision.

## Diagnostics and finite budgets

- Bayesian arms (gp, hmsc2, hmsc6): four chains, 3,000 warm-up and 4,000 retained draws per chain, with one prescribed extension to 8,000 and 16,000 when the screen fails. Screen: rank-normalised Rhat at most 1.01 and bulk and tail ESS at least 400 for every scored species probability mean at the training sites, the species-level average of the held-out probability, the training field average per species, and the spatial SD and range (gp) or the reconstructed species fields' averages (Hmsc), or any native convergence warning. No quantity that uses the truth enters a fitted arm's screen. A fit that fails after extension is `unresolved`.
- gllvm2: three starts, with three more only when the best objective is not replicated within tolerance; a fit whose best solution is not replicated is `unstable`. The env arm under occJSDM follows the Bayesian schedule.
- Statuses are `pass`, `unresolved`, `unstable`, `unsupported` and `no fit`. Every status keeps its best prescribed numerical output for descriptive scoring with the status attached, except `no fit`, which has no fabricated prediction. All statuses remain in every denominator. No diagnostic threshold or attempt ceiling is relaxed after results are seen.
- Attempt ceiling: 200 comparison fits, at most two attempts each (initial plus extension or the extra starts), so 400 attempts.

The first comparison found all twelve occJSDM GP fits unresolved at 100 sites. If the gp arm ends unresolved in a stage 2 cell, that is reported as a finding about the implementation at that scale, with its descriptive scores, and not promoted or dropped.

## Decision rule, fixed before any fit

Within each of the four cells, paired differences in the headroom share across the ten communities, with a 95 percent t interval, for these comparisons: gp against hmsc2; gp against hmsc6; gp against gllvm2; hmsc6 against hmsc2; hmsc2 against gllvm2. A community whose two fits do not both have status `pass` is kept in the descriptive table and excluded from the pre-specified interval, and the number excluded is reported with the interval; an interval from fewer than six communities is reported as inconclusive.

A difference counts as clear only if its estimate is at least 0.10 in absolute value and its interval excludes zero. "Not clearly worse" means the estimate is above minus 0.10 or the interval includes zero.

Development reading, applied to the 900-site cells first and then checked at 300:

- Factor arms clearly better than gp under truth 1, and not clearly worse under truth 2: develop shared spatial factors with species loadings in occJSDM.
- gp not clearly worse than the factor arms under truth 1, and clearly better under truth 2: develop per-species ranges, and the trait-linked prototype is stage 3.
- Factor arms clearly better under truth 1 and gp clearly better under truth 2: the representations trade off, and the development choice depends on which truth Doug expects in the eDNA data; the report states both shortfalls and makes no recommendation.
- No clear difference anywhere, and every arm below 0.5 of the oracle headroom in every cell: representation is not the bottleneck at these scales, and the next study is about inference, not representation.
- hmsc6 clearly better than hmsc2 under truth 2 is read as factor dimension, not representation, driving the two-factor shortfall; hmsc2 clearly different from gllvm2 in either direction is read as a fitting-method effect and qualifies every factor-arm conclusion.

The absolute level of the headroom share is reported for every arm and cell, with between-community SDs, as the stage 3 target.

## Sizing

Ten communities per cell, paired within community. With a paired SD s of the headroom-share difference and the agreed meaningful difference of 0.10, about 7.85 × (s / 0.10)² communities give 80 percent power for a two-sided 5 percent test. The pilot's between-community SDs of headroom at 300 and 900 sites, about 1 percentage point against means of 7 to 10, put s near 0.10 to 0.15 if arm differences vary as much as the oracle does, so ten communities resolve differences of about 0.10 to 0.12. Smaller differences are reported as inconclusive, not absent.

## Verification before comparison fits

- The adapters' certification tests from the first comparison are re-run under the pinned runtime: the Hmsc conditional field prediction, the gllvm conditional-moment adapter's normal integrals and training-location identity, and occJSDM's spatial predictor against a kernel-only basis at 50 lattice points.
- The scoring code is the pilot's, extended for the headroom share and the arm contrast, with tests that reproduce known values on synthetic predictions, including H equal to 1 when the arm's prediction equals the oracle's, 0 when it equals the env arm's, and negative when worse.
- Input hashes match the pilot's manifest; the oracle reference chosen for each dataset is recorded with its file hash and screen status.
- An independent script recomputes every paired interval and the development reading from the saved per-dataset scores.
- The timing pilot: one 900-site dataset per truth for each of the five arms (10 fits, kept as study attempts), including prediction. The all-extension projection prices every Bayesian fit at 4.43 times its initial cost, as the first comparison did. Dispatch proceeds only if the projected remaining wall time at eight single-threaded workers is under 48 hours; otherwise execution stops and Doug decides, and no cell, arm or diagnostic is dropped unilaterally. The gp support-point count and the Hmsc spatial method are frozen from this pilot's measured costs before dispatch and recorded in the config.

## Execution and provenance

- New worktree and branch for stage 2 under `~/src/occJSDM/.worktrees/`, built from the pilot branch, with the first comparison's adapters, certification tests and launcher conventions brought across from `codex/spatial-comparison-plan`. The pinned R 4.5.0 runtime and the first comparison's isolated package library for Hmsc and gllvm; versions and hashes recorded in the run config.
- Run index at `dev/simstudy/results/spatial-stage2-20261011/` (ignored), job files on `/Volumes/occJSDM Research/occJSDM/spatial-stage2-20261011/jobs/`, compact scores and the verification record committed under the study directory in `dev/simstudy/`.
- Single-threaded workers, long jobs detached with `nohup`, no refit of an existing job file, status files per job, and the pilot's lighter provenance: seeds, hashes, conditioning records where the adapters expose them, independent recomputation. Thresholds, the decision rule and the attempt ceiling are written to the config before the timing pilot.
- Report: the headroom-share table by arm and cell, the paired intervals, the development reading, the cross-truth misspecification table (each arm under the truth it does not match), diagnostics by arm, the Hmsc and gllvm range definitions with explicitly descriptive half-correlation distances, and the statement, verbatim from the pilot spec, that a single-season footprint is not evidence of dispersal as a mechanism.

## Risks stated in advance

- Hmsc at 900 sites with the full spatial method may be too slow; the predictive-process fallback changes the Hmsc arm's approximation and is recorded as such. hmsc6 may converge poorly; an unresolved hmsc6 cell is reported with its descriptive scores and weakens, not voids, the factor-dimension reading.
- occJSDM's GP may be unresolved at 900 sites as it was at 100; its descriptive scores are still informative about the representation, and the status is part of the finding.
- The oracle's headroom uses true intercepts and slopes, so H understates what a perfectly specified fitted model could reach. The report says so.
- Ten communities resolve differences of about 0.10 and nothing finer.

## Revision record

11 October 2026, at planning. The convergence screen for fitted arms uses the first comparison's truth-free training-side quantities as implemented in its `diagnostics.R`: per-species training probability mean, per-site training field, realised field SD, intercepts, coefficients, spatial prior SD and range (occJSDM), coefficients and the species spatial correlation at distance 0.5 (Hmsc), objective, gradient and native convergence (gllvm). The held-out probability average named in the design is not screened, because held-out prediction runs only after selection. No threshold changes.

## Future work, out of scope here

- Stage 3: the trait-linked per-species range prototype in occJSDM, with its target set by this study's headroom shares.
- Imperfect detection and eDNA replicate structure, applied to the winning representation.
- A mixed truth combining shared gradients and per-species fields.
- The clustered arrangement and the 0.15 and 0.03 ranges.

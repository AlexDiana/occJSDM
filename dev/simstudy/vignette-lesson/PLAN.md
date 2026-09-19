# First truth-based occJSDM lesson implementation plan

**Goal:** Produce one rendered non-spatial teaching lesson whose data, simulated truth, fitted results and detection examples belong to the same reproducible experiment.

**Architecture:** Use the existing simulator and fitter unchanged. An explicit build script runs three fits and exports a compact teaching bundle; the R Markdown lesson renders that bundle without running MCMC. Keep the original reference walkthrough available and link it to the new lesson.

**Tech stack:** R, occJSDM at main revision b53048a, ggplot2, knitr, rmarkdown, testthat and posterior for diagnostics.

**Spec:** The agreed design is copied to `DESIGN.md` beside this plan. This plan implements its first non-spatial lesson, including the four detection cases and a prior-sensitivity comparison. Spatial examples, the remaining reference-output lessons and Paper2Agent execution are subsequent work.

## Constraints

- Say variation partitioning in prose and figure labels; retain existing API names.
- Use at most six PCR replicates per primer.
- Declare seed and design before examining recovery; do not select cases using fitted probabilities.
- Distinguish underlying occupancy probability, realized site presence, sample presence and PCR outcomes.
- Compare the threshold-adjusted detection truth with fitted p and q.
- Report poor recovery, convergence limitations and the absence of held-out prediction honestly.
- Retain all generating inputs, identifiers, seeds, MCMC settings, source hashes and fit diagnostics.
- Do not alter sampler code, model defaults, shipped sample objects or beta-release criteria.

## Review focus

1. Row and species reordering: joins use declared identities rather than relying on table order.
2. Missing PCR results: shown as missing, never counted as negatives.
3. Laboratory contamination at an occupied site: labelled by sample state, distinct from field contamination.
4. Nominal event rates versus threshold-positive rates: use the rounded-lognormal survival probability.
5. Numerical stability and saved results: expose diagnostics and invalidate mismatched source/input caches.

## Task 1: Matching inputs and trustworthy case labels

Files: `helpers.R`, `test_lesson.R`, `build_lesson.R` in this directory.

Interfaces: `make_lesson()` returns the complete simulation and generating settings. `lesson_observations()` returns a row per species/PCR with IDs, counts, thresholded outcomes and actual z/w. `select_lesson_cases()` selects four truth-defined sample cases before fitting. `effective_detection_rate()` maps event rates to fitted-scale truth.

- [x] Add tests with manually labelled true, laboratory-false and field-false detections, reordered species/rows and a missing PCR. Check threshold-one survival against its analytic expression and a higher threshold.
- [x] Run `Rscript dev/simstudy/vignette-lesson/test_lesson.R`; observe the missing helper failure.
- [x] Implement helpers and pass those tests. Simulate 100 sites, 10 species, two field samples per site, two primers and six PCRs per primer, seed 20260919, two measured environmental covariates, two measured traits and two hidden community factors, without spatial effects.
- [x] Record category counts and selected IDs before fitting. Use all cases when summarizing performance.

## Task 2: Fits, truth comparisons and compact evidence

Files: `build_lesson.R`, `summarise_lesson.R`, `verify_lesson.R`, `vignettes/teaching-data/nonspatial-lesson.rds`.

Interfaces: the build writes the unchanged simulation, binary input made from its actual z, three matching full fits and a source/input manifest to a user-specified scratch directory. The summarizer exports data frames for figures, diagnostics and the selected cases into the compact RDS.

- [x] Fit the same community in binary mode with true z, then in two-stage mode with PCR observations and default priors, then with q and theta0 priors changed from Beta(1,20) to Beta(1,4). Keep the p prior unchanged. The alternative is an explicit stress test, not a recommended default.
- [x] Use four chains, 3,000 burn-in and 6,000 retained iterations initially. If important diagnostic warnings remain, extend the affected fit once to 6,000 burn-in and 12,000 retained draws; preserve the initial output and report remaining limitations.
- [x] Independently reconstruct occupancy probabilities from saved model draws. Compare reconstruction with saved two-stage psi, align fitted site/species/sample identities, and export posterior means and intervals alongside exact generating truth.
- [x] Export signed error and mean absolute error overall and within low/middle/high true-probability bands, rate recovery, selected-case probabilities, all-case classification summaries and parameter/probability diagnostics.
- [x] Run `Rscript dev/simstudy/vignette-lesson/verify_lesson.R SCRATCH_DIRECTORY` to check source hashes, exact-z input identity, threshold truth, case labels, posterior summaries and compact-bundle consistency.

## Task 3: Pedagogical lesson and rendered results

Files: `vignettes/occJSDM-first-lesson.Rmd`, its rendered `.md` and figures; targeted edits to `vignettes/occJSDM.Rmd`; a build README beside these scripts.

- [x] Explain the three stages and show actual site presence separately from its generating probability.
- [x] Plot true versus estimated occupancy for the perfect-observation and two-stage fits, with signed/absolute errors explained in percentage points.
- [x] Explain the priors as starting expectations, not proof of good practice. Plot generating p/q/theta0 beside fitted estimates and show the threshold correction.
- [x] Draw the four selected detection grids, true site/sample states and fitted probabilities. Explain mistakes and uncertainty using the actual computed results. Include the prior stress test on the same cases.
- [x] Add numerical diagnostics, version provenance, reproducible code and explicit limits of a single fitted-site example.
- [x] Link the lesson from the original walkthrough, correct its categorical claim that false positives are always weak, and use variation partitioning in reader-facing prose.
- [x] Render HTML and GitHub Markdown, inspect all figures, run the helper/evidence checks and the ordinary package test suite. Inspect the final diff and obtain an independent scientific/code review before reporting completion.

## Progress

Worktree: `codex/truth-based-vignette`, based on b53048a. The user approved the design and the false-positive teaching additions. Execution is local; publication and Paper2Agent installation are outside this first lesson.

All three tasks completed, 19 September 2026. The baseline and ordinary package tests passed (the opt-in coverage study was skipped). The lesson helper suite passes 24 assertions. The independent evidence checks pass for all three selected fits; both HTML and GitHub Markdown render without warnings. Seven figures were visually checked and clipped/overlapping labels corrected.

Independent review found no critical or scientific-formula defects. Its provenance finding was fixed by exporting each fit's R/package session metadata and summary-script hashes. The binary-input verification gap was fixed by checking the observations stored in the actual fit, with regression cases for changed outcomes, reordered identities and fitter-added columns. Both field samples' posterior DNA probabilities and true/fitted collection probabilities were added to the worked cases. These latter items were treated as required verification/pedagogy completeness, rather than deferred polish. All review findings are addressed.

The alternative-prior fit was extended once. Two named parameters remain above Rhat 1.01 (maximum 1.014343); minimum parameter ESS is 360.9. The lesson reports these limits. No further MCMC or release-validation work is implied by completion of this teaching lesson. The complete raw archive remains in the task workspace at work/vignette-lesson-20260919.

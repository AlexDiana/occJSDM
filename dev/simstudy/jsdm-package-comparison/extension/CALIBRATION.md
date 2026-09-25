# Saved-fit bias and interval coverage

Doug approved this Lesson 4 extension on 25 September 2026, including occJSDM, gllvm, Hmsc and sjSDM. It reanalyses the completed September experiment. No package code, fitting settings, selected starts, datasets or diagnostics are changed, and no model is refitted.

## What is measured

- Test-site probability bias, mean absolute error and RMSE use all 300 held-out sites and every species in each of the 626 scored fits. These predictions already integrate over unknown hidden conditions. Each independent community receives equal weight within its design cell. The aggregate RMSE is the square root of the mean community MSE.
- Marginal probability intervals use five fixed raw environmental settings: `(0, 0)`, `(-1, 0)`, `(1, 0)`, `(0, -1)` and `(0, 1)`. All packages receive the same grid, transformed using the original training centres and scales. These interval results are explicitly distinct from the 300-site point-error results.
- occJSDM and Hmsc intervals use all saved global draws from all four chains, integrating hidden conditions within each draw before calculating the mean and 2.5%/97.5% quantiles. Hmsc's probit-normal integral is analytic. occJSDM's quadrature is refined through 31, 61, 121 and, if required, 241 nodes until the mean and both endpoints change by at most 0.0001 probability. A failure stops extraction. The earlier posterior-mean thinning check is not used to certify interval endpoints.
- Conditional environmental coefficients are compared on original predictor units for occJSDM, gllvm and sjSDM. Only occJSDM has saved uncertainty for this target. The intercept and slopes are transformed jointly within every posterior draw. Hmsc is excluded because probit coefficients do not share the logit truth. As a conservative whole-cell scope, the curved scenario fitted with a linear response is excluded from coefficient recovery, including the unaffected species; its probability results remain included.
- Existing occJSDM and gllvm trait-effect intervals are reused with their saved scale correction. Hmsc's trait coefficients have no exact logit numerical truth, and sjSDM did not enter the trait experiment. Missing intervals are unavailable, never zero coverage. Nonfinite or reversed endpoints are excluded from interval summaries and counted explicitly.

Within a community, coverage and interval width are averaged over eligible elements of the same target and term. The study then gives each community equal weight. Coverage MCSE is the standard deviation of the community coverage fractions divided by the square root of their count. Bias MCSE is calculated analogously. Species, settings and posterior draws are not treated as independent simulation replicates. These are pointwise intervals, not simultaneous bands. With ten communities per design cell, coverage and its MCSE remain exploratory; the varying ecological coefficients also mean this describes the declared distribution of communities rather than repeated samples at one fixed coefficient vector.

All scored fits remain in the main summaries, including unresolved diagnostics. Separate passed-only summaries are a sensitivity analysis with potentially selective missingness. The complete 640-row manifest preserves the 14 failed gllvm fits and all 121 flagged fits. No undefined score is imputed for a failed fit, and no package is selected using truth error.

## Reproduction

Run from the repository root. The original archive is read-only input; use a new output directory for a new analysis. The exporter checks source and input hashes before reusing its own cached job extractions. The original fits and scores are never overwritten.

```sh
study_archive=/Users/douglasyu/src/occJSDM/dev/simstudy/results/lesson-4-extension-20260923
calibration_output=/Users/douglasyu/src/occJSDM/dev/simstudy/results/lesson-4-calibration-20260925/reproduction
calibration_code=dev/simstudy/jsdm-package-comparison/extension

Rscript "$calibration_code/test-calibration.R" "$calibration_code"
Rscript "$calibration_code/export-calibration.R" "$study_archive" . "$calibration_output"
Rscript "$calibration_code/verify-calibration.R" "$calibration_output/calibration.rds" "$study_archive" &&
  cp "$calibration_output/calibration.rds" vignettes/teaching-data/lesson-4-calibration.rds
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-4.Rmd")'
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-4.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE))'
```

Always run the independent verifier successfully before copying, using or distributing a newly exported bundle; stop if any command fails. The exporter writes only to the analysis output directory. A new MD5 merely records what was read; it does not establish that a truth file is correct. The verifier reconstructs all 50 communities using the original frozen simulator, checks their complete truth objects, reconstructs each training input, and compares saved test observations and marginal truth with those communities. It separately recomputes all 626 point scores, fixed-grid truth, coefficient backtransformations and grouped summaries. It checks input/source hashes, unavailable intervals and the numerical refinement records. It never fits a model or changes the archive.

The completed analysis uses an immutable source copy under `dev/simstudy/results/lesson-4-calibration-20260925/source-v1/` and writes resumable per-job extraction records, CSV summaries and logs beneath that analysis directory. The compact `vignettes/teaching-data/lesson-4-calibration.rds` contains the point metrics, grid/coefficient records, community summaries, diagnostics, integration records and provenance required to render Lesson 4 without the 22 GB fitting archive. `calibration-section.Rmd` is the reusable source for sections 14-15; the same text is included in `lesson-section.Rmd` and the full lesson.

The immutable first snapshot also wrote a copy inside its isolated snapshot directory. That copy was not used by the lesson until the independent verifier passed. Review then removed automatic teaching-bundle installation from the reusable exporter; its numerical extraction is unchanged. The reproduction commands above make the verification and installation steps explicit.

## Verification and observed results, 25 September 2026

Extraction processed all 626 scored fits in about 6.3 minutes, followed by independent verification of the 50 simulated communities, 626 point scores and 444 grouped summaries. Every logit probability calculation met the mean/endpoint refinement tolerance; the largest observed change was 0.000002955 probability. The compact bundle is about 1.15 MB. The probability identities, paired-simulator checks, calibration checks, runner-resume checks and lesson-link checks passed. Both HTML and Markdown render successfully, and the two added figures were visually inspected. The review concerns about truth identity and premature teaching-bundle installation are addressed above.

Baseline probability coverage across the five settings was 95.4%/94.8% for occJSDM at 100/300 sites and 96.6%/95.6% for Hmsc. These are exploratory averages over ten communities, with original diagnostic flags retained. Under curved truth at 300 sites, linear-response coverage was 64.2%/65.0% for occJSDM/Hmsc; including the quadratic term gave 93.6%/94.8%. The lesson reads these values from the saved bundle. The lower occJSDM coverage in rare-species and correlated-environment conditions is also shown, without assigning an untested mechanism or declaring a software defect. gllvm/sjSDM probability coverage remains unavailable from these saved records.

## What still requires new work

gllvm/sjSDM probability intervals and ordinary coefficient uncertainty require a declared uncertainty method and further computation, potentially targeted refits. More independent communities are required for a precise coverage assessment. A probit-generating arm is needed to examine link sensitivity fairly. The existing unresolved optimisation and MCMC diagnostics must remain visible when planning further fits. The old July/August occJSDM-only study predates model changes and is not pooled into this extension.

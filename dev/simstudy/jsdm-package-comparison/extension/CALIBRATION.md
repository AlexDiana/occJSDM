# Saved-fit bias and interval coverage

Doug approved this Lesson 4 extension on 25 September 2026, including occJSDM, gllvm, Hmsc and sjSDM. It reanalyses the completed September experiment. The original package code, fitting settings, selected starts, datasets and diagnostics are retained. The initial posterior extraction required no refits. The native coefficient extension subsequently replays selected gllvm fits to recover their fitting objects and computes sjSDM uncertainty directly from its saved weights, as described below.

## What is measured

- Test-site probability bias, mean absolute error and RMSE use all 300 held-out sites and every species in each of the 626 scored fits. These predictions already integrate over unknown hidden conditions. Each independent community receives equal weight within its design cell. The aggregate RMSE is the square root of the mean community MSE.
- Marginal probability intervals use five fixed raw environmental settings: `(0, 0)`, `(-1, 0)`, `(1, 0)`, `(0, -1)` and `(0, 1)`. All packages receive the same grid, transformed using the original training centres and scales. These interval results are explicitly distinct from the 300-site point-error results.
- occJSDM and Hmsc intervals use all saved global draws from all four chains, integrating hidden conditions within each draw before calculating the mean and 2.5%/97.5% quantiles. Hmsc's probit-normal integral is analytic. occJSDM's quadrature is refined through 31, 61, 121 and, if required, 241 nodes until the mean and both endpoints change by at most 0.0001 probability. A failure stops extraction. The earlier posterior-mean thinning check is not used to certify interval endpoints.
- Conditional environmental coefficients are compared on original predictor units for occJSDM, gllvm and sjSDM. occJSDM supplies posterior uncertainty; the native extension adds gllvm and sjSDM Wald intervals for compatible ecological scenarios. The intercept and slopes are transformed jointly within every posterior draw or through the full within-species coefficient covariance matrix. Hmsc is excluded because probit coefficients do not share the logit truth. As a conservative whole-cell scope, the curved scenario fitted with a linear response is excluded from coefficient recovery, including the unaffected species; its probability results remain included.
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

gllvm/sjSDM marginal-probability intervals still require a separate uncertainty method and further computation. Their ordinary ecological coefficient intervals now use the native procedures below. Species-specific slopes in the hierarchical gllvm trait models combine fixed and random effects; this extension does not invent their joint uncertainty. Existing measured trait-effect intervals remain available. More independent communities are required for a precise coverage assessment. A probit-generating arm is needed to examine link sensitivity fairly. The existing unresolved optimisation and MCMC diagnostics must remain visible when planning further fits. The old July/August occJSDM-only study predates model changes and is not pooled into this extension.

## Native gllvm and sjSDM coefficient intervals

Doug approved native coefficient intervals first, with bootstrap probability intervals deferred. The scope is the same compatible ecological coefficient target already present in the bundle: baseline, rare species, correlated environments and correctly specified quadratic responses. There are 80 planned fits per package in these cells; 79 gllvm fits and all 80 sjSDM fits were scored in the original experiment. The original failed gllvm fit remains missing. The trait experiment retains its existing measured trait-effect intervals, while species-specific hierarchical gllvm slope intervals remain unavailable.

For gllvm 2.0.15, the worker replays only the archived selected start, using the original seed, inputs, VA method, optimiser settings and two latent variables. It enables native standard errors and requires coefficients, loadings and log likelihood to agree with the archive within 0.00000001. It extracts each species' full coefficient block from the native covariance matrix. The original point estimates and fit diagnostics are retained. A non-positive-definite joint covariance flags the calculation; a non-positive-definite species block has no interval. Finite native intervals from flagged fits remain in the main summary, with a separate sensitivity summary requiring both original fit and new uncertainty checks to pass.

For sjSDM 1.0.7 from the frozen v0.2.1 source, the worker rebuilds the CPU float64 linear model and restores its exact saved coefficients and loadings. It never optimises the model. A scoped observer records the matrices already inverted by the unmodified native Python `model.se()` function and is removed after each call, including errors. Native standard errors are checked against the recorded covariance diagonal, and both coefficient and loading arrays must remain unchanged. This preserves the native Hessian calculation, including its addition of 0.001 to every Hessian entry. That adjustment is not changed into a diagonal ridge, and the fitting optimiser's weight-decay curvature is not added to the standard-error routine. These are native conditional intervals: associations and the other species' coefficients are held fixed.

The sjSDM routine's supported `sampling` argument is increased from its default 100. Two independent streams are evaluated at 2,000 draws; if any raw-scale standard error differs by more than 5% relative to the larger of the two, both calculations are repeated at 8,000 and then 32,000 draws. The first stream at the final level defines the reported interval; the second is a numerical audit. The rule is fixed before reading truth and does not select a seed based on coverage. Remaining instability is flagged. A calculation error produces missing intervals and a recorded error. Ordinary fit flags, uncertainty flags and failures are distinct; ten original communities remain the replication limit.

Both packages' native coefficient covariance blocks are transformed with the same affine coefficient transformation as the point estimates. Intervals are the transformed estimate plus or minus the 0.975 normal quantile times its standard error. This includes intercept-slope covariance and the already-centred quadratic basis where applicable. No interval is inferred from standard errors alone when their covariance is needed.

The native worker reads training inputs and archived selected fits, never truth files. Its results are saved separately under `dev/simstudy/results/lesson-4-native-intervals-20260925/`. The first immutable source snapshot is `source-v1/`; its original serial log is retained after the run was split into three disjoint resumable partitions. The reusable worker also hashes `pilot-math.R`; for the first snapshot, the augmentation provenance explicitly records that dependency separately. No original archive record is overwritten.

To reproduce, first export and verify an unaugmented base bundle using the earlier instructions, but retain that base at a separate path. Then run:

```sh
native_launcher=/Users/douglasyu/Documents/Codex/2026-09-10/fam/work/lesson-n-environment-20260922/run-r
native_output=/Users/douglasyu/src/occJSDM/dev/simstudy/results/lesson-4-native-intervals-20260925/reproduction
native_base="$calibration_output/calibration.rds"

Rscript "$calibration_code/test-native-intervals.R" "$calibration_code"
"$native_launcher" "$calibration_code/test-native-sjsdm.R" "$calibration_code"
"$native_launcher" "$calibration_code/native-intervals.R" "$study_archive" "$calibration_code" "$native_output"
Rscript "$calibration_code/add-native-calibration.R" "$native_base" "$native_output" "$calibration_code" "$native_output/calibration.rds"
Rscript "$calibration_code/verify-native-calibration.R" "$native_output/calibration.rds" &&
  Rscript "$calibration_code/verify-calibration.R" "$native_output/calibration.rds" "$study_archive" &&
  cp "$native_output/calibration.rds" vignettes/teaching-data/lesson-4-calibration.rds
```

Stop after any failure. Neither exporter installs a teaching bundle. The additional independent verifier checks the native covariance identities, raw-scale interval arithmetic, archive and source hashes, reconstruction differences, and the combined sensitivity summary. It also checks that all original point estimates, probabilities, Bayesian intervals, trait intervals and diagnostic classifications are unchanged. The original verifier then repeats the truth and all grouped-summary checks on the augmented bundle.

## Completed native calculation, 25 September 2026

All 159 eligible calculations completed and supplied 5,170 coefficient intervals. All 79 gllvm selected-start reconstructions matched their archived coefficients, loadings and log likelihood exactly. Twenty gllvm joint covariance matrices failed the positive-definiteness check; their finite species blocks remain in the main results and their jobs are excluded from the combined sensitivity summary. All 80 sjSDM calculations passed the independent-stream check: 77 used 2,000 integration draws and three required 8,000. The largest final relative raw-scale SE difference was 0.04425. No original fit, coefficient point estimate, probability result or trait interval changed.

For the baseline at 100 sites, native gllvm coverage was 94%, 95% and 98% for the first slope, second slope and intercept; sjSDM coverage was 82%, 88% and 86%. At 300 sites the corresponding values were 89%, 98%, 91% and 92%, 97%, 93%. These are averages over the ten communities, with flags retained; they are not a package ranking. The lesson reports MCSE and interval widths, the original diagnostic sensitivity and the stricter native-interval sensitivity. A flagged gllvm rare-species calculation has exceptionally wide slope intervals and remains visible on the full-range coefficient-width plot.

Both independent verifiers passed: native covariance and transformation checks for all 159 jobs, then the original truth checks, all 626 point scores and 444 community-weighted summaries. The native transformation and sjSDM zero-loading logistic-information tests, study probability/simulator/runner regressions, and lesson-link checks passed. HTML and Markdown rendered without warnings. The new coefficient coverage/width figure was visually inspected. The immutable report source is `report-v1/`; the full run, verification and render logs accompany the saved outputs.

## Combined four-package presentation

The lesson now presents baseline probability and coefficient results in shared tables with the same four package columns. Each row identifies one target and measure; values from different targets are never averaged together. Probability bias, RMSE and intervals in these combined tables all refer to the five fixed environments, while section 14 retains the separate 300-test-site prediction results. Bias and coverage retain their community MCSE. The trait comparison uses the same format, and figures retain all four package positions with explicit missing-result labels.

`Unavailable` identifies an interval not calculated; `Different link` identifies the lack of a matching numerical coefficient truth for Hmsc; `Not fitted` identifies sjSDM's omitted trait experiment. Community counts remain explicit. The combined diagnostic table compares the same coefficient-compatible ecological scope for all four packages. The sensitivity display applies original fit checks to posterior intervals and both original fit and native interval checks to the new gllvm/sjSDM coefficient intervals; it does not substitute one package's retained communities for another's.

The presentation helper is `vignettes/lesson-4-calibration-tables.R`. This revision changes no saved fits, numerical summaries or compact-bundle data. Rendering and cell checks use that verified bundle directly. Presentation logs and the previous lesson sources are preserved under `dev/simstudy/results/lesson-4-combined-presentation-20260925/`.

## Presentation following the validation article

Sections 14-15 now use the validation article's report style: a plain-language findings summary, separate compact tables for prediction error, coefficient bias and interval coverage, neutral horizontal plots with zero-bias or 95% reference lines, and a disclosure containing all previous detailed results, widths and methods. Scoped styles in `vignettes/lesson-4-calibration.css` leave the worked lesson unchanged. HTML and GitHub Markdown show the same data; the rendering code stays in the R Markdown source.

The summary plots show one MCSE on either side of the community average, calculated from independent community summaries. They do not reuse the older validation article's binomial error bars based on individual interval decisions. Missing package intervals retain labelled panels. The baseline probability errors still use 300 test sites; probability coverage still uses five fixed environmental settings. Coefficient targets remain separate, and all main summaries retain diagnostic flags.

The presentation archive is `dev/simstudy/results/lesson-4-validation-presentation-20260925/`, including the preceding source files, bundle checksum, render logs and independent display checks. No model fits or saved calibration data are changed.

HTML and Markdown rendered successfully. Independent display checks reconstructed all 64 prediction-summary cells from their community records and checked the compact table values, plot MCSE endpoints and missing-package labels. The HTML structure check found all eight detailed tables and three detailed figures inside the closed disclosure, with unique anchors. The summary PNGs were visually inspected, lesson-link checks passed, and review found no remaining presentation issues. The saved calibration bundle is byte-identical to the preceding presentation.

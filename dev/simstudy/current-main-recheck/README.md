# Current-code comparison: reproduction and evidence

Read the [report](REPORT.md) for the results and the [protocol](protocol.md) for the design fixed before comparing errors. This is a paired rerun of the main non-spatial comparisons in PR #11, not the full alternative-prior grid or the four-package Lesson 4 study.

The historical package revision is `80d449dc8593b272d5fae1f15407356aa41d6c3b`. The current production revision for this comparison is `2a75bf1ed07bf37faa46c2800c3a7e61fa71e5bc`. Production fitting code is unchanged by this research branch. Numeric-trait preprocessing changed between these revisions; the resulting hierarchical prior is not assumed to be equivalent.

The local execution archive is `/Users/douglasyu/src/occJSDM/dev/simstudy/results/pr11-current-20260927`. The three original study archives are under `/Users/douglasyu/Documents/Codex/2026-09-10/fam/work`: `nonspatial-bias-recheck-20260914`, `nonspatial-design-recheck-20260919`, and `jsdm-sample-size-20260919`. Their inputs, random states, fits and results are read-only inputs to this comparison. Full posterior fits are intentionally excluded from Git; compact evidence in `results/` includes paths and hashes. Reproduction needs those input archives and the R dependencies recorded in `results/environment.txt`.

To reproduce in a fresh archive, export the production revision with `git archive`, including `R`, `src`, `DESCRIPTION`, `NAMESPACE`, `LICENSE`, `man` and `data`, into its `source-main/` directory. Write the full revision to `source-revision.txt` and install the exported package into its own `library/` using `R CMD INSTALL --library=... source-main`. Do not install the historical build or invoke the historical `load_frozen()` helpers. The runner verifies the loaded package path and records source and library hashes.

From the repository, run the following sequence, substituting absolute paths for `REPO`, `STUDY` and `ARCHIVES`. `STUDY` must be a fresh archive for a new run; resume is permitted only when source hashes and settings agree.

```sh
Rscript dev/simstudy/current-main-recheck/run.R --repo=REPO --study=STUDY --archives=ARCHIVES --mode=pilot --workers=4
Rscript dev/simstudy/current-main-recheck/verify.R --repo=REPO --study=STUDY --mode=pilot
Rscript dev/simstudy/current-main-recheck/run.R --repo=REPO --study=STUDY --archives=ARCHIVES --mode=initial --workers=4
Rscript dev/simstudy/current-main-recheck/select.R --repo=REPO --study=STUDY --mode=plan
```

The selector writes `long-selection.csv` and `long-keys.txt` before any current-versus-historical error comparison. Supply that comma-separated key list as the `--keys` argument in the next command. It includes the ten historically extended cases and any additional current initial fits flagged by warnings or the prespecified Rhat screen. Retain all initial fits.

```sh
Rscript dev/simstudy/current-main-recheck/run.R --repo=REPO --study=STUDY --archives=ARCHIVES --mode=long --workers=4 --keys=KEYS
Rscript dev/simstudy/current-main-recheck/select.R --repo=REPO --study=STUDY --mode=final
Rscript dev/simstudy/current-main-recheck/verify.R --repo=REPO --study=STUDY --mode=selected
Rscript dev/simstudy/current-main-recheck/summarise.R --repo=REPO --study=STUDY
Rscript dev/simstudy/current-main-recheck/diagnostics.R --repo=REPO --study=STUDY
Rscript dev/simstudy/current-main-recheck/chain-sensitivity.R --repo=REPO --study=STUDY
```

`helpers.R` imports only the named scoring functions and helper definitions from the archived research code; it never invokes their fitting entry points. `run.R` uses the saved input RNG state for every fit. The known-site-factor control is a private clone of the current fitter and does not modify the installed namespace. It is a research control, not an additional public model option.

`verify.R` independently reconstructs every selected fitted-site mean from all retained intercept, environmental, site-factor and loading draws. It binds all 110 actual result jobs to the selection and original input hashes, checks current trait preprocessing, checks the declared and actual MCMC schedule, and checks all four primary probability-band scores per fit. `summarise.R` requires that verification before exporting the paired results and plots.

Each community has equal weight. All primary comparisons use the same original 100 sites. `delta` is current minus archived, so a negative change in MAE means improvement. Confidence limits are 95% t intervals for paired changes across ten communities. They are not posterior intervals or coverage estimates. Both versions' initial fits are exported as a sensitivity comparison. RMSE in the new summaries is the mean of the ten community RMSEs; the historical design report instead used the square root of the mean squared community RMSE, so those aggregate RMSE columns should not be compared without accounting for that definition.

`design-paired-summary.csv` separately compares designs within each code version. Its `reduction_mae` is reference-design error minus alternative-design error, so positive values mean the alternative helps. Rare/common species-prevalence rows and the obsolete nominal-q comparator are retained in the raw community exports but excluded from the fixed-ten aggregate summaries; the available rare-species cases cannot support a general assessment.

The descriptive chain-sensitivity check was added after longer-run convergence flags appeared. It includes every selected fit still flagged by warnings or the Rhat screen, verifies each pooled mean against the independently checked result, and scores the individual-chain posterior means as well. In `flagged-chain-scores.csv`, chain zero denotes the pooled estimate. `chain-sensitivity-summary.csv` keeps all ten communities and takes the lowest and highest MAE and signed bias obtained by using either the pooled or an individual-chain mean in each flagged fit; other fits retain their selected pooled estimate. `design-chain-sensitivity.csv` propagates those choices into design comparisons. Different groups and design arms can attain their extrema using different chains, so the endpoints do not represent one jointly selected set of fits. These ranges describe sensitivity to the observed chains. They are not confidence intervals, convergence guarantees, or bounds on the true posterior; unvisited modes remain possible.

The report does not rerun the historical interval-coverage study, establish a universal error floor, measure held-out prediction, or approve a beta-release target. PRs #13 and #14 remain pending Alex's review.

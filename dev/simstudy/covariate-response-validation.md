# Environmental response output correction

Status: implemented on `codex/fix-covariate-response`, awaiting Alex's review.

`returnCovariateEffect()` and `plotCovariateEffect()` previously returned impossible occupancy probabilities. The numeric calculation converted the environmental term to a probability and then added the intercept: `B0 + plogis(XB)`. The intercept is measured in log-odds, so it must be added **before** the conversion. The stored `X0_psi` also contains standardized numeric predictors, despite its historical description as raw data; standardizing it again gave the wrong x-axis and environmental contribution. The categorical path omitted the intercept and first level. Both public wrappers ignored the requested `confidence`.

## What the corrected curve means

Imagine changing temperature while keeping the rest of the environment fixed. Each point answers: **what occupancy probability does the fitted model predict at this temperature, under those reference conditions?**

- The chosen numeric predictor varies over its observed range, shown in its original input units. A categorical predictor varies over every fitted level, including the first level.
- Other numeric predictors stay at their observed medians. Other categorical predictors stay at their first fitted levels.
- Hidden site effects and spatial contributions are set to zero.
- For each retained posterior draw, all environmental contributions are added to the intercept. Only then is the sum converted to an occupancy probability: `plogis(B0 + X_grid %*% beta)`. Continuous-response models retain the untransformed sum.

This is a conditional environmental response curve. It is not a prediction for any particular sampled site, and it does not average over possible hidden site effects. A comparison against simulation truth must use the same reference conditions and zero hidden contributions. It also does not establish that the environmental relationship is causal.

Numeric plots show the posterior median and an equal-tail credible interval. Categorical plots now show the same summaries as points and bars, replacing boxplots whose quartiles and whiskers did not represent the requested interval. The default is 95%; `confidence = 0.8` displays the 10th to 90th posterior percentiles. Categorical return tables still contain every draw, so changing `confidence` changes the plotted summary, not those raw draws. For compatibility, the numeric table still calls its median column `mean`; the help now explains that historical name.

## Reproduction with the unchanged teaching fit

The validation used the archived Lesson 3 `default-fit.rds`, MD5 `67a364872efeb742de94e705194a7c33`, with four chains and 6,000 retained draws per chain. No new model was fitted.

| Output for `X_psi.EnvCov.1` | Previous output | Corrected output |
| --- | --- | --- |
| OTU_1 median probabilities | 1.772 to 1.815 | approximately 0.625 to 0.912 |
| OTU_10 median probabilities | -0.610 to -0.484 | approximately 0.040 to 0.898 |
| Predictor range displayed | -2.018 to 2.816 | -19.432 to 30.947, matching the original input |
| Medians outside 0 to 1, these two species | 400 of 400 | 0 of 400 |

The corrected medians and 80% interval endpoints agree with an independent calculation at the beginning, middle and end of each curve, to tolerance `1e-12`. That calculation reads original covariates from `data_info` and extracts posterior draws chain by chain, independently of the corrected grid builder and draw reshaping. The plotted data equal the returned data.

These values verify the output arithmetic. They do not show that the model has learned the true ecological response accurately. The fitting algorithm, posterior draws and existing occupancy recovery results have not changed.

To repeat the archived-fit check using an installation of this branch:

```sh
R_LIBS=/path/to/library-containing-this-branch Rscript \
  dev/simstudy/validate_covariate_response.R /path/to/default-fit.rds
```

The full archived fit is a local validation artifact and is not added to the package. The regression tests below are self-contained and need no archived fit.

## Regression coverage

`tests/testthat/test-covariate-response.R` uses known posterior draws and independently computed expected responses. The original implementation fails these tests; the corrected implementation passes. Cases include the full intercept-plus-environment calculation, original units, numeric medians for other predictors, every categorical level, ordered and custom factor coding, unused levels, training spline knots, predictors named `x` and `x2`, non-default intervals, species subsets and labels, singleton dimensions, continuous responses and invalid arguments. A separate test first reproduced the custom full-rank contrast issue found during independent review, then passed after correcting the factor-block width.

```r
devtools::test(filter = "covariate-response")
```

The code deliberately reuses fitted categorical encodings rather than assuming that the current R contrast options are those used during fitting. If a factor level was unused, reconstruction is checked against the observed training rows; incompatible contrast settings produce an explicit error instead of an invented response.

The duplicated plotting calculations, including unreachable copies of the old arithmetic, were removed. Plotting now uses the same returned responses. Existing fits with `X0_psi` and its scaling metadata work without refitting; the existing error for older fits missing that metadata remains.

## Package checks

The final source test suite reports **706 passes, no failures, errors or test warnings, and one skip**, including the new output tests and the existing recovery checks. The opt-in coverage study is skipped. The focused output tests contain 65 passing expectations. Independent code review found the custom factor-width issue described above; its regression test fails before the width correction and passes afterward.

`R CMD check --no-manual --ignore-vignettes` (built with `--no-build-vignettes`) completed with **0 errors, 3 warnings and 3 notes**. Its installed-package tests report 670 passes, no failures or test warnings, and eight expected skips (six CRAN-skipped checks and two requiring source access). Vignettes were not rebuilt for this output-only fix.

The warnings concern unchanged files: the R header's unsupported compiler warning option, the undocumented `predictNewSites(verbose)` argument, and GNU extensions in `src/Makevars`. Notes concern the existing LICENSE declaration and package-wide global-symbol checks, plus the `.git` pointer file included by building in a worktree. None is a probability-calculation failure; this is not a claim of a warning-free package build.

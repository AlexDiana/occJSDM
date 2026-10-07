# Remaining native plotting examples

This deliverable restores checked native `plotOccupancyGradient()` examples for both environmental predictors and separate `plotStage1FPRates()`, `plotStage2FPRates()`, and `plotDetectionRates()` examples. `remaining-plots-examples.Rmd` is mirrored in Lesson 3's sections. It contains the exact bodies executed to make the five exported PNGs. Knitting the lesson only reads `remaining-plots-data.rds` and the PNGs; it does not need the full fit or run MCMC.

The compact bundle contains truth, the original predictor means/standard deviations, native plot summaries, built coordinates, and provenance. It contains no posterior arrays or replacement fit. Every actual native call uses the unchanged `default-fit.rds`, with all 6,000 iterations from each of four chains.

Run from the repository root:

```sh
R_LIBS='/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples/library' Rscript dev/simstudy/vignette-lesson/remaining-plots-export.R /Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
R_LIBS='/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples/library' Rscript dev/simstudy/vignette-lesson/remaining-plots-verify.R /Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
```

`remaining-plots-diagnose-covariate.R` reproduces the pre-PR #13 covariate-response bug against the old archive and package only; it cannot pass on code after PR #13, which removed `plotCovariateEffect_base()` and corrected the calculation, so it is not run in a rebuild from current main.

The verifier also checks that Lesson 3 displays the exact exported plotting bodies. Keep the snippet and lesson synchronized when changing these examples.

Verification independently reconstructs every gradient's 95% interval and median from all 24,000 matched draws. It checks the standardized design against the original site predictors and verifies the generating `B0 + X B + U L` identity before removing site factors for these curves. The native grid has 40 points from the 2nd to 98th percentiles. Other environmental predictors stay at observed medians; site factors equal zero. This target does not marginalize unknown site factors. The checker verifies all 800 truth points, native intervals, medians, rugs, and facet identities against their plotted coordinates.

All 50 native rate intervals are checked against posterior draws. Field false-positive truth is generating `theta0`. Both laboratory helpers preserve species/primer identity; their native confidence intervals do not pool primers. Laboratory truth uses the corresponding read-threshold-adjusted event rate, separately for true and contamination reads. The verifier maps actual native colours to primers and checks that each black cross shares its matching error bar's dodged position.

The exporter and verifier check the full-fit/input hashes, source hashes, installed library location and hashes, and installed function bodies/formals against the current source. Confirmed archive hashes, for the three-sample rebuild of 2 October 2026 in `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples` (source commit eeb1675), are `035f82108bbd94b1230e82e725de1e1d` for `default-fit.rds` and `a788a24de838ee78e19ce3973b8956e5` for `input.rds`. The two-sample build's were `67a364872efeb742de94e705194a7c33` and `a34d5d3282519756b87c1f54a77e6f94`.

## Why `plotCovariateEffect()` was documented as blocked before PR #13

The cheap diagnostic reproduced the source defects against the two-sample full fit and the archived pre-PR #13 package. The line numbers and values below describe that package and fit. In `R/jsdmfun.R:420`, the numeric helper inverse-logit-transforms only the selected covariate contribution; line 426 then adds the log-odds intercept to that probability. The result equals `B0 + plogis(X_sub B)`, not `plogis(B0 + X B)`. The diagnostic independently reconstructs the first plotted interval to prove this is the executed calculation. For `OTU_1`, all 200 median points exceed one, ranging from 1.7716351 to 1.8151802. For `OTU_10`, all 200 are negative, ranging from -0.6104285 to -0.4837591.

There is a separate grid defect: `create_covariates_matrix()` standardizes `df` at `R/jsdmfun.R:115` and returns that same standardized object as `X0` at line 186. The effect helper treats `X0` as raw and standardizes it again. The original first predictor spans -19.432456 to 30.947095; stored `X0` and the plotted horizontal values span -2.018381 to 2.815735. Source inspection also confirms that the numeric effect helper omits other predictor contributions and that the wrapper does not forward its `confidence` argument (`R/output.R:616`). These latter observations are not needed to establish the out-of-range probability failure.

The lesson made this limitation explicit and used the independently checked native `plotOccupancyGradient()` example for response curves, which it still does. This deliverable changes no production R/C++ code, package library, stored fit, or dependencies, and never calls an MCMC fitting function.

The documented defect was corrected by [PR #13](https://github.com/AlexDiana/occJSDM/pull/13), which has merged. It also removed `plotCovariateEffect_base()` and renamed the numeric median column from `mean` to `median`. The archived reproducer and snippet remain unchanged to preserve their recorded provenance (the snippet's md5 is stamped in `remaining-plots-data.rds`, so editing it requires a re-export); Lesson 3 describes the correction outside the exported code. The verified gradient figures do not depend on that correction.

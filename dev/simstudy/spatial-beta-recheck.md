# Spatial fitting: review guide and beta checks

**PR #8 now contains the spatial correction only.** The optional continuous-noise prior is reviewed separately on branch `codex/continuous-noise-prior`, after PR #8. The non-spatial investigation and consolidated TODO are on branch `codex/nonspatial-bias-recheck` and can be reviewed independently. No simulations were rerun to make this split; the saved results and original source fingerprints are unchanged.

## Start the code review here

1. Read the [derivation and original validation](spatial-range-validation.md), including the earlier examples that still had large errors.
2. In `R/jsdmfun.R`, follow `spatialBasis()`, `computeSpatialSummaries()`, `spatial_range_logweights()` and `update_jSDMcoef()`. The range is selected after integrating over the coefficient block, and that entire block is immediately redrawn. The shared basis and this update order belong together.
3. Check the spatial prior means and equivalent matrix arithmetic in `src/jsdm.cpp`, and the same basis in `R/output.R`. `R/runOccJSDM.R` wires in the complete basis and preserves matrix dimensions.
4. Read `test-spatial-range.R`, `test-spatial-support.R`, `test-spatial-dense-algebra.R` and `test-spatial-grouped.R`. They check the mathematics and wiring, including zero residual factors and one support point.
5. Use the results below to assess the tested designs. The CSVs are supporting records, not files that need line-by-line code review.

All current default priors are unchanged. Binary, continuous, occupancy and two-stage models remain available. Previously saved spatial fits should be refitted with the corrected model.

## The agreed target

Doug provisionally chose an average systematic occupancy-probability error of no more than five percentage points in sufficiently informative simulations. Check low true probabilities (below 0.2), medium probabilities (0.2 to 0.8), and high probabilities (above 0.8) separately. For example, when the true probability is 20%, an average estimate between 15% and 25% meets this provisional target. This concerns average error across independently simulated datasets, not every prediction or the width of an uncertainty interval. Report the uncertainty in that average. Rare species need separate consideration because five points can be a large relative error for them.

The threshold was chosen before the stronger-data results below. Increasing the amount of simulated information is a diagnostic design choice; it does not change the target or turn earlier failures into successes.

## What the additional code does

- A user can now request one support point at every unique observed location. Previously the fitter silently capped the request at one fewer location and used clustering to select the points. Full support now uses the observed coordinates directly. This avoids losing a spatial direction solely because of the cap. The existing numerical jitter remains. The default count, 20% of unique locations rounded down, is unchanged; users must check whether increasing it materially changes their results. A small support set can still miss short-range spatial variation.
- Repeated binary observations at the same location now share the spatial part of the range calculation. Their individual environmental covariates, latent factors, outcomes and Polya-Gamma precisions are all retained. The calculation rearranges exactly the same sums. It is about 25 times faster for the range-weight calculation in the 4,800-observation benchmark, with a maximum range-probability difference of 1.5e-13. This is not a 25-fold speed claim for complete fits. Continuous fits retain the existing constant-precision shortcut.
## What the independent binary calculation establishes

The original binary example had six observations at each of 80 locations. occJSDM overestimated low occupancy probabilities by 14.13 percentage points and underestimated high probabilities by 12.60 points. An independent reference calculation, supplied with the true regression coefficients, range, spatial SD and full generating covariance, still had errors of +12.84 and -11.51 points. Most of the attenuation therefore remains even when these parameters are known. This is evidence of limited information, not permission to accept arbitrary package bias.

The reference uses [elliptical slice sampling](https://proceedings.mlr.press/v9/murray10a.html) of the spatial field, without calling occJSDM's sampler or spatial basis routines. A one-dimensional numerical integration checks its implementation. Applying the logistic transform to every saved package draw and then averaging reproduces the original probability estimates to numerical precision, so the earlier errors do not come from transforming a mean linear predictor incorrectly.

At 30 observations per location, the reference's low-probability error was still +5.23 points. At 60 observations per location, its low, medium and high errors were +3.29, -0.052 and -2.32 points, respectively. Their Monte Carlo standard errors were 0.028, 0.011 and 0.022 percentage points. Some individual spatial-field coordinates still mixed slowly, so the 60-repeat result is an information-adequacy pilot, not a general convergence claim or an independent-replicate beta test.

All eight species in those designs have mean prevalence around 41% to 64%. Low-probability observations within those species do not constitute a rare-species validation. The package's spatial-variance prior remains unchanged; the reference calculation alone cannot isolate its contribution to the smaller package-reference discrepancy.

## Complete-model checks

Three independent binary datasets were fitted at interior range-grid indices 4, 6 and 8. Each has 80 unique locations, 60 observations per location, eight species, 80 support points, one environmental covariate, nonzero intercepts and no latent factors or traits. These are independently generated binary observations with distinct site identities sharing coordinates, not repeated PCR detections of one fixed occupancy state. Each fit estimates the regression coefficients, spatial range and spatial SD, using two chains with 300 burn-in and 600 retained iterations. The priors are unchanged. The different range cells use independently seeded datasets; there is one dataset per range, not three replicates of each range.

Giving each dataset equal weight, the probability results are:

- Low probabilities: average error +2.51 percentage points; between-dataset standard error 0.53 points. Individual dataset errors range from +1.93 to +3.57 points.
- Medium probabilities: average error -0.081 points; between-dataset standard error 0.080 points. Individual errors range from -0.236 to +0.031 points.
- High probabilities: average error -2.38 points; between-dataset standard error 0.38 points. Individual errors range from -2.94 to -1.65 points.

All three average errors and all nine dataset-specific group errors lie within the unchanged five-point target. Five narrower probability bins and the very-low/very-high probability subsets also stay within five points in each dataset. The closest subset is true probability above 0.95 in the longest-range dataset: -4.44 points, with Monte Carlo standard error 0.062 points. None of these datasets contains deliberately rare species.

![Average probability errors in three independent informative datasets, with the provisional five-point target shaded](spatial-beta-probability-bias.png)

The between-dataset standard errors use only three independent datasets and include variation between spatial ranges. They do not quantify uncertainty at each range or demonstrate general calibration. Numerical uncertainty from the MCMC draws is much smaller: the combined Monte Carlo errors for the main groups are 0.012, 0.005 and 0.015 points. Their rank-normalized R-hat values are below 1.005, and the smallest effective sample size for a mean is 228. The largest difference between chain means is 0.037 points. Independent reconstruction agrees with the saved probability and field means to numerical precision. No fitting warnings were recorded.

Individual predictions can still be more than five points wrong. The probability RMSE, which measures error size without allowing positive and negative errors to cancel, is 5.63, 5.71 and 4.94 percentage points in the three datasets. The agreed target concerns average systematic error within each probability group, so this is a different measure. [Compact binary results and their column definitions](spatial-beta-binary-results/README.md) include the group errors, parameter estimates, simulation seeds and source fingerprints.

Spatial range means are 0.1063, 0.1676 and 0.2356 against truths 0.1067, 0.1711 and 0.2356. Spatial SD estimates are 0.997, 0.940 and 0.985 against a true value of 1. The posterior mean fields remain attenuated: their estimated-on-true slopes are 0.907, 0.883 and 0.887. This attenuation must remain visible even though the probability target is met. The known-parameter reference also attenuates individual field estimates; agreement with the probability target does not imply perfect recovery of every latent spatial effect.

## Reproduction and remaining review

The original continuous checks with unchanged priors remain in the [original validation](spatial-range-validation.md). The separate [noise-prior report](https://github.com/AlexDiana/occJSDM/blob/codex/continuous-noise-prior/dev/simstudy/continuous-noise-prior-validation.md) contains the small-noise experiments and their reproduction commands. Their improved results rely on the optional noise-prior change and must not be attributed to PR #8 alone.

The binary runner and reference save input data, truth, seeds, settings, source/library fingerprints and complete fits. These commands reproduce the existing work; reorganising the PR did not execute them again.

```sh
Rscript dev/simstudy/validate_spatial_binary.R --source=/path/to/compiled-source --out=/path/to/binary-results/package60-range6-full --grid-index=6 --repeats=60 --knots=80 --burn=300 --iter=600 --chains=2
Rscript dev/simstudy/summarise_spatial_binary.R --base=/path/to/binary-results --out=/path/to/new-binary-analysis
Rscript dev/simstudy/plot_spatial_beta_bias.R
Rscript dev/simstudy/verify_binary_oracle.R --out=/path/to/new-oracle-check
Rscript dev/simstudy/validate_binary_oracle.R --input=/path/to/binary-pilot/data-truth.rds --out=/path/to/new-oracle --burn=2000 --iter=4000 --thin=3
```

Use binary grid indices 4, 6 and 8 in separate sibling output directories named `package60-range4-full`, `package60-range6-full` and `package60-range8-full`. The 60-observation reference pilot used four chains, 2,000 burn-in iterations and 4,000 retained draws with thinning by three. The saved binary source has unchanged binary priors, complete spatial support and grouped arithmetic. The later approved residual-correlation correction does not change these no-factor, no-trait fits. Raw files remain in the task's `work/spatial-binary-resume-20260913` directory; compact CSVs do not replace them.

Alex still needs to review the spatial code, the support-point advice and these results. Rare species and other sampling designs need separate assessment. The small-noise prior/default decision belongs to the separate noise PR. The consolidated [release TODO](https://github.com/AlexDiana/occJSDM/blob/codex/nonspatial-bias-recheck/TODO.md) records all remaining work. Passing these three informative binary cases does not clear every beta gate.

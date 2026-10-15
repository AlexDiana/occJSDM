# Interval calibration of the 15 October 2026 paper study

This is the interval-calibration assessment that the *Paper simulation study* item in `TODO.md` asks for. It uses the R = 300 run recorded in [README.md](README.md): 3,000 fits of the ten production cells on the combined, unmerged code of PRs 31 and 32. No package code changed. Every number below comes from [calibration.R](calibration/calibration.R) or [bisect.R](calibration/bisect.R) and is in one of the CSV files beside them.

## Summary

- **The uncertainty is smaller than the runner suggests.** Coverage standard errors computed across replicates are 0.1 to 0.9 percentage points per block and cell, against the 1.3 points the runner quotes. The elements of one block in one replicate turn out to be only weakly correlated (design effects 0.7 to 1.8), so at R = 300 a deviation of about one point from 0.95 is detectable.
- **Close to nominal:** `q` (0.939 to 0.949) and `G` (0.938 to 0.956, apart from 0.912 in `species_20`), and `B0` in `base` (0.942).
- **Too wide:** `beta_theta` (0.960 to 0.967) and `theta0` (0.974 to 0.983) overcover in every cell.
- **Too narrow:** `B` undercovers in every cell (0.875 to 0.927), `B0` in eight of nine (0.929 to 0.937, and 0.868 in `binary`), and `p` in all seven that estimate it (0.908 to 0.921).
- `sigma_b` and `resid_cor` say nothing about calibration, for reasons already recorded.
- `low_information` is poorly calibrated in every informative block.
- **The comparison with 29 July is paired.** Replicates 1 to 100 have the same simulated data and the same truth in both runs; only the code and the MCMC stream differ.
- **`beta_theta`** rose from about 0.75 to 0.96 because of the September collection-covariate alignment fix, not because of the change of scoring units.
- **`B`, and `B0` in `binary`,** lost coverage because their intervals narrowed by 10 to 29% while the estimates did not improve. Paired refits put the change at one commit, `f1d2e12` (11 September), which by code reading made the spatial part of the sampler coherent: before it, the spatial range was stuck at its grid maximum and the spatial amplitude was inflated about fivefold. So July's near-nominal figures for these blocks were not evidence of calibration. Why the corrected model undercovers is not established; the leading candidate is that the simulated spatial range lies beyond the fitted range grid in every spatial cell.

## How coverage was assessed

For each cell and block, each replicate contributes the fraction of its elements whose 95% interval contains the truth. Coverage is the mean of those 300 fractions and its standard error is their standard deviation over the square root of 300. Replicates are independent, because each seeds itself from its cell and replicate number; the elements within a replicate are not, because they share one dataset and one MCMC run (`dev/simstudy/PLAN.md` section 9). The 95% intervals quoted are the estimate plus or minus 1.96 standard errors.

The runner's figure of about 1.3 points is `sqrt(0.95 * 0.05 / R)`, which treats each replicate as a single check. The design effect, the variance of the replicate-based estimate over the variance that independent elements would give, is 0.7 to 1.8 across the 73 cell-block pairs where it is defined (`sigma_b` always covers, so it has none). So the effective number of checks is close to the number of elements, not the number of replicates, and the standard errors are 0.1 to 0.9 points. Every coverage below carries its replicate-based standard error in brackets.

Alongside coverage the CSVs report the mean interval width, the bias and RMSE of the posterior mean, the slope of the posterior mean on the truth (1 means no shrinkage), and an error-to-width ratio: the root mean square of the error divided by the posterior SD that the interval implies (its width over 3.92). That ratio is close to 1 for a calibrated, roughly normal posterior; above 1 means the intervals are narrow for the errors the estimates actually make.

What this run cannot answer:

- **High contamination.** Every production cell draws the nominal `q` from 0.01 to 0.05 and `theta0` from 0.02 to 0.10, both well inside the central 95% interval of their Beta(1, 20) priors (0.001 to 0.168). The whole grid is therefore a low-contamination assessment, which is what the TODO item asks to keep separate. There is no high-contamination cell; the prior-mismatch stress tests remain those in the [targeted interval follow-up](../interval-calibration/REPORT.md) and the [theta0 prior study](../theta0-prior-mean/REPORT.md).
- **The historical spatial continuous configuration.** The `continuous` arm of `PLAN.md` section 16.5 (`B0` coverage 0.879) is not a production cell, so this run says nothing about it. It stays open.
- **Convergence.** The results file keeps interval endpoints and posterior means but no R-hat or effective sample sizes, so no fit's convergence can be audited here. With 2 chains of 1,000 draws the interval endpoints also carry Monte Carlo noise, which biases coverage slightly downward (`PLAN.md` section 8).
- **Spatial range and amplitude.** `l_s` and `sigma_bs` are not scored, as before, so no cell says anything about spatial range.

## Coverage at R = 300, nine cells

These are the nine cells other than `low_information`. Cells that do not estimate a block are omitted from its line.

- `B0`: `base` 0.942 (0.004), `spatial_isolated` 0.935 (0.005), `traits_isolated` 0.933 (0.005), `primers_3` 0.937 (0.004), `d_underfit` 0.929 (0.005), `d_overfit` 0.933 (0.005), `species_20` 0.930 (0.003), `occupancy` 0.929 (0.005), `binary` 0.868 (0.006). Mean width 1.149 to 3.110.
- `B`: `base` 0.927 (0.004), `spatial_isolated` 0.875 (0.005), `traits_isolated` 0.909 (0.004), `primers_3` 0.917 (0.004), `d_underfit` 0.901 (0.005), `d_overfit` 0.914 (0.004), `species_20` 0.884 (0.004), `occupancy` 0.916 (0.004), `binary` 0.902 (0.004). Mean width 1.357 to 3.202.
- `G`: `base` 0.955 (0.007), `traits_isolated` 0.946 (0.007), `primers_3` 0.951 (0.006), `d_underfit` 0.938 (0.007), `d_overfit` 0.940 (0.007), `species_20` 0.912 (0.009), `occupancy` 0.956 (0.006), `binary` 0.951 (0.007). Mean width 1.285 to 1.920.
- `beta_theta`: `base` 0.963 (0.003), `spatial_isolated` 0.967 (0.002), `traits_isolated` 0.963 (0.003), `primers_3` 0.962 (0.003), `d_underfit` 0.965 (0.002), `d_overfit` 0.963 (0.002), `species_20` 0.960 (0.002), `occupancy` 0.961 (0.002). Mean width 1.491 to 1.821.
- `theta0`: `base` 0.976 (0.003), `spatial_isolated` 0.981 (0.002), `traits_isolated` 0.975 (0.003), `primers_3` 0.978 (0.003), `d_underfit` 0.983 (0.002), `d_overfit` 0.974 (0.003), `species_20` 0.978 (0.002), `occupancy` 0.979 (0.003). Mean width 0.132 to 0.137.
- `p`: `base` 0.911 (0.004), `spatial_isolated` 0.915 (0.003), `traits_isolated` 0.909 (0.004), `primers_3` 0.921 (0.003), `d_underfit` 0.908 (0.004), `d_overfit` 0.908 (0.004), `species_20` 0.911 (0.003). Mean width 0.190 to 0.222.
- `q`: `base` 0.945 (0.003), `spatial_isolated` 0.949 (0.003), `traits_isolated` 0.946 (0.003), `primers_3` 0.946 (0.002), `d_underfit` 0.947 (0.003), `d_overfit` 0.939 (0.003), `species_20` 0.948 (0.002). Mean width 0.031 to 0.036.
- `sigma_b`: `base` 1.000 (0.000), `spatial_isolated` 1.000 (0.000), `traits_isolated` 1.000 (0.000), `primers_3` 1.000 (0.000), `d_underfit` 1.000 (0.000), `d_overfit` 1.000 (0.000), `species_20` 1.000 (0.000), `occupancy` 1.000 (0.000), `binary` 1.000 (0.000). Mean width 0.213 to 0.219.
- `resid_cor`: `base` 0.760 (0.004), `spatial_isolated` 0.761 (0.003), `traits_isolated` 0.756 (0.004), `primers_3` 0.766 (0.004), `d_underfit` 0.979 (0.001), `d_overfit` 0.758 (0.004), `species_20` 0.751 (0.002), `occupancy` 0.757 (0.004), `binary` 0.752 (0.004). Mean width 1.752 to 1.993.

Reading across blocks:

- **`B0`** looks close to nominal in the occupancy and two-stage cells, but that is partly by construction. `B0` is the one block whose truth is drawn from exactly its fitted prior, Normal(0, 1), so its average coverage would be near nominal even if the data said little. They do say little: in those eight cells the slope of the posterior mean on the truth is 0.32 to 0.39, the intervals are 74 to 79% as wide as the prior's central 95% interval (3.92), and in replicates 1 to 100 truths beyond 2 in absolute value are covered about half the time or less. In `binary` the slope is 0.80 and the intervals are under a third of the prior's width, so there the data dominate, and coverage is 0.868.
- **`B`** undercovers everywhere, and the shortfall is concentrated in large true effects. Pooled over the seven cells that fit both traits and a spatial field (replicates 1 to 100), elements whose true value exceeds 3 in absolute value are covered 0.733 of the time, while those below 1 are covered 0.95. Their posterior means are pulled towards zero by 1.27 on average. The error-to-width ratio is 1.12 to 1.39.
- **`beta_theta`** overcovers in both rows: intercepts 0.963 to 0.972 and slopes 0.954 to 0.961, with error-to-width ratios of 0.85 to 0.90. The slope overconfidence that `PLAN.md` sections 13 to 15 chased is gone (see the July comparison below).
- **`p`** undercovers by 3 to 4 points with a bias of +0.03 to +0.04, as in July: the informative detection prior pulls the estimates up.
- **`G`** is within 2 SE of nominal in every cell except `species_20` (0.912, 4.3 SE low).
- **`q`** is within 2 SE of nominal in six of seven cells; `d_overfit` is 0.939 (3.5 SE low).

## `low_information`, reported separately

`low_information` has 40 sites rather than 100, true detection `p` from 0.1 to 0.3 and baseline collection probability from 0.05 to 0.2, and 8 support points. It changes three things at once, so nothing in it is attributable to one cause.

- `B0`: 0.845 (0.006), width 4.220, bias -1.189, RMSE 1.575.
- `B`: 0.799 (0.006), width 4.398, bias +0.009, RMSE 1.623.
- `G`: 0.838 (0.011), width 2.392, bias +0.006, RMSE 0.798.
- `beta_theta`: 0.995 (0.001), width 3.857, bias +0.153, RMSE 0.649.
- `theta0`: 0.604 (0.009), width 0.068, bias -0.041, RMSE 0.047.
- `p`: 0.103 (0.004), width 0.609, bias +0.492, RMSE 0.512.
- `q`: 0.867 (0.004), width 0.049, bias +0.011, RMSE 0.016.
- `sigma_b`: 1.000 (0.000), width 0.206, bias -0.025, RMSE 0.025.
- `resid_cor`: 0.753 (0.004), width 1.993, bias +0.001, RMSE 0.709.

Against the other occupancy and two-stage cells, `B`, `B0` and `G` lose 7 to 13 points, `theta0` loses 37, `p` collapses (its truth sits deep in the tail of the informative detection prior, `PLAN.md` section 12.1), `q` loses 8, and `beta_theta` overcovers further because its intervals are more than twice as wide. Users with small or low-detection datasets remain the most exposed, as recorded in July.

## `theta0` interval widths

- In the eight cells other than `low_information` the mean `theta0` interval width is 0.132 to 0.137, about 80% of the width of the Beta(1, 20) prior's central 95% interval (0.167). Every true `theta0` lies inside that prior interval, so overcoverage of 0.974 to 0.983 is what a mostly prior-driven interval produces here. The bias is small and negative, -0.001 to -0.007.
- These widths agree with the targeted theta0 study's 0.130 for its baseline design with the same prior.
- Paired against July (replicates 1 to 100), the intervals are 6 to 8% narrower (width ratios 0.920 to 0.939). Coverage moved by -0.010 to +0.001, within 2 SE of zero in every cell except `d_overfit` (-0.010, SE 0.005).
- In `low_information` the width is 0.068 with coverage 0.604 and bias -0.041: narrower intervals around estimates that are too low.

## Comparison with the 29 July 2026 baseline

### The comparison is paired

Replicates 1 to 100 of every cell share their simulated data and truth with the July run (`simstudy-20260729-143756.rds`). The evidence:

- The truths of `B0`, `B`, `theta0`, `sigma_b` and `resid_cor` match element by element: exactly for `sigma_b` and `resid_cor`, and to within 2e-15 for the other three, which is rounding from the different BLAS ([pairing-check.csv](calibration/pairing-check.csv)).
- The truths of `beta_theta`, `G`, `p` and `q` differ only by the harness's change of scoring units since July: collection slopes are now expressed per fitted standard deviation, trait coefficients per trait standard deviation, and `p` and `q` as probabilities of a recorded read. The elements those transformations leave unchanged, such as zero slopes, match exactly.
- The simulator is unchanged since July. `simulateOccJSDMData()`, `simulateData()` and every R function they call deparse identically at July's commit `867fc53` and now, and the C++ kernel `K2()` differs only in a commented-out pragma. `draw_truth()` gained an optional `q_range` that the production cells do not set.

So each July-to-October difference below is a difference in code (and MCMC stream) on the same datasets, and its standard error is computed from the per-replicate differences.

The 2 August and 10 August R = 100 runs on the same seeds give the same B and B0 coverage and widths as July to within noise (for example `binary` `B0`: 0.962, 0.958, 0.958 with widths 1.619, 1.616, 1.616). So the changes discussed below entered between 10 August and 15 October.

### `beta_theta`: the collection-covariate alignment fix

Coverage rose by 0.19 to 0.26 in the eight cells other than `low_information` (which rose by 0.01). Rescoring July's own intervals against the corrected truth changes July's coverage by at most 0.005, so the scoring change explains none of it. The data show what did: in July the nonzero slopes, whose true values are about -1 and +1 on the fitted scale, were estimated at -0.106 and +0.104 on average and covered 0.33 of the time; now they are estimated at -1.088 and +1.076 and covered 0.965. Zero slopes covered 0.94 in both runs, and intercepts 0.953 in July and 0.969 now. That is the attenuation that [the collection-alignment validation](../collection-alignment-validation.md) traced to collection covariates being matched to the wrong samples, fixed in September (`90026d6` and `d88f52e`). The July "slope overconfidence" was this bug, not a property of the Polya-Gamma update.

### `B`

What the data show, paired over replicates 1 to 100 ([paired-july-r100.csv](calibration/paired-july-r100.csv), [shrinkage-b-b0.csv](calibration/shrinkage-b-b0.csv)):

- Coverage fell in the seven cells that fit both observed traits and a spatial field: `base` -0.011 (SE 0.006), `primers_3` -0.029 (0.007), `d_underfit` -0.039 (0.008), `d_overfit` -0.028 (0.006), `species_20` -0.045 (0.006), `occupancy` -0.022 (0.007) and `binary` -0.047 (0.006). `low_information`, which also has both, fell by 0.035 (0.011).
- It did not move in the two cells that have only one of the two: `spatial_isolated` (spatial field, traits not fitted) +0.002 (0.006) and `traits_isolated` (traits fitted, no spatial field) -0.003 (0.006).
- In the seven cells the intervals are 10 to 22% narrower (width ratios 0.78 to 0.90), against 9% in `spatial_isolated` and 5% in `traits_isolated`.
- The estimates did not get worse on average: RMSE changed by at most 5.4% either way and the mean error stayed near zero. But they are shrunk harder towards zero. The slope of posterior mean on truth fell in every one of the seven cells, from between 0.74 and 0.98 in July to between 0.66 and 0.86 now (`binary` 0.98 to 0.86, `species_20` 0.82 to 0.67), and did not change in `spatial_isolated` (0.55 to 0.57) or `traits_isolated` (0.66 to 0.66).
- So the error-to-width ratio rose from 0.97 to 1.04 in July to 1.11 to 1.29 now, and, as above, the lost coverage sits in the large true effects.

### `B0`

What the data show, on the same paired replicates:

- In `binary` coverage fell from 0.962 to 0.860 (paired -0.102, SE 0.010). The intervals are 29% narrower (1.62 to 1.16) while the estimates are no better: mean error +0.012 and +0.005, RMSE 0.391 and 0.386, slope on truth 0.84 and 0.81. The error-to-width ratio went from 0.95 to 1.34.
- In the eight occupancy and two-stage cells other than `low_information` the intervals are 8 to 14% narrower, and the estimates did improve: the mean error moved a quarter to a half of the way towards zero (`base` -0.208 to -0.138, `primers_3` -0.151 to -0.088), RMSE fell 4 to 8%, and the slope on truth rose from between 0.25 and 0.31 to between 0.31 and 0.40. Coverage fell slightly, by 0.007 to 0.030; the fall exceeds 2 SE in `d_underfit` (-0.030, SE 0.007), `species_20` (-0.020, 0.005), `primers_3` (-0.015, 0.006) and `occupancy` (-0.018, 0.008).
- `low_information` did not change (-0.005, SE 0.008).

### Other blocks

- `G`: paired changes of -0.020 to +0.010, all within 2 SE of zero. July's `G` intervals are scored against July's raw-unit truth because July's fit used raw traits.
- `p` and `q`: once July is rescored in the same units, the paired changes are -0.004 to +0.006 for `p` and -0.003 to +0.011 for `q`. Three of the sixteen are a little over 2 SE, about what chance gives among this many comparisons. The apparent fall of `q` in `low_information` from 0.911 to 0.862 is the change of units: July rescored is 0.851.
- `theta0`: see above.
- `resid_cor` and `sigma_b`: unchanged, as their construction implies.

## Which change did it: bisection refits

This check was added to the assessment because history alone left three candidates (the collection-alignment fix, trait standardisation in `9a784a9`, and the spatial pull request, PR 8). `base` and `binary` were refitted at replicates 1 to 20 with the study's runner at six commits between the 10 August run and this one ([run-bisect.sh](calibration/run-bisect.sh)). The runner seeds every replicate exactly as in production, so each refit sees the same datasets as July and October; only the code and the MCMC stream differ. At R = 20 a coverage difference of a few points is noise, but interval widths and the slope of posterior mean on truth are paired and precise. Results are in [bisect-base-binary.csv](calibration/bisect-base-binary.csv), all scored against the October truth.

- **`5382460`** (main, 15 September: the collection-alignment fix and the serial random-number steps). `binary` is bit-identical to 10 August (B0 width 1.605, B slope 0.976). In `base`, `beta_theta` is repaired (coverage 0.78 to 0.975, slope on truth 0.41 to 1.09), the `B0` error falls from -0.262 to -0.174 and its width from 3.50 to 3.30, and `B` narrows from 3.59 to 3.42 while its slope improves from 0.73 to 0.77 and its coverage stays at 0.948.
- **`9a784a9`** (main, 16 September: traits standardised before fitting). No material change in either cell (`binary` B0 width 1.597, B slope 0.973).
- **`54b5cca`** (base of the spatial branch, with the alignment fix and serial steps but no spatial change). Identical to `5382460` in both cells.
- **`f1d2e12`** (spatial branch, 11 September: "Align spatial range updates with the fitted coefficient model"). This is where the change happens. In `binary`, the B0 width falls from 1.60 to 1.16 and its coverage from 0.950 to 0.845, and B narrows from 1.53 to 1.38 with its slope falling from 0.98 to 0.87 and coverage from 0.945 to 0.892. In `base`, B narrows from 3.42 to 3.22, its slope falls from 0.77 to 0.68 and coverage from 0.948 to 0.912, and B0 narrows from 3.30 to 3.12.
- **`e9fd2d3`** (blocked range and coefficient draws), **`cc1b4e6`** (PR 8 merged into main) and **this run** (adding everything since, including PRs 31 and 32): the same as `f1d2e12` to within noise (`binary` B0 width 1.14 to 1.16, B slope 0.86 to 0.87; `base` B width 3.09 to 3.21).

So the `B` change in the cells that fit traits and a spatial field, and the `B0` change in `binary`, both entered with one commit, `f1d2e12`. The smaller `B0` narrowing in the two-stage cells came in two steps, about half with the alignment fix and half with `f1d2e12`. Nothing after the PR 8 merge, PRs 31 and 32 included, moved these blocks.

`f1d2e12` changed three things at once (code reading):

1. **The spatial basis.** Each site's field used to be built from the whitened basis columns of only its five nearest support points (`maxPoints = 5`); it now uses all of them, which is what the subset-of-regressors covariance requires.
2. **The prior mean of the spatial coefficients.** The joint draw of `B0`, `B`, `L` and the spatial coefficients `Bs` used a prior mean of zero for `Bs`, while the updates of the spatial amplitude `sigma_bs` and of the spatial trait loadings `Gs`, `As` and `Cs` treated the trait regression `Tr Gs + As Cs` as that mean. The draw now uses the same mean, in both the R and the C++ samplers.
3. **The range update.** A neighbouring-grid Metropolis step that scored a full-GP density of the reconstructed field (with the amplitude fixed at 1) became an exact categorical draw from the conditional likelihood.

A small diagnostic shows what this did to the spatial term ([spatial-diagnostic-binary.csv](calibration/spatial-diagnostic-binary.csv)). Five datasets were simulated with the `binary` cell's settings (not the production replicates, because this script did not reproduce the runner's random-number setup) and fitted at `54b5cca` and at `f1d2e12`:

- Before: the range sat at the top of its grid (0.30) in every retained draw of every fit, and `sigma_bs` was 1.50 to 1.64, about five times its prior mean and far outside the prior's central 95% interval (0.24 to 0.46). The posterior SD of `B0` averaged 0.40 to 0.44.
- After: the range moves (posterior means 0.02 to 0.26), `sigma_bs` sits at 0.31 to 0.34, the prior mean, and the posterior SD of `B0` averages 0.28 to 0.37.

A second check tested whether the amplitude prior is responsible, using the package's experimental half-Cauchy amplitude prior (`listPriors$sigma_bs_prior = "half_cauchy"`) on the study's code, `base` and `binary`, replicates 1 to 20, paired with the production datasets ([halfcauchy.R](calibration/halfcauchy.R)). It changed nothing that matters: B width 3.24 against 3.20 in `base` and 1.38 against 1.37 in `binary`, `binary` B0 width 1.165 against 1.164, and coverage within noise. The refits did not save `sigma_bs`, so whether that prior moved the amplitude is not shown.

## Best-supported explanation, and its limits

What the data establish:

- The July and October runs fit the same datasets, so their differences are code differences.
- `beta_theta`'s recovery is the collection-alignment fix: July's nonzero slopes were estimated near zero.
- The narrower `B` intervals, the stronger shrinkage of large `B` effects and the narrower `binary` `B0` intervals all arrive with `f1d2e12`, and with nothing else between 10 August and this run.
- Before that commit the spatial range was stuck at the top of its grid and the spatial amplitude was inflated about fivefold; after it, both behave differently and the `B0` posterior is about a quarter narrower on the same data.
- The experimental half-Cauchy amplitude prior does not undo the change.
- The simulated spatial range, 0.15 on the raw unit square, is 0.48 to 0.54 on the standardised coordinates the fitter uses (checked on five datasets with the `binary` cell's settings). The fitted range grid stops at 0.30, so no fit in the nine spatial cells can represent the true smoothness. The spatial design sweep chose its range (about 0.10 standardised) to sit inside the grid; the production grid did not.

What is inferred:

- **July's near-nominal `B` and `binary` `B0` coverage was not evidence of calibration.** Before `f1d2e12` the spatial block was not a coherent sampler for the stated model: the coefficient draw and the amplitude and trait-loading draws assumed different priors for the same coefficients, and the range update scored the wrong density. The inflated amplitude gave the spatial field a large, loose prior, which is consistent with the wider intervals of everything the field competes with, `B0` most of all. Those intervals happened to cover.
- **The October intervals are the stated model's.** Their undercoverage of `B` and of `B0` in `binary` therefore describes how the stated model fits this grid, not a new sampler bug introduced since July. This rests on the spatial tests and reports (`spatial-range-validation.md`, `spatial-beta-recheck.md`), which check the mathematics; nothing in this run verifies the sampler.
- **The leading untested candidate is the range mismatch.** With the true range beyond the grid, the fitted field cannot reproduce the simulated spatial variation, and in a logistic model unmodelled variation attenuates the coefficients and makes their intervals too narrow. That matches the pattern: more shrinkage of large `B` effects, narrower intervals, unchanged average error. Before `f1d2e12` the range sat at 0.30, the grid value nearest the truth, with an inflated amplitude.

Limits:

- `f1d2e12` bundles the three changes above, and the refits do not separate them.
- The refits cover `base` and `binary` only, at R = 20; the other cells are attributed by their shared pattern.
- The traits-and-spatial pattern in `B` rests on two contrast cells. `spatial_isolated`, the only spatial cell without fitted traits, also narrowed (9%) but held its coverage while its RMSE improved by 4%; whether fitted traits are needed for the effect, perhaps through the trait-structured spatial mean, is not established.
- No fit here was checked for convergence, and the spatial range and amplitude were not recorded in the production run.

## What Doug needs to decide

- **What the paper may claim.** On this grid, at low contamination: `q` and `G` are close to nominal; `beta_theta` and `theta0` are conservative; `B` is 2 to 8 points short and much worse for large effects; `B0` is 1 to 2 points short in the occupancy and two-stage cells (and mostly prior there) and 8 points short in `binary`; `p` is 3 to 4 points short. Doug decides whether that is reported as it stands, with the caveats above, or whether the `B` question is settled first.
- **Whether to test the range explanation before the paper.** Two ways, each well under an hour on the laptop at R = 50 for `base` and `binary`: simulate with a range inside the fitted grid (for example 0.05 raw, about 0.17 standardised), which needs only a new harness cell; or extend the fitted range grid, which is a package change. The first is cheaper and changes no package code. A third, separating the three parts of `f1d2e12`, needs patched copies of the package and is mainly of interest to Alex.
- **Whether the fixed range grid is a user-facing limitation.** The grid stops at 0.30 in standardised coordinates, where the correlation one coordinate SD apart is about 0.004, so any field smoother than that cannot be represented. That is a modelling question for Alex, and possibly a documentation item.
- **The runner's uncertainty statement.** The runner and `PLAN.md` section 9 quote `sqrt(0.95 * 0.05 / R)`, about 1.3 points at R = 300; the replicate-based standard errors are 0.1 to 0.9 points. Updating the runner's message and section 9 is a small change outside the package.
- **What to save next time.** The production runner keeps no R-hat, effective sample sizes, range or amplitude draws. Saving them would let a run like this audit convergence and the spatial term directly.
- **The two TODO sub-items this run cannot answer.** High-contamination stress tests and the historical spatial continuous configuration need their own cells or targeted runs, or an explicit decision that the targeted reports suffice.
- **The validation article.** The TODO item says to regenerate it from the new results. That is now possible for these numbers, but the article should carry the `B` caveat and the fact that the run is on the unmerged code of PRs 31 and 32. The refits show those two PRs do not move `B`, `B0`, `beta_theta` or `theta0`, but if either changes in review the run's provenance no longer matches merged code.
- **`resid_cor`.** It still cannot be quoted until the statistic is repaired (`PLAN.md` section 17.6).

## Files and reproduction

- [calibration/calibration.R](calibration/calibration.R) writes `coverage-r300.csv` (every cell and block at R = 300), `coverage-r300-nine-cells.csv`, `beta-theta-rows-r300.csv`, `pairing-check.csv`, `paired-july-r100.csv`, `shrinkage-b-b0.csv` and `august-10-b-b0.csv`. Run it from the repository root, giving it the July and 10 August results files from the main checkout's `dev/simstudy/results/`, which are not committed.
- [calibration/run-bisect.sh](calibration/run-bisect.sh) ran the bisection refits with the study's runner, on 15 October 2026 under R 4.6.1 with Accelerate BLAS, 5 replicate processes and one thread per fit, on mains power. [calibration/halfcauchy.R](calibration/halfcauchy.R) ran the half-Cauchy refits and [calibration/spatial-diag.R](calibration/spatial-diag.R) the spatial diagnostic, each from a `git archive` export of the commit concerned. Their paths point at a session scratch directory, and their raw results are not committed.
- [calibration/bisect.R](calibration/bisect.R) summarises the refits into `bisect-base-binary.csv`; the diagnostic's summary is `spatial-diagnostic-binary.csv`.

Both summary scripts run from the repository root. `MAIN` stands for the main checkout:

```sh
Rscript dev/simstudy/paper-study-20261015/calibration/calibration.R MAIN/dev/simstudy/results/simstudy-20260729-143756.rds MAIN/dev/simstudy/results/simstudy-20260810-121954.rds
```

```sh
Rscript dev/simstudy/paper-study-20261015/calibration/bisect.R /path/to/bisect-directory MAIN/dev/simstudy/results/simstudy-20260729-143756.rds MAIN/dev/simstudy/results/simstudy-20260810-121954.rds
```

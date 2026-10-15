# The spatial 2x2: does the simulated spatial field cause B's undercoverage?

The R = 300 paper study found the B intervals too narrow in every cell, worst for the largest true effects ([calibration report](../paper-study-20261015/REPORT.md)). Doug's hypothesis was that the simulated spatial field, which the fitted spatial term cannot fully represent, attenuates B. This run tested that on the base cell's settings by crossing whether a field is simulated with whether a spatial term is fitted. The design and the pairing are in [README.md](README.md).

## Summary

-   **Not at the strength the study simulates.** Leaving a simulated field unmodelled changed B coverage by +0.009 (SE 0.010) and the slope of the estimates on the truth from 0.674 to 0.666. There is no measurable attenuation.
-   **Fitting a spatial term raises B coverage** by 1.3 to 2.1 points, whether or not a field is present, by widening B's intervals about 4%. The difference between the two, the interaction, is +0.008 (SE 0.008) and cannot be told from zero.
-   **B undercovers most when there is no spatial component anywhere:** 0.906 (SE 0.008) in `no_spatial`. The shortfall for true effects larger than 3 in absolute value is present in all four arms, at 0.685 to 0.818. So the main cause of B's undercoverage is not spatial.
-   **B0** shows the same ordering, smaller, and no contrast exceeds two standard errors. Its coverage is 0.928 to 0.939 in all four arms.
-   **For the change since July,** this supports the calibration report's reading. The spatial term pads B's intervals. The sampler before `f1d2e12` padded them more, through its inflated spatial amplitude, and removing that padding exposed a shortfall that does not depend on the spatial field.

## Run

The three new arms ran on 15 October 2026 from 21:38 to 22:22 at revision `01c55a1`, with a clean tree at start and at write: 300 fits, no failures, 5 replicate processes, one thread per fit, on mains power with Low Power Mode off, under R 4.6.1 with Accelerate BLAS. `base` is replicates 1 to 100 of the R = 300 paper study at `de4389c`, whose package code is identical. Every replicate is present in all four arms.

## Coverage by arm

Replicates 1 to 100 in each arm. Standard errors are computed across replicates, as in the calibration report. The error-to-width ratio is near 1 for calibrated intervals and above 1 when they are narrow for the errors the estimates make.

-   `base`, field simulated and fitted: B 0.936 (0.006), width 3.22, RMSE 0.850, slope 0.706, error-to-width 1.11. B0 0.939 (0.008).
-   `field_unfitted`, field simulated, not fitted: B 0.915 (0.007), width 3.08, RMSE 0.850, slope 0.666, error-to-width 1.18. B0 0.930 (0.007).
-   `nofield_fitted`, no field, spatial term fitted: B 0.919 (0.006), width 3.18, RMSE 0.893, slope 0.707, error-to-width 1.15. B0 0.933 (0.008).
-   `no_spatial`, neither: B 0.906 (0.008), width 3.06, RMSE 0.889, slope 0.674, error-to-width 1.20. B0 0.928 (0.008).

## Paired contrasts

Each contrast is the mean of per-replicate differences in coverage, with its standard error. The first two compare fits of identical datasets. The last two compare different datasets drawn from identical truths.

-   **Fitting a spatial term, field present** (`base` minus `field_unfitted`): B +0.021 (0.006), intervals 4% wider. B0 +0.009 (0.006).
-   **Fitting a spatial term, no field** (`nofield_fitted` minus `no_spatial`): B +0.013 (0.006), intervals 4% wider. B0 +0.005 (0.005).
-   **The field, spatial term fitted** (`base` minus `nofield_fitted`): B +0.017 (0.008). B0 +0.006 (0.008).
-   **The field, no spatial term** (`field_unfitted` minus `no_spatial`): B +0.009 (0.010). B0 +0.002 (0.008).
-   **Interaction:** B +0.008 (0.008). B0 +0.004 (0.007).

Fitting a spatial term changes the B estimates themselves very little: their slope on the truth moves by 0.03 to 0.04 and their RMSE not at all. Its effect on coverage comes through the width of the intervals.

## B by size of the true effect

Coverage of B elements grouped by the absolute size of the true effect, 2,000 elements per arm.

-   Under 1 (917 elements): `base` 0.972, `field_unfitted` 0.961, `nofield_fitted` 0.953, `no_spatial` 0.948.
-   1 to 2 (605): 0.940, 0.927, 0.929, 0.917.
-   2 to 3 (297): 0.889, 0.865, 0.899, 0.889.
-   Over 3 (181): 0.818, 0.724, 0.746, 0.685.

The large effects are under-covered in every arm, including the one with no spatial component. Their estimates are pulled towards zero by 1.30 to 1.46 on average, again in every arm.

## Pairing

The truths of B, B0, G, `theta0`, `p` and `q` are identical, element for element, in all four arms. `field_unfitted` also shares `base`'s datasets, and `nofield_fitted` shares `no_spatial`'s. The collection slopes in `beta_theta` differ in the two no-field arms, because the simulator draws them after the field; `beta_theta` is not compared here. The harness tests in `tests/testthat/test-simstudy-spatial-2x2.R` check the pairing at replicate 1, and [pairing-check.csv](pairing-check.csv) checks it at all 100.

## What this does not show

-   **One design and one field.** Only the base cell's settings were run, with the field the study always simulates: 2 spatial factors, scale 0.5, range 0.15 on the raw unit square. A stronger or smoother field might matter more. The `binary` cell, where B0 dropped most since July, was not run.
-   **Small effects.** At 100 replicates a contrast of about 2 points is detectable and one of 1 point is not.
-   **Why an unmodelled field does so little.** One plausible reason, not tested here, is that the model's two latent factors are site-level and can absorb part of a site-level field. The study's `d_underfit` cell, where the truth has 4 factors and the fit 2, loses B coverage by a similar amount.
-   **An unexplained difference in error.** B's RMSE is about 0.04 higher in the two no-field arms than in the two field arms (per-replicate difference +0.041, SE 0.018, and +0.038, SE 0.018), while fitting a spatial term leaves it unchanged. That is about two standard errors in one of several comparisons, so it may be chance.
-   **The range question.** Whether a range inside the fitted grid changes the picture is being tested separately, on branch `claude/range-inside-grid`.

## What this changes

-   **The paper's B caveat stands,** but it should not be put down to the spatial field. B undercovers without any spatial component, and most in that case.
-   **The next candidates are non-spatial:** the latent factors competing with the covariates, the hierarchical prior on B with `sigma_b` held near its prior mean, and Monte Carlo error from 2 chains of 1,000 draws. Each can be tested as an arm paired with `no_spatial`: longer chains, no latent factors in the truth and the fit, and traits left out.

## Files and reproduction

-   [analyse.R](analyse.R) writes [coverage-by-arm.csv](coverage-by-arm.csv), [paired-contrasts.csv](paired-contrasts.csv), [b-by-effect-size.csv](b-by-effect-size.csv) and [pairing-check.csv](pairing-check.csv) from this run's results file and the R = 300 paper study's.
-   `results/simstudy-20261015-222234.rds` holds the three new arms; `run.log` and `provenance.txt` record the run.

Run from the repository root:

```sh
Rscript dev/simstudy/spatial-field-2x2-20261015/analyse.R
```

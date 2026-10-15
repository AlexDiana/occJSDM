# Non-spatial candidates for B's undercoverage

The spatial 2x2 ([its report](../spatial-field-2x2-20261015/REPORT.md)) found B's 95% intervals covering worst when there is no spatial component at all: 0.906 (SE 0.008) in `no_spatial`, and 0.685 for true effects over 3 in absolute value. This run tests the three non-spatial candidates the 2x2 named, each as an arm paired with `no_spatial` on the base cell's settings, with no spatial field simulated or fitted: Monte Carlo error, latent factors competing with the covariates, and B's hierarchical prior. The design and the pairing are in [README.md](README.md).

## Summary

-   **B's hierarchical prior accounts for most of the shortfall.** On identical datasets, fitting 4 latent traits instead of 2, which loosens the prior on B, raised B coverage from 0.906 to 0.936 (paired +0.030, SE 0.006), about two-thirds of the 4.4-point shortfall. The error-to-width ratio fell from 1.20 to 1.02 and the slope of the estimates on the truth rose from 0.674 to 0.769. Coverage of true effects over 3 rose from 0.685 to 0.829. Leaving the traits out, which takes the measured-trait part out of B's prior mean, moved every one of these the other way: coverage 0.871 (paired -0.036, SE 0.007).
-   **Monte Carlo error is a smaller, real part, and it is unconverged sampling, not only noisy endpoints.** Running 5,000 + 5,000 iterations per chain instead of 1,000 + 1,000 raised B coverage to 0.923 (paired +0.017, SE 0.005) and widened the intervals by 2%. At the study's chain length B has not converged: in 10 refitted replicates its median bulk effective sample size is 57 of 2,000 draws, and 91% of B elements have R-hat over 1.01. At 10,000 draws the median is 230 and 41% still exceed 1.01, so this measures only part of what chain length does.
-   **Latent factors competing with the covariates: no.** Fitting 2 factors to data simulated with none changed B coverage by +0.002 (SE 0.005). Taking the factors out of the truth as well gave 0.923 (+0.017, SE 0.010), but those are different datasets, and the difference is under two standard errors.
-   **B0** covers at 0.928 in `no_spatial`. Only removing the latent factors from the truth moves it by more than 1.4 points, to 0.950 (+0.022, SE 0.009, different datasets), and fitting factors that are absent costs nothing (-0.001, SE 0.005).
-   **So no single candidate accounts for all of it.** The prior is the main lever, unconverged chains add a smaller independent part, and the two were not run together.

## Run

The six cells ran from 15 October 2026 at 23:04 to 16 October at 00:41 at revision `c40f385`, with a clean tree at start and at write: 600 fits in 97 minutes, no failures, 5 replicate processes, one thread per fit, on mains power with Low Power Mode off, under R 4.6.1 with Accelerate BLAS. `no_spatial_long` ran its own 5,000 + 5,000 iterations per chain, and the runner recorded that in the results file; every other cell ran 1,000 + 1,000, with 2 chains throughout. Every replicate is present in all six arms.

## Coverage by arm

Replicates 1 to 100 in each arm. Standard errors are computed across replicates, as in the calibration report. The error-to-width ratio is near 1 for calibrated intervals and above 1 when they are narrow for the errors the estimates make.

-   `no_spatial`, the control: B 0.906 (0.008), width 3.07, RMSE 0.889, slope 0.674, error-to-width 1.20. B0 0.928 (0.008).
-   `no_spatial_long`, 5,000 + 5,000 iterations: B 0.923 (0.007), width 3.14, RMSE 0.876, slope 0.680, error-to-width 1.15. B0 0.934 (0.008).
-   `no_factors`, no latent factors in the truth or the fit: B 0.923 (0.007), width 3.22, RMSE 0.860, slope 0.756, error-to-width 1.13. B0 0.950 (0.007).
-   `no_factors_fit2`, the same datasets fitted with 2 factors: B 0.924 (0.007), width 3.18, RMSE 0.851, slope 0.751, error-to-width 1.11. B0 0.949 (0.008).
-   `no_spatial_notraits`, traits left out of the fit: B 0.871 (0.009), width 2.92, RMSE 0.983, slope 0.561, error-to-width 1.42. B0 0.914 (0.009).
-   `no_spatial_gt4`, 4 latent traits fitted instead of 2: B 0.936 (0.006), width 3.64, RMSE 0.885, slope 0.769, error-to-width 1.02. B0 0.930 (0.009).

## Paired contrasts

Each contrast is the mean of per-replicate differences in coverage, with its standard error. All but two compare fits of identical datasets; the two that involve factors in the truth compare different datasets drawn from identical B truths.

-   **Longer chains** (`no_spatial_long` minus `no_spatial`, same data): B +0.017 (0.005), intervals 2% wider, mean squared error down 0.023 (0.009). B0 +0.006 (0.005).
-   **4 latent traits instead of 2** (`no_spatial_gt4` minus `no_spatial`, same data): B +0.030 (0.006), intervals 19% wider, mean squared error unchanged (-0.007, SE 0.017). B0 +0.002 (0.006).
-   **Traits left out** (`no_spatial_notraits` minus `no_spatial`, same data): B -0.036 (0.007), intervals 5% narrower, mean squared error up 0.176 (0.026). B0 -0.014 (0.005).
-   **Fitting 2 factors that are not there** (`no_factors_fit2` minus `no_factors`, same data): B +0.002 (0.005), intervals 1% narrower. B0 -0.001 (0.005).
-   **No factors in the truth or the fit** (`no_factors` minus `no_spatial`, different data): B +0.017 (0.010). B0 +0.022 (0.009).
-   **Factors in the truth, 2 fitted** (`no_spatial` minus `no_factors_fit2`, different data): B -0.018 (0.010), slope 0.751 down to 0.674. B0 -0.021 (0.010).

The two prior arms change the estimates, not only the widths. With 4 latent traits the large effects are pulled towards zero by 1.16 on average instead of 1.46, and small effects are pushed away from zero a little more, so the RMSE is unchanged while coverage rises. Without the traits the large effects are pulled in by 1.89. Longer chains barely move the estimates (slope 0.674 to 0.680) and act mainly on the intervals.

## B by size of the true effect

Coverage of B elements grouped by the absolute size of the true effect, 2,000 elements per arm, with the paired difference from `no_spatial` and its standard error, treating replicates as clusters.

-   Under 1 (917 elements): `no_spatial` 0.948. Longer chains 0.958 (+0.010, SE 0.005), 4 latent traits 0.952 (+0.004, SE 0.005), no traits 0.954 (+0.007, SE 0.007).
-   1 to 2 (605): 0.917. Longer chains 0.927 (+0.010, SE 0.008), 4 latent traits 0.934 (+0.017, SE 0.011), no traits 0.886 (-0.031, SE 0.011).
-   2 to 3 (297): 0.889. Longer chains 0.902 (+0.013, SE 0.014), 4 latent traits 0.953 (+0.064, SE 0.018), no traits 0.768 (-0.121, SE 0.026).
-   Over 3 (181): 0.685. Longer chains 0.762 (+0.077, SE 0.024), 4 latent traits 0.829 (+0.144, SE 0.033), no traits 0.564 (-0.122, SE 0.039).

The two arms without factors in the truth cover the largest effects at 0.796 and 0.790. Every candidate acts mostly on the large effects, and with 4 latent traits coverage is close to flat up to effects of 3.

## Chain length and mixing

On identical datasets, going from 2,000 to 10,000 draws moved B's interval endpoints by 9% of the interval width on average, and its posterior means by 4%. The covered flag changed for 4.3% of B elements, 3.0% gaining coverage and 1.3% losing it. That asymmetry alone does not separate noise from bias, since noisy endpoints also lose coverage on balance. What does is that the intervals widen systematically with longer chains (2% on average, 4% at the median element) and that the short chains have not converged, below.

[mixing.R](mixing.R) refitted replicates 1 to 10 of `no_spatial` and `no_spatial_gt4`, and 1 to 5 of `no_spatial_long`, in one process with the draws kept. Each refit reproduced the run's posterior means exactly.

-   `no_spatial`, 2,000 draws: B bulk effective sample size 22, 57 and 117 at the 10th, 50th and 90th percentiles, tail effective sample size 43, 95 and 197. R-hat over 1.01 for 91% of B elements, maximum 1.27. B0 is similar: median bulk 68, R-hat over 1.01 for 83%.
-   `no_spatial_long`, 10,000 draws: B bulk 115, 230 and 464, tail 176, 352 and 758. R-hat over 1.01 for 41%, maximum 1.055.
-   `no_spatial_gt4`, 2,000 draws: B bulk 14, 45 and 103, R-hat over 1.01 for 92%, maximum 1.41. It mixes no better than 2 latent traits, so its gain in coverage is not a mixing effect.

## The prior on B

In the fit, each species' B is normal around `Tr G + A C` with standard deviation `sigma_b`, and the elements of G, A and C have N(0, 1) priors. `sigma_b`'s posterior sits at its prior in every arm: posterior mean 0.33, interval about 0.24 to 0.46, against a truth of 0.35. So B is held close to the trait and latent-trait regression, and how far the estimates can spread across species is set by G and the latent traits. With 2 latent traits the prior's expected variance for a B element is about 4.1, with 4 about 6.1; the truths' mean square is 2.96. The prior is therefore not too narrow on average. What 4 latent traits change is the shape of the prior on how spread out B is across species, which puts less weight on small spreads. That is a reading of the result, not something this run tests.

In `no_spatial`'s truth the latent factors carry 42% of the linear predictor's variance, averaged over species. With them in the data the estimates of B are more shrunk (slope 0.674 against 0.751 on the factor-free data, both fitted with 2 factors). That fits the same picture: anything that makes the data less informative about B lets the prior's shrinkage bite harder, and fitting factors as such costs nothing.

## Pairing

The truths of B, B0, G, `sigma_b`, `theta0`, `p` and `q` are identical, element for element, in all six arms ([pairing-check.csv](pairing-check.csv)). The collection slopes in `beta_theta` differ in the two no-factor arms, because the simulator draws them after the factors; `beta_theta` is not compared here. Re-simulating all 100 replicates confirmed that `no_spatial_long`, `no_spatial_notraits` and `no_spatial_gt4` fit exactly `no_spatial`'s datasets, that `no_factors_fit2` fits exactly `no_factors`', and that the no-factor datasets differ ([dataset-pairing.csv](dataset-pairing.csv)). This run's `no_spatial` reproduces the 2x2's bit for bit in all 14,262 rows. The harness tests in `tests/testthat/test-simstudy-b-candidates.R` check the pairing at replicate 1.

## What this does not show

-   **4 latent traits is a diagnostic, not a fix.** It fits more latent traits than the simulator uses, G then over-covers (0.973 against 0.945), and it was tested in one design. It shows that the prior's shrinkage is the lever, not that this setting should be the default.
-   **The `sigma_b` layer was not tested.** `runOccJSDM()` hard-codes its InvGamma(10, 1) prior and `listPriors` cannot change it, so a looser prior there needs a package change.
-   **The two main effects were not combined.** Whether 4 latent traits at 5,000 + 5,000 iterations reaches 0.95 is not known. With 4 latent traits B still covers at 0.936 (SE 0.006).
-   **Chain length is only partly measured.** Even at 10,000 draws 41% of B elements have R-hat over 1.01, so the +0.017 is a lower bound on what converged chains would give. The mixing figures come from 10 replicates, and 5 for the long chains.
-   **These truths are not draws from the fit's prior.** The simulator draws G and C from -1, 0 and 1 with C's diagonal fixed at 1, which the fit's normal priors do not reproduce. Coverage here is measured against the simulator's truths, so the prior result is a mismatch between prior and truth at least as much as a defect. Simulation-based calibration (PLAN.md 7), which draws the truth from the prior, would separate the two.
-   **The factor-in-the-truth contrasts are suggestive only.** They compare different datasets and sit at 1.7 to 2.1 standard errors: B -0.018 (0.010), B0 -0.021 (0.010).
-   **One design.** Only the base cell's settings without a spatial field were run: 10 species, 2 covariates, 2 traits. At 100 replicates a paired contrast of about 1.5 points on the same data is detectable and one of 1 point is not.

## What this changes

-   **The paper's B caveat** can be attributed: mainly to shrinkage under B's hierarchical prior with 2 latent traits, which over-pulls large effects for the study's truths, and secondly to chains too short to converge. Neither the spatial field nor latent factors competing with the covariates is a cause.
-   **The study runs chains a fifth as long as the package default.** `runOccJSDM()` defaults to 5,000 + 5,000 iterations, and the study uses 1,000 + 1,000. In this cell that costs B 1.7 points of coverage, so the study understates the calibration a user gets at the defaults, and B is far from converged at the study's length.
-   **Next tests, in order of cost:** 4 latent traits at 5,000 + 5,000 iterations, to see whether the two together reach nominal; the same pair in one spatial cell, to check that the result carries over; a `listPriors` option for `sigma_b`'s prior, which is a package change and Alex's call; and simulation-based calibration for B. B's slow mixing, with median ESS 57 of 2,000 draws, also makes the sampler's separate updates of B and of G, A and C worth a look; that is a hypothesis.

## Files and reproduction

-   [analyse.R](analyse.R) writes [coverage-by-arm.csv](coverage-by-arm.csv), [paired-contrasts.csv](paired-contrasts.csv), [b-by-effect-size.csv](b-by-effect-size.csv), [b-by-effect-size-contrasts.csv](b-by-effect-size-contrasts.csv), [chain-length.csv](chain-length.csv) and [pairing-check.csv](pairing-check.csv) from this run's results file and the 2x2's.
-   [mixing.R](mixing.R) writes [mixing.csv](mixing.csv) and [dataset-pairing.csv](dataset-pairing.csv). It refits 25 replicates and takes about 15 minutes.
-   `results/simstudy-20261016-004116.rds` holds all six arms; `run.log` and `provenance.txt` record the run, and [launch.sh](launch.sh) launched it.

Run from the repository root:

```sh
Rscript dev/simstudy/b-nonspatial-candidates-20261015/analyse.R
Rscript dev/simstudy/b-nonspatial-candidates-20261015/mixing.R
```

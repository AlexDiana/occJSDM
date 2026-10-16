# 4 latent traits with long chains

The non-spatial candidates run ([its report](../b-nonspatial-candidates-20261015/REPORT.md)) found two levers on B's 95% coverage in `no_spatial`, each tested alone: 4 latent traits instead of 2 (0.906 to 0.936) and 5,000 + 5,000 iterations per chain instead of 1,000 + 1,000 (0.906 to 0.923). This run fits `no_spatial`'s 100 datasets with both, in the cell `no_spatial_gt4_long`, and pairs it against `no_spatial`, `no_spatial_gt4` and `no_spatial_long`. The design is in [README.md](README.md).

## Summary

-   **B reaches nominal on average.** Coverage is 0.956 (SE 0.005, 95% interval 0.946 to 0.966), against 0.906 in `no_spatial` on the same datasets, a paired gain of +0.050 (SE 0.006). The error-to-width ratio is 0.96, against 1.20, and the slope of the posterior means on the truth is 0.769, against 0.674.
-   **The two effects add.** Alone, 4 latent traits gave +0.030 (SE 0.006) and long chains +0.017 (SE 0.005); the additive prediction is 0.952 and the observed coverage is 0.956. The interaction, (`gt4_long` minus `gt4`) minus (`long` minus `no_spatial`) per replicate, is +0.004 (SE 0.007), a 95% interval of -0.009 to +0.017. Long chains add +0.021 (SE 0.004) given 4 latent traits, and 4 latent traits add +0.034 (SE 0.005) given long chains.
-   **Not uniformly across effect sizes.** True effects over 3 in absolute value now cover at 0.895 (SE 0.025), up from 0.685, while effects under 3 cover at 0.961 to 0.966. The average of 0.956 mixes slight over-coverage of small and medium effects with remaining under-coverage of the largest.
-   **B0** covers at 0.941 (SE 0.008), +0.013 (SE 0.005) over `no_spatial`; neither lever's share is clearly separated from zero. **G** over-covers at 0.983 (SE 0.006), against 0.945.
-   **B is still not converged.** At 10,000 draws B's median bulk effective sample size is 207 and 41% of its elements have R-hat over 1.01, no better than with 2 latent traits at the same length.

## Run

The cell ran on 16 October 2026 from 07:18 to 07:47 at revision `815ba6f`: 100 fits in 29.4 minutes, about 88 seconds each, no failures, 5 replicate processes, one thread per fit, 2 chains of 5,000 + 5,000, on mains power with Low Power Mode off, under R 4.6.1 with Accelerate BLAS. A stuck system process (`appstoreagent`) was holding one core at about 97% at launch. The tree was clean at launch. The runner records it as dirty at write because `analyse.R` and `mixing.R` were written, untracked, in this directory while the run was in flight; nothing the run reads changed. The package code (`R/`, `src/`, `DESCRIPTION`, `NAMESPACE`) and `run_study.R` are identical to `c40f385`, the revision that produced the three comparison arms.

## Pairing

Checked before anything else. The truths of every block (B, B0, G, `sigma_b`, `beta_theta`, `theta0`, `p`, `q` and the residual correlations) are identical element for element between the new cell and each of the three comparison arms, with the same rows present in both: 42,786 rows compared, maximum difference 0 ([pairing-check.csv](pairing-check.csv)). Re-simulating all 100 replicates confirmed the cell hands the fitter exactly `no_spatial`'s datasets ([dataset-pairing.csv](dataset-pairing.csv)), and the five refits in [mixing.R](mixing.R) reproduced the run's posterior means exactly.

## Coverage by arm

100 replicates per arm, standard errors across replicates. The error-to-width ratio is near 1 for calibrated intervals and above 1 when the intervals are narrow for the errors.

-   `no_spatial`, 2 latent traits, 1,000 + 1,000: B 0.906 (0.008), width 3.07, RMSE 0.889, slope 0.674, error-to-width 1.20. B0 0.928 (0.008), width 3.01, error-to-width 1.08.
-   `no_spatial_long`, 2 latent traits, 5,000 + 5,000: B 0.923 (0.007), width 3.14, RMSE 0.876, slope 0.680, error-to-width 1.15. B0 0.934 (0.008), width 3.05, error-to-width 1.05.
-   `no_spatial_gt4`, 4 latent traits, 1,000 + 1,000: B 0.936 (0.006), width 3.64, RMSE 0.885, slope 0.769, error-to-width 1.02. B0 0.930 (0.009), width 3.07, error-to-width 1.05.
-   `no_spatial_gt4_long`, 4 latent traits, 5,000 + 5,000: B 0.956 (0.005), width 3.71, RMSE 0.864, slope 0.769, error-to-width 0.96. B0 0.941 (0.008), width 3.12, error-to-width 1.02.

## Paired contrasts

Means of per-replicate differences, all on identical datasets ([paired-contrasts.csv](paired-contrasts.csv)).

-   **Both against neither** (`gt4_long` minus `no_spatial`): B +0.050 (0.006), intervals 21% wider, mean squared error down 0.044 (0.017). B0 +0.013 (0.005), intervals 4% wider.
-   **Long chains, given 4 latent traits** (`gt4_long` minus `gt4`): B +0.021 (0.004), intervals 2% wider, mean squared error down 0.037 (0.013). B0 +0.011 (0.005).
-   **4 latent traits, given long chains** (`gt4_long` minus `long`): B +0.034 (0.005), intervals 18% wider, mean squared error down 0.021 (0.013). B0 +0.007 (0.004).
-   **Interaction**: B +0.004 (0.007), B0 +0.005 (0.008). The width ratios multiply almost exactly (ratio of ratios 0.997).

The two levers act on different things. 4 latent traits change the estimates: large effects are pulled towards zero by 1.18 on average instead of 1.46, and the slope rises from 0.674 to 0.769 whatever the chain length. Long chains leave the estimates' shrinkage alone (slope 0.769 with and without them at 4 latent traits) and act on the intervals and on Monte Carlo error.

## B by size of the true effect

Coverage by the absolute size of the true effect, with cluster standard errors over replicates, and the paired gain of `gt4_long` over `no_spatial` ([b-by-effect-size.csv](b-by-effect-size.csv), [b-by-effect-size-contrasts.csv](b-by-effect-size-contrasts.csv)).

-   Under 1 (917 elements): `no_spatial` 0.948, `long` 0.957, `gt4` 0.952, `gt4_long` 0.961 (0.006). Gain +0.013 (0.006).
-   1 to 2 (605): 0.917, 0.927, 0.934, 0.962 (0.008). Gain +0.045 (0.010).
-   2 to 3 (297): 0.889, 0.902, 0.953, 0.966 (0.010). Gain +0.077 (0.017).
-   Over 3 (181): 0.685, 0.762, 0.829, 0.895 (0.025). Gain +0.210 (0.035).

The interaction is within 1.5 standard errors of zero in every bin. Over 3, long chains add +0.066 (0.020) at 4 latent traits, about what they added at 2.

## Chain length and mixing

At 4 latent traits, going from 2,000 to 10,000 draws moved B's interval endpoints by 10% of the interval width on average and its posterior means by 4%, and changed the covered flag for 3.4% of B elements, 2.7% gaining coverage and 0.7% losing it ([chain-length.csv](chain-length.csv)). At 2 latent traits the figures were 9%, 4% and 4.3%.

[mixing.R](mixing.R) refitted replicates 1 to 5 of the new cell with the draws kept ([mixing.csv](mixing.csv)). At 10,000 draws B's bulk effective sample size is 85, 207 and 401 at the 10th, 50th and 90th percentiles, tail 122, 287 and 657, with R-hat over 1.01 for 41% of B elements (maximum 1.056). `no_spatial_long` at the same replicates gave 115, 230 and 464, with 41% over 1.01. B0's median bulk is 295, with 22% over 1.01. At this length 4 latent traits mix about as well as 2, slightly worse at the low end (10th percentile 85 against 115).

## Limits

-   **Nominal on average, not conditionally.** The largest true effects still undercover (0.895, SE 0.025) and effects under 3 slightly over-cover. A marginal 0.95 here is partly a balance of the two.
-   **4 latent traits remains a diagnostic.** The simulator uses 2, G over-covers at 0.983, and the setting was tested in one design. The result locates the problem in B's hierarchical prior; it does not make 4 latent traits a recommended default.
-   **Not converged.** 41% of B elements have R-hat over 1.01 even at 10,000 draws, so the long-chain part of the gain may still be a lower bound. The mixing figures come from 5 replicates.
-   **One design.** Only the base cell's settings without a spatial field: 10 species, 2 covariates, 2 traits. No spatial cell has been run with both levers.
-   **Truths not drawn from the fit's prior.** The simulator draws G and C from -1, 0 and 1, so coverage is measured against truths the fit's normal priors do not generate. Simulation-based calibration (PLAN.md 7) would separate prior mismatch from sampler error.
-   **The `sigma_b` layer is untested.** Its InvGamma(10, 1) prior is hard-coded in `runOccJSDM()`.

## What this changes

-   **B's shortfall in `no_spatial` is accounted for.** The two levers together close the 4.4-point gap to 0.95, additively: about 3 points from the shrinkage of B's hierarchical prior and about 2 from chains too short to converge, depending a little on which is added first.
-   **The study's chain length matters for the paper.** At the package's default 5,000 + 5,000 iterations, B gains 1.7 to 2.1 points of coverage in this cell whichever prior is fitted, so the study's 1,000 + 1,000 understates the calibration a user gets.
-   **Next tests:** both levers in one spatial cell, to check the result carries over; a `listPriors` option for `sigma_b`'s prior, a package change and Alex's call; simulation-based calibration for B.

## Files and reproduction

-   [analyse.R](analyse.R) writes [pairing-check.csv](pairing-check.csv), [coverage-by-arm.csv](coverage-by-arm.csv), [paired-contrasts.csv](paired-contrasts.csv), [b-by-effect-size.csv](b-by-effect-size.csv), [b-by-effect-size-contrasts.csv](b-by-effect-size-contrasts.csv) and [chain-length.csv](chain-length.csv) from this run's results and the candidates run's.
-   [mixing.R](mixing.R) writes [mixing.csv](mixing.csv) and [dataset-pairing.csv](dataset-pairing.csv). It refits 5 replicates and takes about 5 minutes.
-   `results/simstudy-20261016-074749.rds` holds the cell; `run.log` and `provenance.txt` record the run, and [launch.sh](launch.sh) launched it.

Run from the repository root:

```sh
Rscript dev/simstudy/b-gt4-long-20261016/analyse.R
Rscript dev/simstudy/b-gt4-long-20261016/mixing.R
```

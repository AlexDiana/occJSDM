# The prior on sigma_b

The two B reports ([non-spatial candidates](../b-nonspatial-candidates-20261015/REPORT.md), [4 latent traits with long chains](../b-gt4-long-20261016/REPORT.md)) traced B's undercoverage in `no_spatial` to shrinkage under B's hierarchical prior and to unconverged chains, and left one lever untested: the hard-coded InvGamma(10, 1) prior on `sigma_b^2`, under which `sigma_b`'s posterior sat at the prior's 0.33 in every arm. This run adds a `listPriors` option for that prior and fits `no_spatial`'s 100 datasets under a half-Cauchy(1) on `sigma_b` and under a weak InvGamma(2, 1/9) on `sigma_b^2`, each at 1,000 + 1,000 and at 5,000 + 5,000 iterations per chain, plus the half-Cauchy with 4 latent traits at 5,000 + 5,000, and pairs each arm with the earlier arm it changes one thing from. The design, the priors and the success criteria were written before the run, in [README.md](README.md).

## Summary

-   **A freer `sigma_b` partly restores B's coverage and does not reach nominal.** The half-Cauchy raises B's 95% coverage by +0.015 (SE 0.006) at 1,000 + 1,000 and +0.016 (SE 0.004) at 5,000 + 5,000, on identical datasets. At the package's default chain length B covers at 0.938 (SE 0.006, 95% interval 0.927 to 0.949) against 0.922 under the default prior. By the pre-registered criteria that is "partly restores": the gain is clear, G stays at 0.963 and B0 is unchanged, but the interval excludes 0.95.
-   **It is the form of the prior that matters, not only its weight.** The weak inverse-gamma with the default's prior mean changes nothing: +0.002 (SE 0.006) short, +0.002 (SE 0.004) long. Under either inverse-gamma prior `sigma_b`'s posterior is its prior; under the half-Cauchy it spreads over 0.03 to 1.6 and its mean lands at 0.57, above the truth of 0.35 in 88% of replicates.
-   **The gain is in the large effects and is smaller than the latent-trait probe's.** True effects over 3 cover at 0.829 (SE 0.032) under the half-Cauchy at 5,000 + 5,000, up from 0.762 but short of the 0.90 threshold and of the 0.895 that 4 latent traits reached. Against `no_spatial_gt4_long` the half-Cauchy falls 0.018 (SE 0.004) short, and large effects are still pulled towards zero by 1.35 on average against 1.18.
-   **It does not replace the probe, and it adds nothing to it.** Given 4 latent traits and long chains, the half-Cauchy changes B coverage by -0.003 (SE 0.003); given the half-Cauchy, 4 latent traits still add +0.016 (SE 0.005). The two act on the same thing and the latent-trait count does more of it.
-   **Chain length is still needed, and it adds.** Going to 5,000 + 5,000 adds +0.017 (SE 0.004) under the half-Cauchy, the same as under the default prior, with an interaction of +0.001 (SE 0.006).
-   **Cost.** B's intervals are 7% wider, their mean squared error slightly lower (-0.016, SE 0.006, at the long length). B0 is unaffected (+0.006, SE 0.004). G's coverage rises by +0.008 (SE 0.006) to 0.963 with 6% wider intervals, within the 0.97 limit. `sigma_b` itself becomes essentially unidentified: its interval is 7 times wider and its posterior mean ranges from 0.24 to 1.5 across replicates. And the chains mix worse: at 10,000 draws B's median bulk effective sample size falls from 230 to 169 and the share of B elements with R-hat over 1.01 rises from 41% to 72%, while `sigma_b` itself has a bulk effective sample size of 34 to 89.
-   **Mixing is a second cost.** B's chains mix worse under the half-Cauchy than under the default at either length, and `sigma_b` itself moves slowly once it is free: bulk effective sample sizes of 5 to 31 in 2,000 draws and 34 to 89 in 10,000.

## Run

The five cells ran on 17 October 2026 from 16:14 to 18:11 at revision `6d1d0c2` on branch `claude/sigma-b-prior`: 500 fits in 117 minutes, no failures, 5 replicate processes, one thread per fit, 2 chains, on mains power with Low Power Mode off, under R 4.6.1 with Accelerate BLAS. The three long cells ran their own 5,000 + 5,000 iterations and the runner recorded that in the results file; the two short cells ran the runner's 1,000 + 1,000. The tree was clean at launch; the runner records it as dirty at write because `analyse.R` and `mixing.R` were written, untracked, in this directory while the run was in flight. Nothing the run reads changed. The package code is `c40f385`, which produced the comparison arms, plus the option, whose default path is bit-identical ([test-sigma-b-prior.R](../../../tests/testthat/test-sigma-b-prior.R)).

## Pairing

Checked before anything else. The truths of every block (B, B0, G, `sigma_b`, `beta_theta`, `theta0`, `p`, `q`) are identical element for element between each new arm and each of the four comparison arms, with the same rows present: 210,000 rows compared, maximum difference 0, and no posterior mean shared between a new arm and an old one ([pairing-check.csv](pairing-check.csv)). Re-simulating all 100 replicates confirmed that every new cell hands the fitter exactly `no_spatial`'s datasets ([dataset-pairing.csv](dataset-pairing.csv)).

## Coverage by arm

100 replicates per arm, standard errors across replicates. The error-to-width ratio is near 1 for calibrated intervals and above 1 when the intervals are narrow for the errors.

-   `no_spatial`, default prior, 1,000 + 1,000: B 0.906 (0.008), width 3.07, RMSE 0.889, slope 0.674, error-to-width 1.20. B0 0.928 (0.008). G 0.945 (0.013), width 1.84.
-   `no_spatial_sbhc`, half-Cauchy(1), 1,000 + 1,000: B 0.921 (0.007), width 3.27, RMSE 0.888, slope 0.702, error-to-width 1.14. B0 0.934 (0.008). G 0.945 (0.013), width 1.94.
-   `no_spatial_sbig2`, InvGamma(2, 1/9), 1,000 + 1,000: B 0.908 (0.008), width 3.05, RMSE 0.893, slope 0.669, error-to-width 1.23. B0 0.926 (0.008). G 0.940 (0.012), width 1.83.
-   `no_spatial_long`, default prior, 5,000 + 5,000: B 0.922 (0.007), width 3.14, RMSE 0.876, slope 0.680, error-to-width 1.15. B0 0.934 (0.008). G 0.955 (0.013), width 1.89.
-   `no_spatial_sbhc_long`, half-Cauchy(1), 5,000 + 5,000: B 0.938 (0.006), width 3.36, RMSE 0.867, slope 0.712, error-to-width 1.07. B0 0.940 (0.008). G 0.963 (0.011), width 1.99.
-   `no_spatial_sbig2_long`, InvGamma(2, 1/9), 5,000 + 5,000: B 0.924 (0.007), width 3.12, RMSE 0.881, slope 0.677, error-to-width 1.16. B0 0.936 (0.007). G 0.953 (0.011), width 1.88.
-   `no_spatial_gt4_long`, default prior, 4 latent traits, 5,000 + 5,000: B 0.956 (0.005), width 3.71, RMSE 0.864, slope 0.769, error-to-width 0.96. B0 0.941 (0.008). G 0.983 (0.006), width 2.18.
-   `no_spatial_sbhc_gt4_long`, half-Cauchy(1), 4 latent traits, 5,000 + 5,000: B 0.954 (0.006), width 3.93, RMSE 0.872, slope 0.792, error-to-width 0.93. B0 0.944 (0.007). G 0.985 (0.006), width 2.27.

## Where sigma_b goes

The simulated truth is 0.35 in every replicate ([sigma-b-posterior.csv](sigma-b-posterior.csv)). Posterior means are averaged across replicates, with their spread across replicates in parentheses, and the interval is the mean 95% posterior interval.

-   **Default InvGamma(10, 1)**: 0.329 (spread 0.004), interval 0.24 to 0.46, width 0.21. That interval is the prior's own 95% interval (0.24 to 0.46) and the spread across replicates is a hundredth of its width: the posterior is the prior, whatever the data.
-   **InvGamma(2, 1/9)**: 0.293 (spread 0.02), interval 0.14 to 0.64, width 0.50, below the truth in every replicate. Again the prior's own interval (0.14 to 0.68). The lighter prior moved `sigma_b` down, to its own median, not towards the truth; shape 2 against 10 B-element pairs still leaves the data without a say.
-   **Half-Cauchy(1)**: 0.573 (spread 0.23, from 0.24 to 1.5), interval 0.05 to 1.58, width 1.53, 7 times the default's. The posterior mean is above the truth in 88% of replicates and the interval always covers it. The spread across replicates and the width within each say that 20 B elements do not identify `sigma_b` once the prior lets go, and the mean is pulled up by the posterior's heavy right tail (in the five refitted replicates the posterior medians are 0.27 to 0.70, about 0.1 below the means, and above the truth in four of the five; see the mixing section).
-   **Half-Cauchy(1) with 4 latent traits**: 0.516 (spread 0.13), interval 0.03 to 1.58. With more of B's spread carried by the latent traits, `sigma_b` settles a little lower.

Under the half-Cauchy the chain length does not move `sigma_b` (0.573 at both lengths), and its interval widens by 5% with longer chains.

## Paired contrasts

Means of per-replicate differences, all on identical datasets ([paired-contrasts.csv](paired-contrasts.csv)). Standard errors in parentheses.

-   **Half-Cauchy, 1,000 + 1,000** (`sbhc` minus `no_spatial`): B +0.015 (0.006), intervals 7% wider, mean squared error unchanged (-0.002, SE 0.014). B0 +0.006 (0.005). G 0.000 (0.010), intervals 6% wider.
-   **Weak InvGamma, 1,000 + 1,000** (`sbig2` minus `no_spatial`): B +0.002 (0.006), intervals 1% narrower. B0 -0.002 (0.006). G -0.005 (0.009).
-   **Half-Cauchy against weak InvGamma, 1,000 + 1,000**: B +0.013 (0.005).
-   **Half-Cauchy, 5,000 + 5,000, the primary contrast** (`sbhc_long` minus `no_spatial_long`): B +0.016 (0.004), intervals 7% wider, mean squared error down 0.016 (0.006), slope 0.680 to 0.712. B0 +0.006 (0.004), intervals 1% wider. G +0.008 (0.006), intervals 6% wider.
-   **Weak InvGamma, 5,000 + 5,000** (`sbig2_long` minus `no_spatial_long`): B +0.002 (0.004). B0 +0.002 (0.004). G -0.003 (0.008).
-   **Half-Cauchy against weak InvGamma, 5,000 + 5,000**: B +0.014 (0.004).
-   **Long chains, given the half-Cauchy** (`sbhc_long` minus `sbhc`): B +0.017 (0.004), intervals 3% wider, mean squared error down 0.037 (0.011). B0 +0.006 (0.005). G +0.018 (0.006). Given the weak InvGamma, +0.016 (0.004); given the default prior, +0.017 (0.005).
-   **Interaction of prior and chain length**, (`sbhc_long` minus `sbhc`) minus (`long` minus `no_spatial`): B +0.001 (0.006), a 95% interval of -0.012 to +0.013. The width ratios multiply exactly (ratio of ratios 1.003). The two effects add.
-   **Half-Cauchy against the latent-trait probe** (`sbhc_long` minus `no_spatial_gt4_long`): B -0.018 (0.004), intervals 10% narrower. B0 -0.001 (0.005). G -0.020 (0.008).
-   **4 latent traits, given the half-Cauchy and long chains** (`sbhc_gt4_long` minus `sbhc_long`): B +0.016 (0.005), intervals 17% wider. B0 +0.004 (0.005). G +0.023 (0.008), intervals 14% wider.
-   **Half-Cauchy, given 4 latent traits and long chains** (`sbhc_gt4_long` minus `no_spatial_gt4_long`): B -0.003 (0.003), intervals 6% wider, mean squared error up 0.014 (0.008). B0 +0.003 (0.003). G +0.003 (0.003).
-   **Interaction of prior and latent traits at 5,000 + 5,000**, (`sbhc_gt4_long` minus `sbhc_long`) minus (`gt4_long` minus `long`): B -0.018 (0.005). The two do not add: the half-Cauchy's gain is a part of what 4 latent traits already give.

## B by size of the true effect

Coverage by the absolute size of the true effect, with cluster standard errors over replicates ([b-by-effect-size.csv](b-by-effect-size.csv), [b-by-effect-size-contrasts.csv](b-by-effect-size-contrasts.csv)). The pull is the mean of the true absolute effect minus the estimated absolute effect, so a positive pull is shrinkage towards zero.

-   Under 1 (917 elements): `no_spatial` 0.948, `long` 0.957, `sbhc` 0.956, `sbhc_long` 0.955, `sbig2_long` 0.957, `gt4_long` 0.961, `sbhc_gt4_long` 0.955. The half-Cauchy changes nothing here (`sbhc_long` minus `long` -0.002, SE 0.004).
-   1 to 2 (605): 0.917, 0.927, 0.931, 0.947 (0.011), 0.927, 0.962, 0.950. Gain of the half-Cauchy at the long length +0.020 (0.006).
-   2 to 3 (297): 0.889, 0.902, 0.899, 0.933 (0.015), 0.912, 0.966, 0.980. Gain +0.030 (0.012).
-   Over 3 (181): 0.685, 0.762, 0.746, 0.829 (0.032), 0.762, 0.895, 0.912 (0.025). Gain +0.066 (0.021); against `gt4_long`, -0.066 (0.023). Pull 1.46, 1.45, 1.39, 1.35, 1.46, 1.18, 1.11.

The half-Cauchy acts only on effects above 1, most on the largest, and there it does about half of what 4 latent traits do. With both, effects over 3 reach 0.912 and effects from 2 to 3 over-cover at 0.980.

## Chain length and mixing

Under the half-Cauchy, going from 2,000 to 10,000 draws moved B's interval endpoints by 10% of the interval width on average and its posterior means by 4%, and changed the covered flag for 3.9% of B elements, 2.8% gaining coverage and 1.1% losing it ([chain-length.csv](chain-length.csv)). Under the default prior the figures were 9%, 4% and 4.3%, so the prior does not change what longer chains do. `sigma_b`'s own interval under the half-Cauchy moves its upper end by 14% of its width with longer chains and its lower end by 2%.

[mixing.R](mixing.R) refitted replicates 1 to 5 of each new arm, and of `no_spatial` and `no_spatial_long` for `sigma_b`'s figures under the default prior, with the draws kept ([mixing.csv](mixing.csv)); every refit reproduced the run's posterior means exactly. Effective sample sizes are bulk, at the 10th, 50th and 90th percentiles over B's elements, and R-hat is the share of elements over 1.01.

-   **B mixes worse under the half-Cauchy.** At 10,000 draws `no_spatial_sbhc_long` gives B a bulk ESS of 69, 169 and 292, with 72% of elements over R-hat 1.01 (maximum 1.09), against 115, 230 and 464 and 41% in `no_spatial_long` at the same replicates. At 2,000 draws it is 11, 50 and 105 with 94% over 1.01, against 26, 56 and 110 and 92% under the default. The weak inverse-gamma is between the two (197 median, 56% over 1.01 at 10,000 draws), and the half-Cauchy with 4 latent traits is the worst (128 median, 70%). B0 is similar across priors (median bulk ESS 263 to 297 at 10,000 draws).
-   **`sigma_b` itself mixes slowly once it is free.** Under the default prior its bulk ESS is 570 to 730 of 2,000 draws and 2,760 to 3,280 of 10,000: it is the prior, and the prior mixes. Under the half-Cauchy it is 5 to 31 of 2,000 draws (R-hat up to 1.31) and 34 to 89 of 10,000 (R-hat up to 1.04); under the weak inverse-gamma 79 to 245 and 550 to 850. The auxiliary-variable Gibbs step for the half-Cauchy moves `sigma_b` in small steps when B's residuals around the regression are themselves uncertain, and B and `sigma_b` then move together slowly.
-   **Posterior median against mean.** In the five refitted replicates of `no_spatial_sbhc_long` the posterior medians of `sigma_b` are 0.38, 0.27, 0.61, 0.70 and 0.44, about 0.1 below the posterior means (0.47, 0.34, 0.65, 0.80, 0.52), with 2.5% quantiles near zero and 97.5% quantiles of 1.1 to 2.0. The heavy right tail lifts the mean; the median is above the truth of 0.35 in four of the five, so the shift upward is not only the tail.

## The prior on B, read from these arms

In the fit each species' B is normal around `Tr G + A C` with standard deviation `sigma_b`, with N(0, 1) priors on the elements of G, A and C. The results separate two layers of that prior.

-   **The `sigma_b` layer is not where most of the shortfall is.** Letting `sigma_b` go free buys 1.5 points of B coverage and leaves large effects pulled in by 1.35 instead of 1.45. Even with `sigma_b`'s posterior spread over 0.03 to 1.6, B stays close to its trait regression.
-   **The regression's own shrinkage is the larger part.** 4 latent traits, which widen what `Tr G + A C` can reach, buy 3.4 points and bring the pull on large effects to 1.18, and the half-Cauchy adds nothing on top of them. The N(0, 1) priors on A and C, and on G, set how far the regression mean can move, and `sigma_b` cannot compensate for a mean that is too shrunk because the data do not identify it (20 B elements per fit).
-   **Within the inverse-gamma family the data have no say.** Both inverse-gamma arms return the prior as the posterior, with the weak one pulling `sigma_b` down to its own median. A looser inverse-gamma is therefore not a route; the half-Cauchy's heavy tail is what lets `sigma_b` leave 0.33.

That reading is consistent with every contrast here, but the prior on A, C and G was not varied in this run and it is a hypothesis.

## Limits

-   **Partly, not restored.** B reaches 0.938 (0.927 to 0.949) at the package's default chain length, not 0.95, and the largest effects still cover at 0.829.
-   **`sigma_b` is unidentified under the half-Cauchy.** Its posterior mean varies from 0.24 to 1.5 across replicates of the same design, with intervals from near 0 to 1.6. A user reading `sigma_b` from a half-Cauchy fit of a small community learns little about it; the option changes B's intervals, not the estimate of `sigma_b`.
-   **One scale, one design.** Only half-Cauchy(1) and InvGamma(2, 1/9) were run, in the base cell's settings without a spatial field: 10 species, 2 covariates, 2 traits, 2 latent traits. A different scale, a spatial cell or a larger community were not tested.
-   **Truths not drawn from the fit's prior.** As in the earlier reports, the simulator draws G and C from -1, 0 and 1, so coverage is measured against truths the fit's normal priors do not generate. Simulation-based calibration (PLAN.md 7) would separate prior mismatch from sampler error.
-   **Not converged, and the half-Cauchy makes it worse.** 72% of B elements have R-hat over 1.01 at 10,000 draws under the half-Cauchy, against 41% under the default, and `sigma_b`'s own draws have a bulk effective sample size of 34 to 89 in 10,000, so the long-chain gain may still be a lower bound, and a user of this option needs longer chains than the default to reach the same precision.
-   **The mixing figures come from 5 replicates per arm.**

## What this changes

-   **The `sigma_b` prior is a real but second-order lever.** The half-Cauchy is a cheap partial improvement with no cost to B0 or G: +1.5 points on B at either chain length, B at 0.938 at the default length. Whether to offer it as an option, make it the default, or leave the default alone is Alex's call; the default is unchanged in this branch and nothing already run is affected by the option.
-   **The attribution of B's shortfall sharpens.** It is mostly shrinkage of B's trait regression, the G, A and C layer, not the `sigma_b` layer, plus chains too short to converge. The two B reports' "shrinkage under B's hierarchical prior" is now located one layer up.
-   **Next tests:** a wider or heavy-tailed prior on the latent-trait loadings A and C (and on G), which is the layer these results point at; the half-Cauchy at a larger scale, to check 1 is not the limit; the half-Cauchy in one spatial cell; simulation-based calibration for B.

## Files and reproduction

-   [analyse.R](analyse.R) writes [pairing-check.csv](pairing-check.csv), [coverage-by-arm.csv](coverage-by-arm.csv), [sigma-b-posterior.csv](sigma-b-posterior.csv), [paired-contrasts.csv](paired-contrasts.csv), [b-by-effect-size.csv](b-by-effect-size.csv), [b-by-effect-size-contrasts.csv](b-by-effect-size-contrasts.csv), [chain-length.csv](chain-length.csv) and [criteria.csv](criteria.csv), the pre-registered criteria evaluated mechanically, from this run's results and the two earlier runs'.
-   [mixing.R](mixing.R) writes [mixing.csv](mixing.csv) and [dataset-pairing.csv](dataset-pairing.csv). It refits 5 replicates of each of the five new arms and of `no_spatial` and `no_spatial_long` on 5 workers and takes about 15 minutes.
-   `results/simstudy-20261017-181148.rds` holds the five arms; `run.log` and `provenance.txt` record the run, and [launch.sh](launch.sh) launched it.

Run from the repository root:

```sh
Rscript dev/simstudy/sigma-b-prior-20261017/analyse.R
Rscript dev/simstudy/sigma-b-prior-20261017/mixing.R
```

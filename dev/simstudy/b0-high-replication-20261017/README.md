# Why B0 is miscalibrated at high field replication, 17 October 2026

The community-mean ladder ([REPORT](../community-mean-ladder-20261016/REPORT.md)) found that at 20 field samples per site occupancy is estimated well while B0 is not. In the rare-dominated community with B0's community mean (`base_rare_cm_m20`), the community occupancy error is +0.003 and the mean absolute error 0.018, yet B0 is overestimated by 0.38 on average and by 1.22 for true B0 below -4.25, B0 coverage is 0.82, and `mu_b0` sits 0.43 above the realised mean of the true B0. In the ordinary community at 20 field samples per site (`M20`, `M20_cm`), under either prior, true B0 below -1.5 is overestimated by about 0.5 and above 1.5 underestimated by about 0.5, with coverage 0.68 to 0.74, and overall B0 coverage is 0.87.

This run asks why, and tests the leading explanation.

## Candidate explanations

Each species' occupancy linear predictor is `eta = B0 + r`, where `r` is the sum of the covariate term (X times B), the latent-factor term (U times L) and the spatial term. Write `r` for that spread around B0.

-   **Spread (the leading hypothesis).** The fit underestimates the spread of `r` across sites, through the shrinkage of B (sigma_b's hard-coded InvGamma(10, 1) prior), of the latent factors, or of the spatial term. A species' share of occupied sites is the site average of plogis(B0 + r), and plogis is convex below 0 and concave above it. So with less spread a rare species needs a higher B0, and a common species a lower one, to keep the share the data fix (Jensen's inequality). That is the direction of both errors.
-   **B0's prior.** B0 has a Normal(centre, 1) prior, with the centre at 0 (default) or at the estimated `mu_b0` (community mean). Its pull towards the centre is largest where the data say least, which is the rarest species, and it gives errors of the observed sign at both extremes. In a correctly specified model this pull is the right Bayesian answer and leaves the marginal coverage nominal, so it can explain the pattern by true B0 but not the overall coverage of 0.82 to 0.87 by itself.
-   **Location.** The site averages of the fitted latent-factor and spatial terms need not match those of the truth's factors and field, and B0 absorbs the difference. That adds error the intervals may not carry, but it has no reason to depend on true B0.

The three are not exclusive. An accounting identity (below) splits each species' B0 error into terms, but it cannot say which way the causation runs: if the prior pulls B0 up, the fit may shrink the spread to keep the occupied sites right, and the identity would show a spread shortfall either way. So the run also changes each suspect on the same datasets.

## Design

All fits are two-stage, the `base` cell family (100 sites, 10 species, 2 primers, 3 PCR replicates per field sample, a spatial field, traits), at 20 field samples per site, with 2 chains of 1,000 draws after 1,000 burn-in, the study's settings.

-   **Diagnosis.** Replicates 1 to 40 of `M20` (ordinary, default prior), `M20_cm` (ordinary, community mean) and `base_rare_cm_m20` (rare-dominated, community mean), refitted by [refit_spread.R](refit_spread.R), which keeps what the runner discards. Each refit must reproduce its arm's stored rows exactly (every block, every element): `M20` from `../field-sample-ladder/results/simstudy-20261016-031826.rds`, the other two from `../community-mean-ladder-20261016/results/simstudy-20261017-025051.rds`. That shows the refits are the study's own fits. Before launch, replicate 1 of each arm was refitted this way as a check of the script, and its parameter and occupancy rows matched the stored ones exactly; its spread numbers were not looked at.
-   **B0-prior control (`sb0`).** The same 120 datasets refitted with `listPriors$sigma_b0 = 2`, a prior four times as wide in variance. If B0's prior drives the pattern, the pull should fall by about three quarters.
-   **Spread remedy.** Run only if the spread criterion below holds, on the same 120 datasets. The remedy is picked by the rule below: `d4`, 4 latent factors instead of the simulated 2 (`n_factors = 4`), or `gt4`, 4 latent traits instead of 2 (`n_lattrait = 4`).

Every variant keeps its arm's seed, so it fits exactly the diagnosis datasets, and every comparison between variants is paired by species-fit.

## What is measured

For every posterior draw and species, [refit_spread.R](refit_spread.R) rebuilds the fitted linear predictor from the stored draws of B0, B, U, L, the spatial coefficients and the range index, on the fit's own scaled covariates and spatial basis. The site average over draws of plogis of the rebuilt predictor must match the fit's posterior mean occupancy probability to within 1e-10. (Checked on one short fit of `base_rare_cm_m20` before writing this: largest difference 7e-16.) From the simulator it takes the true B0, the true covariate term, U times L and the spatial field.

Per draw and species, over the 100 sites:

-   the site mean and standard deviation of `r` and of each of its parts: covariate, latent-factor, spatial, and the site-level part (latent-factor plus spatial, compared with the truth's U times L plus field, since the fit can move variation between its factors and its spatial term);
-   the variance of the covariate part and of the site-level part, and their covariance;
-   the fitted expected occupancy, the site average of plogis(eta);
-   the correlation of `r` with the true `r`;
-   two solved intercepts, defined next.

The decomposition. Let `m` be the site mean of `r` and `a = B0 + m` the intercept at the site-average of `r`; the true expected occupancy is the site average of plogis(true eta). Then:

-   `aJ` is the intercept that, with this draw's centred `r`, gives the true expected occupancy;
-   `aS` is the intercept that gives it with the true centred `r` rescaled to this draw's standard deviation.

Averaging over draws, the posterior-mean B0 error `e` splits exactly as `e = J + V - D`:

-   `J = aJ - a_true`, the shape term: what the fitted spread and shape of `r` alone require of the intercept;
-   `V = a_fit - aJ`, the level term: the fitted expected occupancy differing from the truth;
-   `D = m_fit - m_true`, the location term: the fitted factors and spatial term absorbing a different site average.

The spread-only term `P = aS - a_true` is the part of `J` that the ratio of the fitted to the true standard deviation of `r` produces with the true shape kept. If `P` reproduces `e`, the spread shortfall predicts the B0 error.

The prior-pull term `Q` approximates the pull of B0's Normal prior: with `w` the posterior variance of B0 over sigma_b0 squared, `Q = w / (1 - w)` times (the prior centre minus the posterior mean of B0), with the centre at 0 under the default prior and at the posterior mean of `mu_b0` under the community mean. It is exact for a Normal likelihood and approximate here.

Also kept per fit: the posterior means and 95% intervals of sigma_b, the spatial amplitude sigma_bs, the factor-score scale sigma_h, the spatial range and `mu_b0`; the posterior means of `r` and its parts at each site; R-hat and bulk ESS of B0, the standard deviation of `r`, `aJ` and the site mean of `r`; and the study's own parameter and occupancy rows.

## Bands of true B0

As in the ladder report. Ordinary community: at most -1.5, -1.5 to -0.5, -0.5 to 0.5, 0.5 to 1.5, above 1.5. Rare-dominated: at most -4.25, -4.25 to -3.25, -3.25 to -2.25, -2.25 to -1.25, above -1.25. The extreme bands are, in the ordinary community, at most -1.5 and above 1.5; in the rare-dominated community, at most -4.25 and -4.25 to -3.25, where its error is largest.

## Does the spread shortfall predict the B0 error?

Written on 17 October 2026 before any fit ran. In an arm, the spread criterion holds if both:

-   **(a)** The least-squares slope of `e` on `P` across the arm's 400 species-fits is between 0.5 and 1.5, with a 95% interval, from standard errors clustered by replicate, that excludes 0.
-   **(b)** In each extreme band, the band mean of `P` has the sign of the band mean of `e` and at least half its size.

The answer is "yes" if the criterion holds in at least two of the three arms, "partly" if in one, "no" if in none. The REPORT also gives, per arm and band, the means of `e`, `J`, `V`, `D`, `P` and `Q`, the ratio of the fitted to the true standard deviation of `r` and of each part, and the regression of `e` on `P` and `Q` together.

## Choosing the remedy

Written before any fit ran, and applied by [choose.R](choose.R) after the diagnosis step, before the intervention step. If the spread answer is "no", no spread remedy runs and only the `sb0` control does. Otherwise, over the species-fits of the extreme bands of all three arms:

-   the covariate shortfall is the mean of the true variance of the covariate part minus its posterior-mean fitted variance;
-   the site-level shortfall is the same for the site-level part.

If the site-level shortfall is at least the covariate shortfall, the remedy is `d4` (4 latent factors): the truth's site-level variation has rank 4 across species (2 factors plus a field built from 2 underlying patterns), and the fit's 2 factors and smooth 20-knot spatial term may not hold it all. Otherwise it is `gt4` (4 latent traits), which frees each species' covariate effects from the measured-trait prediction. `choose.R` writes the choice to `draws/interventions.txt`, which the intervention step reads.

## When does an intervention work?

Written before any fit ran. Each is judged per arm, paired with the diagnosis refits of the same species-fits.

The spread remedy works in an arm if all four hold: B0 coverage rises by more than two paired standard errors (by replicate); the mean B0 error in each extreme band shrinks by at least a third; the fitted-to-true standard deviation ratio of `r` in the extreme bands moves towards 1; and, as a guard, the mean absolute error of species occupancy worsens by no more than 0.005 and the community occupancy error moves by no more than 0.01. It works overall if it works in at least two arms, one of them `base_rare_cm_m20`.

B0's prior is the main cause of an arm's pattern by true B0 if, under `sb0`, the mean B0 error in each extreme band falls by at least half; a minor contributor if it falls by less than a quarter in both; partial otherwise. The REPORT gives B0 coverage under `sb0` alongside.

Both interventions change what is estimated, so the REPORT also gives occupancy, `mu_b0`, B and the collection intercepts under each, against the diagnosis refits.

## Mixing

Decided before any fit ran. If, in any arm and variant, the median bulk ESS (of 2,000 draws) of B0 or of the standard deviation of `r` is under 100, or more than 10% of species-fits have R-hat above 1.05 for either, then replicates 1 to 5 of that arm and variant are refitted at 5,000 + 5,000 by `launch.sh long`, and the REPORT says how far each extreme band's mean `e`, `P` and standard deviation ratio move. The ladder's refits put B0's median bulk ESS at 144 to 216 of 2,000 in these arms.

## Running

Five fitting workers, one thread per fit, launched with `nohup` and `setsid` from [launch.sh](launch.sh) only when no other study, refit or R worker is running, on mains power with Low Power Mode off, under R 4.6.1 with Accelerate BLAS. Each M = 20 fit took about 115 s on 5 workers in the ladder's refits, so the diagnosis (120 fits) should take about 45 minutes and each intervention about the same. Steps, one at a time:

-   `launch.sh diag`: the 120 diagnosis refits;
-   `launch.sh choose`: `choose.R` applies the spread criterion and the remedy rule;
-   `launch.sh interventions`: `sb0` and the chosen remedy on the same 120 datasets;
-   `launch.sh long`: only if the mixing rule fires.

Per-fit files go to `draws/` (ignored); [analyse.R](analyse.R) combines them into the committed tables. Logs and provenance are kept beside them.

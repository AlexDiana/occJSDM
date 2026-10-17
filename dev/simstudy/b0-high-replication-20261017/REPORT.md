# Why B0 is miscalibrated at high field replication: a narrow linear predictor and B0's own prior

At 20 field samples per site occupancy is estimated well while B0 is not: the rarest species' B0 is overestimated, the commonest's underestimated, and B0's 95% intervals cover 0.87 (ordinary community) and 0.82 (rare-dominated community with B0's community mean) of the truths ([ladder report](../community-mean-ladder-20261016/REPORT.md)). This run asks why. It refits replicates 1 to 40 of `M20`, `M20_cm` and `base_rare_cm_m20`, rebuilds every posterior draw's occupancy linear predictor, splits each species' B0 error into exact terms, and then changes two suspects on the same 120 datasets. The candidate explanations, the decomposition, and the rules for every verdict below were written before any fit ran and are in [README.md](README.md).

## Answer

-   **The leading hypothesis holds: the fitted linear predictor is too narrow, and B0 compensates.** Rebuilt draw by draw, the spread of each species' predictor around B0 is 0.81 to 0.83 of the truth's (standard deviation 2.0 against 2.5) in all three arms, and it is short in every band of true B0. Because plogis is convex below 0 and concave above, a narrower predictor needs B0 nearer 0 to give the occupancy the data fix. That spread shortfall alone (the term `P`) predicts errors of the observed sign and about two thirds of their size in the extreme bands: +0.27 of +0.44 and -0.22 of -0.44 in `M20`, +0.87 of +1.30 for true B0 below -4.25 in the rare-dominated community. The slope of the error on `P` across species-fits is 0.97, 1.00 and 0.87. The README's spread criterion holds in all three arms.
-   **In the rare-dominated community the narrow predictor is most of the community-wide shift.** The mean B0 error there is +0.40, of which the shape term is +0.34. So most of `mu_b0`'s overshoot (+0.46 above the realised mean of the true B0 in these 40 replicates) is this compensation, not a failure of the community mean.
-   **The rest of the extreme-band error is B0's own prior.** The Normal(centre, 1) prior holds the expected occupancy of the rarest species 0.02 to 0.03 above the truth and of the commonest 0.02 to 0.03 below it (the level term `V`). Widening that prior to sd 2 (`sb0`) removes the level term and cuts the extreme-band errors by 44% to 70%. By the README's rule B0's prior is the main cause of the pattern in `M20_cm` and a partial cause in `M20` and `base_rare_cm_m20`. It barely changes coverage in the ordinary community (0.882 to 0.887, 0.875 to 0.882), lifts it in the rare-dominated one (0.818 to 0.868), widens the intervals by 11% to 19%, and leaves the mean absolute B0 error where it was or slightly worse.
-   **The pre-registered spread remedy, 4 latent factors instead of 2 (`d4`), does not work.** It widens the predictor (spread ratio 0.83 to 0.90 in the ordinary community, 0.81 to 0.86 in the rare-dominated one) and shrinks the shape term in the extreme bands by 9% to 37%, as the hypothesis predicts, but it leaves the level term, so the extreme-band errors fall by only 6% to 16%. B0 coverage rises by more than two standard errors only in the rare-dominated community (0.818 to 0.875).
-   **Why the predictor is narrow: the scale hyperparameters sit at their prior means.** The latent-factor term carries 0.40 to 0.49 of the truth's spread (0.60 to 0.71 with 4 factors), the spatial term 0.58 to 0.60 and the covariate term 0.89 to 0.90. The factor-score scale `sigma_h` is 0.34 to 0.39 in every arm and variant against a simulated 1, under an InvGamma(10, 1) prior on its square that puts it near 0.33; the spatial amplitude `sigma_bs` and `sigma_b` sit at 0.33 too. The priors of `sigma_h` and `sigma_b` are hard-coded; `sigma_bs` has an experimental half-Cauchy option, not tried here. Adding factors does not move `sigma_h`.
-   **Exploratory, not set out in the README: in the rare-dominated community the narrow predictor accounts for all of B0's undercoverage.** With the shape term subtracted, the errors are 0.99 posterior standard deviations in root mean square, and 95% of them fall within 1.96 posterior standard deviations. In the ordinary community it accounts for little (coverage 0.88 to 0.91 net of it). There the extra scatter is in the level and location terms, which the intervals also do not carry.
-   **Not a mixing artefact.** At 1,000 + 1,000 the spread's R-hat exceeds 1.05 in 20% to 35% of species-fits, so by the README's rule replicates 1 to 5 of every arm and variant were refitted at 5,000 + 5,000. Between the two lengths the B0 error and the shape term move by at most 0.03, and `sigma_h` stays at 0.34 to 0.38 with bulk ESS of 125 to 391 of 10,000. The short-chain results stand.

## Run

-   The 120 diagnosis refits ran on 17 October 2026 from 07:25 to 08:09 on 5 workers, at revision `f8b0698`. `choose.R` ran at 08:09, at `480f3ce`, which adds only `analyse.R`. The 240 intervention refits ran from 08:09 to 09:40 at `480f3ce`. The 45 long-chain refits ran from 09:40 to 11:02 at `480f3ce`.
-   No fit failed. Each took a median 100 to 118 s, or 537 s at 5,000 + 5,000. All runs used mains power with Low Power Mode off, R 4.6.1 with Accelerate BLAS, and one thread per fit. The diagnosis and choose steps launched on a clean tree; the later steps launched with one untracked file, `choose.csv`, the choose step's own output.
-   All 120 diagnosis refits reproduce their arms' stored rows exactly, every block and element (17,096 rows), and their occupancy rows ([reproduction.csv](reproduction.csv)). So these are the study's own fits, and the 40 replicates look like the 100: B0 coverage 0.882, 0.875 and 0.818 here against 0.873, 0.870 and 0.821.
-   Every rebuilt predictor reproduces its fit's posterior mean occupancy probabilities to within 5e-15.

## The predictor is narrow everywhere

From [decomposition.csv](decomposition.csv). For each draw and species, `r` is the predictor minus B0 across the 100 sites; the spread ratio is the posterior mean of its standard deviation over the truth's, averaged over species-fits. Parts are compared as ratios of their mean standard deviations, because species with no true latent-factor loading have a true latent-factor term of 0.

-   Spread ratio of `r`: 0.827 (`M20`), 0.824 (`M20_cm`), 0.806 (`base_rare_cm_m20`); variance ratio 0.77, 0.77, 0.75. In a correctly specified model, averaged over truths drawn from its priors, posterior draws have the truth's spread, so a shortfall in every arm and band points to misspecified priors, not to the usual shrinkage of posterior means.
-   By band it is shortest at the extremes: 0.78 and 0.77 in the ordinary community's outer bands against 0.86 in the middle one, and 0.70 for true B0 below -4.25 in the rare-dominated community against 0.86 in its top band.
-   Parts: covariate term 0.90, 0.90, 0.89; latent-factor term 0.49, 0.48, 0.40; spatial term 0.60, 0.58, 0.60; latent-factor and spatial terms together (the site-level part) 0.56, 0.55, 0.50.
-   The fitted predictor tracks the truth's pattern only partly: the correlation of `r` with the true `r` is 0.69, 0.69 and 0.65.

## Where each species' B0 error comes from

The README's exact split is `e = J + V - D`: the shape term `J` (what the fitted spread and shape of `r` require of B0 to give the true expected occupancy), the level term `V` (the fitted expected occupancy differing from the truth) and the location term `D` (the fitted latent-factor and spatial terms having a different site average from the truth's). `P` is the part of `J` produced by the standard deviation ratio alone. Values are means over the band's species-fits; the standard error of `e` treats each replicate as one cluster.

`M20` (ordinary community, default prior):

-   True B0 at most -1.5 (25 species-fits, mean -1.89): `e` +0.436 (SE 0.061) = `J` +0.267 + `V` +0.187 - `D` +0.019; `P` +0.268; coverage 0.76.
-   -0.5 to 0.5 (159): `e` +0.026 = `J` 0.000 + `V` +0.009 - `D` -0.017; coverage 0.92.
-   Above 1.5 (27, mean 1.89): `e` -0.443 (0.096) = `J` -0.227 + `V` -0.177 - `D` +0.039; `P` -0.225; coverage 0.82.

`M20_cm` (ordinary community, community mean):

-   At most -1.5: `e` +0.414 (0.060) = `J` +0.284 + `V` +0.175 - `D` +0.045; `P` +0.285; coverage 0.84.
-   Above 1.5: `e` -0.456 (0.093) = `J` -0.238 + `V` -0.169 - `D` +0.049; `P` -0.234; coverage 0.74.

`base_rare_cm_m20` (rare-dominated community, community mean):

-   True B0 at most -4.25 (27, mean -4.65): `e` +1.296 (0.099) = `J` +0.885 + `V` +0.345 - `D` -0.066; `P` +0.865; coverage 0.56.
-   -4.25 to -3.25 (99, mean -3.70): `e` +0.646 (0.080) = `J` +0.482 + `V` +0.132 - `D` -0.032; `P` +0.468; coverage 0.76.
-   -3.25 to -2.25 (145): `e` +0.390 = `J` +0.309 + `V` +0.076 - `D` -0.005.
-   -2.25 to -1.25 (94): `e` +0.122 = `J` +0.202 + `V` -0.066 - `D` +0.013.
-   Above -1.25 (35): `e` -0.156 = `J` +0.068 + `V` -0.194 - `D` +0.030.
-   All bands: `e` +0.404 = `J` +0.344 + `V` +0.051 - `D` -0.009.

Read together:

-   `P` and `J` agree to within 0.02 in every band, so the fitted shape of `r` matters only through its standard deviation.
-   `J` is the larger term at the extremes. In the middle of the ordinary community it vanishes, although the predictor is as narrow there, because B0 sits near 0 where the curvature cannot act.
-   `V` grows towards the extremes and changes sign with them. It is the fitted expected occupancy being 0.021 to 0.026 too high in the lowest bands and 0.023 to 0.026 too low in the highest; on the share of occupied sites, the occupancy error there is +0.011 to +0.014 and -0.006. That is what is left of the occupancy squeeze at 20 field samples per site.
-   `D` averages near 0 in every band, but it is not small species by species: in the ordinary community it accounts for 0.38 to 0.40 of the between-species variance of `e`, against 0.15 for `J` and 0.45 to 0.47 for `V`. In the rare-dominated community `J` accounts for 0.43, `V` for 0.43 and `D` for 0.14 ([variance-shares.csv](variance-shares.csv)).

## The spread criterion

From [spread-criterion.csv](spread-criterion.csv), the same numbers [choose.R](choose.R) used ([choose.csv](choose.csv)).

-   (a) Slope of `e` on `P`: `M20` 0.97 (95% interval 0.66 to 1.29), `M20_cm` 1.00 (0.72 to 1.28), `base_rare_cm_m20` 0.87 (0.72 to 1.01). R-squared 0.14, 0.15 and 0.37: the spread predicts the systematic part, and the rest is scatter.
-   (b) `P` as a share of `e` in the extreme bands: 0.61 and 0.51 (`M20`), 0.69 and 0.51 (`M20_cm`), 0.67 and 0.72 (`base_rare_cm_m20`), all of the same sign.
-   Both parts hold in all three arms, so the answer is "yes".
-   The approximate pull of B0's prior, `Q`, has band means of the same size as `P`'s (for example +0.31 and -0.27 in `M20`'s outer bands, +0.74 in the rare-dominated lowest band), but it does not predict `e` species by species: correlation 0.00, -0.05 and 0.04, and with `P` in the same regression its coefficient is -0.32 (SE 0.19), -0.55 (0.09) and 0.01 (0.02) ([regression.csv](regression.csv)). `Q` is computed from the estimated B0, so it is partly a function of the error itself; the intervention below is the cleaner test.

## The remedy rule

Over the 230 species-fits in the extreme bands of the three arms, the covariate part's variance falls 1.21 short of the truth and the site-level part's 1.29 short. The site-level shortfall is the larger, narrowly, so the rule picked `d4`. In the rare-dominated lowest band the covariate shortfall is the larger (2.48 against 1.35); the ordinary community's outer bands decided it.

## B0's prior widened (`sb0`)

`listPriors$sigma_b0 = 2` on the same 120 datasets, paired by species-fit with the diagnosis refits ([interventions.csv](interventions.csv)).

-   Extreme-band error, diagnosis to `sb0`: `M20` +0.436 to +0.181 (58% smaller) and -0.443 to -0.235 (47%); `M20_cm` +0.414 to +0.160 (61%) and -0.456 to -0.223 (51%); `base_rare_cm_m20` +1.296 to +0.724 (44%) and +0.646 to +0.192 (70%). By the README's rule (at least half in both extreme bands), B0's prior is the main cause in `M20_cm` and a partial cause in the other two.
-   What moves is the level term: in `M20`'s outer bands `V` goes from +0.187 to +0.031 and from -0.177 to -0.039, and in the rare-dominated lowest band from +0.345 to +0.007 (in the next band it overshoots, +0.132 to -0.189). The shape term falls by about a quarter (+0.267 to +0.195, -0.227 to -0.177, +0.885 to +0.675), because the predictor widens where B0 is freed: the extreme-band spread ratio goes from 0.78 to 0.83 in `M20` and from 0.70 to 0.77 in the rare-dominated lowest band.
-   B0 coverage: 0.882 to 0.887 (paired difference +0.005, SE 0.013), 0.875 to 0.882 (+0.007, 0.012) and 0.818 to 0.868 (+0.050, 0.016). Interval width rises from 1.40 to 1.57, 1.40 to 1.56 and 1.96 to 2.35. The mean absolute B0 error goes from 0.357 to 0.374, 0.359 to 0.366 and 0.583 to 0.587: less bias at the extremes, more variance elsewhere.
-   In the rare-dominated community the mean B0 error falls from +0.40 to +0.24, but `mu_b0`'s error does not move (+0.46 to +0.47): with a wider prior the species' B0 spread further from the community mean.
-   Occupancy barely moves: the mean absolute error of species occupancy changes by at most 0.0014 and the community occupancy error by at most 0.0021.

## Four latent factors (`d4`)

`n_factors = 4` on the same 120 datasets, against the simulated 2.

-   The predictor widens: spread ratio 0.827 to 0.899, 0.824 to 0.905 and 0.806 to 0.864; the latent-factor part from 0.49 to 0.70, 0.48 to 0.71 and 0.40 to 0.60. `sigma_h` stays at 0.37 to 0.38 (0.34 in the rare-dominated community), and its median bulk ESS halves (66 to 31).
-   The shape term shrinks as the hypothesis predicts (`M20` outer bands +0.267 to +0.222 and -0.227 to -0.175; rare-dominated lowest band +0.885 to +0.801), but the level term does not (+0.187 to +0.200, -0.177 to -0.196, +0.345 to +0.348).
-   So the extreme-band errors fall only a little: `M20` 9% and 8%, `M20_cm` 16% and 14%, `base_rare_cm_m20` 6% and 14%. By the README's rule (at least a third in both extreme bands, with coverage up by more than two standard errors), the remedy fails in every arm.
-   B0 coverage: 0.882 to 0.890 (+0.007, SE 0.008), 0.875 to 0.892 (+0.018, 0.010) and 0.818 to 0.875 (+0.058, 0.013). Intervals widen by 3% to 5%. B coverage rises from 0.894, 0.890 and 0.880 to 0.905, 0.910 and 0.897. Occupancy does not move. In the rare-dominated community the mean B0 error falls from +0.40 to +0.33 and `mu_b0`'s error from +0.46 to +0.39.

## Coverage, term by term (exploratory)

Not set out in the README. From [coverage-attribution.csv](coverage-attribution.csv): each term is subtracted from the error and the result compared with 1.96 posterior standard deviations (a Normal approximation of the intervals, which reproduces the study's coverage to within 0.015). `J` is systematic given the fit, so subtracting it reads as an attribution; `V` and `D` also carry ordinary estimation noise that a calibrated interval should cover, so their rows do not.

-   `base_rare_cm_m20`: coverage 0.805 (approximation), 0.950 net of `J`; root mean square error 1.54 posterior standard deviations, 0.99 net of `J`. The narrow predictor accounts for all of the undercoverage.
-   `M20` and `M20_cm`: 0.880 and 0.875, 0.910 and 0.900 net of `J`; root mean square 1.29 and 1.31, 1.21 and 1.22 net of `J`, 0.96 and 0.97 net of `J` and `V`. Here most of the excess is in `V` and `D` species by species, not in the curvature term.

## Mixing

From [mixing-check.csv](mixing-check.csv), [long-chains-summary.csv](long-chains-summary.csv), [long-chains-hyper.csv](long-chains-hyper.csv) and [long-chains.csv](long-chains.csv).

-   At 1,000 + 1,000, B0's median bulk ESS is 120 to 286 of 2,000 and the spread's 84 to 131, but the spread's R-hat exceeds 1.05 in 20% to 35% of species-fits and B0's in 6% to 22%. The README's rule fired for all nine arms and variants, so replicates 1 to 5 of each were refitted at 5,000 + 5,000 (45 fits, 50 species-fits per arm and variant).
-   At 10,000 draws, B0's R-hat exceeds 1.05 in at most 2% of species-fits, with median bulk ESS 433 to 739. The spread mixes more slowly: its R-hat exceeds 1.05 in 0% to 36% of species-fits (most in `M20_cm` with `sb0`), with median bulk ESS 155 to 477.
-   Between the two lengths, each arm and variant's mean B0 error moves by -0.023 to +0.009, the shape term by -0.027 to +0.007 and the level term by -0.002 to +0.007. The spread ratio moves by -0.01 to +0.08, mostly upwards in the ordinary arms, to 0.88 to 0.96 of the truth at 5,000 + 5,000. These 5 replicates are less short than the 40 (0.82 to 0.93 at 1,000 + 1,000), and B0 does not move with the spread.
-   `sigma_h` is 0.34 to 0.39 at 1,000 + 1,000 and 0.34 to 0.38 at 5,000 + 5,000, with median bulk ESS rising from 21 to 68 to 125 to 391 and R-hat at most 1.10. So its prior, not slow mixing, holds it there. `sigma_bs` is 0.33 at both lengths. The spatial range does not converge even at 10,000 draws (R-hat up to 1.78, median bulk ESS 6 to 51), as earlier runs found.
-   The short-chain results stand.

## What this does not show

-   **One design.** Two-stage model, 100 sites, 10 species, 20 field samples per site, a spatial field built from 2 underlying patterns and 2 latent factors, traits. The narrow predictor was not measured at fewer field samples, where the occupancy squeeze dominates.
-   **A fix.** Neither intervention restores B0's calibration in the ordinary community, and their combination was not run. The evidence points at the scale hyperpriors. Those of `sigma_h` and `sigma_b` cannot be changed without a package change, and `sigma_bs`'s half-Cauchy option was not tried; neither was part of this run.
-   **Which scale matters most.** The latent-factor term is the most short, but the fit can move variation between its factors and its spatial term, so the split between them is not identified. The covariate term is short mostly for the rarest species (0.73 of the truth below -4.25 in the rare-dominated community).
-   **Forty replicates.** The extreme bands hold 25 to 27 species-fits in the ordinary arms and 27 in the rare-dominated lowest band, so their means carry standard errors of 0.06 to 0.10.
-   **`Q` is an approximation.** It assumes a Normal likelihood for B0 and uses the estimated B0, so it is not an independent predictor; the `sb0` arm is the direct test.

## What this changes

-   **Reading B0.** At high replication B0 is the intercept of a predictor that is about 17% too narrow, so it sits too close to 0: rare species' B0 too high, common species' too low, and in a rare-dominated community the community mean too high. Occupancy itself is right. Users should read occupancy, not B0, and the documentation of B0 and of `b0_prior = "community"` should say so.
-   **The community mean is not the problem.** In the rare-dominated community most of `mu_b0`'s overshoot is the narrow predictor, which the default prior shares.
-   **`sigma_b0` is not a fix.** Widening B0's prior halves the extreme-band bias but trades it for variance, and does not restore coverage in an ordinary community.
-   **The next lever is the scale priors.** `sigma_h`, `sigma_bs` and `sigma_b` sit at their InvGamma(10, 1) prior means. Loosening them, `sigma_h` first (which needs a package option) and `sigma_bs` through its half-Cauchy option, alone and with `sigma_b0 = 2`, on these same 120 datasets, would test whether the predictor and B0 recover.

## Files and reproduction

-   [refit_spread.R](refit_spread.R) refits and rebuilds the predictor; [choose.R](choose.R) and [mixing_check.R](mixing_check.R) applied the README's rules; [launch.sh](launch.sh) ran every step; [common.R](common.R) holds what they share.
-   [analyse.R](analyse.R) writes [reproduction.csv](reproduction.csv), [species-fits.csv](species-fits.csv) (one row per species-fit and variant), [decomposition.csv](decomposition.csv), [variance-shares.csv](variance-shares.csv), [regression.csv](regression.csv), [spread-criterion.csv](spread-criterion.csv), [interventions.csv](interventions.csv), [hyper.csv](hyper.csv), [coverage-attribution.csv](coverage-attribution.csv), [long-chains.csv](long-chains.csv), [long-chains-hyper.csv](long-chains-hyper.csv) and [long-chains-summary.csv](long-chains-summary.csv).
-   The logs and provenance files record the run. The per-fit files (`draws/`) are not committed; `launch.sh diag`, `launch.sh interventions` and `launch.sh long` regenerate them exactly.

From the repository root:

```sh
Rscript dev/simstudy/b0-high-replication-20261017/analyse.R
```

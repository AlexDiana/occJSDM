# Community means for both intercepts: the overshoot goes, an undershoot takes its place

The community-mean run ([its report](../b0-community-mean-20261016/REPORT.md)) estimated the centre of B0's prior, `mu_b0`. In rare-dominated eDNA communities that removed most of the occupancy shift, but the collection intercepts' error reversed from -0.47 to +0.28 (two-stage) and from -0.37 to +0.22 (one-stage). The occupancy-collection split at community scale was then decided by the collection intercepts' own fixed Normal(0, 1) prior, centred above 98 to 99% of their true values. This run tests the companion change that report suggested: an experimental, off-by-default option, `listPriors = list(theta_intercept_prior = "community")`, that centres the collection intercepts' prior on an estimated community mean `mu_theta`. Six arms fitted with both community means at R = 100 are compared with default-prior and B0-only fits of identical datasets. The design and the success criteria, written before any fit ran, are in [README.md](README.md).

## Summary

-   **The community's occupancy comes out close to right, though no better than with B0's mean alone.** In the rare-dominated arms the community's mean occupancy error is +0.033 (SE 0.006) in `base_rare_cm_both` (two-stage) and +0.040 (SE 0.006) in `occupancy_rare_cm_both` (one-stage). That is against +0.167 and +0.171 with the default prior, and -0.050 and -0.013 with B0's mean alone. B0 coverage is 0.94 and 0.92 (default 0.42 and 0.39), and `mu_theta`'s 95% interval contains the realised mean collection intercept in 89% and 95% of replicates.
-   **The collection intercepts' overshoot goes, but an undershoot takes its place, so the success definition is not met.** Their error is -0.52 (SE 0.05) in `base_rare_cm_both`, as large as the default prior's -0.47, and -0.24 (SE 0.03) in `occupancy_rare_cm_both`, against +0.28 and +0.22 with B0's mean alone (C5 fails in both). `mu_theta` is estimated at -1.13 and -0.87 against a true mean of -0.64.
-   **The undershoot sits in the rarest species, and they are again estimated too common.** Species truly under 10% occupancy have collection intercepts 1.07 logit too low in the two-stage arm and 0.59 in the one-stage arm, and their B0 1.1 too high in both. Species truly at up to 5% of sites are overestimated by 0.14 in both arms, against 0.30 and 0.32 with the default prior and 0.03 and 0.06 with B0's mean alone (C3 fails in the two lowest bands). Each species' collection intercept has a Normal(`mu_theta`, 1) prior, while the true intercepts spread with sd about 0.4. With both centres free, nothing holds a rare species' collection rate near its community's, and B0's spread prior pulls its occupancy up towards the community mean.
-   **At community scale the data at 2 field samples per site do identify the split, though not tightly.** Widening both hyperpriors from sd 2 to sd 5 widens the intervals for `mu_b0` and `mu_theta` by only 4 to 12%, against up to 150% if the hyperpriors decided the split. It changes the community occupancy error by -0.027 and -0.012 and the collection-intercept error by +0.13 and +0.04 (C15 passes in both rare arms, narrowly in the two-stage one). At 5,000 + 5,000, on five replicates, the two-stage arm's changes shrink to -0.006 and -0.02, so much of that sensitivity was slow mixing. `mu_theta`'s posterior sd is 0.43 to 0.51 in the rare arms, against 2 for its hyperprior and 0.31 if every species' intercept were known. The sd-2 hyperprior still pulls the two-stage arm's `mu_b0` up, from -2.74 to -2.26 (truth -2.75; by 0.21 on the five long-chain replicates), about the ordinary shrinkage expected of a mean whose posterior sd is about 0.8.
-   **The parents lose a little.** Their community occupancy errors are small (+0.007 in `base_cm_both`, +0.018 in `occupancy_cm_both`). But species under 30% occupancy are overestimated by 0.03 to 0.07 more than with the default prior (C10 fails in both). Collection intercepts move by -0.19 and -0.12, with coverage 0.943 against 0.970 and 0.966 (C11 fails in both).
-   **Mixing gets worse.** At the study's 1,000 + 1,000, `mu_b0`'s median bulk ESS is 8 to 21 of 2,000 draws, against 21 to 29 with B0's mean alone, and its R-hat is under 1.05 in only 1 or 2 of 10 replicates per arm. `mu_theta` mixes well in the parents (median ESS 216 and 225) but not in the rare arms (23 and 75). B0 has R-hat over 1.05 in 47 to 82% of elements, against 30 to 42% with the default prior (C12 and C13 fail). At 5,000 + 5,000 both community means converge (`mu_b0`'s R-hat 1.02 to 1.05 on average, worst 1.12; `mu_theta`'s 1.01 to 1.02), and the averaged estimates move by at most 0.01 in occupancy and 0.07 in the collection intercepts, so the study's averages stand.

In short, the second community mean does what it was meant to do at community scale: no fixed prior centre decides the occupancy-collection split, and the data, not the hyperpriors, set its spread. The failure has moved down a level. Within a community, a rare species' own split is still decided by the fixed spreads of the species-level priors, and with both centres free that split goes the default prior's way: too common, collected too rarely.

## The model change

Each species' collection intercept, `beta_theta[1, s]`, gets a `Normal(mu_theta, 1)` prior in place of the fixed `Normal(0, 1)`, and `mu_theta` gets a `Normal(0, mu_theta_sd^2)` hyperprior with `mu_theta_sd` defaulting to 2, as `mu_b0_sd` does. After the collection coefficients are drawn, each iteration draws `mu_theta` from its conjugate normal conditional given the species' collection intercepts, and the next iteration's collection draw uses it as their prior mean. The collection intercepts' prior enters the sampler in one place, the prior-mean vector `b_betatheta` handed to the collection sampler, and no step integrates the collection coefficients out, so unlike B0's mean there is no second site to keep in agreement. The option is `listPriors = list(theta_intercept_prior = "community")`; this run uses it together with `b0_prior = "community"`. The community mean is called `mu_theta` because `theta0` is already the field-contamination rate.

One detail differs from the default sampler. `sample_betatheta_cpp_parallel()` leaves a species' collection coefficients unchanged in an iteration where the species is absent from every site, although with no collection data their conditional is their prior. Under the option a stale intercept would feed `mu_theta`'s draw, so such species are drawn from the prior the sampler used. Default fits do not take this path; whether they should is a separate question.

Draws are returned in `results_output$mu_theta_output`, a niter by nchain array present only when the option is on; `returnConvergenceDiagnostics()` reports them, `thinOutput()` thins them, and `infos$collection_intercept_prior` records the prior. Binary and continuous fits, which have no collection stage, reject the option. `?runOccJSDM` documents it as experimental and not a recommended setting.

## Tests

`tests/testthat/test-theta-community-prior.R`, all passing:

-   The option's parser returns nothing by default, so default fits record nothing new, and rejects malformed values, `mu_theta_sd` without the option, and the option in binary and continuous fits.
-   The conjugate `mu_theta` draw matches its closed form, to rounding error from a fixed seed and in mean and variance over 20,000 draws.
-   Both C++ collection samplers centre the intercept on the prior mean they are given: with one occupied site and two samples, repeated Polya-Gamma steps match the posterior mean from numerical integration at prior means of 0, -2.2 and 1.4. A sampler ignoring the mean could not pass.
-   Species absent from every site are drawn from their prior, other species are untouched, and with every species present no random number is drawn.
-   Within a fit, each iteration's collection draw receives the `mu_theta` drawn in the iteration before, verified against the saved draws; `mu_theta` is drawn from that iteration's intercepts with the sampler's own prior sd; the slopes keep a zero mean; and by default the sampler receives 0 throughout and nothing is drawn or saved.
-   Default fits are identical with and without `theta_intercept_prior = "fixed"`, two-stage and one-stage, spatial and not; with the option the output gains only `mu_theta_output`, and it combines with `b0_prior = "community"`.
-   A 30-species one-stage simulation with true collection intercepts of -2.2 to -1.1 recovers their community mean, and the community fit's intercepts lie below the fixed prior's.

That last test asserts less than the matching B0 test does, for a reason that bears on this report. On that dataset the community fit's average intercept error is -0.14 against the fixed prior's +0.14. Over 8 further seeds at 2 and at 4 samples per site the community fit was always the lower, its error averaging -0.25, and its interval for `mu_theta` always covered the realised mean. Species by species, the cause is B0's fixed Normal(0, 1) prior: it pulls each rare species' occupancy up, and that species' collection intercept falls to match. The fixed collection prior at 0 used to offset part of that; the community mean does not. When the intercepts are drawn from the model's own prior, `mu_theta` is close to calibrated (mean error -0.03 over 10 seeds), so this is the model, not the sampler.

The harness setting is tested in `tests/testthat/test-simstudy-theta-community.R`. The full suite passes: 2,058 expectations, 0 failures, and the same 5 skips as on the branch base (the env-gated coverage study and four parallel-chain identity tests that need a declared single-threaded BLAS).

**Fits that do not use the option are bit-identical to before, checked twice.** [check_default_identity.R](check_default_identity.R) fits twelve seeded models under the branch base (`cf9c133`) and under this change: two-stage, occupancy, binary and continuous, spatial and not, one with a non-default `sigma_b0`, one with an explicit `theta_intercept_prior = "fixed"`, and three B0-only community fits. All twelve whole fit objects are identical apart from `infos$execution`, which records timings and process ids, and a control fit that uses the option differs ([default-identity.csv](default-identity.csv)). At study scale, refitting replicates 1 to 3 of the eight comparison arms on this branch reproduced, element for element, the rows the paper study, the rare-dominated run and the community-mean run had stored ([reproduction.csv](reproduction.csv)).

## Run

The six arms ran from 12:06 to 13:30 on 16 October 2026 at revision `05ea03c`, with a clean tree at start and at write and `HEAD` unmoved: 600 fits, no failures, 5 replicate processes, one thread per fit, 2 chains of 1,000 draws after 1,000 burn-in, under R 4.6.1 with Accelerate BLAS and Low Power Mode off. The laptop was on battery for about 15 minutes around 13:00, with all five workers still at about 95% CPU. Each replicate seeds itself, so neither affects the results.

Every arm is paired with a default-prior arm and a B0-only arm on identical datasets: all truths and every species' true occupancy agree at all 100 replicates in all ten pairs ([pairing-check.csv](pairing-check.csv)). The default-prior arms are replicates 1 to 100 of `base` and `occupancy` from the R = 300 paper study and `base_rare` and `occupancy_rare` from the rare-dominated run; the B0-only arms are `base_cm`, `occupancy_cm`, `base_rare_cm` and `occupancy_rare_cm` from the community-mean run.

Two refits kept the draws that the runner discards:

-   Replicates 1 to 10 of the six arms at the study's settings, and replicates 1 to 3 of the eight comparison arms (11:56 to 12:06, 84 fits, no failures). Every refit reproduces its arm's stored rows exactly ([reproduction.csv](reproduction.csv)), so the mixing numbers describe the study's own fits.
-   Replicates 1 to 5 of the six arms at 5,000 + 5,000 (13:31 to 13:48, 30 fits, no failures). It runs because C12 failed. The README planned it for the four main arms; the two `_wide` arms were added after the first refits, because C15 compares them with the main rare arms and both mix slowly.

## Rare-dominated eDNA communities

Community-level numbers average each replicate's 10 species first ([community.csv](community.csv)); standard errors treat replicates as independent. For each arm the default prior comes first, then B0's mean alone, then both means:

-   `base_rare_cm_both` (two-stage). Mean occupancy 0.199 true; estimated 0.366, 0.149 and 0.232, error +0.167, -0.050 and +0.033 (SE 0.006). B0 error +1.89, -0.46 and +0.40. `mu_b0` -2.26 against -2.75 (coverage 0.91). B0 coverage 0.422, 0.917 and 0.936. Share of species estimated under 0.1: 0.000, 0.336 and 0.043, against 0.234 truly; under 0.2: 0.029, 0.753 and 0.406, against 0.544.
-   `occupancy_rare_cm_both` (one-stage). Mean occupancy 0.197 true; estimated 0.367, 0.184 and 0.237, error +0.171, -0.013 and +0.040 (SE 0.006). B0 error +1.87, -0.01 and +0.47. `mu_b0` -2.17 (coverage 0.84). B0 coverage 0.393, 0.940 and 0.919. Under 0.1: 0.002, 0.184 and 0.054, against 0.226; under 0.2: 0.032, 0.601 and 0.377, against 0.556.

By true occupancy band, paired species by species ([occupancy-by-band.csv](occupancy-by-band.csv)), bias with the default prior, with B0's mean alone, and with both means:

-   `base_rare_cm_both`. Up to 0.05: +0.300, +0.030 and +0.144. 0.05 to 0.1: +0.262, +0.004 and +0.106. 0.1 to 0.2: +0.191, -0.033 and +0.053. 0.2 to 0.3: +0.122, -0.071 and -0.001. 0.3 to 0.5: +0.054, -0.131 and -0.064.
-   `occupancy_rare_cm_both`. Up to 0.05: +0.317, +0.059 and +0.136. 0.05 to 0.1: +0.261, +0.040 and +0.104. 0.1 to 0.2: +0.200, +0.000 and +0.057. 0.2 to 0.3: +0.120, -0.031 and +0.008. 0.3 to 0.5: +0.046, -0.090 and -0.045.

Both means land between the other two priors. The community's average is about as close as with B0's mean alone, but the rarest species are pulled back up: half way from B0's mean alone towards the default in the lowest band. The slope of estimated on true occupancy is 0.40 and 0.47, against 0.27 and 0.22 with the default prior and 0.52 and 0.56 with B0's mean alone.

B0 by its true value ([b0-by-truth.csv](b0-by-truth.csv)) shows the same thing on the logit scale. In `base_rare_cm_both`, truths below -5 carry a bias of +2.98 (coverage 0.18), -5 to -4 +1.66 (0.79), -4 to -3 +0.91 (0.93), -3 to -2 +0.26 (0.98), -2 to -1 -0.45 and -1 to 0 -1.09 (0.86). The slope of estimate on truth is 0.26, against 0.34 with B0's mean alone and 0.14 with the default prior.

## The occupancy and collection trade-off

The collection-side truths are identical in each pair, so their change measures the fit's response ([detection-paired.csv](detection-paired.csv)):

-   Collection intercept: bias -0.471, +0.284 and -0.518 in `base_rare_cm_both`; -0.372, +0.216 and -0.243 in `occupancy_rare_cm_both`. Coverage 0.921, 0.981 and 0.932; 0.912, 0.969 and 0.956.
-   Collection slope: unbiased under all three priors; its intervals are 45% wider than the default's in both arms, a little narrower than with B0's mean alone.
-   `theta0`: bias -0.023, -0.005 and -0.007 in `base_rare_cm_both`, -0.018, 0.000 and -0.001 in `occupancy_rare_cm_both`; coverage 0.92 to 0.97 throughout.
-   `p` (`base_rare_cm_both` only): bias +0.066, +0.081 and +0.083, coverage 0.888, 0.851 and 0.849. The lab stage does not take up the change.

Species by species ([tradeoff-by-band.csv](tradeoff-by-band.csv)), collection-intercept bias with the default prior, with B0's mean alone and with both means:

-   `base_rare_cm_both`. True occupancy up to 0.1: -0.97, +0.15 and -1.07. 0.1 to 0.2: -0.49, +0.31 and -0.51. 0.2 to 0.3: -0.25, +0.36 and -0.27. 0.3 to 0.5: -0.13, +0.30 and -0.19.
-   `occupancy_rare_cm_both`. Up to 0.1: -0.78, +0.11 and -0.59. 0.1 to 0.2: -0.43, +0.22 and -0.28. 0.2 to 0.3: -0.22, +0.21 and -0.11. 0.3 to 0.5: +0.04, +0.36 and +0.09.

With B0's mean alone, the collection intercepts' error was a community-wide tilt, the same in every band. With both means it grows with rarity again, as it did under the default prior: in the two-stage arm it has the default's size in every band, though B0's error is far smaller (+1.13 against +2.81 in the lowest band). Each rare species trades its own few detections between occupancy and collection, and with no fixed centre on either side, what decides its trade is the two spread priors. B0 has sd 1 around `mu_b0`, which matches the simulator's spread and pulls a rare species' B0 towards the community's centre. The collection intercept has sd 1 around `mu_theta`, more than twice the true intercepts' spread of 0.39. Nothing pulls back as hard, so the species' collection rate falls to fit its detections at a higher occupancy. The same mechanism appeared in the package's own recovery test (see Tests above), with B0's fixed prior in place of its community mean.

On the probability scale, from the 10 refitted replicates of each arm ([collection-rate.csv](collection-rate.csv)), the community's mean collection rate (true 0.357 and 0.346) is estimated at 0.261, 0.413 and 0.240 in the two-stage arm and 0.308, 0.426 and 0.358 in the one-stage arm. The one-stage arm is close to right on this scale. Its logit-scale error comes from the rarest species, whose intercepts are poorly determined and estimated far down the logit scale, where a large logit error is a small error in the rate.

## Do the data at 2 field samples per site identify the split?

At community scale, yes, though not tightly. With both centres free, the only things that could decide where the community sits along the occupancy-collection trade-off are the data and the two hyperpriors. Both hyperpriors are centred at 0, so widening both by the same factor would leave an unidentified posterior centred near the same place but spread along the trade-off in proportion to the hyperprior sd. From the `_wide` arms against the sd-2 arms, paired by replicate ([wide-paired.csv](wide-paired.csv)):

-   `base_rare_cm_both_wide`. Interval widths grow by 10% for `mu_b0` and 12% for `mu_theta`, against 150% for the hyperprior. The community occupancy error moves from +0.033 to +0.006 (paired change -0.027, SE 0.005), the collection-intercept error from -0.52 to -0.39 (+0.13, SE 0.03), and `mu_b0` from -2.26 to -2.74 (truth -2.75; B0 error +0.40 to -0.03).
-   `occupancy_rare_cm_both_wide`. Widths grow by 6% and 4%. The community occupancy error moves from +0.040 to +0.028 (-0.012, SE 0.003), the collection-intercept error from -0.24 to -0.20 (+0.04, SE 0.01), and `mu_b0` from -2.17 to -2.40.

The intervals barely widen, so the data, not the hyperpriors, set the spread of the split, and C15 passes in both arms. The centre of `mu_b0` still moves, by half a logit in the two-stage arm and 0.23 in the one-stage arm. A sd-2 prior centred at 0 is expected to shift a mean that the data pin only to a posterior sd of about 0.8 by about 0.4, so this reads as ordinary shrinkage of a weakly determined mean, not as a prior choosing between equally good fits. The long chains (see Mixing) show that at convergence the two-stage arm's sensitivity is smaller still: on the five replicates refitted at 5,000 + 5,000, widening changes its community occupancy error by -0.006, its collection-intercept error by -0.02 and `mu_b0` by -0.21. The one-stage model, which observes each field sample directly, is a third to a half as sensitive to the hyperpriors as the two-stage model, whose field samples are seen only through PCR replicates.

From the refits ([identification.csv](identification.csv)), median over replicates: `mu_theta`'s posterior sd is 0.36 and 0.34 in the parents and 0.51 and 0.43 in the rare arms, against 2 for its hyperprior and about 0.31 if every species' intercept were known exactly. So the data remove 93 to 97% of the hyperprior's variance. Of what remains, the floor that 10 species set accounts for about three quarters in the parents and a third to a half in the rare arms. `mu_b0`'s posterior sd is 0.54, 0.50, 0.76 and 0.59. Within a fit, `mu_b0` and `mu_theta` are only weakly correlated (-0.21 to -0.42), while the community's mean B0 and mean collection intercept, draw by draw, correlate at -0.47 to -0.65: the trade-off is carried mainly by the species' deviations, not by the two community means. These sds come from chains that mix slowly, but the long chains give the same values to within 0.06.

At species scale, no. A rare species' own split is decided by its priors, as the previous section shows, and estimating the community means does not change that.

## The parents

Against `base` and `occupancy`, and their B0-only arms, on identical datasets ([community.csv](community.csv), [blocks.csv](blocks.csv)), default first, then B0's mean alone, then both means:

-   `base_cm_both`: B0 coverage 0.939, 0.927 and 0.937; community B0 error -0.14, -0.32 and +0.07; `mu_b0` +0.10 against 0 (coverage 0.95); `mu_theta` -0.81 against -0.64 (coverage 0.99); community occupancy error -0.023, -0.047 and +0.007; collection-intercept bias -0.002, +0.052 and -0.190, coverage 0.970, 0.971 and 0.943.
-   `occupancy_cm_both`: B0 coverage 0.935, 0.909 and 0.934; community occupancy error -0.011, -0.015 and +0.018; `mu_b0` +0.18 (coverage 0.92); `mu_theta` -0.71 (coverage 1.00); collection-intercept bias +0.042, +0.056 and -0.080, coverage 0.966, 0.965 and 0.943.
-   Occupancy by band: every band moves up. Species truly at 10 to 30% of sites are overestimated by 0.03 to 0.07 more than with the default prior (+0.206 against +0.139 for 0.1 to 0.2 in `base_cm_both`, from 18 species), and the commonest species are underestimated by 0.03 less.
-   B, G, `p` and `theta0` coverage stay within about 0.01 of the default, or improve (B 0.936 to 0.949 in `base_cm_both`).

In the rare arms B coverage rises too: 0.889 to 0.950 and 0.876 to 0.935.

## Mixing

From the refits of replicates 1 to 10 at the study's settings ([mixing-summary.csv](mixing-summary.csv), [mixing.csv](mixing.csv)), rank-normalised R-hat and bulk ESS out of 2,000 draws, with the default and B0-only figures from the community-mean run's refits of the same datasets:

-   `mu_b0`: median bulk ESS 12, 21, 8 and 15 in `base_cm_both`, `occupancy_cm_both`, `base_rare_cm_both` and `occupancy_rare_cm_both`, against 27, 29, 21 and 25 with B0's mean alone. R-hat under 1.05 in 1, 2, 1 and 1 of 10 replicates (B0 only: 5, 4, 2 and 3).
-   `mu_theta`: median bulk ESS 216, 225, 23 and 75; R-hat under 1.05 in 8, 10, 4 and 6 of 10.
-   B0, share of elements with R-hat over 1.05: 0.67, 0.47, 0.82 and 0.54, against 0.30, 0.31, 0.42 and 0.34 with the default prior and 0.50, 0.45, 0.50 and 0.45 with B0's mean alone. Median bulk ESS 30, 43, 14 and 37, against 61, 62, 44 and 55.
-   Collection intercepts, share with R-hat over 1.05: 0.27, 0.15, 0.51 and 0.22, against 0.10, 0.12, 0.31 and 0.21 with the default prior. Median bulk ESS 87, 106, 41 and 88, against 132, 135, 70 and 80.
-   The `_wide` arms mix about as the sd-2 arms do (`mu_b0` median ESS 12 and 21, `mu_theta` 27 and 73).

The second community mean adds a second slow direction. In the two-stage rare arm, `mu_b0` and B0 mix worst of any arm so far.

**Longer chains** ([long-chains-summary.csv](long-chains-summary.csv), [long-chains.csv](long-chains.csv)), replicates 1 to 5 at 5,000 + 5,000 against the same datasets at 1,000 + 1,000:

-   `mu_b0`'s R-hat falls to 1.02 to 1.05 on average (worst 1.12, in `base_rare_cm_both`), and its bulk ESS rises to 59 to 125 of 10,000 draws on average (minimum 14). `mu_theta`'s R-hat is 1.01 to 1.02 on average (worst 1.05), its ESS 140 to 1,170. The worst B0 R-hat in a fit averages 1.03 to 1.05. Both community means converge, `mu_b0` most slowly in the two-stage rare arm (mean ESS 59, against 100 to 132 with B0's mean alone at the same length).
-   Averaged over the five replicates, the community occupancy error changes by at most 0.010, the collection-intercept error by at most 0.07 and `mu_b0` by at most 0.11 in the four main arms, so the study's averages stand. Single replicates move more: occupancy by up to 0.04, `mu_b0` and the community's B0 error by up to 0.51.
-   In the two-stage `_wide` arm the long chains move `mu_b0` by +0.36 on average, and by 1.0 in one replicate, more than in any sd-2 arm. Paired with the sd-2 arm at 5,000 + 5,000, widening the hyperpriors changes the community occupancy error by -0.006, the collection-intercept error by -0.02 and `mu_b0` by -0.21, against -0.030, +0.17 and -0.50 for the same five replicates at 1,000 + 1,000. In the one-stage arm the two chain lengths agree (-0.011 against -0.013 for occupancy, -0.30 against -0.35 for `mu_b0`). So much of the two-stage arm's sensitivity to the hyperpriors at the study's chain length is slow mixing, and C15, measured there, overstates it.

So the conclusions do not rest on the short chains, but a user fitting one dataset with both options needs chains of at least the package default, 5,000 + 5,000, and should check `mu_b0`'s convergence.

## Against the success criteria

From [criteria.csv](criteria.csv) and [criteria-mixing.csv](criteria-mixing.csv):

-   C1, community occupancy error within ±0.05: passes in both rare arms (+0.033, +0.040).
-   C2, community B0 error within ±0.5: passes in both (+0.40, +0.47).
-   C3, low-band occupancy bias under +0.10: fails in the two lowest bands of both arms (+0.144 and +0.106; +0.136 and +0.104), passes in the 0.1 to 0.2 band.
-   C4, B0 coverage at least 0.85: passes (0.936, 0.919).
-   C5, collection-intercept bias within ±0.15: fails in both (-0.52, -0.24), and is not halved relative to B0's mean alone.
-   C6, `mu_b0` covers -2.75 in at least 85%: passes in `base_rare_cm_both` (0.91), fails narrowly in `occupancy_rare_cm_both` (0.84).
-   C14, `mu_theta` covers the realised mean in at least 85%: passes in all four arms (0.99, 1.00, 0.89, 0.95).
-   C8, parent B0 coverage no more than 0.02 lower: passes (0.937, 0.934).
-   C9, parent community occupancy error no more than 0.02 larger: passes (+0.007, +0.018).
-   C10, no occupancy band more than 0.03 worse: fails in both (+0.067 and +0.034 in their worst bands).
-   C11, detection side and B: fails on the collection intercept in both parents (bias moves by -0.19 and -0.12, coverage 0.943 against 0.970 and 0.966); `p`, `theta0` and B pass.
-   C12, community-mean mixing: fails in all four arms on `mu_b0`; `mu_theta` passes only in `occupancy_cm_both`.
-   C13, B0 and collection-intercept mixing relative to the default: fails on B0 in all four arms, and on the collection intercepts in the two two-stage arms.
-   C15, hyperprior sensitivity: passes in both rare arms, narrowly in `base_rare_cm_both` (-0.027 against 0.03, +0.13 against 0.15, +12% against 25%).
-   C16, descriptive: above.

By the rule written beforehand the shift is not removed, because C5 fails, and not reduced relative to B0's mean alone, because the collection-intercept error is not halved. The parents show small losses (C10, C11) and mixing fails (C12, C13). The data identify the split at community scale (C15).

## What this does not show

-   **One design.** 10 species, 100 sites, 2 field samples per site, 3 PCR replicates, true collection rates of 20 to 50%, and one B0 shift. With 10 species the community means themselves are uncertain (a floor of 0.31 on `mu_theta`'s sd); with more field samples per site each species' split is better identified, which is where the remaining error sits.
-   **No collection-spread change.** The collection intercepts' prior sd stays fixed at 1, and so does `sigma_b0`. The species-level result suggests the collection spread is now the prior that matters, but no fit here varies it.
-   **No mixed community, no continuous or binary model.** Binary fits have no collection stage.
-   **Point estimates for occupancy.** As before, the occupancy rows are posterior means with no interval.

## What this suggests

-   **For Alex: estimating both community means fixes the community-scale problem, not the species-scale one.** The occupancy-collection split at community scale is now set by the data, which pin it at 2 field samples per site to a posterior sd of about 0.4 to 0.8 logit in the rare arms. But within a community, each rare species' split is set by the two fixed spread priors. The collection intercepts' sd of 1 is too loose for these communities, so rare species are again estimated too common and collected too rarely: the default prior's error, at about half its size in occupancy.
-   **The natural next step is the collection intercepts' spread.** Estimating it (for example a half-normal or inverse-gamma prior on the collection intercepts' sd, conjugate given the intercepts) would let a community whose collection rates are similar hold each rare species' rate near the community's, which is what fixes its occupancy. A quick diagnostic before building that: refit the rare arms with the collection sd fixed at 0.5, which needs only the `prior_beta_theta_sd` constant made a setting.
-   **Mixing needs attention before any default.** Both community means add slow directions, and B0 mixes worst of any arm so far in the two-stage rare arm. A joint move along the B0 and collection-intercept trade-off, or non-centred intercepts, is worth trying before longer chains.
-   **Not a default.** Both options stay experimental and off; default fits are unchanged.

Should the next test estimate the collection intercepts' spread, and if so under what prior?

## Files and reproduction

-   [README.md](README.md): design and success criteria, written before the run.
-   [analyse.R](analyse.R) writes [pairing-check.csv](pairing-check.csv), [community.csv](community.csv), [occupancy-by-band.csv](occupancy-by-band.csv), [b0-by-truth.csv](b0-by-truth.csv), [detection-paired.csv](detection-paired.csv), [blocks.csv](blocks.csv), [tradeoff-by-band.csv](tradeoff-by-band.csv), [wide-paired.csv](wide-paired.csv) and [criteria.csv](criteria.csv) from `results/simstudy-20261016-132957.rds` and the three comparison runs' results.
-   [mixing.R](mixing.R) writes [reproduction.csv](reproduction.csv), [mixing.csv](mixing.csv), [mixing-summary.csv](mixing-summary.csv), [identification.csv](identification.csv), [collection-rate.csv](collection-rate.csv) and [criteria-mixing.csv](criteria-mixing.csv) from the refits' draws, and [long.R](long.R) writes [long-chains.csv](long-chains.csv) and [long-chains-summary.csv](long-chains-summary.csv). The draws themselves are not committed; [refit.R](refit.R), run by [launch.sh](launch.sh), regenerates them exactly, as [reproduction.csv](reproduction.csv) shows. The default and B0-only arms' probability-scale rows in [collection-rate.csv](collection-rate.csv) need the community-mean run's draws, named by `OLD_DRAWS`; without them `mixing.R` writes the both-means rows only.
-   [check_default_identity.R](check_default_identity.R) writes [default-identity.csv](default-identity.csv) from exports of the branch base and of the change.
-   `draws.log`, `main.log`, `long.log` and the three `*-provenance.txt` files record the three steps.

Run from the repository root:

```sh
Rscript dev/simstudy/theta-community-mean-20261016/analyse.R
OLD_DRAWS=<community-mean run>/draws Rscript dev/simstudy/theta-community-mean-20261016/mixing.R
Rscript dev/simstudy/theta-community-mean-20261016/long.R
```

# Rare-dominated communities at 3 field samples per site: B0's community mean alone is enough

The both-means run ([its report](../theta-community-mean-20261016/REPORT.md)) and the community-mean run before it ([its report](../b0-community-mean-20261016/REPORT.md)) fitted rare-dominated eDNA communities with 2 field samples per site. At M = 2, B0's community mean alone removed most of the occupancy shift but left the collection intercepts 0.28 and 0.22 logit too high. Adding a community mean for the collection intercepts turned that into an undershoot of 0.52 and 0.24. This run repeats the comparison at 3 field samples per site, a more realistic design, under all three priors: the default, B0's community mean alone (`b0_prior = "community"`), and both community means (`b0_prior` and `theta_intercept_prior` both `"community"`). The design and the success criteria, written before any fit ran, are in [README.md](README.md).

## Summary

-   **At M = 3, B0's community mean alone meets every criterion set for it.** In the two-stage arm `base_rare_cm_m3` the community's mean occupancy error is -0.022 (SE 0.003), its B0 error -0.11 and its collection-intercept bias +0.141 (SE 0.018). In the one-stage arm `occupancy_rare_cm_m3` they are -0.007, +0.07 and +0.135 (SE 0.016). C1, C2 and C5 hold in both, so by the rule written beforehand B0's mean alone suffices. C3, C4 and C6 hold too: species truly at up to 5% of sites are overestimated by only 0.06, B0 coverage is 0.94, and `mu_b0` is estimated at -2.75 and -2.56 against a generating -2.75. C5 holds with little room: the limit is 0.15.
-   **Both community means fail at M = 3, as at M = 2.** The community's occupancy is close (+0.031 in both arms), but the collection intercepts still undershoot (-0.38, SE 0.03, and -0.22, SE 0.03; C5 fails), and the undershoot does not halve (C17 fails). Species truly at up to 5% of sites are overestimated by 0.13 (C3 fails).
-   **The default prior's bias barely moves with a third field sample.** The community occupancy error is +0.158 and +0.154, against +0.167 and +0.171 at M = 2, and species at up to 5% of sites are still estimated about 0.31 to 0.34 too high. B0 coverage stays at 0.40 and 0.43. A third sample does not rescue the fixed B0 prior.
-   **The third sample helps where the collection prior's centre was deciding the split.** Paired with the M = 2 fits of the same occupancy truth, B0's mean alone gains most: its collection overshoot halves (+0.28 to +0.14 two-stage, +0.22 to +0.13 one-stage) and its B0 error falls from -0.46 to -0.11 in the two-stage arm. Under both means the third sample changes little (collection -0.52 to -0.38, -0.24 to -0.22), because there the error comes from each rare species' own split. The collection intercepts' sd of 1 still lets a rare species' collection rate fall while its occupancy rises.
-   **Mixing improves with the third sample but still fails at the study's chain length.** With B0's mean alone, `mu_b0` has R-hat under 1.05 in 6 and 5 of 10 replicates (2 and 3 at M = 2) and median bulk ESS 36 and 53 of 2,000 draws (C12 fails). B0 mixes about as under the default prior in the one-stage arm and somewhat worse in the two-stage arm (C13 fails on R-hat there). At 5,000 + 5,000 `mu_b0` converges in all four community-mean arms (R-hat 1.01 to 1.03 on average, worst 1.05), and the averaged estimates move by under 0.007 in occupancy and at most 0.05 in the collection intercepts, so the study's averages stand.

In short, at a realistic 3 field samples per site the simpler option is the one that works. Estimating B0's community mean moves a rare community's occupancy to the right level, and the third field sample gives the collection intercepts enough information to stay near their truth despite their fixed prior. Estimating the collection intercepts' community mean as well frees each rare species' collection rate to fall, and the rarest species come out too common again.

## Design and run

Six arms at R = 100: `base_rare_m3` and `occupancy_rare_m3` (default prior), `base_rare_cm_m3` and `occupancy_rare_cm_m3` (B0's mean alone), and `base_rare_cm_both_m3` and `occupancy_rare_cm_both_m3` (both means), with true B0 drawn around -2.75 and 3 field samples per site. Within each model the three priors fit identical datasets. Each arm takes the seed label of its M = 2 counterpart, so B0, the occupancy predictor, the true occupancy states, `p`, `q`, `theta0` and the collection intercepts are identical to the M = 2 datasets; only the field samples, the collection covariate and its slopes differ ([pairing-check.csv](pairing-check.csv)). Comparisons across M are therefore paired on occupancy and B0, while each side's collection intercepts are scored against their own truth on the fit's centred scale.

The main run took 95 minutes, from 14:10 to 15:46 on 16 October 2026, at revision `168ec9e`, with a clean tree at start and at write and `HEAD` unmoved: 600 fits, no failures, 5 replicate processes, one thread per fit, 2 chains of 1,000 draws after 1,000 burn-in, under R 4.6.1 with Accelerate BLAS, on AC power with Low Power Mode off. Two refits kept the draws:

-   Replicates 1 to 10 of the six arms at the study's settings (14:01 to 14:10, 60 fits, no failures). Every refit reproduces its arm's stored rows exactly ([reproduction.csv](reproduction.csv)).
-   Replicates 1 to 5 of the four community-mean arms at 5,000 + 5,000 (15:46 to 15:59, 20 fits, no failures), because C12 failed.

## B0's community mean alone

For each model, default prior first, then B0's mean alone, at M = 3 ([community.csv](community.csv)):

-   `base_rare_cm_m3` (two-stage). Mean occupancy 0.199 true; estimated 0.358 and 0.177. B0 error +1.82 and -0.11, coverage 0.401 and 0.944. `mu_b0` -2.75 (coverage 0.93). Collection-intercept bias -0.445 and +0.141. Share of species estimated under 0.1: 0.001 and 0.225, against 0.234 truly; under 0.2: 0.051 and 0.636, against 0.544.
-   `occupancy_rare_cm_m3` (one-stage). Mean occupancy 0.197 true; estimated 0.351 and 0.190. B0 error +1.76 and +0.07, coverage 0.425 and 0.937. `mu_b0` -2.56 (coverage 0.92). Collection-intercept bias -0.364 and +0.135. Under 0.1: 0.001 and 0.175, against 0.226; under 0.2: 0.062 and 0.591, against 0.556.

By true occupancy band, on the same species ([occupancy-by-band.csv](occupancy-by-band.csv)), bias with the default prior and with B0's mean alone:

-   `base_rare_cm_m3`. Up to 0.05: +0.310 and +0.056. 0.05 to 0.1: +0.245 and +0.024. 0.1 to 0.2: +0.183 and -0.013. 0.2 to 0.3: +0.110 and -0.044. 0.3 to 0.5: +0.048 and -0.080.
-   `occupancy_rare_cm_m3`. Up to 0.05: +0.341 and +0.064. 0.05 to 0.1: +0.251 and +0.034. 0.1 to 0.2: +0.164 and -0.001. 0.2 to 0.3: +0.096 and -0.021. 0.3 to 0.5: +0.056 and -0.063.

The commoner species in these communities are pulled down a little towards the community's centre, as at M = 2, but less: -0.080 in the 0.3 to 0.5 band against -0.131 at M = 2. The slope of estimated on true occupancy is 0.62 and 0.64, against 0.52 and 0.56 at M = 2 and 0.26 and 0.25 with the default prior.

The collection side ([detection-paired.csv](detection-paired.csv), [tradeoff-by-band.csv](tradeoff-by-band.csv)): the collection intercepts' coverage is 0.973 and 0.974, against 0.870 and 0.887 with the default prior, with intervals 26% and 27% wider. By true occupancy band their bias is -0.03 for species under 10% occupancy and +0.18 to +0.23 above it in the two-stage arm (+0.15 to +0.36 across the bands at M = 2), and -0.03 and +0.18 to +0.20 in the one-stage arm. What remains is the community-wide tilt the M = 2 run described: with B0's centre free, the collection intercepts' fixed Normal(0, 1) prior, centred above nearly all of their true values, still pulls the community's collection rate up. A third field sample gives the data enough weight to hold it to about +0.14. `theta0` and `p` move by under 0.02 (`p` bias +0.053 against +0.044), and B coverage rises from 0.882 to 0.917 and from 0.884 to 0.928.

On the probability scale, from the 10 refitted replicates of each arm ([collection-rate.csv](collection-rate.csv)), the community's mean collection rate (true 0.356 and 0.348) is estimated at 0.274 with the default prior, 0.378 with B0's mean alone and 0.278 with both means in the two-stage arm, and at 0.304, 0.395 and 0.345 in the one-stage arm.

## Both community means

Against the default prior and B0's mean alone on the same M = 3 datasets:

-   `base_rare_cm_both_m3`. Community occupancy error +0.031 (SE 0.005); B0 error +0.40, coverage 0.916; `mu_b0` -2.26 (coverage 0.86); `mu_theta` -1.00 against -0.64 (coverage 0.94); collection-intercept bias -0.377. Occupancy bias up to 0.05: +0.132; 0.05 to 0.1: +0.089; 0.1 to 0.2: +0.043.
-   `occupancy_rare_cm_both_m3`. Community occupancy error +0.031 (SE 0.004); B0 error +0.42, coverage 0.921; `mu_b0` -2.22 (coverage 0.83); `mu_theta` -0.84 (coverage 0.99); collection-intercept bias -0.215. Occupancy bias up to 0.05: +0.131; 0.05 to 0.1: +0.080; 0.1 to 0.2: +0.041.

The pattern is the M = 2 one, slightly smaller. Species under 10% occupancy have collection intercepts 0.83 and 0.62 logit too low and B0 1.0 too high (at M = 2: 1.07 and 0.59, and 1.1), so the error still grows with rarity. `mu_b0` sits 0.5 too high in both arms, as at M = 2. Against B0's mean alone on the same datasets, both means put species under 10% occupancy 0.05 to 0.08 higher in occupancy and pull every collection intercept down: an extra fitted parameter that the data cannot pin species by species gives the wrong answer for the rarest species.

## The third field sample

Each M = 3 arm against its M = 2 counterpart, paired by replicate on the same B0 and true occupancy ([m-effect.csv](m-effect.csv), [occupancy-by-band.csv](occupancy-by-band.csv)), M = 2 first:

-   Default prior: community occupancy error +0.167 and +0.158 (paired change -0.008, SE 0.003) two-stage, +0.171 and +0.154 (-0.017, SE 0.003) one-stage; collection intercepts -0.47 and -0.44, -0.37 and -0.36; B0 coverage 0.42 and 0.40, 0.39 and 0.43. Every occupancy band moves by under 0.04.
-   B0's mean alone: community occupancy error -0.050 and -0.022 (+0.028, SE 0.004), -0.013 and -0.007 (+0.006, SE 0.004); collection intercepts +0.28 and +0.14 (-0.14, SE 0.02), +0.22 and +0.13 (-0.08, SE 0.02); B0 error -0.46 and -0.11, -0.01 and +0.07; `mu_b0` -3.10 and -2.75, -2.64 and -2.56.
-   Both means: community occupancy error +0.033 and +0.031, +0.040 and +0.031; collection intercepts -0.52 and -0.38 (+0.14, SE 0.05), -0.24 and -0.22 (+0.03, SE 0.04); B0 error +0.40 and +0.40, +0.47 and +0.42.

So a third field sample does most where a prior centre on the collection side was deciding the community-level split, which is the B0-only fit. It does little for the default prior, whose error comes from B0's fixed centre, and little for both means, whose error comes from the species-level spreads. The field-sample ladder (branch `claude/field-sample-ladder`) found that in ordinary communities the default prior's squeeze is undone mostly by 10 samples per site; going from 2 to 3 here is a small step on that ladder.

## Mixing

From the refits of replicates 1 to 10 at the study's settings ([mixing-summary.csv](mixing-summary.csv), [mixing.csv](mixing.csv)), rank-normalised R-hat and bulk ESS out of 2,000 draws, M = 2 figures from the earlier runs' refits of the paired datasets:

-   `mu_b0`, B0's mean alone: R-hat under 1.05 in 6 and 5 of 10 replicates (M = 2: 2 and 3); median bulk ESS 36 and 53 (M = 2: 21 and 25).
-   `mu_b0`, both means: 1 and 7 of 10 (M = 2: 1 and 1); median ESS 15 and 34 (8 and 15). `mu_theta`: 6 and 10 of 10 (4 and 6), median ESS 53 and 165 (23 and 75).
-   B0, share of elements with R-hat over 1.05: 0.31 and 0.27 with the default prior, 0.43 and 0.31 with B0's mean alone, 0.58 and 0.35 with both means. Median bulk ESS 47 and 75, 49 and 62, 36 and 49.
-   Collection intercepts: B0's mean alone mixes them better than the default prior does (median ESS 146 and 161 against 84 and 107; share with R-hat over 1.05 0.07 and 0.11 against 0.23 and 0.12).

**Longer chains** ([long-chains-summary.csv](long-chains-summary.csv), [long-chains.csv](long-chains.csv)), replicates 1 to 5 of the four community-mean arms at 5,000 + 5,000 against the same datasets at 1,000 + 1,000:

-   `mu_b0`'s R-hat falls to 1.01 to 1.03 on average (worst 1.05, in `occupancy_rare_cm_m3`) and its bulk ESS rises to 120 to 138 of 10,000 draws on average (minimum 45). `mu_theta`'s R-hat is 1.00 and 1.02 on average (worst 1.06). The worst B0 R-hat in a fit averages 1.03 to 1.04.
-   Averaged over the five replicates, the community occupancy error changes by under 0.007, the collection-intercept error by at most 0.05 and `mu_b0` by at most 0.08 in every arm. Single replicates move more: occupancy by up to 0.03, `mu_b0` by up to 0.53.

So the conclusions do not rest on the short chains, but a user fitting one dataset with `b0_prior = "community"` needs chains of at least the package default, 5,000 + 5,000, and should check `mu_b0`'s convergence.

## Against the success criteria

From [criteria.csv](criteria.csv) and [criteria-mixing.csv](criteria-mixing.csv).

Both means (`base_rare_cm_both_m3`, `occupancy_rare_cm_both_m3`):

-   C1, community occupancy error within ±0.05: passes (+0.031, +0.031).
-   C2, community B0 error within ±0.5: passes (+0.40, +0.42).
-   C3, low-band occupancy bias under +0.10: fails in the lowest band of both (+0.132, +0.131); passes in the next two.
-   C4, B0 coverage at least 0.85: passes (0.916, 0.921).
-   C5, collection-intercept bias within ±0.15: fails (-0.38, -0.22).
-   C6, `mu_b0` covers -2.75 in at least 85%: passes in the two-stage arm (0.86), fails in the one-stage arm (0.83).
-   C14, `mu_theta` covers the realised mean in at least 85%: passes (0.94, 0.99).
-   C17, collection-intercept error at most half its M = 2 value: fails (0.38 against 0.26; 0.22 against 0.12).

B0's mean alone (`base_rare_cm_m3`, `occupancy_rare_cm_m3`):

-   C1: passes (-0.022, -0.007). C2: passes (-0.11, +0.07). C3: passes in all three bands of both (at most +0.064). C4: passes (0.944, 0.937). C5: passes (+0.141, +0.135). C6: passes (0.93, 0.92).

Mixing:

-   C12: fails in all four community-mean arms at 1,000 + 1,000; `mu_theta` passes in `occupancy_rare_cm_both_m3`.
-   C13: B0's mean alone passes on B0's ESS and on the collection intercepts in both arms, and on B0's R-hat in the one-stage arm; it fails on B0's R-hat in the two-stage arm (0.43 against 0.31). Both means fail on B0's R-hat in both arms and on B0's ESS in the one-stage arm.
-   C18, descriptive: above.

By the rules written beforehand, both means do not remove the shift at M = 3 (C5 fails), and B0's mean alone suffices (C1, C2 and C5 hold in both arms). Neither passes the mixing criteria at the study's chain length.

## What this does not show

-   **No parents at M = 3.** At M = 2, B0's mean alone cost ordinary communities a little (B0 coverage 0.935 to 0.909 in `occupancy_cm`, community occupancy error -0.023 to -0.047 in `base_cm`). Whether that cost remains at M = 3 is untested, and it decides whether the option could become a default.
-   **One design otherwise.** 10 species, 100 sites, 3 PCR replicates on 2 primers, true collection rates of 20 to 50%, and one B0 shift. C5 passes with little room (+0.141 and +0.135 against 0.15): communities whose collection rates sit further below 50% would feel the collection prior's pull more.
-   **No mixed community.** A few common species among many rare ones would show how far one community mean drags the common species down; the 0.3 to 0.5 band (-0.08) suggests a little.
-   **Point estimates for occupancy.** As before, the occupancy rows are posterior means with no interval.

## What this suggests

-   **For Alex: at a realistic 3 field samples per site, B0's community mean alone is the candidate.** It removes the rare-community occupancy shift by every criterion written beforehand, and adding the collection community mean makes the rarest species worse. The collection-mean option can stay experimental, or be dropped.
-   **What is left before it could be a default:** the parents at M = 3 (does the small cost to ordinary communities remain?), and mixing (`mu_b0` needs at least 5,000 + 5,000 and a convergence check).
-   **The collection intercepts' fixed prior centre still tilts the community's collection rate up by about 0.14 logit.** A weakly informative collection prior centred nearer typical eDNA collection rates, or a wider fixed sd, would reduce that without the instability of a free community mean. That is a smaller question than the ones above.

Should the next test run the parents at M = 3 with B0's community mean, to see whether it could become the default?

## Files and reproduction

-   [README.md](README.md): design and success criteria, written before the run.
-   [analyse.R](analyse.R) writes [pairing-check.csv](pairing-check.csv), [community.csv](community.csv), [occupancy-by-band.csv](occupancy-by-band.csv), [b0-by-truth.csv](b0-by-truth.csv), [detection-paired.csv](detection-paired.csv), [blocks.csv](blocks.csv), [tradeoff-by-band.csv](tradeoff-by-band.csv), [m-effect.csv](m-effect.csv) and [criteria.csv](criteria.csv) from `results/simstudy-20261016-154536.rds` and the three M = 2 runs' results.
-   [mixing.R](mixing.R) writes [reproduction.csv](reproduction.csv), [mixing.csv](mixing.csv), [mixing-summary.csv](mixing-summary.csv), [identification.csv](identification.csv), [collection-rate.csv](collection-rate.csv) and [criteria-mixing.csv](criteria-mixing.csv) from the refits' draws, and [long.R](long.R) writes [long-chains.csv](long-chains.csv) and [long-chains-summary.csv](long-chains-summary.csv). The draws are not committed; [launch.sh](launch.sh) regenerates them with the both-means run's [refit.R](../theta-community-mean-20261016/refit.R).
-   `draws.log`, `main.log`, `long.log` and the three `*-provenance.txt` files record the three steps.

Run from the repository root:

```sh
Rscript dev/simstudy/rare-m3-20261016/analyse.R
Rscript dev/simstudy/rare-m3-20261016/mixing.R
Rscript dev/simstudy/rare-m3-20261016/long.R
```

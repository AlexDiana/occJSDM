# Ordinary communities at 3 field samples per site: B0's community mean costs nothing measurable

At 3 field samples per site, the rare-M3 run ([its report](../rare-m3-20261016/REPORT.md)) found that `b0_prior = "community"` alone removes the rare-dominated community's occupancy shift by every criterion set for it. At 2 field samples per site the community-mean run ([its report](../b0-community-mean-20261016/REPORT.md)) had found that the same option cost ordinary communities a little: three of its no-harm criteria failed narrowly. This run refits the two eDNA parents, `base` (two-stage) and `occupancy` (one-stage), with true B0 drawn from the simulator's Normal(0, 1), at 3 field samples per site under the default prior and under B0's community mean, to decide whether the option could become the default. The design and the success criteria, written before any fit ran, are in [README.md](README.md).

## Summary

-   **At M = 3, B0's community mean passes every no-harm criterion in both parents, with wide margins.** Against the default prior on identical datasets, B0 coverage is 0.937 against 0.938 (two-stage) and 0.931 against 0.930 (one-stage). The community occupancy error is -0.015 against -0.009, and +0.000 against -0.002. No occupancy band is more than 0.006 worse, against a limit of 0.03. The collection-intercept bias moves by +0.017 and +0.001, against a limit of 0.05. `p` and `theta0` move by at most 0.003 in coverage and 0.001 in bias, and B coverage rises slightly. No criterion uses more than about a third of its allowance; the closest are the two-stage community occupancy error (0.006 against 0.02) and collection-intercept bias (0.017 against 0.05). C8 to C11 hold in both arms.
-   **The small costs seen at M = 2 go away with the third field sample.** On the same replicates, the option's cost to the two-stage community occupancy error falls from -0.024 to -0.006 (paired change +0.018, SE 0.005). Its cost to the collection-intercept bias falls from +0.053 to +0.017. In the one-stage arm, its cost to B0 coverage falls from -0.026 to +0.001 (paired change +0.027, SE 0.011). At M = 2 the two-stage fit placed `mu_b0` at -0.29 on average, against a realised mean of the true B0 of +0.03, which pulled every species down. At M = 3 it places `mu_b0` at -0.05, and its interval covers the generating mean of 0 in 96% of replicates (88% at M = 2).
-   **What remains is spread, not bias, and the intervals account for it.** The default prior is told the community's centre, which here is the truth; the community mean has to estimate it from 10 species (posterior sd about 0.45). So the community-level estimates vary more between datasets: the per-replicate community occupancy error has sd 0.044 against 0.028 (two-stage) and 0.042 against 0.030 (one-stage). B0's 95% intervals are 6 to 7% wider (3.02 against 2.82, 2.82 against 2.65) and their coverage is unchanged.
-   **Mixing still fails at the study's chain length.** `mu_b0` has R-hat under 1.05 in 5 and 8 of 10 replicates, with median bulk ESS 38 and 56 of 2,000 draws (C12 fails in both arms). In the two-stage arm B0 mixes worse than under the default prior (median ESS 48 against 71; share with R-hat over 1.05 0.43 against 0.25), so C13 fails there; in the one-stage arm it passes. At 5,000 + 5,000 `mu_b0` converges (R-hat at most 1.03), and the averaged errors and the community-minus-default differences move by under 0.007 in occupancy, so the results stand.

So `b0_prior = "community"` is harmless to ordinary communities at 3 field samples per site. With the rare-M3 result, that makes it a candidate for the default prior on B0. It needs longer chains than the study's 1,000 + 1,000 for `mu_b0` to converge; at the package default of 5,000 + 5,000 it does.

## Design and run

Four arms at R = 100: `base_m3` and `occupancy_m3` (default prior), and `base_cm_m3` and `occupancy_cm_m3` (`b0_prior = "community"`), at 3 field samples per site. Within each model the two priors fit identical datasets. Each arm takes its parent's seed label, so B0, B, G, the true occupancy states, `p`, `q` and `theta0` are identical to those of the M = 2 parents; only the field samples, the collection covariate and its slopes differ ([pairing-check.csv](pairing-check.csv)). The M = 2 comparisons use replicates 1 to 100 of `base` and `occupancy` from the R = 300 paper study and `base_cm` and `occupancy_cm` from the community-mean run. Recomputed on these replicates, the M = 2 criteria reproduce the community-mean run's figures exactly.

The main run went from 16:28 to 18:35 on 16 October 2026 at revision `c7b6db5`, with a clean tree at launch: 400 fits, no failures, 5 replicate processes, one thread per fit, 2 chains of 1,000 draws after 1,000 burn-in, under R 4.6.1 with Accelerate BLAS and Low Power Mode off. The laptop was on battery at launch, its lid was closed at 16:30, and it slept until the lid was opened at 17:48, after which it ran on AC power. So of the run's 126 minutes, about 48 were spent fitting. Sleep pauses the fits without changing them: every refit below reproduces its stored rows exactly. Two refits kept the draws:

-   Replicates 1 to 10 of the four arms at the study's settings (16:23 to 16:28, 40 fits, no failures). Every refit reproduces its arm's stored rows exactly ([reproduction.csv](reproduction.csv)).
-   Replicates 1 to 5 of all four arms at 5,000 + 5,000 (18:35 to 18:48, 20 fits, no failures), because C12 failed. The default arms are included, so the community-minus-default comparison can be repeated at the longer chain length.

## No harm to ordinary communities

For each model, default prior first, then B0's community mean, at M = 3 ([community.csv](community.csv), [detection-paired.csv](detection-paired.csv), [blocks.csv](blocks.csv)):

-   Two-stage (`base_m3`, `base_cm_m3`). Mean occupancy 0.505 true; estimated 0.496 and 0.490. Community occupancy error -0.009 and -0.015 (paired difference -0.006, SE 0.003). B0 error -0.05 and -0.09, coverage 0.938 and 0.937, interval width 2.82 and 3.02. Collection-intercept bias -0.016 and +0.001, coverage 0.961 and 0.962. `p` bias +0.023 in both, coverage 0.937 in both. `theta0` coverage 0.973 and 0.972. B coverage 0.912 and 0.919.
-   One-stage (`occupancy_m3`, `occupancy_cm_m3`). Mean occupancy 0.504 true; estimated 0.502 and 0.504. Community occupancy error -0.002 and +0.000 (paired difference +0.002, SE 0.003). B0 error -0.04 and -0.03, coverage 0.930 and 0.931, interval width 2.65 and 2.82. Collection-intercept bias +0.007 and +0.008, coverage 0.965 and 0.962. `theta0` coverage 0.976 and 0.973. B coverage 0.915 and 0.919.

By true occupancy band, on the same species ([occupancy-by-band.csv](occupancy-by-band.csv)), bias with the default prior and with B0's community mean:

-   Two-stage. 0.1 to 0.2: +0.155 and +0.146. 0.2 to 0.3: +0.105 and +0.097. 0.3 to 0.5: +0.022 and +0.010. 0.5 to 0.7: -0.041 and -0.041. 0.7 to 0.9: -0.126 and -0.127.
-   One-stage. 0.1 to 0.2: +0.172 and +0.178. 0.2 to 0.3: +0.107 and +0.096. 0.3 to 0.5: +0.025 and +0.026. 0.5 to 0.7: -0.035 and -0.031. 0.7 to 0.9: -0.096 and -0.088.

Both priors leave the same within-community squeeze: species at 0.1 to 0.2 are estimated 0.15 to 0.18 too common and species at 0.7 to 0.9 0.09 to 0.13 too rare. The slope of estimated on true occupancy is 0.57 and 0.60 (two-stage) and 0.61 and 0.63 (one-stage), so B0's community mean neither causes the squeeze nor removes it. These parents have almost no species under 10% or over 90% occupancy (one species-fit in each in the one-stage arm, none in the two-stage arm), so C10 is judged on the five bands from 0.1 to 0.9. The rare end is where the rare-M3 run speaks.

## Why the M = 2 cost goes away

Each arm against its M = 2 counterpart, paired by replicate on the same B0 and true occupancy ([m-effect.csv](m-effect.csv), [mu-b0.csv](mu-b0.csv)), M = 2 first:

-   Default prior: community occupancy error -0.023 and -0.009 (two-stage), -0.011 and -0.002 (one-stage); B0 error -0.14 and -0.05, -0.09 and -0.04; B0 coverage 0.939 and 0.938, 0.935 and 0.930.
-   B0's community mean: community occupancy error -0.047 and -0.015, -0.015 and +0.000; B0 error -0.32 and -0.09, -0.12 and -0.03; B0 coverage 0.927 and 0.937, 0.909 and 0.931; collection-intercept bias +0.052 and +0.001, +0.056 and +0.008.
-   `mu_b0`, against a realised mean of the true B0 of +0.03 (two-stage) and +0.05 (one-stage): -0.29 and -0.05, -0.06 and +0.03. Its 95% interval covers the generating mean of 0 in 88% and 96% of replicates, 94% and 93%; it is 12% and 9% narrower at M = 3.

At M = 2 the two-stage fit tends to put the community's centre too low: `mu_b0` averages 0.32 logit below the realised mean of the species' true B0, B0 drops with it by 0.19 more than under the default prior, and the collection intercepts rise by 0.05 to compensate. This is the same occupancy-collection trade-off that the rare-dominated runs describe, here in a community where the fixed Normal(0, 1) prior is exactly the truth's distribution and so costs nothing. A third field sample pins each species' collection rate, so the split no longer has room to drift and `mu_b0` comes to within 0.09 of the realised mean.

In the one-stage arm `mu_b0` is close on average at M = 2 already (-0.06), and the B0 coverage loss there (0.935 to 0.909) comes from spread instead. The community-level B0 error varies more between replicates under the community mean (sd 0.47 against 0.26), so B0's RMS error rises from 0.78 to 0.90 while its intervals widen only from 2.92 to 3.16. At M = 3 the spread is smaller (sd 0.33 against 0.23, RMS error 0.76 against 0.72) and the wider intervals (2.82 against 2.65) cover it.

## Mixing

From the refits of replicates 1 to 10 at the study's settings ([mixing-summary.csv](mixing-summary.csv), [mixing.csv](mixing.csv), [identification.csv](identification.csv)), rank-normalised R-hat and bulk ESS out of 2,000 draws, M = 2 figures from the community-mean run's refits of the paired datasets:

-   `mu_b0`: R-hat under 1.05 in 5 and 8 of 10 replicates (M = 2: 5 and 4); worst 1.12 and 1.07. Median bulk ESS 38 and 56 (M = 2: 27 and 29); in three two-stage replicates it is under 20.
-   B0, default prior and community mean: median bulk ESS 71 and 48 (two-stage), 78 and 79 (one-stage); share of elements with R-hat over 1.05, 0.25 and 0.43, 0.20 and 0.20.
-   Collection intercepts: median ESS 141 and 135, 190 and 202; share with R-hat over 1.05, 0.14 and 0.10, 0.04 and 0.09.
-   Within a fit, `mu_b0`'s draws correlate at about 0.71 with the mean of the species' B0 draws, and that mean correlates at about -0.5 with the mean collection intercept. The chains move slowly along that direction: the community's centre and its collection rate trading against each other.

**Longer chains** ([long-chains-summary.csv](long-chains-summary.csv), [long-chains.csv](long-chains.csv), [long-cost.csv](long-cost.csv)), replicates 1 to 5 of all four arms at 5,000 + 5,000 against the same datasets at 1,000 + 1,000:

-   `mu_b0`'s R-hat falls to 1.02 and 1.01 on average (worst 1.03 and 1.01), and its bulk ESS rises to 194 and 367 of 10,000 draws on average (minimum 115). The worst B0 R-hat in a fit averages 1.03 and 1.02, the same as under the default prior; at 1,000 + 1,000 it averaged 1.13 and 1.08, against 1.13 and 1.07 under the default prior. So B0 is not fully converged under either prior at the study's chain length; what the community mean adds is `mu_b0`.
-   Averaged over the five replicates, the community occupancy error moves by under 0.007 in every arm, the collection-intercept error by at most 0.011, B0's error by at most 0.05 and `mu_b0` by at most 0.05. Single replicates move more: occupancy by up to 0.014, `mu_b0` by up to 0.12.
-   The community-minus-default difference on these five replicates barely changes with the longer chains: in community occupancy error it is +0.022 and +0.017 (two-stage) and -0.018 and -0.018 (one-stage), and in B0 coverage 0 and 0, and 0 and -0.02 (one species-fit in 50). Five replicates are too few to estimate the difference itself (its 100-replicate values are -0.006 and +0.002); what they show is that the chain length does not change it.

So the conclusions do not rest on the short chains. A user fitting one dataset with `b0_prior = "community"` should run at least the package default, 5,000 + 5,000, and check `mu_b0`'s R-hat.

## Against the success criteria

From [criteria.csv](criteria.csv) and [criteria-mixing.csv](criteria-mixing.csv). Each line gives the worsening against the limit, with the paired standard error.

No harm, `base_cm_m3` against `base_m3`:

-   C8, B0 coverage no more than 0.02 lower: passes (0.001 lower, SE 0.007).
-   C9, community occupancy error no more than 0.02 larger: passes (0.006 larger, SE 0.003).
-   C10, no band more than 0.03 worse: passes (worst band 0.5 to 0.7, 0.001 worse, SE 0.003).
-   C11: collection-intercept bias changed by 0.017 (SE 0.006), coverage 0.001 higher; `p` bias changed by 0.0003, coverage unchanged; `theta0` bias changed by 0.001, coverage 0.001 lower; B coverage 0.008 higher. Passes.

No harm, `occupancy_cm_m3` against `occupancy_m3`:

-   C8: passes (0.001 higher). C9: passes (0.002 smaller). C10: passes (worst band 0.1 to 0.2, 0.006 worse, SE 0.012). C11: collection-intercept bias changed by 0.001, coverage 0.003 lower; `theta0` coverage 0.003 lower; B coverage 0.004 higher. Passes.

The same criteria at M = 2 on the same replicates fail where the community-mean run said: C9 (0.024 against 0.02) and C11 (collection-intercept bias 0.053 against 0.05) in `base_cm`, and C8 (0.026 against 0.02) in `occupancy_cm`. All three failures fall within two paired standard errors of their limits.

Mixing:

-   C12: fails in both arms (5 and 8 of 10 replicates; median ESS 38 and 56).
-   C13: fails in the two-stage arm on B0 (ESS 48 against 71, and R-hat share 0.43 against 0.25); passes on the collection intercepts there, and on everything in the one-stage arm.
-   D3, descriptive: at 5,000 + 5,000 `mu_b0` converges (R-hat at most 1.03), the averaged errors move by under 0.007 in occupancy and at most 0.05 in B0, and the community-minus-default differences by under 0.006 in occupancy.

By the rule written beforehand, `b0_prior = "community"` is harmless to ordinary communities at M = 3: C8 to C11 hold in both arms. Mixing fails at the study's chain length, but the long chains leave the averaged results and the comparison with the default prior in place (D3). That makes it, in the README's terms, a candidate default that needs longer chains to mix; the default prior's B0 is not fully converged at 1,000 + 1,000 either.

## What this does not show

-   **Only communities centred where the fixed prior is.** True B0 here is drawn from Normal(0, 1), the default prior itself, so the default is at its best and the community mean can only match or cost. This run shows that the cost is negligible at M = 3. A community centred elsewhere on the common side (most species above 70% occupancy) is not tested, though it is the rare-dominated case mirrored.
-   **One design otherwise.** 10 species, 100 sites, 3 PCR replicates on 2 primers, true collection rates of 20 to 50%. With 2 field samples per site the option's cost is small but measurable, so designs with 2 field samples per site, or with fewer species to estimate `mu_b0` from, would feel it more.
-   **No binary parent.** It has no field samples to vary; at M = 2 the community mean passed C8 and C11 in it.
-   **The squeeze.** Both priors leave species at 0.1 to 0.2 occupancy estimated 0.15 to 0.18 too common in these communities. Neither prior addresses that.
-   **Point estimates for occupancy.** As before, the occupancy rows are posterior means with no interval.

An aside, the same under both priors: the third field sample narrows the intervals of B by 12 to 13% but its coverage falls, from 0.936 to 0.912 in the two-stage default arm and from 0.925 to 0.915 in the one-stage arm. That is the B undercoverage found by the R = 300 paper study, which other runs are investigating; it does not bear on the choice of B0 prior.

## What this suggests

-   **For Alex: B0's community mean is a candidate default at a realistic 3 field samples per site.** It removes the rare-community occupancy shift (rare-M3 run) and costs ordinary communities nothing measurable (this run). At 2 field samples per site its costs are small: three no-harm criteria fail, each within two standard errors of its limit.
-   **What it needs before it could be a default is mixing.** `mu_b0` mixes slowly because it moves together with the community's mean B0, which trades against the mean collection intercept. Options are longer default chains with a convergence check on `mu_b0`, or a sampler move that updates `mu_b0` and the species' B0 together (a joint shift, or a non-centred parametrisation), which would target exactly the direction in which the chains are slow.
-   **The collection-intercept community mean can stay experimental or be dropped**, as the rare-M3 run concluded.

## Files and reproduction

-   [README.md](README.md): design and success criteria, written before the run.
-   [analyse.R](analyse.R) writes [pairing-check.csv](pairing-check.csv), [community.csv](community.csv), [occupancy-by-band.csv](occupancy-by-band.csv), [b0-by-truth.csv](b0-by-truth.csv), [detection-paired.csv](detection-paired.csv), [blocks.csv](blocks.csv), [tradeoff-by-band.csv](tradeoff-by-band.csv), [m-effect.csv](m-effect.csv), [mu-b0.csv](mu-b0.csv) and [criteria.csv](criteria.csv) from `results/simstudy-20261016-183451.rds` and the two M = 2 runs' results.
-   [mixing.R](mixing.R) writes [reproduction.csv](reproduction.csv), [mixing.csv](mixing.csv), [mixing-summary.csv](mixing-summary.csv), [identification.csv](identification.csv), [collection-rate.csv](collection-rate.csv) and [criteria-mixing.csv](criteria-mixing.csv) from the refits' draws, and [long.R](long.R) writes [long-chains.csv](long-chains.csv), [long-chains-summary.csv](long-chains-summary.csv) and [long-cost.csv](long-cost.csv). The draws are not committed; [launch.sh](launch.sh) regenerates them with the both-means run's [refit.R](../theta-community-mean-20261016/refit.R).
-   `draws.log`, `main.log`, `long.log` and the three `*-provenance.txt` files record the three steps.

Run from the repository root:

```sh
Rscript dev/simstudy/ordinary-m3-20261016/analyse.R
Rscript dev/simstudy/ordinary-m3-20261016/mixing.R
Rscript dev/simstudy/ordinary-m3-20261016/long.R
```

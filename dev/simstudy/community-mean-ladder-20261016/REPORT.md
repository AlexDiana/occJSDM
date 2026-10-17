# The field-sample ladder with B0's community mean: how many field samples remove the squeeze?

B0's community mean (`b0_prior = "community"`) removes the community-wide occupancy shift in rare-dominated communities, but at 2 and 3 field samples per site it left each species' occupancy squeezed towards the middle, with a slope of estimated on true occupancy of about 0.6. This run asks how many field samples per site remove that squeeze, in an ordinary community (true B0 drawn from Normal(0, 1)) and in a rare-dominated one (`b0_mean = -2.75`), and compares the default prior, whose ladder had never been run on a rare-dominated community. The design, the reused arms and the rules for "mostly removed" and "removed", written before any fit ran, are in [README.md](README.md).

## Answer

-   **With the community mean, both communities need the same field effort: the squeeze is mostly removed at 10 field samples per site and removed at 20.** That is also what the ordinary community needs under the default prior. The slope of estimated on true occupancy at 2, 3, 5, 10 and 20 field samples per site is 0.52, 0.62, 0.73, 0.85 and 0.94 in the rare-dominated community and 0.50, 0.60, 0.69, 0.85 and 0.95 in the ordinary one.
-   **Under the default prior, 20 field samples per site are not enough for a rare-dominated community.** Its slope is 0.27, 0.26, 0.34, 0.52 and 0.80. At 20 it fails "mostly removed" narrowly on both counts (slope 0.799, and species under 0.05 true occupancy overestimated by 0.103), each within two standard errors of its threshold. Species at a true occupancy of 0.036 are estimated at 0.34 with 2 field samples per site and still at 0.14 with 20.
-   **In the rare-dominated community, the community mean already does better at 5 field samples per site than the default prior does at 20.** The mean absolute error of a species' occupancy is 0.050 and the community average occupancy error -0.004, against 0.035 and +0.028 under the default prior at 20. Every band of true occupancy is within 0.06 of the truth at 5 samples. What holds the rule back at 5 is the slope, 0.73, because the community's true occupancies span a narrow range.
-   **In the ordinary community the community mean changes almost nothing from 5 field samples per site up.** On identical datasets the paired differences from the default prior in community occupancy error, mean absolute error and B0 coverage are at most 0.003 at 5, 10 and 20. The mean absolute error is 0.003 higher at 5 and 10, which is detectable (2.6 and 4.4 standard errors) but negligible. This extends the M = 3 finding that the option costs ordinary communities nothing that matters.
-   **B0's intervals undercover as field samples accumulate, and most in the rare-dominated community with the community mean.** B0 coverage there is 0.917, 0.944, 0.927, 0.867 and 0.821 at 2 to 20 field samples per site (ordinary community: 0.94 falling to 0.87 under both priors; rare-dominated default: 0.40 to 0.45 throughout). At 20 samples occupancy is right (community error +0.003, mean absolute error 0.018), yet B0 is overestimated by 0.38 on average and by 1.2 for the rarest species, and `mu_b0` sits 0.43 above the realised mean of the true B0. The ordinary community shows the same pattern at its extremes. See *B0 and mu_b0* below.
-   **`mu_b0` mixes better as field samples rise.** Its median bulk ESS in 2,000 draws is 27, 38, 109, 440 and 971 in the ordinary community and 21, 36, 45, 113 and 339 in the rare-dominated one. It converges by the README's rule at 10 and 20 in both communities. At 5 it fails narrowly; at 5,000 + 5,000 it converges and the averaged results stay where they were.

## Run

-   The 60 mixing refits ran on 16 October 2026 from 22:53 to 23:10 on 5 workers.
-   The main run of 900 fits started at 23:10 on 5 workers. Doug asked for 6 workers for the overnight run, so it was stopped at 23:31, just after the checkpoint held 80 fits. It resumed at 23:31 on 6 workers and finished on 17 October at 02:50, with no failures.
-   The long chains (replicates 1 to 5 of `M5_cm` and `base_rare_cm_m5` at 5,000 + 5,000) ran on 6 workers from 02:50 to 02:58. `mu_b0_check.R` launched them because `mu_b0` failed the convergence rule in those two arms.
-   Revisions: `5db9419` for the refits and the first 80 fits, `0dbb93f` for the rest. They differ only in the launcher and the README. The tree was clean at every launch. All runs used mains power with Low Power Mode off, R 4.6.1 with Accelerate BLAS, and one thread per fit.
-   All 60 refits, made on 5 workers, reproduce their stored rows exactly, including the 56 whose stored fits ran on 6 workers after the resume. So neither the resume nor the worker count changed any fit.
-   The pairing holds in all 234 checks of [pairing-check.csv](pairing-check.csv).
    -   Within each M, the two priors' truths are identical in every block and in each species' true occupancy.
    -   Across M within a seed label, B0, B, G, theta0, `p`, `q` and the true occupancies are identical.
    -   The stored collection-coefficient truths differ across M, as expected: the collection covariate is redrawn, and the intercepts are scored at its mean. The raw intercepts are identical, which `tests/testthat/test-simstudy-cm-ladder.R` checks.
-   Throughput: 5 workers fitted 4.0 fits per minute (80 in 20 minutes) and 6 workers 4.1 (820 in 200 minutes), with the CPU-to-wall ratio rising from 4.5 to 5.1. The sixth worker added about 3%.

## By M

Values are at 2, 3, 5, 10 and 20 field samples per site, in that order, over replicates 1 to 100 of each arm. Errors are estimate minus truth. The full table is [by-m.csv](by-m.csv).

Ordinary community, default prior (the ladder's `M2`, `M5`, `M10` and `M20`, with `base_m3` at 3):

-   Occupancy slope 0.44, 0.57, 0.66, 0.86, 0.95; mean absolute error 0.085, 0.073, 0.065, 0.040, 0.020; community occupancy error -0.020, -0.009, +0.000, +0.003, -0.002.
-   B0 coverage 0.942, 0.938, 0.916, 0.889, 0.873, with mean interval widths 3.11, 2.82, 2.42, 1.81, 1.40.
-   Collection-intercept bias -0.019, -0.016, -0.041, -0.027, -0.003; `p` bias +0.038, +0.023, +0.015, +0.007, +0.004.

Ordinary community, community mean (`base_cm`, `base_cm_m3`, `M5_cm`, `M10_cm`, `M20_cm`):

-   Occupancy slope 0.50, 0.60, 0.69, 0.85, 0.95; mean absolute error 0.094, 0.077, 0.068, 0.043, 0.020; community occupancy error -0.047, -0.015, -0.000, +0.002, -0.002.
-   B0 coverage 0.927, 0.937, 0.918, 0.890, 0.870, with widths 3.39, 3.02, 2.53, 1.85, 1.40.
-   Collection-intercept bias +0.052, +0.001, -0.038, -0.027, -0.002; `p` bias +0.040, +0.024, +0.015, +0.007, +0.004.
-   `mu_b0`: R-hat under 1.05 in 5, 5, 8, 9 and 10 of 10 refitted replicates; median bulk ESS 27, 38, 109, 440, 971.

Rare-dominated community, default prior (`base_rare`, `base_rare_m3`, `base_rare_m5`, `base_rare_m10`, `base_rare_m20`):

-   Occupancy slope 0.27, 0.26, 0.34, 0.52, 0.80; mean absolute error 0.170, 0.164, 0.138, 0.083, 0.035; community occupancy error +0.167, +0.158, +0.131, +0.075, +0.028.
-   52% of species-fits have true occupancy under 0.2; the share estimated under 0.2 is 0.03, 0.05, 0.12, 0.28, 0.44.
-   B0 coverage 0.422, 0.401, 0.413, 0.440, 0.454, with widths 3.51, 3.27, 2.91, 2.30, 1.73; B0 bias +1.89, +1.82, +1.65, +1.30, +0.98.
-   Collection-intercept bias -0.471, -0.445, -0.375, -0.231, -0.087; `p` bias +0.066, +0.044, +0.027, +0.015, +0.007.

Rare-dominated community, community mean (`base_rare_cm`, `base_rare_cm_m3`, `base_rare_cm_m5`, `base_rare_cm_m10`, `base_rare_cm_m20`):

-   Occupancy slope 0.52, 0.62, 0.73, 0.85, 0.94; mean absolute error 0.073, 0.059, 0.050, 0.034, 0.018; community occupancy error -0.050, -0.022, -0.004, +0.004, +0.003.
-   Share estimated under 0.2 occupancy 0.75, 0.64, 0.56, 0.53, 0.53, against 0.52 true.
-   B0 coverage 0.917, 0.944, 0.927, 0.867, 0.821, with widths 3.92, 3.57, 3.16, 2.46, 1.98; B0 bias -0.46, -0.11, +0.15, +0.34, +0.38.
-   Collection-intercept bias +0.284, +0.141, +0.035, -0.025, -0.017; `p` bias +0.081, +0.053, +0.031, +0.016, +0.007.
-   `mu_b0`: R-hat under 1.05 in 2, 6, 6, 10 and 9 of 10 refitted replicates; median bulk ESS 21, 36, 45, 113, 339.

The ordinary community's 2 and 3 field-sample arms with the community mean carry a different set of 100 truths (seed label `base`) from its 5 to 20 arms (`mladder`). Under the default prior both sets were fitted at 2 field samples per site, and the change of truths moves little: slope 0.443 against 0.448, mean absolute error 0.085 against 0.081, B0 coverage 0.942 against 0.939.

## Species occupancy by true occupancy

Mean occupancy error by band of true occupancy, at 2, 3, 5, 10 and 20 field samples per site. The full table, with clustered standard errors, is [occupancy-by-band.csv](occupancy-by-band.csv).

Rare-dominated community (species-fits per arm and mean true occupancy in brackets), default prior:

-   0 to 0.05 (87, 0.036): +0.300, +0.310, +0.292, +0.233, +0.103.
-   0.05 to 0.1 (147, 0.082): +0.262, +0.245, +0.215, +0.123, +0.049.
-   0.1 to 0.2 (310, 0.156): +0.191, +0.183, +0.142, +0.072, +0.023.
-   0.2 to 0.3 (274, 0.251): +0.122, +0.110, +0.082, +0.040, +0.012.
-   0.3 to 0.5 (181, 0.367): +0.054, +0.048, +0.044, +0.021, +0.006.

Rare-dominated community, community mean:

-   0 to 0.05: +0.030, +0.056, +0.056, +0.043, +0.021.
-   0.05 to 0.1: +0.004, +0.024, +0.031, +0.018, +0.007.
-   0.1 to 0.2: -0.033, -0.013, +0.001, +0.005, +0.004.
-   0.2 to 0.3: -0.071, -0.044, -0.020, -0.004, -0.002.
-   0.3 to 0.5: -0.131, -0.080, -0.041, -0.015, -0.006.

Under the default prior the whole rare community sits on the squeeze's high side, so every band is overestimated and the error shrinks only slowly with field samples. With the community mean the centre moves down to the community, and what is left is a squeeze around that centre: commoner species underestimated at 2 field samples per site, and the rarest species slightly overestimated until 20.

In the ordinary community the two priors agree from 5 field samples per site up. In the band 0.1 to 0.2 (18 species-fits per arm), the error at 5, 10 and 20 field samples per site is +0.142, +0.061 and -0.000 under the default prior and +0.133, +0.060 and -0.000 with the community mean. In the band 0.7 to 0.9 (72 species-fits) it is -0.088, -0.025 and -0.025, against -0.084, -0.027 and -0.023.

## The rules

From [squeeze.csv](squeeze.csv) and [squeeze-answer.csv](squeeze-answer.csv). "Mostly removed" needs a slope of at least 0.80 and every band holding at least 10 species-fits within 0.10 of zero; "removed" needs 0.90 and 0.05.

-   Ordinary, default prior: mostly removed at 10 (slope 0.86, worst band +0.061), removed at 20 (0.95, -0.025).
-   Ordinary, community mean: mostly removed at 10 (0.85, +0.060), removed at 20 (0.95, -0.023).
-   Rare-dominated, default prior: not mostly removed by 20 (0.80, +0.103; both within two standard errors, so a narrow failure).
-   Rare-dominated, community mean: mostly removed at 10 (0.85, +0.043), removed at 20 (0.94, +0.021). At 5 the worst band is +0.056 but the slope is 0.73.

At the rung each rule names, B0 coverage is 0.889 and 0.873 for the ordinary community under the default prior, 0.890 and 0.870 with the community mean, and 0.867 and 0.821 for the rare-dominated community with the community mean. In every community the slope near 1 comes with B0 intervals that are too narrow.

## The community mean against the default prior

Community mean minus default prior on identical datasets, paired by replicate, at 2, 3, 5, 10 and 20 field samples per site. The full table, with standard errors, is [prior-paired.csv](prior-paired.csv).

-   Rare-dominated: community occupancy error -0.217, -0.180, -0.135, -0.071, -0.025; mean absolute error -0.098, -0.105, -0.088, -0.049, -0.017; B0 coverage +0.50, +0.54, +0.51, +0.43, +0.37; collection-intercept error +0.75, +0.59, +0.41, +0.21, +0.07; slope +0.25, +0.36, +0.39, +0.33, +0.14. Standard errors are at most 0.005 for the occupancy errors and 0.023 for B0 coverage, so every difference is many standard errors from zero.
-   Ordinary: community occupancy error -0.024, -0.006, -0.001, -0.000, +0.000; mean absolute error +0.013, +0.005, +0.003, +0.003, +0.000; B0 coverage -0.012, -0.001, +0.002, +0.001, -0.003; collection-intercept error +0.053, +0.017, +0.003, -0.000, +0.000. B0 intervals are wider with the community mean by 0.31, 0.20, 0.10, 0.05 and 0.00.

## B0 and mu_b0

From [b0-by-truth.csv](b0-by-truth.csv) and [mu-b0.csv](mu-b0.csv).

-   With the community mean in the rare-dominated community, B0 moves from underestimated to overestimated as field samples rise: the mean error is -0.46, -0.11, +0.15, +0.34 and +0.38. At 20 field samples per site the error by true B0 is +1.22 below -4.25 (coverage 0.58), +0.63 from -4.25 to -3.25 (0.79), +0.36 from -3.25 to -2.25 (0.84) and +0.10 from -2.25 to -1.25 (0.88).
-   `mu_b0` follows the B0 estimates, not the truth. It is -3.10, -2.75, -2.50, -2.32 and -2.28 against a realised mean of the true B0 of -2.72. It covers the generating mean of -2.75 in 0.92, 0.93, 0.90, 0.82 and 0.72 of replicates, and the realised mean in 0.93, 0.98, 0.98, 0.91 and 0.91. In the ordinary community it stays within 0.03 of the realised mean from 5 field samples per site up.
-   The ordinary community has the same pattern at its extremes under both priors. At 20 field samples per site, true B0 below -1.5 is overestimated by 0.50 under the default prior and 0.48 with the community mean, and above 1.5 underestimated by 0.53 and 0.52, with coverage 0.68 to 0.74. A rare-dominated community puts most of its species at the low extreme, so its overall B0 coverage falls further.
-   Occupancy is nonetheless right at 20 field samples per site. That suggests B0 is trading off against something else in the occupancy predictor, not against the data. One candidate: the fit underestimates the spread of the spatial field and latent factors around B0. For a rare species, less spread must be offset by a higher B0 to give the same share of occupied sites, and for a common species by a lower one, which is the direction of both errors. This run did not test that.

## Mixing

From [mixing-summary.csv](mixing-summary.csv), [mu-b0-by-m.csv](mu-b0-by-m.csv) and [long-chains-summary.csv](long-chains-summary.csv). The M = 2 and M = 3 values come from the community-mean, ordinary M = 3 and rare M = 3 runs' refits at the same chain length.

-   `mu_b0` converges by the README's rule (R-hat under 1.05 in at least 9 of 10 replicates, median bulk ESS at least 100) at 10 and 20 field samples per site in both communities. It fails at 5: in the ordinary community 8 of 10 replicates are under 1.05, with a median ESS of 109; in the rare-dominated one 6 of 10, with 45.
-   B0 and the collection intercepts also mix better with more field samples. B0's median bulk ESS in the refitted community-mean arms is 85, 145 and 216 at 5, 10 and 20 in the ordinary community (40 and 48 at 2 and 3) and 58, 94 and 144 in the rare-dominated one (37 and 49); the collection intercepts' rises from 108 at 2 to 438 at 20 in the ordinary community and from 150 to 436 in the rare-dominated one.
-   At 5,000 + 5,000, replicates 1 to 5 of both 5-sample arms converge: `mu_b0` has R-hat at most 1.02, with median bulk ESS 238 (ordinary) and 143 (rare-dominated) of 10,000. Averaged over those 5 replicates, the results barely move.
    -   Community occupancy error: -0.043 to -0.041 (ordinary) and -0.008 to -0.010 (rare-dominated).
    -   Mean absolute error: 0.070 to 0.068 and 0.051 to 0.048.
    -   B0 coverage: 0.94 to 0.92 and 0.90 to 0.90.
    -   Collection-intercept error: +0.023 to +0.013 and +0.009 to +0.027.
    -   The short-chain results stand.

## Other blocks

B coverage falls with field samples in the rare-dominated community: from 0.889 to 0.785 under the default prior and from 0.920 to 0.884 with the community mean. It holds at 0.89 to 0.94 in the ordinary community. `theta0` coverage is 0.89 to 0.98 in every arm, and the collection intercepts cover 0.85 to 0.98. `p` is overestimated at few field samples in every arm, and least at 20 (+0.004 to +0.007).

## What this does not show

-   **One design.** The answer is for 100 sites, 10 species, 3 PCR replicates per field sample, a spatial field and traits, and true sample-capture rates of 0.2 to 0.5. Lower capture rates would probably need more field samples. Only the two-stage model was run.
-   **Rungs, not a curve.** Nothing between 5 and 10 or between 10 and 20 field samples per site was fitted, so "10" means "somewhere above 5 and at most 10".
-   **Mixed truths in the ordinary community-mean ladder.** Its 2 and 3 rungs carry different truths from its 5 to 20 rungs. The default prior's two sets at M = 2 differ little, which suggests the effect is small.
-   **B0's mechanism.** The explanation offered for B0's overestimate at high replication is a hypothesis. Checking it needs the fitted spread of the spatial field and latent factors, which this run did not store.
-   **Convergence elsewhere.** Only `mu_b0`, B0 and the collection intercepts were checked, on 10 replicates per arm. Long chains ran only for the two 5-sample arms, on 5 replicates each.

## What this changes

-   **Survey design for rare-dominated communities.** Under the default prior even 20 field samples per site leave rare species overestimated and the community average 0.03 too high. With the community mean, rare-dominated and ordinary communities need the same effort: about 10 field samples per site for most of the squeeze and 20 for nearly all of it. At 5 the community mean already brings every band within 0.06.
-   **The default.** With the M = 3 run, this makes `b0_prior = "community"` close to costless for ordinary communities at every M tested from 3 up (mean absolute error at most 0.005 higher), and decisively better for rare-dominated ones at every M. Short chains mix `mu_b0` poorly at 5 or fewer field samples per site. The package default of 5,000 + 5,000 mixes it at 5.
-   **B0's calibration is a separate problem.** As field samples rise, B0's intervals narrow faster than its errors at the extremes, under either prior. In a rare-dominated community with the community mean, overall B0 coverage falls to 0.82 at 20 field samples per site while occupancy is right. Users reading B0 rather than occupancy should know this. It bears on the question for Alex in `TODO.md`.

## Files and reproduction

-   [analyse.R](analyse.R) writes [pairing-check.csv](pairing-check.csv), [by-m.csv](by-m.csv), [occupancy-by-band.csv](occupancy-by-band.csv), [squeeze.csv](squeeze.csv), [squeeze-answer.csv](squeeze-answer.csv), [prior-paired.csv](prior-paired.csv), [b0-by-truth.csv](b0-by-truth.csv) and [mu-b0.csv](mu-b0.csv). It reads this run's results and replicates 1 to 100 of the reused arms in the field-sample ladder, paper study, community-mean, rare-dominated, ordinary M = 3 and rare M = 3 runs.
-   [mixing.R](mixing.R) writes [reproduction.csv](reproduction.csv), [mixing.csv](mixing.csv), [mixing-summary.csv](mixing-summary.csv), [mu-b0-by-m.csv](mu-b0-by-m.csv), [long-chains.csv](long-chains.csv) and [long-chains-summary.csv](long-chains-summary.csv) from the kept draws.
-   [mu_b0_check.R](mu_b0_check.R) decided the long step. [launch.sh](launch.sh) ran every step.
-   `results/simstudy-20261017-025051.rds` holds the run. The logs and provenance files record it.
-   The kept draws (`draws/`) and the checkpoint are not committed; `launch.sh draws` and `launch.sh long` regenerate the draws exactly.

From the repository root:

```sh
Rscript dev/simstudy/community-mean-ladder-20261016/analyse.R
Rscript dev/simstudy/community-mean-ladder-20261016/mixing.R
```

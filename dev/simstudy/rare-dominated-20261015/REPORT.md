# Rare-dominated communities: how far is rare species' occupancy overestimated?

The R = 300 paper study found species-level occupancy compressed towards the middle in the two-stage cells, but its simulated communities contain almost no rare species. The true B0 values are drawn from the fitter's own Normal(0, 1) prior, so only about 2% of species have true occupancy under 0.2, and over a whole community the upward and downward errors cancel. This run moved the mean of the true B0 draw to -2.75, giving communities with a median true occupancy of 0.19, and refitted the `base`, `occupancy` and `binary` cells at R = 100, each paired with its parent replicate by replicate. The design and the pairing are in [README.md](README.md).

## Summary

-   **Rare species are estimated at about a third of sites, whatever their true occupancy.** In `base_rare` (two-stage eDNA), species truly present at up to 5% of sites are estimated at 0.34 on average (true mean 0.036, 9.4 times too high), at 5 to 10% at 0.34 (true 0.082, 4.2 times), at 10 to 20% at 0.35 (true 0.156, 2.2 times), and at 20 to 30% at 0.37 (true 0.251, 1.5 times). Every species under 0.1 is overestimated, and over 99% of those under 0.2. The one-stage `occupancy_rare` cell agrees to within 0.02 in every band.
-   **The community's average occupancy is overestimated by 84%.** In `base_rare` it is estimated at 0.366 against a true 0.199 (error +0.167, SE 0.004); in `base` the error was -0.023 (SE 0.003). 54% of `base_rare`'s species have true occupancy under 0.2 but 3% are estimated there, and 23% are truly under 0.1 but none is estimated there. `occupancy_rare` gives 0.367 against 0.197, 87% too high.
-   **Is the shrinkage worse than in the paper study's cells? Per species, no; per community, much.** In narrow bands of true occupancy where parent and rare arm overlap, the bias agrees to within 0.02 in 11 of 12 comparisons: for example +0.171 (SE 0.009) in `base` and +0.168 (0.006) in `base_rare` for species at 15 to 20%. The bias belongs to the species' own truth, and a rare community puts most species where it is large. The paper study never showed the bottom of the curve, which is flat: below 0.2 the estimates sit near 0.35 whatever the truth, so within a rare community they barely rank the species (correlation of estimate with truth 0.36, against 0.67 in `base`; slope 0.27, against 0.45).
-   **B0 recovers about a quarter of the 2.75 shift in the eDNA cells.** Each species' posterior mean of B0 moved by -0.72 (SE 0.02) from its parent in `base_rare` and by -0.79 (0.02) in `occupancy_rare`: 26% and 29% of the true shift. B0 coverage fell from 0.939 to 0.422 and from 0.935 to 0.393, and true values below -3 are covered at most 7% of the time. In `binary`, where presence is observed, 72% of the shift is recovered, but coverage still falls from 0.860 to 0.449, and there the rare community does worsen B0 at a given true value.
-   **The fit trades occupancy against the collection rate.** The collection-side truths are identical in each pair, yet the rare arms' fits put each species' collection intercept 0.47 (SE 0.02) logit lower in `base_rare` and 0.41 (0.02) lower in `occupancy_rare`, and `theta0` about 0.017 lower. The data mainly fix how often a species is detected. The Normal(0, 1) prior penalises a B0 near -2.75 heavily, so the fit raises occupancy and lowers the chance that a sample catches the species, even though the collection intercept moves away from its own prior mean to do so.
-   **B suffers too.** B coverage falls from 0.936 to 0.889 in `base_rare`, from 0.925 to 0.876 in `occupancy_rare` and from 0.905 to 0.795 in `binary_rare`, with RMSE up 18 to 40% and intervals at most 10% wider.

## Run

The three arms ran from 23:43 on 15 October to 00:38 on 16 October 2026 at revision `b0d5811`, with a clean tree at start and at write and `HEAD` unmoved: 300 fits, no failures, 5 replicate processes, one thread per fit, 2 chains of 1,000 draws after 1,000 burn-in, on mains power with Low Power Mode off, under R 4.6.1 with Accelerate BLAS. Another study run with 5 workers shared the laptop throughout, so the 55 minutes is not a speed measurement; each replicate seeds itself, so the results do not depend on the load. The run used the worktree `.worktrees/rare-dominated`, as `launch.sh` and `provenance.txt` record; its outputs were then copied, checksums matching, to the worktree this report was committed from. The parents are replicates 1 to 100 of the R = 300 paper study at `de4389c`, whose package code is identical. Every replicate is present in both arms of every pair.

## Pairing

At all 100 replicates the rare arms' true B0 equals the parent's minus 2.75 exactly, and B, G, `beta_theta`, `theta0`, `p` and `q` equal the parent's exactly ([pairing-check.csv](pairing-check.csv)). The datasets differ, because the occupancy states do.

## How many species fall in each band

True occupancy is the share of the 100 sites a species occupies. Over replicates 1 to 100 (1,000 species per cell):

-   `base_rare`: 87 at up to 0.05, 147 at 0.05 to 0.1, 310 at 0.1 to 0.2, 274 at 0.2 to 0.3, 181 at 0.3 to 0.5, 1 above 0.5. Median 0.19.
-   `occupancy_rare`: 78, 148, 330, 282, 159 and 3. Median 0.19.
-   `base`: none under 0.1, 18 at 0.1 to 0.2, 68 at 0.2 to 0.3, 398 at 0.3 to 0.5, 440 at 0.5 to 0.7, 76 at 0.7 to 0.9. Median 0.51.
-   `occupancy`: 1 under 0.1, 19 at 0.1 to 0.2, then 60, 416, 431 and 72, and 1 above 0.9. Median 0.51.

The binary cells have no occupancy rows: presence is observed there, so the occupancy of the surveyed sites is known rather than estimated. Before fitting, the shift put 247 of `binary_rare`'s 1,000 species at 0.1 or less, against 2 in `binary` ([README.md](README.md)).

## Species-level occupancy by true band

Estimated occupancy is the posterior mean of the occupancy state averaged over sites. Bias is estimated minus true occupancy; standard errors treat each replicate as one cluster. Bands with fewer than 10 species are left out here and kept in [occupancy-by-band.csv](occupancy-by-band.csv).

-   `base_rare`. Up to 0.05: true 0.036, estimated 0.336, bias +0.300 (SE 0.008). 0.05 to 0.1: 0.082, 0.344, +0.262 (0.007). 0.1 to 0.2: 0.156, 0.347, +0.191 (0.005). 0.2 to 0.3: 0.251, 0.372, +0.122 (0.004). 0.3 to 0.5: 0.367, 0.421, +0.054 (0.005).
-   `occupancy_rare`. Up to 0.05: 0.035, 0.352, +0.317 (0.011). 0.05 to 0.1: 0.081, 0.342, +0.261 (0.009). 0.1 to 0.2: 0.156, 0.356, +0.200 (0.006). 0.2 to 0.3: 0.249, 0.369, +0.120 (0.004). 0.3 to 0.5: 0.369, 0.415, +0.046 (0.006).
-   `base`, replicates 1 to 100. 0.1 to 0.2: 0.175, 0.314, +0.139 (0.012). 0.2 to 0.3: 0.266, 0.391, +0.124 (0.010). 0.3 to 0.5: 0.422, 0.441, +0.019 (0.004). 0.5 to 0.7: 0.584, 0.519, -0.065 (0.004). 0.7 to 0.9: 0.769, 0.600, -0.169 (0.010).
-   `occupancy`, replicates 1 to 100. 0.1 to 0.2: 0.165, 0.329, +0.165 (0.023). 0.2 to 0.3: 0.261, 0.377, +0.116 (0.011). 0.3 to 0.5: 0.421, 0.447, +0.026 (0.003). 0.5 to 0.7: 0.590, 0.540, -0.050 (0.004). 0.7 to 0.9: 0.764, 0.621, -0.143 (0.009).

Of the 234 `base_rare` species at 0.1 or less, the middle 90% are estimated between 0.21 and 0.46. Over all 300 replicates the parents give the same picture with more rare species: `base` at 0.1 to 0.2 (49 species) is estimated at 0.350 against 0.170 ([occupancy-by-band-parents-r300.csv](occupancy-by-band-parents-r300.csv)).

**Like for like.** Banding by 0.05 between 0.1 and 0.4, with the parents over all 300 replicates, each line gives `base` and `base_rare`, then `occupancy` and `occupancy_rare` ([matched-truth.csv](matched-truth.csv)):

-   0.1 to 0.15: +0.209 (0.036, 12 species) and +0.213 (0.007); +0.167 (0.034, 12 species) and +0.226 (0.007).
-   0.15 to 0.2: +0.171 (0.009) and +0.168 (0.006); +0.187 (0.011) and +0.176 (0.007).
-   0.2 to 0.25: +0.135 (0.011) and +0.133 (0.006); +0.150 (0.010) and +0.132 (0.006).
-   0.25 to 0.3: +0.109 (0.006) and +0.106 (0.006); +0.092 (0.007) and +0.103 (0.007).
-   0.3 to 0.35: +0.072 (0.006) and +0.081 (0.007); +0.084 (0.006) and +0.065 (0.009).
-   0.35 to 0.4: +0.042 (0.004) and +0.038 (0.008); +0.044 (0.004) and +0.048 (0.008).

The one difference over 0.02, in `occupancy` at 0.1 to 0.15, rests on 12 parent species and is 1.7 standard errors.

## B0 by true value

Bias is the posterior mean minus the truth ([b0-by-truth.csv](b0-by-truth.csv)).

-   `base_rare`. Truth -4 to -3 (288 species): bias +2.51 (SE 0.04), coverage 0.07. -3 to -2 (385): +1.69 (0.02), coverage 0.49. -2 to -1 (188): +0.93 (0.03), coverage 0.92. Below -4 (95): +3.4 to +4.3, coverage 0. The posterior means range only from -2.5 to 0.9 while the truths run from -6.0 to 0.2.
-   `occupancy_rare`. -4 to -3 (288): +2.47 (0.03), coverage 0.06. -3 to -2 (380): +1.65 (0.03), coverage 0.44. -2 to -1 (195): +0.92 (0.03), coverage 0.87. Below -4 (94): coverage 0.
-   `binary_rare`. -4 to -3 (276): +1.02 (0.02), coverage 0.30. -3 to -2 (373): +0.62 (0.02), coverage 0.55. -2 to -1 (202): +0.40 (0.03), coverage 0.64. Below -4 (116): +1.5 to +2.1, coverage 0.06 or less.
-   The parents, for comparison: in `base` truths of -3 to -2 (19 species) carry +1.48 and coverage 0.63, and truths over 2 (28) carry -1.82 and coverage 0.29. The pull is symmetric; the parents simply have few species far from 0.

At matched true B0 the bias in the two eDNA cells agrees with their parents' to within about two standard errors. `binary` does not: at truths of -1.5 to -1 the bias is +0.23 (0.02) in `binary` and +0.38 (0.04) in `binary_rare`, and coverage 0.84 against 0.57. In the binary model the other species' rarity costs a species precision, which it does not in the eDNA models.

## The community as a whole

Per replicate, the mean over its 10 species, then the mean of those over replicates ([community.csv](community.csv)):

-   `base`: occupancy 0.505 true, 0.482 estimated, error -0.023 (0.003). Mean B0 0.03 true, -0.11 estimated. B0 coverage 0.939, slope of B0 estimate on truth 0.35.
-   `base_rare`: occupancy 0.199 true, 0.366 estimated, error +0.167 (0.004). Mean B0 -2.72 true, -0.82 estimated. B0 coverage 0.422, slope 0.14.
-   `occupancy`: occupancy 0.504 true, 0.493 estimated, error -0.011 (0.003). B0 coverage 0.935, slope 0.40.
-   `occupancy_rare`: 0.197 true, 0.367 estimated, error +0.171 (0.003). Mean B0 -2.70 true, -0.82 estimated. B0 coverage 0.393, slope 0.14.
-   `binary`: mean B0 -0.01 true, 0.00 estimated. B0 coverage 0.860, slope 0.81.
-   `binary_rare`: mean B0 -2.76 true, -1.98 estimated; on the probability scale 0.12 against 0.06. B0 coverage 0.449, slope 0.62.

A user reading a rare-dominated eDNA community off the default fit would see a community of mid-frequency species: half its species are truly under 0.2 and 3% are estimated there.

## How much of the shift the fits recover

Each rare-arm species against the same species in its parent, whose true B0 is 2.75 higher ([b0-shift-recovered.csv](b0-shift-recovered.csv)): the posterior mean moved by -0.72 (SE 0.02) in `base_rare`, -0.79 (0.02) in `occupancy_rare` and -1.98 (0.02) in `binary_rare`, that is 26%, 29% and 72% of the shift.

## The detection side

The collection-side truths are identical in each pair, so any change in their estimates is the fit's response to the rarer occupancy ([detection-paired.csv](detection-paired.csv)):

-   Collection intercept: bias -0.002 in `base` and -0.471 (SE 0.022) in `base_rare`; +0.042 in `occupancy` and -0.372 (0.017) in `occupancy_rare`. Coverage falls from 0.97 to 0.92 and from 0.97 to 0.91.
-   Collection slope: unchanged within a standard error.
-   `theta0`: 0.017 lower in both rare arms; coverage 0.98 to 0.92 in `base_rare`.
-   `p`: 0.028 higher in `base_rare` (bias +0.066 against +0.039). `q`: unchanged.

So the fit explains a rare species' few detections partly as a commoner species that samples seldom catch. With 2 samples per site the data separate the two only weakly, which is the trade-off the [occupancy-intercept prior study](../occupancy-intercept-prior/REPORT.md) described.

## Other blocks

From [blocks.csv](blocks.csv): B coverage falls from 0.936 to 0.889 (`base`), 0.925 to 0.876 (`occupancy`) and 0.905 to 0.795 (`binary`); B RMSE rises from 0.85 to 1.00, 0.81 to 0.98 and 0.43 to 0.60; mean width goes from 3.22 to 3.41, 2.92 to 3.20 and 1.37 to 1.37. G coverage holds in the eDNA cells (0.96 to 0.94) and falls in `binary` (0.94 to 0.90).

## What this does not show

-   **One shift and one design.** One B0 mean, -2.75, with 10 species, 100 sites, 2 samples per site, 3 PCR replicates and the paper study's detection ranges. More samples per site should separate occupancy from the collection rate better; that is not tested here.
-   **Point estimates for occupancy.** The occupancy rows hold the posterior mean of the occupancy state averaged over sites, with no interval, so there is no coverage statement for a species' occupancy.
-   **No convergence check.** The runner stores no R-hat. The pattern, estimates near the prior's centre with collection intercepts moved to compensate, is what a dominant prior produces, but slow mixing along the same trade-off would look similar and has not been ruled out at 2 chains of 1,000.
-   **The default prior only.** The occupancy-intercept prior study widened the fixed prior to standard deviations of 2, 3 and 5; in two-stage fits that cut the low-band error only from about 15 to 12 points and worsened mixing and mid-range bias. It is not repeated here.
-   **The binary cells' occupancy.** They have no occupancy rows; their B0 bias matters for prediction at new sites.

## What this changes

-   **For the paper.** With the default prior, the eDNA models' species-level occupancy estimates for species under about 0.2 carry little information: they sit near 0.35 whatever the truth. The paper study's B0 coverage of about 0.93 holds because its truths are drawn from the prior; in a rare-dominated community it is about 0.4. Community summaries, such as mean occupancy or the share of rare species, are badly biased in exactly the communities eDNA surveys usually find.
-   **The paper study's per-band numbers carry over.** In the eDNA cells the bias at a given true occupancy does not depend on the rest of the community, so the paper study's band results describe a rare species' bias correctly; what they understate is how many species a real community puts in those bands.
-   **The cause is the prior's centre and the trade-off, not the sampler.** A fix has to move where B0 is pulled to, or give the data more to separate occupancy from the collection rate.

## For Alex: what would estimating the community mean of B0 need?

The question: should B0's prior become `Normal(mu_b0, sigma_b0^2)` with `mu_b0` estimated, say under a Normal(0, 4) hyperprior, instead of fixed at 0? In a rare community `mu_b0` would move down and release the rare species from the pull towards 0.5. As far as this run can see, testing it would need:

-   **Every place the B0 prior mean enters, kept in agreement.** The zero mean is implicit in the joint B0/B/L/Bs draw (`sample_BBsL_cpp()` and `sample_BBsL_parallel()` in `src/jsdm.cpp`) and in the spatial range step's integrated likelihood (`spatial_range_logweights()` and its grouped variant in `R/jsdmfun.R`, whose prior-mean row for B0 is zero). If the draw and the range step used different means, the range moves would target the wrong posterior. Tier-1 tests should check the log weights against direct integration at a non-zero mean, and check that a mean fixed at 0 reproduces the current sampler bit for bit.
-   **A conjugate step for `mu_b0`** given the B0 values, tested against its closed form, and its draws carried through the results, summaries, posterior-draw extractors and convergence diagnostics.
-   **The trade-off at community scale.** This run shows the fit lowering the collection intercepts to keep the B0 penalty down. With `mu_b0` free, the community means of B0 and of the collection intercepts can trade off in the same way, anchored only by the collection intercepts' fixed Normal(0, 1) prior. Is that anchor enough, or should the collection intercepts' mean be estimated too, or deliberately kept fixed?
-   **`sigma_b0`.** Fixed at 1 with a free mean, or estimated as well? The occupancy-intercept prior study found a wider fixed `sigma_b0` worsened mixing in two-stage fits.
-   **Recovery and calibration.** The six paired cells here, three parents and three rare arms, can be rerun with the change at no harness cost. Success would mean low-band occupancy bias well below the current +0.19 to +0.30, B0 coverage near nominal across its true range, unbiased collection intercepts, and no loss in the parents. Mixing of `mu_b0` and B0 needs checking at the runner's 2 chains of 1,000, with longer chains if R-hat is high. A mixed community, a few common species among many rare ones, would check that one community mean does not drag the common species down.
-   **Defaults.** Whether it becomes the default or an option; as the default it would change the lessons' outputs and `sampleresults`.

Is estimating the community mean of B0 the change you would make, and if so, would you estimate `sigma_b0` and the collection intercepts' community mean as well?

## Files and reproduction

-   [analyse.R](analyse.R) writes every CSV here from this run's results file and the R = 300 paper study's: [occupancy-by-band.csv](occupancy-by-band.csv), [occupancy-by-band-parents-r300.csv](occupancy-by-band-parents-r300.csv), [matched-truth.csv](matched-truth.csv), [b0-by-truth.csv](b0-by-truth.csv), [community.csv](community.csv), [b0-shift-recovered.csv](b0-shift-recovered.csv), [detection-paired.csv](detection-paired.csv), [blocks.csv](blocks.csv) and [pairing-check.csv](pairing-check.csv).
-   `results/simstudy-20261016-003851.rds` holds the three arms; `run.log` and `provenance.txt` record the run.
-   The harness setting is `b0_mean` in `tests/testthat/helper-simstudy.R`, tested in `tests/testthat/test-simstudy-b0-shift.R`.

Run from the repository root:

```sh
Rscript dev/simstudy/rare-dominated-20261015/analyse.R
```

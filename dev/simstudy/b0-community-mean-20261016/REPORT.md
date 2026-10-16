# A community mean for the occupancy intercepts: does it remove the rare-community shift?

The rare-dominated run ([its report](../rare-dominated-20261015/REPORT.md)) found that in eDNA fits of a community of mostly rare species, every rare species is estimated at about a third of sites, because B0's fixed Normal(0, 1) prior, centred on 50% occupancy, decides the weakly identified split between how often a species is present and how often a sample catches it. This run tests the change that report put to Alex: an experimental, off-by-default option, `listPriors = list(b0_prior = "community")`, that centres B0's prior on an estimated community mean `mu_b0` instead of 0. Six arms fitted with it at R = 100 are compared with default-prior fits of identical datasets. The design and the success criteria, written before any fit ran, are in [README.md](README.md).

## Summary

-   **The option fixes most of the occupancy and B0 shift in rare-dominated eDNA communities.** The community's mean occupancy error falls from +0.167 to -0.050 (SE 0.004) in `base_rare_cm` (two-stage) and from +0.171 to -0.013 (SE 0.003) in `occupancy_rare_cm` (one-stage). Species truly at up to 5% of sites are now estimated at 0.07 and 0.09, against 0.34 and 0.35 before. B0 coverage rises from 0.42 to 0.92 and from 0.39 to 0.94, and `mu_b0`'s 95% interval contains the generating mean, -2.75, in 92% and 97% of replicates.
-   **It does not meet the success definition written beforehand.** The collection intercepts' error does not shrink: it reverses, from -0.47 to +0.28 (SE 0.02) in `base_rare_cm` and from -0.37 to +0.22 (SE 0.02) in `occupancy_rare_cm` (C5 fails in both). In the two-stage arm the fit now overshoots, putting `mu_b0` at -3.10 and the community's occupancy 25% too low. The binary rare arm improves but falls short of its targets: B0 coverage 0.45 to 0.83, community B0 error +0.78 to +0.35 (C7 fails).
-   **The cause has moved to the collection intercepts' own prior.** With B0's centre free, the only thing deciding the occupancy-collection split at community scale is the collection intercepts' fixed Normal(0, 1) prior. 98 to 99% of the true collection intercepts lie below its centre of 0, and it now pulls the split the other way: collection up, occupancy down. Across replicates the community's mean collection-intercept error and `mu_b0`'s error are negatively correlated (-0.50 in `base_rare_cm`, -0.58 in `occupancy_rare_cm`).
-   **The parents lose a little.** In `base_cm` the community occupancy error goes from -0.023 to -0.047 and `mu_b0` settles at -0.29 against a true 0 (C9 and C11 fail narrowly); in `occupancy_cm` B0 coverage falls from 0.935 to 0.909 (C8 fails). The binary parent is unchanged.
-   **`mu_b0` mixes slowly in the eDNA models.** At the study's 2 chains of 1,000 + 1,000, its median bulk ESS is 21 to 29 of 2,000 draws and its R-hat is under 1.05 in only 2 to 5 of 10 replicates per eDNA arm (C12 fails in all four), and B0 mixes somewhat worse than under the default prior (C13 fails on R-hat in all four). In the binary arms `mu_b0` mixes well. Refits at 5,000 + 5,000 bring `mu_b0`'s R-hat to 1.07 or less, and averaged over replicates they leave the community occupancy, B0, `mu_b0` and collection-intercept errors within 0.05 of their 1,000 + 1,000 values, so the conclusions do not rest on the short chains.
-   **The species within a community are still compressed.** In `base_rare_cm`, averaged by true value, the B0 estimates run only from -3.7 for truths below -5 to -2.2 for truths above 0 (slope of estimate on truth 0.34, against 0.14 with the default prior), and species truly at 30 to 50% of sites are now underestimated by 0.13. The community mean moves the whole community to the right place; it does not spread the species out. That is ordinary partial pooling, not a defect of the option: `sigma_b0` = 1 is the simulator's own spread of B0, and each species' data say little, so its posterior mean regresses most of the way to the community mean. The parents show the same slope (0.35 in `base`).

In short, estimating the community mean of B0 removes the block shift in occupancy and B0 but leaves the occupancy-collection split decided by a prior, now the collection intercepts' prior, so it needs a companion change before it could become a default. The natural next step is to estimate the collection intercepts' community mean as well.

## The model change

Each species' intercept B0 gets a `Normal(mu_b0, sigma_b0^2)` prior, with `sigma_b0` still fixed at 1 by default, and `mu_b0` gets a `Normal(0, mu_b0_sd^2)` hyperprior with `mu_b0_sd` defaulting to 2. After the coefficient block is drawn, each iteration draws `mu_b0` from its conjugate normal conditional given the species' B0 draws. The next iteration's spatial range step and coefficient draw then use it as B0's prior mean. Those are the five places that encode B0's prior: both C++ coefficient samplers (`sample_BBsL_cpp()`, `sample_BBsL_parallel()`), the R reference sampler (`sample_BBsL()`), and both range log-weight functions (`spatial_range_logweights()` and its grouped variant), which integrate B0 out before the range is chosen. Keeping all five in agreement matters because commit `f1d2e12` showed that steps disagreeing about a prior break the sampler.

The draws are returned in `results_output$jsdm_output$mu_b0_output`, a niter by nchain array present only when the option is on. `returnConvergenceDiagnostics()` reports them, and `infos$intercept_prior` records the applied prior. The default, `b0_prior = "fixed"`, passes a zero mean, which is the arithmetic and random stream the samplers had before; it adds no output element and leaves `infos$intercept_prior` as it was. The option is documented in `?runOccJSDM` as experimental and not a recommended setting.

## Tests

`tests/testthat/test-b0-community-prior.R`, all passing:

-   The option's parser keeps the default prior list unchanged and rejects malformed values.
-   The conjugate `mu_b0` draw matches its closed form, to rounding error from a fixed seed and in mean and variance over 20,000 draws.
-   In all three coefficient samplers, B0's conditional draws centre on `mu_b0` at three values of `sigma_b0` and `mu_b0`, with the binary and the continuous likelihood. A sampler ignoring `mu_b0` could not pass, because the designs separate means by more than 10 Monte Carlo standard errors. With no data, only B0's draws move with `mu_b0`; B, L and Bs keep their means.
-   Both range log-weight functions, and the public function's route to the grouped one, match an independent observation-space integration at non-zero means, and a one-site binary case matches numerical integration.
-   Within a fit, both sampling sites receive the `mu_b0` drawn in the iteration before, verified against the saved draws, and receive 0 throughout by default.
-   Default fits are identical with and without `b0_prior = "fixed"`, spatial and not. Community fits have the expected output shape, and `thinOutput()` and `returnConvergenceDiagnostics()` handle `mu_b0`.
-   A 30-species binary simulation with true intercepts around -2 recovers the community mean: posterior mean -1.75 (95% interval -2.13 to -1.36) against a realised -1.82, with the species' average intercept error 0.05 against 0.22 under the fixed prior.

The harness setting is tested in `tests/testthat/test-simstudy-b0-community.R`, and the existing intercept-prior, parallel-chains, RNG-safety, diagnostics, posterior-draws and spatial suites pass unchanged.

**Default fits are bit-identical to before, checked twice.** [check_default_identity.R](check_default_identity.R) fits seven seeded models under the branch base (`34eac47`) and under this branch: two-stage, occupancy, binary and continuous, spatial and not, and one with a non-default `sigma_b0`. All seven whole fit objects are identical, apart from `infos$execution`, which records timings and process ids ([default-identity.csv](default-identity.csv)). And at study scale, refitting replicates 1 to 10 of the six default-prior arms on this branch reproduced, element for element, the rows that the paper study (at `de4389c`) and the rare-dominated run (at `b0d5811`) had stored ([reproduction.csv](reproduction.csv)).

## Run

The six community arms ran from 08:48 to 11:06 on 16 October 2026 at revision `62b3138`, with a clean tree at start and at write and `HEAD` unmoved: 600 fits, no failures, 5 replicate processes, one thread per fit, 2 chains of 1,000 draws after 1,000 burn-in, under R 4.6.1 with Accelerate BLAS. The laptop slept for about 50 minutes mid-run and was on battery for part of it, and until 10:07 the two refits below ran alongside it on 5 more processes each. So the 138 minutes is not a speed measurement. Each replicate seeds itself, so none of this affects the results.

Every arm is paired with a default-prior arm on identical datasets: all truths and every species' true occupancy agree at all 100 replicates in all six pairs ([pairing-check.csv](pairing-check.csv)). The default-prior arms are replicates 1 to 100 of `base`, `occupancy` and `binary` from the R = 300 paper study and the three arms of the rare-dominated run, whose package code is identical to this branch's base.

Two refits kept the draws that the runner discards:

-   [refit_draws.R](refit_draws.R) refitted replicates 1 to 10 of all twelve arms at the study's settings (08:48 to 10:05, 120 fits, no failures). Every refit reproduces its arm's stored rows exactly, the 60 community refits included, so the mixing numbers describe the study's own fits.
-   [refit_long.R](refit_long.R) refitted replicates 1 to 5 of the four eDNA community arms and the two default rare arms at 5,000 + 5,000 (30 fits, no failures). It was added after the 1,000 + 1,000 refits showed slow mixing, so it is not one of the success criteria. From 10:07 to 11:06 its processes were paused (SIGSTOP, then SIGCONT once the main run finished) to keep the laptop at 5 workers in total. The pause line the launcher appended to `long.log` was overwritten by the R process writing the same file.

## Rare-dominated eDNA communities

Community-level numbers average each replicate's 10 species first ([community.csv](community.csv)); standard errors treat replicates as independent. For each arm the default-prior value comes first and the community-mean value second:

-   `base_rare_cm`. Mean occupancy 0.199 true; estimated 0.366 and 0.149, error +0.167 and -0.050 (SE 0.004). Mean B0 -2.72 true; estimated -0.82 and -3.18. `mu_b0` -3.10 (coverage of -2.75: 0.92). B0 coverage 0.422 and 0.917. Share of species estimated under 0.1: 0.000 and 0.336, against 0.234 truly; under 0.2: 0.029 and 0.753, against 0.544.
-   `occupancy_rare_cm`. Mean occupancy 0.197 true; estimated 0.367 and 0.184, error +0.171 and -0.013 (SE 0.003). Mean B0 -2.70 true; estimated -0.82 and -2.71. `mu_b0` -2.64 (coverage 0.97). B0 coverage 0.393 and 0.940. Under 0.1: 0.002 and 0.184, against 0.226; under 0.2: 0.032 and 0.601, against 0.556.

By true occupancy band, paired species by species ([occupancy-by-band.csv](occupancy-by-band.csv)), bias with the default prior and then with the community mean:

-   `base_rare_cm`. Up to 0.05: +0.300 and +0.030. 0.05 to 0.1: +0.262 and +0.004. 0.1 to 0.2: +0.191 and -0.033. 0.2 to 0.3: +0.122 and -0.071. 0.3 to 0.5: +0.054 and -0.131.
-   `occupancy_rare_cm`. Up to 0.05: +0.317 and +0.059. 0.05 to 0.1: +0.261 and +0.040. 0.1 to 0.2: +0.200 and +0.000. 0.2 to 0.3: +0.120 and -0.031. 0.3 to 0.5: +0.046 and -0.090.

The bias no longer depends on rarity the way it did. With the default prior every species was pulled towards 0.35; with the community mean the rare species are close to right and the community's commonest species are pulled down towards the community's centre. The slope of estimated on true occupancy rises from 0.27 to 0.52 and from 0.22 to 0.56.

B0 by its true value ([b0-by-truth.csv](b0-by-truth.csv)): in `base_rare_cm`, truths of -4 to -3 now carry a bias of -0.05 with coverage 0.99 (default +2.51, 0.07), but truths of -1 to 0 carry -1.85 with coverage 0.44 (default +0.17, 0.98). `occupancy_rare_cm` is similar, centred higher: -4 to -3 at +0.45 (coverage 0.96), -1 to 0 at -1.13 (0.79). This is the within-community compression of the summary, and it is what a posterior mean does when each species' data are weak. The 95% interval for a species' B0 is 3.92 wide on average in `base_rare_cm`, against about 4.67 for the community prior alone (`sigma_b0` = 1 widened by `mu_b0`'s own uncertainty), so the data remove only about 30% of the prior variance and the estimate moves only about a third of the way from `mu_b0` towards the truth. The slope of estimate on truth matches: 0.34. The same calculation gives 0.36 for `occupancy_rare_cm` (slope 0.37) and 0.84 for `binary_rare_cm` (slope 0.76), where presence is observed. In the eDNA models a species detected at a few sites could be rare and easily caught or commoner and seldom caught, so its own data barely pin its B0. The compression is not a sign of overestimating rare species. Truths above the community mean are pulled down as much as those below are pulled up, and on the occupancy scale the pull on the rarest species hardly shows, because occupancy changes little with B0 down there: B0 truths of -5 to -4 are still overestimated by 0.71, yet species at up to 5% of sites are only 0.03 too high. And the centre itself sits 0.35 too low in this arm, which offsets the pull at the bottom and adds to it at the top.

## The occupancy and collection trade-off

The collection-side truths are identical in each pair, so their change measures the fit's response ([detection-paired.csv](detection-paired.csv)):

-   Collection intercept: bias -0.471 and +0.284 in `base_rare_cm` (paired change +0.75, SE 0.02); -0.372 and +0.216 in `occupancy_rare_cm` (+0.59, SE 0.02). Coverage rises (0.92 to 0.98, 0.91 to 0.97) because the intervals widen by 25%.
-   Collection slope: unbiased in both, but intervals 70% wider (2.16 to 3.65 in `base_rare_cm`), because the fit now thinks fewer sites are occupied.
-   `theta0`: bias -0.023 and -0.005 in `base_rare_cm`, -0.018 and 0.000 in `occupancy_rare_cm`; coverage 0.92 to 0.93 and 0.95 to 0.96.
-   `p` (`base_rare_cm` only): bias +0.066 and +0.081, coverage 0.888 and 0.851.

Species by species ([tradeoff-by-band.csv](tradeoff-by-band.csv)), the default prior's collection-intercept error grew with rarity (-0.97 for species under 0.1 in `base_rare`), as each species traded its own detections between the two parts. With the community mean the error is +0.15 to +0.36 in every band of `base_rare_cm`, a community-wide tilt rather than a species-by-species one. The true collection intercepts average -0.64, and 98 to 99% lie below 0, so the collection intercepts' Normal(0, 1) prior pulls every species' collection rate up. With B0's centre free, nothing pulls back, and the fit moves the community's occupancy down to match.

This is the possibility the rare-dominated report raised: with `mu_b0` free, the community means of B0 and of the collection intercepts can trade off, anchored only by the collection intercepts' fixed prior. The anchor holds the split, but in the wrong place for these communities, whose collection rates are below 50%.

## The binary rare arm

Presence is observed in `binary_rare_cm`, so there is no collection stage. The community mean helps but not fully: community B0 error +0.78 to +0.35 (SE 0.02), B0 coverage 0.449 to 0.830, and `mu_b0` at -2.35 with coverage 0.76. Species truly below -4 are still overestimated by 0.8 to 1.3. Most of the remaining error belongs to the rarest species: their few presences cap how common they can be but say little about how rare they are, so their B0 sits near the community mean. B coverage rises from 0.795 to 0.897.

## The parents

Against `base`, `occupancy` and `binary` on identical datasets ([community.csv](community.csv), [blocks.csv](blocks.csv)):

-   `base_cm`: B0 coverage 0.939 and 0.927 (SE 0.009 for the difference); community B0 error -0.14 and -0.32; `mu_b0` -0.29 against 0 (coverage 0.88); community occupancy error -0.023 and -0.047; collection-intercept bias -0.002 and +0.052. The same tilt as in the rare arm, smaller, because common species identify their own collection rates better.
-   `occupancy_cm`: B0 coverage 0.935 and 0.909 (difference -0.026, SE 0.009); community occupancy error -0.011 and -0.015; `mu_b0` -0.06 (coverage 0.94).
-   `binary_cm`: unchanged to within a standard error in every block; `mu_b0` -0.003 (coverage 0.97).
-   Occupancy by band changes by at most 0.054 in any band with 10 or more species. The bands under 0.3 move towards the truth, by 0.04 to 0.05 in `base_cm` and 0.02 to 0.04 in `occupancy_cm`, and in `base_cm` the bands above 0.5 move away from it by about 0.02.
-   B, G, `p` and `theta0` coverage stay within about 0.01 of the default; B coverage moves from 0.936 to 0.930, 0.925 to 0.919 and 0.905 to 0.912.

In the rare arms B coverage rises: 0.889 to 0.920, 0.876 to 0.917 and 0.795 to 0.897.

## Mixing

From the refits of replicates 1 to 10 at the study's settings ([mixing-summary.csv](mixing-summary.csv), [mixing.csv](mixing.csv)), rank-normalised R-hat and bulk ESS out of 2,000 draws:

-   `mu_b0`, eDNA arms: median R-hat 1.06 to 1.08, maximum 1.20 to 1.37; median bulk ESS 21 to 29, minimum 5 to 10. R-hat under 1.05 in 5, 4, 2 and 3 of 10 replicates of `base_cm`, `occupancy_cm`, `base_rare_cm` and `occupancy_rare_cm`.
-   `mu_b0`, binary arms: R-hat at most 1.02; median bulk ESS 1,474 and 587.
-   B0, share of elements with R-hat over 1.05, default and community: `base` 0.30 and 0.50, `occupancy` 0.31 and 0.45, `base_rare` 0.42 and 0.50, `occupancy_rare` 0.34 and 0.45, `binary` 0.07 and 0.03, `binary_rare` 0.00 and 0.01. Median bulk ESS: 61 and 40, 62 and 52, 44 and 37, 55 and 51, 454 and 526, 413 and 267.
-   Collection intercepts mix better in the rare arms with the community mean (median bulk ESS 70 to 150 and 80 to 163), because B0's prior no longer opposes them.

The default-prior eDNA fits do not mix well either: B0 has R-hat over 1.05 in 30 to 42% of elements at 1,000 + 1,000. The community mean adds a slow direction, the community-wide trade-off, on top.

**Longer chains** ([long-chains-summary.csv](long-chains-summary.csv), [long-chains.csv](long-chains.csv)), replicates 1 to 5 at 5,000 + 5,000 against the same datasets at 1,000 + 1,000:

-   `mu_b0`'s R-hat falls to 1.01 to 1.03 on average (maximum 1.07), and its bulk ESS rises to 100 to 132 of 10,000 draws on average (minimum 31). It converges, but about 60 times more slowly per draw than in the binary model.
-   Averaged over the five replicates, community occupancy error changes by under 0.01 and community B0 error by under 0.05 in every arm. For single replicates, occupancy moves by at most 0.03 (0.05 in one `occupancy_cm` replicate), while `mu_b0` and the mean of B0 move by up to 0.47.
-   The default rare arms keep their bias at 5,000 + 5,000 (community B0 error +1.72 and +1.73 in `base_rare`, +1.82 and +1.87 in `occupancy_rare`). This answers the rare-dominated report's open question: the default prior's shift is not slow mixing.

So the study's averages stand, but a user fitting one dataset with this option needs chains of at least the package default, 5,000 + 5,000, and should check `mu_b0`'s convergence.

## Against the success criteria

From [criteria.csv](criteria.csv) and [criteria-mixing.csv](criteria-mixing.csv):

-   C1, community occupancy error within ±0.05: `base_rare_cm` -0.0503, fails at the boundary; `occupancy_rare_cm` -0.013, passes.
-   C2, community B0 error within ±0.5: passes in both (-0.46, -0.01).
-   C3, low-band occupancy bias under +0.10: passes in all six bands.
-   C4, B0 coverage at least 0.85: passes (0.917, 0.940).
-   C5, collection-intercept bias within ±0.15: fails in both (+0.28, +0.22), and is not halved.
-   C6, `mu_b0` covers -2.75 in at least 85%: passes in the eDNA arms (0.92, 0.97), fails in `binary_rare_cm` (0.76).
-   C7, binary rare B0 coverage at least 0.85 and community error within ±0.3: fails (0.830, +0.35).
-   C8, parent B0 coverage no more than 0.02 lower: passes in `base_cm` (-0.012) and `binary_cm`, fails in `occupancy_cm` (-0.026).
-   C9, parent community occupancy error no more than 0.02 larger: fails in `base_cm` (0.047 against 0.023), passes in `occupancy_cm`.
-   C10, no occupancy band more than 0.03 worse: passes.
-   C11, detection side and B: fails only on `base_cm`'s collection-intercept bias, which moves by 0.053 against a limit of 0.05; every coverage passes.
-   C12, `mu_b0` mixing: fails in all four eDNA arms, passes in both binary arms.
-   C13, B0 mixing relative to the default: fails on R-hat in all four eDNA arms, and on ESS in `base_cm` and `binary_rare_cm`.

By the rule written beforehand, the shift is removed only if C1, C2 and C5 hold in both rare eDNA arms, and reduced only if each of those errors is at least halved. Neither holds, because of C5: the occupancy and B0 errors fall by 70 to 99%, but the collection-intercept error reverses rather than shrinking. The parents show small losses (C8, C9 and C11 each fail in one arm by a small margin), and mixing fails at the study's chain length.

## What this does not show

-   **One design.** 10 species, 100 sites, 2 samples per site, 3 PCR replicates, true collection rates of 20 to 50%, and one B0 shift. With collection rates around 50% the collection intercepts' prior would pull less; with more samples per site the split is better identified.
-   **No mixed community.** A few common species among many rare ones would show how far one community mean drags the common species down; the 30 to 50% band in `base_rare_cm` (-0.13) suggests it would.
-   **`sigma_b0` stays fixed at 1.** That is the simulator's own spread of B0, so estimating it would not be expected to loosen the within-community compression, which comes from how little each species' data say. More field samples per site are what add that information: the field-sample ladder (branch `claude/field-sample-ladder`) found that they undo the default prior's squeeze and PCR replicates do not. With real communities, whose spread is unknown, estimating `sigma_b0` is a separate question.
-   **No continuous model.** The option applies there too but is not tested.
-   **Point estimates for occupancy.** As in the rare-dominated run, the occupancy rows are posterior means with no interval.

## What this suggests

-   **For Alex: the community mean is the right direction but not a complete fix.** It moves the community to the right level in the one-stage model and to within about 25% in the two-stage model, where the default prior's error was +84%. What remains is the same problem one level up. The occupancy-collection split at community scale is still set by a fixed prior, now the collection intercepts' Normal(0, 1).
-   **The natural companion is a community mean for the collection intercepts.** If the collection intercepts' prior were centred on an estimated mean too, neither side's fixed centre would decide the split, and the split would be left to whatever the data say, with the two hyperpriors as the only anchors. Whether the data then identify it at M = 2 is the question. It might need an informative hyperprior on the collection side, which users can reason about (how often does a sample catch a species where it is present?) more easily than about occupancy.
-   **Mixing needs attention before any default.** The community-wide trade-off is a slow direction for the Gibbs sampler: `mu_b0` and the collection intercepts move together only through many small steps. A joint move along that direction, or a non-centred parameterisation of B0, could be worth trying.
-   **Not a default yet.** The option stays experimental and off; default fits are unchanged.

Should the next test estimate the collection intercepts' community mean alongside `mu_b0`, and if so, under what hyperprior?

## Files and reproduction

-   [README.md](README.md): design and success criteria, written before the run.
-   [analyse.R](analyse.R) writes [pairing-check.csv](pairing-check.csv), [community.csv](community.csv), [occupancy-by-band.csv](occupancy-by-band.csv), [b0-by-truth.csv](b0-by-truth.csv), [detection-paired.csv](detection-paired.csv), [blocks.csv](blocks.csv), [tradeoff-by-band.csv](tradeoff-by-band.csv) and [criteria.csv](criteria.csv) from `results/simstudy-20261016-110650.rds` and the two comparison runs' results.
-   [mixing.R](mixing.R) writes [reproduction.csv](reproduction.csv), [mixing.csv](mixing.csv), [mixing-summary.csv](mixing-summary.csv) and [criteria-mixing.csv](criteria-mixing.csv) from the refits' draws, and [long.R](long.R) writes [long-chains.csv](long-chains.csv) and [long-chains-summary.csv](long-chains-summary.csv). The draws themselves (142 MB) are not committed; [refit_draws.R](refit_draws.R) and [refit_long.R](refit_long.R) regenerate them exactly, as [reproduction.csv](reproduction.csv) shows for the first.
-   [check_default_identity.R](check_default_identity.R) writes [default-identity.csv](default-identity.csv) from an export of the branch base.
-   `run.log`, `provenance.txt`, `draws.log`, `draws-provenance.txt`, `long.log` and `long-provenance.txt` record the three runs; [launch.sh](launch.sh), [launch-draws.sh](launch-draws.sh) and [launch-long.sh](launch-long.sh) launched them.

Run from the repository root:

```sh
Rscript dev/simstudy/b0-community-mean-20261016/analyse.R
Rscript dev/simstudy/b0-community-mean-20261016/mixing.R
Rscript dev/simstudy/b0-community-mean-20261016/long.R
```

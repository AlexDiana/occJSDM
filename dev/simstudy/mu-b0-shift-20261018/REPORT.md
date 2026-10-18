# The collapsed shift of the community's occupancy intercepts: report, 18 October 2026

**Verdict: the step removes the failure the check found, and the rule as it stands passes for the ordinary community and on its median reading for the rare one; one rare fit stays under the strict ESS bar.** With the shift step on, `mu_b0` has R-hat under 1.05 in all 10 fits at 100 species and 5,000 + 5,000 (7 of 10 without it) and bulk ESS of at least 100 in 9 of 10 (7 of 10): the three fits that failed before, ordinary replicate 2 and rare replicates 3 and 4, now have R-hat 1.008, 1.003 and 1.003 and bulk ESS 273, 379 and 383 (1.060, 1.128 and 1.189; 32, 13 and 8). The one fit under the bar is rare replicate 5, R-hat 1.037 and bulk ESS 64 (1.025 and 120 before), a slow wander of the mean intercept within one chain that is not the common-shift direction. In the diagnostic's three fits the shift between the chains is gone: the chain difference in the mean B0 is a tenth of a posterior sd where it was one. The posterior itself is unchanged where the earlier chains had mixed, the step accepts at its 0.44 target with an adapted proposal sd of 0.10 to 0.16, and it costs about 7% of an iteration at 100 species.

## The step

A collapsed Metropolis shift of `mu_b0` and every species' occupancy intercept B0 with the latent occupancy states z integrated out, as the drift diagnostic recommended ([mu-b0-drift-20261018](../mu-b0-drift-20261018/REPORT.md)). It is in `R/b0_shift.R`, runs in `run_one_chain()` immediately before the z draw, and is on only with `b0_prior = "community"` in an occupancy or two-stage fit. Each iteration proposes, 5 times by default, a common delta from Normal(0, sd^2) for every B0 and `mu_b0`, and accepts it under the product over sites and species of (psi' L1 + (1 - psi') L0) over (psi L1 + (1 - psi) L0) times the hyperprior ratio of `mu_b0`, where psi' = logistic(eta + delta) and L1 and L0 are the probabilities of the site's field-sample outcomes given z = 1 and z = 0, the quantities the z draw forms. The proposal sd starts at 2 / sqrt(S) and is adapted during burn-in towards an acceptance rate of 0.44, or is set by `b0_shift_sd`; `infos$b0_shift` records each chain's final sd and acceptance rates. The settings and the tests are described in [README.md](README.md).

## What ran

Both steps ran on 18 October 2026 from the clean tree at 04afeef, R 4.6.1 with Accelerate's `libBLAS.dylib`, one RcppParallel and one BLAS thread per fit (`VECLIB_MAXIMUM_THREADS=1`), mains power with Low Power Mode off (`*-provenance.txt`), the shift step at its defaults and nothing else changed from the earlier studies: the same datasets and the same chain seeds, so the chains differ from the earlier ones only because the step draws random numbers. No fit failed.

-   `drift`: the diagnostic's three fits, 14:34 to 15:01 on 3 workers, 26 to 28 minutes per fit at 5,000 + 500 x 10 with the latent states kept (the diagnostic's own run took 20 to 21 minutes per fit; see Cost).
-   `main`: the 100-species check's 10 fits, 15:02 to 16:05 on 5 workers, 34.7 to 34.9 minutes per ordinary fit and 27.6 to 28.1 per rare fit at 5,000 + 5,000 (the check's own run: 33.9 to 34.3 and 27.0 to 27.5).

Throughout both steps an App Store background process (`appstoreagent`) held about 85 to 95% of a core, so the laptop ran 4 and then 6 busy processes on its 4 performance cores; the elapsed times below carry that.

## The diagnostic's three fits with the step on

The shift between the chains is gone (`drift.csv`, `traces-drift.png`). Chain 2 minus chain 1 of the mean B0, with the diagnostic's value for the same chains without the step in brackets: -0.016 (-0.025) in ordinary replicate 2, -0.018 (+0.155) in rare replicate 3, +0.018 (+0.181) in rare replicate 4; as fractions of the posterior sd, -0.15 (-0.22), -0.13 (+0.97) and +0.13 (+1.11). The same for `mu_b0`: -0.022 (-0.026), -0.018 (+0.165) and +0.020 (+0.184). Mean psi and the occupied count follow: in rare replicate 3 the chains now differ by 0.6 occupied pairs of 10,000 (43 before), in rare replicate 4 by 19 (55), in ordinary replicate 2 by 24 (46).

Mixing of the 1,000 stored draws (rank-normalised R-hat, then bulk ESS; the diagnostic's values in brackets):

-   `mu_b0`: ordinary replicate 2, 1.015 and 292 (1.056 and 39); rare replicate 3, 1.000 and 317 (1.145 and 11); rare replicate 4, 1.002 and 318 (1.183 and 8). Tail ESS 754, 757 and 596.
-   Mean B0: 1.01 and 160 (1.10 and 19); 1.00 and 262 (1.19 and 10); 1.00 and 251 (1.27 and 6).
-   Mean psi: 1.01 and 145 (1.10 and 18); 1.00 and 336 (1.05 and 82); 1.02 and 258 (1.05 and 44). The occupied count the same to within a few.
-   The mean collection intercept: 1.00 and 281 (1.10 and 18); 1.00 and 534 (1.01 and 265); 1.02 and 137 (1.00 and 279).

Species by species (`drift-species.csv`), the sd across species of the chain difference in B0 is 0.12, 0.10 and 0.09 (0.12, 0.11 and 0.13 before) with its mean at -0.016, -0.018 and +0.018: the common shift is gone and what remains is each species' own Monte Carlo noise, as it should be with 500 stored draws per chain and a per-species posterior sd of 0.7 to 0.8.

The step itself: acceptance over the kept iterations 0.42 to 0.46 in every chain (target 0.44; 0.44 over burn-in), and the adapted proposal sd 0.107 to 0.132, next to `mu_b0`'s posterior sd of 0.146 to 0.168 in these fits. That is the scale the diagnostic anticipated (0.1 to 0.2 at 100 species).

## The 100-species check with the step on

`mu_b0`'s R-hat, then bulk ESS and tail ESS of 10,000 draws, replicates 1 to 5 (`fits.csv`; the check's values without the step are in `before-after.csv` and in its [report](../b0-community-s100-20261018/REPORT.md)):

-   Ordinary: R-hat 1.003, 1.008, 1.002, 1.002, 1.006 (before 1.006, 1.060, 1.011, 1.009, 1.006); bulk ESS 408, 273, 395, 687, 428 (217, 32, 187, 257, 165); tail ESS 928, 1207, 1111, 2378, 1268 (449, 406, 1021, 1078, 441).
-   Rare: R-hat 1.005, 1.003, 1.003, 1.003, 1.037 (1.003, 1.008, 1.128, 1.189, 1.025); bulk ESS 201, 312, 379, 383, 64 (146, 149, 13, 8, 120); tail ESS 383, 996, 1105, 1348, 207 (440, 373, 183, 62, 279).

By the rule (`rule.csv`): the ordinary arm has R-hat under 1.05 and bulk ESS of at least 100 in 5 of 5 replicates and passes both readings (median ESS 409). The rare arm has R-hat under 1.05 in 5 of 5 and bulk ESS of at least 100 in 4 of 5, so it passes the median reading (median ESS 312) and fails the strict one on replicate 5. Without the step both arms failed both readings, since both readings require R-hat under 1.05 in every replicate. For scale, the 10-species rare arm at the same settings failed the strict reading too in the check's own refits (bulk ESS 78 and 71 in replicates 2 and 3).

The species' mean B0 moves with `mu_b0`: its bulk ESS is 155 to 364 in the ordinary fits (21 to 129 before) and 158 to 254 in rare replicates 1 to 4 (6 to 104), and 44 in rare replicate 5 (93). The slow wander of `mu_b0`, the spread of 500-draw batch means as a fraction of its posterior sd, falls from 0.26 to 0.41 to 0.18 to 0.25 in the ordinary fits and from 0.36 to 0.61 to 0.22 to 0.27 in rare replicates 1 to 4; in rare replicate 5 it is 0.42 (0.37 before). The within-chain correlation of `mu_b0` with the mean B0 (0.68 to 0.76 ordinary, 0.80 to 0.87 rare) and of the mean B0 with the mean collection intercept (-0.53 to -0.64, -0.34 to -0.48) are what they were: the step changes how the posterior is explored, not the posterior. B0's share of species with R-hat over 1.05 is 0 to 1% in every fit (0 to 6% before), its median bulk ESS 316 to 412 (303 to 408), and the collection intercepts' 714 to 967 (671 to 1006).

Rare replicate 5 (`traces-mu-b0.png`, bottom right): the chains overlap but each wanders; their means differ by 0.085 in the mean B0, 0.4 posterior sds, and the mean B0 has R-hat 1.05 and bulk ESS 44. The step is not what is slow there: it accepts 41 to 43% of its proposals at an adapted sd of 0.15, as in the other rare fits. Nor is it clearly the detection-side ridge the diagnostic named as the next candidate: the correlation of the 250-draw batch means of the mean B0 and the mean collection intercept is -0.36 in one chain and -0.04 in the other (-0.54 and -0.40 in rare replicate 4, which mixes well), and the mean collection intercept itself has bulk ESS 522. The least-mixed species are rare ones (true B0 -3.3 to -2.1, bulk ESS 54 to 91). The refit kept no spatial or factor draws, so whether the wander rides on the spatial range (which sat at different grid points per chain in every fit of the diagnostic) cannot be read from it; a refit with `refit-drift.R`'s storage would say. With ESS near 100 the estimate itself is rough: before the step this fit scored 120 and after it 64, and the check's rare arm sat near the bar at 10 species as well.

The posterior is unchanged where the earlier chains had mixed (`before-after.csv`): `mu_b0`'s posterior mean moves by at most 0.03 in the seven fits that passed before and by 0.07 to 0.10 in the three that had been stuck (where the earlier mean was the average of two separated chains); its posterior sd by at most 0.026; the community occupancy error by at most 0.0035 in every fit (ordinary -0.027 to -0.001, rare -0.050 to -0.019); B0 coverage by at most 0.02. `mu_b0`'s 95% interval covers the realised mean of the true B0 in 9 of 10 fits (10 of 10 before): rare replicate 4's interval, now (-2.62, -1.97), misses -2.64 by 0.02 where before the lower chain had stretched it.

The step itself (`fits.csv`): acceptance over the kept iterations 0.41 to 0.46 in every chain (0.44 over burn-in, the target); the adapted proposal sd 0.096 to 0.111 in the ordinary fits and 0.116 to 0.157 in the rare ones, against `mu_b0`'s posterior sd of 0.136 to 0.155 and 0.166 to 0.201.

## Cost

The main step's fits took 34.7 to 34.9 minutes (ordinary) and 27.6 to 28.1 (rare) against the check's 33.9 to 34.3 and 27.0 to 27.5 on the same 5 workers: 2 to 3% more. The drift step's fits took 26 to 28 minutes against the diagnostic's 20 to 21 on 3 workers; that run shared the laptop's four performance cores with the App Store process, and the controlled timing below, not the wall clock, is the measure of the step. Timed after the runs on the 100-species rare arm, replicate 3, one chain at 300 + 300 in one process, the step off, on and off again: 53.4, 56.8 and 52.7 ms per iteration, so the step costs 3.8 ms per iteration, 7% of an iteration at 100 species and 100 sites (the isolated timing before the runs, 0.5 ms for the per-pair log-likelihoods and 2 ms for 5 proposals, had put it at 3%). On the validity test's small dataset (8 species, 60 sites, occupancy model) the step adds 29% to a 2 x 2,300 fit, 8.6 to 11.0 s, because the R overhead of a proposal is fixed while the iteration is cheap; there it raises `mu_b0`'s bulk ESS from 139 to 253 and the mean B0's from 81 to 152.

## What it means for PR #36

The decision the check was for, whether `b0_prior = "community"` could become the default, turned on `mu_b0` converging at the default chain length at a realistic number of species. With the step it does by R-hat in every fit of this check and by bulk ESS in all but one, where it is 64 of a required 100; the chain-to-chain shift of every intercept that the diagnostic traced to the latent occupancy states is gone. What remains in that one fit is a slow within-chain wander whose partner the kept draws cannot name. A user who opts in sees `mu_b0`'s R-hat and ESS in `returnConvergenceDiagnostics()` as before. Suggested wording for decision 2: the mixing check at 100 species passes its R-hat criterion in 10 of 10 fits and its ESS criterion in 9 of 10 once `mu_b0` and the intercepts are shifted jointly (this folder), against 7 of 10 on both before; the option's default is a judgement on that one fit.

## Files

-   `README.md`: the step, the design and the rule, written before the runs.
-   `launch.sh`, `refit-drift.R`, `refit-s100.R`, `*-provenance.txt`, `*.log`: the two steps as run.
-   `analyse.R`: writes everything below from the kept draws.
-   `drift.csv`: per diagnostic fit and quantity, the chain difference, posterior sd, R-hat and bulk ESS before and after. `drift-species.csv`: the across-species sd of the B0 shift before and after, the step's acceptance rates and proposal sd, elapsed time.
-   `fits.csv`, `before-after.csv`, `rule.csv`: the 100-species check with the step on, per fit, against the check's values, and the rule.
-   `traces.rds`: per-draw `mu_b0`, mean B0 and mean collection intercept of every fit; `traces-drift.png` and `traces-mu-b0.png`: the trace plots before and after.
-   The full draws (300 MB for `drift`, about 420 MB for `main`) are in the ignored `draws/` folder of this worktree only; `launch.sh` regenerates them exactly.

# The cap on the embedding scales: evidence steps 1 and 2, 25 November 2026

This records evidence steps 1 and 2 of `dev/superpowers/specs/2026-11-25-axis-slab-design.md`, run as Task 6 of `dev/superpowers/plans/2026-11-25-axis-slab.md`, and the cost measurement made in Task 5. The reading rules are the design's, fixed before either run. Step 3 (the five simulated data sets with and without the cap) is not run yet.

## What was run

- Library: built from commit c10ee23 ("Make a cap of 2 the default for the embedding scales", the head of `claude/axis-slab` at the time) with `git archive c10ee23` and `R CMD INSTALL --no-docs` into a separate library, `$SP/lib-evidence`, where `$SP` is the session's scratch directory (`/private/tmp/claude-501/-Users-douglasyu-src-occJSDM--worktrees-axis-slab/d420624c-3ca4-4bb1-ade6-94fff9a195ac/scratchpad/sp`). Both scripts loaded the package from that library and checked that they had. The library was not reinstalled during the runs.
- R 4.6.1, on Doug's laptop. Power: on AC power with Low Power Mode off (`lowpowermode 0`) before and after both steps, battery 95% to 96% and charging. Background load was not recorded.
- Step 1: `prior_check.R`, started and finished 21:56:17 GMT (under a second of computing).
- Step 2: `separation.R`, both arms started detached with `nohup` at 21:57:13 GMT, one R process each and chains one after another in each (`ncores = 1`), so 2 fitting workers in all. No other R fitting job was running. The uncapped arm finished at 21:57:46 GMT (0.55 minutes), the capped arm at 21:58:02 GMT (0.82 minutes).
- Output, which git ignores: `dev/simstudy/results/axis-slab-20261125/` holds `prior_check.out`, `prior_check.rds`, `separation-default.rds`, `separation-Inf.rds`, both arms' logs, and the power state before and after step 2.

## Step 1: the prior check

250 sites, 8 latent traits, 2,000 prior draws of one species' coefficients, the site by column matrix of the earlier check (`dev/simstudy/grouped-prior-checks-20261111/prior_check.R`), seed 20261111 at 64 columns and 20261112 at 128. The prior is read by the package's `read_coef_prior()` and its standard deviations come from the package's `coef_prior_sds()`, so the cap is the package's. The capped arm passes no `slab_scale` and gets the default, checked to be 2. Both arms at one column count use the same site by column matrix and the same random numbers. The number is the standard deviation across sites of the embedding's contribution to the species' linear predictor, in logits. `axis_scale` is the default, `0.7 / K`.

At 64 embedding columns (`axis_scale` 0.01094):

- With the cap (2): 5% point 0.300, median 1.294, 95% point 5.720, 99% point 9.241, largest 20.636.
- Without the cap (`Inf`): 5% point 0.304, median 1.425, 95% point 28.489, 99% point 157.132, largest 18,349.378.

At 128 embedding columns (`axis_scale` 0.00547):

- With the cap (2): 5% point 0.322, median 1.301, 95% point 5.733, 99% point 9.914, largest 23.769.
- Without the cap (`Inf`): 5% point 0.326, median 1.485, 95% point 27.775, 99% point 211.846, largest 820,917.563.

Reading: **passes**. With the cap the median is between 1 and 2 (1.294 and 1.301) and the 95% point under 8 (5.720 and 5.733) at both column counts.

Beside it: without the cap the 95% points are 28.5 and 27.8, near the 27 the design expected. The capped rows agree with the scratch simulation the design quotes (median 1.3 and 95% point 5.7 at 64 columns, 1.26 and 5.7 at 128).

## Step 2: a rare species that one column separates

The data set, built by `separation.R` with seed 3001: 150 sites and 10 embedding columns of independent standard Normals, standardised; species 1 present at exactly the 6 sites with the largest values of column 1 and absent at the other 144; species 2 to 6 present independently with probability 0.3 at each site (44, 43, 50, 49 and 45 presences). The fit: model type binary (one row per site, inferred by the package as one-stage binary data), all 10 columns in `embeddingVars` and `occCovariates`, no measured covariate, no traits, no spatial covariates, `n_lattrait = 2`, `n_factors = 2`, 4 chains of 5,000 burn-in and 5,000 kept sweeps, fit seed 3101. The capped arm passed no `slab_scale` (`listPriors = list()`), so it ran at the default cap, which the script checked to be 2. The other arm passed `listPriors = list(slab_scale = Inf)`.

The number read is the 99% point, over the 20,000 kept draws of all four chains, of the absolute value of species 1's coefficient on column 1.

With the cap (2):

- 99% point 2.349, largest absolute draw 3.214. Per chain, the 99% points are 2.338, 2.322, 2.413 and 2.311. Posterior mean of the coefficient 1.267.
- Reading: **passes** (under 14).
- The parts of that coefficient, 99% point of the absolute value and largest absolute draw: community average 0.866 and 2.001, latent-trait part 2.071 and 3.866, leftover 1.925 and 4.092.
- Effective axis scale of column 1: median 0.418, 5% point 0.158, 95% point 1.115, 99% point 1.569, largest 2.000 (the cap). 0.15% of kept draws are above 1.9.
- `sigma_emb`: median 0.586.
- `returnGroupedConvergence()`: verdict "ok". Every chain in the shared-structure state (median `sigma_emb` by chain 0.489, 0.555, 0.594, 0.698). R-hat of the two shared scales 1.005. R-hat of the axis scales, largest 1.122 (1 of 10 over 1.1), which is not part of the verdict; `runOccJSDM()` also warned that the largest R-hat of log lambda is over 1.1. Measured coefficients, gradients and disagreement not assessed (no measured covariate, and the embedding's effect is small).
- 0.82 minutes.

Without the cap (`Inf`):

- 99% point 2.358, largest absolute draw 3.363. Per chain, the 99% points are 2.435, 2.330, 2.397 and 2.240. Posterior mean of the coefficient 1.274.
- Reading: **no harm seen in this data set** (under 14).
- The parts: community average 0.805 and 1.839, latent-trait part 2.134 and 4.551, leftover 1.981 and 3.534.
- Effective axis scale of column 1 (here the axis scale itself): median 0.396, 5% point 0.156, 95% point 1.032, 99% point 1.623, largest 4.248. 0.57% of kept draws are above 1.9.
- `sigma_emb`: median 0.723.
- `returnGroupedConvergence()`: verdict "ok". Every chain in the shared-structure state (median `sigma_emb` by chain 0.783, 0.598, 0.780, 0.717). R-hat of the two shared scales 1.009, of the axis scales at most 1.005. Measured coefficients, gradients and disagreement not assessed, as above.
- 0.55 minutes.

The two arms differ little. The uncapped prior did not run away on this data set: species 1's coefficient on column 1 stays near 1.3 with a 99% point of about 2.4 in both arms, and only the rare largest draws of the axis scale (over 2 in the uncapped arm) are cut by the cap. As the design said might happen, the case for the cap so far rests on the prior check alone.

## Task 5: the cost of the cap

Measured in Task 5 (its report, `.superpowers/sdd/2026-11-25-axis-slab/task-5-report.md`, outside the repository), not in this task. Data set 1 of the `sigma_b` move study (seed 1001: 256 sites, 48 species, 5 measured covariates, 64 embedding columns), simulated with the call in `dev/simstudy/sigma-b-move-20261122/repro.R`. Each fit is one chain of 300 burn-in and 0 kept sweeps, fit seed 2001, one thread, run as a separate `Rscript --vanilla` process, one at a time, alternating `slab_scale` 2 and `Inf`. The time is the elapsed time of the whole `runOccJSDM()` call. On AC power with Low Power Mode off before and after every set, battery 80% to 88% and charging. No other `Rscript` job was running. The laptop was never idle: system processes (spotlightknowledged, PerfPowerServicesSignpostReader, ReportCrashService, a WebKit process, Backblaze) were near 100% CPU at times during the sets.

With the axis-scale draws in vector form (as committed), two sets of three alternating fits per arm, in seconds:

- Set 1: with the cap 12.749, 12.624, 13.452 (mean 12.942); without 12.302, 12.239, 12.519 (mean 12.353). The cap costs +4.8%.
- Set 2: with the cap 12.860, 12.857, 13.084 (mean 12.934); without 12.381, 12.482, 12.480 (mean 12.448). The cap costs +3.9%.
- Pooled: 12.938 against 12.401, +4.3%.

The loop form first built in Task 2, one slice draw per column, cost +16.7% (with the cap 14.218, 14.181, 14.239; without 12.092, 12.188, 12.251; means 14.213 and 12.177). It was replaced by the vector form for that reason.

Micro-benchmark of the scale draws alone, per sweep: with the loop form, `sample_coef_scales()` took 5.41 ms with the cap against 0.21 ms without, of which the 64 axis-scale draws took 5.24 ms. With the vector form it takes 0.43 ms with the cap against 0.23 ms without, about 0.5% of a sweep of about 41 ms. So the design's requirement, that the scale draws add under 5% to the time of a sweep at 64 columns, holds by a wide margin. The rest of the whole-fit cost is mostly the capped forms of the two scale moves. The whole-fit times include the fixed setup of `runOccJSDM()`, so the cost per sweep is a little higher than the whole-fit percentage, and set 1 may be at or just over 5% per sweep.

# The cap on the embedding scales: evidence steps 1 to 3, 25 and 26 November 2026

This records evidence steps 1 and 2 of `dev/superpowers/specs/2026-11-25-axis-slab-design.md`, run as Task 6 of `dev/superpowers/plans/2026-11-25-axis-slab.md`, the cost measurement made in Task 5, and evidence step 3 (the five simulated data sets with and without the cap), run as Task 7 and recorded in its own section below. The reading rules are the design's, fixed before any run.

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

The capped fit took 0.82 minutes against 0.55 without, about 49% longer. This was not remarked on when step 2 was first recorded. On a fit this small (150 sites, 6 species, 10 columns) the fixed cost in R of each sweep's capped scale draws and moves probably dominates, though that was not measured. It is outside the design's rule that the scale draws add under 5% to a sweep, which is about 64 columns and was measured in Task 5 below.

## Task 5: the cost of the cap

Measured in Task 5 (its report, `.superpowers/sdd/2026-11-25-axis-slab/task-5-report.md`, outside the repository), not in this task. Data set 1 of the `sigma_b` move study (seed 1001: 256 sites, 48 species, 5 measured covariates, 64 embedding columns), simulated with the call in `dev/simstudy/sigma-b-move-20261122/repro.R`. Each fit is one chain of 300 burn-in and 0 kept sweeps, fit seed 2001, one thread, run as a separate `Rscript --vanilla` process, one at a time, alternating `slab_scale` 2 and `Inf`. The time is the elapsed time of the whole `runOccJSDM()` call. On AC power with Low Power Mode off before and after every set, battery 80% to 88% and charging. No other `Rscript` job was running. The laptop was never idle: system processes (spotlightknowledged, PerfPowerServicesSignpostReader, ReportCrashService, a WebKit process, Backblaze) were near 100% CPU at times during the sets.

With the axis-scale draws in vector form (as committed), two sets of three alternating fits per arm, in seconds:

- Set 1: with the cap 12.749, 12.624, 13.452 (mean 12.942); without 12.302, 12.239, 12.519 (mean 12.353). The cap costs +4.8%.
- Set 2: with the cap 12.860, 12.857, 13.084 (mean 12.934); without 12.381, 12.482, 12.480 (mean 12.448). The cap costs +3.9%.
- Pooled: 12.938 against 12.401, +4.3%.

The loop form first built in Task 2, one slice draw per column, cost +16.7% (with the cap 14.218, 14.181, 14.239; without 12.092, 12.188, 12.251; means 14.213 and 12.177). It was replaced by the vector form for that reason.

Micro-benchmark of the scale draws alone, per sweep: with the loop form, `sample_coef_scales()` took 5.41 ms with the cap against 0.21 ms without, of which the 64 axis-scale draws took 5.24 ms. With the vector form it takes 0.43 ms with the cap against 0.23 ms without, about 0.5% of a sweep of about 41 ms. So the design's requirement, that the scale draws add under 5% to the time of a sweep at 64 columns, holds by a wide margin. The rest of the whole-fit cost is mostly the capped forms of the two scale moves. The whole-fit times include the fixed setup of `runOccJSDM()`, so the cost per sweep is a little higher than the whole-fit percentage, and set 1 may be at or just over 5% per sweep.

## Step 3: the same simulated data with and without the cap

### What was run

- Scripts: `compare.R` fits one arm over the five data sets and writes one file per fit; `summarise.R` reads the ten files and prints the numbers and the four readings. `compare.R` was not edited after the runs began.
- Data: the five data sets of the `sigma_b` move study, simulated as in `dev/simstudy/sigma-b-move-20261122/repro.R` (256 sites, 48 species, 5 measured covariates, 64 embedding columns, 2 traits, 4 latent traits, binary), data seed 1000 + k and fit seed 2000 + k for data set k. Fits: 4 chains of 5,000 burn-in and 5,000 kept sweeps, `n_factors = 2`, `n_lattrait = 4`, chains one after another in each fit (`ncores = 1`). One arm passed `listPriors = list(slab_scale = 2)`, the other `listPriors = list(slab_scale = Inf)`, and `compare.R` checked that each fit's prior record held the value asked for.
- Library: `$SP/lib-evidence`, the same one as steps 1 and 2, installed from commit c10ee23 and not reinstalled. The branch's R code is the same at c10ee23 and at 5859874: the two commits between (0af646c and 5859874) change only help, `NEWS.md`, tests and files under `dev/`. The worktree's HEAD, which each file stores, was 0af646c when the uncapped arm ran and 5859874 when the capped arm ran.
- R 4.6.1 on Doug's laptop. Power: on AC power with Low Power Mode off (`lowpowermode 0`) at the start and end of both arms and at the start and end of every fit. `caffeinate -i` held the laptop awake from about 22:10 GMT.
- Times: the uncapped arm first, from 22:06:58 to 23:37:02 GMT on 25 November; then the capped arm, from 23:37:03 GMT on 25 November to 00:38:52 GMT on 26 November. Each arm ran 3 fits at once (data sets 1 to 3, then 4 and 5), detached with `nohup`, and the two arms ran one after the other.
- Other load: while the uncapped arm ran (about 22:10 to 22:45 GMT), other work ran on the laptop: code reviews with short R checks, a reinstall of a scratch library, and two filtered test runs (`axis-slab` and `grouped`). Nothing comparable ran during the capped arm.
- The truth's two numbers per data set were computed from the stored truth after fitting, not printed before fitting as the design says. They depend only on the data seed (1000 + k), so they are the same as they would have been before fitting, and `summarise.R` checks that both arms of each data set have the same truth.
- Output, which git ignores: `dev/simstudy/results/axis-slab-20261125/` holds `slabInf-d1.rds` to `slabInf-d5.rds`, `slab2-d1.rds` to `slab2-d5.rds`, both arms' logs (`step3-slabInf.log`, `step3-slab2.log`) and the power state before and after (`step3-power-before.txt`, `step3-power-after.txt`). Each file also holds the whole result of `returnGroupedConvergence()`, the per-chain posterior means of the embedding coefficients, the axis scales' pooled medians, the draws of `sigma_emb` and `sigma_b`, and every warning of the fit.

### The numbers, per data set

The truth is described by the largest absolute true coefficient on an embedding column and the largest root mean square across species of one column's true coefficients (written "column RMS"). Recovery is the root mean squared error (RMSE), over all 256 sites and 48 species, of the posterior mean of the embedding's contribution to the linear predictor against the true contribution. The shrinkage slope is `sum(est * truth) / sum(truth^2)` over the 48 species of the 5 columns with the largest true column RMS. A fit fails the convergence rule if `returnGroupedConvergence()` flags `rhat_measured`, `gradients` or `disagreement`. Ratios are with the cap over without.

Data set 1 (seeds 1001 and 2001):

- Truth: largest absolute coefficient 3.370; largest column RMS 1.404. The 5 strongest columns are `X_emb.26`, `X_emb.28`, `X_emb.42`, `X_emb.33` and `X_emb.4`, with column RMS 1.404, 1.366, 1.335, 1.141 and 0.987.
- With the cap: RMSE 0.6806; slope 0.7525; verdict "review", flagged `rhat_scales` (1.173), so the fit passes the rule; `rhat_measured` 1.023, `gradients` 0.999, `disagreement` 0.057; 32.1 minutes.
- Without the cap: RMSE 0.6777; slope 0.7425; verdict "ok", nothing flagged, so the fit passes; `rhat_measured` 1.021, `gradients` 0.999, `disagreement` 0.047; 45.3 minutes.
- Ratios: RMSE 1.0043; slope 1.0135.

Data set 2 (seeds 1002 and 2002):

- Truth: largest absolute coefficient 3.966; largest column RMS 1.338. The 5 strongest columns are `X_emb.52`, `X_emb.58`, `X_emb.4`, `X_emb.28` and `X_emb.50`, with column RMS 1.338, 1.106, 1.037, 1.025 and 1.012.
- With the cap: RMSE 0.6404; slope 0.8783; verdict "review", flagged `rhat_scales` (1.107), so the fit passes; `rhat_measured` 1.017, `gradients` 0.999, `disagreement` 0.043; 32.1 minutes.
- Without the cap: RMSE 0.6566; slope 0.8741; verdict "ok", nothing flagged, so the fit passes; `rhat_measured` 1.017, `gradients` 0.999, `disagreement` 0.044; 45.2 minutes.
- Ratios: RMSE 0.9753; slope 1.0048.

Data set 3 (seeds 1003 and 2003):

- Truth: largest absolute coefficient 3.314; largest column RMS 1.114. The 5 strongest columns are `X_emb.26`, `X_emb.12`, `X_emb.28`, `X_emb.20` and `X_emb.1`, with column RMS 1.114, 1.103, 0.875, 0.840 and 0.823.
- With the cap: RMSE 0.7620; slope 0.8518; verdict "review", flagged `state`, `rhat_measured` and `rhat_scales` (9.355), so the fit **fails** the rule; `rhat_measured` 1.408, `gradients` 0.995, `disagreement` 0.154; 32.1 minutes.
- Without the cap: RMSE 0.7237; slope 0.8713; verdict "review", flagged `rhat_scales` (1.163), so the fit passes; `rhat_measured` 1.044, `gradients` 1.000, `disagreement` 0.062; 45.5 minutes.
- Ratios: RMSE 1.0529; slope 0.9777.

Data set 4 (seeds 1004 and 2004):

- Truth: largest absolute coefficient 2.911; largest column RMS 1.173. The 5 strongest columns are `X_emb.10`, `X_emb.22`, `X_emb.16`, `X_emb.58` and `X_emb.14`, with column RMS 1.173, 1.080, 1.026, 0.983 and 0.735.
- With the cap: RMSE 0.7502; slope 0.8060; verdict "review", flagged `state`, `rhat_measured` and `rhat_scales` (6.252), so the fit **fails**; `rhat_measured` 1.743, `gradients` 0.997, `disagreement` 0.118; 29.6 minutes.
- Without the cap: RMSE 0.7480; slope 0.8121; verdict "review", flagged `state`, `rhat_measured` and `rhat_scales` (1.821), so the fit **fails**; `rhat_measured` 1.272, `gradients` 0.998, `disagreement` 0.093; 43.7 minutes.
- Ratios: RMSE 1.0029; slope 0.9925.

Data set 5 (seeds 1005 and 2005):

- Truth: largest absolute coefficient 2.947; largest column RMS 1.109. The 5 strongest columns are `X_emb.12`, `X_emb.43`, `X_emb.5`, `X_emb.29` and `X_emb.40`, with column RMS 1.109, 1.040, 0.982, 0.908 and 0.900.
- With the cap: RMSE 0.7027; slope 0.9177; verdict "review", flagged `rhat_measured` only, so the fit **fails**; `rhat_measured` 1.308, `gradients` 0.992, `disagreement` 0.163; 29.6 minutes.
- Without the cap: RMSE 0.7104; slope 0.8970; verdict "review", flagged `state`, `rhat_measured` and `rhat_scales` (4.870), so the fit **fails**; `rhat_measured` 1.261, `gradients` 0.992, `disagreement` 0.178; 44.7 minutes.
- Ratios: RMSE 0.9892; slope 1.0231.

In every data set the strongest true column has a column RMS of 1.1 to 1.4, 55% to 70% of the cap of 2, and the largest true coefficient is 2.9 to 4.0.

### The readings

- *Recovery*: **no harm**. The RMSE ratios are 1.0043, 0.9753, 1.0529, 1.0029 and 0.9892: at most 1.10 in 5 of 5 data sets (the rule asks for at least 4), and the largest, 1.0529, is under 1.25.
- *Shrinkage of the strongest columns*: **no extra shrinkage**. The slopes with the cap are 0.7525, 0.8783, 0.8518, 0.8060 and 0.9177; without, 0.7425, 0.8741, 0.8713, 0.8121 and 0.8970. The ratios, 1.0135, 1.0048, 0.9777, 0.9925 and 1.0231, are at least 0.95 in 5 of 5 data sets (the rule asks for at least 4). Both priors shrink the strongest columns, by 8% to 26% in the slope, and the cap adds no measurable shrinkage to that.
- *Convergence*: **no worse**. 3 fits fail with the cap (data sets 3, 4 and 5) and 2 without (data sets 4 and 5). The rule allows the count without plus 1, which is 3.
- *Time*: **not slower** by the rule. The mean is 31.1 minutes per fit with the cap (32.1, 32.1, 32.1, 29.6, 29.6) and 44.9 without (45.3, 45.2, 45.5, 43.7, 44.7), a ratio of 0.693. This does not show that the cap is faster. The arms ran one after the other, the uncapped arm first, and the uncapped arm shared the laptop with other work for about 35 of its first 40 minutes. Run order and that extra load are confounded with the arm, and both could favour the capped arm. The Task 5 measurement, made under controlled alternation, found the cap about 4% slower per fit. The time reading here says only that the cap was not more than 10% slower in this run.

None of the four readings fails, so by the design's rule step 3 does not stop the cap of 2 being kept as the default. The design asks for a look at the chains only for a "harm" or "extra shrinkage" reading; neither was read, but the chains were looked at anyway, because the convergence count went up by one.

### What the chains show

From `returnGroupedConvergence()$chains`, the median `sigma_emb` of each chain and its state, with the cap and then without:

- Data set 1: with the cap 0.304, 0.509, 0.455, 0.461, all "shared structure"; without 0.485, 0.559, 0.489, 0.468, all "shared structure".
- Data set 2: with the cap 0.235, 0.296, 0.277, 0.235, all "shared structure"; without 0.310, 0.328, 0.323, 0.310, all "shared structure".
- Data set 3: with the cap 0.621, 0.123, 1.678, 0.189, chain 3 "leftovers dominate" and the rest "shared structure"; without 0.254, 0.168, 0.191, 0.207, all "shared structure".
- Data set 4: with the cap 0.675, 5.668, 0.690, 0.705, chain 2 "leftovers dominate"; without 0.778, 0.689, 0.687, 5.320, chain 4 "leftovers dominate".
- Data set 5: with the cap 0.032, 0.049, 0.041, 0.030, all "shared structure"; without 0.043, 1.303, 0.056, 0.052, chain 2 "leftovers dominate".

So in each arm one chain settled in the "leftovers dominate" state in 2 of the 5 data sets: data sets 3 and 4 with the cap, 4 and 5 without. The extra failing fit with the cap (data set 3) is that one chain; the data set that passes with the cap but not without (data set 5) is the mirror case, though with the cap it still fails on `rhat_measured` (1.308) with all four chains in the same state. This is the noise the design's *Risks* names (chains of this model settle in different states), and with one such chain in 2 of 5 fits in each arm, five data sets cannot tell whether the cap makes it more or less likely.

Beside the rule, and not part of any reading:

- The package's own verdict is "ok" in 2 of 5 fits without the cap (data sets 1 and 2) and in none with it. In those two data sets the capped fit is flagged only for `rhat_scales`, the R-hat of `sigma_emb` and `sigma_b`, at 1.173 and 1.107, with every chain in the same state.
- The largest R-hat of the 64 axis scales (`rhat_axis`, reported but not part of the verdict) is higher with the cap in every data set: 4.264 against 1.443, 1.287 against 1.153, 2.032 against 1.232, 4.053 against 1.732, and 3.608 against 2.452. The number of axis scales with an R-hat over 1.1 is 7, 2, 17, 8 and 19 with the cap and 4, 2, 2, 7 and 21 without. Under the cap the axis scale stored is the effective one, so the two arms' R-hats are of the same quantity only for columns whose scale is well under the cap (the largest pooled median of any axis scale was 1.337 with the cap). The draws of each axis scale were not stored, so which columns drive the larger R-hats with the cap was not looked at. Whether the cap slows the mixing of the axis scales was not part of this step and is not settled by it.

### `summarise.R` output

```
Directory: /Users/douglasyu/src/occJSDM/.worktrees/axis-slab/dev/simstudy/results/axis-slab-20261125 
All 10 files present.

d1 truth: largest |true coefficient| 3.370; largest column RMS 1.404; strongest 5 columns X_emb.26, X_emb.28, X_emb.42, X_emb.33, X_emb.4 (RMS 1.404, 1.366, 1.335, 1.141, 0.987)
d1 slab2   slab_scale 2: RMSE 0.6806; slope 0.7525; verdict review, flagged: rhat_scales -> passes; rhat_measured 1.023, gradients 0.999, disagreement 0.057; 32.1 minutes; power Now drawing from 'AC Power';  lowpowermode         0
d1 slabInf slab_scale Inf: RMSE 0.6777; slope 0.7425; verdict ok, flagged: none -> passes; rhat_measured 1.021, gradients 0.999, disagreement 0.047; 45.3 minutes; power Now drawing from 'AC Power';  lowpowermode         0
d1 ratios cap / none: RMSE 1.0043; slope 1.0135

d2 truth: largest |true coefficient| 3.966; largest column RMS 1.338; strongest 5 columns X_emb.52, X_emb.58, X_emb.4, X_emb.28, X_emb.50 (RMS 1.338, 1.106, 1.037, 1.025, 1.012)
d2 slab2   slab_scale 2: RMSE 0.6404; slope 0.8783; verdict review, flagged: rhat_scales -> passes; rhat_measured 1.017, gradients 0.999, disagreement 0.043; 32.1 minutes; power Now drawing from 'AC Power';  lowpowermode         0
d2 slabInf slab_scale Inf: RMSE 0.6566; slope 0.8741; verdict ok, flagged: none -> passes; rhat_measured 1.017, gradients 0.999, disagreement 0.044; 45.2 minutes; power Now drawing from 'AC Power';  lowpowermode         0
d2 ratios cap / none: RMSE 0.9753; slope 1.0048

d3 truth: largest |true coefficient| 3.314; largest column RMS 1.114; strongest 5 columns X_emb.26, X_emb.12, X_emb.28, X_emb.20, X_emb.1 (RMS 1.114, 1.103, 0.875, 0.840, 0.823)
d3 slab2   slab_scale 2: RMSE 0.7620; slope 0.8518; verdict review, flagged: state, rhat_measured, rhat_scales -> fails; rhat_measured 1.408, gradients 0.995, disagreement 0.154; 32.1 minutes; power Now drawing from 'AC Power';  lowpowermode         0
d3 slabInf slab_scale Inf: RMSE 0.7237; slope 0.8713; verdict review, flagged: rhat_scales -> passes; rhat_measured 1.044, gradients 1.000, disagreement 0.062; 45.5 minutes; power Now drawing from 'AC Power';  lowpowermode         0
d3 ratios cap / none: RMSE 1.0529; slope 0.9777

d4 truth: largest |true coefficient| 2.911; largest column RMS 1.173; strongest 5 columns X_emb.10, X_emb.22, X_emb.16, X_emb.58, X_emb.14 (RMS 1.173, 1.080, 1.026, 0.983, 0.735)
d4 slab2   slab_scale 2: RMSE 0.7502; slope 0.8060; verdict review, flagged: state, rhat_measured, rhat_scales -> fails; rhat_measured 1.743, gradients 0.997, disagreement 0.118; 29.6 minutes; power Now drawing from 'AC Power';  lowpowermode         0
d4 slabInf slab_scale Inf: RMSE 0.7480; slope 0.8121; verdict review, flagged: state, rhat_measured, rhat_scales -> fails; rhat_measured 1.272, gradients 0.998, disagreement 0.093; 43.7 minutes; power Now drawing from 'AC Power';  lowpowermode         0
d4 ratios cap / none: RMSE 1.0029; slope 0.9925

d5 truth: largest |true coefficient| 2.947; largest column RMS 1.109; strongest 5 columns X_emb.12, X_emb.43, X_emb.5, X_emb.29, X_emb.40 (RMS 1.109, 1.040, 0.982, 0.908, 0.900)
d5 slab2   slab_scale 2: RMSE 0.7027; slope 0.9177; verdict review, flagged: rhat_measured -> fails; rhat_measured 1.308, gradients 0.992, disagreement 0.163; 29.6 minutes; power Now drawing from 'AC Power';  lowpowermode         0
d5 slabInf slab_scale Inf: RMSE 0.7104; slope 0.8970; verdict review, flagged: state, rhat_measured, rhat_scales -> fails; rhat_measured 1.261, gradients 0.992, disagreement 0.178; 44.7 minutes; power Now drawing from 'AC Power';  lowpowermode         0
d5 ratios cap / none: RMSE 0.9892; slope 1.0231

Readings
Recovery: RMSE ratio cap / none 1.0043, 0.9753, 1.0529, 1.0029, 0.9892; at most 1.10 in 5 of 5, largest 1.0529 -> no harm
Shrinkage: slopes with the cap 0.7525, 0.8783, 0.8518, 0.8060, 0.9177; without 0.7425, 0.8741, 0.8713, 0.8121, 0.8970; ratio 1.0135, 1.0048, 0.9777, 0.9925, 1.0231; at least 0.95 in 5 of 5 -> no extra shrinkage
Convergence: failing fits with the cap 3 (d3,d4,d5), without 2 (d4,d5) -> no worse
Time: mean minutes with the cap 31.1, without 44.9, ratio 0.693 -> not slower
```

# Free spatial amplitude (15 October 2026)

This follows the [range-inside-grid diagnostic](../range-inside-grid-20261015/README.md). That run found that the fitted spatial amplitude `sigma_bs` sat at the mean of its default prior, 0.33, in every fit, and left one lead open: a spatial term held small by that prior might leave spatial variation unmodelled, narrowing the `B` intervals in the cells that fit traits and a spatial field, and the `B0` intervals in `binary`. Here the amplitude is freed and the same datasets are refitted. No package code changed. The run used the same `R/` and `src/` as the production run of the [paper study](../paper-study-20261015/REPORT.md).

## Answer

- **The pinned amplitude is not the cause.** Under the half-Cauchy prior the amplitude's posterior mean averages 0.82 to 2.68 across the four cells, against 0.33 under the default, although single fits range from 0.23 to 13.5. Yet `B` and `B0` coverage, interval widths and shrinkage do not change. Width ratios against the default prior on the same datasets are 0.99 to 1.01, and the slope of posterior mean on truth moves by at most 0.012. Seven of the eight paired `B` and `B0` coverage changes are within 1.8 SE of zero.
- **The fitted spatial share does not change either.** In `base` the fitted spatial term already carries more of the variation in occupancy probability than the simulated field does (0.27 against 0.18), and it does in both priors. In `binary` it carries less (0.13 against 0.20), and freeing the amplitude does not enlarge it.
- **The spatial part of these fits has not converged, under either prior.** The range's median split R-hat is 1.16 to 1.68 per cell, with a median bulk effective sample size of 3 to 13 out of 2,000 draws, in every cell, including the production cells under the default prior. Under the half-Cauchy prior the amplitude mixes as badly (median R-hat 1.32 to 1.72, effective sample size 3 to 5). Under the default prior, which holds it in place, it mixes well (R-hat 1.01 to 1.03).

So neither the range grid nor the amplitude prior explains the `B` undercoverage, or `binary` `B0`'s. What the fits do show is a spatial term whose range, and whose amplitude once freed, have not mixed within the study's 2 chains of 1,000 draws. Whether the `B` and `B0` intervals are too narrow because of that is the next thing to test. It needs only longer chains, no package change.

## Why only the half-Cauchy prior

The default amplitude prior is InvGamma(10, 1) on the variance: `runOccJSDM()` sets `a_sigmabs = 10` and `b_sigmabs = 1`, and `read_spatial_sd_prior()` returns shape 10 and rate 1 with no option to change them. So that prior cannot be loosened through `listPriors` without a package change. The only alternative `listPriors` offers is `sigma_bs_prior = "half_cauchy"`, with `sigma_bs_scale` defaulting to 1. That prior's median is 1 and its 95th percentile is about 12.7, against a central 95% interval of 0.24 to 0.46 for the default, so it leaves the amplitude essentially free. Its sampler, `sample_spatial_sd_half_cauchy()`, is the standard auxiliary-variable form. The calibration report had already refitted `base` and `binary` with it at replicates 1 to 20 and found no change in widths or coverage, but it did not record the amplitude. The fits here at those replicates are the same fits as those, since the data, seeds and settings match.

## Design

Four harness-only cells at the end of `simstudy_scenarios()` in `tests/testthat/helper-simstudy.R`, each with `listPriors = list(sigma_bs_prior = "half_cauchy")` and its parent's `seed_label`:

- `base_hc` and `binary_hc`: `base` and `binary` at the production range (0.15 raw, about 0.52 standardised).
- `base_ls075_hc` and `binary_ls075_hc`: the same at the in-grid range of 0.075 raw (about 0.26 standardised), as in `base_ls075` and `binary_ls075`.

Each cell sees exactly the datasets of a default-prior cell already on disk. So the prior's effect is a paired comparison on identical data: `base_hc` against production `base` (replicates 1 to 50 of the R = 300 run), and `base_ls075_hc` against `base_ls075` from the range-inside-grid run, and likewise for `binary`. Comparing `base_ls075_hc` with production `base` gives the effect of changing range and prior together.

Before the run, a check confirmed three things: all 36 earlier cells are unchanged, `draw_truth()` returns the same output for `base` and `binary` at replicates 1 to 50, and a short fit returns the half-Cauchy prior and the variance partition in the expected shape. After the run, the `B0`, `B` and `G` truths of every pair matched exactly ([pairing-check.csv](pairing-check.csv)).

## Run

- Launched on 15 October 2026 at 21:47 with [launch.sh](launch.sh), detached, through `dev/simstudy/run_study.R` at revision `8a2c628` on `claude/amplitude-prior-check`, with a clean tree. Its `R/` and `src/` are identical to the production run's `f3668e0`.
- R = 50 per cell, 5 replicate processes, one thread per fit, and 2 chains of 1,000 burn-in and 1,000 retained draws, as in production. R 4.6.1 with Accelerate BLAS, on mains power with Low Power Mode off.
- The 200 fits took 36.1 minutes, and none failed. Another session's study (`spatial-field-2x2`, 5 workers) shared the machine for most of that time, so the elapsed time is not a timing; the results are unaffected, because every fit is single-threaded and seeded.
- The results file records `git_dirty_at_write` as TRUE because the comparison scripts were written, untracked, while the run was in progress ([provenance.txt](provenance.txt), [run.log](run.log)).
- [amp-diag.R](amp-diag.R) then ran from the same launcher, from 22:24 to 22:35 (80 refits).

## Coverage

Coverage is the mean of the per-replicate coverage fractions, and its standard error is their SD over the square root of 50, as in the calibration report's `calibration.R`. Each figure below gives coverage (SE), the paired change against the default prior on the same datasets (SE), the width ratio, and the slope of posterior mean on truth. All rows are in [coverage-r50.csv](coverage-r50.csv) and [paired-vs-default-prior.csv](paired-vs-default-prior.csv).

### `B`

- `base_hc`: 0.930 (0.008) against 0.924; paired +0.006 (0.007); width ratio 1.004; slope 0.683 against 0.685.
- `base_ls075_hc`: 0.930 (0.009) against 0.928; paired +0.002 (0.009); width ratio 1.006; slope 0.679 against 0.682.
- `binary_hc`: 0.917 (0.009) against 0.900; paired +0.017 (0.006); width ratio 1.004; slope 0.875 against 0.869.
- `binary_ls075_hc`: 0.911 (0.008) against 0.916; paired -0.005 (0.005); width ratio 0.993; slope 0.856 against 0.860.

The `binary_hc` change is the only one of the eight `B` and `B0` contrasts beyond 2 SE (2.7 SE). It comes with no wider intervals, and it reverses at the in-grid range, so it gives no sign that the amplitude limits coverage. Error-to-width ratios stay at 1.12 to 1.22. Coverage of true effects of 3 or more in absolute value moves in both directions (0.76 to 0.75 and 0.72 to 0.72 in `base`, 0.78 to 0.83 and 0.86 to 0.82 in `binary`, 72 to 87 elements each; [b-b0-by-effect-size.csv](b-b0-by-effect-size.csv)).

### `B0`

- `base_hc`: 0.928 (0.009) against 0.932; paired -0.004 (0.008); width ratio 1.000.
- `base_ls075_hc`: 0.934 (0.013) against 0.946; paired -0.012 (0.007); width ratio 1.006.
- `binary_hc`: 0.868 (0.015) against 0.864; paired +0.004 (0.008); width ratio 0.996; error-to-width ratio 1.347 against 1.336.
- `binary_ls075_hc`: 0.896 (0.014) against 0.902; paired -0.006 (0.008); width ratio 0.996.

Against production `binary`, `binary_ls075_hc` gains +0.032 (0.015) in `B0`, the same as the range change alone gave (+0.038), and +0.011 (0.011) in `B`.

### Other blocks

`G`, `beta_theta`, `theta0`, `p` and `q` moved by -0.025 to +0.010 in every contrast, all within 2 SE of zero, with width ratios of 0.99 to 1.02.

## Amplitude, range and spatial share

The runner saves none of these, so [amp-diag.R](amp-diag.R) refitted replicates 1 to 10 of eight cells exactly as the runner does. The eight cells are the four half-Cauchy cells and their four default-prior partners. Each refit is the same fit as the run's: its `B` and `B0` posterior means equal the stored rows exactly. The spatial share is the package's variance partition (`computeVariancePartitioning_stdevs()`), which the simulator also computes for the truth: the share of each species' variation in occupancy probability across sites that is attributed to the spatial term, averaged over species. For each pair the list gives the default-prior value first, then the half-Cauchy value ([amp-diag.csv](amp-diag.csv)):

- `base` (true range 0.52): amplitude 0.33 then 1.47 (replicate means 0.37 to 4.8); fitted range 0.148 then 0.126; spatial share 0.27 then 0.27, against a true 0.18.
- `base_ls075` (0.26): amplitude 0.33 then 2.68 (0.37 to 13.5); range 0.143 then 0.140; share 0.26 then 0.29, against 0.20.
- `binary` (0.52): amplitude 0.33 then 0.82 (0.23 to 1.96); range 0.077 then 0.067; share 0.13 then 0.14, against 0.20.
- `binary_ls075` (0.26): amplitude 0.33 then 1.62 (0.39 to 4.0); range 0.065 then 0.051; share 0.12 then 0.11, against 0.21.

Mixing, from split R-hat and bulk effective sample size over the 2 chains of each fit (medians over the 10 replicates):

- Range, default prior: R-hat 1.68, 1.34, 1.36 and 1.26 in `base`, `base_ls075`, `binary` and `binary_ls075`, with effective sample sizes of 3 to 7.
- Range, half-Cauchy prior: R-hat 1.49, 1.16, 1.16 and 1.28, with effective sample sizes of 4 to 13.
- Amplitude, default prior: R-hat 1.01 to 1.03, with effective sample sizes of 92 to 101.
- Amplitude, half-Cauchy prior: R-hat 1.32 to 1.72, with effective sample sizes of 3 to 5.

Freed, the amplitude wanders over a wide range without settling, while the fitted spatial share, and with it the coefficients, stay where they were. The amplitude values above are therefore not reliable estimates, but on average they sit far above the default prior's central 95% interval of 0.24 to 0.46. The range has not converged in any fit, under either prior, so the fitted-range figures here and in the range-inside-grid record are indicative only; a note to that effect has been added to that record.

## What this leaves

- Neither the range lying beyond the fitted grid nor the amplitude prior explains why `B` undercovers, or why `B0` does in `binary`. The intervals, the shrinkage of large effects and the fitted spatial share all stay put when either is changed.
- The study's fits do not mix in their spatial part. The production cells, with their default prior, have a range R-hat of 1.3 to 1.7. The calibration report already noted that no fit had been checked for convergence. If the spatial field is slow to explore, the `B` and `B0` posteriors in each fit may be narrower than the model's true posterior. The direct test is to refit `base` and `binary` at replicates 1 to 20 with much longer and more chains, keep R-hat for `B` and `B0` as well as for the range and amplitude, and see whether the intervals widen and coverage recovers. That needs no package change.

## Limits

- Only the half-Cauchy prior at scale 1 was tried, because the default prior's parameters cannot be changed without a package change.
- R = 50 and two parents. The other cells with fitted traits and a spatial field were not rerun.
- The mixing figures rest on 2 chains per fit and 10 replicates per cell. With effective sample sizes this small, R-hat itself is rough, but values above 1.1 in every cell leave no doubt that the range has not converged.

## Files and reproduction

- [launch.sh](launch.sh), [provenance.txt](provenance.txt) and [run.log](run.log) record the run and the refits. `results/` holds the runner's output (`simstudy-20261015-222358.rds` and its CSV summaries) and the refits' `amp-diag.rds`; the runner's checkpoint is not committed. `launch.sh` has this worktree's absolute paths and refuses to start if `results/` exists. It was started from this directory, in a non-interactive shell, as shown below.
- [compare.R](compare.R) writes `coverage-r50.csv`, `paired-vs-default-prior.csv`, `b-b0-by-effect-size.csv` and `pairing-check.csv`. It reads the default-prior rows from the R = 300 production run and from the range-inside-grid run.
- [amp-diag.R](amp-diag.R) ran the refits; it refuses to overwrite its output, and its log, `amp-diag.log`, is empty. [amp-diag-summary.R](amp-diag-summary.R) writes `amp-diag.csv` from the refits, with R-hat and effective sample sizes from the `posterior` package, and checks each refit against the run's rows.

From this directory:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh > launch.log 2>&1 < /dev/null &
```

From the repository root:

```sh
Rscript dev/simstudy/amplitude-prior-20261015/compare.R
Rscript dev/simstudy/amplitude-prior-20261015/amp-diag-summary.R dev/simstudy/amplitude-prior-20261015/results/amp-diag.rds
```

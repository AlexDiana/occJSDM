# Spatial range inside the fitted grid (15 October 2026)

This tests the leading untested explanation in the [calibration report](../paper-study-20261015/REPORT.md) of the 15 October 2026 paper study. That report found `B` undercovering in the cells that fit traits and a spatial field, and `B0` undercovering in `binary`, and suggested the cause might be the simulated spatial range: 0.15 on the raw unit square, which is 0.47 to 0.58 on the fitter's standardised coordinates, beyond the fitted range grid's top of 0.30. No package code changed. The run used the same `R/` and `src/` as the production run.

## Answer

- **`B` does not recover.** With the range inside the grid, `B` coverage in `base` is 0.928 and 0.914 against 0.924 for the same production replicates, and in `binary` 0.916 and 0.896 against 0.900. Every paired change is within 1.6 SE of zero. The intervals are as wide as before (width ratios 0.98 to 1.00), the estimates are shrunk as hard as before (slope of posterior mean on truth 0.67 to 0.69 in `base`, 0.85 to 0.87 in `binary`), and the error-to-width ratio stays at 1.13 to 1.26.
- **`binary` `B0` recovers only partly, and not through wider intervals.** Coverage rises from 0.864 to 0.902 at 0.075 raw (paired +0.038, SE 0.014) and to 0.890 at 0.05 raw (+0.026, SE 0.016), which is still 3.5 and 4.7 SE below 0.95. The intervals are about 3% narrower, and the gain comes from estimates about 9% closer to the truth (RMSE 0.391 against 0.355 and 0.356). The gain is larger at the range near the grid's top than at the one deeper inside it, so it does not track the grid.
- **The fitted range does not follow the simulated one.** Refits of replicates 1 to 10 that keep the range draws put the posterior mean of the fitted range at 0.15 in `base` and 0.08 in `binary` when the truth is 0.52. That is well inside the grid, not pressed against its top: only 9% and 1% of draws are at 0.30. Moving the truth down to 0.26 and 0.17 moves the fitted range only to 0.14 and 0.12 in `base`, and to 0.07 and 0.05 in `binary`. The spatial amplitude stays at about its prior mean, 0.33, in all six cells. These range figures come from chains that have not converged for the range (median R-hat 1.14 to 1.68), so they are indicative only; see the note under *Fitted range and amplitude*.

So the range mismatch is not why `B` undercovers, and it explains at most a small part of `binary` `B0`'s shortfall. The explanation's premise does not hold either: in these fits the range is not trying to reach beyond the grid, and the fitted spatial term barely responds to the true range, although the range draws have not converged.

## Design

Four harness-only cells, defined at the end of `simstudy_scenarios()` in `tests/testthat/helper-simstudy.R`:

- `base_ls05` and `binary_ls05`: copies of `base` and `binary` with `l_s = 0.05` raw, which is 0.157 to 0.194 standardised over both axes and replicates 1 to 50, mid-grid.
- `base_ls075` and `binary_ls075`: the same with `l_s = 0.075` raw, 0.236 to 0.290 standardised, near the grid's top. A first choice of 0.08 reached 0.31 on some replicates, outside the grid, and was dropped before any fit.

There are two ranges so that a recovery at the shorter one could not be put down to the field merely being rougher.

`draw_truth()` now reads an optional scenario `l_s`, defaulting to 0.15, and adds no random draws. Each new cell sets `seed_label` to its parent's label, so it uses the parent's seed at every replicate. The simulator's random-number use does not depend on `l_s`, which only rescales the kernel before the same `chol()` and `rnorm()` calls. So the site coordinates and every JSDM coefficient (`B0`, `B`, `G`, `A`, `C`, `L`, `U`) are bit-identical to the parent's, and only the spatial field, and therefore the data, differ. The field's spread does not change: its SD is about 0.69 in `base` and 0.71 in `binary` at all three ranges.

[check-pairing.R](check-pairing.R) verified three things before the run, with no fitting. For all ten production cells at replicates 1 to 50, `draw_truth()` returns output identical to the pre-change harness. Under the runner's L'Ecuyer-CMRG generator, the `base` and `binary` truths at replicates 1 to 3 match the production rows. And the pairing above holds at replicates 1 to 50. [standardised-range.csv](standardised-range.csv) has the per-replicate standardised ranges. After the run, the `B0`, `B` and `G` truths of each cell matched its parent's production rows exactly ([pairing-check.csv](pairing-check.csv)).

## Run

- Launched on 15 October 2026 at 20:59 with [launch.sh](launch.sh), detached, through `dev/simstudy/run_study.R` at revision `f4d7a64` on `claude/range-inside-grid`. Its `R/` and `src/` are identical to those of the production run's `f3668e0`.
- R = 50 per cell, 5 replicate processes, one thread per fit, and 2 chains of 1,000 burn-in and 1,000 retained draws, as in production.
- R 4.6.1 with Accelerate BLAS, on mains power with Low Power Mode off. The 200 fits took 23.7 minutes, and none failed.
- The tree was clean when the run started. The results file records `git_dirty_at_write` as TRUE because the comparison scripts were written, untracked, while the run was in progress ([provenance.txt](provenance.txt), [run.log](run.log)).

## Comparison

Each cell is compared with replicates 1 to 50 of its parent in the R = 300 production run (`dev/simstudy/paper-study-20261015/results/`), which share its seeds and coefficients. As in the calibration report's `calibration.R`, coverage is the mean of the per-replicate coverage fractions, and its standard error is their SD over the square root of 50. Paired differences are computed per replicate. Replicates 1 to 50 are representative of the R = 300 run: `base` `B` 0.924 against 0.927, `binary` `B` 0.900 against 0.902, and `binary` `B0` 0.864 against 0.868.

Each figure below gives coverage (SE), mean interval width, the slope of posterior mean on truth, and the error-to-width ratio. That ratio is the RMSE of the error over the posterior SD that the interval implies, and it is near 1 when intervals are calibrated. Paired differences are against the parent.

### `B`

- `base` (production, true range 0.52 standardised): 0.924 (0.008), width 3.215, slope 0.685, ratio 1.144.
- `base_ls075` (0.26): 0.928 (0.010), width 3.227, slope 0.682, ratio 1.132. Paired +0.004 (0.010), width ratio 1.004.
- `base_ls05` (0.17): 0.914 (0.010), width 3.155, slope 0.669, ratio 1.175. Paired -0.010 (0.011), width ratio 0.981.
- `binary` (production, 0.52): 0.900 (0.010), width 1.390, slope 0.869, ratio 1.206.
- `binary_ls075` (0.26): 0.916 (0.008), width 1.374, slope 0.860, ratio 1.217. Paired +0.016 (0.010), width ratio 0.989.
- `binary_ls05` (0.17): 0.896 (0.010), width 1.357, slope 0.849, ratio 1.255. Paired -0.004 (0.010), width ratio 0.977.

The shortfall is still concentrated in large effects ([b-b0-by-effect-size.csv](b-b0-by-effect-size.csv), pooled elements, no SE). In `base` and its 0.075 and 0.05 arms, true effects of 3 or more in absolute value are covered 0.76, 0.72 and 0.67 of the time (87 elements each). In `binary` the figures are 0.78, 0.86 and 0.82 (72 elements). The average pull towards zero on those effects is unchanged: 1.44 to 1.51 in `base`, 0.54 to 0.61 in `binary`.

### `B0`

- `binary` (production): 0.864 (0.015), width 1.167, slope 0.806, ratio 1.336, RMSE 0.391.
- `binary_ls075`: 0.902 (0.014), width 1.139, slope 0.802, ratio 1.224, RMSE 0.355. Paired +0.038 (0.014), width ratio 0.975.
- `binary_ls05`: 0.890 (0.013), width 1.125, slope 0.787, ratio 1.233, RMSE 0.356. Paired +0.026 (0.016), width ratio 0.964.
- `base`, for completeness: 0.932 (0.010) in production, and 0.946 (0.012) and 0.946 (0.010) in the two arms, each a paired +0.014 (0.009), with widths unchanged. `B0` in `base` is mostly prior-driven (slope 0.33 to 0.35), as the calibration report explains.

### Other blocks

`G`, `beta_theta`, `theta0`, `p` and `q` moved by -0.005 to +0.015, all within 1.4 SE of zero, with width ratios of 0.98 to 1.01 ([paired-vs-production.csv](paired-vs-production.csv)).

## Fitted range and amplitude

The runner saves neither, so [range-diag.R](range-diag.R) refitted replicates 1 to 10 of the two parents and the four new cells (60 fits) exactly as the runner does, keeping the range and amplitude draws. Each refit is the same fit as the run's: its `B` and `B0` posterior means equal the production or diagnostic rows exactly. For each cell the list gives the true standardised range; the posterior mean of the fitted range, averaged over replicates, with its range across replicates; the share of draws at the grid's top (0.30) and bottom (0.01); and the mean amplitude `sigma_bs` ([range-diag.csv](range-diag.csv)):

- `base`: true 0.52, fitted 0.148 (0.079 to 0.223), 9% at top, 11% at bottom, `sigma_bs` 0.328.
- `base_ls075`: true 0.26, fitted 0.143 (0.042 to 0.259), 11% at top, 16% at bottom, `sigma_bs` 0.329.
- `base_ls05`: true 0.17, fitted 0.124 (0.039 to 0.248), 7% at top, 18% at bottom, `sigma_bs` 0.328.
- `binary`: true 0.52, fitted 0.077 (0.019 to 0.182), 1% at top, 35% at bottom, `sigma_bs` 0.330.
- `binary_ls075`: true 0.26, fitted 0.065 (0.014 to 0.143), none at top, 34% at bottom, `sigma_bs` 0.325.
- `binary_ls05`: true 0.17, fitted 0.045 (0.011 to 0.141), none at top, 45% at bottom, `sigma_bs` 0.330.

The fitted range falls a little as the truth falls, by 16% in `base` and 42% in `binary`, while the truth falls threefold. In every cell it sits well below the truth. Only 16 to 17% of `base`'s draws, and under 3% of `binary`'s, fall within one grid step of an in-grid truth. The amplitude does not move off its prior mean. On this design the fitted spatial term does not track the true range, so moving that range inside the grid cannot change the fit much, which is what the coverage shows.

**Note added later on 15 October 2026:** the range draws have not converged. Split R-hat and bulk effective sample size, computed afterwards from the same draws by [range-diag-summary.R](range-diag-summary.R), put the range's median R-hat at 1.14 to 1.68 in the six cells, with a median effective sample size of 3 to 13 out of 2,000 draws. So each fit's range summary is unreliable, and the averages above are only indicative. The amplitude mixes well under its default prior (median R-hat 1.01 to 1.03). The coverage comparison still describes what the study's fits give. But slow mixing of the spatial part is itself a candidate cause of the narrow `B` intervals, which the [amplitude follow-up](../amplitude-prior-20261015/README.md) on `claude/amplitude-prior-check` takes up.

## Limits

- R = 50 and two parents. A change of about 2 points in `B` coverage would have been detectable, and none was seen. The other cells with fitted traits and a spatial field were not rerun; the calibration report found that they share `base`'s pattern.
- The 100 random sites and 20 support points of production are kept. On the standardised scale those support points are about 0.8 apart, farther than any range in the grid, so the fitted spatial term may be unable to resolve a short range whatever the grid. The 0.075 arm, whose field is smoother, gives the same `B` result as the 0.05 arm.
- The range and amplitude figures rest on 10 replicates per cell, and the range draws have not converged (see the note above).
- This run does not establish what does cause the `B` undercoverage. It leaves one lead open: the amplitude stays at its prior mean in every cell, so a spatial term held small by that prior could leave spatial variation unmodelled. The half-Cauchy refits in the calibration report changed neither widths nor coverage, but they did not record the amplitude, so they do not settle this.

## Files and reproduction

- [check-pairing.R](check-pairing.R) runs the pre-run checks and writes `standardised-range.csv`.
- [launch.sh](launch.sh), [provenance.txt](provenance.txt) and [run.log](run.log) record the run. `results/` holds the runner's output (`simstudy-20261015-212334.rds` and its CSV summaries); the runner's checkpoint is not committed. `launch.sh` has this worktree's absolute paths and refuses to start if `results/` exists. It was started from this directory, in a non-interactive shell, as shown below.
- [compare.R](compare.R) writes `coverage-r50.csv` (every cell and block, with replicate-based SEs), `paired-vs-production.csv`, `b-b0-by-effect-size.csv` and `pairing-check.csv`.
- [range-diag.R](range-diag.R) ran the range refits (`results/range-diag.rds`, logged in `range-diag.log`; it refuses to overwrite its output), and [range-diag-summary.R](range-diag-summary.R) writes `range-diag.csv` from them. The log's start line counts two untracked files: `compare.R` and `range-diag.R` itself, committed afterwards unchanged.

From this directory:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh > launch.log 2>&1 < /dev/null &
```

From the repository root:

```sh
Rscript dev/simstudy/range-inside-grid-20261015/check-pairing.R
Rscript dev/simstudy/range-inside-grid-20261015/compare.R
Rscript dev/simstudy/range-inside-grid-20261015/range-diag.R dev/simstudy/range-inside-grid-20261015/results/range-diag.rds
Rscript dev/simstudy/range-inside-grid-20261015/range-diag-summary.R dev/simstudy/range-inside-grid-20261015/results/range-diag.rds
```

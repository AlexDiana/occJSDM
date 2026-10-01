# Site-arrangement sweep for the spatial field (Lesson 2)

This study asks how much of a shared spatial field occJSDM can recover when the survey budget is fixed at 100 sites and only the arrangement of those sites changes: spread at random, spread with close pairs, clustered, or a regular grid. It compares an oracle sampler that knows the true parameters, occJSDM fitted to true occupancy states, and occJSDM fitted to a simulated two-stage eDNA survey, in three independent communities with eight species each. The protocol is in `PLAN.md`; the study was completed on 1 October 2026.

## Results

All numbers below are read from `results/` and are means over the three communities unless a range is given. Three communities cannot support confidence intervals, so none are claimed.

Reading labels from the prespecified rules, applied to the binary-control fits: 11 of the 12 arrangement and species-group cells are uninformative. The one exception is the clustered design for the 25% species group, which is intermediate: the RMSE reductions in the three communities are 9.3% to 11.8%, below the 20% needed in all three to call it informative, but not below 10% in all three, so it is not uninformative either (the correlations are 0.64 or more). The clustered design for the 75% group is uninformative because the reductions in all three communities are below 10% (7.9% to 8.7%), although its correlations are 0.52 or more. Spread, pairs and grid are uninformative for all three groups.

Field recovery at the surveyed sites, as the RMSE reduction relative to a flat field (binary control; oracle in brackets; eDNA survey after the slash), for the 5%, 25% and 75% groups. Spread: 0.3% (-0.2%) / -0.1%, 1.9% (8.8%) / 0.3%, 2.1% (10.2%) / 0.6%. Pairs: 0.4% (-0.5%) / -0.1%, 3.0% (13.2%) / 1.0%, 2.7% (13.4%) / 0.7%. Clustered: 2.9% (10.6%) / 0.4%, 10.5% (27.9%) / 2.6%, 8.4% (19.7%) / 2.7%. Grid: 0.4% (1.4%) / 0.1%, 1.6% (7.9%) / 0.4%, 1.6% (7.5%) / 0.3%. The binary-control median-field correlations with the truth are 0.11, 0.38 and 0.40 for spread, 0.11, 0.46 and 0.46 for pairs, 0.46, 0.69 and 0.61 for clustered, and 0.18, 0.40 and 0.36 for the grid. Even the oracle, which knows the true range and amplitude, recovers at most 28% of the field RMSE (clustered, 25% group) and about 10% for a plain spread design.

Range recovery fails in every arrangement and both arms. The posterior mass within one grid step of the true range (standardised value about 0.10) never reaches 0.5 in any of the 24 fits: it ranges from 0.18 to 0.37, and the arrangement means are 0.30 to 0.33 for both arms. Posterior mass at the two ends of the range grid has arrangement means of 0.14 to 0.20. The spatial amplitude is also pulled below its true value of 1: the posterior median is 0.31 to 0.36 and the mean 95% interval runs from about 0.24 to 0.53 across arrangements and arms.

Estimation cost, the centred field RMSE of the binary control minus the oracle's, averages 0.070 on the log-odds scale and is positive in 30 of the 36 community, arrangement and group cells. By arrangement it is 0.049 for spread, 0.069 for pairs, 0.117 for clustered and 0.044 for the grid, and by group 0.017, 0.102 and 0.089 for the 5%, 25% and 75% species. Detection cost, the two-stage RMSE minus the binary control's, averages 0.022 and is positive in 35 of 36 cells: 0.011 for spread, 0.015 for pairs, 0.051 for clustered and 0.009 for the grid. Imperfect parameter estimation therefore costs about three times as much field accuracy as imperfect detection, and the eDNA survey recovers at most 2.7% of the field RMSE in any cell, with median-field correlations of 0.41 or lower.

Prediction at the 1,600 unsurveyed lattice locations, mean absolute error in probability points by distance to the nearest surveyed site (up to 0.02, 0.02 to 0.05, 0.05 to 0.1, above 0.1), using the spatial term. Binary control, spread: 11.68, 11.65, 11.69, 11.98. Pairs: 11.78, 11.71, 11.86, 12.01. Clustered: 9.96, 11.08, 11.40, 11.51. Grid: 11.75, 11.69, 11.61, and no cells above 0.1 because the 10 by 10 grid leaves no location further than 0.071 from a site. Removing the spatial term changes the error by at most 0.17 points for spread, pairs and grid; for the clustered design it adds 0.73 points in the nearest bin and 0.11 in the next, and nothing beyond. Error does not grow with distance from the survey in any design, because the fitted field is too weak to be informative at any distance. The two-stage fits have mean absolute errors of 18.7 to 21.6 points in every bin and a positive bias of 5.7 to 8.5 points; removing their spatial term changes the error by under 0.1 points.

## Convergence

Three of the 24 initial fits met the prespecified rule for a longer run: `rep01-spread-two_stage` (group and element Rhat above 1.05 and a native convergence warning), `rep01-grid-two_stage` (a native convergence warning) and `rep03-spread-two_stage` (occupancy species mean ESS below 100). Each was rerun with four chains, 6,000 burn-in and 12,000 retained draws, and the longer fit replaces the initial one in every summary. After the longer runs no fit carries a flag: the maximum group Rhat is 1.012 for the three longer fits and 1.033 for the 21 initial fits, there are no native warnings, and no occupancy group or species mean ESS falls below 100 (the minima are 140.5 and 119.6). The spatial amplitude is the weak point: ten of the 24 selected fits have a spatial amplitude mean ESS below 100 in `results/groups.csv`, with a minimum of 49.5 for `rep03-clustered-binary` (then 68.5 for `rep02-pairs-binary`, 72.8 for `rep02-clustered-two_stage` and 80.7 for `rep01-clustered-binary`). These do not trigger a longer run, because the rule covers occupancy groups and species only, and the amplitude quantiles quoted above come from fits with such ESS values, so they are imprecise. The longer runs moved the occupancy group bias and mean absolute error by at most 0.013 in probability, and by at most 0.036 across all scored metrics, the largest being the logit-scale intercept (`results/long-run-sensitivity.csv`). The audit passed for all 24 selected fits and the 96 oracle results: probabilities, range tables, and the range mass within one step and at the grid ends, rebuilt from the saved draws, agree with the scored values exactly, the field medians and lattice means reproduce through the scoring functions exactly, and a basis for 50 lattice cells built from the kernel alone agrees with the package's prediction matrix to within 6e-12.

## Reproduction

From the repository root, with `STUDY` an absolute path to a fresh directory under `dev/simstudy/results/`. The long runs, the oracle and the initial fits are launched detached with `nohup` so that they survive the session.

```sh
STUDY=$PWD/dev/simstudy/results/spatial-design-20261001
Rscript dev/simstudy/spatial-design-sweep/test-generator.R
Rscript dev/simstudy/spatial-design-sweep/test-oracle.R
Rscript dev/simstudy/spatial-design-sweep/test-score.R
Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=$STUDY --mode=prepare
nohup Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=$STUDY --mode=oracle --workers=2 > $STUDY/oracle.log 2>&1 &
Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=$STUDY --mode=freeze
nohup Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=$STUDY --mode=initial --workers=8 > $STUDY/initial.log 2>&1 &
Rscript dev/simstudy/spatial-design-sweep/summarise.R --repo=. --study=$STUDY --mode=select
KEYS=$(paste -sd, $STUDY/long-keys.txt)
nohup Rscript dev/simstudy/spatial-design-sweep/run.R --repo=. --study=$STUDY --mode=long --workers=8 --keys=$KEYS > $STUDY/long.log 2>&1 &
Rscript dev/simstudy/spatial-design-sweep/summarise.R --repo=. --study=$STUDY --mode=final
Rscript dev/simstudy/spatial-design-sweep/verify.R --repo=. --study=$STUDY
Rscript dev/simstudy/spatial-design-sweep/plot.R --repo=. --study=$STUDY
Rscript dev/simstudy/spatial-design-sweep/export-teaching.R --repo=. --study=$STUDY
```

Each `nohup` step must finish before the next command starts (the logs end with a completion line). `export-teaching.R` is listed in `PLAN.md` as the last step and is added by the lesson task, not by this study. `run.R` never overwrites an existing input, result or fit (it checks that the saved one matches and reuses it), and `summarise.R --mode=final` refuses a non-empty output directory; `verify.R` and `plot.R` rewrite `audit/` and `results/` on each run.

## Files

Scripts:

- `PLAN.md`: the protocol fixed before any simulation, with dated amendments and the decision log.
- `generator.R`: simulates the landscapes, the four arrangements, the survey data and the design statistics.
- `oracle.R`: the oracle sampler wrapper, with the true states and parameters supplied.
- `score.R`: scoring of every full fit, including the independent basis reconstruction and lattice prediction.
- `run.R`: the driver, with modes prepare, oracle, freeze, pilot, initial and long.
- `summarise.R`: selects the fits needing longer runs (`--mode=select`) and builds the final tables (`--mode=final`); also defines `diagnostic_reasons()`, `read_sweep_selection()` and `reading_labels()`.
- `verify.R`: the independent audit of every selected fit and of the oracle result hashes, including the fully independent lattice-basis check.
- `plot.R`: copies the compact results into `results/` and draws the three figures.
- `test-generator.R`, `test-oracle.R`, `test-score.R`: the numerical validation tests listed in `PLAN.md`.

Compact results in `results/`:

- `reading.csv`: the prespecified reading label for each arrangement and species group.
- `range-reading.csv`: whether the range was recovered, by arrangement and arm.
- `aggregate.csv`: field RMSE and correlation by arrangement, arm and species group, with the minimum and maximum over the three communities.
- `paired.csv`: the estimation cost and detection cost for every community, arrangement and group.
- `groups.csv`, `species.csv`: occupancy accuracy and convergence diagnostics for each fit, by group and by species.
- `field.csv`, `oracle.csv`: field recovery per species for the full fits and for the oracle.
- `range.csv`, `amplitude.csv`: posterior range and amplitude summaries for each fit.
- `lattice.csv`: prediction error at the unsurveyed lattice, overall and by distance bin, with and without the spatial term.
- `selected-fits.csv`, `long-run-sensitivity.csv`: which result file represents each fit, the reasons for any longer run, and how the longer runs changed the group results.
- `design-statistics.csv`: nearest-neighbour distance, neighbour counts, effective rank and standardised ranges for each arrangement.
- `audit.csv`: the audit differences for each selected fit and whether it passed.
- `field-recovery.png`, `lattice-prediction.png`, `range-recovery.png`: the three figures.

The raw archive (fits, draws, inputs, the frozen package library and logs) lives at `dev/simstudy/results/spatial-design-20261001/`. It is gitignored and is needed only for re-auditing the draws; everything quoted above is in the compact results.

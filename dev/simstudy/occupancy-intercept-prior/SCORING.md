# Scoring phase A

How the controller scores and audits the new arms after the scorer has passed review. The definitions are those of [README.md](README.md), [AMENDMENT-1.md](AMENDMENT-1.md) and [AMENDMENT-2.md](AMENDMENT-2.md); departures are recorded in [DEVIATIONS.md](DEVIATIONS.md). Substitute `REPO`, `STUDY` and `INPUTS` as given in the README. Run every step from a non-interactive shell in the detached form of AMENDMENT-1.md, and start a step only when the previous one has written its `.exit` file holding 0. The timings below were measured on the control arm while eight long fits were running.

## Scripts

- `analysis.R` scores each selected fit of the requested arms from its saved draws and writes one record per fit to `STUDY/summary/<phase>/scores`. Arm 1 is the control. A record of the same fit and the same scorer code is resumed, so the control records made during validation are reused. It refuses new arms before `STUDY/selection/<phase>/selected-fits.csv` exists, and before the other sub-phase's selection exists.
- `summarise.R` reads only those records and the selection tables, and writes `band-error.csv`, `coverage.csv`, `mae.csv`, `b0.csv`, `convergence.csv`, `summary-means.csv`, `b0-means.csv`, `provenance.csv` and, with the control and a new arm, `gate.csv`, `gate-detail.csv` and `schedule-matched.csv` to `STUDY/summary/<phase>`. Once both sub-phase gates exist it also writes `STUDY/summary/A/gate.csv`, with one phase A row per SD. It prints no verdict.
- `verify.R` recomputes, with its own code and from the saved draws, every per-community summary, the stored convergence diagnostics, and from those every across-community table: `summary-means.csv`, `b0-means.csv`, `convergence.csv`, `gate-detail.csv`, `gate.csv` and `schedule-matched.csv`. It exits 1 on any disagreement. It writes `verify-arms-<arms>.csv` (one row per fit), `tables-arms-<arms>.csv` (one row per table check), `gate-audit-arms-<arms>.csv` (its own gate) and `audit-source.csv` (the md5 of `verify.R`, `jobs.R` and every table it read) to `STUDY/verify/<phase>`. With `--phase=A` it recomputes the phase A rows from its own A1 and A2 gate results.
- `supplement.R` is used only in the R27 case below.
- `crosscheck-controls.R` compares the control outcomes with the archived control results; it was run during validation and need not be rerun.
- A verdict is read from `gate.csv` only after `verify.R` has exited 0 for the same phase and arms, and the phase A verdict only after `verify.R --phase=A` has exited 0.

## Ordering rule

AMENDMENT-1.md says that no occupancy error is computed before both `STUDY/selection/A1/selected-fits.csv` and `STUDY/selection/A2/selected-fits.csv` exist. Under controller ruling R29 (post-look, see DEVIATIONS.md), the official scoring of the new arms waits for both selections, exactly as that rule requires. `analysis.R` and `supplement.R` refuse to score a new arm before then, and the `--allow-before-both-selections` option of `analysis.R` is not used. Only the control arm was scored earlier, for validation (AMENDMENT-2.md, and DEVIATIONS.md for its basis).

## The gate's tolerance

Criteria 2 and 3 are "at most" comparisons: band `abs(E)` worse by at most 1 point, MAE by at most 0.5 point, band coverage lower by at most 0.03. They are evaluated as the difference being at most the threshold plus `GATE_TOL = 1e-9`, so that a difference equal to the threshold in exact arithmetic is not failed by binary rounding: for example coverage .53 minus .50 is 0.030000000000000027 in doubles. The allowance is far below any attainable nonzero difference in these outcomes. Strict comparisons (criterion 1, and criterion 4's counts) use none. Controller ruling R30 accepts it, subject to the task review. `gate-detail.csv` has a `tolerance_decisive` column, TRUE where the allowance decided a comparison, that is, where the difference lies above the threshold but within `GATE_TOL` of it.

A criterion evaluated over no community (a band with no cells in any community, or an arm without valid fits) is NA, never PASS. A sub-phase result is FAIL if any criterion fails, UNDETERMINED if none fails but one is NA, and PASS otherwise. Only PASS passes.

## Phase A1

Run these only after the A2 selection (AMENDMENT-1.md Step 6) has finished.

Step A1-1, score the 60 selected new-arm fits; the 20 control records are resumed. The 20 controls took 11 minutes at two workers, so allow about 35 minutes:

```sh
mkdir -p STUDY/summary/logs STUDY/verify/logs
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/analysis.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phase=A1 --arms=1,2,3,5 --workers=2; echo $? > STUDY/summary/logs/analysis-a1.exit' < /dev/null > STUDY/summary/logs/analysis-a1.out 2>&1 &
echo $! > STUDY/summary/logs/analysis-a1.pid
```

Step A1-2, the tables and the A1 gate (a few seconds):

```sh
Rscript REPO/dev/simstudy/occupancy-intercept-prior/summarise.R --repo=REPO --study=STUDY --phase=A1 --arms=1,2,3,5 > STUDY/summary/logs/summarise-a1.out 2>&1; echo $? > STUDY/summary/logs/summarise-a1.exit
```

Step A1-3, the independent audit of all 80 fits and of every table. The 20 controls took 3 minutes at two workers, so allow about 15 minutes:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/verify.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phase=A1 --arms=1,2,3,5 --workers=2; echo $? > STUDY/verify/logs/verify-a1.exit' < /dev/null > STUDY/verify/logs/verify-a1.out 2>&1 &
echo $! > STUDY/verify/logs/verify-a1.pid
```

## Phase A2

After `STUDY/logs/select-a2.exit` holds 0, run the same three steps with `--phase=A2` in place of `--phase=A1` and the log names `analysis-a2`, `summarise-a2` and `verify-a2`. The 9 controls took 9 minutes to score at two workers, and 20 minutes to audit at two workers, or about 40 at one. Allow about 30 minutes for the analysis of the 27 new fits and about 1.5 hours for the audit of all 36 at two workers; use one worker if memory is short. The A2 summarise step also writes `STUDY/summary/A/gate.csv`, since `STUDY/summary/A1/gate.csv` exists by then; `summarise.R --repo=REPO --study=STUDY --phase=A` rewrites it from the two sub-phase gates.

## Phase A

Once both sub-phase audits have exited 0 with `--arms=1,2,3,5`, audit the phase A rows:

```sh
Rscript REPO/dev/simstudy/occupancy-intercept-prior/verify.R --repo=REPO --study=STUDY --phase=A > STUDY/verify/logs/verify-a.out 2>&1; echo $? > STUDY/verify/logs/verify-a.exit
```

## A community without a valid selected fit (R27)

If a first fit of some SD fails for good (R25), `select.R` refuses to record any selection that includes that SD, because its `check_first_fits` refuses while a first fit is missing. Frozen `select.R` and `flags.R` are not edited. Instead:

- Write `STUDY/selection/<phase>/manual-missing.csv` with the columns `phase`, `sd`, `key` and `reason`, one row per community without a valid selected fit, naming the log or quarantine file that shows the cause.
- Record the selection of the other SDs with `select.R` and `--sds` limited to them.
- Write `STUDY/selection/<phase>/manual-selected.csv` with the columns `phase`, `sd`, `key`, `schedule`, `fit` and `fit_md5`: one row for every other community of the affected SD, with the fit the selection rule selects. `fit` is relative to `STUDY`, for example `fits/A1/sd2/initial/jsdm-n0300-05-fit.rds`.
- Once both selections exist, run `supplement.R` (below). It applies `flags.R`'s `evaluate_item` with the frozen flag rule to each listed fit, and to the first fit of each longer repeat, exactly as `select.R` does. It checks the selection rule, writes `STUDY/summary/<phase>/manual-flags.csv` and `manual-selection.csv` once, and scores the listed fits.
- Then run the analysis, summarise and audit steps as above with `--arms=1,2,3,5`. The affected SD fails criterion 4 of its sub-phase, and its criteria 1 to 3 are reported, descriptively, over its valid communities. An SD with a `manual-missing.csv` row but no `manual-selected.csv` rows is reported with criteria 1 to 3 NA. `verify.R` reads the same manual record and audits the affected SD's valid fits and its gate.

```sh
Rscript REPO/dev/simstudy/occupancy-intercept-prior/supplement.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phase=A1 --workers=2 > STUDY/summary/logs/supplement-a1.out 2>&1; echo $? > STUDY/summary/logs/supplement-a1.exit
```

## Reading the result

`STUDY/summary/A/gate.csv` holds, per SD, the A1 row, the A2 row and the phase A row. The phase A result is PASS only when both sub-phase rows are PASS, and `smallest_passing` marks the smallest SD that passes. Each sub-phase row gives every number the gate used:

- `gate.csv` columns are named `c1_<group>_` for criterion 1, `c2_<band>_<stratum>_` and `c3_<band>_<stratum>_` for criteria 2 and 3, and `c4_<stratum>_` for criterion 4, followed by `c1_pass` to `c4_pass` and `result`.
- Criterion 1: the number of communities, the number improved, the number needed and both across-community means. For A1 these are pooled over the n100 and n300 fits of each generating community (R21).
- Criteria 2 and 3: both arms' means and their difference in each band and stratum.
- Criterion 4: both arms' flagged selected fits, and the new arm's communities without a valid selected fit (R27).
- `gate-detail.csv` holds the same numbers one criterion per row, with `tolerance_decisive`.
- `schedule-matched.csv` holds criteria 1 to 3 with each new arm's first fits. In phase A these are the selected fits, unless an A1 first fit was replaced by its longer repeat.

# Scoring phase A

How the controller scores and audits the new arms after the scorer has passed review. The definitions are those of [README.md](README.md), [AMENDMENT-1.md](AMENDMENT-1.md) and [AMENDMENT-2.md](AMENDMENT-2.md). Substitute `REPO`, `STUDY` and `INPUTS` as given in the README. Run every step from a non-interactive shell in the detached form of AMENDMENT-1.md, and start a step only when the previous one has written its `.exit` file holding 0. The timings below were measured on the control arm at two workers while eight long fits were running.

## Scripts

- `analysis.R` scores each selected fit of the requested arms from its saved draws and writes one record per fit to `STUDY/summary/<phase>/scores`. Arm 1 is the control. A record of the same fit and the same scorer code is resumed, so the control records made during validation are reused. It refuses new arms before `STUDY/selection/<phase>/selected-fits.csv` exists.
- `summarise.R` reads only those records and the selection tables, and writes `band-error.csv`, `coverage.csv`, `mae.csv`, `b0.csv`, `convergence.csv`, `summary-means.csv`, `b0-means.csv`, `provenance.csv` and, with the control and a new arm, `gate.csv`, `gate-detail.csv` and `schedule-matched.csv` to `STUDY/summary/<phase>`. Once both sub-phase gates exist it also writes `STUDY/summary/A/gate.csv`, with one phase A row per SD.
- `verify.R` recomputes every per-community summary and the stored convergence diagnostics from the saved draws with its own code, and exits 1 on any disagreement. Its result is `STUDY/verify/<phase>/verify-arms-<arms>.csv`.
- `crosscheck-controls.R` compares the control outcomes with the archived control results; it was run during validation and need not be rerun.
- `gate.csv` is read only after `verify.R` has exited 0 for the same phase and arms.

## Ordering rule

AMENDMENT-1.md says that no occupancy error is computed before both `STUDY/selection/A1/selected-fits.csv` and `STUDY/selection/A2/selected-fits.csv` exist. AMENDMENT-2.md relaxes this for the control arm only. `analysis.R` therefore refuses to score a new arm while the other sub-phase's selection is missing. To score the A1 new arms before the A2 selection exists, first record a dated amendment that permits it, and only then add `--allow-before-both-selections=TRUE` to the A1 analysis command. Without such an amendment, run the A1 commands below after the A2 selection (AMENDMENT-1.md Step 6) has finished.

## Phase A1

Step A1-1, score the 60 selected new-arm fits; the 20 control records are resumed. The 20 controls took 11 minutes, so allow about 35 minutes:

```sh
mkdir -p STUDY/summary/logs STUDY/verify/logs
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/analysis.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phase=A1 --arms=1,2,3,5 --workers=2; echo $? > STUDY/summary/logs/analysis-a1.exit' < /dev/null > STUDY/summary/logs/analysis-a1.out 2>&1 &
echo $! > STUDY/summary/logs/analysis-a1.pid
```

Step A1-2, the tables and the A1 gate (a few seconds):

```sh
Rscript REPO/dev/simstudy/occupancy-intercept-prior/summarise.R --repo=REPO --study=STUDY --phase=A1 --arms=1,2,3,5 > STUDY/summary/logs/summarise-a1.out 2>&1; echo $? > STUDY/summary/logs/summarise-a1.exit
```

Step A1-3, the independent audit of all 80 fits. The 20 controls took 2.5 minutes, so allow about 15 minutes:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/verify.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phase=A1 --arms=1,2,3,5 --workers=2; echo $? > STUDY/verify/logs/verify-a1.exit' < /dev/null > STUDY/verify/logs/verify-a1.out 2>&1 &
echo $! > STUDY/verify/logs/verify-a1.pid
```

## Phase A2

After `STUDY/logs/select-a2.exit` holds 0, run the same three steps with `--phase=A2` in place of `--phase=A1` and the log names `analysis-a2`, `summarise-a2` and `verify-a2`. The 9 controls took 9 minutes to score and 20 minutes to audit, so allow about 30 minutes for the analysis of the 27 new fits and about 1.5 hours for the audit of all 36. The A2 summarise step also writes `STUDY/summary/A/gate.csv` when `STUDY/summary/A1/gate.csv` exists; `summarise.R --repo=REPO --study=STUDY --phase=A` rewrites it from the two sub-phase gates.

## Reading the result

`STUDY/summary/A/gate.csv` holds, per SD, the A1 row, the A2 row and the phase A row. The phase A result is PASS only when both sub-phase rows are PASS, and `smallest_passing` marks the smallest SD that passes. Each sub-phase row gives every number the gate used: for criterion 1 the number of communities, the number improved, the number needed and both across-community means (for A1 pooled over the n100 and n300 fits of each generating community, R21); for criteria 2 and 3 both arms' means and their difference in each band and stratum; for criterion 4 both arms' flagged selected fits and the new arm's communities without a valid selected fit (R27). `gate-detail.csv` holds the same numbers one criterion per row, and `schedule-matched.csv` the criteria 1 to 3 with each new arm's first fits (in phase A these are the selected fits unless an A1 first fit was replaced by its longer repeat).

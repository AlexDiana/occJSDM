# Scoring phases A and B

How the controller scores and audits the new arms after the scorer has passed review, for phase A and then phase B. The definitions are those of [README.md](README.md), [AMENDMENT-1.md](AMENDMENT-1.md) and [AMENDMENT-2.md](AMENDMENT-2.md); departures are recorded in [DEVIATIONS.md](DEVIATIONS.md). Substitute `REPO`, `STUDY` and `INPUTS` as given in the README. Run every step from a non-interactive shell in the detached form of AMENDMENT-1.md, and start a step only when the previous one has written its `.exit` file holding 0. The timings below were measured on the control arm while eight long fits were running.

## Scripts

- `analysis.R` scores each selected fit of the requested arms from its saved draws and writes one record per fit to `STUDY/summary/<phase>/scores`. Arm 1 is the control. A record of the same fit and the same scorer code is resumed, so the control records made during validation are reused. It refuses new arms before `STUDY/selection/<phase>/selected-fits.csv` exists, and before the other sub-phase's selection exists.
- `summarise.R` reads only those records and the selection tables, and writes `band-error.csv`, `coverage.csv`, `mae.csv`, `b0.csv`, `convergence.csv`, `summary-means.csv`, `b0-means.csv`, `provenance.csv` and, with the control and a new arm, `gate.csv`, `gate-detail.csv` and `schedule-matched.csv` to `STUDY/summary/<phase>`. Once both sub-phase gates exist it also writes `STUDY/summary/A/gate.csv`, with one phase A row per SD. It prints no verdict.
- `verify.R` recomputes, with its own code and from the saved draws, every per-community summary, the stored convergence diagnostics, and from those every across-community table: `summary-means.csv`, `b0-means.csv`, `convergence.csv`, `gate-detail.csv`, `gate.csv` and `schedule-matched.csv`. It exits 1 on any disagreement. It writes `verify-arms-<arms>.csv` (one row per fit), `tables-arms-<arms>.csv` (one row per table check), `gate-audit-arms-<arms>.csv` (its own gate) and `audit-source.csv` (the md5 of `verify.R`, `jobs.R` and every table it read) to `STUDY/verify/<phase>`. With `--phase=A` it recomputes the phase A rows from its own A1 and A2 gate results.
- `analysis-b.R` scores phase B in the same way, into `STUDY/summary/B/scores`. It is a separate file so that `analysis.R` stays byte-identical: the md5 of `analysis.R` is one of the scorer hashes of every phase A record, so editing it would invalidate all phase A records and force them to be rescored. Phase B records carry their own version and hash set, which covers both files and the archived phase B scorers. It refuses new arms before `STUDY/selection/B/selected-fits.csv` exists.
- `supplement.R` is used only in the R27 case below, in phase A or phase B.
- `crosscheck-controls.R` compares the control outcomes with the archived control results; it was run during validation of each phase and need not be rerun.
- A verdict is read from `gate.csv` only after `verify.R` has exited 0 for the same phase and arms, the phase A verdict only after `verify.R --phase=A` has exited 0, the phase B verdict only after `verify.R --phase=B` has exited 0, and the final decision only after `verify.R --phase=final` has exited 0.

## Ordering rule

AMENDMENT-1.md says that no occupancy error is computed before both `STUDY/selection/A1/selected-fits.csv` and `STUDY/selection/A2/selected-fits.csv` exist. Under controller ruling R29 (post-look, see DEVIATIONS.md), the official scoring of the new arms waits for both selections, exactly as that rule requires. `analysis.R` and `supplement.R` refuse to score a new arm before then, and the `--allow-before-both-selections` option of `analysis.R` is not used. Only the control arm was scored earlier, for validation (AMENDMENT-2.md, and DEVIATIONS.md for its basis).

In phase B no new-arm occupancy error is computed before `select.R --phases=B --mode=final` has written `STUDY/selection/B/selected-fits.csv` (README); `analysis-b.R` and `supplement.R` refuse before then. The 20 phase B controls were scored earlier for validation, as the phase A controls were (DEVIATIONS.md).

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

## Phase B

Phase B runs only for the SDs that pass phase A (README), and its gate is separate from the phase A gate. After the phase A verdict that is SD 2, 3 and 5; if fewer SDs pass, replace `--sds=2,3,5` by those SDs and `--arms=1,2,3,5` by `1` and those SDs. The definitions differ from phase A only as the protocol says:

- Scored cells: the original 100 sites of all 10 species, and every fitted site of a rare species (mean truth over the 300 fitted sites below 0.2, replicate 07 only). There is no all-site scope.
- Criterion 1 (R24): the community is the replicate; its value is the mean of `abs(E)` in the low band over its qnear and qfar fits; at least 7 of 10 communities improved and a strictly lower across-community mean. The low band is the only gated group; the rare group is descriptive (R22).
- Criteria 2 to 4 hold at qnear and at qfar separately, each over its 10 fits. `select.R` writes `convergence.csv` with one stratum `all` for phase B, so `summarise.R` counts the flags per contamination level from `selected-fits.csv` and checks that they add up to the `select.R` totals; the control flags 1 qnear and 3 qfar selected fits.
- Secondary outcomes: `beta-theta.csv` and `beta-theta-means.csv` hold the collection intercept bias, absolute bias and coverage on the standardised collection-covariate scale, and the posterior correlation of `B0` and the collection intercept, per community and averaged over communities.

Step B-1, after the phase B launcher has finished (its `.exit` file holds 0 and its `DONE` file reports `exit_status: 0`), the selection plan, written to `STUDY/selection/B`:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/select.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=B --sds=2,3,5 --mode=plan --workers=2; echo $? > STUDY/logs/select-b-plan.exit' < /dev/null > STUDY/logs/select-b-plan.out 2>&1 &
echo $! > STUDY/logs/select-b-plan.pid
```

Step B-2, only if `STUDY/selection/B/long-keys.txt` lists any fit, the single longer repeats of the flagged initial-schedule new-arm fits (with an empty list this command reports an empty queue and exits 0). Fits that start at the longer schedule are never repeated:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/launch.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=B --sds=2,3,5 --max-procs=8 --only=STUDY/selection/B/long-keys.txt --schedule-override=long; echo $? > STUDY/logs/launch-b-long.exit' < /dev/null > STUDY/logs/launch-b-long.out 2>&1 &
echo $! > STUDY/logs/launch-b-long.pid
```

Step B-3, the final selection, after Step B-2 has finished (its `.exit` file holds 0 and its `DONE` file reports `exit_status: 0`), or directly after Step B-1 when `long-keys.txt` is empty:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/select.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=B --sds=2,3,5 --mode=final --workers=2; echo $? > STUDY/logs/select-b-final.exit' < /dev/null > STUDY/logs/select-b-final.out 2>&1 &
echo $! > STUDY/logs/select-b-final.pid
```

Step B-4, score the selected new-arm fits and the first fits of any longer repeats; the 20 control records are resumed. The 20 controls took 30 minutes at one worker while eight fits were running: about 30 seconds per initial-schedule fit and 2 to 3 minutes per longer-schedule fit. The archived design scorer needs about 8.4 GB at its peak for a longer-schedule fit, so use two workers only when no fits are running, and one otherwise. For the 60 new fits (24 at the initial and 36 at the longer schedule) and the first fits of any longer repeats, allow about 1 hour at two workers or 2 hours at one:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/analysis-b.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phase=B --arms=1,2,3,5 --workers=2; echo $? > STUDY/summary/logs/analysis-b.exit' < /dev/null > STUDY/summary/logs/analysis-b.out 2>&1 &
echo $! > STUDY/summary/logs/analysis-b.pid
```

Step B-5, the tables and the phase B gate, written to `STUDY/summary/B` (a few seconds):

```sh
Rscript REPO/dev/simstudy/occupancy-intercept-prior/summarise.R --repo=REPO --study=STUDY --phase=B --arms=1,2,3,5 > STUDY/summary/logs/summarise-b.out 2>&1; echo $? > STUDY/summary/logs/summarise-b.exit
```

Step B-6, the independent audit of every phase B fit and table. The 20 controls took 12 minutes at one worker while eight fits were running (about 1 minute per longer-schedule fit, at about 2 GB), so allow about 50 minutes at one worker, or 25 at two, for all 80 selected fits and the first fits of any longer repeats:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/verify.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phase=B --arms=1,2,3,5 --workers=2; echo $? > STUDY/verify/logs/verify-b.exit' < /dev/null > STUDY/verify/logs/verify-b.out 2>&1 &
echo $! > STUDY/verify/logs/verify-b.pid
```

The phase B verdict, `STUDY/summary/B/gate.csv`, is read only after `STUDY/verify/logs/verify-b.exit` holds 0.

## Final decision

Once `STUDY/verify/logs/verify-a.exit` and `STUDY/verify/logs/verify-b.exit` both hold 0, write the final decision table `STUDY/summary/final/decision.csv` and audit it (a few seconds each):

```sh
Rscript REPO/dev/simstudy/occupancy-intercept-prior/summarise.R --repo=REPO --study=STUDY --phase=final > STUDY/summary/logs/summarise-final.out 2>&1; echo $? > STUDY/summary/logs/summarise-final.exit
Rscript REPO/dev/simstudy/occupancy-intercept-prior/verify.R --repo=REPO --study=STUDY --phase=final > STUDY/verify/logs/verify-final.out 2>&1; echo $? > STUDY/verify/logs/verify-final.exit
```

`decision.csv` has one row per SD: `phase_a_result`, `phase_b_result` (`NOT RUN` for an SD that did not pass phase A), `passes_both`, `recommended`, `consequence` and `study_decision`. The consequences are the README's: an SD that passes both phases passes; one that passes phase A and fails phase B is a binary-only improvement and does not change the default; an UNDETERMINED phase B is not a pass, so that SD does not change the default either; the smallest SD that passes both phases is recommended, and any default change is made in a separate reviewed pull request; if no SD passes both, the default stays at 1. `verify.R --phase=final` recomputes the table from its own phase A and phase B gate audits, which must have passed, and it is read only after `STUDY/verify/logs/verify-final.exit` holds 0.

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

In phase B the procedure is the same with `--phase=B` and `STUDY/selection/B`; `supplement.R` then scores with `analysis-b.R`'s definitions, and the scoring, summary and audit steps are B-4 to B-6.

## Reading the result

`STUDY/summary/B/gate.csv` holds one phase B row per SD, with columns named as below; its strata are `qnear` and `qfar`, and criterion 1 is pooled over them (R24). `STUDY/summary/final/decision.csv` holds the final decision per SD.

`STUDY/summary/A/gate.csv` holds, per SD, the A1 row, the A2 row and the phase A row. The phase A result is PASS only when both sub-phase rows are PASS, and `smallest_passing` marks the smallest SD that passes. Each sub-phase row gives every number the gate used:

- `gate.csv` columns are named `c1_<group>_` for criterion 1, `c2_<band>_<stratum>_` and `c3_<band>_<stratum>_` for criteria 2 and 3, and `c4_<stratum>_` for criterion 4, followed by `c1_pass` to `c4_pass` and `result`.
- Criterion 1: the number of communities, the number improved, the number needed and both across-community means. For A1 these are pooled over the n100 and n300 fits of each generating community (R21).
- Criteria 2 and 3: both arms' means and their difference in each band and stratum.
- Criterion 4: both arms' flagged selected fits, and the new arm's communities without a valid selected fit (R27).
- `gate-detail.csv` holds the same numbers one criterion per row, with `tolerance_decisive`.
- `schedule-matched.csv` holds criteria 1 to 3 with each new arm's first fits. In phase A these are the selected fits, unless an A1 first fit was replaced by its longer repeat.

# Amendment 1 to the frozen protocol

Dated 29 September 2026. This amendment was made after the protocol freeze (commit `31d5e9a`) and before any scientific fit of a new arm was started or any outcome was computed: at this date the study archive holds only the three pilot fits, which are never scored, and no selection has been recorded. It was revised in place the same day after the task review, still before any fit, because it had not yet governed any run. Where this amendment and [README.md](README.md) differ, this amendment applies; the rest of the README is unchanged.

## Status of the rulings

Doug approved in chat, on 29 September 2026, the candidate SDs 2, 3 and 5, the controller rulings R1 to R15, and the longer schedule for every phase A2 fit. The README's "Approvals" section says that he approved the rulings R1 to R18; that is inaccurate. R16, R17 and R18 were recorded after his approval of R1 to R15, and they, like R19 to R23 and R25 below, are controller rulings made under the subagent-driven execution Doug requested. They are recorded here before any fit or outcome and are pending his confirmation.

The same review confirmed three points of the frozen text without change: A2 longer fits are flagged by the A2 archive's robust-v1 rule (R19); the convergence criterion stays as frozen, so in phase A every selected new-arm fit must be unflagged (R20); and the phase B rare group is reported descriptively and not gated (R22).

## R21: phase A1 is judged over its 10 generating communities

Reason: each 300-site A1 dataset contains the 100-site dataset of the same replicate, with identical observations and generating probabilities at sites 1 to 100. The 20 A1 fits are therefore 10 generating communities, each fitted at two sample sizes, not 20 independent communities. Counting 20 in the first gate criterion would count each community twice.

- Communities: replicate r, for r from 01 to 10, is one community with two fits, `jsdm-n0100-r` (100 sites) and `jsdm-n0300-r` (300 sites). Both fits are scored on the same original 100 sites with the same generating probabilities, so each band holds the same cells in both.
- Criterion 1 in A1 (improvement). The gated group is the below-20% band only; the rare group is not gated in A1 (R18). For each community r and for each arm (the new arm and the control), the community value is the mean of the absolute signed error of its two fits, (`abs(E(n100))` + `abs(E(n300))`) / 2, with E defined as in the README. Criterion 1 passes only if the new arm's community value is strictly smaller than the control's in at least 7 of the 10 communities (ceil of two thirds of 10) and the across-community mean of the new arm's community values is strictly smaller than that of the control's.
- No-harm criteria in A1 (criteria 2, 3 and 4 of the README) apply at 100 sites and at 300 sites separately, each over its 10 communities: no band's across-community mean absolute signed error worse by more than 1 point, overall MAE not worse by more than 0.5 point, no band's across-community mean coverage lower by more than 0.03, and the number of flagged selected fits of the new arm no more than the control's. Each must hold at both sizes.
- Phase A1 passes only if criterion 1, pooled as above, passes and all no-harm criteria pass at both 100 and 300 sites.
- Phase A2 is unchanged: 9 communities, criterion 1 needs at least 6 of 9 (ceil of two thirds of 9) in each gated group, and criteria 2 to 4 apply over the 9. Phase A passes only if A1 and A2 both pass. Phase B is unchanged by this amendment.
- `select.R` writes `convergence.csv` with a `stratum` column: `n100` and `n300` for A1, `all` for the other phases, so that criterion 4 is counted at each A1 size. No gate or scoring code exists yet; the Task 4 scorer implements R21 as stated here.

## R25: failed or quarantined fits

A fit that the launcher reports as failed or quarantined is diagnosed from its log and, if present, its quarantine file (`<key>-fit.QUARANTINE-<time>.rds`, which holds the complete fit and the reason).

- If the cause is infrastructural (a killed process, a full disk, a stale lock or another failure outside the model and its data), the fit is rerun once with identical settings and seed, by rerunning the same step's command after the recovery described below; `run.R` resumes every completed fit and refits only the missing ones.
- If the fit fails a post-fit invariant, or errors reproducibly, it is not rerun again: that arm counts the community as flagged for the convergence criterion, and the case is reported.

## Correction: what the selection reads

This replaces the README's statement, in "Convergence flags and selection", that `flags.R` "computes no estimation error: generating truth enters only to decide which cells belong to a band or group". The accurate description is as follows. No occupancy error is computed or read by the selection, and none is used. Generating truth enters the flag computation in two ways. In every phase it decides which cells belong to a band or group. In A2 the archived robust rule also screens the `field_projection` traces, which project each draw of a species' spatial field onto that species' centred true field, and the archived function that computes the spatial screens (`robust_field_draws`) also returns a field-recovery summary against the true field, which the flag code discards unread. The recorded control selection is read from archived files. For A2 that is the spatial-amplitude `fits.csv`, which also holds the controls' archived occupancy bias, MAE and coverage columns: the selection reads that file but uses only its schedule, fit path and md5 columns. The cross-check `check-flags.R` likewise uses only the reason strings of that file.

## R23: phase A1 is launched, selected and repeated before phase A2

Reason: the 60 A1 fits take seconds each, while the 27 A2 fits take about 10 hours. Launching A1 first lets its selection and any single longer repeats finish within about an hour, instead of waiting for the A2 batch: the launcher refuses to run while another launcher holds its lock, so a combined launch would hold every A1 longer repeat back until all A2 fits were done. The fits, schedules and selection rule are unchanged; only the order of launches changes. The two launches are the two halves of `results/phase-a-queue.csv`, recorded by dry runs as `results/phase-a1-queue.csv` (60 initial fits) and `results/phase-a2-queue.csv` (27 longer fits).

- `select.R` now defaults its output to `STUDY/selection/<phases>`, for example `STUDY/selection/A1`, so the A1 and A2 selections are recorded separately and cannot overwrite each other. Both are recorded before any occupancy error is computed.
- The commands below replace the "Running phase A" section of the README. Substitute `REPO`, `STUDY` and `INPUTS` as given in the README. Each command runs in its own session (perl's `setsid`, so that it survives the end of the session that started it, not only a hangup), writes its output to `STUDY/logs/NAME.out`, its process id to `STUDY/logs/NAME.pid` and, when it ends, its exit status to `STUDY/logs/NAME.exit`. Start each step only when the previous one has finished: its `.exit` file exists and holds 0, and, for a launcher, the `DONE` file of its log directory `STUDY/logs/launch-<timestamp>/` also reports `exit_status: 0`.

Step 1, the A1 initial fits:

```sh
mkdir -p STUDY/logs
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/launch.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A1 --sds=2,3,5 --max-procs=8; echo $? > STUDY/logs/launch-a1.exit' < /dev/null > STUDY/logs/launch-a1.out 2>&1 &
echo $! > STUDY/logs/launch-a1.pid
```

Step 2, the A1 selection, written to `STUDY/selection/A1`:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/select.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A1 --sds=2,3,5 --mode=plan --workers=4; echo $? > STUDY/logs/select-a1-plan.exit' < /dev/null > STUDY/logs/select-a1-plan.out 2>&1 &
echo $! > STUDY/logs/select-a1-plan.pid
```

Step 3, only if `STUDY/selection/A1/long-keys.txt` lists any fit, the single longer repeats (with an empty list this command reports an empty queue and exits 0):

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/launch.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A1 --sds=2,3,5 --max-procs=8 --only=STUDY/selection/A1/long-keys.txt --schedule-override=long; echo $? > STUDY/logs/launch-a1-long.exit' < /dev/null > STUDY/logs/launch-a1-long.out 2>&1 &
echo $! > STUDY/logs/launch-a1-long.pid
```

Step 4, the A1 final selection:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/select.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A1 --sds=2,3,5 --mode=final --workers=4; echo $? > STUDY/logs/select-a1-final.exit' < /dev/null > STUDY/logs/select-a1-final.out 2>&1 &
echo $! > STUDY/logs/select-a1-final.pid
```

Step 5, the A2 longer fits:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/launch.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A2 --sds=2,3,5 --max-procs=8; echo $? > STUDY/logs/launch-a2.exit' < /dev/null > STUDY/logs/launch-a2.out 2>&1 &
echo $! > STUDY/logs/launch-a2.pid
```

Step 6, the A2 selection, written to `STUDY/selection/A2`. A2 fits are never escalated, so the final mode follows the plan directly:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/select.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A2 --sds=2,3,5 --mode=plan --workers=4 && Rscript REPO/dev/simstudy/occupancy-intercept-prior/select.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A2 --sds=2,3,5 --mode=final --workers=4; echo $? > STUDY/logs/select-a2.exit' < /dev/null > STUDY/logs/select-a2.out 2>&1 &
echo $! > STUDY/logs/select-a2.pid
```

No occupancy error is computed before both `STUDY/selection/A1/selected-fits.csv` and `STUDY/selection/A2/selected-fits.csv` exist.

## Watching a step and recovering a dead launcher

- A step is running while `kill -0 "$(cat STUDY/logs/NAME.pid)"` succeeds; the pid is that of the step's shell, whose child is the R process. A launcher's own R process id is also in `launcher.txt` in its log directory and in `STUDY/logs/launcher.lock/owner`.
- The `.exit` file appears when the step ends. For a launcher, a non-zero value with a `DONE` file means some fits failed or were quarantined (apply R25 to the fits `summary.csv` lists); a non-zero value without `DONE` means the launcher itself died (128 plus the signal number if it was killed). If the pid is gone and no `.exit` file exists, the step's shell was killed.
- A launcher that dies leaves no `DONE` and a stale `STUDY/logs/launcher.lock`, and the fits it had started keep running. To recover without disturbing them: list the fit processes still running with `pgrep -fl occupancy-intercept-prior/run.R` and wait until none remains (each finishes, saves its fit and releases its key lock); check that the launcher pid in `STUDY/logs/launcher.lock/owner` no longer runs (`kill -0` fails), then remove that directory by hand; remove any key lock `STUDY/fits/<phase>/sd<sd>/<schedule>/<key>.lock` whose owner pid no longer runs, after checking it; then rerun the same step's command. `run.R` resumes every completed fit with identical settings and fits the others with the same settings and seeds.

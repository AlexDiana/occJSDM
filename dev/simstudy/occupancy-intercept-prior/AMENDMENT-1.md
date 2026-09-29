# Amendment 1 to the frozen protocol

Dated 29 September 2026. This amendment was made after the protocol freeze (commit `31d5e9a`) and before any scientific fit of a new arm was started or any outcome was computed: at this date the study archive holds only the three pilot fits, which are never scored, and no selection has been recorded. It records two controller rulings, R21 and R23, made on review of the frozen protocol and approved under Doug's standing approval of controller rulings. Where this amendment and [README.md](README.md) differ, this amendment applies; the rest of the README is unchanged.

The same review confirmed three points of the frozen text without change: A2 longer fits are flagged by the A2 archive's robust-v1 rule (R19); the convergence criterion stays as frozen, so in phase A every selected new-arm fit must be unflagged (R20); and the phase B rare group is reported descriptively and not gated (R22).

## R21: phase A1 is judged over its 10 generating communities

Reason: each 300-site A1 dataset contains the 100-site dataset of the same replicate, with identical observations and generating probabilities at sites 1 to 100. The 20 A1 fits are therefore 10 generating communities, each fitted at two sample sizes, not 20 independent communities. Counting 20 in the first gate criterion would count each community twice.

- Communities: replicate r, for r from 01 to 10, is one community with two fits, `jsdm-n0100-r` (100 sites) and `jsdm-n0300-r` (300 sites). Both fits are scored on the same original 100 sites with the same generating probabilities, so each band holds the same cells in both.
- Criterion 1 in A1 (improvement). The gated group is the below-20% band only; the rare group is not gated in A1 (R18). For each community r and for each arm (the new arm and the control), the community value is the mean of the absolute signed error of its two fits, (`abs(E(n100))` + `abs(E(n300))`) / 2, with E defined as in the README. Criterion 1 passes only if the new arm's community value is strictly smaller than the control's in at least 7 of the 10 communities (ceil of two thirds of 10) and the across-community mean of the new arm's community values is strictly smaller than that of the control's.
- No-harm criteria in A1 (criteria 2, 3 and 4 of the README) apply at 100 sites and at 300 sites separately, each over its 10 communities: no band's across-community mean absolute signed error worse by more than 1 point, overall MAE not worse by more than 0.5 point, no band's across-community mean coverage lower by more than 0.03, and the number of flagged selected fits of the new arm no more than the control's. Each must hold at both sizes.
- Phase A1 passes only if criterion 1, pooled as above, passes and all no-harm criteria pass at both 100 and 300 sites.
- Phase A2 is unchanged: 9 communities, criterion 1 needs at least 6 of 9 (ceil of two thirds of 9) in each gated group, and criteria 2 to 4 apply over the 9. Phase A passes only if A1 and A2 both pass. Phase B is unchanged by this amendment.
- `select.R` writes `convergence.csv` with a `stratum` column: `n100` and `n300` for A1, `all` for the other phases, so that criterion 4 is counted at each A1 size. No gate or scoring code exists yet; the Task 4 scorer implements R21 as stated here.

## R23: phase A1 is launched, selected and repeated before phase A2

Reason: the 60 A1 fits take seconds each, while the 27 A2 fits take about 10 hours. Launching A1 first lets its selection and any single longer repeats finish within about an hour, instead of waiting for the A2 batch: the launcher refuses to run while another launcher holds its lock, so a combined launch would hold every A1 longer repeat back until all A2 fits were done. The fits, schedules and selection rule are unchanged; only the order of launches changes. The two launches are the two halves of `results/phase-a-queue.csv`, recorded by dry runs as `results/phase-a1-queue.csv` (60 initial fits) and `results/phase-a2-queue.csv` (27 longer fits).

- `select.R` now defaults its output to `STUDY/selection/<phases>`, for example `STUDY/selection/A1`, so the A1 and A2 selections are recorded separately and cannot overwrite each other. Both are recorded before any occupancy error is computed.
- The commands below replace the "Running phase A" section of the README. Run them in this order, each detached, and start each step only when the previous one has finished; for a launcher that means its `DONE` file in `STUDY/logs/launch-<timestamp>/` reports `exit_status: 0`. Substitute `REPO`, `STUDY` and `INPUTS` as given in the README.

Step 1, the A1 initial fits:

```sh
mkdir -p STUDY/logs
nohup Rscript REPO/dev/simstudy/occupancy-intercept-prior/launch.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A1 --sds=2,3,5 --max-procs=8 < /dev/null > STUDY/logs/launch-a1.out 2>&1 &
echo $! > STUDY/logs/launch-a1.pid
```

Step 2, the A1 selection, written to `STUDY/selection/A1`:

```sh
nohup Rscript REPO/dev/simstudy/occupancy-intercept-prior/select.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A1 --sds=2,3,5 --mode=plan --workers=4 < /dev/null > STUDY/logs/select-a1-plan.out 2>&1 &
```

Step 3, only if `STUDY/selection/A1/long-keys.txt` lists any fit, the single longer repeats:

```sh
nohup Rscript REPO/dev/simstudy/occupancy-intercept-prior/launch.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A1 --sds=2,3,5 --max-procs=8 --only=STUDY/selection/A1/long-keys.txt --schedule-override=long < /dev/null > STUDY/logs/launch-a1-long.out 2>&1 &
echo $! > STUDY/logs/launch-a1-long.pid
```

Step 4, the A1 final selection:

```sh
nohup Rscript REPO/dev/simstudy/occupancy-intercept-prior/select.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A1 --sds=2,3,5 --mode=final --workers=4 < /dev/null > STUDY/logs/select-a1-final.out 2>&1 &
```

Step 5, the A2 longer fits:

```sh
nohup Rscript REPO/dev/simstudy/occupancy-intercept-prior/launch.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A2 --sds=2,3,5 --max-procs=8 < /dev/null > STUDY/logs/launch-a2.out 2>&1 &
echo $! > STUDY/logs/launch-a2.pid
```

Step 6, the A2 selection, written to `STUDY/selection/A2`. A2 fits are never escalated, so the final mode follows directly:

```sh
nohup sh -c 'Rscript REPO/dev/simstudy/occupancy-intercept-prior/select.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A2 --sds=2,3,5 --mode=plan --workers=4 && Rscript REPO/dev/simstudy/occupancy-intercept-prior/select.R --repo=REPO --study=STUDY --inputs-root=INPUTS --phases=A2 --sds=2,3,5 --mode=final --workers=4' < /dev/null > STUDY/logs/select-a2.out 2>&1 &
```

No occupancy error is computed before both `STUDY/selection/A1/selected-fits.csv` and `STUDY/selection/A2/selected-fits.csv` exist.

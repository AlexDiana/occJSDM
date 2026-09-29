# Amendment 2 to the frozen protocol

Dated 29 September 2026. This amendment was made before any scoring of any arm in this study: no occupancy error, interval coverage or `B0` outcome of any fit, new arm or control, had been computed under this protocol, and no scoring code existed. At this date the phase A1 fits and their selection are complete (`STUDY/selection/A1`) and the phase A2 fits are running. Where this amendment differs from [README.md](README.md) or [AMENDMENT-1.md](AMENDMENT-1.md), this amendment applies; the rest of both files is unchanged.

## Status of the rulings

On 29 September 2026 Doug confirmed in chat the controller rulings R16 to R26, which had been listed to him as pending, and approved R27. The statements in AMENDMENT-1.md that R16 to R25 are pending his confirmation are therefore resolved: all rulings R1 to R27 are approved by Doug.

## R24: phase B is judged over its 10 generating communities

Reason: `design-qnear_K6-sites300-r` and `design-qfar_K6-sites300-r`, for the same replicate number r, share their generating occupancy probabilities; only detection and contamination differ. This is the dependence that R21 corrects in A1. Phase B therefore has 10 generating communities, each fitted at two contamination levels, not 20 independent communities.

- Criterion 1 in B (improvement). The gated group is the below-20% band only; the rare group is not gated in B (R22). For each replicate r and for each arm (the new arm and the control), the community value is the mean over the two contamination levels of the absolute signed error, (`abs(E(qnear))` + `abs(E(qfar))`) / 2, with E defined as in the README. Criterion 1 passes only if the new arm's community value is strictly smaller than the control's in at least 7 of the 10 communities (ceil of two thirds of 10) and the across-community mean of the new arm's community values is strictly smaller than that of the control's.
- No-harm criteria in B (criteria 2, 3 and 4 of the README) apply to each contamination level separately, `qnear` and `qfar`, each over its 10 communities: no band's across-community mean absolute signed error worse by more than 1 point, overall MAE not worse by more than 0.5 point, no band's across-community mean coverage lower by more than 0.03, and the number of flagged selected fits of the new arm no more than the control's. Each must hold at both levels.
- Phase B passes only if criterion 1, pooled as above, passes and all no-harm criteria pass at both contamination levels. This mirrors R21. It is recorded before any phase B fit exists.

## R27: a community without a valid selected fit

A community lacks a valid selected fit for an arm when no selected fit can be recorded for it, for example a fit that failed permanently or was quarantined and, under R25, is not rerun again.

- Such an arm fails the convergence criterion (criterion 4) for that sub-phase, whatever its flag counts, so the SD fails that sub-phase, and therefore fails phase A (or phase B). This extends R20 and R25.
- Criteria 1 to 3 are still computed and reported, descriptively, over the communities that have valid selected fits in both arms. For a pooled criterion 1 (R21 in A1, R24 in B) a generating community enters only if both of its fits are valid in both arms. These descriptive values cannot make the sub-phase pass.
- Interaction with `select.R`: its `check_first_fits` refuses to record any selection while a first fit of a requested SD is missing, so the frozen selection code cannot itself record such a case. If it ever occurs, it is handled by a documented manual record, not by editing frozen code. The controller writes `STUDY/selection/<phase>/manual-missing.csv` with the columns `phase`, `sd`, `key` and `reason`, one row per community lacking a valid selected fit, where the reason names the log or quarantine file that shows the cause, and reports the case in the study report. The selection of the other SDs is recorded by running `select.R` with `--sds` limited to them. The scorer treats every community listed in that file, and every community absent from `selected-fits.csv` for a scored SD, as lacking a valid selected fit.

## Correction: the detached launch wrapper

The wrapper used in AMENDMENT-1.md, `perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die'`, gives a step its own session only when it is started from a non-interactive shell. Under an interactive zsh, job control makes the wrapper a process-group leader, so `POSIX::setsid` fails with EPERM and returns -1; -1 is true in perl, so `or die` does not fire and the step runs on without its own session. The controller launches every step from a non-interactive shell, where the wrapper works as documented (the detached-wrapper test in `test-launch.R` checks this). The commands are unchanged.

## Scoring

The phase A scorer (Task 4: `analysis.R`, `summarise.R`, `verify.R`) implements the README's outcome and gate definitions with R18, R20, R21, R25 and R27 as stated in the README, AMENDMENT-1.md and this amendment. It is developed and validated without computing any new-arm outcome (controller ruling R28): its tests use hand-computed synthetic examples, and it is validated by scoring the control arm and reproducing the archived control results. The new arms are scored only after the scorer has passed review.

AMENDMENT-1.md ends its run procedure with the rule that no occupancy error is computed before both `STUDY/selection/A1/selected-fits.csv` and `STUDY/selection/A2/selected-fits.csv` exist. Under R28 the control arm is scored for this validation while the A2 fits are still running, so that rule is relaxed for the control arm only: its outcomes are already recorded in the archives, and no selection uses them, since every selection is computed by frozen code from convergence flags alone. For the new arms the rule is unchanged by this amendment.

# Protocol deviations

Departures from the frozen protocol ([README.md](README.md), [AMENDMENT-1.md](AMENDMENT-1.md), [AMENDMENT-2.md](AMENDMENT-2.md)) are recorded here as they are found, with what was done, when, by whom and what was seen. No outcome value is recorded in this file.

## Deviation 1: provisional look at phase A1 outcomes before the A2 selection

Recorded 29 September 2026.

What happened: from 12:21 to 12:28 BST on 29 September 2026, at Doug's request, the controller computed a provisional, unaudited summary of the phase A1 occupancy outcomes of the new arms and of the control from the selected A1 fits: signed error by band, MAE, interval coverage and `B0` bias. It was computed by a scratch script outside the repository, not by the study's scorer (`analysis.R`, `summarise.R`), and before `STUDY/selection/A2/selected-fits.csv` existed. This contravenes `AMENDMENT-1.md:83`, which says that no occupancy error is computed before both `STUDY/selection/A1/selected-fits.csv` and `STUDY/selection/A2/selected-fits.csv` exist.

Who saw the outcomes: the controller and Doug only. The scorer agent did not see them.

Facts that bound the effect of the look:

- AMENDMENT-2.md was committed before the look (`f02224e` at 11:58 and `04a34c5` at 12:00 BST).
- The scorer commit `f822020` (12:24 BST) was written by the blind scorer agent, which did not see the interim results.
- `select.R`, `flags.R` and the flag rules were frozen at `31d5e9a` and `077baff` before the look. The A2 selection they perform is mechanical, with no discretion, so the look cannot steer it.
- No protocol definition or gate criterion has been changed since the look.
- Every controller ruling made after the look is labelled post-look. R29 and R30 are the first: under R29, official new-arm scoring waits for the A2 selection, exactly as `AMENDMENT-1.md:83` requires, with no relaxation amendment; under R30, the gate's `GATE_TOL` of 1e-9 and the `manual-missing.csv` record of R27 are accepted subject to the task review.
- The interim numbers are not part of the study evidence. The evidence comes only from the official scoring in [SCORING.md](SCORING.md), after the A2 selection.

## Note: the basis of AMENDMENT-2's control-arm relaxation

Recorded 29 September 2026. AMENDMENT-2.md relaxes the rule of `AMENDMENT-1.md:83` for the control arm only, so that the scorer could be validated on the controls before the A2 selection exists, and cites R28 for it. R28 is a controller ruling, not among the rulings Doug approved (R1 to R27). The relaxation rests on controller ruling R28, because the control outcomes were already published in the archived studies, so scoring them early reveals nothing new. It is confirmed by controller ruling R31 (post-look), and is pending Doug's confirmation. AMENDMENT-2.md's body is unchanged.

## Note: post-look changes to the scorer's reporting

Recorded 29 September 2026. The task review of the scorer, after the look of Deviation 1, found gaps that fix round 1 corrects. Neither correction changes a gate criterion, threshold, band or group.

- A criterion evaluated over no community, which the scorer had reported as a pass for criteria 2 and 3, is now reported as NA, and a sub-phase with such a criterion and no failing one is UNDETERMINED rather than PASS. With the phase A truths every band and gated group has cells in every community, as the control tables show (`STUDY/summary/A1` and `A2`, `summary-means.csv`), and band membership depends only on the truth. The change can therefore matter only in the R27 case of a community without a valid selected fit.
- The manual record of R27 now includes `manual-selected.csv` and `supplement.R` (SCORING.md), so that the case which frozen `select.R` cannot record can be scored and reported as R27 requires.

The scorer agent made these changes without seeing any new-arm outcome.

## Note: phase B control scoring before the phase B selection

First written at 23:53 BST on 29 September 2026, before the scoring it describes, and reworded on 30 September 2026 to describe what then happened. At the controller's instruction for Task 5, the phase B scorer (`analysis-b.R`, with the phase B parts of `summarise.R` and `verify.R`) was validated by scoring the 20 phase B control fits (SD 1, the archived pr11 fits) from 00:05 to 00:35 BST on 30 September and auditing them from 00:35 to 00:47. Phase B new-arm fits were still running, and `STUDY/selection/B/selected-fits.csv` was written only at 01:11. The controls' occupancy estimates and errors and their `B0` and collection-intercept posterior means are recorded in the pr11 archive and were reproduced (`results/control-validation-B.csv`). Their occupancy interval coverage, their `B0` and collection-intercept interval coverage, and the posterior correlation of `B0` and the collection intercept were newly computed, because the archive records no intervals and no correlation. They are control outcomes only, and the phase B selection is computed by frozen code (`select.R`, `flags.R`) from convergence flags alone, so scoring the controls early could not steer it. This is controller ruling R35, made after scoring began on the same basis as R28 and R31, and not yet confirmed by Doug. No new-arm phase B outcome was computed before the phase B selection: `analysis-b.R` refuses new arms until it exists. The scorer agent did not load or inspect any phase B new-arm fit.

The same qualification applies to the phase A1 controls scored for validation before the A2 selection: the A1 archive scorer computes no intervals (README, "Community outcomes"), so their interval coverage was also newly computed, which the note above on AMENDMENT-2's relaxation does not say.

## Note: status of the controller rulings

Recorded 30 September 2026. Doug confirmed the post-look controller rulings R28 to R32 in chat on 29 September 2026. This resolves the statement above that the control-arm relaxation, which rests on R28 and R31, is pending his confirmation. Controller rulings R33 (phase B scoring in a separate `analysis-b.R`, so that `analysis.R` and the phase A score records stay unchanged), R34 (an UNDETERMINED phase B result counts as not adopted) and R35 (the phase B control scoring above) were made after scoring began and are not yet confirmed by Doug. None of them changes a gate criterion, threshold, band or group.

## Note: idle sleep during the A2 fits

Recorded 30 September 2026. This is a note, not a deviation from the protocol. The Mac idle-slept from about 16:39 to 17:49 BST on 29 September 2026 while the phase A2 fits were running. The fits paused and then resumed: none failed, was lost or was restarted (all 27 completed, `results/launch-A2-summary.csv`), and a pause cannot change the draws of a seeded fit, so the results are unaffected; it only lengthens the recorded elapsed times of the fits that were running. The controller's session was later restarted and resumed at 19:01; the detached launcher survived, and from 19:05 `caffeinate` kept the machine awake until the launcher exited.

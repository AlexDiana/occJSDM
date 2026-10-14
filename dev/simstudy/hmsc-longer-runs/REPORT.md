# Hmsc longer runs: results

Run on 14 October 2026 under the frozen protocol in [README.md](README.md) (commit `c34d411`); no amendment was needed. The 46 refits started at 01:53, when the mirror-frequency check finished, four at once, and the comparison finished at 04:03; the chain exited 0. A refit took 5.8 to 16.8 minutes (mean 9.8), well under the planning estimate, which had scaled the archive's own process times. Hmsc 3.3-7 from the check's library under R 4.6.1. The tables are under `results/`; the refits' draws and diagnostics (about 3.6 GB) stay in the archive `dev/simstudy/results/hmsc-longer-runs-20261013`.

## Answer

- **Running four times longer does not change what Lesson 6 reports for Hmsc.** In the lesson's main population ("All scored fits") no Hmsc cell moved by more than 0.2 percentage points of coverage, 0.003 points of bias or 0.004 points of RMSE at the five fixed environments, and none of the 18 prediction cells (300 test sites) moved by more than 0.003 points of RMSE or 0.005 of bias. The protocol's thresholds are 1 point of coverage and 0.5 points of error.
- **`decision.txt` nevertheless says "yes"**, because 7 of the 18 cells of the sensitivity population ("Passed diagnostics only") crossed a threshold. That is a change of composition, not of numbers: 33 of the 46 refits now pass the diagnostics, so those cells gained communities. For example the baseline at 300 sites goes from 5 to 8 passing communities and its coverage from 96.4% to 95.2% (MCSE 1.5), and the rare-species cell at 300 sites from 1 to 7 communities and from 98.0% to 94.6%. Every flagged cell is in that population; the cells' values for any one community barely moved.
- **Per fit, the saved and refit posteriors agree.** The largest change in any held-out site's probability was 0.91 percentage points, and in any estimate at the five environments 0.54 points; each job's coverage at the five environments is the same to two decimals (mean 0.908 before and after). The probit coefficients' posterior means moved by at most 0.055 of a saved posterior standard deviation, and every saved and refit 95% interval overlaps (46 jobs, three coefficients per species). The Hmsc trait rows of the trait jobs moved by at most 0.005. So the flagged fits' numbers are robust to the schedule, and the flags were about mixing, not about where the posterior is.
- **Diagnostics after the longer run.** 33 of 46 pass every check; 8 keep at least one failing coefficient quantity, 6 of them badly unmixed 30-species trait jobs (the largest coefficient Rhat falls from 1.293 to 1.170, and the smallest coefficient effective sample size among them is 16); 5 are flagged only in the residual species covariance block, which Lesson 6 does not report. Of the 14 badly unmixed jobs, 4 now pass. At 100 sites 5 of 6 pass; at 300 sites 28 of 40. The 30-species trait communities are where even 96,000 iterations per chain do not mix the coefficients and the covariance.

## What this means for Lesson 6

- Nothing in the main tables needs updating: the flagged fits stand, as the TODO item provided for when nothing changes.
- The lesson's counts and its passed-only sensitivity summaries would change if the refits replaced the saved fits: 53 Hmsc fits would be flagged instead of 86. Whether to do that is Doug's decision. It is a separate task, because the lesson's bundle is built by the calibration exporter from the extension archive and verified against it; this check wrote nothing into that archive or the bundle. It would mean a second archive or a documented replacement of those 46 `result.rds` files, a new export and both verifiers.
- The TODO item's premise was corrected in README.md: Lesson 6 reports no Hmsc coefficient recovery or coefficient coverage, so "coefficient errors and interval coverage" for Hmsc have no truth to compare against; the coefficient comparison above is saved against refit.

## Limits

- R 4.6.1 against the archive's R 4.5.0, the same Hmsc version but unpinned dependency versions, and a new seed: the refits are fresh runs, not replays, so no bitwise agreement is claimed or expected.
- The chain is thinned by four to keep 16,000 retained draws per chain rather than the TODO item's 64,000; the diagnostics' effective-sample-size threshold applies to the same retained count as the saved attempt.
- The decision rule flags any cell that crosses the thresholds, including by composition; the report says which kind each change is.

## Files

- `results/jobs.csv`: the 46 jobs and why each was selected.
- `results/comparison-by-job.csv`: one row per job (diagnostics before and after, 300-site errors, five-environment coverage and width, coefficient and trait changes, timings).
- `results/lesson6-probability-cells.csv` and `results/lesson6-prediction-cells.csv`: Lesson 6's Hmsc cells before and after substitution, with the change flags; `results/summary.csv` and `results/decision.txt`.
- Archive: `jobs/<job>/attempt-3-parameters.rds`, `attempt-3-diagnostics.csv`, `result.rds`; `logs/chain.log`; `library/`; `CHAIN-DONE`.

# Hmsc longer runs: do they change Lesson 6's Hmsc results?

This is the protocol for the TODO.md item "Check whether longer Hmsc runs change Lesson 6's coefficient results". The commit that adds this file is the freeze: it was written after the scripts and tests below were finished and before any refit was launched. Nothing in it may be edited after the launch except by a dated amendment file in this directory (`AMENDMENT-YYYY-MM-DD.md`) that says what changed and why. Doug asked for the check on 13 October 2026, to run after the mirror-frequency check finishes; it was built the same day without a Codex review, as a low-priority check.

`ROOT` below is the Lesson 6 extension's run archive, `dev/simstudy/results/lesson-4-extension-20260923` (git-ignored, read only here); `ARCHIVE` is this check's archive, `dev/simstudy/results/hmsc-longer-runs-20261013` (git-ignored). No production code, lesson, teaching bundle or earlier archive is changed.

## The question, and a correction to its premise

Of Lesson 6's 180 Hmsc fits, 86 stayed flagged after their one longer attempt (8,000 warm-up and 16,000 retained iterations per chain, four chains). Their failures are mostly in the residual species covariance, which Lesson 6 does not report, but 46 of them also fail at least one environmental coefficient (rank-normalised Rhat above 1.01 or an effective sample size below 400). Sixteen of the 86 are badly unmixed (more than 20 failing quantities), 14 of them among the 46; the largest coefficient Rhat is 1.293. The TODO item asks whether running those 46 four times longer changes "Lesson 6's coefficient results".

The premise needs one correction, found while building this check: Lesson 6 reports no Hmsc coefficient recovery or coefficient interval coverage. The [calibration protocol](../jsdm-package-comparison/extension/CALIBRATION.md) excludes Hmsc from the coefficient target because its probit coefficients have no logit truth, and the lesson labels its coefficient cells "Different link". What Lesson 6 reports for Hmsc is: prediction error at the 300 held-out sites, marginal-probability bias, RMSE and 95% interval coverage at five fixed environments, trait-effect rows for the trait experiment, and the counts of passed and flagged fits. So the question this check answers is whether a four-times-longer run changes those reported Hmsc numbers, and whether the coefficient draws themselves move; the latter is reported as a change between the saved and refit posteriors, with no truth.

## Which fits

`results/jobs.csv`, written once by `select-jobs.R` from the archive's `status.csv` and each flagged Hmsc job's `attempt-2-diagnostics.csv`: the Hmsc jobs with `fit_ok` true and `diagnostic_pass` false whose longer attempt has at least one failing row in the `coefficient` block. There are 46 (6 at 100 sites, 40 at 300; by scenario 3 baseline, 5 correlated, 12 curved, 6 rare, 20 traits). The file records each job's failing-quantity counts by block, its largest coefficient Rhat, whether it is badly unmixed (more than 20 failing quantities; 14 of the 46) and the saved run time of its two attempts. `test-refit.R` checks that the file reproduces the rule.

## What is fitted

Each job is refitted with `fit_hmsc()` in `common.R`, which copies the Hmsc statements of the frozen worker (`ROOT/worker-final/extension/fit-job.R`, md5 `67b62869071142a9426fe4467f342d53`, checked before any fit) verbatim: the same input (md5-checked against the manifest, test predictors dropped before fitting), `XScale = FALSE`, the probit link, one site-level random level with exactly two latent factors (`nfMin = nfMax = 2`), and the trait formula where traits were supplied. Only the sampler settings and the seed differ.

- **Schedule**: `transient = 32000`, `samples = 16000`, `thin = 4`, four chains, `nParallel = 1`: 32,000 warm-up and 64,000 post-warm-up iterations per chain, four times the longer attempt's 8,000 and 16,000, with the same 16,000 retained draws per chain. The TODO item said "64,000 retained draws"; retaining every draw would quadruple the saved posterior and the memory each worker holds (the saved attempt holds about 78 MB of draws per 300-site trait job), so the retained count is kept and the chain is thinned by four instead. The diagnostics' effective-sample-size threshold of 400 then applies to the same number of retained draws as before.
- **Seed**: the frozen worker's rule, `26100000 + 100 * job number + attempt`, with attempt 3; attempts 1 and 2 are the archive's.
- **Software**: Hmsc 3.3-7, the version the archive records, installed on 13 October 2026 from CRAN into `ARCHIVE/library` together with any dependency the system library lacked, under system R 4.6.1 (the archive's fits ran under R 4.5.0); `common.R` refuses any other Hmsc version or location. Each process pins the BLAS thread variables to 1. Dependency versions other than Hmsc's are recorded in each result's session information, not pinned.
- **Extraction and diagnostics**: the extension's own `extract_bayesian()` and `parameter_diagnostics()` from `dev/simstudy/jsdm-package-comparison/` (`pilot-math.R`, `extension/math.R`, `extension/study.R`, `extension/calibration-math.R`, refused unless their md5s are the ones the calibration bundle records), so a refit is judged by exactly the rule that flagged the saved fit: every coefficient, residual covariance, new-site probability and trait quantity must have rank-normalised Rhat at most 1.01 and bulk and tail effective sample sizes of at least 400, with the integration check.
- **Output**: `ARCHIVE/jobs/<job>/attempt-3-parameters.rds`, `attempt-3-diagnostics.csv` and `result.rds` (settings, seed, diagnostics summary, warnings, timings, md5s, session), written under a per-job lock. A job already saved with identical settings is reported as resumed; one with different settings is never overwritten. `refit.R` never opens `ROOT/truth/`.

## Comparison and decision rules

Fixed before any refit. `compare.R` runs after every refit has finished and is the first script here to read truth.

- **Per job**, saved longer attempt against refit: whether the refit passes the diagnostics and its failing-quantity counts by block; the 300-test-site bias, mean absolute error and RMSE (percentage points, the saved values recomputed and checked against the archived score) and the largest change in any site-by-species probability; at the five fixed environments, interval coverage, mean interval width, bias and RMSE, and the largest change in any estimate (the saved values recomputed and checked against the calibration bundle's rows); the probit coefficients on the fitted scale, as the largest and root-mean-square change in posterior means, that change in units of the saved posterior standard deviation, the share of elements whose saved and refit 95% intervals overlap, and the mean interval widths; and for jobs with traits supplied, the four Hmsc trait rows as the archived scorer computed them.
- **Lesson 6's Hmsc cells**: the calibration bundle's Hmsc rows (marginal probability at the five environments) and per-job predictions are copied with the 46 jobs' values replaced by the refits' and their diagnostic status updated, and the lesson's grouped summaries recomputed with the extension's `calibration_group_summary()` for both populations ("All scored fits" and "Passed diagnostics only"), and cell means of the 300-site bias, MAE and RMSE with communities weighted equally. The saved-side recomputation must match the bundle's summary before any comparison.
- **A Lesson 6 cell counts as changed** when its coverage moves by at least 1 percentage point, or its bias or RMSE (at the five environments, or at the 300 sites) by at least 0.5 points of probability. These thresholds are at the resolution the lesson prints (coverage to one decimal, errors to one decimal) and below the coverage MCSE of every baseline cell. `decision.txt` says `lesson6_hmsc_numbers_change: yes` if any cell is changed, else `no`, with the counts.
- **What follows**: if no cell changes, the flagged fits stand and the TODO item closes with this record. If a cell changes, the report says which, and updating Lesson 6 is a separate task for Doug to decide, because the lesson's bundle is built by the calibration exporter from the archive and this check writes nothing there.
- Every per-job quantity is reported whatever the cell decision, including how many refits now pass and whether the coefficient failures resolved, so the result also says whether a four-times-longer run is enough for these fits.

## Cost and launch

Planning estimate: the 46 jobs' two attempts (30,000 iterations per chain in all) took 5.9 hours of process time together (300-site jobs 7.1 minutes on average, 100-site jobs 12.1); the refit's 96,000 iterations per chain is 3.2 times that, about 19 hours of process time, or roughly five to seven hours with four at once. Each worker holds one job; four at once took about the same memory as the original study's four workers. Disk: about 78 MB of draws per 300-site trait job, under 4 GB in all.

The check starts automatically when the mirror-frequency check writes its `CHAIN-DONE` file, through `wait-and-run.sh`, which then runs the launcher and the comparison and writes `ARCHIVE/CHAIN-DONE` with its exit status. Keep the Mac on mains power with the lid open. Run from a non-interactive shell. The first command is a dry run that starts nothing.

```sh
REPO=/Users/douglasyu/src/occJSDM/.worktrees/mirror-frequency
ROOT=/Users/douglasyu/src/occJSDM/dev/simstudy/results/lesson-4-extension-20260923
STUDY=/Users/douglasyu/src/occJSDM/dev/simstudy/results/hmsc-longer-runs-20261013
MIRROR=/Users/douglasyu/src/occJSDM/dev/simstudy/results/mirror-frequency-20261013
S=$REPO/dev/simstudy/hmsc-longer-runs
Rscript "$S/launch.R" --repo="$REPO" --root="$ROOT" --study="$STUDY" --max-procs=4 --dry-run
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh -c "caffeinate -i sh $S/wait-and-run.sh $REPO $ROOT $STUDY $MIRROR 4" < /dev/null > "$STUDY/logs/chain.out" 2>&1 &
```

`ARCHIVE/logs/chain.log` is the chain's log; each launch writes `ARCHIVE/logs/launch-<timestamp>/` with `progress.csv`, one `.log` and `.exit` file per job, and `DONE`. A failed or killed job is relaunched with `--jobs=N` (completed jobs resume), after removing any lock a killed process left, once no `refit.R` process holds it; then run `compare.R` by hand.

## Files

- `common.R`: constants, the md5-checked code loader, the Hmsc loader and `fit_hmsc()`.
- `select-jobs.R`: writes `results/jobs.csv` once.
- `refit.R`: one job per process; `launch.R`: the launcher; `compare.R`: the comparison; `wait-and-run.sh`: the detached chain.
- `test-refit.R`: the tests (the schedule and seeds; the job selection against the archive; the code and worker md5s; a tiny four-chain Hmsc fit through the frozen statements that is reproducible from its seed, differs under another seed, extracts, passes through the diagnostics and the comparison helpers).
- `results/jobs.csv`: the selection; after the comparison, `comparison-by-job.csv`, `lesson6-probability-cells.csv`, `lesson6-prediction-cells.csv`, `summary.csv` and `decision.txt`.

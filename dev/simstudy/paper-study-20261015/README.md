# Paper simulation study on the combined branch of PRs 31 and 32, 15 October 2026

The full paper simulation grid at R = 300, run on the combined, unmerged code of [PR 31](https://github.com/AlexDiana/occJSDM/pull/31) (parallel chains with per-chain seeds) and [PR 32](https://github.com/AlexDiana/occJSDM/pull/32) (symmetric Gram matrices). This is the rerun that the *Paper simulation study* item in `TODO.md` asks for, done ahead of the merges. If either PR changes in review, the code that produced these results is not the code that merged, and the run has to be repeated.

## Code

Branch `claude/paper-study-combined`, pushed to origin, in the worktree `.worktrees/paper-study-combined`. It is main at `d156492` with the two PR heads merged in order (`7a47d1b`, then `375b387`); the only conflict was NEWS.md, resolved by keeping both entries. The full local test suite passed on it: 1386 passed, 0 failed, 5 skipped. This run directory and the timing pre-study's live on the same branch, committed before launch, so the revision in `provenance.txt` and in the results file includes them. The `.gitignore` here keeps the log, pid files, provenance and results out of the index while the run is in flight, so the runner's dirty-tree check stays clean. Do not commit or edit anything in the worktree until the run finishes: a commit moves `HEAD` under the run, and an edit to the runner risks a torn read.

## Run

The production ten cells, R = 300 per cell (3,000 fits), 5 replicate processes, one thread per fit, the runner's default 1000 burn-in and 1000 sampling iterations with 2 chains. The runner is `dev/simstudy/run_study.R` from the worktree, unchanged by either PR, pointed at the worktree with `--pkg`. There is no time cap.

Launched detached from this directory with:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh > launch.log 2>&1 < /dev/null &
```

Monitor with `tail -f run.log`; the progress line's cpu/wall figure should stay near 4x, and a fall towards zero means the machine slept or throttled. To stop it, kill the R master named in `rscript.pid` and then its five PSOCK workers, which are left running their current fit when the master dies. Resume an interrupted run with `RESUME=1` set in front of the same command; the checkpoint under `results/` holds every completed chunk. A closed lid on battery still sleeps the machine, which `caffeinate` cannot prevent.

## First attempt, stopped

The first launch, at 09:38, used a run directory created by mistake in the main checkout rather than on this branch. It was stopped after three minutes, before any chunk completed, and the directory was moved here. Its log and provenance are kept, uncommitted, in `aborted-0938/`.

## Expected duration

The [timing pre-study](../paper-study-prestudy-20261015/README.md) earlier the same morning did 70 fits in 10.4 minutes on this code and machine, on mains power with Low Power Mode off: 6.7 fits per minute. At that rate 3,000 fits take about 7.4 hours. Low Power Mode on battery halved the sampler speed when it was measured on 14 October 2026, so a run on battery with it on takes about twice as long. The power state at launch is in `provenance.txt`.

## Results

The run finished at 17:28 with exit status 0: all 3,000 fits, 300 in each of the ten cells, and no failed replicate. Results are in `results/`: the full results file `simstudy-20261015-172834.rds`, the coverage summary `simstudy-20261015-172834.csv` and the occupancy-recovery summary `simstudy-occupancy-20261015-172834.csv`. The checkpoint is not committed; it duplicates the rows of the results file.

### Timeline and worker count

The run took 7 hours 45 minutes of wall time, from 09:43 to 17:28, in three invocations joined by `--resume`. Each stop came right after a chunk was checkpointed, so no finished fit was lost or repeated.

-   09:43 to 11:39, 5 workers: 620 fits. The laptop was on battery with Low Power Mode on until about 10:38, which held it to about 3.2 fits per minute; with Low Power Mode off it ran at about 7.3 fits per minute, and it was on mains by 11:35.
-   11:39 to 13:07, 7 workers, at Doug's request: 644 fits, also about 7.3 fits per minute overall. It started at about 8.6 fits per minute and fell to about 5.9 by the end, while the CPU time per fit rose from about 0.70 to about 0.95 CPU-minutes, against about 0.55 with 5 workers. That pattern fits a fanless chip slowing its clock under sustained heat, with the two extra workers on efficiency cores; macOS does not report throttling on this chip without admin rights, so it is inferred. A Time Machine backup and a Teams call also competed for part of that stretch.
-   13:07 to 17:28, 5 workers again: 1,736 fits at about 6.6 fits per minute, recovering to about 7 to 7.7 once the machine was quiet.

On this laptop, 5 workers is the better setting for a run of hours; 6 and 4 were not tried. The launchers for the two resumed invocations are `resume-cores7.sh` and `resume-cores5.sh`; they ran from outside the worktree to keep its tree clean, and are copied here for the record. `provenance.txt` and `run.log` record every stop, resume and power change.

### Provenance

The results file holds all 3,000 replicates, but its provenance block describes only the last invocation: revision `de4389c`, clean tree at start and at write, `HEAD` unmoved, 5 workers, one thread per fit, started 13:07, elapsed 261.5 minutes. The 7-worker invocation started with one untracked file in the tree, an empty log, so its own dirty flag was set; that invocation wrote no results file, and the cause is recorded in `provenance.txt`. Each replicate seeds itself from its scenario and replicate number before simulating and fitting, so its result does not depend on which worker ran it or how many workers there were. R used Apple's Accelerate BLAS throughout, so a rerun reproduces these numbers to rounding error, not bit for bit.

### Coverage at a glance

This is a summary for orientation, not the assessment the *Paper simulation study* item in `TODO.md` asks for. The runner gives the coverage standard error at R = 300 as about 1.3 percentage points for one block in one cell. The comparison is with the 29 July 2026 R = 100 run (`dev/simstudy/results/simstudy-20260729-143756.csv`, uncommitted in the main checkout), which ran on much earlier code; differences larger than a few points are real but not yet explained.

-   `beta_theta`: 0.960 to 0.967 in the eight cells with it apart from `low_information` (0.995), against 0.71 to 0.77 in July.
-   `B`: 0.875 (`spatial_isolated`) to 0.927 (`base`), and 0.799 in `low_information`; below nominal in every cell, and lower than July on average (0.895 against 0.922).
-   `B0`: 0.929 to 0.942 in most cells, but 0.868 in `binary` (0.962 in July) and 0.845 in `low_information`.
-   `p`: 0.908 to 0.921, unchanged from July on average; 0.103 in `low_information`.
-   `theta0`: 0.974 to 0.983, and 0.604 in `low_information`.
-   `q` and `G`: within about two standard errors of nominal, except in `low_information` and for `G` in `species_20` (0.912).
-   `resid_cor`: 0.751 to 0.766, and 0.979 in `d_underfit`, as in July. `dev/simstudy/PLAN.md` section 17 explains why this block's coverage says nothing about the estimator.
-   `sigma_b`: 1.000 everywhere; the runner notes it is prior-dominated.
-   Occupancy recovery: the mean correlation between true and estimated occupancy probabilities is 0.58 to 0.68 and the mean AUC 0.77 to 0.80 in eight of the nine cells that report it; `low_information` is far weaker (0.13 and 0.56). The `binary` cell has no occupancy summary.

`low_information` is the small, low-detection cell: 40 sites rather than 100, with `p` drawn from 0.1 to 0.3 and the baseline `theta` from 0.05 to 0.2. `dev/simstudy/PLAN.md` already recorded it as compromised across the board in the July run.

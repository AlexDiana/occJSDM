# Paper study timing pre-study, 15 October 2026

A short run of the paper simulation grid to estimate how long the full study would take on the combined code of [PR 31](https://github.com/AlexDiana/occJSDM/pull/31) (parallel chains with per-chain seeds) and [PR 32](https://github.com/AlexDiana/occJSDM/pull/32) (symmetric Gram matrices), before either is merged. It is a timing run only. Its coverage numbers are at R = 7 per cell and say nothing about calibration.

## Code

Branch `claude/paper-study-combined`, pushed to origin, in the worktree `.worktrees/paper-study-combined`. It is main at `d156492` with the two PR heads merged in order (`7a47d1b`, then `375b387`); the only conflict was NEWS.md, resolved by keeping both entries. The full local test suite passed on it before launch: 1386 passed, 0 failed, 5 skipped. The exact revision is in `provenance.txt`.

## Run

The production ten cells, R = 7 per cell (70 fits), 5 replicate processes, one thread per fit, the runner's default 1000 burn-in and 1000 sampling iterations with 2 chains. The runner is `dev/simstudy/run_study.R` from the worktree, unchanged by either PR, pointed at the worktree with `--pkg`. The run is capped at 30 minutes by `timeout`; exit status 124 in `provenance.txt` means the cap fired, in which case the checkpoint under `results/` holds the completed chunks and `run.log` holds the per-chunk progress lines.

Launched detached from the run directory with:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh > launch.log 2>&1 < /dev/null &
```

Monitor with `tail -f run.log`. Machine on mains power, Low Power Mode off (recorded in `provenance.txt`).

This directory was created under `dev/simstudy/` in the main checkout and moved into the worktree after the run, so `provenance.txt` names the old path; `launch.sh` now names the new one. The `.gitignore` here keeps the checkpoint and the launcher's pid files out of the index; `run.log`, `provenance.txt`, the results tables and the results file that `estimate.R` reads are committed with `git add -f`; the checkpoint is not.

## Sizing basis

The 29 July 2026 full grid (`dev/simstudy/results/full-run-20260729-095230.log`) did 1000 fits in 285 minutes on 5 cores, about 3.5 fits per minute, on the pre-Accelerate build. At that rate 70 fits take about 20 minutes. The spatial cell should be faster now under Accelerate and PR 32; the other nine cells are not expected to change much.

## Estimate

The run finished in 10.4 minutes: 70 of 70 fits, 0 failures, at 6.7 fits per minute on 5 cores with cpu/wall between 3.7x and 4.2x throughout (`run.log`). That is about 1.9 times the 29 July rate of 3.5 fits per minute. The results file is `results/simstudy-20261015-082406.rds`, and `estimate.R` reproduces the figures below from it.

Scaling linearly, since the runner mixes all ten cells across the workers and its throughput was flat across the whole run, the full grid on this machine at these settings would take about:

-   R = 100: 1,000 fits, about 2.5 hours.
-   R = 200: 2,000 fits, about 5 hours.
-   R = 300: 3,000 fits, about 7.5 hours.
-   R = 500: 5,000 fits, about 12.5 hours.

These are slightly conservative. At 70 jobs the runner used chunks of 5, so each chunk waited for its slowest fit; the full grid uses chunks of 20 with load balancing inside each, which wastes less of the tail. They assume the same conditions as this run: 5 replicate processes, one thread per fit, mains power, Low Power Mode off, Accelerate BLAS under R 4.6.1, and nothing else heavy running. A sleeping or throttled laptop stretches the wall time without changing the compute; the cpu/wall figure in the progress line is the way to tell.

To run the full study on the same code, launch from this directory or a sibling run directory with the same command as `launch.sh` but `--R=` set to the chosen replicate count, no `timeout` wrapper, and a fresh `--out`.

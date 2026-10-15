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

Filled in after the run.

# Stage 2: spatial representation comparison

This study fits five existing spatial representations on 40 of the recoverability pilot's communities and measures, against the pilot's oracle, how much held-out headroom each recovers under each truth. The design is `dev/superpowers/specs/2026-10-11-spatial-stage2-comparison-design.md` and the plan is `dev/superpowers/plans/2026-10-11-spatial-stage2-comparison.md`. **Status: shelved on 12 October 2026, before any comparison fit, until server time is available.** `SHELVED.md` explains what was built, what was measured, why it cannot run on the laptop, and how to resume; `STATUS.md` is the short status. The compact evidence is in `results/`. There is no `REPORT.md` because there are no results.

Everything runs through `run-r`, the pinned R 4.5.0 launcher, never system `Rscript`. Hmsc 3.3.7 and gllvm 2.0.15 live in the study library `dev/simstudy/results/spatial-stage2-20261011/library`, and occJSDM 0.1.0 is the first comparison's pinned build. Nothing in `R/`, `src/`, priors, defaults or lessons changes. The pilot's inputs and job files are read-only.

## Reproducing the run

From the worktree root, with the research volume mounted:

```sh
P=dev/simstudy/spatial-stage2-comparison
STUDY=/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-stage2-20261011
PILOT=/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010
JOBS="/Volumes/occJSDM Research/occJSDM/spatial-stage2-20261011/jobs"
for t in runtime data adapters score jobs run verify plot microtiming; do $P/run-r -f $P/test-$t.R || exit 1; done
$P/run-r -f $P/run.R --args --mode=prepare --study=$STUDY --pilot=$PILOT --jobs="$JOBS"
$P/run-r -f $P/run.R --args --mode=microtiming --study=$STUDY
nohup caffeinate -i $P/run-r -f $P/run.R --args --mode=timing --study=$STUDY --step=1 > $STUDY/logs/timing-step1.log 2>&1 &
$P/run-r -f $P/run.R --args --mode=freeze --study=$STUDY --step=1
nohup caffeinate -i $P/run-r -f $P/run.R --args --mode=fit --study=$STUDY --workers=8 > $STUDY/logs/fit.log 2>&1 &
$P/run-r -f $P/run.R --args --mode=score --study=$STUDY
$P/run-r -f $P/verify.R --args --study=$STUDY
$P/run-r -f $P/plot.R --args --study=$STUDY
```

Only the commands down to the micro-timing were run (the study was shelved at that point, see `SHELVED.md`). The timing ladder steps down through `--step=2`, `--step=3` and so on until a step passes the 48-hour rule, and `freeze` takes the number of the step that passed. `fit` can be relaunched after any interruption: completed attempts are reused and never refitted.

## Files

- `config.R` is the frozen protocol. `data.R` bridges the pilot's communities to the adapter contract and loads the oracle reference.
- `upstream/` holds the first comparison's adapters, basis, diagnostics, selection code and certification tests, copied with recorded hashes in `sources.csv`. `modified` marks the files edited for stage 2.
- `score.R` computes the headroom share, pattern, slope, aggregation, paired comparisons and the development reading.
- `jobs.R` drives one job, one dataset and arm: attempts, selection, prediction and storage. `run.R` drives the study: prepare, micro-timing, timing, freeze, fit and score.
- `probe-gllvm.R`, `probe-gllvm.sh` and `probe-costs.R` are small watched probes of cost and memory (not study attempts) whose output is in `results/`. `probe-accelerate-gp.R`, `probe-accelerate-hmsc.R` and `probe-accelerate-projection.R` are the 14 October 2026 re-check of the gp and Hmsc costs under Apple's Accelerate BLAS; they run under system R 4.6.1 rather than `run-r`, and their output is `results/probe-accelerate.txt`. `SHELVED.md` is the report on the shelved study, with the re-check as its section 9.
- `verify.R` recomputes the intervals, the reading, the statuses and the study bookkeeping independently. `plot.R` copies the compact results and draws the figures.
- `test-*.R` are the tests, one per file above. `results/` holds the committed compact evidence and `STATUS.md` the short status.

# Handover: running the spatial recoverability pilot

Written 10 October 2026 by the session that designed and built the study. Doug decided that the design and critique were the expensive work and that the run itself should be carried out by a cheaper model. This file is the complete state at handover and the exact remaining steps. Read the spec (`dev/superpowers/specs/2026-10-10-spatial-recoverability-pilot-design.md`) and the plan (`dev/superpowers/plans/2026-10-10-spatial-recoverability-pilot.md`) first; the plan's Tasks 1 to 6 and the scripts of Task 8 are complete and tested, and Task 7 is in progress.

## State at handover

- Branch `claude/spatial-recoverability-pilot` in worktree `~/src/occJSDM/.worktrees/spatial-recoverability-pilot`, nine commits ahead of `main`, working tree clean after this file is committed. Run every command from that worktree root.
- All six test files pass under the pinned launcher: `test-seeds.R`, `test-kernel.R`, `test-generator.R`, `test-oracle.R`, `test-score.R`, `test-run.R`. The generator and score tests take two to four minutes each.
- The `prepare` step is running detached (started 14:27 BST, PID 18073, eight workers) and writes to `~/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010/`. Its log is `logs/prepare.log` in that directory and ends with `Completed mode prepare` on success. At handover 40 of 120 inputs existed; the 900-site communities are the slow ones, so expect another 20 to 40 minutes. The process is not attached to any session and will finish on its own.
- The job directory on the research volume, `/Volumes/occJSDM Research/occJSDM/spatial-recoverability-20261010/jobs`, is empty. The volume must stay mounted for the whole run.
- An earlier `prepare` attempt failed because clipped clustered sites duplicated a corner point; the generator was fixed (commit `8b75843`), the half-built directory was deleted, and `prepare` was relaunched from scratch. No fit has ever been run.
- The thresholds and budget are frozen in `config.rds`, written before any community was generated. Do not edit them, and do not re-run `prepare` on this directory; it refuses to.

## Remaining steps, in order

Set the shell variables used below from the worktree root:

```sh
P=dev/simstudy/spatial-recoverability-pilot
STUDY=/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-recoverability-20261010
```

1. Wait for `prepare`. Check with `tail -2 $STUDY/logs/prepare.log` and `ls $STUDY/inputs | wc -l` (target 120). Success means `input-manifest.csv` exists with 120 rows and the log ends with the completion line. If the log instead shows `Execution halted`, read the error, fix the cause in the generator with a test, delete the study directory and the empty jobs directory, and relaunch `prepare` exactly as the README shows.
2. Run the timing probe: `$P/run-r -f $P/run.R --args --mode=timing --study=$STUDY`. It fits one 900-site truth 2 job (twelve per-species oracles) and one 900-site, range 0.15 truth 1 joint job at the full initial budget, including lattice prediction, and writes `timing.json` with an upper-bound CPU-hours projection. These two jobs are kept and reused by `fit`.
3. Apply the fixed rule from the spec: proceed only if `projected_cpu_hours_upper / 8` is under 12 wall-clock hours. If it is not, stop and report the numbers to Doug. Do not reduce the grid, the budget or the number of communities yourself. A crude earlier measurement put a single 900-site per-species oracle at about 8 seconds for 7,000 iterations, so the rule is expected to pass comfortably.
4. Launch truth 2 detached: `nohup $P/run-r -f $P/run.R --args --mode=fit --study=$STUDY --truth=t2 --workers=8 > $STUDY/logs/fit-t2.log 2>&1 &`. Wait until the log ends with `Completed mode fit`. Then launch truth 1 the same way with `--truth=t1` into `fit-t1.log`. Progress: `ls $STUDY/status/*.json | wc -l` reaches 60 after truth 2 and 240 after truth 1; `ls $STUDY/status/*.error` must stay empty. An interrupted `fit` can simply be relaunched; existing job files are reused, never refitted.
5. Score: `$P/run-r -f $P/run.R --args --mode=score --study=$STUDY`. Expect `scores/species.csv` with 2,880 rows, `gate.csv` with 18 rows (nine truth 1 cells, three classes by three site counts for truth 2) and `map.json` naming one of `both`, `truth1_only`, `truth2_only`, `stop`.
6. Verify independently: `$P/run-r -f $P/verify.R --args --study=$STUDY`. It recomputes every gate decision and the map with separate code, checks input hashes, seeds, status files and conditioning tolerances, writes `results/verification.csv`, and exits non-zero on any failure. A failure here is a defect to diagnose, not a result to report.
7. Plot and copy compact results: `$P/run-r -f $P/plot.R --args --study=$STUDY`. Open `results/recoverability-map.png` and check for clipped or overlapping labels.
8. Write `REPORT.md` following the plan's Task 8 Step 4 outline. Every number comes from `results/gate.csv`, `results/aggregate.csv`, `results/jobs.csv` and `results/timing.json`, never from memory. The report must state, verbatim from the spec, that the single-season footprint is not evidence of dispersal as a mechanism, and that truth 1's advantage over truth 2 follows from twelve species pooling onto two fields. It must also state the posterior mean, not median, is the field summary, and list the factorisation routes used by cell from `results/factors.csv` and `results/conditionals.csv`.
9. Add one line to the spec's "Revision record" naming the study directory, the timing numbers and the launch and completion times. Add one bullet to `TODO.md` in the spatial section pointing at the report. Commit.
10. Request a whole-branch review against `main` (the plan's Task 8 Step 7), fix findings, and hand Doug the report and the map. Doug reviews rendered output, not diffs.

## Conventions that bind the run

- Always use `$P/run-r`, never system `Rscript` (system R is 4.6.1; the study is pinned to 4.5.0 through the private R home of the first comparison).
- Markdown in this repo: one line per paragraph, no pipe tables, no em-dashes, escape `-\>` and `\~` in prose and never inside code. Re-read a markdown file before editing it, since Doug may have it open in RStudio.
- Commit messages end with the attribution line used by the earlier commits on this branch.
- Long jobs are launched with `nohup` so that an app restart does not kill them. Workers are single-threaded through the launcher's environment variables.
- Nothing in the occJSDM package, its priors, defaults or lessons changes.

## Expected result and how to read it

The spec expects truth 1 (common gradients) to pass at a smaller site count than truth 2 (per-species dispersal fields), because the joint oracle pools twelve species' observations onto two fields. The map-reading rule in the spec's section 3 covers all four branches. If truth 2 passes while truth 1 fails, inspect the truth 1 oracle for defects before writing anything else. If nothing passes even at 900 sites and range 0.15, that is the finding, and no stage 2 is proposed.

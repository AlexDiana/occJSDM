# Spatial comparison handover to Claude

Continue Doug's authorised spatial JSDM comparison and its Lesson 5 and Lesson 7 teaching deliverables. The study runs in detached, supervised processes on this Mac, independently of the originating chat. Changing assistants does not require restarting the fits. This handover was prepared on 8 October 2026 at approximately 19:46 BST for a possible transfer on Saturday 10 October. Read the live records below when taking over; the dated status is not a forecast of Saturday's state.

## Start here

Read the repository's `AGENTS.md`, the current specification and implementation plan, and the most recent entries at the bottom of the execution ledger. Its opening task checklist is historical and is less current than its dated entries.

| Purpose | Location |
|---|---|
| Working checkout | `/Users/douglasyu/src/occJSDM/.worktrees/spatial-comparison-plan` |
| Branch | `codex/spatial-comparison-plan` |
| Specification | `dev/superpowers/specs/2026-10-08-spatial-package-comparison-design.md` in that checkout |
| Implementation plan | `dev/superpowers/plans/2026-10-08-spatial-package-comparison.md` |
| Execution ledger | `.superpowers/sdd/2026-10-08-spatial-package-comparison/progress.md` |
| Persistent run, called RUN below | `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-package-comparison-20261008` |
| Large fit archive | `/Volumes/occJSDM Research/occJSDM/spatial-package-comparison-20261008/jobs` |

RUN is outside the disposable worktree. Its `jobs` symlink points to the external research volume. Keep that volume mounted and the Mac on mains power with the lid open. The supervisor uses `caffeinate`; it does not make shutdown or unplugging the archive safe. The separate Time Machine backup volume is not part of this study.

## Live execution

As of preparation, the stage was `finish_pilot`, the two extended occJSDM GP fits had completed in approximately 110 minutes each, and all 16 Hmsc attempts had been recovered successfully without refitting. Predictions and remaining pilot extensions were proceeding. No scientific recovery scores had been opened by the assistant. The supervisor's stderr was empty.

Use these files as the current source of execution status:

- `RUN/continuation-status.json`: current stage, or a stage requiring attention.
- `RUN/launch-runtime-repair-01/owner.json`, `child.json`, `stdout.log`, `stderr.log`, and eventual `complete.json`: the active continuation supervisor. Its recorded supervisor PID at preparation was 52764; read the owner file and verify the command before relying on a PID.
- `RUN/control/dispatch.json`: the current batch. Each task directory has `owner.json`, logs, and an `exit.json` when its process finishes. A task in this batch can already be finished; inspect its files when counting remaining jobs.
- Each job's `attempts/NN/summary.rds` and `complete.rds`: actual fit/diagnostic status and hashes. Process exit zero alone does not mean a statistically acceptable fit.
- Each job's `selection.rds` and `prediction/summary.rds`: immutable diagnostic selection and prediction status.

Read existing records for status. Calling the CLI's `queue` or `select` modes can create requests and selections and is not a read-only status check. Do not launch a second dispatcher while the existing supervisor/controller is alive. Verify ownership through the recorded command before signalling any process. The original pilot dispatcher was deliberately held and drained for the E12 repair; older continuation logs can contain intentional exit 143 records. Those are not current study failures.

The active driver is `RUN/runtime-repair-v1/resume-runtime-repair.py`. It finishes the pilot with four workers, creates the runtime/storage checkpoint, checks that it passes, dispatches the remaining study with six workers, independently verifies full results, and scores them. Its terminal success stage is `scientific_results_ready`, with `lesson_export_complete=false`. It does not write the final teaching bundle, figures, report or lessons. A `needs_attention` stage needs diagnosis, not an automatic restart or weaker criteria.

## Scientific and resource commitments

The study has 192 comparison cells: 16 arms across three communities and four arrangements. There are also 96 species-specific oracle references. The pilot is community 1's spread and clustered arrangements, covering all arms and 16 oracle references. Pilot fits and predictions belong to the full study and are reused.

Fit occJSDM, sjSDM, Hmsc and gllvm only. Doug explicitly deferred spOccupancy's `sfJSDM()`. The Lesson 7 Backblaze archive was unavailable; audited reconstruction of the original seeded communities was accepted. Preserve the `reconstructed` provenance, the 61 audit checks and the uniform replay evidence. Later restoration is an independent comparison, not permission to replace inputs or previous fits.

Initial Bayesian budgets are four chains, 3,000 burn-in and 4,000 retained draws. The single permitted extension uses four chains, 8,000 burn-in and 16,000 retained draws. Point models have three starts, with three more only under the declared agreement rule. Diagnostics, prediction targets, links, uncertainty limitations and range comparability are frozen in the specification and adapters. Keep failed, unstable and unresolved cells in the accounting. Keep truth out of fitting and selection, and do not open recovery scores before the timing checkpoint.

Doug authorised six concurrent single-threaded R jobs for the main run, while noting that workers five and six would be slower. The pilot remains at four. The conservative remaining-time forecast still assumes four-worker capacity and must be at most 48 hours. The separate per-attempt limit is 24 hours. Attempt ceilings remain 768 comparison plus 192 oracle attempts. A failed time/storage gate requires a documented amendment or another host; it does not justify dropping difficult cells or weakening diagnostics.

## Runtime and source identity

Use the pinned `run-r` launcher, not the system `Rscript`. System R is 4.6.1; these compiled dependencies use R 4.5.0 and were observed to segfault under the newer runtime. The launcher pins the private R home, library, Python environment and one-thread settings. It does not change global installations.

The original scientific snapshot is `RUN/snapshot`, frozen at commit `c448cb2bdc7808b8e6443661999aa0182876725b`. Execution correction E12 is committed at `0601865243f552fa4d28f651b18a36c21b1553d4`. Amendments E11 and E12 are in the current specification; the original frozen protocol predates them. Read `RUN/execution-policy.json`, `execution-policy-runtime-repair.json`, and `runtime-repair-v1/manifest.json` together.

Do not edit the executing snapshot, control overlay, config, selections or checksum-bound artifacts. New fixes need a separately recorded correction with tests and preserved provenance. Development and reporting changes belong in the working checkout.

E12 fixes an administrative comparison between Hmsc namespace version `3.3-7` and the identical numeric package version `3.3.7`. The overlay checks that its runner differs from the original by exactly one normalization line, and retains installed/loaded path and installation-hash checks. It preserves original administrative summary/completion bytes under each recovered attempt's `version-recovery/` directory. `RUN/runtime-recovery.rds` records the completed recovery. No Hmsc selection was revised. The eight Hmsc pilot extensions were triggered by this administrative failure; do not describe them as an observed rate of poor mixing. The four held workers' original coordinator exit records and independently observed finish-time bounds are also preserved.

The E12 repair passed two loaded-version assertions, 17 recovery assertions, 44 runner assertions, six Python tests, a recovery test on a copy of a real Hmsc attempt, and independent review with no remaining Critical or Important findings. The earlier full infrastructure suite passed 215 top-level assertions plus native subprocess checks. These are infrastructure checks, not evidence of good scientific recovery or convergence.

## Remaining deliverables

1. Let the authorised execution finish, or investigate its recorded failure. Inspect `runtime-checkpoint.rds`, `verification-full.rds` and `scores-full.rds` when available. Preserve every failure and the distinction between diagnostics, numerical verification and recovery of truth.
2. Finish Task 6's independent checks. The current verifier reconstructs saved training predictions and Bayesian audit draws/means; supported interval summaries and Monte Carlo SE checks still need the broader audit described in the plan. Complete score, attribution, group-weight and denominator verification without changing frozen estimands after inspecting results.
3. Finish Task 7's compact teaching export and figure generation. Five untracked drafts exist in the working checkout: `export-teaching.R`, `verify-teaching.R`, `plot.R`, `tests/export-fixture.R`, and `tests/test-export.R`. They are unfinished and are not used by executing workers. The pure export tests previously passed 18 assertions; plots have not had visual QA, and no scientific bundle or plots have been generated. Add E11/E12 execution provenance and administrative-extension accounting to the export. Persist explicit Moran-basis distance and basis-amplitude summaries, the deferred minor review item. Finish CLI integration and the archive-free arithmetic audit.
4. Write the comparison report and extend Lesson 5, with a link from Lesson 7's spatial-effects explanation. Teach pattern versus amplitude recovery, held-out prediction, uncertainty limitations, and rare-species per-site attribution against the oracle. State unsupported features and compare ranges only when definitions and scales match. Use the existing specification for exact teaching outputs.
5. Render and inspect figures and lessons, complete the final whole-branch review, and report the evidence and remaining limitations. Doug reviews lesson prose before publication. Preserve publication flags and the quickstart's spatial advice; this task does not authorise publishing or a production sampler/API change.

The five reporting drafts have a verified persistent copy at `RUN/handoff/reporting-drafts-20261008`, with a SHA256 manifest. Preserve newer edits if those files have changed by handover time. The current handover itself also has a persistent copy under `RUN/handoff/`. No pull request has been created or branch pushed for this work.

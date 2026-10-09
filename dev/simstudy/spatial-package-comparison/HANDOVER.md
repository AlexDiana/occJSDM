# Spatial comparison handover to Claude

The spatial study is paused after E14's authorised window. The controller and supervisor exited normally at 03:04:37 BST on Friday 9 October 2026, after stopping new starts at 03:00. No unfinished attempt was interrupted. Leave the study paused until Doug explicitly authorises resumption. The pilot remains complete, and the full comparison remains unfinished. This handover supports a possible transfer to Claude on Saturday 10 October; read the current pause receipt and records below when taking over.

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

RUN is outside the disposable worktree. Its `jobs` symlink points to the external research volume. While paused, Doug can eject the physical drive normally before disconnecting it. During future execution, keep that volume mounted and the Mac on mains power with the lid open. The supervisor uses `caffeinate`; it does not make shutdown or unplugging an active archive safe. The separate Time Machine backup volume is not part of this study.

## Live execution

**E14 window finished at 03:04:37 BST on 9 October:** `RUN/window-20261008-2224/coordinator/complete.json` records `paused_window`, 252 completed tasks, zero interrupted tasks, 142 deferred tasks from the latest staged queue, no queue interruption and zero owned processes remaining. The supervisor completion record reports exit zero. All 252 worker exit records also report zero, and every recorded start preceded 03:00. The 240 fit attempts comprise 40 occJSDM initial fits, 120 sjSDM starts and 80 oracle initial fits. Twelve additional occJSDM prediction jobs completed. These process counts do not establish diagnostic acceptance or scientific recovery. The 142 deferred tasks are only the latest staged queue, not the complete remaining study; Hmsc/gllvm work, subsequent diagnostic attempts and predictions still remain.

`RUN/window-20261008-2224/pause-complete.json` is the current final pause receipt. Its supporting integrity audit passed for all 252 completed tasks and 2,072 artifact files (4.14 GB), checking completion identities and saved hashes without fitting, selecting or scoring. A fresh process check found no owned controller, supervisor, R worker, queue or descendant process, and `lsof` found no open files anywhere under the research archive at 03:10:59 BST. Numerical prediction/interval verification and all recovery scoring remain pending. The historical E13 pause receipt and earlier E14 launch handover are preserved separately. Recorded supervisor PID 10898, coordinator PID 10899 and all worker PIDs are historical; never signal or resume them without fresh identity checks.

The deadline heartbeat `pause-spatial-study-before-laptop-use` must remain paused after this confirmed shutdown. It must never restart the study. Read `RUN/window-20261008-2224/heartbeat-pause-confirmation.json` for the final automation confirmation. The persistent copy of this updated handover is `RUN/handoff/HANDOVER-20261009-reporting-drafts.md`. The E14 pause receipt continues to bind the preserved historical copy `RUN/handoff/HANDOVER-20261009-paused-E14.md`.

**Historical E13 pause completed on 8 October at 20:28 BST:** Doug needs the laptop at 09:00 BST on Friday 9 October. Read `RUN/laptop-pause-complete.json` and `RUN/laptop-deadline-hold.json` before interpreting older records. The pilot coordinator completed with an empty queue after verifying completion identities and artifact MD5s. Only then was supervisor PID 52764 terminated; its cleanup stopped held driver PID 52765 without sending SIGCONT. All owned study processes exited, and `lsof` found no open files under the study job archive. No unfinished fit or prediction was interrupted. The supervisor's exit 143 records this intentional pause. The thread heartbeat was PAUSED at that point; E14 subsequently reactivated it with a new monitoring prompt. Do not revive any recorded PID or run the old continuation without new authorisation.

All 48 pilot cells have immutable selections, and all 47 expected prediction jobs completed. One sjSDM cell has no usable fit and therefore no prediction. Diagnostic selections are: occJSDM six selected and two unresolved; Hmsc five selected and three unresolved; gllvm six selected and two unstable; sjSDM six selected, one unstable and one failed; all 16 oracle references selected. Completion does not mean every fit passed diagnostics. All 16 Hmsc attempts were recovered without refitting, and the eight administrative extensions remain labelled as such. No scientific recovery scores have been opened.

The immutable `RUN/runtime-checkpoint.rds` passed both gates: 31.23 conservative remaining wall hours at four-worker capacity, against the 48-hour limit; 90.88 GB estimated additional archive storage plus the 10 GiB reserve, with approximately 1.147 TB available when checked. It contains 109 completed-fit timing records and 47 prediction timings. This permits later execution under the frozen budget but does not override the user pause. The checkpoint is a resource projection, not a recovery result or guaranteed finish time.

Use these files as the current source of execution status:

- `RUN/continuation-status.json`: current stage, or a stage requiring attention.
- `RUN/window-20261008-2224/pause-complete.json`, `pause-process-check.json`, `pause-drive-check.json` and `pause-integrity-audit.json`: current E14 pause receipt, verified process exit evidence, archive-open-file check and completion-integrity audit. The receipt binds supporting records by SHA256. Root `laptop-pause-*` files refer to the historical E13 pilot pause.
- `RUN/launch-runtime-repair-01/owner.json`, `child.json`, `stdout.log`, `stderr.log`, and `complete.json`: the finished continuation supervisor. Its intentional exit is 143, its final pilot queue completed, and stderr is empty. Recorded PIDs are historical and must never be signalled without fresh identity checks.
- `RUN/control/dispatch.json`: the current batch. Each task directory has `owner.json`, logs, and an `exit.json` when its process finishes. A task in this batch can already be finished; inspect its files when counting remaining jobs.
- Each job's `attempts/NN/summary.rds` and `complete.rds`: actual fit/diagnostic status and hashes. Process exit zero alone does not mean a statistically acceptable fit.
- Each job's `selection.rds` and `prediction/summary.rds`: immutable diagnostic selection and prediction status.

Read existing records for status. Calling the CLI's `queue` or `select` modes can create requests and selections and is not a read-only status check. Do not launch a second dispatcher while the existing supervisor/controller is alive. Verify ownership through the recorded command before signalling any process. The original pilot dispatcher was deliberately held and drained for the E12 repair; older continuation logs can contain intentional exit 143 records. Those are not current study failures.

The archived driver is `RUN/runtime-repair-v1/resume-runtime-repair.py`. When explicitly restarted, it checks the recovered pilot with four workers, reuses the existing snapshot-bound runtime/storage checkpoint, checks that it passes, dispatches the remaining study with six workers, independently verifies full results, and scores them. Its terminal success stage is `scientific_results_ready`, with `lesson_export_complete=false`. It does not write the final teaching bundle, figures, report or lessons. A `needs_attention` stage needs diagnosis, not an automatic restart or weaker criteria.

## Restart only after Doug authorises

1. Confirm Doug authorises resumption and that his next laptop/drive availability window is adequate. Record that authorisation beside the pause receipt; preserve the historical hold and completed records.
2. Remount the research volume. Verify `RUN/jobs`, the storage sentinel and volume UUID `8667ADC0-D4C2-4B79-85DF-C44B181E8C24`; check current free space against the checkpoint's remaining bytes plus reserve. Recheck the frozen snapshot, input/runtime identities and every file in the E12 control manifest. A previously passing capacity check does not establish current free space.
3. Confirm there is no live owned supervisor, coordinator or R worker. Do not rerun the pilot fits or add attempts to cells with exhausted diagnostics. The driver's pilot queue rechecks completion checksums and reuses selections and predictions.
4. Prefer a new bounded window using the reviewed `coordinate-window.py` and `launch.py`, with a fresh policy and record directory reflecting Doug's newly authorised availability. Preserve the completed E14 window and all existing requests, selections and attempts. Recheck deferred requests as never started before dispatch. The older `resume-runtime-repair.py` proceeds into unrestricted fitting, verification and scoring; do not launch it unless Doug explicitly authorises that scope and duration. Never reuse a supervisor log directory or send SIGCONT to historical PIDs. Preserve the existing immutable checkpoint; its 31.23-hour figure is the pre-E14 projection, not a current countdown.
5. Confirm the new owner/child records and live commands, then update the handover with the new supervisor directory and status. Any failure needs investigation under the frozen scientific criteria. The paused deadline heartbeat must not restart the study.

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

1. Resume only after fresh authorisation and the restart checks above. E14 finished 240 additional fit attempts and 12 additional prediction jobs; all started work is saved, and the pilot is retained as part of the full study. Complete remaining Hmsc/gllvm work, diagnostic extensions or extra starts where prescribed, and remaining predictions. Do not mistake the 142 currently deferred tasks for the entire remainder. Inspect `verification-full.rds` and `scores-full.rds` when available; neither has been produced yet. Preserve every failure and the distinction between diagnostics, numerical verification and recovery of truth.
2. Finish Task 6's independent checks. The current verifier reconstructs saved training predictions and Bayesian audit draws/means; supported interval summaries and Monte Carlo SE checks still need the broader audit described in the plan. Complete score, attribution, group-weight and denominator verification without changing frozen estimands after inspecting results.
3. Finish Task 7's compact teaching export and figure generation. Five unfinished drafts are now tracked in a draft checkpoint on this branch: `export-teaching.R`, `verify-teaching.R`, `plot.R`, `tests/export-fixture.R`, and `tests/test-export.R`. They are unfinished and are not used by executing workers. The pure export tests previously passed 18 assertions; plots have not had visual QA, and no scientific bundle or plots have been generated. Add E11/E12 execution provenance and administrative-extension accounting to the export. Persist explicit Moran-basis distance and basis-amplitude summaries, the deferred minor review item. Finish CLI integration and the archive-free arithmetic audit.
4. Write the comparison report and extend Lesson 5, with a link from Lesson 7's spatial-effects explanation. Teach pattern versus amplitude recovery, held-out prediction, uncertainty limitations, and rare-species per-site attribution against the oracle. State unsupported features and compare ranges only when definitions and scales match. Use the existing specification for exact teaching outputs.
5. Render and inspect figures and lessons, complete the final whole-branch review, and report the evidence and remaining limitations. Doug reviews lesson prose before publication. Preserve publication flags and the quickstart's spatial advice; this task does not authorise publishing or a production sampler/API change.

The five reporting drafts have a verified persistent copy at `RUN/handoff/reporting-drafts-20261008`, with a SHA256 manifest. Preserve newer edits if those files have changed by handover time. The current handover itself also has a persistent copy under `RUN/handoff/`. No pull request has been created or branch pushed for this work.

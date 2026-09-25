# Lesson 4 replicated extension

Status: all 640 planned combinations were attempted, exported and verified on 24 September 2026. There are 626 successful scored fits: 505 passed the declared diagnostics and 121 remain flagged. The 14 failed gllvm fits are preserved. Completion means the fixed experiment has finished, not that every fit passed its numerical checks. See [PLAN.md](PLAN.md) for the approved scenarios, statistical targets and fitting budgets.

## Where everything lives

- Working tree: `/Users/douglasyu/src/occJSDM/.worktrees/lesson-4-extension/`, branch `codex/lesson-4-extension`.
- Run archive: `/Users/douglasyu/src/occJSDM/dev/simstudy/results/lesson-4-extension-20260923/`.
- Source lesson: `vignettes/occJSDM-lesson-4.Rmd` in the worktree. Its extension reads `vignettes/teaching-data/lesson-4-extension.rds`.
- The existing frozen R/Python installation is reused read-only through `/Users/douglasyu/Documents/Codex/2026-09-10/fam/work/lesson-n-environment-20260922/run-r`. No packages or Mojo versions were upgraded. All newly created source, output and logs are within `~/src/occJSDM/`.

The run manifest declares 640 package/dataset/model combinations before fitting, excluding restarts. `jobs/*/status.json` distinguishes a completed attempt from a successful fit and from passed diagnostics. A failed completed fit is preserved, not silently rerun. `score.rds` exists only after fitted parameters and selection were saved. `completion.json` is written by the supervisor when the complete pipeline finishes or needs attention. It is absent while the supervisor is running. Read `logs/full-study.log` for current activity.

## Resume and reproduce

The completed run used the immutable `worker-normal-integral-v1/` source snapshot. The previous `worker-final/` snapshot is retained unchanged. Do not edit either snapshot. Source corrections belong in the worktree, followed by a new snapshot and an explicitly documented retry; preserve affected old records. The command below records the resume procedure; the finished archive does not need further fitting.

```sh
study_root="$HOME/src/occJSDM/dev/simstudy/results/lesson-4-extension-20260923"
study_repo="$HOME/src/occJSDM/.worktrees/lesson-4-extension"
study_launcher="$HOME/Documents/Codex/2026-09-10/fam/work/lesson-n-environment-20260922/run-r"

python3 "$study_root/worker-normal-integral-v1/extension/finish-study.py" \
  "$study_root" "$study_root/worker-normal-integral-v1" "$study_repo" \
  --launcher "$study_launcher" --workers 4
```

A live `runner.lock` prevents two supervisors from fitting the same rows. Do not remove it just to force a resume. After an interrupted process, verify that the recorded PID is no longer running before removing a stale lock. The runner checks hashes before accepting completed results and stops scheduling when disk space drops below 8 GiB. It leaves recorded numerical failures visible. A scoring or export error stops the supervisor and records `needs_attention`, rather than publishing a misleading completed report.

The supervisor fits and scores the fixed manifest, exports the compact bundle, independently verifies identifiers/truth/errors, then renders both HTML and GitHub Markdown. It does not commit, push or merge. Source edits and generated figures remain reviewable in the worktree.

The 23 September scoring recovery starts from 79 verified completed fits and 78 scores. Job 39 (`curved-r01-n300-s10-quadratic-traits0-sjSDM`) failed only during scoring. The corrected adaptive helper uses the equivalent normal-density integral for ordinary residual SDs and retains the probability-domain identity for SD >= 10, where it avoids a nearly discontinuous normal-domain integrand. Six earlier point-model scores also exceeded the integration tolerance because the old formula silently missed rare-event probability mass. All 35 existing point-model scores are preserved under `implementation-checks/scoring-integral-recovery/scores-before/` and regenerated from their saved fits, along with job 39's missing score. Bayesian scores and all fit records are retained. `source-manifest-normal-integral-v1.json` records the replacement snapshot; the recovery directory contains original failure logs, the pre-resume file inventory and numerical checks. No completed fit is repeated.

## Checks

```sh
Rscript dev/simstudy/jsdm-package-comparison/test-pilot-math.R dev/simstudy/jsdm-package-comparison
Rscript dev/simstudy/jsdm-package-comparison/extension/test-study.R dev/simstudy/jsdm-package-comparison/extension
Rscript dev/simstudy/jsdm-package-comparison/extension/test-math.R dev/simstudy/jsdm-package-comparison/extension
PYTHONDONTWRITEBYTECODE=1 python3 dev/simstudy/jsdm-package-comparison/extension/test-runner.py
```

`verify.R RUN_ROOT REPOSITORY_ROOT` checks available scored results, including an independent normal integral for true probability and recomputation of signed and absolute errors from prediction cells. Rendering never fits models. The initial source-interface problems and their original records are retained in `implementation-checks/` and described in PLAN.md. In particular, the initial sjSDM driver incorrectly assigned through `reticulate::py`; the corrected binding was tested in an initialised Python session and the complete native sjSDM job was rerun. Three independent starts then passed the declared stability checks on that first dataset. This is an interface check, not evidence that every later dataset will fit stably.

Pre-commit verification on 25 September 2026: all four checks above passed, and the archive verifier independently checked all 626 scored results and the 640-job manifest. The compact results and six teaching figures are included with Lesson 4; full fits, failed attempts and recovery records remain in the ignored run archive. No completed fit was repeated during scoring recovery or these checks.

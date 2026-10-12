# First spatial package comparison: record

Written 12 October 2026. On that day Doug decided that this study is not worth keeping: it was underpowered, so its fits, inputs and scores are being discarded and only this record is kept. The spatial model work begins with the recoverability pilot (`dev/simstudy/spatial-recoverability-pilot/REPORT.md`), which is the positive control this study lacked.

## What it was

A comparison of whether a different way of representing space would let a model recover spatial structure better than occJSDM's own Gaussian-process field, which at 100 sites in the Lesson 7 site-arrangement sweep learned almost nothing. It was built and run by Codex on 8 to 10 October 2026, on the branch `codex/spatial-comparison-plan`.

- **Data.** The sweep's three seeded binary communities, eight species each, with regional occupancy targets of 5, 25 and 75 percent. Each was surveyed under four site arrangements (spread, pairs, clustered, grid) of 100 training sites, and scored on a held-out lattice of 1,600 sites. The generating spatial field had a squared-exponential covariance with range 0.03 of the landscape width and a logit-scale standard deviation of 1. The original input archive had been lost, so the inputs were regenerated from the original generator and seeds and checked against surviving evidence.
- **Arms.** Four packages (occJSDM, sjSDM, Hmsc, gllvm), each fitted four ways: environment only, a quadratic trend surface in the coordinates, Moran eigenvectors, and the package's own spatial term. That is 16 arms and 192 comparison cells, plus 96 per-species oracle references.
- **Budget.** Four chains of 3,000 warm-up and 4,000 retained draws, with one prescribed extension when the convergence screen failed (Rhat at most 1.01, effective sample size at least 400).

## What happened

All 288 jobs ran and the numerical reconstruction of every fit passed its checks, in 462 comparison attempts and 165 oracle attempts. Of the 192 comparison cells, 146 passed selection, 29 were unresolved (they failed the convergence screen after the extension), 6 were unstable point fits and 11 had no acceptable fit. All twelve occJSDM GP fits and all twelve Hmsc two-factor fits were unresolved, and sjSDM's coordinate network produced no accepted fit.

## What it found

Nothing that could decide anything. Even the oracle, which knew every generating parameter except the hidden field, barely recovered held-out spatial structure: the correlation between its predicted spatial contrast and the truth was 0.074, 0.200 and 0.207 for the 5, 25 and 75 percent occupancy groups, and its field attenuation slopes were 0.011 to 0.043. The occJSDM GP improved held-out probability error over an environment-only fit by only 0.04 to 0.12 percentage points, and every GP fit was unresolved. The Hmsc spatial factors improved it by even less. The gllvm spatial factors slightly worsened it. The fixed trend surfaces and Moran bases worsened held-out prediction in all four packages.

## Why it was underpowered

The design had no positive control: no setting in which the oracle was shown to recover the structure at held-out sites. It used the hardest case, 100 sites with range 0.03. The recoverability pilot later showed that range 0.03 could not be certified at any site count it tried (100, 300 or 900), and that the short-range class of the per-species fields fails at 100 sites. So the study could not separate package differences from a lack of recoverable information, and its own report said so and recommended establishing a positive control first. That recommendation became the recoverability pilot.

## What was kept and what was discarded

- **Kept.** This record. The adapters, basis, diagnostics and selection code that the pilot and the shelved stage 2 reused, which are copied into `dev/simstudy/spatial-stage2-comparison/upstream/` with their hashes in `sources.csv`. The pinned R 4.5.0 home and the occJSDM library (the `runtime/` and `library/` folders of the laptop's git-ignored folder `dev/simstudy/results/spatial-package-comparison-20261008/`), because the pilot's and stage 2's launchers (`run-r`) point at them. The full report, design, code and a Lesson 5 draft section with 19 figures are reachable through the local tag `archive/first-spatial-comparison` (commit `f7ec484`). That tag has never been pushed, and the branch was never merged or pushed.
- **Discarded on 12 October 2026.** The laptop's run index apart from the two folders above (reconstructed inputs, truth, bases, scores and verification files, about 2.1 GB). The branch `codex/spatial-comparison-plan` and its worktree. The Lesson 5 draft section and its figures were only ever on that branch, so no lesson on `main` depended on this study.
- **To be discarded when the research drive is reconnected** (it was unavailable on 12 October 2026): the 288 fit directories and the 342 MB teaching bundle `spatial-package-comparison.rds` (358,611,608 bytes) in `/Volumes/occJSDM Research/occJSDM/spatial-package-comparison-20261008/`. This record will be updated when that is done.

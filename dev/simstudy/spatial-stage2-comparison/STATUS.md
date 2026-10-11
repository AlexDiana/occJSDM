# Stage 2 status

This file is the live state of the stage 2 spatial representation comparison, for Doug to read on return. It is updated and committed at every phase change. The newest entry is at the top.

Branch `claude/spatial-stage2-comparison`, worktree `~/src/occJSDM/.worktrees/spatial-stage2-comparison`. Plan `dev/superpowers/plans/2026-10-11-spatial-stage2-comparison.md`, spec `dev/superpowers/specs/2026-10-11-spatial-stage2-comparison-design.md`. Study directory (ignored) `dev/simstudy/results/spatial-stage2-20261011/`. Job files on the research volume `/Volumes/occJSDM Research/occJSDM/spatial-stage2-20261011/jobs/`.

## Current phase

- 12 October 2026, 00:05: micro-timing exposed a feasibility blocker for gllvm. The first micro-timing run (5 workers, 10 micro fits at 900 sites) was stopped after 2.5 hours because the two `gllvm/spatial2` fits at 900 sites each grew to about 38 GB (top's MEM column, of which 34 to 37 GB compressed) on a 32 GB machine, filled 28 of 28.6 GB of swap, and ran at 12 to 22 percent CPU. gllvm's spatial variational fit stores dense 900 by 900 by 2 structures and an automatic differentiation tape. The gllvm settings are not part of the timing ladder, so no ladder step can change this. The micro fits are discarded by design and no study attempt exists. The micro-timing is being redone so that each micro fit saves its own result, with the gllvm 900-site fit watched under a memory cap. Details follow in the next entry.
- 11 October 2026, 21:30: built and prepared. Every test file passes under `run-r`, the code is committed, and `prepare` has written the frozen config and the 40-dataset table (hashes checked against the pilot manifest). Nothing has been fitted. Next: the micro-timing, then the timing ladder from step 1.

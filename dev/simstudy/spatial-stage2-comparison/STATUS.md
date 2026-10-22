# Stage 2 status

**Shelved on 12 October 2026, before any comparison fit, until server time is available.** No study attempt exists. `SHELVED.md` in this folder is the full account: what was built, what was measured, why the study cannot run on the laptop it was built on, and how to resume.

- **Why.** The gllvm arm needs more than 33 GB of memory at 300 sites on a 32 GB machine, and the Hmsc and gp arms at 900 sites project to thousands of wall-hours against the 48-hour limit. Hmsc with the GPP method also fails at 200 or more knots.
- **Re-check, 14 October 2026.** With the laptop's R on Apple's Accelerate BLAS, the gp and Hmsc arms run 5 to 60 times faster per iteration, and the study projects to about 26 wall-hours at eight workers if no job needs the extension and 116 if every job does, at the cheapest ladder step (`SHELVED.md` section 9). Still outside the 48-hour rule; gllvm's memory problem is unchanged, and the pinned R 4.5.0 runtime still loads the reference BLAS, so the runtime must be switched or replaced before the study sees the gain.
- **State of the work.** The code, tests, frozen configuration and dataset table are complete and committed on branch `claude/spatial-stage2-comparison`, so the study can resume without rebuilding anything. The compact evidence is in `results/`.
- **To resume.** Decide the gllvm arm and the Hmsc and gp settings (each is a spec amendment), then follow `README.md` and section 6 of `SHELVED.md`.
- **Files.** The laptop's study directory is `dev/simstudy/results/spatial-stage2-20261011/` (git-ignored). The research drive's folder for this study is an empty skeleton, and `/Volumes/occJSDM Research/occJSDM/README.md` explains every folder on the drive.

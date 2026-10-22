# v0.1.1-beta release checks, 8 October 2026

Release branch `release/v0.1.1-beta`. The checked source is commit `195bb7d` (version bump to 0.1.1); the only later commit on the branch, `ea12937`, changes `AGENTS.md`, which `.Rbuildignore` excludes, so the package contents are identical. Changes since `v0.1.0-beta` are listed in `NEWS.md` under 0.1.1. Model code, priors and the shipped example (`data/sampleresults.rda`, `data/sampledata.rda`) are unchanged, so the example was not refitted.

## Environment

All checks ran under R 4.6.1 on macOS 27.0.1, Apple silicon, with Apple clang 21.0.0; the last release was checked under R 4.5.0. R 4.6 became the system default on 7 October 2026, and its packages were installed fresh from CRAN into the user library on 8 October, including ggtern 4.0.0, so the ternary-plot compatibility fix is exercised. BLAS and OpenMP thread variables were set to 1. Full session details are in [session-info.txt](session-info.txt).

One build difference: under this machine's R 4.6 configuration the compiled code links Homebrew's `libomp` (`-lomp`). No OpenMP compile flag is used and the package's threading comes from RcppParallel, so this does not change the sampler. It is a property of the local R setup, not of the package.

## Results

-   **Source tests** ([source-tests.log](source-tests.log)): 1,145 expectations, 0 failures, 0 errors, 0 warnings, 1 skip (the opt-in coverage study). The last release had 1,050; the new tests cover single-species fits, `returnPosteriorDraws()` and the ggtern compatibility fix.
-   **Package build** ([package-build.log](package-build.log)): both vignettes (quickstart and simulator guide) built. Archive `occJSDM_0.1.1.tar.gz`, 52,585,460 bytes, SHA-256 `d7fa080490d6a934b6c9414d858e2ebf3cdcffd8fa0dc9ddd54c9548a23e1c41`. It contains `NEWS.md` and excludes Lessons 0-7, teaching bundles and development records.
-   **`R CMD check --no-manual`** ([package-check.log](package-check.log)): 0 errors, 2 warnings, 2 notes. Installed-package tests pass 1,141 expectations with 0 failures and 3 standard skips ([installed-tests.log](installed-tests.log)); examples run and vignettes rebuild.

The warnings and notes are the same categories recorded for v0.1.0-beta: the undocumented `predictNewSites(verbose)` argument, GNU extensions in `src/Makevars`, `LICENSE` not mentioned in `DESCRIPTION`, and existing R-code findings (two partial argument matches and undefined globals in internal helpers). The previous release's third warning, about the R 4.5 header's unsupported clang option, does not occur under R 4.6. No new diagnostic category appears.

## Not checked

Lesson rendering, the lesson-site tests and the MCP pilot were not rerun: none of them changed since `v0.1.0-beta` and none ships in the package. Hmsc, gllvm and sjSDM remain installed only under R 4.5. These are macOS checks only, not CRAN acceptance or cross-platform equivalence.

## Release contents

The proposed GitHub prerelease description is [release-notes.md](release-notes.md). Proposed assets: the checked `occJSDM_0.1.1.tar.gz` with a `SHA256SUMS` file, and the unchanged `occJSDM-cloud-pilot.zip` carried over from `v0.1.0-beta`. **Published 8 October 2026** as the [GitHub prerelease](https://github.com/AlexDiana/occJSDM/releases/tag/v0.1.1-beta), tag `v0.1.1-beta` on merge commit `cdb29da`, whose package contents are identical to the checked commit. All three assets were downloaded back from the release and match `SHA256SUMS`; the cloud-pilot zip is byte-identical to the v0.1.0-beta asset.

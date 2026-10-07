# Beta release refresh implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task by task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Refresh the shipped example and beta vignettes on merged main, verify the installed release package, and prepare reviewable beta release materials.

**Architecture:** Fit the existing quickstart model from an isolated installation of main `a56f5548e68e09d244513d09fcc04ad52f0caa19`, retaining its defaults and two-chain schedule. Save explicit provenance beside the development report, rebuild only the two beta vignettes, and check a clean source archive with vignette building enabled. Keep the historical scientific and lesson evidence intact.

**Tech stack:** R 4.5.0, Rcpp, testthat, rmarkdown/knitr, Git and GitHub CLI.

**Spec:** User's 7 October 2026 instruction to commence the beta release preparation described in this chat; `TODO.md`, release preparation items 4 and 5.

## Global constraints

- Preserve the primary checkout's unrelated MCP planning edits and all historical simulation archives.
- Retain the quickstart's spatial specification, two factors, threshold 1, default priors, two chains, 5,000 burn-in iterations and 5,000 retained iterations per chain, with no thinning.
- Set the fitting seed to 20261007 and record one requested RcppParallel thread and one BLAS/OpenMP thread.
- Keep Lessons 0-7 withheld. Their scientific rebuilding and publication remain deferred.
- Write paragraphs in planning documents on one line and use no em-dashes.
- Changes to model behavior require diagnosis and a regression test; this refresh initially changes data and documentation only.

## Review focus

1. The shipped fit must be produced by the current installed package and match the quickstart's input and settings, with convergence concerns disclosed.
2. Rebuilding must use the refreshed fit, retain threshold and prior metadata, and include both beta vignettes in the installed archive.
3. Remaining package-check warnings must be classified honestly; no claim of CRAN readiness is made.
4. The announcement must match the README's limitations and say that teaching lessons follow later.
5. Historical teaching bundles, simulation studies and the primary checkout's unrelated changes must remain untouched.

## Task 1: Refit and verify the shipped model

**Files:** Create `dev/release/2026-10-07-beta/refit-sampleresults.R`, `provenance.rds`, `diagnostics.csv` and `session-info.txt`; modify `data/sampleresults.rda`.

**Interfaces:** The refit script consumes the repository root, a private installed-package library and an evidence output directory. It produces `sampleresults` in the documented six-element fit format and records source/data hashes, arguments, seeds, threads, warnings and elapsed time.

- [ ] Install main into a private library and run the existing source tests as the baseline. Expected: zero test failures and no new test warnings.
- [ ] Run the exact quickstart model with the fixed seed and thread settings; validate model identity, dimensions, threshold, default prior metadata, and finite summaries. Expected: a complete two-chain fit and a readable provenance record.
- [ ] Export convergence diagnostics and inspect flagged parameters before accepting the saved example. Expected: findings disclosed, with a longer run considered only if the example cannot responsibly support the quickstart.

## Task 2: Refresh beta documentation and release materials

**Files:** Modify `vignettes/occJSDM.Rmd`, its generated Markdown/figures, `vignettes/simulateOccJSDMData.md` and any regenerated simulator figures, `TODO.md`; create `dev/release/2026-10-07-beta/REPORT.md`, `release-notes.md` and `announcement.md`.

**Interfaces:** Rendering consumes the refreshed installed data. Release notes identify version 0.1.0 and the GitHub prerelease tag `v0.1.0-beta`; the announcement describes the available quickstart and simulator guide and links current limitations.

- [ ] Show the refit seed in the quickstart and document the example's provenance and any meaningful convergence flags without changing the taught model specification. Expected: teaching remains complete and accurate.
- [ ] Render both beta vignettes to HTML and GitHub Markdown in fresh R sessions, with unwrapped Markdown; visually inspect regenerated plots. Expected: successful standalone rendering and current output.
- [ ] Mark PR #14 merged and record actual release-check progress in TODO; reconcile the announcement's lesson-publication wording. Expected: no premature release or announcement completion claim.

## Task 3: Check the distributable package and prepare integration

**Files:** Update the release report with final verification; create a clean source archive and check logs in a temporary execution directory.

**Interfaces:** The release source archive includes the refreshed data and the quickstart/simulator vignettes, excluding all withheld lessons and internal development material.

- [ ] Build and run `R CMD check --no-manual` with vignette construction and rebuilding enabled. Expected: zero errors; inspect and record every warning/note and installed test result.
- [ ] Verify both installed vignettes and shipped data; run source tests, lesson-link and publication-flag checks on the final tree when relevant changes justify them. Expected: current fit metadata and only the two beta vignettes shipped.
- [ ] Obtain one independent whole-change review, address material findings, and commit the refresh on a `codex/` branch. Expected: a reviewable release refresh with evidence attached and unrelated work preserved.
- [ ] Prepare the GitHub release description and announcement for the reviewed release revision. Any final shared-branch integration, release publication or outbound announcement must use established authorization and actual destinations; do not invent announcement recipients.

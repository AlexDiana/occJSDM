# Lesson verifiers after the plotting-colour change: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the lesson verifiers pass on main again without refitting, by giving them a code fingerprint that ignores comments and allows declared plot-only changes; re-export the figure bundles with the current package; drop Lesson 3's own colour scales now that the package sets them; and make the coloured figures in Lessons 2, 5, 6 and 7 colour-blind-safe.

**Architecture:** The verifiers compare the source hashes recorded when the fits were made with `lesson_source_hashes()` on the current tree. A new helper in `dev/simstudy/vignette-lesson/helpers.R` replaces those `identical()` checks: for each R file whose raw md5 differs, it finds the recorded version in git history, parses both versions, and compares them function by function with comments dropped; differences are allowed only in a declared list of plot-only functions. The archive's package library is reinstalled from current main (the old one kept), and every bundle that records the library's files is re-exported from the existing fits. Lessons 2, 5, 6 and 7 draw their figures live, so their colours change in the lesson code alone.

**Tech Stack:** R 4.5, rmarkdown, knitr, ggplot2, testthat, git; the full-fit archive at `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples`.

**Spec:** Doug's decisions of 4 October 2026: verifiers may ignore comment lines; for #23's plot-only code changes, use a function-level fingerprint with a named allowlist and re-export the figure bundles from the existing fits (no refitting); drop Lesson 3's manual colour scales; colour-blind fixes for Lessons 2, 5, 6 and 7. Writing rules for lesson prose: `dev/simstudy/vignette-lesson/LESSON-STYLE-SPEC.md`, "The reader" and "The criteria, applied".

## Facts checked while planning (4 October 2026, main at 471611b)

- Since the fits' source commit (fbe3ed6), `R/` changed in `R/simulateData.R` (roxygen only, #22), `R/output.R` (`plotFPTPStage2Rates()`, `plotDetectionRates()`, `plotStage2FPRates()`: colour scales), `R/diagnostics.R` (`plotTraceplot()`: colour scale) and a new `R/colours.R` (`okabe_ito()`), all #23. `src/`, `DESCRIPTION` and `NAMESPACE` are unchanged.
- On main, all 11 fit-reading lesson verifiers stop at the fingerprint (`identical(... source_hashes, lesson_source_hashes())`, or in `prediction-verify.R` `identical(b$source_hashes, tools::md5sum(names(b$source_hashes)))`). `verify-teaching.R` and the spatial sweep's `verify-lesson.R` pass.
- Files that check the installed library (`find.package` or `package_files_md5`): `export_native_plots.R`, `verify_native_plots.R`, `native-traits-export.R`, `native-traits-verify.R`, `ordination-export.R`, `ordination-verify.R`, `remaining-plots-export.R`, `remaining-plots-verify.R`, `covariate-effect-export.R`, `covariate-effect-verify.R`, `prediction-build.R`, `prediction-verify.R`, `remaining-plots-diagnose-covariate.R`.

## Global Constraints

- Worktree `/Users/douglasyu/src/occJSDM-worktrees/lesson-verify`, branch `codex/lesson-verify` from main at 471611b. Never work in `/Users/douglasyu/src/occJSDM`; never commit to main; never push (the controller pushes); never force anything.
- No change to `R/`, `src/`, `DESCRIPTION` or `NAMESPACE`. No refitting: every fit in the archive stays byte-identical (check `md5` of every `*-fit.rds` in the archive before and after).
- Commit messages end with a blank line and `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`, exactly, whatever model writes them.
- Markdown rules from `/Users/douglasyu/src/occJSDM/CLAUDE.md`: one line per paragraph; no em-dashes anywhere (code comments and commit messages included; use `--`); no pipe tables; escape a literal `>` `<` `~` `|` in prose, never inside code; re-read a file before and after editing it.
- Shell is zsh: `/usr/bin/grep`, `/bin/ls`.
- Step runner: `zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-VERIFY-PLAN/run-step.zsh STEP command args` (sets `R_LIBS` to the archive library, runs from the worktree, logs to the archive's `logs/` and `BUILD-LOG.md`). Step names start with `v` and are new.
- Render commands, from the worktree root:

```sh
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/<stem>.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/<stem>.Rmd", output_format="rmarkdown::html_vignette")'
```

- PNG policy: keep a changed PNG only when its chunk or its exported figure changed on purpose; restore render noise (`vignettes/occJSDM-lesson-0_files/figure-gfm/environmental-maps-1.png`, `vignettes/occJSDM-lesson-2_files/figure-gfm/probability-and-state-1.png`, `vignettes/occJSDM-lesson-3_files/figure-gfm/mirror-labelling-chains-1.png`, `vignettes/occJSDM-lesson-3_files/figure-gfm/trait-cancellation-1.png`) with `git checkout --`. Use `zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-VERIFY-PLAN/revert-unchanged-png.zsh PATHSPEC` to revert pixel-identical PNGs.
- Checks at the end of every task: `Rscript dev/simstudy/vignette-lesson/test_lesson_links.R`; `node --test dev/simstudy/lesson-site/test_lessons.js`; `Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/check-lesson-links.R vignettes/occJSDM-lesson-[0-7].Rmd vignettes/occJSDM-lesson-[0-7].md vignettes/occJSDM.Rmd vignettes/LESSON-PLAN.md TODO.md`; em-dash counts on `git diff 471611b` and on `git log 471611b..HEAD --format=%B` are 0.
- Stop and report if: a fit would need refitting; a change to `R/` would be needed; a verifier fails for a reason this plan does not expect.

---

### Task 1: A comment-blind, function-level source check

**Files:** Modify `dev/simstudy/vignette-lesson/helpers.R`; modify every lesson verifier and exporter that compares source hashes (find them with `/usr/bin/grep -rln "source_hashes" dev/simstudy/vignette-lesson/*.R`); create `dev/simstudy/vignette-lesson/test_source_check.R`.

**Produces:** `lesson_source_check(recorded, allow = lesson_plot_only_changes)` in `helpers.R`, returning `TRUE` invisibly or stopping with a message naming each file and function that differs; and `lesson_plot_only_changes`, a character vector of function names with a comment above it giving the date, the pull request and why each cannot affect a fit.

- [ ] **Step 1: Write the test first.** `test_source_check.R` (run with `Rscript`, `stopifnot()` style like `test_lesson.R`) builds a temporary git repository with an `R/` file, records `tools::md5sum()` of it, then checks that `lesson_source_check()`: passes when nothing changed; passes after a comment-only edit (`#` and `#'` lines, and a trailing comment); fails after a code change in a function not on the allow list, naming it; passes after a code change in a function on the allow list; fails when a recorded file is deleted; fails when a new file defines a function not on the allow list, and passes when the new file defines only allow-listed functions; and that non-R files (`DESCRIPTION`, `NAMESPACE`, `src/*.cpp`) are still compared by raw md5. Run it and see it fail (the function does not exist).

- [ ] **Step 2: Implement.** For each recorded path: equal raw md5 passes. Otherwise, for an `R/*.R` file, find the recorded version: search `git log --format=%H -- <path>` (from the repository root) for the first commit whose `git show <commit>:<path>` has the recorded md5; stop if none matches. Parse both versions with `parse(keep.source = FALSE)`, map each top-level assignment `name <- function(...)` (or `name = function`) to its deparsed definition, and compare: names added, removed or changed outside `allow` fail; other top-level expressions must be identical after deparsing. New `R/` files not in the record are parsed the same way: only allow-listed functions may appear. Non-R files compare by raw md5. Cache `git show` results within a run.

- [ ] **Step 3: The allow list.** `lesson_plot_only_changes <- c("plotFPTPStage2Rates", "plotDetectionRates", "plotStage2FPRates", "plotTraceplot", "okabe_ito")`, with the comment: added 4 October 2026 for PR #23 (Okabe-Ito colours), colour scales only, these functions draw plots from a finished fit and are never called while fitting. Before committing, check with the new function that these are exactly the functions whose code differs between fbe3ed6 and HEAD.

- [ ] **Step 4: Use it everywhere.** Replace every `identical(<recorded>$source_hashes, lesson_source_hashes())` and the `prediction-verify.R` form with `lesson_source_check(<recorded>$source_hashes)`, in verifiers and exporters alike. Leave the recording side (`source_hashes = lesson_source_hashes()`) unchanged, so new bundles record raw md5s of the current tree. Do not weaken any other check.

- [ ] **Step 5: Run.** The test passes. Run the verifiers through the step runner (`v1x-` names): they should now pass the source check and, for the bundles that record the library, stop later at the library or API checks; record exactly where each stops. Those are Task 2's job.

- [ ] **Step 6: Commit** `helpers.R`, the edited verifiers and exporters, and the test: "Make the lesson verifiers' source check ignore comments and allow plot-only changes".

### Task 2: Reinstall the archive library and re-export the figure bundles, dropping Lesson 3's manual scales

**Files:** the archive library (outside the repository); `dev/simstudy/vignette-lesson/remaining-plots-examples.Rmd` and `native-plot-examples.Rmd` (the four recoloured chunks); `vignettes/occJSDM-lesson-3.Rmd` and `.md`; every re-exported bundle and PNG under `vignettes/teaching-data/`.

- [ ] **Step 1: Record the fits.** md5 of every `*.rds` in the archive and its subfolders, to compare at the end.
- [ ] **Step 2: Reinstall.** Move `ARCHIVE/library` to `ARCHIVE/library-fbe3ed6` (keep it; do not delete), create a new `ARCHIVE/library`, copy into it every package from the old library except occJSDM (so dependency versions are unchanged), then `R CMD INSTALL --preclean --library=ARCHIVE/library .` from the worktree through the step runner (`v20-install`). Append a short paragraph to the archive's `BUILD-LOG.md` saying why.
- [ ] **Step 3: Drop Lesson 3's manual colour scales.** Remove the `scale_colour_manual()` lines Task 6 of the port added to `native-primer-1-example`, `native-primer-2-example`, `remaining-plots-stage2-fp` and `remaining-plots-detection`, identically in the snippets and in Lesson 3, so the package's own Okabe-Ito colours draw them. Keep the manual scale in `combine-rate-plots` only if the package no longer covers it (check); keep the invented figure's scale only if it is still needed (`plotDetectionRates()` now colours primers A, B, C in the same order, so it is probably redundant: remove it if the colours match). Remove any `message=FALSE` that only existed for the scale message.
- [ ] **Step 4: Re-export.** Through the step runner, in dependency order: `export_native_plots.R`, `native-traits-export.R`, `ordination-export.R`, `remaining-plots-export.R`, `covariate-effect-export.R`, `prediction-export.R` (read each header for its arguments; `prediction-build.R` fits and is NOT run). Then every verifier: `verify_lesson.R`, `verify_outputs.R`, `verify_diagnostics.R`, `verify_latent_tables.R`, `verify_native_plots.R`, `native-traits-verify.R`, `ordination-verify.R`, `remaining-plots-verify.R`, `covariate-effect-verify.R`, `prediction-verify.R`, `unbalanced-verify.R`, `test_lesson.R`, `test_source_check.R`, `dev/simstudy/jsdm-package-comparison/verify-teaching.R archive dev/simstudy/jsdm-package-comparison`, `dev/simstudy/jsdm-package-comparison/extension/verify.R /Users/douglasyu/src/occJSDM/dev/simstudy/results/lesson-4-extension-20260923 .`, `dev/simstudy/spatial-design-sweep/verify-lesson.R .`. All exit 0. If a summary bundle (for example `output-lesson.rds`) must be rebuilt, rebuild it only with its existing summarise script and say why.
- [ ] **Step 5: Compare.** Every fit md5 unchanged. Look at every re-exported PNG (Read tool) against its previous version: only the trace plots' chain colours (Okabe-Ito now) and nothing else may differ; colours in the four primer figures must be identical to before (the package now draws the same Okabe-Ito colours the manual scales did).
- [ ] **Step 6: Render Lessons 3 and 4** (both commands), PNG policy, and confirm "Scale for colour is already present" no longer appears in `vignettes/occJSDM-lesson-3.md`. Update the controller's frozen-chunk tool if any chunk it pins changed (`/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-VERIFY-PLAN/check-frozen-chunks.R`).
- [ ] **Step 7: Commit** "Re-export the lesson figures with the package's Okabe-Ito colours".

### Task 3: Colour-blind-safe figures in Lessons 2, 5, 6 and 7

**Files:** `vignettes/occJSDM-lesson-2.Rmd`, `-5.Rmd`, `-6.Rmd`, `-7.Rmd`, their `.md`, and their figures (Lessons 5 to 7 write to `vignettes/teaching-data/lesson-N-*`).

The audit from the port (Task 6 report, `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-PORT-PLAN/task-6-report.md`): Lessons 2 and 7 use purple `#7B3294` beside Okabe-Ito blue and vermillion (Lesson 7's fit-recovery figure is most at risk); Lesson 5 uses a custom 4-colour palette, and its errors-by-species and two gradient figures rely on colour alone; Lesson 6 has six figures on ggplot2's default palette with 2 or 3 groups.

- [ ] **Step 1:** Re-audit each lesson's figures (colour and fill aesthetics, manual values, defaults). Replace every categorical palette with Okabe-Ito in a consistent order (`#0072B2`, `#E69F00`, `#009E73`, `#CC79A7`, `#56B4E9`, `#D55E00`, `#F0E442`, `#000000`), keeping any meaning a lesson attaches to a colour consistent across that lesson. Where a figure relies on colour alone to separate groups that matter, also add a second cue (shape or linetype) if it fits without clutter. Continuous scales: use viridis if a scale is red-green; leave others.
- [ ] **Step 2:** Before editing a chunk, check whether a verifier compares or runs it (`/usr/bin/grep` the chunk label in `dev/simstudy/`); if so, stop and report rather than edit it.
- [ ] **Step 3:** Update every prose sentence and caption that names a colour.
- [ ] **Step 4:** Render the four lessons (both commands), PNG policy, look at every changed figure, run the lesson-specific verifiers (`verify-teaching.R`, `extension/verify.R`, the sweep's `verify-lesson.R`) and the global checks.
- [ ] **Step 5: Commit** "Use colour-blind-safe Okabe-Ito colours in Lessons 2, 5, 6 and 7".

### Task 4: Final verification, HTML and PR body

- [ ] Run every verifier and test from Task 2 Step 4 again, plus `Rscript -e 'devtools::test()'`. Confirm every archive fit md5 is unchanged and `git diff --stat 471611b -- R src DESCRIPTION NAMESPACE` is empty.
- [ ] Update `vignettes/LESSON-PLAN.md` (dated log entry) and `TODO.md`. Remove or update any item this branch resolves, for example notes that the verifiers fail on main, and say in the log entry that the archive library was reinstalled and why.
- [ ] Render HTML for Lessons 2 to 7 into `vignettes/` (gitignored).
- [ ] Draft `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-VERIFY-PLAN/pr-body.md`. It covers:
  - the fingerprint design and the allow list;
  - that no fit changed;
  - the library reinstall;
  - the re-exported bundles;
  - the colour changes per lesson;
  - the checks run.

  It ends with `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.

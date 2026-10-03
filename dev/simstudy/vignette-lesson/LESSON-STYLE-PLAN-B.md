# Lesson Rewrite Phase B Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restructure Lessons 3 and 4, move new-site prediction and the ten-community extension into lessons of their own, make the package's own plot with truth the main figure of each Lesson 3 section, and renumber the eight lessons 0 to 7 once, with every link, registration, figure prefix and living document updated, in one pull request.

**Architecture:** Every lesson is an R Markdown file under `vignettes/` that renders from committed teaching bundles without fitting anything. Task 1 adds four groups of package plots to Lesson 3's evidence by extending four existing exporters, re-running them on the three-sample full-fit archive and adding verifier checks of the existing standard. Tasks 2 to 5 split and restructure Lessons 3 and 4 under provisional file names, moving blocks unchanged wherever a verifier compares displayed code with exported code. Tasks 6 and 7 renumber once: first the files, figure prefixes, titles, registrations and script paths, then every live cross-reference and living document. Task 8 runs every check and drafts the pull request for the controller.

**Tech Stack:** R 4.5.0 with rmarkdown and knitr, dplyr, tidyr, tibble, ggplot2, ggtern and png; the occJSDM library installed in the full-fit archive; the committed bundles in `vignettes/teaching-data/`; the verifiers in `dev/simstudy/vignette-lesson/` and `dev/simstudy/jsdm-package-comparison/`; Node for the site test; perl for one mechanical link rewrite; git.

**Spec:** `dev/simstudy/vignette-lesson/LESSON-STYLE-SPEC.md`: "Sequence" (Phase B), "Lesson 3", "Lesson 4", "The intuition lesson", "Decisions taken on 2 October 2026", "What does not change" and "Verification for every lesson pull request". The owner's decisions of 3 October 2026 in `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/notes.md` override the spec where they differ: Phase B starts now, before PR #14 merges and before the beta tag; eight lessons numbered 0 to 7; path-only edits to verifiers and `finish-study.py` for the renumbering are approved; Lesson 3's four figure groups (environmental gradient 2, perfect-observation versions, `plotVariancePartitioning()`, the unaligned ordination plot) are re-exported now with new exporter code and matching verifier checks, no existing check loosened or removed. Research maps in the same directory: `lesson-3-map.md`, `lesson-4-map.md`, `renumber-inventory.md`. Format model: `dev/simstudy/vignette-lesson/LESSON-STYLE-PLAN-A.md`.

## Final numbering

The owner confirmed this mapping on 3 October 2026. Titles marked new are set by this plan.

- Lesson 0, `occJSDM-lesson-0` (unchanged): "Lesson 0 (optional): Create and explore a simulated survey".
- Lesson 1, `occJSDM-lesson-1`, from `occJSDM-lesson-intuition`: "Lesson 1: What occupancy models and joint species distribution models do" (new prefix).
- Lesson 2, `occJSDM-lesson-2`, from `occJSDM-lesson-1`: "Lesson 2: Fit the model and compare its answers with truth".
- Lesson 3, `occJSDM-lesson-3` (stem unchanged, prediction removed): "Lesson 3: Understand the model's outputs by comparing them with truth".
- Lesson 4, `occJSDM-lesson-4`, from the provisional `occJSDM-lesson-prediction` (split from old Lesson 3): "Lesson 4: Predict occupancy at new sites and compare models" (new).
- Lesson 5, `occJSDM-lesson-5`, from `occJSDM-lesson-4` sections 1 to 11: "Lesson 5: Compare four JSDMs with a community whose truth we know".
- Lesson 6, `occJSDM-lesson-6`, from the provisional `occJSDM-lesson-extension` (old Lesson 4 sections 12 to 15 and the reduced full results): "Lesson 6: Repeat the four-JSDM comparison across ten communities, with traits" (new).
- Lesson 7, `occJSDM-lesson-7`, from `occJSDM-lesson-2`: "Lesson 7: Spatial landscapes and survey design".

## Decisions this plan takes

Each is recorded here so the reviewer and the owner can see it; the final report lists the ones the owner may want to revisit.

- **Perfect-observation versions.** The package plots re-exported for the perfect-observation fit are `plotOccupancyCovariates()` for both gradients and `plotTraitsCoefficients()` for both gradients, because the environmental and trait sections compare the two fits and their custom two-fit figures are then duplicates. `plotVariancePartitioning()` is exported for the PCR fit only: the triangle does not label species and, in a non-spatial model, puts every point on one edge, so the custom per-species scatter `variation-truth` stays as the two-fit comparison.
- **Custom figures left in Lesson 3's main text:** `ordination-contribution` (no package plot draws the rotation-invariant combined contribution) and `variation-truth` (above). Each gets one sentence saying why it is custom. Removed as duplicates of a package plot with truth: `environmental-coefficients`, `response-profile-1`, `response-profile-2`, `trait-effects`, `correlation-truth`, `collection-effects`, `detection-effort`. The `baseline-probabilities` table stays (a table, not a figure).
- **Lesson 4 section 3** keeps its heading and number with the two-sentence reminder (Option A in `lesson-4-map.md`), so sections 4 to 9 keep their numbers; old section 11 becomes section 10 once section 10 moves to the appendix.
- **The extension lesson** numbers its sections 1 to 4 (old 12 to 15) and takes its final figure prefix `teaching-data/lesson-6-` when it is created in Task 4, so its twelve PNGs are renamed once.
- **Provisional stems** `occJSDM-lesson-prediction` and `occJSDM-lesson-extension` are registered unpublished in Tasks 2 and 4, as the intuition lesson was in Phase A, and take their numbers in Task 6.
- **Data file names keep `lesson-4-`** (`lesson-4-species-errors.csv`, `lesson-4-extension.rds`, `lesson-4-calibration.rds`, `lesson-4-calibration-tables.R`), per the planning ruling; only figure prefixes change.
- **The dev mirror** `dev/simstudy/vignette-lesson/unbalanced-lesson-1.Rmd` is renamed `unbalanced-lesson-2.Rmd` with the path-only verifier edit, so its name matches the lesson it mirrors. `unbalanced-lesson-0.Rmd` keeps its name.
- **Historical documents keep their numbers:** completed plans and reports (`LESSON-STYLE-PLAN-A.md`, `style-diagnostic/`, `dev/simstudy/vignette-lesson/PLAN.md` and `DESIGN.md`, `dev/simstudy/lesson-site/PLAN.md` and `DESIGN.md`, `dev/simstudy/jsdm-package-comparison/` reports and plans other than the two READMEs and `CALIBRATION.md`, `dev/simstudy/convergence-flag-diagnosis/REPORT.md`), the dated decision log in `vignettes/LESSON-PLAN.md`, struck-through DONE items in `TODO.md`, dated records inside `dev/simstudy/vignette-lesson/README.md`, `LESSON-STYLE-SPEC.md` itself, and all of `dev/simstudy/spatial-design-sweep/`. The one exception is a link target: `FITTING-REPORT.md` line 3 links to `occJSDM-lesson-4.md`, which would silently point at the prediction lesson, so its target (not its prose) changes to `occJSDM-lesson-5.md`.
- **Snippet files are frozen after Task 1.** The four canonical snippet files' md5 sums are stamped in their bundles, so no later task edits them, even where their prose names an old lesson number (`ordination-examples.Rmd` line 23 says "Lesson 1"). Only the lesson copy of such prose is renumbered.
- **Deferred to Phase C** (not done here): the teaching and copyedit prose passes, the prerequisites block, R-gloss cuts, the `ggtern::theme_bw()` explanation, the Lesson 4 sjSDM v0.1.0 citation and `SJSDM_MOJO_BACKEND` line, the line 765 causal claim, renaming `calibration_targets` in `calibration-combined-targets`, and trimming Lesson 3 to the spec's 1,150-line target. Expected Phase B lengths: Lesson 3 about 1,650 to 1,750 lines (about 1,250 before the appendix), the prediction lesson about 380, Lesson 5 about 860, Lesson 6 about 700.
- **Memory notes** are updated by the controller after the pull request opens, not by an implementer (Task 8, controller steps).

## Global Constraints

- Work only in the worktree `/Users/douglasyu/src/occJSDM-worktrees/lesson-phase-b`, on branch `codex/lesson-phase-b` (created from main at 0f48dae). Never work in `/Users/douglasyu/src/occJSDM`, never commit to main, never push (the controller pushes), never force anything (no `push --force`, no `reset --hard` of committed work, no `-f` on `git mv`). Reading the Lesson 4 extension archive under `/Users/douglasyu/src/occJSDM/dev/simstudy/results/` for its verifier is read-only use; nothing is written there.
- The full-fit archive is `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples`. Every exporter and archive verifier runs from the worktree root with `R_LIBS=/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples/library`, through the step runner below so that each run is logged in the archive's `BUILD-LOG.md`. Agent shells keep no variables between calls, so commands below spell paths out in full.
- Commit messages end with a blank line and `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. The pull request body (drafted in Task 8, opened by the controller) ends with `🤖 Generated with [Claude Code](https://claude.com/claude-code)` and carries a "needs Alex" list containing only package facts the lessons newly state that the spec marks unconfirmed.
- Markdown rules for every `.md` and `.Rmd` touched, from `CLAUDE.md` and `AGENTS.md` lines 89 to 124: one line per paragraph, no hard wrapping; no em-dashes anywhere, including commit messages (use `--`); no pipe tables in lessons or in `vignettes/LESSON-PLAN.md` (bullet lists instead); escape a literal `>` `<` `~` `|` in prose as `\>` `\<` `\~` `\|`, never inside code spans or fenced blocks; keep inline code spans short; never add an `editor_options` block. Re-read a file immediately before and after editing it (it may be open in RStudio); if an edit looks corrupted or duplicated, rewrite the whole file with Write.
- Out of scope unless a task names it: package code (`R/`, `src/`, `DESCRIPTION`, `NAMESPACE`); everything under `dev/simstudy/spatial-design-sweep/`; the Lesson 2 (spatial) and Lesson 4 (comparison and extension) bundles and their exporters; any verifier change other than the approved path edits (Tasks 2, 4, 6) and the Task 1 additions.
- Line numbers in this plan are from main at 0f48dae and from the research maps. Earlier tasks shift them. Executors find the quoted heading, chunk label or sentence, and never trust a line number alone.
- Shell is zsh and `grep` and `ls` are aliased: use `/usr/bin/grep` and `/bin/ls` in commands.
- Render commands, run from the worktree root (replace `<stem>`). The four `library()` calls keep R 4.5.0's "built under" warnings out of the rendered `.md`:

```sh
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/<stem>.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/<stem>.Rmd", output_format="rmarkdown::html_vignette")'
```

  Both must finish without a new warning in the console or in the `.md`. The HTML output is gitignored and never committed.

- PNG policy (as Plan A, made mechanical). A render rewrites every figure it draws. Keep a changed PNG only when its content changed; revert the others with the shared revert script below, which compares pixels with the index copy. Use `git mv` for every moved or renamed figure so history follows. Never apply the revert script to exporter-written PNGs in `vignettes/teaching-data/` (`native-plot-*`, `native-traits-*`, `remaining-plots-*`, `ordination-*`): their md5 sums are stamped in their bundles and must be exactly what the exporter wrote. Check with `git status --short vignettes/` before every commit.
- Checks that must pass at the end of every task from Task 2 onward:

```sh
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
node --test dev/simstudy/lesson-site/test_lessons.js
git diff --name-only 0f48dae -- ':!*.png' ':!*.rds' | xargs /usr/bin/grep -n $'\u2014' | wc -l
git grep -n "^|" -- 'vignettes/occJSDM*.Rmd' | wc -l
```

  The `$'\u2014'` argument is zsh's spelling of the em-dash character, so this plan never contains it. Expected: the link test ends "Shared lesson-link regression checks passed."; the site test reports `fail 0`; the em-dash count is 0; the pipe-table count falls as tasks convert tables (Lesson 3's 32 lines in Task 3, Lesson 4's in Tasks 4 and 5) and is 0 from Task 5 on.

- Stop and report to the controller, without working around it, if: a baseline verifier fails before this plan has changed anything; a package plotting call errors on the perfect-observation fit or the truth layer cannot be drawn on the ternary plot; an exporter or verifier would need a change outside its task's list; a verifier check would have to be loosened or removed; a render needs a bundle change; or the spec and this plan conflict.

---

## Shared tools (created in Task 1 Step 1, kept outside the repository)

These three files live in `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/`, the controller's working directory. They are not committed.

`run-step.zsh`, the logged step runner (the archive's own `run-step.sh` points at a removed worktree):

````zsh
#!/bin/zsh
# Usage: zsh run-step.zsh STEPNAME command args...
# Runs from the Phase B worktree with the archive library first in R_LIBS,
# logs to the archive's logs/STEPNAME.log and appends a record to BUILD-LOG.md.
ARCHIVE=/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
WT=/Users/douglasyu/src/occJSDM-worktrees/lesson-phase-b
export R_LIBS=$ARCHIVE/library
step=$1; shift
cd $WT || exit 1
log=$ARCHIVE/logs/$step.log
start=$(date +%s)
"$@" > $log 2>&1
rc=$?
end=$(date +%s)
{
  echo ""
  echo "## $step"
  echo ""
  echo "Started $(date -r $start '+%Y-%m-%d %H:%M:%S'), exit $rc, $((end-start)) s wall, HEAD $(git rev-parse --short HEAD), worktree lesson-phase-b"
  echo ""
  echo '```'
  echo "R_LIBS=$R_LIBS $*"
  echo '```'
  echo ""
  echo '```'
  tail -n 25 $log
  echo '```'
} >> $ARCHIVE/BUILD-LOG.md
echo "EXIT $rc" >> $log
tail -n 6 $log
exit $rc
````

`revert-unchanged-png.zsh`, the PNG policy (one pathspec argument, never `vignettes/teaching-data/` exporter outputs):

```zsh
#!/bin/zsh
# Usage: zsh revert-unchanged-png.zsh PATHSPEC
# Reverts each modified tracked PNG under PATHSPEC whose pixels equal the index copy.
cd /Users/douglasyu/src/occJSDM-worktrees/lesson-phase-b || exit 1
tmp=$(mktemp -d)
for f in $(git diff --name-only -- "$1"); do
  [[ $f == *.png ]] || continue
  git show ":$f" > $tmp/old.png 2>/dev/null || continue
  if Rscript -e 'a <- png::readPNG(commandArgs(TRUE)[1]); b <- png::readPNG(commandArgs(TRUE)[2]); quit(status = if (identical(a, b)) 0 else 1)' $tmp/old.png $f; then
    git checkout -- $f; echo "reverted, pixels unchanged: $f"
  else
    echo "kept, content changed: $f"
  fi
done
rm -rf $tmp
```

`check-lesson-links.R`, which checks that every lesson link or file name in the files given names an existing lesson source and, if it has a fragment, an existing heading; for lesson sources it also checks same-file `](#anchor)` links. Historical text is skipped: everything from `### Decisions made after this update` on (the dated log of `vignettes/LESSON-PLAN.md`) and any line containing a struck-through `~~` item (the DONE items of `TODO.md`). Run from the worktree root: `Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/check-lesson-links.R FILE...`.

```r
slug <- function(heading) {
  text <- sub("^#+\\s+", "", heading)
  text <- sub("\\s*\\{[^}]*\\}\\s*$", "", text)
  text <- tolower(gsub("`", "", text))
  text <- gsub("[^a-z0-9 _-]", "", text)
  gsub(" ", "-", text)
}
anchors <- function(rmd) {
  lines <- readLines(rmd, warn = FALSE)
  inside <- cumsum(grepl("^```", lines)) %% 2 == 1
  vapply(lines[!inside & grepl("^#{1,6} ", lines)], slug, character(1), USE.NAMES = FALSE)
}
bad <- character()
files <- commandArgs(TRUE)
for (file in files) {
  lines <- readLines(file, warn = FALSE)
  log_start <- grep("^### Decisions made after this update", lines)
  if (length(log_start)) lines <- lines[seq_len(log_start[1] - 1L)]
  for (i in seq_along(lines)) {
    if (grepl("~~", lines[i], fixed = TRUE)) next
    hits <- regmatches(lines[i], gregexpr("occJSDM(-lesson-[a-z0-9]+)?\\.(md|Rmd)(#[A-Za-z0-9_-]+)?", lines[i]))[[1]]
    for (hit in hits) {
      parts <- strsplit(hit, "#", fixed = TRUE)[[1]]
      target <- file.path("vignettes", sub("\\.md$", ".Rmd", parts[1]))
      if (!file.exists(target)) {
        bad <- c(bad, sprintf("%s:%d no such lesson: %s", file, i, parts[1])); next
      }
      if (length(parts) == 2 && !(parts[2] %in% anchors(target)))
        bad <- c(bad, sprintf("%s:%d no such heading: %s", file, i, hit))
    }
    if (grepl("^vignettes/occJSDM.*\\.Rmd$", file)) {
      local <- regmatches(lines[i], gregexpr("\\]\\(#[A-Za-z0-9_-]+\\)", lines[i]))[[1]]
      for (hit in sub("^\\]\\(#(.*)\\)$", "\\1", local))
        if (!(hit %in% anchors(file))) bad <- c(bad, sprintf("%s:%d no such heading here: #%s", file, i, hit))
    }
  }
}
if (length(bad)) { writeLines(bad); quit(status = 1) }
cat("Every lesson link in", length(files), "files names an existing lesson and heading.\n")
```

---

## File structure

- `dev/simstudy/vignette-lesson/{native-plot-examples.Rmd,export_native_plots.R,verify_native_plots.R}`, `{native-traits-examples.Rmd,native-traits-export.R,native-traits-verify.R}`, `{remaining-plots-examples.Rmd,remaining-plots-export.R,remaining-plots-verify.R}`, `{ordination-examples.Rmd,ordination-export.R,ordination-verify.R}`: Task 1 extends them; their bundles and PNGs in `vignettes/teaching-data/` are regenerated.
- `vignettes/occJSDM-lesson-3.Rmd`, `.md`, `occJSDM-lesson-3_files/`: Task 1 (new chunks in the old gallery), Task 2 (prediction out), Task 3 (restructure).
- `vignettes/occJSDM-lesson-prediction.Rmd` (Task 2), renamed `occJSDM-lesson-4.Rmd` in Task 6.
- `vignettes/occJSDM-lesson-4.Rmd` (Tasks 4 and 5), renamed `occJSDM-lesson-5.Rmd` in Task 6.
- `vignettes/occJSDM-lesson-extension.Rmd` (Task 4), renamed `occJSDM-lesson-6.Rmd` in Task 6.
- Registrations, changed in Tasks 2, 4 and 6: `vignettes/lesson-links.R`, `dev/simstudy/vignette-lesson/test_lesson_links.R`, `_config.yml`, `assets/js/lessons.js`, `dev/simstudy/lesson-site/test_lessons.js`, `.Rbuildignore`, and in Task 6 `dev/simstudy/lesson-site/build_preview.R`.
- Path-only script edits: `dev/simstudy/vignette-lesson/prediction-verify.R` (Tasks 2 and 6), `dev/simstudy/vignette-lesson/unbalanced-verify.R`, `dev/simstudy/jsdm-package-comparison/verify-teaching.R`, `dev/simstudy/jsdm-package-comparison/extension/finish-study.py` (Task 6).
- Living documents, Task 7: `vignettes/LESSON-PLAN.md`, `TODO.md` (also Task 3), `AGENTS.md`, `dev/simstudy/vignette-lesson/README.md` (also Task 1), `dev/simstudy/vignette-lesson/remaining-plots-README.md` (Task 3), `dev/simstudy/jsdm-package-comparison/README.md`, `extension/README.md`, `extension/CALIBRATION.md`, `FITTING-REPORT.md` line 3. `README.md` and the quickstart `vignettes/occJSDM.Rmd` carry no lesson number or link (inventory 1.1 and 1.5) and are checked, not edited.

## Known collision with PR #14 (open, `codex/site-waic`)

PR #14 edits Lesson 3 lines 1066 to 1083 plus 16 new lines after 1083 (all inside `### Check the additional fit and understand the WAIC limitation`, which Task 2 moves to the prediction lesson), Lesson 3 line 1601 (the WAIC paragraph in `### Summaries and what to do next`, which Task 2 edits), `prediction-export.R` line 57 and `prediction-verify.R` line 46 (Tasks 2 and 6 edit line 263 of the same file), `vignettes/LESSON-PLAN.md` line 64 and `TODO.md` (Task 7). Git does not follow a cross-file move, so whichever branch merges second re-applies the other's WAIC hunks by hand in `vignettes/occJSDM-lesson-4.Rmd` (the final prediction lesson) and in Lesson 3's diagnostics summary, re-renders Lessons 3 and 4, and re-runs `prediction-verify.R`. The pull request body records this (Task 8).

---

### Task 1: Re-export Lesson 3's package plots for gradient 2, the perfect-observation fit, variation partitioning and the ordinary biplot

**Files:**
- Modify: `dev/simstudy/vignette-lesson/native-plot-examples.Rmd`, `export_native_plots.R`, `verify_native_plots.R`
- Modify: `dev/simstudy/vignette-lesson/native-traits-examples.Rmd`, `native-traits-export.R`, `native-traits-verify.R`
- Modify: `dev/simstudy/vignette-lesson/remaining-plots-examples.Rmd`, `remaining-plots-export.R`, `remaining-plots-verify.R`
- Modify: `dev/simstudy/vignette-lesson/ordination-examples.Rmd`, `ordination-export.R`, `ordination-verify.R`
- Regenerate (exporter output, committed exactly as written): `vignettes/teaching-data/native-plots.rds`, `native-traits-data.rds`, `remaining-plots-data.rds`, `ordination-examples.rds` and the existing `native-plot-*.png` (9), `native-traits-gradient-*.png` (2), `remaining-plots-*.png` (5), `ordination-{sites,loadings,biplot}.png` (3)
- Create (exporter output): `vignettes/teaching-data/native-plot-environment-2.png`, `native-plot-environment-perfect-1.png`, `native-plot-environment-perfect-2.png`, `native-traits-perfect-gradient-1.png`, `native-traits-perfect-gradient-2.png`, `remaining-plots-variation.png`, `ordination-ordinary-biplot.png`
- Modify: `vignettes/occJSDM-lesson-3.Rmd` and `.md` (new chunks inserted in the existing gallery; `load-full-fit`), `dev/simstudy/vignette-lesson/README.md` (one dated record)
- Read: `lesson-3-map.md` section 3; the four exporters and verifiers in full before editing

**Interfaces:**
- Consumes: the archive fits `default-fit.rds` and `perfect-fit.rds` (manifests `lesson$manifests$default` and `$perfect` in `vignettes/teaching-data/nonspatial-lesson.rds`, md5-checked), `input.rds`, and `output-lesson.rds` (`coefficients`, `variation`), all unchanged.
- Produces for Task 3: locked chunks `native-environment-2-example`, `native-environment-perfect-1-example`, `native-environment-perfect-2-example` (objects `native_environment_2`, `native_environment_perfect_1`, `native_environment_perfect_2`; the perfect chunks use `environment_truth` from `native-environment-example` and `environment_2_truth`); `native-traits-perfect-gradient-1`, `native-traits-perfect-gradient-2` (objects `native_traits_perfect_1`, `native_traits_perfect_2`); `remaining-plots-variation` (object `remaining_variation`, truth table `remaining_examples$truth$variation` with columns `Species`, `Env`, `Biotic`, `Spatial`); unlocked image chunks `native-environment-2-image`, `native-environment-perfect-1-image`, `native-environment-perfect-2-image`, `native-traits-perfect-gradient-1-image`, `native-traits-perfect-gradient-2-image`, `remaining-plots-variation-image`, `ordination-ordinary-biplot-image`. The reader's perfect-observation fit object is `fitmodel_perfect`. Bundle additions: `provenance$perfect_fit_manifest` in `native-plots.rds` and `native-traits-data.rds`; `plots$environment_2`, `plots$environment_perfect_1`, `plots$environment_perfect_2`; `plots$perfect_gradient_1`, `plots$perfect_gradient_2`; `plots$variation`, `truth$variation`; `plots$ordinary_biplot`.
- Verifier chunk counts after this task: native 13 (was 10), traits 5 (was 3), remaining 7 (was 6), ordination 5 (unchanged).

- [ ] **Step 1: Confirm the branch and create the shared tools**

```sh
cd /Users/douglasyu/src/occJSDM-worktrees/lesson-phase-b
git status --short && git branch --show-current && git log --oneline -1
```

Expected: clean tree, `codex/lesson-phase-b`, HEAD at or after 0f48dae (the plan commit). Write the three files of "Shared tools" with the Write tool, then mark the start of Phase B in the archive log:

```sh
printf '\n## Phase B begins\n\nFrom here, steps run from worktree lesson-phase-b on branch codex/lesson-phase-b through .sdd/LESSON-STYLE-PLAN-B/run-step.zsh.\n' >> /Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples/BUILD-LOG.md
```

- [ ] **Step 2: Run the four existing verifiers as a baseline**

```sh
S=/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/run-step.zsh; A=/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
zsh $S 40-baseline-verify_native_plots Rscript dev/simstudy/vignette-lesson/verify_native_plots.R $A
zsh $S 41-baseline-native-traits-verify Rscript dev/simstudy/vignette-lesson/native-traits-verify.R $A
zsh $S 42-baseline-ordination-verify Rscript dev/simstudy/vignette-lesson/ordination-verify.R $A
zsh $S 43-baseline-remaining-plots-verify Rscript dev/simstudy/vignette-lesson/remaining-plots-verify.R $A
```

Expected final lines, in order: "Nine PNGs, compact export, full-fit hash, source hashes and exact teaching-code provenance verified."; "generating trait-scale conversion, plotted truth coordinates, code and PNG evidence."; the ordination `print(maximum_error)` after "marginal quantiles, exact teaching code, native plotted quantities and truth overlays."; "Five PNGs, compact bundle, full-fit/library/source hashes, and exact displayed code verified." Each log also contains "Lesson 3 displays the exact exported ...". If any fails, stop (Global Constraints).

- [ ] **Step 3: Trial the new calls without writing anything**

```sh
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/run-step.zsh 44-trial-new-plots Rscript -e 'suppressMessages({library(occJSDM); library(ggplot2)}); A <- "/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples"; fp <- readRDS(file.path(A, "perfect-fit.rds"))$fit; fd <- readRDS(file.path(A, "default-fit.rds"))$fit; stopifnot(identical(fp$X_psi, fd$X_psi)); print(dim(plotOccupancyCovariates(fp, covName = "X_psi.EnvCov.2")$data)); print(dim(plotTraitsCoefficients(fp, covName = "X_psi.EnvCov.1")$data)); tr <- data.frame(Env = c(.2, .5), Biotic = c(.8, .5), Spatial = 0); b <- ggplot_build(plotVariancePartitioning(fd) + geom_point(data = tr, aes(x = Env, y = Biotic, z = Spatial), inherit.aes = FALSE, shape = 4)); print(b$data[[2]][, c("x", "y", "z")]); print(length(ggplot_build(plotBiplot(fd))$data))'
```

Expected: two data dimensions printed without error, the truth layer printed with columns `x`, `y`, `z` holding 0.2/0.8/0 and 0.5/0.5/0, and `3` for the biplot layers. A warning "Ignoring unknown aesthetics: z" is expected and harmless (the ternary coordinates still read `z`, as the printed layer shows). If any call errors, stop.

- [ ] **Step 4: Add the three environmental chunks to `native-plot-examples.Rmd`**

Insert directly after the `native-environment-image` chunk and before the `native-collection-example` chunk, exactly:

````markdown
The second environmental gradient uses the same call with its own name.

```{r native-environment-2-example, eval=FALSE, purl=TRUE}
environment_2_truth <- native_truth$environment |>
  filter(covariate == "X_psi.EnvCov.2")

native_environment_2 <- plotOccupancyCovariates(
  fitmodel, covName = "X_psi.EnvCov.2"
) +
  native_theme +
  geom_point(
    data = environment_2_truth, aes(x = species, y = truth),
    inherit.aes = FALSE, shape = 4, size = 3, stroke = 1
  ) +
  labs(
    title = "Environmental gradient 2",
    y = "Effect on log-odds per standard deviation",
    caption = "Black cross: truth. Bar: native 95% interval. Red line: zero effect."
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

native_environment_2
```

```{r native-environment-2-image, echo=FALSE, purl=FALSE, out.width="100%"}
knitr::include_graphics("teaching-data/native-plot-environment-2.png")
```

The same calls on the perfect-observation fit, loaded as `fitmodel_perfect`, show what the occupancy model recovers when every presence and absence is known. The truth crosses are the same; only the fit changes.

```{r native-environment-perfect-1-example, eval=FALSE, purl=TRUE}
native_environment_perfect_1 <- plotOccupancyCovariates(
  fitmodel_perfect, covName = "X_psi.EnvCov.1"
) +
  native_theme +
  geom_point(
    data = environment_truth, aes(x = species, y = truth),
    inherit.aes = FALSE, shape = 4, size = 3, stroke = 1
  ) +
  labs(
    title = "Environmental gradient 1, perfect-observation fit",
    y = "Effect on log-odds per standard deviation",
    caption = "Black cross: truth. Bar: native 95% interval. Red line: zero effect."
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

native_environment_perfect_1
```

```{r native-environment-perfect-1-image, echo=FALSE, purl=FALSE, out.width="100%"}
knitr::include_graphics("teaching-data/native-plot-environment-perfect-1.png")
```

```{r native-environment-perfect-2-example, eval=FALSE, purl=TRUE}
native_environment_perfect_2 <- plotOccupancyCovariates(
  fitmodel_perfect, covName = "X_psi.EnvCov.2"
) +
  native_theme +
  geom_point(
    data = environment_2_truth, aes(x = species, y = truth),
    inherit.aes = FALSE, shape = 4, size = 3, stroke = 1
  ) +
  labs(
    title = "Environmental gradient 2, perfect-observation fit",
    y = "Effect on log-odds per standard deviation",
    caption = "Black cross: truth. Bar: native 95% interval. Red line: zero effect."
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

native_environment_perfect_2
```

```{r native-environment-perfect-2-image, echo=FALSE, purl=FALSE, out.width="100%"}
knitr::include_graphics("teaching-data/native-plot-environment-perfect-2.png")
```
````

- [ ] **Step 5: Extend `export_native_plots.R`**

After the `stopifnot(identical(saved$source_hashes, input$source_hashes), ...)` block that ends with `c(6000L, 4L)))`, insert:

```r
perfect_manifest <- lesson$manifests$perfect
perfect_path <- file.path(archive, perfect_manifest$file)
stopifnot(identical(perfect_manifest, outputs$fit_manifests$perfect),
          identical(perfect_manifest$md5, md5(perfect_path)))
perfect_saved <- readRDS(perfect_path)
fitmodel_perfect <- perfect_saved$fit
validate_lesson_fit_identity(fitmodel_perfect, input)
stopifnot(identical(perfect_saved$source_hashes, input$source_hashes),
          identical(perfect_saved$input_md5, lesson$input_md5),
          identical(perfect_saved$mcmc, perfect_manifest$mcmc),
          identical(fitmodel_perfect$infos$speciesNames, fitmodel$infos$speciesNames),
          identical(fitmodel_perfect$X_psi, fitmodel$X_psi),
          identical(dim(fitmodel_perfect$results_output$jsdm_output$B_output),
                    c(2L, 10L, 6000L, 4L)))
```

After `environment$fitmodel <- fitmodel` add `environment$fitmodel_perfect <- fitmodel_perfect`. Replace the `plot_names` vector and the `height` vector so that the three new plots come last:

```r
plot_names <- c("environment", "collection", "occupancy_rates", "collection_rates",
                "primer_1", "primer_2", "correlations", "effort_k", "effort_m",
                "environment_2", "environment_perfect_1", "environment_perfect_2")
```

```r
  width = 9, height = c(4.8, 4.8, 5.2, 5.2, 5.5, 5.5, 7, 7.5, 9, 4.8, 4.8, 4.8)
```

In `provenance`, add `perfect_fit_manifest = perfect_manifest,` directly after `fit_manifest = manifest,`, and replace the note's first string with `"Native calls used the unchanged complete 6000 x 4 default fit; the two perfect-observation plots used the unchanged complete 6000 x 4 perfect fit."`.

- [ ] **Step 6: Add the two perfect-observation trait chunks and extend `native-traits-export.R`**

In `native-traits-examples.Rmd`, insert after the `native-traits-gradient-2-image` chunk and before the paragraph "Read each bar against both references.":

````markdown
The same two calls on the perfect-observation fit, `fitmodel_perfect`, use the same generating crosses.

```{r native-traits-perfect-gradient-1, eval=FALSE, purl=TRUE}
native_traits_perfect_1 <- occJSDM::plotTraitsCoefficients(
  fitmodel_perfect, covName = "X_psi.EnvCov.1"
) +
  native_trait_theme +
  geom_point(
    data = filter(native_trait_truth, covariate == "X_psi.EnvCov.1"),
    aes(x = trait, y = truth), inherit.aes = FALSE,
    shape = 4, size = 3, stroke = 1
  ) +
  labs(
    title = "Trait effects on the response to gradient 1, perfect-observation fit",
    x = "Measured trait",
    y = "Change in environmental coefficient\nper trait standard deviation",
    caption = "Black cross: generating effect on the fitted scale. Bar: native 95% interval."
  )

native_traits_perfect_1
```

```{r native-traits-perfect-gradient-1-image, echo=FALSE, purl=FALSE, out.width="100%"}
knitr::include_graphics("teaching-data/native-traits-perfect-gradient-1.png")
```

```{r native-traits-perfect-gradient-2, eval=FALSE, purl=TRUE}
native_traits_perfect_2 <- occJSDM::plotTraitsCoefficients(
  fitmodel_perfect, covName = "X_psi.EnvCov.2"
) +
  native_trait_theme +
  geom_point(
    data = filter(native_trait_truth, covariate == "X_psi.EnvCov.2"),
    aes(x = trait, y = truth), inherit.aes = FALSE,
    shape = 4, size = 3, stroke = 1
  ) +
  labs(
    title = "Trait effects on the response to gradient 2, perfect-observation fit",
    x = "Measured trait",
    y = "Change in environmental coefficient\nper trait standard deviation",
    caption = "Black cross: generating effect on the fitted scale. Bar: native 95% interval."
  )

native_traits_perfect_2
```

```{r native-traits-perfect-gradient-2-image, echo=FALSE, purl=FALSE, out.width="100%"}
knitr::include_graphics("teaching-data/native-traits-perfect-gradient-2.png")
```
````

In `native-traits-export.R`, after the `stopifnot(identical(saved$source_hashes, ...` block that ends `< 1e-12)`, insert:

```r
perfect_manifest <- lesson$manifests$perfect
perfect_path <- file.path(archive, perfect_manifest$file)
stopifnot(identical(perfect_manifest, outputs$fit_manifests$perfect),
          identical(perfect_manifest$md5, md5(perfect_path)))
perfect_saved <- readRDS(perfect_path)
fitmodel_perfect <- perfect_saved$fit
validate_lesson_fit_identity(fitmodel_perfect, input)
stopifnot(identical(perfect_saved$source_hashes, input$source_hashes),
          identical(perfect_saved$input_md5, lesson$input_md5),
          identical(perfect_saved$mcmc, perfect_manifest$mcmc),
          identical(fitmodel_perfect$Tr, fitmodel$Tr),
          identical(fitmodel_perfect$X_psi, fitmodel$X_psi),
          identical(dim(fitmodel_perfect$results_output$jsdm_output$G_output),
                    c(2L, 2L, 6000L, 4L)))
```

After `student_environment$outputs <- outputs` add `student_environment$fitmodel_perfect <- fitmodel_perfect`; after `stopifnot(identical(student_environment$fitmodel, saved$fit))` add `stopifnot(identical(student_environment$fitmodel_perfect, perfect_saved$fit))`. Replace the `plots` list and the `figures` file vector:

```r
plots <- list(gradient_1 = student_environment$native_traits_1,
              gradient_2 = student_environment$native_traits_2,
              perfect_gradient_1 = student_environment$native_traits_perfect_1,
              perfect_gradient_2 = student_environment$native_traits_perfect_2)
figures <- data.frame(plot = names(plots),
                      file = c("native-traits-gradient-1.png", "native-traits-gradient-2.png",
                               "native-traits-perfect-gradient-1.png",
                               "native-traits-perfect-gradient-2.png"),
                      width = 8, height = 4.8)
```

In `provenance`, add `perfect_fit_manifest = perfect_manifest,` after `fit_manifest = manifest,`; change the note to `"Native trait intervals use all 6000 x 4 draws from the unchanged default fit, and from the unchanged perfect fit for the two perfect-observation plots."`; change the final `cat()` to `cat("Exported four native trait plots from all 24000 retained draws of each fit.\n")`.

- [ ] **Step 7: Add the variation-partitioning chunk and extend `remaining-plots-export.R`**

In `remaining-plots-examples.Rmd`, append at the end of the file:

````markdown
### Variation partitioning with truth

`plotVariancePartitioning()` draws each species as one point in a triangle whose corners are the environmental, spatial and residual (`Biotic`) shares, each a posterior mean. In this non-spatial model every spatial share is zero, so every point lies on the edge between `Env` and `Biotic`. Black crosses mark the true shares under the package's own definition. The plot keeps its own ternary theme, so no theme is added.

```{r remaining-plots-variation, eval=FALSE, purl=TRUE}
# ggplot2 may warn that it ignores z; the ternary coordinates still use it.
remaining_variation <- plotVariancePartitioning(fitmodel) +
  geom_point(
    data = remaining_truth$variation,
    aes(x = Env, y = Biotic, z = Spatial),
    inherit.aes = FALSE, shape = 4, size = 3, stroke = 1
  ) +
  labs(
    title = "Variation partitioning, PCR-observation fit",
    caption = "Dot: posterior mean shares for one species. Black cross: true shares, same definition."
  )

remaining_variation
```

```{r remaining-plots-variation-image, echo=FALSE, purl=FALSE, out.width="100%"}
knitr::include_graphics("teaching-data/remaining-plots-variation.png")
```
````

In `remaining-plots-export.R`: append `"returnVariancePartitioning", "plotVariancePartitioning", "returnVariancePartitioningMatrix", "plotVarPart"` to `api_names`; before `remaining_examples <- list(truth = list(` insert

```r
variation_truth <- outputs$variation |>
  filter(arm == "default") |>
  select(species, component, truth) |>
  pivot_wider(names_from = component, values_from = truth) |>
  transmute(Species = species, Env = Environmental, Biotic = Residual, Spatial)
```

and add `variation = variation_truth` as the last element of that `truth` list. Replace `plot_names` and the `figures` tibble:

```r
plot_names <- c("gradient_1", "gradient_2", "stage1_fp", "stage2_fp", "detection", "variation")
```

```r
figures <- tibble(plot = plot_names,
                  file = paste0("remaining-plots-", gsub("_", "-", plot_names), ".png"),
                  width = c(10, 10, 10, 10, 10, 8), height = c(8, 8, 5.5, 6, 6, 7))
```

Change the last `cat()` to `cat("Exported six native plots, matching truth, scaling, and exact-code provenance.\n")`.

- [ ] **Step 8: Add the ordinary-biplot image and extend `ordination-export.R`**

In `ordination-examples.Rmd`, insert directly after the `ordination-standard` chunk:

````markdown
The ordinary biplot uses the package's stored orientation. It carries no truth overlay: on unaligned axes a correct configuration can be rotated or reflected away from the generating one, so a mismatch would not show an error.

```{r ordination-ordinary-biplot-image, echo=FALSE, purl=FALSE, out.width="100%"}
knitr::include_graphics("teaching-data/ordination-ordinary-biplot.png")
```
````

In `ordination-export.R`, replace the `plots` and `figures` assignments:

```r
plots <- c(setNames(lapply(paste0("native_", plot_names), get,
                           envir = student_environment), plot_names),
           list(ordinary_biplot = student_environment$ordinary_biplot))
figures <- tibble(plot = names(plots),
                  file = paste0("ordination-", gsub("_", "-", names(plots)), ".png"),
                  width = c(12, 12, 9, 9), height = c(7, 7, 7, 7))
```

Add the string `"The ordinary biplot uses the unaligned fit and carries no truth overlay."` as the last element of the note's `paste()`, and change the final `cat()` to `cat("Exported three aligned native ordination plots from 24000 draws and the ordinary biplot.\n")`.

- [ ] **Step 9: Commit the snippets and exporters before running them**

The exporters stamp the md5 of the snippet and exporter files into their bundles, so these files must be final and committed before Step 10. If an exporter needs a fix after this commit, commit the fix first and rerun it.

```bash
git add dev/simstudy/vignette-lesson/native-plot-examples.Rmd dev/simstudy/vignette-lesson/export_native_plots.R dev/simstudy/vignette-lesson/native-traits-examples.Rmd dev/simstudy/vignette-lesson/native-traits-export.R dev/simstudy/vignette-lesson/remaining-plots-examples.Rmd dev/simstudy/vignette-lesson/remaining-plots-export.R dev/simstudy/vignette-lesson/ordination-examples.Rmd dev/simstudy/vignette-lesson/ordination-export.R
git commit -m "Extend the Lesson 3 plot exporters: gradient 2, perfect-observation fit, variation partitioning, ordinary biplot

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 10: Run the four exporters**

```sh
S=/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/run-step.zsh; A=/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
zsh $S 45-export_native_plots Rscript dev/simstudy/vignette-lesson/export_native_plots.R $A
zsh $S 46-native-traits-export Rscript dev/simstudy/vignette-lesson/native-traits-export.R $A
zsh $S 47-ordination-export Rscript dev/simstudy/vignette-lesson/ordination-export.R $A
zsh $S 48-remaining-plots-export Rscript dev/simstudy/vignette-lesson/remaining-plots-export.R $A
git status --short vignettes/teaching-data/
```

Expected: "Exported 12 native plots and compact truth/provenance tables."; "Exported four native trait plots from all 24000 retained draws of each fit."; "Exported three aligned native ordination plots from 24000 draws and the ordinary biplot."; "Exported six native plots, matching truth, scaling, and exact-code provenance."; seven new PNGs untracked and the four `.rds` bundles modified. Open each new PNG with the Read tool: the gradient plots show ten species bars with crosses; the trait plots show two traits; the triangle shows dots and crosses along the `Env` to `Biotic` edge; the biplot shows grey sites, blue arrows and labels.

- [ ] **Step 11: Add the new verifier checks**

Change only what is listed; every existing `stopifnot()` stays as it is, apart from the three count literals, which rise.

`verify_native_plots.R`: after the `stopifnot(identical(saved$source_hashes, provenance$source_hashes), ...)` block, insert

```r
stopifnot(identical(provenance$perfect_fit_manifest, lesson$manifests$perfect))
perfect_path <- file.path(archive, provenance$perfect_fit_manifest$file)
stopifnot(identical(md5(perfect_path), provenance$perfect_fit_manifest$md5))
perfect_saved <- readRDS(perfect_path)
fp <- perfect_saved$fit
validate_lesson_fit_identity(fp, input)
stopifnot(identical(perfect_saved$source_hashes, provenance$source_hashes),
          identical(perfect_saved$input_md5, provenance$input_md5),
          identical(perfect_saved$mcmc, provenance$perfect_fit_manifest$mcmc),
          identical(dim(fp$results_output$jsdm_output$B_output), c(2L, 10L, 6000L, 4L)))
```

change `stopifnot(length(expected_chunks) == 10L)` to `13L`; after `cat("All coefficient/rate intervals, standardized truth, species and primer axes verified.\n")` insert

```r
# Gradient 2 and both perfect-observation gradients, from their own draws.
close(fp$X_psi, f$X_psi)
stopifnot(identical(fp$infos$speciesNames, sp))
covariates <- colnames(f$X_psi)
for (name in c("environment_2", "environment_perfect_1", "environment_perfect_2")) {
  rec <- x$plots[[name]]
  k <- if (name == "environment_perfect_1") 1L else 2L
  draws <- if (name == "environment_2") j$B_output else fp$results_output$jsdm_output$B_output
  table <- x$truth$environment[x$truth$environment$covariate == covariates[k], ]
  stopifnot(nrow(rec$data) == length(sp), nrow(table) == length(sp))
  check_crosses(rec, table, 3)
  for (i in seq_len(nrow(rec$data))) {
    s <- match(rec$data$Output[i], sp)
    close(unlist(rec$data[i, c("2.5%", "97.5%")]), quantile(draws[k, s, , ], c(.025, .975)))
    close(table$truth[table$species == sp[s]], b$B[k, s])
  }
}
cat("Gradient 2 and both perfect-observation gradients: 30 intervals from all 24000 draws and 30 truth crosses verified.\n")
```

and change the last message's "Nine PNGs" to "Twelve PNGs".

`native-traits-verify.R`: after the `stopifnot(... G_output), c(2L, 2L, 6000L, 4L)))` block, insert

```r
stopifnot(identical(provenance$perfect_fit_manifest, lesson$manifests$perfect))
perfect_path <- file.path(archive, provenance$perfect_fit_manifest$file)
stopifnot(identical(md5(perfect_path), provenance$perfect_fit_manifest$md5))
perfect_saved <- readRDS(perfect_path)
fp <- perfect_saved$fit
validate_lesson_fit_identity(fp, input)
stopifnot(identical(perfect_saved$source_hashes, provenance$source_hashes),
          identical(perfect_saved$input_md5, provenance$input_md5),
          identical(perfect_saved$mcmc, provenance$perfect_fit_manifest$mcmc),
          identical(dim(fp$results_output$jsdm_output$G_output), c(2L, 2L, 6000L, 4L)))
close(fp$Tr, fit$Tr)
```

Change `stopifnot(length(expected) == 3L)` to `5L`. After the closing brace of the existing `for (environment in seq_along(environments)) {` loop, insert

```r
perfect_draws <- fp$results_output$jsdm_output$G_output
for (environment in seq_along(environments)) {
  record <- examples$plots[[paste0("perfect_gradient_", environment)]]
  bars <- record$layers[[1]]
  crosses <- record$layers[[3]]
  stopifnot(nrow(bars) == length(traits), nrow(crosses) == length(traits),
            all(crosses$shape == 4), all(record$layers[[2]]$yintercept == 0),
            identical(sort(as.character(record$x_labels)), sort(traits)))
  for (row in seq_len(nrow(bars))) {
    trait <- match(record$x_labels[as.integer(bars$x[row])], traits)
    stopifnot(!is.na(trait))
    close(c(bars$ymin[row], bars$ymax[row]),
          quantile(perfect_draws[trait, environment, , ], c(.025, .975)))
  }
  for (row in seq_len(nrow(crosses))) {
    trait <- match(record$x_labels[as.integer(crosses$x[row])], traits)
    close(crosses$y[row], factors$G[trait, environment] * trait_sd[trait])
  }
}
cat("Perfect-observation trait plots: intervals from all 24000 perfect-fit draws and generating crosses verified.\n")
```

`remaining-plots-verify.R`: after the `for (name in p$api_names) {` loop add `stopifnot(all(c("plotVariancePartitioning", "plotVarPart") %in% p$api_names))`; change `stopifnot(length(expected_chunks) == 6L)` to `7L`; after `cat("All 50 rate intervals and truth crosses ...")` insert

```r
rec <- x$plots$variation
vp <- j$varPart_output
stopifnot(identical(dim(vp), c(length(sp), 4L, 6000L, 4L)),
          identical(as.character(rec$data$Species), sp),
          nrow(rec$layers[[1]]) == length(sp), nrow(rec$layers[[2]]) == length(sp),
          all(rec$layers[[2]]$shape == 4))
share <- function(k) vapply(seq_along(sp), function(s) mean(vp[s, k, , ]), numeric(1))
close(rec$data$Env, share(1)); close(rec$data$Spatial, share(2)); close(rec$data$Biotic, share(3))
close(rec$layers[[1]]$x, share(1)); close(rec$layers[[1]]$y, share(3)); close(rec$layers[[1]]$z, share(2))
# Independent truth: in this non-spatial model the package partition reduces to two shares.
true_env <- vapply(seq_along(sp), function(s) {
  measured <- b$B0[s] + drop(f$X_psi %*% b$B[, s])
  hidden <- drop(b$U %*% b$L[, s])
  ve <- sd(plogis(measured)); vh <- sd(plogis(hidden)); vb <- sd(plogis(measured + hidden))
  ce <- ve + max(vb - vh, 0); ch <- vh + max(vb - ve, 0)
  ce / (ce + ch)
}, numeric(1))
stopifnot(identical(x$truth$variation$Species, sp))
close(x$truth$variation$Env, true_env); close(x$truth$variation$Biotic, 1 - true_env)
close(x$truth$variation$Spatial, rep(0, length(sp)))
close(rec$layers[[2]]$x, true_env); close(rec$layers[[2]]$y, 1 - true_env)
close(rec$layers[[2]]$z, rep(0, length(sp)))
cat("Variation partitioning: ten posterior-mean shares from all 24000 draws and ten independent true shares verified in plotted ternary coordinates.\n")
```

and change the last message's "Five PNGs" to "Six PNGs".

`ordination-verify.R`: before the final `for (index in seq_len(nrow(examples$figures))) {` loop, insert

```r
ordinary <- examples$plots$ordinary_biplot
ordinary_sites <- apply(original$U_output, c(1, 2), median)
ordinary_loadings <- apply(original$L_output, c(1, 2), median)
ordinary_multiplier <- .8 * max(sqrt(rowSums(ordinary_sites^2))) /
  max(sqrt(colSums(ordinary_loadings^2)))
stopifnot(length(ordinary$layers) == 3L,
          identical(as.character(ordinary$layers[[3]]$label), species))
close(ordinary$layers[[1]]$x, ordinary_sites[, 1])
close(ordinary$layers[[1]]$y, ordinary_sites[, 2])
close(ordinary$layers[[2]]$x, rep(0, length(species)))
close(ordinary$layers[[2]]$xend, ordinary_loadings[1, ] * ordinary_multiplier)
close(ordinary$layers[[2]]$yend, ordinary_loadings[2, ] * ordinary_multiplier)
close(ordinary$layers[[3]]$x, ordinary_loadings[1, ] * ordinary_multiplier)
close(ordinary$layers[[3]]$y, ordinary_loadings[2, ] * ordinary_multiplier)
close(examples$ordinary_quantiles$sites[2, , ], ordinary_sites)
cat("Ordinary biplot: 100 unaligned site medians and ten loading arrows and labels verified from all 24000 draws; no truth layer.\n")
```

- [ ] **Step 12: Put the new chunks into Lesson 3's existing gallery**

The four verifiers require every locked chunk of a snippet file to be present in Lesson 3 with an identical body, so copy, do not retype. In `vignettes/occJSDM-lesson-3.Rmd`:
  - After the `native-environment-image` chunk, paste Step 4's inserted text exactly.
  - After the `native-traits-gradient-2-image` chunk, paste Step 6's inserted Rmd text exactly.
  - Before `## Reproduce the extraction or find a function`, paste Step 7's inserted Rmd text exactly (its `###` heading included).
  - After the `ordination-standard` chunk, paste Step 8's inserted Rmd text exactly.
  - In the `load-full-fit` chunk, after `fitmodel <- saved_fit$fit`, add the line `fitmodel_perfect <- readRDS("/path/to/full-fits/perfect-fit.rds")$fit`.

- [ ] **Step 13: Run the four verifiers**

```sh
S=/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/run-step.zsh; A=/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
zsh $S 49-verify_native_plots Rscript dev/simstudy/vignette-lesson/verify_native_plots.R $A
zsh $S 50-native-traits-verify Rscript dev/simstudy/vignette-lesson/native-traits-verify.R $A
zsh $S 51-ordination-verify Rscript dev/simstudy/vignette-lesson/ordination-verify.R $A
zsh $S 52-remaining-plots-verify Rscript dev/simstudy/vignette-lesson/remaining-plots-verify.R $A
/usr/bin/grep -h "Lesson 3 displays\|verified" /Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples/logs/{49,50,51,52}-*.log
```

Expected: all four exit 0; each log has its "Lesson 3 displays the exact exported ..." line, the new lines of Step 11, "Twelve PNGs, ..." and "Six PNGs, ...". A failure in a new check means a wrong layer index or a wrong draws array; fix the verifier (never the bundle) and rerun.

- [ ] **Step 14: Render Lesson 3 and record the build**

Render Lesson 3 with both commands. Open the rendered `.md` at the seven new images and confirm each shows. Apply the PNG policy to the lesson's own figures only: `zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/revert-unchanged-png.zsh vignettes/occJSDM-lesson-3_files` (expected: every changed PNG reverted, since no knit-time figure changed). In `dev/simstudy/vignette-lesson/README.md`, add this record as the last section of the file:

```markdown
### Lesson 3: package plots for gradient 2, the perfect-observation fit, variation partitioning and the ordinary biplot, added 3 October 2026

Phase B of the lesson rewrite makes the package's own plot with truth the main figure of each Lesson 3 section. Four exporters were extended and re-run on `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples` without refitting: `export_native_plots.R` adds `plotOccupancyCovariates()` for gradient 2 and for both gradients of the perfect-observation fit; `native-traits-export.R` adds `plotTraitsCoefficients()` for both gradients of the perfect-observation fit; `remaining-plots-export.R` adds `plotVariancePartitioning()` with the true shares; `ordination-export.R` saves the ordinary, unaligned biplot without a truth layer. The matching verifiers re-derive every new plotted interval or point from all 24,000 draws of the relevant fit, re-derive the truth independently, and still require the lesson's displayed code to equal the exported code. The runs are logged as steps 40 to 52 in that archive's `BUILD-LOG.md`.
```

- [ ] **Step 15: Commit the bundles, figures, verifiers and lesson**

```bash
git add vignettes/teaching-data/native-plots.rds vignettes/teaching-data/native-traits-data.rds vignettes/teaching-data/remaining-plots-data.rds vignettes/teaching-data/ordination-examples.rds vignettes/teaching-data/native-plot-*.png vignettes/teaching-data/native-traits-*.png vignettes/teaching-data/remaining-plots-*.png vignettes/teaching-data/ordination-*.png dev/simstudy/vignette-lesson/verify_native_plots.R dev/simstudy/vignette-lesson/native-traits-verify.R dev/simstudy/vignette-lesson/remaining-plots-verify.R dev/simstudy/vignette-lesson/ordination-verify.R dev/simstudy/vignette-lesson/README.md vignettes/occJSDM-lesson-3.Rmd vignettes/occJSDM-lesson-3.md vignettes/occJSDM-lesson-3_files
git status --short
git commit -m "Re-export Lesson 3's package plots with gradient 2, the perfect-observation fit, variation partitioning and the ordinary biplot

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

Expected `git status --short` after the commit: clean.

---

### Task 2: Move new-site prediction and model comparison into a provisional lesson

**Files:**
- Create: `vignettes/occJSDM-lesson-prediction.Rmd`, `.md`, `vignettes/occJSDM-lesson-prediction_files/figure-gfm/` (two PNGs moved with `git mv`)
- Modify: `vignettes/occJSDM-lesson-3.Rmd`, `.md`; `vignettes/occJSDM-lesson-2.Rmd` line 46 (one link target; re-render its `.md`)
- Modify (registration): `vignettes/lesson-links.R`, `dev/simstudy/vignette-lesson/test_lesson_links.R`, `_config.yml`, `assets/js/lessons.js`, `dev/simstudy/lesson-site/test_lessons.js`, `.Rbuildignore`
- Modify (path only): `dev/simstudy/vignette-lesson/prediction-verify.R` line 263
- Read: `lesson-3-map.md` sections 2d and 4

**Interfaces:**
- Consumes: Task 1's Lesson 3.
- Produces: `vignettes/occJSDM-lesson-prediction.Rmd` with title "Predict occupancy at new sites and compare models", headings `## What this lesson answers`, `## Predict occupancy at genuinely new sites` (anchor `#predict-occupancy-at-genuinely-new-sites`, linked from the spatial lesson) and its five `###` subsections unchanged, `## Where to go next`, `## Appendix: evidence and reproduction`; chunk labels `prediction-load`, `prediction-one-site-meaning`, `prediction-native-call`, `prediction-native-intervals`, `prediction-marginal-recovery`, `prediction-probability-errors`, `prediction-observed-scores`, `prediction-paired-scores`, `prediction-fit-one-factor`, `prediction-diagnostics`, `prediction-extract-waic`, `prediction-waic-values`, `prediction-reproduction`, all bodies unchanged; site entry label "Prediction", `published: false`.

- [ ] **Step 1: Move the two prediction figures**

```bash
mkdir -p vignettes/occJSDM-lesson-prediction_files/figure-gfm
git mv vignettes/occJSDM-lesson-3_files/figure-gfm/prediction-native-intervals-1.png vignettes/occJSDM-lesson-prediction_files/figure-gfm/
git mv vignettes/occJSDM-lesson-3_files/figure-gfm/prediction-marginal-recovery-1.png vignettes/occJSDM-lesson-prediction_files/figure-gfm/
```

- [ ] **Step 2: Create the new file's opening**

Create `vignettes/occJSDM-lesson-prediction.Rmd` with this content (the `theme_set` comment and line are copied from Lesson 3's `load-results` chunk):

````markdown
---
title: "Predict occupancy at new sites and compare models"
output:
  rmarkdown::html_vignette:
    toc: true
    css: teaching.css
    fig_width: 9
    fig_height: 5
vignette: >
  %\VignetteIndexEntry{Predict occupancy at new sites and compare models}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE, collapse = FALSE, comment = "#>",
  fig.width = 9, fig.height = 5, dpi = 150, dev = "png"
)

source("lesson-links.R")
```

## What this lesson answers

Genuine prediction at an unsurveyed site requires keeping its observations out of fitting and averaging appropriately over its unknown conditions. Reusing the fitting sites' covariates is not an independent prediction test. This lesson supplies 300 independently generated sites that neither model has seen, predicts occupancy there with the package's new-site function, and compares two models by how well they predict what actually occurred. Spatial prediction is in [the spatial lesson](occJSDM-lesson-2.md).

It continues [Lesson 3](occJSDM-lesson-3.md), which reads the same fits' outputs against truth, and uses the same non-spatial community: 100 sites, 10 species, two measured environmental covariates, two measured traits, three field samples per site, two primers and six PCR replicates per primer.

All teaching code is visible. Figures and tables use compact saved summaries; chunks labelled `eval=FALSE` show how to obtain the underlying outputs from a full fit, where `fitmodel` is the full object returned by `runOccJSDM()`. The [appendix](#appendix-evidence-and-reproduction) shows how to rebuild the new-site comparison.

```{r load-results, message=FALSE, warning=FALSE}
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

lesson <- readRDS("teaching-data/nonspatial-lesson.rds")
known_truth <- lesson$input$sim$true_params
species_order <- colnames(lesson$input$sim$data_list$OTU)

# This theme also works after occJSDM loads its ternary-plot dependency.
theme_set(ggtern::theme_bw(base_size = 12))
```

In the figures, **black crosses or lines show truth**. Blue points or lines show estimates, and blue intervals show posterior uncertainty. An interval is not a measurement of how far the estimate actually is from truth; simulation lets us check both separately.

````

- [ ] **Step 3: Move the prediction block**

Cut from Lesson 3 everything from the line `## Predict occupancy at genuinely new sites` up to, but not including, `## Put the observations, inferred states and truth in one table` (main lines 786 to 1085), and paste it at the end of the new file. Remove one blank line from the double blank line at its end. Do not change any chunk body or header. Then make only these edits inside the moved text:
  - The sentence at the end of `### Which true probability should a new-site prediction recover?` that refers to "the earlier response profiles" (main line 862): turn "the earlier response profiles" into "[Lesson 3's response profiles](occJSDM-lesson-3.md#what-does-an-effect-mean-for-a-species-distribution)".
  - Directly before the `prediction-diagnostics` chunk, add the sentence: "The screen below is the one [Lesson 3's diagnostics](occJSDM-lesson-3.md#check-computation-as-well-as-ecological-recovery) introduces: flag a parameter whose Rhat exceeds 1.01 or whose effective sample size is below 400."
  - Any other `](#...)` link in the moved text whose heading stayed in Lesson 3: point it at `occJSDM-lesson-3.md#...`. Find them with `/usr/bin/grep -n "](#" vignettes/occJSDM-lesson-prediction.Rmd`.

- [ ] **Step 4: Close the new lesson and give it its appendix**

Append to the new file:

````markdown
## Where to go next

[Lesson 3](occJSDM-lesson-3.md) reads the fitted outputs at the surveyed sites. [The four-JSDM comparison](occJSDM-lesson-4.md) asks the same two prediction questions of four packages on a perfectly observed community. The appendix below holds the commands that rebuild this lesson's comparison.

## Appendix: evidence and reproduction

This appendix is for readers who want to reproduce the lesson; its conclusions do not depend on reading it.

````

Then cut from Lesson 3 the paragraph beginning "The new-site comparison needs the additional one-factor fit once." and the `prediction-reproduction` chunk after it, paste them at the end of the new file, and in that paragraph replace "Follow the README" with "Follow [the lesson build README](https://github.com/AlexDiana/occJSDM/blob/main/dev/simstudy/vignette-lesson/README.md)".

- [ ] **Step 5: Bridge what Lesson 3 says about the moved material**

In `vignettes/occJSDM-lesson-3.Rmd`:
  - In `## What this lesson answers` (main line 33), replace the sentence "The new-site comparison below adds one fit to the same PCR observations, changing only the number of hidden site factors from two to one." with "[The prediction lesson](occJSDM-lesson-prediction.md) adds one fit to the same PCR observations, changing only the number of hidden site factors from two to one."
  - Replace the last paragraph of `## Distinguish fitted probabilities, occupancy states and new-site predictions` (main line 784, beginning "Genuine prediction at an unsurveyed site") with: "Prediction at unsurveyed sites, and comparing models by what they predict there, are the subject of [the prediction lesson](occJSDM-lesson-prediction.md)."
  - In the WAIC paragraph at the end of `### Summaries and what to do next` (main line 1601), replace "The new-site section above demonstrates extraction" with "[The prediction lesson](occJSDM-lesson-prediction.md#check-the-additional-fit-and-understand-the-waic-limitation) demonstrates extraction". This is the paragraph PR #14 rewrites.
  - In the paragraph after the `load-full-fit` chunk, delete the sentence "The separate prediction builder is the command that fits the one-factor comparison model."
  - In the function-finder table row "How well does it predict unsurveyed sites?", replace its third cell with "[The prediction lesson](occJSDM-lesson-prediction.md): 300 independent non-spatial sites, with clearly distinguished probability targets". (Task 3 converts the table to bullets.)
  - In `vignettes/occJSDM-lesson-2.Rmd` (spatial, main line 46), change the link target `occJSDM-lesson-3.md#predict-occupancy-at-genuinely-new-sites` to `occJSDM-lesson-prediction.md#predict-occupancy-at-genuinely-new-sites`; leave its link text "Lesson 3" for Task 7.

- [ ] **Step 6: Check that nothing was lost in the move**

```sh
diff <(git show HEAD:vignettes/occJSDM-lesson-3.Rmd | sort) <(cat vignettes/occJSDM-lesson-3.Rmd vignettes/occJSDM-lesson-prediction.Rmd | sort) | /usr/bin/grep "^<" | head -40
```

Expected: only the lines replaced in Steps 3 to 5 (the old sentence at line 33, the old line 784 paragraph, the old WAIC paragraph, the old load-full-fit paragraph, the old finder row, the old line 862 sentence, the old "Follow the README" paragraph, and one empty line). Any other line means text was dropped.

- [ ] **Step 7: Register the provisional lesson, unpublished**

  - `vignettes/lesson-links.R` and `dev/simstudy/vignette-lesson/test_lesson_links.R`: add `"occJSDM-lesson-prediction"` to the `lessons` vector after `"occJSDM-lesson-3"`.
  - `_config.yml`: add `  - vignettes/occJSDM-lesson-prediction.md` after the `occJSDM-lesson-4.md` line and `  - vignettes/occJSDM-lesson-prediction_files/` after the `occJSDM-lesson-3_files/` line.
  - `.Rbuildignore`: after `^vignettes/occJSDM-lesson-intuition_files$` add `^vignettes/occJSDM-lesson-prediction\.(Rmd|md)$` and `^vignettes/occJSDM-lesson-prediction_files$`.
  - `assets/js/lessons.js`: insert after the `occJSDM-lesson-3` entry `    { stem: "occJSDM-lesson-prediction", label: "Prediction", title: "Predict occupancy at new sites and compare models", published: false },`.
  - `dev/simstudy/lesson-site/test_lessons.js`: change the `occJSDM-lesson-4.html` index assertion from `6` to `7`; change the stem filter to `/^occJSDM(-lesson-(\d+|intuition|prediction))?\.Rmd$/`; change the label array to `["Quickstart", "Intuition", "Lesson 0", "Lesson 1", "Lesson 2", "Lesson 3", "Prediction", "Lesson 4"]`; change the last neighbours assertion to `assert.deepEqual(nav.neighbours(7, nav.LESSONS), { previous: nav.LESSONS[6], next: null });`. The published-label guard and the position assertions are unchanged.
  - `dev/simstudy/vignette-lesson/prediction-verify.R` line 263: `"vignettes/occJSDM-lesson-3.Rmd"` becomes `"vignettes/occJSDM-lesson-prediction.Rmd"`. Nothing else in the file changes.

- [ ] **Step 8: Render, verify and commit**

Render `occJSDM-lesson-prediction`, `occJSDM-lesson-3` and `occJSDM-lesson-2` with both commands. Apply the PNG policy: `zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/revert-unchanged-png.zsh vignettes/occJSDM-lesson-prediction_files`, then `vignettes/occJSDM-lesson-3_files`, then `'vignettes/teaching-data/lesson-2-*'` (a quoted glob, so git matches the file names; expected: all reverted). Then:

```sh
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/run-step.zsh 53-prediction-verify Rscript dev/simstudy/vignette-lesson/prediction-verify.R /Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples /Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples/prediction
Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/check-lesson-links.R vignettes/occJSDM-lesson-prediction.Rmd vignettes/occJSDM-lesson-3.Rmd vignettes/occJSDM-lesson-2.Rmd
```

Expected: "All prediction verification checks passed."; the link checker's success line. Run the Global Constraints checks (link test now "All 8 lesson sources ..."; site test `fail 0`). Then:

```bash
git add -A vignettes/occJSDM-lesson-prediction.Rmd vignettes/occJSDM-lesson-prediction.md vignettes/occJSDM-lesson-prediction_files vignettes/occJSDM-lesson-3.Rmd vignettes/occJSDM-lesson-3.md vignettes/occJSDM-lesson-3_files vignettes/occJSDM-lesson-2.Rmd vignettes/occJSDM-lesson-2.md vignettes/lesson-links.R dev/simstudy/vignette-lesson/test_lesson_links.R _config.yml .Rbuildignore assets/js/lessons.js dev/simstudy/lesson-site/test_lessons.js dev/simstudy/vignette-lesson/prediction-verify.R
git status --short
git commit -m "Move new-site prediction and model comparison out of Lesson 3 into a provisional lesson

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Restructure Lesson 3 in place

**Files:**
- Modify: `vignettes/occJSDM-lesson-3.Rmd`, `.md`, `occJSDM-lesson-3_files/` (seven PNGs of removed custom figures deleted with `git rm`)
- Modify: `TODO.md` (one new item), `dev/simstudy/vignette-lesson/remaining-plots-README.md` (line 3)
- Read: `lesson-3-map.md` sections 2a to 2l and 3, the spec's "Lesson 3" bullets

**Interfaces:**
- Consumes: Task 1's chunks and Task 2's Lesson 3.
- Produces: the final Lesson 3 outline below. Headings kept verbatim because other files link to them: `## Check computation as well as ecological recovery`, `### When chains settle on two different explanations`, `### Put the true coefficient beside its estimate`, `## Residual species associations: did we recover what was put in?`, `## Ordination: compare the combined effect before naming the axes`, `## Variation partitioning: an allocation within the model`, `## References and further reading`. New headings: `## How many effects does each fit resolve?`, `### Inspect a primer's detection rate`, `## Find the function for your question`, `## Where to go next`, `## Appendix: evidence and reproduction` and its subsections `### A real cancellation inside this simulated community`, `### Read native ordination plots after aligning their axes`, `### The mirror-labelling study`, `### Read trace draws from the fit object`, `### The corrected response-curve helper`, `### Reproduce the teaching figures` (anchor `#reproduce-the-teaching-figures`).

Target outline, in reading order (chunk labels in brackets; "gallery" means the old `## Appendix: use the package's plotting functions`, which disappears):

1. YAML without `editor_options`; `setup`.
2. `## What this lesson answers`: opening prose, [load-results], the truth convention sentence.
3. `## Check computation as well as ecological recovery`: intro; `### Get the diagnostics for your own fit` [extract-full-fit-diagnostics], parameter-block bullets; `### Find the parameters that need attention` [flag-default-diagnostics, flag-alternative-diagnostics]; `### Read a traceplot for a collection covariate` [collection-trace]; `### Inspect a primer's detection rate` [primer-trace]; `### Follow up an actual diagnostic warning` [flagged-field-trace]; `### When chains settle on two different explanations` (short warning, [extract-chain-means], warning sign and action); `### Summaries and what to do next` [fitting-diagnostics, coefficient-diagnostics].
4. `## How many effects does each fit resolve?` [summarize-effect-evidence]; `### Put the true coefficient beside its estimate` [extract-environment-coefficients, native-read-examples, native-plot-setup, native-environment-example + image, native-environment-2-example + image, native-environment-perfect-1-example + image, native-environment-perfect-2-example + image].
5. `## What does an effect mean for a species' distribution?` [extract-response-profile, remaining-plots-read, remaining-plots-setup, remaining-plots-gradient-1 + image, remaining-plots-gradient-2 + image, extract-baseline-probabilities, native-baseline-occupancy-example + image, baseline-probabilities].
6. `## Traits ask a harder, different question` [extract-trait-coefficients, native-traits-setup, native-traits-gradient-1 + image, native-traits-gradient-2 + image, native-traits-perfect-gradient-1 + image, native-traits-perfect-gradient-2 + image, summarise-trait-array].
7. `## Residual species associations: did we recover what was put in?` [extract-correlations, native-correlation-example + image].
8. `## Ordination: compare the combined effect before naming the axes` [extract-ordination, ordination-contribution (custom), ordination-standard, ordination-ordinary-biplot-image].
9. `## Variation partitioning: an allocation within the model` [extract-variation, remaining-plots-variation + image, variation-truth (custom)].
10. `## Collection effects and detection probabilities` [native-collection-example + image, extract-collection-coefficients, native-baseline-collection-example + image]; `### Laboratory true-positive and false-positive rates by primer` [native-primer-1-example + image, native-primer-2-example + image]; `### Separate false-positive and detection-rate plots` [remaining-plots-stage1-fp, -stage2-fp, -detection + images, combine-rate-plots]; `### Would more field samples or PCR replicates help detection?` [expected-detections-formula, native-effort-k-example + image, native-effort-m-example + image].
11. `## Distinguish fitted probabilities, occupancy states and new-site predictions` [extract-conditional-occupancy].
12. `## Put the observations, inferred states and truth in one table` and its subsection, unchanged.
13. `## Find the function for your question` (bullets).
14. `## Where to go next`.
15. `## References and further reading` (Ji et al. only).
16. `## Appendix: evidence and reproduction` with the six subsections listed under Interfaces.

- [ ] **Step 1: Record the before state**

Keep the rendered markdown committed at the end of Task 2 for comparison: `cp vignettes/occJSDM-lesson-3.md /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/lesson-3-before.md`. Note `wc -l vignettes/occJSDM-lesson-3.Rmd`.

- [ ] **Step 2: Move commit, part 1: diagnostics to the front**

Cut `## Check computation as well as ecological recovery` through the end of `### Summaries and what to do next` (up to, not including, `## Appendix: use the package's plotting functions`) and paste it directly after the truth-convention paragraph ("In the figures, **black crosses or lines show truth**...") at the end of `## What this lesson answers`. Then, still within the moved block, cut these to a holding area at the end of the file under a temporary line `## Appendix: evidence and reproduction` (the heading is final; Step 7 orders its contents):
  - From `### Read a traceplot for a collection covariate`: the second and third sentences of its first paragraph ("First locate the collection covariate by name. The array has dimensions ...") and the `extract-collection-trace` chunk, plus the sentence "That call plots every species." from the paragraph after it.
  - From `### Extract a primer and inspect its detection rate`: its first paragraph ("The same method applies to `p` and `q` ...") and the `extract-primer-trace` chunk.
  - From `### Follow up an actual diagnostic warning`: the first two sentences of its first paragraph ("For a species-only parameter such as `theta0` ... Do not apply the four-index expression above to it.") and the `extract-field-trace` chunk.
  - From `### When chains settle on two different explanations`: the paragraphs from "One of our simulation studies" (the second sentence of the first paragraph onward) through the paragraph after the `mirror-labelling-chains` chunk ("Two chains of the original fit sit far above the rest ..."), the `mirror-labelling-chains` chunk itself, the array-layout sentences of the paragraph before `extract-chain-means` (from "`theta0_output` and `jsdm_output$B0_output` are species by iteration by chain." to the end of that paragraph), and, from the last paragraph, the clause "and in the study the counts differed between the two ways of running chains for reasons not established" and every sentence from "The study ran 8 chains" to the end, including the study-report link.
  - The whole of `### Where to find other parameter draws` with its table and the `direct-coefficient-diagnostics` chunk.

- [ ] **Step 3: Move commit, part 2: fold the gallery into the sections**

Move these gallery blocks, chunks and their prose unchanged, to the positions in the target outline:
  - The gallery's first paragraph ("These are the package's own plots ...") and the chunks `native-read-examples` and `native-plot-setup`: into `### Put the true coefficient beside its estimate`, after `extract-environment-coefficients` and the paragraph after it.
  - `native-environment-example` and its image, the Task 1 gradient 2 and perfect-observation chunks with their images and lead sentences: directly after `native-plot-setup`. The custom chunk `environmental-coefficients` goes to a holding list for deletion in Step 5. The gallery's guidance paragraph "A cross outside its bar ..." follows these figures.
  - `### Occupancy response curves with the package helper`: its two opening paragraphs, `remaining-plots-read`, `remaining-plots-setup`, both gradient chunks and images and the paragraph between them, into `## What does an effect mean for a species' distribution?` after `extract-response-profile`. Its two closing paragraphs (the PR #13 history, beginning "The original walkthrough used `plotCovariateEffect()`" in Lesson 3's corrected form) go to the appendix holding area. The custom `response-profile-1` and `response-profile-2` go to the deletion list.
  - `### Baseline occupancy and collection probabilities`: its first paragraph's occupancy sentences and `native-baseline-occupancy-example` with its image into the response section after `extract-baseline-probabilities`; its collection sentences ("Baseline collection is conditional on presence ...") and `native-baseline-collection-example` with its image into `## Collection effects and detection probabilities` after the paragraph that begins "Two related helpers answer different questions about collection."
  - `### Use the native trait-coefficient plots`: its chunks, images and paragraphs stay in the traits section, now including the Task 1 perfect-observation trait chunks; its heading is removed so the plots are the section's main figures. The custom `trait-effects` goes to the deletion list. `### A real cancellation inside this simulated community` moves whole to the appendix holding area.
  - `### Residual correlations and their uncertainty`: its paragraph, `native-correlation-example` and image into the associations section after `extract-correlations`. The custom `correlation-truth` goes to the deletion list.
  - Ordination: `ordination-standard`, its lead-in sentence ("For an ordinary analysis, the same five native functions ...") and the Task 1 ordinary-biplot paragraph and image move into the main section after `ordination-contribution` and its reading paragraph. Everything else under `### Read native ordination plots after aligning their axes` (its paragraphs, `ordination-load-examples`, `ordination-align`, `ordination-score-table`, and the three `####` subsections with their chunks and images) moves to the appendix holding area under that same `###` heading. Then copy, so it appears in both places, the paragraph "A loading is a species' response to a unit change in a hidden site score ..." and the first two sentences of the biplot reading paragraph ("A site lying farther in a species-arrow direction ... The complete probability also includes the intercept and measured environmental effects.") into the main section after the ordinary biplot.
  - Task 1's `### Variation partitioning with truth`: its paragraph, `remaining-plots-variation` and image into the variation section after `extract-variation`, its heading removed.
  - `native-collection-example` and its image into the collection section directly after the section's first paragraph; the gallery's environmental and collection guidance paragraph is already placed with the environmental figures. The custom `collection-effects` goes to the deletion list.
  - `### Laboratory true-positive and false-positive rates by primer` and `### Separate false-positive and detection-rate plots`, each whole, into the collection section, before `### Would more field samples or PCR replicates help detection?`.
  - `### Cumulative detections: compare survey outcomes with survey outcomes`: its paragraphs, both effort chunks and images, into the detection-effort subsection directly after the paragraph that begins "The exporter applies that formula" (so before the custom chunk), its heading removed. The custom `detection-effort` goes to the deletion list.
  - Delete the now-empty `## Appendix: use the package's plotting functions` heading.
  - Split `## Reproduce the extraction or find a function`: its reproduction paragraphs, the unnamed bash chunk and `load-full-fit` go to the appendix holding area; its pipe table stays where it is for Step 6.

- [ ] **Step 4: Check the moves and commit them**

```sh
diff <(git show HEAD:vignettes/occJSDM-lesson-3.Rmd | sort) <(sort vignettes/occJSDM-lesson-3.Rmd)
```

Expected: only the removed gallery heading, the removed `### Use the native trait-coefficient plots`, `### Occupancy response curves with the package helper`, `### Residual correlations and their uncertainty`, `### Cumulative detections: ...`, `### Variation partitioning with truth` and `## Reproduce the extraction or find a function` headings, the added `## Appendix: evidence and reproduction` heading, the copied ordination sentences, and lines split at sentence boundaries in Steps 2 and 3 (each split paragraph appears as two lines whose text together equals the old one). Then:

```bash
git add vignettes/occJSDM-lesson-3.Rmd
git commit -m "Lesson 3: move diagnostics to the front, fold the plot gallery into the sections, gather the appendix

Moves only; the next commit edits.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 5: Remove the duplicated custom figures and bridge the prose**

  - Delete the seven custom chunks on the deletion list: `environmental-coefficients`, `response-profile-1`, `response-profile-2`, `trait-effects`, `correlation-truth`, `collection-effects`, `detection-effort`.
  - Delete their figure files: `git rm vignettes/occJSDM-lesson-3_files/figure-gfm/{environmental-coefficients,response-profile-1,response-profile-2,trait-effects,correlation-truth,collection-effects,detection-effort}-1.png`.
  - Change only the words of the surrounding prose that described a deleted figure's marks, so they describe the package plot now in its place: in the environmental section replace "For the complete comparison we use the same calculation for every coefficient." with "The package's plots below put each coefficient beside its generating value.", and in "Read each species in three steps" replace "how far is the blue point from it, and how wide is the blue interval?" with "is it inside the bar, and how wide is the bar?" (the package plots draw no posterior mean); after the correlation figure replace "We show those true cells in grey and exclude them from" with "The package plot labels those cells `NA` and we exclude them from" (keep the rest of that sentence); in the detection-effort subsection replace the paragraph beginning "The exporter applies that formula to the known probabilities" with "The package's `plotCumulativeSpeciesDetections()` below simulates whole surveys rather than computing this expectation, so its intervals include the random outcome of a single survey; the orange truth is the exact distribution of that outcome under the generating values.", and in the paragraph after the effort figures replace "the one-sample curve levels off" with "the one-sample panel levels off". In the moved native guidance, replace "not an interval around the analytic expectation taught earlier" with "not an interval around the expected count from the formula above".
  - Add the why-sentence for each remaining custom figure. Before `ordination-contribution`: "This figure is custom because no package plot draws the combined contribution, the one comparison with truth that does not depend on how the hidden axes are rotated." Before `variation-truth`: "The scatter below is custom because the package's triangle does not label species and, in a non-spatial model, puts every point on one edge; plotting each species' true share against its estimate, for both fits, is the comparison the triangle cannot draw."
  - In `load-results`, delete the `environment_labels` definition only if `/usr/bin/grep -c environment_labels vignettes/occJSDM-lesson-3.Rmd` shows no use outside it.

- [ ] **Step 6: Make the diagnostics and the appendix read in their new order**

  - `### Read a traceplot for a collection covariate`: after its first sentence add: "`returnConvergenceDiagnostics()` above is the public route to the numbers; `plotTraceplot()` draws the traces from an array of draws that keeps chains separate. The package has no public function yet that returns such an array from a fit: the [appendix](#read-trace-draws-from-the-fit-object) shows how to take one from `fitmodel$results_output`, and a per-chain accessor is planned (`TODO.md`)." In the next paragraph, change "using the same scale conversion as the coefficient figure above" to "using the same scale conversion as the collection-coefficient plot below".
  - Rename `### Extract a primer and inspect its detection rate` to `### Inspect a primer's detection rate`.
  - `### When chains settle on two different explanations`: after its kept first sentence ("In that example the chains explore the same range of values, only slowly.") add: "A different failure is possible: chains that each settle on a different explanation of the same observations, each looking stable on its own, so that the pooled summary averages the two. Rhat compares chains, so it flags this only when chains actually land in different explanations. One of our simulation studies found such a mirror explanation when laboratory contamination was set far above what the default priors assume; the [appendix](#the-mirror-labelling-study) describes it." Keep "To check your own fit, summarise each chain separately for every species." and the sentence after it, `extract-chain-means`, the warning-sign paragraph and the shortened action paragraph.
  - `### Summaries and what to do next` keeps the Task 2 WAIC paragraph.
  - In `## Collection effects and detection probabilities`, the diagnostics section now comes first: change "use the separate chain arrays in the [diagnostics section](#check-computation-as-well-as-ecological-recovery) for convergence checks" to "use the separate chain arrays, as the [diagnostics section](#check-computation-as-well-as-ecological-recovery) at the start of this lesson does, for convergence checks".
  - In `## What this lesson answers`, change the link `[reproduction section](#reproduce-the-extraction-or-find-a-function)` to `[reproduction instructions](#reproduce-the-teaching-figures)`, and change "You do not need the spatial Lesson 2 first." to "You do not need [the spatial lesson](occJSDM-lesson-2.md) first."
  - Delete the paragraph beginning "The older `sampleresults` object currently shipped with the package" and rename `## Do environmental effects really go undetected?` to `## How many effects does each fit resolve?`.
  - Delete the three lines 13 to 15 of the YAML (`editor_options:`, `  markdown:`, `    wrap: none`).
  - In the traits section, after the paragraph that begins "Three of the four generating trait effects are nonzero.", add: "An oracle decomposition in the [appendix](#a-real-cancellation-inside-this-simulated-community) shows why one of them stays unresolved."

- [ ] **Step 7: Order and introduce the appendix**

Arrange the holding area under `## Appendix: evidence and reproduction` in this order, each block unchanged except as stated:
  - The opening sentence: "This appendix is for readers who want the evidence behind the lesson or need to reproduce it; the lesson's conclusions do not depend on reading it."
  - `### A real cancellation inside this simulated community` (whole).
  - `### Read native ordination plots after aligning their axes` (whole, with its `####` subsections).
  - `### The mirror-labelling study`: the moved study paragraphs, `mirror-labelling-chains`, the reading paragraph, the array-layout sentences, the moved clause (as the sentence "In the study, the counts differed between the two ways of running chains for reasons not established.") and the moved Beta(1, 100) sentences with the report link.
  - `### Read trace draws from the fit object`: first the sentence "These are the extraction steps behind the three traceplots in the diagnostics section. They read internal slots of the fit object, which may change between package versions; a per-chain accessor is planned (`TODO.md`)."; then `#### A collection covariate` (its two moved sentences, `extract-collection-trace`, "That call plots every species."), `#### A primer` (its moved paragraph and `extract-primer-trace`), `#### The field-contamination parameter` (its two moved sentences and `extract-field-trace`), `#### Where to find other parameter draws` (the old subsection's content; its table becomes bullets, one per row: "`jsdm_output$B0_output`: occupancy intercept; species, iteration, chain.", and so on for the six rows).
  - `### The corrected response-curve helper`: the two moved PR #13 paragraphs.
  - `### Reproduce the teaching figures`: the moved reproduction paragraphs, the bash chunk and `load-full-fit`. Replace "follow `dev/simstudy/vignette-lesson/README.md`" with "follow [the lesson build README](https://github.com/AlexDiana/occJSDM/blob/main/dev/simstudy/vignette-lesson/README.md)".

- [ ] **Step 8: Tables, finder, closing and references**

  - The `returnConvergenceDiagnostics()` parameter-block table becomes six bullets of the form "- `beta0_psi`: baseline occupancy on the log-odds scale; `label1` species, `label2` placeholder."
  - Replace the function-finder pipe table with a new section `## Find the function for your question`, placed after the latent-table section, holding one bullet per row in the form "- **How does each species respond to the environment?** `returnOccupancyCovariates()`, ... Truth check: ..." with these truth-check texts: data and fitting, "Lessons 0 and 1"; environment, "the coefficient and response-curve plots above, for both fits"; baseline, "the baseline table and package plot above"; traits, "the package trait plots for both fits above, with standardized truth; the cancellation diagnostic in the appendix"; associations, "the package heat map with true correlations above"; ordination, "the combined contribution and the ordinary biplot above; truth-aligned scores, loadings and biplot in the appendix"; variation, "the package triangle with true shares and the matching fractions above"; collection, "the collection plots above; observation process in Lesson 1"; PCR failures, "combined and separate rate plots above; actual cases in Lesson 1"; sampling effort, "package survey-outcome intervals with exact true ranges above"; a particular site or sample, "package tables above, with matching states and probabilities"; computation, "diagnostics at the start of this lesson; these have no single simulated true value"; unsurveyed sites, the Task 2 cell text; spatial prediction, "Spatial model outputs. [The spatial lesson](occJSDM-lesson-2.md) works through a site-arrangement sweep; this lesson does not validate spatial outputs." Keep every function list of the table as it was, and drop the phrase "corrected conditional response curves are also available".
  - Add, after the finder: `## Where to go next` with: "[The prediction lesson](occJSDM-lesson-prediction.md) predicts occupancy at 300 new sites and compares two models by what actually occurred there. The appendix below holds the evidence that needs the simulation's truth or the full fits: the trait cancellation, the truth-aligned ordination plots, the mirror-labelling study, reading trace draws from the fit object, the corrected response-curve helper, and the commands that reproduce the teaching figures."
  - In `## References and further reading`, delete the Cai et al. (2025), Leibold et al. (2021) and Pichler et al. (2025) entries; keep the heading and Ji et al. (2025).
  - Check: `git grep -n "^|" -- vignettes/occJSDM-lesson-3.Rmd | wc -l` gives 0.

- [ ] **Step 9: The two documents that describe Lesson 3**

  - `TODO.md`: directly after the struck-through DONE item titled "Guidance on mirror labellings under high contamination", add the item: "- **Per-chain draws accessor (added 3 October 2026, lesson rewrite Phase B):** add an exported function that returns one parameter's posterior draws with chains kept separate, in the array shape `plotTraceplot()` takes, so that the outputs lesson can show traceplots and per-chain checks without reading `fitmodel$results_output`. Until then, Lesson 3's appendix shows the internal slots and its diagnostics section says the accessor is planned."
  - `dev/simstudy/vignette-lesson/remaining-plots-README.md` line 3: change "is mirrored in the Lesson 3 native plotting appendix" to "is mirrored in Lesson 3's sections".

- [ ] **Step 10: Render, check and commit**

Render Lesson 3 with both commands. Read the rendered `.md` once from top to bottom against the outline and `lesson-3-before.md`: every figure present, no "above" or "below" pointing the wrong way in the sentences touched, no dangling reference to a removed figure. Apply the PNG policy to `vignettes/occJSDM-lesson-3_files` (expected: moved knit-time figures reverted as unchanged). Run:

```sh
S=/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/run-step.zsh; A=/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
zsh $S 54-verify_native_plots Rscript dev/simstudy/vignette-lesson/verify_native_plots.R $A
zsh $S 55-native-traits-verify Rscript dev/simstudy/vignette-lesson/native-traits-verify.R $A
zsh $S 56-ordination-verify Rscript dev/simstudy/vignette-lesson/ordination-verify.R $A
zsh $S 57-remaining-plots-verify Rscript dev/simstudy/vignette-lesson/remaining-plots-verify.R $A
Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/check-lesson-links.R vignettes/occJSDM-lesson-3.Rmd vignettes/occJSDM-lesson-prediction.Rmd vignettes/occJSDM-lesson-1.Rmd vignettes/occJSDM-lesson-intuition.Rmd
wc -l vignettes/occJSDM-lesson-3.Rmd
```

Expected: the four verifiers pass with their "Lesson 3 displays the exact exported ..." lines (moves keep bodies identical); the link checker passes; about 1,650 to 1,750 lines. Run the Global Constraints checks. Then:

```bash
git add vignettes/occJSDM-lesson-3.Rmd vignettes/occJSDM-lesson-3.md vignettes/occJSDM-lesson-3_files TODO.md dev/simstudy/vignette-lesson/remaining-plots-README.md
git status --short
git commit -m "Lesson 3: package plots with truth as main figures, appendix, bullets, references

Removes the seven custom figures that duplicated a package plot, keeps the
ordination contribution and the variation scatter with the reason each is
custom, separates the function finder from reproduction, and fixes the
Lesson 2 statements and the opening sampleresults paragraph.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Move the ten-community extension out of Lesson 4 into a provisional lesson

**Files:**
- Create: `vignettes/occJSDM-lesson-extension.Rmd`, `.md`
- Rename (git mv): the twelve `vignettes/teaching-data/lesson-4-extension-*-1.png` and `lesson-4-calibration-*-1.png` files to `lesson-6-...`
- Modify: `vignettes/occJSDM-lesson-4.Rmd`, `.md`; `vignettes/LESSON-PLAN.md` line 78 (one link target)
- Modify (registration): `vignettes/lesson-links.R`, `dev/simstudy/vignette-lesson/test_lesson_links.R`, `_config.yml`, `assets/js/lessons.js`, `dev/simstudy/lesson-site/test_lessons.js`, `.Rbuildignore`
- Read: `lesson-4-map.md` sections 2 ("Sections 12 to 15", "Full-results appendix", "The echo leak") and 3

**Interfaces:**
- Consumes: Lesson 4 at main.
- Produces: `vignettes/occJSDM-lesson-extension.Rmd`, title "Repeat the four-JSDM comparison across ten communities, with traits", `fig.path` prefix `teaching-data/lesson-6-` (markdown) and `occJSDM-lesson-extension_files/figure-html/` (HTML); headings `## What this lesson adds`, `## 1. Does more data help when ecology becomes harder?` (with `### Which fits passed, remained flagged or failed?`, anchor linked from `LESSON-PLAN.md`), `## 2. Do traits explain species responses, and do they improve prediction?`, `## 3. Across communities, are estimates biased?`, `## 4. Do 95% intervals contain the truth?`, `## What this lesson establishes`, `## Appendix: evidence and reproduction`; every `extension-*`, `calibration-report-*` and kept `calibration-*` chunk label unchanged; site entry label "Extension", `published: false`.

- [ ] **Step 1: Rename the extension's figures**

```sh
for f in vignettes/teaching-data/lesson-4-extension-*-1.png vignettes/teaching-data/lesson-4-calibration-*-1.png; do git mv $f ${f/lesson-4-/lesson-6-}; done
git status --short | /usr/bin/grep -c "^R"
```

Expected: 12 renames (six `extension-*`, three `calibration-report-*`, three other `calibration-*`). The data files `lesson-4-extension.rds`, `lesson-4-calibration.rds` and `lesson-4-species-errors.csv` do not match the `-1.png` pattern and stay.

- [ ] **Step 2: Create the new file's opening**

Create `vignettes/occJSDM-lesson-extension.Rmd` with this content (the last two lines of the load chunk are Lesson 4's `load-comparison` lines 87 and 88):

````markdown
---
title: "Repeat the four-JSDM comparison across ten communities, with traits"
output:
  rmarkdown::html_vignette:
    toc: true
    css: teaching.css
    fig_width: 9
    fig_height: 6
vignette: >
  %\VignetteIndexEntry{Repeat the four-JSDM comparison across ten communities, with traits}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(
  echo = TRUE, collapse = FALSE, comment = "#>",
  fig.width = 9, fig.height = 6, dpi = 150, dev = "png",
  # R's vignette cleanup removes generated figure directories after HTML builds.
  # Keep temporary figures separate from the saved teaching inputs.
  fig.path = if (knitr::is_html_output(excludes = c("markdown", "gfm"))) {
    "occJSDM-lesson-extension_files/figure-html/"
  } else {
    "teaching-data/lesson-6-"
  }
)

source("lesson-links.R")
```

## What this lesson adds

[The four-JSDM comparison](occJSDM-lesson-4.md) gave occJSDM, gllvm, sjSDM and Hmsc the same perfectly observed community and compared their probabilities with truth. This lesson repeats that experiment across ten independent communities, compares ecological and trait scenarios, and examines bias and interval coverage. These experiments do not establish a general ranking of the packages. All numbers below come from saved fits.

You will learn to distinguish bias from error size, and to assess interval coverage alongside interval width.

Sections 1 and 2 show their code. Sections 3 and 4 present the saved results as a report and hide the code that renders them; the appendix shows the code for its tables and figures, and every chunk is in this R Markdown source. Knit this file, or run its chunks with `vignettes` as the working directory. Rendering reads two compact results bundles and does **not** fit any models.

```{r load-extension-setup, message=FALSE, warning=FALSE}
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

package_order <- c("occJSDM", "gllvm", "sjSDM", "Hmsc")

# Keep the registered ternary theme elements valid when vignettes share a session.
theme_set(ggtern::theme_bw(base_size = 12))
```

````

- [ ] **Step 3: Move sections 12 to 15 and the full results**

Cut from Lesson 4 everything from `## 12. Does more data help when ecology becomes harder?` to the end of the file (main lines 857 to 1532) and paste it at the end of the new file. Then, in the new file:
  - Renumber the four headings: `## 12.` to `## 1.`, `## 13.` to `## 2.`, `## 14.` to `## 3.`, `## 15.` to `## 4.`; in the text, "section 12" becomes "section 1" and "section 14" becomes "section 3" (find with `/usr/bin/grep -n "[Ss]ection 1[2-5]" vignettes/occJSDM-lesson-extension.Rmd`).
  - In the first paragraph of section 1, replace "The original example above used one small community." with "[The four-JSDM comparison](occJSDM-lesson-4.md) used one small community."
  - Convert the scenario pipe table (header, separator and four rows) to four bullets of the form "- **Baseline:** neither of the complications below. What we want to learn: does collecting more sites improve estimates?"
  - Fix the echo leak: in the `calibration-report-setup` chunk delete the line `knitr::opts_chunk$set(echo = FALSE)`, and add `echo=FALSE` to the headers of exactly these six chunks: `calibration-report-prediction-table`, `calibration-report-bias`, `calibration-report-coefficient-bias`, `calibration-report-coverage-table`, `calibration-report-probability-coverage`, `calibration-report-coefficient-coverage` (for example `{r calibration-report-bias, echo=FALSE, ...}` keeping their other options). Do not add a set-and-reset pair.

- [ ] **Step 4: Reduce the full results to what sections 3 and 4 do not show**

Replace the heading `## Full results, interval widths and methods` with `## Appendix: evidence and reproduction`, and make its first line "This appendix is for readers who want the full tables, interval widths and methods behind sections 3 and 4; the lesson's conclusions do not depend on reading it." Make `### Detailed prediction errors` and `### Detailed interval comparisons` its subsections (headings kept with their attributes). Then:
  - Delete these, each of which repeats what sections 3 and 4 already show: the paragraph before the old heading that begins "The full results below keep probability bias"; the paragraph "These detailed tables and community-level plots retain ..."; the paragraph beginning "The September fits can answer more questions"; the paragraph beginning "An estimate's **signed error** is" (section 3 defines bias and RMSE); in the paragraph after the `calibration-data` chunk, beginning "Here the target is each species'", every sentence except the one about RMSE being the root of the average community mean squared error; the paragraph beginning "The straight-response fits to curved truth illustrate" (main line 1415, whose numbers are section 4's), except its sentence "More data can narrow intervals around an unsuitable response shape; it does not supply the missing curve.", which stays as a paragraph of its own; the paragraph beginning "The baseline is also not a substitute for checking difficult conditions" (main line 1417, the point of section 4's "occJSDM's baseline is not the whole story"); and, in the paragraph beginning "A flagged gllvm rare-species fit produces" (main line 1484, the point of section 4's figure note), its first sentence, keeping "Hmsc keeps its place on the axis, with its coefficient comparison labelled `Different link`."
  - Keep the `calibration-data` chunk: its communities, flagged and bias MCSE columns are not in section 3's table.
  - From the last paragraph of the file, move the sentence "A larger coverage study should add independent communities, assess gllvm/sjSDM marginal-probability intervals, and include a probit-generating arm under a declared fitting and diagnostic protocol." into the closing section of Step 5; the rest of that paragraph stays.
  - Move up into section 4, directly after its first paragraph: the worked coverage paragraph beginning "**Coverage** asks how often an interval includes the generating truth" and the paragraph beginning "Bayesian credible intervals and approximate frequentist confidence intervals"; and the first sentence of the paragraph before `calibration-grid` ("To keep the interval calculation inspectable, we evaluate five fixed settings ...") after section 4's sentence that names the five fixed environmental settings.
  - Keep every other chunk and paragraph in order, including `calibration-error-comparison`, the "Which saved results support this calculation?" bullets, `calibration-grid`, `calibration-combined-targets` (unchanged), the two combined tables, `calibration-probability-intervals`, `calibration-combined-status`, `calibration-native-coefficient-intervals`, `calibration-combined-sensitivity`, `calibration-trait-effects` and `calibration-inspect`.

- [ ] **Step 5: Close the new lesson**

Insert before `## Appendix: evidence and reproduction`:

```markdown
## What this lesson establishes

Across ten independent communities, more training sites reduce prediction error, average bias is small in the baseline for all four packages, and interval coverage differs more among packages than point accuracy does. Rare species, correlated predictors and an unsuitable response shape change the answer, and the diagnostic warnings and the small number of communities prevent a confident overall ranking.

A larger coverage study should add independent communities, assess gllvm/sjSDM marginal-probability intervals, and include a probit-generating arm under a declared fitting and diagnostic protocol.

The appendix below holds the detailed prediction errors, the interval widths and coverage by community, the diagnostic sensitivity, the trait relationships, and how to filter the saved results yourself.
```

The first paragraph restates only claims already made in section 3's "What the comparison says" bullets and its opening paragraph; if a reviewer finds one unsupported by the rendered output, cut that clause.

- [ ] **Step 6: Leave Lesson 4 pointing at the new lesson**

In `vignettes/occJSDM-lesson-4.Rmd`:
  - Main line 40: replace "Sections 1-11 introduce a **worked pilot using one community**. Sections 12-15 extend it to ten independent communities, compare ecological and trait scenarios, and examine bias and interval coverage." with "This lesson works through **one community**. [The ten-community lesson](occJSDM-lesson-extension.md) repeats the experiment across ten independent communities, compares ecological and trait scenarios, and examines bias and interval coverage."
  - Delete objective 7 ("Distinguish bias from error size, and assess interval coverage alongside interval width.").
  - Main line 52: replace "The code is visible in the worked lesson; sections 14-15 present the saved results as a report, with the rendering code in this R Markdown source." with "All teaching code is visible."
  - Main line 839, section 11's second paragraph: replace "Sections 12-15 add independent simulated communities and sample-size comparisons." with "[The ten-community lesson](occJSDM-lesson-extension.md) adds independent simulated communities, sample-size comparisons and a matched trait example." and delete "Section 13 adds a matched trait example."
  - Remove the trailing blank lines so the file ends after the "Return to ..." line.
  - In `vignettes/LESSON-PLAN.md` line 78, change the link target `occJSDM-lesson-4.md#which-fits-passed-remained-flagged-or-failed` to `occJSDM-lesson-extension.md#which-fits-passed-remained-flagged-or-failed`.

- [ ] **Step 7: Check that nothing was lost**

```sh
diff <(git show HEAD:vignettes/occJSDM-lesson-4.Rmd | sort) <(cat vignettes/occJSDM-lesson-4.Rmd vignettes/occJSDM-lesson-extension.Rmd | sort) | /usr/bin/grep "^<"
```

Expected: only the lines Steps 3, 4 and 6 replace or delete (the four old section headings, the old table lines, the `opts_chunk$set(echo = FALSE)` line, the six old chunk headers, the deleted paragraphs and sentences of Step 4, the old lines 40, 50, 52 and 839, the old "The original example above" sentence, the old full-results heading, and blank lines).

- [ ] **Step 8: Register the provisional lesson, unpublished**

  - `vignettes/lesson-links.R` and `test_lesson_links.R`: append `"occJSDM-lesson-extension"` as the last element of `lessons`.
  - `_config.yml`: add `  - vignettes/occJSDM-lesson-extension.md` after the `occJSDM-lesson-prediction.md` line. No `_files/` line: its figures go to `teaching-data/`, which is already excluded, and the HTML figure folder is removed after each render.
  - `.Rbuildignore`: after the prediction lines add `^vignettes/occJSDM-lesson-extension\.(Rmd|md)$`.
  - `assets/js/lessons.js`: append after the `occJSDM-lesson-4` entry (adding a comma to it) `    { stem: "occJSDM-lesson-extension", label: "Extension", title: "Repeat the four-JSDM comparison across ten communities, with traits", published: false }`.
  - `test_lessons.js`: stem filter `/^occJSDM(-lesson-(\d+|intuition|prediction|extension))?\.Rmd$/`; label array `["Quickstart", "Intuition", "Lesson 0", "Lesson 1", "Lesson 2", "Lesson 3", "Prediction", "Lesson 4", "Extension"]`; last neighbours assertion `assert.deepEqual(nav.neighbours(8, nav.LESSONS), { previous: nav.LESSONS[7], next: null });`.

- [ ] **Step 9: Render, verify and commit**

Render `occJSDM-lesson-extension` and `occJSDM-lesson-4` with both commands. In the extension's `.md`, confirm that the six `calibration-report-*` chunks show no code and that the appendix chunks (including `calibration-inspect`, the "For example" chunk) show theirs. Apply the PNG policy to `'vignettes/teaching-data/lesson-6-*'` and `'vignettes/teaching-data/lesson-4-*'` (quoted globs; expected: all reverted). Run:

```sh
Rscript dev/simstudy/jsdm-package-comparison/verify-teaching.R archive dev/simstudy/jsdm-package-comparison
Rscript dev/simstudy/jsdm-package-comparison/extension/verify.R /Users/douglasyu/src/occJSDM/dev/simstudy/results/lesson-4-extension-20260923 .
Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/check-lesson-links.R vignettes/occJSDM-lesson-4.Rmd vignettes/occJSDM-lesson-extension.Rmd vignettes/LESSON-PLAN.md
```

Expected: "Teaching numerical verification passed for every check the committed archive supports." (its `reproduce-community` and `marginal-extraction-example` chunks are still in Lesson 4); "Verified input identity, held-out states, independent probability truth and error summaries for 626 scored results; 640 of 640 fits attempted."; the link checker's success line. Run the Global Constraints checks. Then:

```bash
git add -A vignettes/occJSDM-lesson-extension.Rmd vignettes/occJSDM-lesson-extension.md vignettes/teaching-data vignettes/occJSDM-lesson-4.Rmd vignettes/occJSDM-lesson-4.md vignettes/LESSON-PLAN.md vignettes/lesson-links.R dev/simstudy/vignette-lesson/test_lesson_links.R _config.yml .Rbuildignore assets/js/lessons.js dev/simstudy/lesson-site/test_lessons.js
git status --short
git commit -m "Move the ten-community extension out of Lesson 4 into a provisional lesson

Fixes the echo leak (echo=FALSE on the six report chunks instead of a
global option), reduces the full results to what the report sections do
not show, and renames the extension's figures to the lesson-6 prefix.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Restructure Lesson 4 in place

**Files:**
- Modify: `vignettes/occJSDM-lesson-4.Rmd`, `.md`
- Read: `lesson-4-map.md` sections "Section 3", "Section 10 and the reproduction record", "The two-optima teaching point", "The three pipe tables"; the intuition lesson's `## What a joint model adds to a factorisation`

**Interfaces:**
- Consumes: Task 4's Lesson 4.
- Produces: section 3 heading kept (`## 3. What a joint model adds to factoring a matrix`) with a two-sentence reminder linking `occJSDM-lesson-intuition.md#what-a-joint-model-adds-to-a-factorisation`; `### Exercise: predict one species given another` (anchor `#exercise-predict-one-species-given-another`, linked from the intuition lesson) as the last subsection of section 4; sections 4 to 9 numbered as before; `## 10. What this lesson establishes, and what remains open` (was 11) as the closing section; `## Appendix: evidence and reproduction` with `### Fit the same configurations yourself` (its `#### What did these runs cost?`), `### The two sjSDM optima: evidence`, `### Reproduction record`. Chunk labels `reproduce-community`, `marginal-extraction-example` and `conditional-prediction` unchanged in name and body.

- [ ] **Step 1: Move commit**

  - Cut the exercise (heading `### Exercise: predict one species given another`, its two paragraphs, the `conditional-prediction` chunk and its two closing paragraphs; main lines 231 to 285) and paste it at the end of section 4, after the paragraph that follows `marginal-extraction-example` and before `## 5. Check the computation before interpreting the ecology`.
  - Add `## Appendix: evidence and reproduction` at the end of the file. Move `## 10. Optional: fit the same configurations yourself` there, renamed `### Fit the same configurations yourself`, with its subsection renamed `#### What did these runs cost?`. Move `## Reproduction record` there as `### Reproduction record`, without its last line ("Return to the [Quickstart and lesson guide] ..."), which moves to the end of the closing section.
  - From `### Two local optima in the sjSDM fit`, move to a new appendix subsection `### The two sjSDM optima: evidence`, placed between the two above: the paragraph beginning "We restored the optimiser's usual weak penalty", the sentence-paragraph "We then ran nine more native sjSDM starts ...", the `sjsdm-basins` chunk, the paragraph beginning "Within each hill the starts agree closely", the line "What distinguishes the two solutions? ...", the `sjsdm-two-solutions` chunk and the `sjsdm-solution-errors` chunk.
  - Check with `diff <(git show HEAD:vignettes/occJSDM-lesson-4.Rmd | sort) <(sort vignettes/occJSDM-lesson-4.Rmd)`: only the renamed headings, the new appendix heading and blank lines differ. Commit:

```bash
git add vignettes/occJSDM-lesson-4.Rmd
git commit -m "Lesson 4: exercise after section 4, fitting calls, sjSDM evidence and reproduction record to an appendix

Moves only; the next commit edits.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 2: Reduce section 3 to its reminder**

Replace everything under `## 3. What a joint model adds to factoring a matrix` (the paragraphs that were main lines 221 to 229) with these two sentences as one paragraph: "A JSDM explains what the measured environment leaves over with hidden site factors and species loadings, so it defines the joint probability of the whole species list at a site; [the intuition lesson](occJSDM-lesson-intuition.md#what-a-joint-model-adds-to-a-factorisation) builds this up from a recommender's viewer-by-film table. Two consequences matter here: recording one species at a site is evidence about the others, and the same model answers two different prediction questions, which the next section separates." Nothing of the deleted text is lost: it is in the intuition lesson (`lesson-4-map.md`, "Section 3").

- [ ] **Step 3: Keep the two-optima point in section 5 and point to its evidence**

Under `### Two local optima in the sjSDM fit`, the text now reads, in order: a new first paragraph "Restoring the optimiser's usual weak penalty and refining the saved fits with an exact optimiser showed why: the penalised fitting problem has two genuine local maxima, two hills, and independent starts climb one or the other. The declared selection rule picked the better hill, reached by four of twelve starts, before any truth was read."; the kept paragraph beginning "In the selected solution, species 2 has a large residual spread", with its last sentence changed from "The next table shows that the choice barely matters" to "The appendix tables show that the choice barely matters"; the kept paragraph beginning "The average errors of the two local maxima are almost identical"; the kept paragraph beginning "We select starts by the training fitting criterion"; and a closing sentence "[The appendix](#the-two-sjsdm-optima-evidence) gives the evidence: the deterministic refinement and curvature check, the twelve-start classification and the tables behind these statements."

- [ ] **Step 4: Close the lesson, tables, YAML and pointers**

  - Rename `## 11. What this lesson establishes, and what remains open` to `## 10. What this lesson establishes, and what remains open`. At its end, after the exercises, add "The appendix below holds each package's fitting call and what the runs cost, the evidence for the two sjSDM optima, and the reproduction record." followed by the moved "Return to ..." line.
  - In section 4's subsection on averaging, change "The full checked extraction scripts are recorded in the reproduction notes below." to "The full checked extraction scripts are named in the [reproduction record](#reproduction-record) in the appendix."
  - Convert the section 2 pipe table (four packages) to four bullets, one per package, "- **occJSDM:** ..." keeping each row's cells as a sentence; convert the section 4 pipe table (two questions) to two bullets the same way.
  - Delete the YAML's `editor_options:`, `  markdown:` and `    wrap: none` lines.
  - Check `git grep -n "^|" -- vignettes/occJSDM-lesson-4.Rmd | wc -l` gives 0 and `/usr/bin/grep -n "section 1[0-5]\|sections 1[0-5]\|[Ss]ection 11" vignettes/occJSDM-lesson-4.Rmd` finds nothing stale.

- [ ] **Step 5: Render, verify and commit**

Render Lesson 4 with both commands; read the rendered `.md` from section 3 to the end. Apply the PNG policy to `'vignettes/teaching-data/lesson-4-*'` (a quoted glob; expected: all seven reverted, since no figure's content changed). Run `Rscript dev/simstudy/jsdm-package-comparison/verify-teaching.R archive dev/simstudy/jsdm-package-comparison` (expected pass line as in Task 4) and the link checker on `vignettes/occJSDM-lesson-4.Rmd vignettes/occJSDM-lesson-intuition.Rmd`. Run the Global Constraints checks; the pipe-table count is now 0. Then:

```bash
git add vignettes/occJSDM-lesson-4.Rmd vignettes/occJSDM-lesson-4.md vignettes/teaching-data
git status --short
git commit -m "Lesson 4: section 3 becomes a reminder of the intuition lesson, two-optima point kept with its evidence in the appendix

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Renumber the files, figure prefixes, titles, registrations and script paths

**Files:**
- Rename (git mv): see Step 1
- Modify: the YAML `title` and `VignetteIndexEntry` of the six lessons whose number changes; `fig.path` in `occJSDM-lesson-5.Rmd`, `occJSDM-lesson-6.Rmd`, `occJSDM-lesson-7.Rmd`
- Modify (registration): `vignettes/lesson-links.R`, `dev/simstudy/vignette-lesson/test_lesson_links.R`, `_config.yml`, `assets/js/lessons.js`, `dev/simstudy/lesson-site/test_lessons.js`, `.Rbuildignore`, `dev/simstudy/lesson-site/build_preview.R` line 62
- Modify (path only): `dev/simstudy/vignette-lesson/unbalanced-verify.R` lines 142, 143, 148, 149; `prediction-verify.R` line 263; `dev/simstudy/jsdm-package-comparison/verify-teaching.R` line 44; `extension/finish-study.py` lines 29 and 30
- Regenerate: every lesson's `.md`
- Read: `renumber-inventory.md` sections 2, 3, 4 and 6 (items 6 and 7)

**Interfaces:**
- Consumes: Tasks 2 to 5.
- Produces: stems `occJSDM-lesson-0` to `occJSDM-lesson-7` with the titles of "Final numbering"; figure folders `occJSDM-lesson-{0,1,2,3,4}_files`; figure prefixes `teaching-data/lesson-5-`, `lesson-6-`, `lesson-7-`; site labels "Lesson 0" to "Lesson 7". Cross-references inside the lessons still use old numbers until Task 7.

- [ ] **Step 1: Rename commit (renames only, in collision-safe order)**

```sh
cd /Users/douglasyu/src/occJSDM-worktrees/lesson-phase-b
git mv vignettes/occJSDM-lesson-2.Rmd vignettes/occJSDM-lesson-7.Rmd
git mv vignettes/occJSDM-lesson-2.md vignettes/occJSDM-lesson-7.md
for f in vignettes/teaching-data/lesson-2-*.png; do git mv $f ${f/lesson-2-/lesson-7-}; done
git mv vignettes/occJSDM-lesson-1.Rmd vignettes/occJSDM-lesson-2.Rmd
git mv vignettes/occJSDM-lesson-1.md vignettes/occJSDM-lesson-2.md
git mv vignettes/occJSDM-lesson-1_files vignettes/occJSDM-lesson-2_files
git mv vignettes/occJSDM-lesson-intuition.Rmd vignettes/occJSDM-lesson-1.Rmd
git mv vignettes/occJSDM-lesson-intuition.md vignettes/occJSDM-lesson-1.md
git mv vignettes/occJSDM-lesson-intuition_files vignettes/occJSDM-lesson-1_files
git mv vignettes/occJSDM-lesson-4.Rmd vignettes/occJSDM-lesson-5.Rmd
git mv vignettes/occJSDM-lesson-4.md vignettes/occJSDM-lesson-5.md
for f in vignettes/teaching-data/lesson-4-*-1.png; do git mv $f ${f/lesson-4-/lesson-5-}; done
git mv vignettes/occJSDM-lesson-prediction.Rmd vignettes/occJSDM-lesson-4.Rmd
git mv vignettes/occJSDM-lesson-prediction.md vignettes/occJSDM-lesson-4.md
git mv vignettes/occJSDM-lesson-prediction_files vignettes/occJSDM-lesson-4_files
git mv vignettes/occJSDM-lesson-extension.Rmd vignettes/occJSDM-lesson-6.Rmd
git mv vignettes/occJSDM-lesson-extension.md vignettes/occJSDM-lesson-6.md
git mv dev/simstudy/vignette-lesson/unbalanced-lesson-1.Rmd dev/simstudy/vignette-lesson/unbalanced-lesson-2.Rmd
/bin/ls vignettes/teaching-data/lesson-*-1.png | sed 's/.*lesson-\([0-9]\)-.*/\1/' | sort | uniq -c
git commit -m "Rename the lessons to the final 0 to 7 numbering

Renames only, so history follows each file; the next commit updates
titles, figure prefixes, registrations and script paths. Old 2 became 7
before old 1 became 2, and old 4 became 5 before the prediction lesson
became 4.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

Expected counts before the commit: 7 PNGs with prefix 5, 12 with prefix 6, 8 with prefix 7, none with 2 or 4. The tests fail at this commit; Step 6 makes them pass.

- [ ] **Step 2: Titles and figure prefixes**

  - Set `title:` and `%\VignetteIndexEntry{...}` to the "Final numbering" titles in `occJSDM-lesson-1.Rmd` (prefix "Lesson 1: " to the intuition title), `-2.Rmd` ("Lesson 1:" to "Lesson 2:"), `-4.Rmd` ("Lesson 4: Predict occupancy at new sites and compare models"), `-5.Rmd` ("Lesson 4:" to "Lesson 5:"), `-6.Rmd` ("Lesson 6: Repeat the four-JSDM comparison across ten communities, with traits"), `-7.Rmd` ("Lesson 2:" to "Lesson 7:").
  - `occJSDM-lesson-7.Rmd` setup: `"occJSDM-lesson-2_files/figure-html/"` to `"occJSDM-lesson-7_files/figure-html/"` and `"teaching-data/lesson-2-"` to `"teaching-data/lesson-7-"`.
  - `occJSDM-lesson-5.Rmd` setup: `"occJSDM-lesson-4_files/figure-html/"` to `"occJSDM-lesson-5_files/figure-html/"` and `"teaching-data/lesson-4-"` to `"teaching-data/lesson-5-"`.
  - `occJSDM-lesson-6.Rmd` setup: `"occJSDM-lesson-extension_files/figure-html/"` to `"occJSDM-lesson-6_files/figure-html/"`.

- [ ] **Step 3: Registrations**

  - `vignettes/lesson-links.R` and `dev/simstudy/vignette-lesson/test_lesson_links.R`, the `lessons` vector:

```r
lessons <- c("occJSDM", "occJSDM-lesson-0", "occJSDM-lesson-1",
             "occJSDM-lesson-2", "occJSDM-lesson-3", "occJSDM-lesson-4",
             "occJSDM-lesson-5", "occJSDM-lesson-6", "occJSDM-lesson-7")
```

  (indent as each file already does). The test's fixture links at its lines 21 to 58 use stems 0, 1, 3 and 4 as test material and the anchor `#check-computation-as-well-as-ecological-recovery`, which still exists; they stay.
  - `_config.yml`: the comment "the five teaching lessons" becomes "the eight teaching lessons"; the lesson lines of `exclude` become `vignettes/occJSDM-lesson-0.md` to `vignettes/occJSDM-lesson-7.md` (eight lines) and `vignettes/occJSDM-lesson-{0,1,2,3,4}_files/` (five lines, in that order), then `vignettes/teaching-data/`; the intuition, prediction and extension lines go.
  - `.Rbuildignore`: `[0-4]` becomes `[0-7]` in both lesson lines; delete the intuition, prediction and extension lines. Line `^vignettes/lesson-4-calibration-tables\.R$` stays (the file keeps its name).
  - `assets/js/lessons.js`: the `LESSONS` array becomes the quickstart entry followed by eight entries, `stem` `occJSDM-lesson-0` to `-7`, `label` "Lesson 0" to "Lesson 7", `title` from "Final numbering", all `published: false`; in the comment above `positionText`, "Lesson 1 of 0-4" becomes "Lesson 1 of 0-7".
  - `dev/simstudy/lesson-site/test_lessons.js`: the first test becomes indexes `occJSDM.html` 0, `occJSDM-lesson-0.html` 1, `occJSDM-lesson-1.html` 2, `occJSDM-lesson-1` 2, `occJSDM-lesson-4.html` 5, `occJSDM-lesson-7.html` 8; add `"/occJSDM/vignettes/occJSDM-lesson-intuition.html"` to the ignored-pages list; in the published-label guard keep the assertion and replace its three comment lines with `// positionText() reads the range from the first and last labels, so a`, `// published page without a "Lesson N" label would break it, as the`, `// unnumbered intuition lesson would have done ("Lesson 1 of Intuition-4").`; the position assertions become `positionText(1, pages)` "Lesson 0 of 0–7", `positionText(2, pages)` "Lesson 1 of 0–7", `positionText(8, pages)` "Lesson 7 of 0–7" (en dash, as in the file); neighbours become `(0)` next `LESSONS[1]`, `(3)` previous `[2]` next `[4]`, `(8)` previous `[7]` next `null`; stem filter `/^occJSDM(-lesson-\d+)?\.Rmd$/`; labels `["Quickstart", "Lesson 0", "Lesson 1", "Lesson 2", "Lesson 3", "Lesson 4", "Lesson 5", "Lesson 6", "Lesson 7"]`.
  - `dev/simstudy/lesson-site/build_preview.R` line 62: `stopifnot(length(lessons) == 6)` becomes `== 9`. Lines 9 and 28 name `occJSDM-lesson-1` as an example page and template; leave them.

- [ ] **Step 4: Script paths (approved path-only edits)**

  - `unbalanced-verify.R`: lines 142 and 148, `"vignettes/occJSDM-lesson-1.Rmd"` to `"vignettes/occJSDM-lesson-2.Rmd"`; lines 143 and 149, `unbalanced-lesson-1.Rmd` to `unbalanced-lesson-2.Rmd`. Lines 132 and 133 (Lesson 0) stay.
  - `prediction-verify.R` line 263: `"vignettes/occJSDM-lesson-prediction.Rmd"` to `"vignettes/occJSDM-lesson-4.Rmd"`.
  - `verify-teaching.R` line 44: `"vignettes/occJSDM-lesson-4.Rmd"` to `"vignettes/occJSDM-lesson-5.Rmd"`.
  - `extension/finish-study.py` lines 29 and 30: both `vignettes/occJSDM-lesson-4.Rmd` to `vignettes/occJSDM-lesson-6.Rmd`.
  - The four Lesson 3 verifiers read `vignettes/occJSDM-lesson-3.Rmd`, which keeps its name; they do not change.
  - Check: `git diff -U0 -- dev/simstudy/vignette-lesson/unbalanced-verify.R dev/simstudy/vignette-lesson/prediction-verify.R dev/simstudy/jsdm-package-comparison/verify-teaching.R dev/simstudy/jsdm-package-comparison/extension/finish-study.py` shows only path strings changing.

- [ ] **Step 5: Re-render every lesson**

Render, with the markdown command, `occJSDM-lesson-0` to `occJSDM-lesson-7`. Apply the PNG policy to `vignettes/occJSDM-lesson-0_files`, `-1_files`, `-2_files`, `-3_files`, `-4_files`, `'vignettes/teaching-data/lesson-5-*'`, `'vignettes/teaching-data/lesson-6-*'` and `'vignettes/teaching-data/lesson-7-*'` (quoted globs; expected: all reverted; a "kept" line means a figure's content changed and must be explained before committing). Confirm `git status --short vignettes/ | /usr/bin/grep "^??"` lists nothing (no stray old-prefix figure was written).

- [ ] **Step 6: Verify and commit**

```sh
S=/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/run-step.zsh; A=/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
zsh $S 58-unbalanced-verify Rscript dev/simstudy/vignette-lesson/unbalanced-verify.R $A/unbalanced
zsh $S 59-prediction-verify Rscript dev/simstudy/vignette-lesson/prediction-verify.R $A $A/prediction
Rscript dev/simstudy/jsdm-package-comparison/verify-teaching.R archive dev/simstudy/jsdm-package-comparison
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
node --test dev/simstudy/lesson-site/test_lessons.js
```

Expected: the unbalanced verifier's "Verified paired deletion, ... and displayed code." lines; "All prediction verification checks passed."; the comparison pass line; "All 9 lesson sources are free of dynamic or encoded link destinations.", "All 9 rendered lessons keep each lesson link on one line for GitHub Pages.", "Shared lesson-link regression checks passed."; `pass 10`, `fail 0`. Then:

```bash
git add -A vignettes dev/simstudy/vignette-lesson/test_lesson_links.R dev/simstudy/vignette-lesson/unbalanced-verify.R dev/simstudy/vignette-lesson/prediction-verify.R dev/simstudy/jsdm-package-comparison/verify-teaching.R dev/simstudy/jsdm-package-comparison/extension/finish-study.py dev/simstudy/lesson-site/test_lessons.js dev/simstudy/lesson-site/build_preview.R _config.yml .Rbuildignore assets/js/lessons.js
git status --short
git commit -m "Renumber the lessons 0 to 7: titles, figure prefixes, registrations and script paths

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Update every live cross-reference and living document

**Files:**
- Modify: `vignettes/occJSDM-lesson-{0..7}.Rmd` and their `.md`; `vignettes/LESSON-PLAN.md`; `TODO.md`; `AGENTS.md` line 317; `dev/simstudy/vignette-lesson/README.md`; `dev/simstudy/jsdm-package-comparison/README.md`; `dev/simstudy/jsdm-package-comparison/extension/README.md`; `dev/simstudy/jsdm-package-comparison/extension/CALIBRATION.md`; `dev/simstudy/jsdm-package-comparison/FITTING-REPORT.md` line 3 (link target only)
- Check, not edit: `README.md`, `vignettes/occJSDM.Rmd`
- Read: `renumber-inventory.md` section 1 (1.2 to 1.6.1 and the [hist] rules), section 6 items 2 to 5

**Interfaces:**
- Consumes: Task 6's file names.
- Produces: a tree in which every live link names its final file and an existing heading, and every live "Lesson N" means the final lesson N; the Phase B entry in the decision log of `vignettes/LESSON-PLAN.md`.

Old-to-new rules for this task (old numbers on the left):
- Files: `occJSDM-lesson-intuition` to `-1`; `-1` to `-2`; `-2` to `-7`; `-3` stays; `-prediction` to `-4`; `-4` to `-5`; `-extension` to `-6`; `-0` stays.
- Prose: "Lesson 1" to "Lesson 2"; "Lesson 2" to "Lesson 7"; "Lesson 3" stays unless the sentence is about new-site prediction, WAIC model comparison or the one-factor fit, which becomes "Lesson 4"; "Lesson 4" becomes "Lesson 5", or "Lesson 6" where the sentence is about the ten communities, sections 12 to 15, bias and coverage across communities, or the calibration results; "the intuition lesson" becomes "Lesson 1"; "the prediction lesson" becomes "Lesson 4"; "the ten-community lesson" becomes "Lesson 6"; "the four-JSDM comparison" (as a lesson name) becomes "Lesson 5"; "the spatial lesson" stays as a role name. Ranges and lists follow the same map: "Lessons 0 and 1" to "Lessons 0 and 2"; "Lessons 0, 1 and 3" to "Lessons 0, 2 and 3"; "Lessons 0 to 4" to "Lessons 0 to 7".
- Historical text keeps its numbers (see "Decisions this plan takes").

- [ ] **Step 1: Rewrite file links in one mechanical pass**

Run once only; a second run would map 1 to 2 to 7. If it must be redone, `git checkout` the files first.

```sh
perl -pi -e 's{occJSDM-lesson-(intuition|prediction|extension|1|2|4)(?=\.md|\.Rmd|\.html|_files)}{"occJSDM-lesson-" . {intuition => 1, prediction => 4, extension => 6, 1 => 2, 2 => 7, 4 => 5}->{$1}}ge' vignettes/occJSDM-lesson-[0-7].Rmd dev/simstudy/vignette-lesson/README.md dev/simstudy/jsdm-package-comparison/README.md dev/simstudy/jsdm-package-comparison/FITTING-REPORT.md
perl -pi -e '$done = 1 if /^### Decisions made after this update/; s{occJSDM-lesson-(intuition|prediction|extension|1|2|4)(?=\.md|\.Rmd|\.html|_files)}{"occJSDM-lesson-" . {intuition => 1, prediction => 4, extension => 6, 1 => 2, 2 => 7, 4 => 5}->{$1}}ge unless $done' vignettes/LESSON-PLAN.md
git diff --stat
```

Then by hand: in `dev/simstudy/jsdm-package-comparison/extension/README.md` line 11 and `extension/CALIBRATION.md` lines 30 and 31, `vignettes/occJSDM-lesson-4.Rmd` becomes `vignettes/occJSDM-lesson-6.Rmd` (these describe the extension, so the perl map's 5 would be wrong; they were left out of the pass). In `dev/simstudy/vignette-lesson/README.md`, check that the render commands now name the right lessons (Lesson 1's old commands now name `occJSDM-lesson-2.Rmd`, the old Lesson 3 commands still `occJSDM-lesson-3.Rmd`) and that dated records (sections whose heading carries a date) are unchanged: revert any perl change inside them with an Edit.

- [ ] **Step 2: Check every link target and heading**

```sh
Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/check-lesson-links.R vignettes/occJSDM-lesson-[0-7].Rmd vignettes/occJSDM.Rmd vignettes/LESSON-PLAN.md TODO.md dev/simstudy/vignette-lesson/README.md dev/simstudy/jsdm-package-comparison/README.md dev/simstudy/jsdm-package-comparison/extension/README.md dev/simstudy/jsdm-package-comparison/extension/CALIBRATION.md dev/simstudy/jsdm-package-comparison/FITTING-REPORT.md
```

Expected: the success line. A "no such heading" line names a link whose content moved: point it at the lesson that now holds the heading. Known cases already handled: the spatial lesson's link to `#predict-occupancy-at-genuinely-new-sites` (now Lesson 4), `LESSON-PLAN.md`'s link to `#which-fits-passed-remained-flagged-or-failed` (now Lesson 6), and the intuition lesson's links to `#exercise-predict-one-species-given-another` (now Lesson 5). A "no such lesson" line naming `occJSDM-lesson-intuition`, `-prediction` or `-extension` is a missed live reference: fix it by the file rule.

- [ ] **Step 3: Renumber the prose in each lesson**

For each of `vignettes/occJSDM-lesson-0.Rmd` to `-7.Rmd`, list the candidates and edit each line by the rules above:

```sh
/usr/bin/grep -n -E "Lessons? [0-9]|intuition lesson|prediction lesson|ten-community lesson|four-JSDM comparison|[Ss]ections? 1[0-5]" vignettes/occJSDM-lesson-0.Rmd
```

Known lines from `renumber-inventory.md` 1.2 (old file, old line): lesson-0 lines 26, 28, 69, 95, 125, 161, 163, 167, 187, 280, 325, 327, 366, 406, 506, 529, 568, 570, 572 (heading "Take the right objects into Lesson 1" becomes "... Lesson 2"), 574, 580, 586; intuition lines 20, 48, 56, 58, 60, 122, 180, 324, 326, 403, 405, 407, 417 to 422; lesson-1 lines 26 ("This is **Lesson 1**" becomes "**Lesson 2**"; "The spatial lesson, [Lesson 2]" becomes "[Lesson 7]"), 200, 334, 519, 566, 923, 982, 986, 1032; lesson-2 lines 46 ("[Lesson 3]" becomes "[Lesson 4]", since it links the prediction section), 62, 74, 92, 221, 427, 574, 600; lesson-3 lines 31, 33, 37, 480, 721, 778, 780, 1107, 1174, 1294 and the finder bullets; lesson-4 (now 5) lines 38, 668, 855. Also the new lessons' own references: Lesson 4's "Lesson 3" links stay, its "[the four-JSDM comparison](occJSDM-lesson-5.md)" text becomes "[Lesson 5](occJSDM-lesson-5.md)"; Lesson 6's "[The four-JSDM comparison](occJSDM-lesson-5.md)" becomes "[Lesson 5](occJSDM-lesson-5.md)"; Lesson 3's "[the prediction lesson](occJSDM-lesson-4.md)" becomes "[Lesson 4](occJSDM-lesson-4.md)"; Lesson 5's "[The ten-community lesson](occJSDM-lesson-6.md)" becomes "[Lesson 6](occJSDM-lesson-6.md)"; Lesson 5's section 3 "[the intuition lesson](occJSDM-lesson-1.md#...)" becomes "[Lesson 1](occJSDM-lesson-1.md#...)". A heading whose text changes (only lesson-0's "Take the right objects into ...") has no inbound link (inventory 1.3); confirm with Step 2's checker after editing.

- [ ] **Step 4: `vignettes/LESSON-PLAN.md`**

Only lines above `### Decisions made after this update` change.
  - Replace the status pipe table (`| Lesson | Status | What it teaches | Source |` and its rows) with nine bullets in reading order, each bullet giving the lesson title (or "Quickstart") in bold, its file name in a code span in brackets, then its status and what it teaches as sentences, taken from the old row's cells, with new bullets for Lesson 1 ("Built in Phase A of the lesson rewrite, merged in PR #19; Phase B numbered it." and what it teaches from the intuition lesson's opening), Lesson 4 ("Split from Lesson 3 in Phase B." and the old Lesson 3 row's prediction clause) and Lesson 6 ("Split from the four-JSDM lesson in Phase B." and the old Lesson 4 row's extension clause); Lessons 3 and 5 lose the clauses that moved.
  - Lines 28, 34, 36: renumber by the rules; line 28's "Lesson N was renamed Lesson 4 on 23 September 2026" is history and stays.
  - `### Lesson 2: spatial landscapes and survey design` becomes `### Lesson 7: spatial landscapes and survey design`; `### Lesson 4: occJSDM, gllvm, sjSDM and Hmsc` becomes `### Lessons 5 and 6: occJSDM, gllvm, sjSDM and Hmsc`; lines 54, 62, 64 ("Lesson 4 sections 14-15" becomes "Lesson 6 sections 3 and 4"; PR #14 also edits this line), 68, 76, 82 to 85 by the rules. Archive paths that contain `lesson-4-` keep it.
  - Replace the migration pipe table (`| Original content | Current home and status | Remaining work |`) with one bullet per row giving the original content in bold, then its current home and status, then "Remaining work:" and the remaining work, renumbering lesson references by the rules (the prediction row's "Lesson 3" becomes "Lesson 4").
  - Append to the decision log, as its last entry, with DATE replaced by the output of `date "+%-d %B %Y"`: "- **DATE, Phase B restructure and renumbering:** combined Phase B pull request (link added when opened). The lessons are now numbered 0 to 7: 0 simulate (unchanged), 1 what occupancy models and JSDMs do (was `occJSDM-lesson-intuition`), 2 fit the model (was 1), 3 outputs (was 3, without prediction), 4 predict at new sites and compare models (split from 3), 5 compare four JSDMs (was 4, sections 1 to 11), 6 ten communities and traits (was 4, sections 12 to 15 and the reduced full results), 7 spatial landscapes and survey design (was 2). Lesson 3 now starts with diagnostics and shows the package's own plot with truth as each section's main figure; gradient 2, the perfect-observation fit, variation partitioning and the ordinary biplot were re-exported from the three-sample archive with new verifier checks, and the ordination contribution and the variation scatter stay custom because the package plots cannot draw those comparisons. Lesson 3's simulation-only evidence, mirror study, internal-slot trace reading, package history and reproduction commands are in its appendix. Lesson 5's matrix-factoring section is a reminder pointing at Lesson 1, its exercise follows section 4, and its fitting calls, sjSDM two-optima evidence and reproduction record are in its appendix. Lesson 6 no longer hides every later chunk's code. Dated entries above keep the numbers they were written with; the prose passes on Lessons 3 to 7 follow in Phase C."
  - Check: `git grep -n "^|" -- vignettes/LESSON-PLAN.md | wc -l` gives 0.

- [ ] **Step 5: `TODO.md` and `AGENTS.md`**

  - `TODO.md` line 114: "currently lessons 0-4" becomes "currently lessons 0-7". Line 115: "Lessons 0 to 4" becomes "Lessons 0 to 7" and "for Lesson 4, `lesson-4-calibration-tables.R`" becomes "for Lesson 6, `lesson-4-calibration-tables.R`".
  - Line 116, "**Move Lesson 2 to the end of the sequence ...**": strike the bold title with `~~`, insert after it (DATE as in Step 4) "**DONE DATE in Phase B of the lesson rewrite: the spatial lesson is Lesson 7, in the 0 to 7 numbering recorded in `vignettes/LESSON-PLAN.md`.**", and wrap the rest of the item's original text in `~~` so it reads as history (the format of the item two lines below it).
  - Line 117, "**Lesson 4's section order:**": the same treatment, with "**DONE DATE in Phase B: the ten-community extension is its own lesson (Lesson 6), and the comparison lesson (Lesson 5) ends with its closing section and an appendix.**"
  - Open items that name lessons (lines 71, 119, 120, 142, 146, 147 by the inventory): renumber by the rules; "the Lesson 2 sweep" becomes "the site-arrangement sweep (Lesson 7)", since the sweep's own files keep their old name. Struck-through DONE items stay as written.
  - `AGENTS.md` line 317: "Lessons 0 to 4" becomes "Lessons 0 to 7". Line 956 is in a dated section and stays.

- [ ] **Step 6: The dev documents that describe the lessons**

  - `dev/simstudy/vignette-lesson/README.md`: lines 3, 5 and 7 (the opening) describe the sequence as it now is (Lesson 0, Lesson 1 on what the models do, Lesson 2 fitting, Lesson 3 outputs, Lesson 4 prediction, Lessons 5 and 6 the four-package comparison, Lesson 7 the spatial sweep), replacing the stale "future outline" wording; line 48's render advice names `occJSDM-lesson-0.Rmd` and `occJSDM-lesson-7.Rmd`; line 66's counts become "the nine lesson sources" and "the nine committed `.md` files". Dated records stay.
  - `dev/simstudy/jsdm-package-comparison/README.md` lines 5 and 65: "Lesson 4" becomes "Lesson 5" (the link targets were set in Step 1).
  - `extension/README.md` line 1: `# Lesson 4 replicated extension` becomes `# Lesson 6: replicated extension`; line 5: "Lesson 4 sections 14-15" becomes "Lesson 6 sections 3 and 4". Lines 9, 10, 21 and 22 are historical worktree and archive paths and stay.
  - `extension/CALIBRATION.md` line 3: "this Lesson 4 extension" becomes "this extension (now Lesson 6)". Archive paths stay.
  - `FITTING-REPORT.md`: only the link target changed (Step 1); its prose stays.
  - `README.md` and `vignettes/occJSDM.Rmd`: confirm with `/usr/bin/grep -n -E "Lessons? [0-9]|occJSDM-lesson" README.md vignettes/occJSDM.Rmd` that they still name no lesson number or file (expected: no output).

- [ ] **Step 7: Re-render, check and commit**

Render all eight lessons with the markdown command and the HTML command; apply the PNG policy to the folders and prefixes of Task 6 Step 5 (expected: all reverted). Run Step 2's checker again, `Rscript dev/simstudy/vignette-lesson/test_lesson_links.R` and the site test, and the Global Constraints counts. Then:

```bash
git add vignettes TODO.md AGENTS.md dev/simstudy/vignette-lesson/README.md dev/simstudy/jsdm-package-comparison/README.md dev/simstudy/jsdm-package-comparison/FITTING-REPORT.md dev/simstudy/jsdm-package-comparison/extension/README.md dev/simstudy/jsdm-package-comparison/extension/CALIBRATION.md
git status --short
git commit -m "Point every live cross-reference and living document at the 0 to 7 numbering

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Final verification and the pull request draft

**Files:**
- Create (not committed): `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/pr-body.md`
- No repository file changes unless a check fails, in which case the fix is committed with a message naming the check.

**Interfaces:**
- Consumes: the branch after Task 7.
- Produces: a branch ready to push, every check's output, rendered HTML of all eight lessons in `vignettes/` (gitignored) for Doug's read, and the drafted pull request body.

- [ ] **Step 1: Every verifier**

```sh
S=/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/run-step.zsh; A=/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
zsh $S 60-test_lesson Rscript dev/simstudy/vignette-lesson/test_lesson.R
zsh $S 61-verify_lesson Rscript dev/simstudy/vignette-lesson/verify_lesson.R $A
zsh $S 62-verify_outputs Rscript dev/simstudy/vignette-lesson/verify_outputs.R $A
zsh $S 63-verify_diagnostics Rscript dev/simstudy/vignette-lesson/verify_diagnostics.R $A
zsh $S 64-verify_latent_tables Rscript dev/simstudy/vignette-lesson/verify_latent_tables.R $A
zsh $S 65-verify_native_plots Rscript dev/simstudy/vignette-lesson/verify_native_plots.R $A
zsh $S 66-native-traits-verify Rscript dev/simstudy/vignette-lesson/native-traits-verify.R $A
zsh $S 67-ordination-verify Rscript dev/simstudy/vignette-lesson/ordination-verify.R $A
zsh $S 68-remaining-plots-verify Rscript dev/simstudy/vignette-lesson/remaining-plots-verify.R $A
zsh $S 69-prediction-verify Rscript dev/simstudy/vignette-lesson/prediction-verify.R $A $A/prediction
zsh $S 70-unbalanced-verify Rscript dev/simstudy/vignette-lesson/unbalanced-verify.R $A/unbalanced
Rscript dev/simstudy/jsdm-package-comparison/verify-teaching.R archive dev/simstudy/jsdm-package-comparison
Rscript dev/simstudy/jsdm-package-comparison/extension/verify.R /Users/douglasyu/src/occJSDM/dev/simstudy/results/lesson-4-extension-20260923 .
Rscript dev/simstudy/spatial-design-sweep/verify-lesson.R .
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
node --test dev/simstudy/lesson-site/test_lessons.js
```

Expected last lines, in order: `test_lesson.R` passes; "All teaching-evidence checks passed."; "Trait decomposition, sampling-effort truth, source hashes and fit hashes verified."; "Diagnostic selection, fit hashes, source hashes and plotting identities verified."; "Truth joins survive row permutations and reject duplicate identities; native table builds."; "Twelve PNGs, ..."; the traits PASS lines; the ordination PASS lines; "Six PNGs, ..."; "All prediction verification checks passed."; the unbalanced "... and displayed code." line; "Teaching numerical verification passed for every check the committed archive supports."; "... 626 scored results; 640 of 640 fits attempted."; "Lesson 2 bundle verified against the committed results." (the sweep verifier's message keeps its old wording; the sweep is out of scope); "Shared lesson-link regression checks passed." after "All 9 ..."; `fail 0`. The four Lesson 3 logs each contain "Lesson 3 displays the exact exported ...".

- [ ] **Step 2: Counts and greps**

```sh
git diff --name-only 0f48dae -- ':!*.png' ':!*.rds' | xargs /usr/bin/grep -n $'\u2014' | wc -l
git log 0f48dae..HEAD --format=%B | /usr/bin/grep -c $'\u2014'
git grep -n "^|" -- 'vignettes/occJSDM*.Rmd' vignettes/LESSON-PLAN.md | wc -l
git grep -n "editor_options" -- 'vignettes/occJSDM*.Rmd' | wc -l
git grep -n -E "teaching-data/lesson-[24]-[a-z0-9-]+-1\.png" -- vignettes | wc -l
git grep -n -E "occJSDM-lesson-(intuition|prediction|extension)" -- . ':!dev/simstudy/vignette-lesson/LESSON-STYLE-PLAN-A.md' ':!dev/simstudy/vignette-lesson/LESSON-STYLE-PLAN-B.md' ':!dev/simstudy/vignette-lesson/LESSON-STYLE-SPEC.md' ':!dev/simstudy/vignette-lesson/style-diagnostic' ':!dev/simstudy/lesson-site/PLAN.md' ':!dev/simstudy/lesson-site/DESIGN.md'
git grep -n -E "Lesson (1: Fit|2: Spatial|4: Compare)" -- 'vignettes/occJSDM*' assets _config.yml dev/simstudy/lesson-site/test_lessons.js
Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/check-lesson-links.R vignettes/occJSDM-lesson-[0-7].Rmd vignettes/occJSDM.Rmd vignettes/LESSON-PLAN.md TODO.md dev/simstudy/vignette-lesson/README.md dev/simstudy/jsdm-package-comparison/README.md dev/simstudy/jsdm-package-comparison/extension/README.md dev/simstudy/jsdm-package-comparison/extension/CALIBRATION.md dev/simstudy/jsdm-package-comparison/FITTING-REPORT.md
/usr/bin/grep -n "Phase B restructure and renumbering" vignettes/LESSON-PLAN.md
wc -l vignettes/occJSDM-lesson-[0-7].Rmd
```

Expected: 0, 0, 0, 0, 0; the stem grep lists only lines inside the dated decision log of `vignettes/LESSON-PLAN.md`, struck-through DONE items of `TODO.md` and dated records of `dev/simstudy/vignette-lesson/README.md` (any other hit is a missed live reference: fix it); the title grep prints nothing (the old titles "Lesson 1: Fit", "Lesson 2: Spatial" and "Lesson 4: Compare" survive nowhere in the lessons, the site list, the site configuration or its test); the link checker's success line; one decision-log line; line counts close to the expectations in "Decisions this plan takes".

- [ ] **Step 3: HTML for Doug's read**

Render each of the eight lessons with the HTML command. Confirm `git status --short` is clean (HTML is gitignored) and list the eight `vignettes/occJSDM-lesson-*.html` paths for the controller.

- [ ] **Step 4: Draft the pull request body**

Write `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/pr-body.md` (markdown rules apply; one line per paragraph; no em-dashes) with these sections:
  - Summary: the final numbering (the "Final numbering" bullets) and one paragraph on what Phase B is (structure, moves, splits, figures, renumbering; no teaching prose, which is Phase C).
  - Per task: what changed, in two to four bullets each, naming the new figures, the two custom figures kept and why, the seven removed, the appendices, the echo leak fix, the reduced full results, and the registrations.
  - Verifier changes: the four extended verifiers with their new checks and raised chunk counts, stating that no existing check was loosened or removed; the path-only edits with file and line.
  - Deliberately left for Phase C: the items under "Deferred to Phase C".
  - Collision with PR #14: the paragraph "Known collision with PR #14" of this plan, with the final file names.
  - Needs Alex: "None. This pull request states no package fact that the spec marks unconfirmed; the `ggtern::theme_bw()` explanation, the one Lesson 3 item marked needs Alex, is Phase C." If any task added such a statement, list it here instead.
  - Checks run: the output lines of Steps 1 and 2.
  - A note that Doug reviews line by line, with the eight HTML paths.
  - The last line: `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.

- [ ] **Step 5: Hand back to the controller**

Report: the commit list (`git log --oneline 0f48dae..HEAD`), the check outputs, the PR body path, and any check that needed a fix. The branch is ready to push.

Controller steps, not for an implementer:
  - Push `codex/lesson-phase-b` and open the pull request with `gh pr create --base main --title "Lesson rewrite Phase B: restructure Lessons 3 and 4, split prediction and the extension, renumber 0 to 7" --body-file /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/pr-body.md`; then add the pull request link to the Phase B decision-log entry in `vignettes/LESSON-PLAN.md` in a follow-up commit.
  - Update the memory notes once (`/Users/douglasyu/.claude/projects/-Users-douglasyu-src-occJSDM/memory/`): add a note recording the 0 to 7 numbering and the old-to-new mapping, so later sessions do not read "Lesson 2" or "Lesson 4" with their old meanings, and index it in `MEMORY.md`; the existing `beta-waits-on-pr14.md` lines that name Lessons 2 and 4 are history and stay.

---

## Self-review notes

- **Spec coverage, Phase B items.** Intuition lesson takes its place as Lesson 1: Task 6 (file, title, site entry) and Task 7 (back-references from Lessons 0, 2, 3 and 5). Lesson 2 to the end: Task 6. Prediction out of Lesson 3: Task 2. Extension out of Lesson 4: Task 4. Lesson 4 section 3 to a two-sentence reminder pointing at the intuition lesson, keeping its exercise, which moves after section 4: Task 5 Steps 1 and 2. Appendices created in 3 and 4: Task 3 Step 7, Task 5 Step 1, and in the two new lessons Tasks 2 and 4. Lesson 3's figures switched to package plots with truth: Task 1 (re-export) and Task 3 Steps 3 and 5. Diagnostics to the front of Lesson 3: Task 3 Step 2. Lesson 4 echo leak: Task 4 Step 3, exactly as `lesson-4-map.md` says. Every cross-link, title, `VignetteIndexEntry`, PNG prefix, the site list and its test: Tasks 6 and 7. `LESSON-PLAN.md` (status table to bullets, decision-log line): Task 7 Step 4. Quickstart and README: Task 7 Step 6 confirms neither names a lesson. Memory notes: Task 8 controller steps. One pull request verified by `test_lesson_links.R`, the site test and every lesson verifier: Task 8.
- **Spec coverage, Lesson 3 and 4 bullets placed in Phase B.** Lesson 3: trait cancellation, Procrustes alignment, mirror study (short warning plus the per-chain check in the text) and PR #13 history to the appendix (Task 3 Steps 2, 3, 6, 7); traces through `returnConvergenceDiagnostics()` and `plotTraceplot()` in the text with internal slots in the appendix and the planned accessor (Task 3 Steps 2, 6, 9); reproduction separated from the function finder (Steps 3, 7, 8); three pipe tables to bullets (Steps 7, 8); three uncited references removed (Step 8); three Lesson 2 statements (line 31 in Step 6, line 784 in Task 2 Step 5 and the prediction lesson's opening, the finder row in Step 8); the `sampleresults` paragraph (Step 6); `editor_options` (Step 6); the appendix's unique guidance folded into the sections (Step 3). Lesson 4: section 10 and the reproduction record to the appendix, two-optima point kept with its evidence in the appendix, pipe tables, `editor_options` (Task 5); full results reduced (Task 4 Step 4). The remaining Lesson 3 and 4 bullets are teaching-prose items and belong to Phase C (listed under "Deferred to Phase C").
- **Owner rulings.** Phase B now, before PR #14 and the beta tag: the collision section and Task 8's PR body. Eight lessons 0 to 7: "Final numbering". Path-only verifier and `finish-study.py` edits: Task 2 Step 7, Task 6 Step 4. Four figure groups re-exported with new checks of the same standard: Task 1 (each new check re-derives plotted numbers from all draws of the relevant fit and truth independently; displayed code stays compared byte for byte; count literals only rise). Historical documents keep old numbers and data file names keep `lesson-4-`: "Decisions this plan takes". `editor_options` removed from Lessons 3 and 4: Task 3 Step 6, Task 5 Step 4. Extension dev copies left unchanged: not touched by any task.
- **Decomposition.** The suggested five tasks became eight: each split and its restructure are separate tasks (2 and 3, 4 and 5) because a reviewer can approve a pure move while rejecting a restructure, and the split commits are checkable mechanically with the sorted-line diff; renumbering is two tasks (6 and 7) because file renames and registrations are mechanical and testable by the existing tests, while the cross-reference pass needs judgement per sentence and its own link checker. Tasks 3 and 5 each make a moves-only commit first, so `git diff --color-moved` shows the restructure as moves.
- **Hash ordering.** Task 1 commits the snippets and exporters (Step 9) before running the exporters (Step 10), and commits the verifiers with the bundles afterwards (Step 15); no later task edits a snippet or exporter, and the PNG revert script is never pointed at exporter output.
- **Placeholder scan.** No "TBD", "similar to" or bracketed fill-ins: the perfect-fit blocks are written out in full for each exporter and verifier, and the extension lesson's theme lines and future-work sentence are quoted. The `<stem>` in the render commands and DATE in the TODO and decision-log text are substitutions (DATE from `date`), not open choices.
- **Name consistency.** Provisional stems `occJSDM-lesson-prediction` and `occJSDM-lesson-extension` are used identically in Tasks 2, 4, 6 and 7 and in the perl map; site labels "Prediction" and "Extension" match `test_lessons.js` in both tasks; `fitmodel_perfect`, `environment_2_truth`, `remaining_truth$variation`, `plots$ordinary_biplot` and `provenance$perfect_fit_manifest` are defined in Task 1 before any verifier or lesson uses them; the anchors `#reproduce-the-teaching-figures`, `#read-trace-draws-from-the-fit-object`, `#the-mirror-labelling-study`, `#a-real-cancellation-inside-this-simulated-community`, `#the-two-sjsdm-optima-evidence` and `#reproduction-record` are the slugs of headings created in the same task.
- **Choices a reviewer may question.** Keeping `calibration-data` in Lesson 6's appendix departs from `lesson-4-map.md`, which recommended cutting it; it stays because three of its columns appear nowhere else, which is the spec's own test for the reduced appendix. Removing the custom `detection-effort` figure departs from the map's suggestion to keep it; the package's effort plots already carry exact truth for the same question, and the formula chunk keeps the analytic expectation in the text.

# Port the old vignettes into the lessons: Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Move everything useful that `vignettes/ORIG_occJSDM.Rmd` and `vignettes/simulateOccJSDMData.Rmd` still teach, and Lessons 0 to 7 do not, into the lessons and the simulator's help page; delete `ORIG_occJSDM.Rmd` now, and record that `simulateOccJSDMData.Rmd` is deleted when the lessons are published.

**Architecture:** Lesson 0 gains a section on how traits build each species' environmental responses (from the saved truth) and an appendix that runs the simulator live with its spatial field switched on. The simulator's roxygen help is expanded and its Rd regenerated. Lesson 3 gains a `plotCovariateEffect()` figure through a new exporter and verifier pair modelled on `remaining-plots-export.R` and `remaining-plots-verify.R`, which writes one new PNG and one new compact bundle and leaves every existing bundle unchanged. Lessons 2, 3 and 4 gain one to three sentences each. One branch, one pull request.

**Tech Stack:** R 4.5 with rmarkdown, knitr, dplyr, tidyr, tibble, ggplot2, roxygen2 8.0.0 (matches `RoxygenNote` in `DESCRIPTION`), testthat; the occJSDM library in the full-fit archive for the exporter and verifiers; git.

**Spec:** `dev/simstudy/vignette-lesson/LESSON-STYLE-SPEC.md` for how lessons are written ("The reader", "The criteria, applied", "What does not change", "Verification for every lesson pull request"). Doug's decisions of 4 October 2026, recorded in "Decisions" below, set the scope and override the spec where they differ. The audit behind them is in this conversation's record and summarised per task below.

## Decisions (Doug, 4 October 2026)

- **Goal.** Port everything useful from the two old vignettes so they can be deleted without loss.
- **Keeper 1, spatial simulation.** Into a Lesson 0 appendix and the simulator's help page.
- **Keeper 2, variation partitioning.** Omitted; it stays a `TODO.md` item from Phase B.
- **Keeper 3, threshold advice.** One sentence in Lesson 2.
- **Keeper 4, the trait-generating formula.** Taught, with the coefficient-against-trait scatter (minor gap 3) integrated.
- **Minor gap 1, `plotCovariateEffect()`.** Checked on 4 October 2026 (unit tests pass; medians match a hand calculation from the default fit's draws to 6e-17; the generating curve lies inside the 95% band everywhere for 6 of 10 species). Shown in Lesson 3 through a new exported figure: Doug approved extending the export and verification machinery and committing new PNGs, with existing bundles unchanged.
- **Minor gap 2.** Two or three sentences beside Lesson 3's primer plots on reading low laboratory detection rates (primer mismatch, species that shed little DNA).
- **Minor gap 4.** Half a sentence in Lesson 4 on `useEnvCov = FALSE`.
- **Minor gap 5, `plotBiplot()` arguments.** Skipped; the help page covers them.
- **Minor gap 6, `sigma_bs` and `sigma_ts`.** Lesson 0 says they have no effect in the simulator; the pull request's "needs Alex" list asks whether to drop or document the arguments.
- **Deletion.** `ORIG_occJSDM.Rmd` is deleted in this pull request. `simulateOccJSDMData.Rmd` stays, because during the beta it is the only simulation guide that ships (the lessons are withheld from the package and the site); `TODO.md`'s publishing item gains the step of deleting it and pointing the quickstart and README to Lesson 0 when Lesson 0 is published.

## Rulings this plan takes

- **Live simulation in Lesson 0.** The spatial appendix runs `occJSDM::simulateOccJSDMData()` while knitting. Checked while planning: with the lesson's settings and `useSpatField = TRUE` it takes about 0.2 seconds, and re-seeding with the saved seed and the saved settings reproduces the saved `data_list` exactly with the current package. Every number the prose states about the live run is computed inline, because random draws can change with the package or R version.
- **Where the trait formula goes.** In Lesson 0's main body, in "Specify the ecological model", because that section already says responses are built from traits; Lesson 3's trait section links to it rather than repeating it. Checked while planning: `t(Tr %*% G + A %*% C + Bt)` equals the saved `B` (`all.equal` TRUE), with `Tr` the raw `traits` table, `G` = rows (-1, 0) and (-1, 1), and `C` = (1, 0).
- **A new exporter pair, not an edit of an old one.** Adding the figure to `remaining-plots-export.R` would change `remaining-plots-data.rds` and the hashes its verifier checks. A new pair (`covariate-effect-examples.Rmd`, `covariate-effect-export.R`, `covariate-effect-verify.R`) leaves every existing bundle and verifier byte-identical.
- **Help page.** Roxygen comments in `R/simulateData.R` only; no code changes. The Rd is regenerated with `roxygen2::roxygenise()`, and only `man/simulateOccJSDMData.Rd` may change.

## Global Constraints

- Work only in the worktree `/Users/douglasyu/src/occJSDM-worktrees/lesson-port`, on branch `codex/lesson-port` (created from main at fbe3ed6). Never work in `/Users/douglasyu/src/occJSDM`, never commit to main, never push (the controller pushes), never force anything.
- The full-fit archive is `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples`. The exporter and verifiers that read fits run from the worktree root through the step runner below, so each run is logged in the archive's `BUILD-LOG.md`. Agent shells keep no variables between calls, so commands spell paths out in full.
- Commit messages end with a blank line and `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.
- Markdown rules for every `.md` and `.Rmd` touched, from `CLAUDE.md` and `AGENTS.md`: one line per paragraph, no hard wrapping; no em-dashes anywhere, including commit messages and roxygen (use `--`); no pipe tables; escape a literal `>` `<` `~` `|` in prose as `\>` `\<` `\~` `\|` and write `-\>` for an arrow, never inside code spans or fenced blocks; keep inline code spans short; never add an `editor_options` block; no double blank lines. Re-read a file immediately before and after editing it (it may be open in RStudio); if an edit looks corrupted or duplicated, rewrite the whole file with Write.
- Writing rules, from the spec's "The reader" and "The criteria, applied": the reader is an ecologist who runs eDNA or presence/absence surveys, knows basic R and the pipe, and has not fitted an occupancy model or a JSDM. Name each concept in the reader's terms before using it, read every printed output in the reader's terms, spell out the reasoning behind every choice, and give practical advice for the reader's own survey. Do not gloss routine dplyr verbs. Match the surrounding lesson's voice and density; additions are short.
- Numbers: every number the prose states that a bundle or a live chunk holds is computed inline with `` `r ...` ``. A number nothing holds is stated with its source (file and line), or not stated.
- Out of scope (stop and ask if one would need a change): package code other than the roxygen block of `simulateOccJSDMData()` (`R/`, `src/`, `DESCRIPTION`, `NAMESPACE`); every existing teaching bundle and exporter-written PNG under `vignettes/teaching-data/`; every existing exporter, verifier and frozen snippet under `dev/simstudy/vignette-lesson/`; everything under `dev/simstudy/spatial-design-sweep/` and the hash-stamped sources under `dev/simstudy/jsdm-package-comparison/`; historical documents (completed plans, `style-diagnostic/`, dated decision-log entries, struck-through `TODO.md` items).
- The 36 frozen chunks listed in `check-frozen-chunks.R` (below) stay byte-identical, header line included. After Task 3, the new `covariate-effect-*` chunks in Lesson 3 must also stay identical to `covariate-effect-examples.Rmd`, which `covariate-effect-verify.R` checks.
- Headings that other files link to keep their exact text; `check-lesson-links.R` must pass after every task.
- Line numbers in this plan are from main at fbe3ed6. Earlier tasks shift lines; find the quoted heading, chunk label or sentence.
- Shell is zsh and `grep` and `ls` are aliased: use `/usr/bin/grep` and `/bin/ls`.
- Render commands, run from the worktree root (replace `<stem>`):

```sh
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/<stem>.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/<stem>.Rmd", output_format="rmarkdown::html_vignette")'
```

Both must finish without a new warning (`/usr/bin/grep -c "built under" vignettes/<stem>.md` is 0). HTML output is gitignored and never committed.

- PNG policy. A render rewrites every figure it draws. Keep a changed PNG only when its chunk changed on purpose or is new; revert the others with `revert-unchanged-png.zsh`. Four figures vary between renders for reasons unrelated to content: `vignettes/occJSDM-lesson-0_files/figure-gfm/environmental-maps-1.png`, `vignettes/occJSDM-lesson-2_files/figure-gfm/probability-and-state-1.png`, `vignettes/occJSDM-lesson-3_files/figure-gfm/mirror-labelling-chains-1.png` and `vignettes/occJSDM-lesson-3_files/figure-gfm/trait-cancellation-1.png`; if one is reported "kept, content changed" and its chunk did not change, restore it with `git checkout -- <path>`. Check `git status --short vignettes/` before every commit.
- Checks that must pass at the end of every task (`<files>` are the Rmd and md files the task changed):

```sh
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
node --test dev/simstudy/lesson-site/test_lessons.js
Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-PORT-PLAN/check-frozen-chunks.R
Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/check-lesson-links.R vignettes/occJSDM-lesson-[0-7].Rmd vignettes/occJSDM-lesson-[0-7].md vignettes/occJSDM.Rmd vignettes/LESSON-PLAN.md TODO.md
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-PORT-PLAN/prose-checks.zsh <files>
git diff --name-only fbe3ed6 -- ':!*.png' ':!*.rds' | xargs /usr/bin/grep -n $'\u2014' | wc -l
git log fbe3ed6..HEAD --format=%B | /usr/bin/grep -c $'\u2014'
```

Expected: "Shared lesson-link regression checks passed."; `fail 0`; "All 36 frozen chunks match main at fbe3ed6."; the link checker's success line; every "must be 0" count is 0, and every long sentence the task added is split or kept for a stated reason; the last two counts are 0.

- Stop and report to the controller, without working around it, if: a baseline check fails before the task changed anything; a verifier fails for any reason other than an accidental change to a frozen chunk; an out-of-scope file would need a change; or the spec, the Decisions and this plan conflict.

---

## Shared tools (controller creates them before Task 1; kept outside the repository)

All live in `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-PORT-PLAN/`.

- `run-step.zsh` and `revert-unchanged-png.zsh`: copies of Phase C's (`/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/`) with every `lesson-phase-c` replaced by `lesson-port`. Usage: `zsh run-step.zsh STEPNAME command args...` logs to the archive's `logs/STEPNAME.log` and appends a record to its `BUILD-LOG.md`; `zsh revert-unchanged-png.zsh PATHSPEC` reverts every modified tracked PNG under the pathspec whose pixels equal the index copy. Step names start with `p` (for example `p31-covariate-effect-export`).
- `check-frozen-chunks.R`: Phase C's, with every `7662f9a` replaced by `fbe3ed6`.
- `prose-checks.zsh`: Phase C's, unchanged.
- `check-lesson-links.R` is used from `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/` unchanged.

## File structure

- `vignettes/occJSDM-lesson-0.Rmd`, `.md`, `occJSDM-lesson-0_files/`: Task 1.
- `R/simulateData.R` (roxygen block only), `man/simulateOccJSDMData.Rd`: Task 2.
- `dev/simstudy/vignette-lesson/covariate-effect-examples.Rmd`, `covariate-effect-export.R`, `covariate-effect-verify.R` (new); `vignettes/teaching-data/covariate-effect-data.rds` and `covariate-effect-1.png` (new); `vignettes/occJSDM-lesson-3.Rmd` and `.md`; `dev/simstudy/vignette-lesson/README.md` (reproduction list): Task 3.
- `vignettes/occJSDM-lesson-2.Rmd`, `-3.Rmd`, `-4.Rmd` and their `.md`: Task 4.
- `vignettes/ORIG_occJSDM.Rmd` (deleted), `.Rbuildignore`, `vignettes/LESSON-PLAN.md`, `TODO.md`, `AGENTS.md`: Task 5.

---

### Task 1: Lesson 0, the trait formula and the spatial-simulation appendix

**Files:**
- Modify: `vignettes/occJSDM-lesson-0.Rmd`; render `vignettes/occJSDM-lesson-0.md` and its new figures under `vignettes/occJSDM-lesson-0_files/figure-gfm/`.

**Interfaces:**
- Produces: a Lesson 0 heading for the trait formula, `## How traits build each species' responses` (Task 3 and Task 4 may link to it as `occJSDM-lesson-0.md#how-traits-build-each-species-responses`), and the appendix heading `### Simulate a survey with a spatial field`.

**Source material.** `vignettes/simulateOccJSDMData.Rmd` lines 36 to 41 (the `useSpatField` switch; coordinates drawn either way), 85 to 120 (`ds`, `sigma_bs`, `sigma_ts`, `sigma_s`, `l_s` comments), 179 to 280 (the on and off runs, `varPart`, the occupancy and `eta` maps, the tip to raise `l_s` or lower `sigma_h` and `sigma_b`), and 282 to 325 (the trait section). The simulator code is `simulateData()` in `R/jsdmfun.R`: coefficients at lines 733 to 752 (`sampleEffects()` at 675 draws from -1, 0 and 1), the spatial field at 783 to 806, `eta <- eta + spatField` at 827. Read `varPart`'s computation in the simulator before explaining its columns; do not repeat the old vignette's description of its "Total" column without checking it.

- [ ] **Step 1: Baseline checks.** Run the Global Constraints checks on the untouched branch and record the output in the report. Stop if any fails.

- [ ] **Step 2: Correct the settings comments.** In the `ecological-settings` chunk (not frozen), change the comments on `sigma_bs` and `sigma_ts` to say they have no effect in the simulator (for example `# No effect: the simulator never uses it.`), and on `ds`, `sigma_s` and `l_s` to say they are spatial settings, off here, explained in the appendix. In the paragraph that begins "The most important choice for this lesson is", say in one sentence that `sigma_bs` and `sigma_ts` are required arguments that change nothing (the simulator sets the residual spatial coefficients to zero and never reads `sigma_ts`), and link to the new appendix for the spatial settings. Also say, in a clause, that the site coordinates are drawn whether or not the spatial field is on.

- [ ] **Step 3: Add `## How traits build each species' responses`.** Place it after the paragraph about `useSpatField` and before "The simulator draws the environmental values" (or move that paragraph after the new section if the flow reads better; say which in the report). Content, in the lesson's voice:
  - The formula in words and in one display line: each species' coefficients for the environmental covariates are the sum of three parts, the measured traits times a trait-effect matrix `G`, the unmeasured traits `A` times their effect matrix `C`, and a species' own remainder `Bt`; in matrix form `B = t(Tr %*% G + A %*% C + Bt)`. Name each matrix's shape in this survey, with inline numbers.
  - How each part is drawn: `G`'s entries are -1, 0 or 1 at random, so a measured trait either raises, lowers or does not change a response; `A` is standard normal with `gt` columns; `C` has 1 on its diagonal and 0 below it, so the first unmeasured trait always acts on the first covariate; `Bt` is normal with standard deviation `sigma_b`. Read the saved `G` and `C` for this survey from `known_truth$jsdmParams_true` in a run chunk and state what they mean in words (which trait raises or lowers which response, which trait has no effect on which covariate).
  - A run chunk that rebuilds `B` from its parts with the saved truth and `survey_data$traits` and checks it equals the saved `B` (`all.equal()`), so the reader sees the formula is exactly what generated the data.
  - The scatter from `simulateOccJSDMData.Rmd` lines 313 to 325, rebuilt from the saved truth: each species' true coefficient for one covariate against one measured trait whose `G` entry is nonzero, chosen in code (not hard-coded), with the line the trait part alone predicts (`Tr %*% G` for that column) drawn beside the points. Read it: the slope's sign follows the `G` entry, and the points scatter around the line because of the unmeasured-trait part and the remainder. Add the old file's point that with a zero `G` entry no trend would be expected.
  - One sentence of practical meaning: this is the structure `runOccJSDM()` estimates when given traits and `n_lattrait`, and Lesson 3's trait section asks how well ten species can recover `G`; link to it (`occJSDM-lesson-3.md#traits-ask-a-harder-different-question`).

- [ ] **Step 4: Add `### Simulate a survey with a spatial field` to the appendix.** Place it as the last subsection of "Appendix: evidence and reproduction". Content:
  - Why it is here: this lesson's survey has no spatial field; this subsection shows how to switch one on, for a reader who wants spatial test data. Fitting spatial models is taught in the spatial lesson (link to Lesson 7 by what it teaches).
  - What the field is, in the reader's terms: a smooth random surface over the site coordinates added to each species' occupancy score, so nearby sites share more of it. Then the settings: `useSpatField = TRUE`; `sigma_s`, the field's variance at a site (the kernel equals `sigma_s` at zero distance, so it is a variance, not a standard deviation); `l_s`, the distance over which sites stay similar, in coordinate units, with coordinates drawn uniformly on the unit square; `ds`, with 0 giving one field shared by every species and a positive value giving each species its own field, correlated across species at that rank.
  - A run chunk that reruns the simulator twice with the lesson's own `survey_settings`, `observation_settings` and `ecological_settings` and the saved seed, once as saved and once with `useSpatField = TRUE` (`modifyList()`), calling `occJSDM::simulateOccJSDMData()` with `model = "two_stage"`. Show and read both `varPart` tables in the reader's terms, with inline numbers (for example the range of the spatial share across species with the field on, and that it is zero with it off). Explain every column you show, checked against the simulator code.
  - One map of the continuous occupancy score (`true_params$jsdmParams_true$eta`) for one species with the field on, plotted at its site coordinates, and one sentence on why the 0/1 occupancy map shows the pattern less clearly (it is thresholded, and the environment and hidden factors add variation that is not spatial).
  - The practical tip: to make the spatial signal dominate, raise `l_s` for broader patches, or lower `sigma_h` and `sigma_b` to reduce the non-spatial variation.
  - The chunk must not change any object the main lesson uses later (it is in the appendix, after everything else; use new object names).

- [ ] **Step 5: Render and check.** Render Lesson 0 with both commands. Apply the PNG policy (new figures are kept; `environmental-maps-1.png` is render noise). Run the Global Constraints checks with `<files>` = `vignettes/occJSDM-lesson-0.Rmd`, plus `Rscript dev/simstudy/vignette-lesson/unbalanced-verify.R /Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples/unbalanced` through `run-step.zsh` as step `p15-unbalanced-verify` (it runs two Lesson 0 chunks). Read the rendered `.md` sections you added and check every sentence against the printed output.

- [ ] **Step 6: Commit.**

```bash
git add vignettes/occJSDM-lesson-0.Rmd vignettes/occJSDM-lesson-0.md vignettes/occJSDM-lesson-0_files
git commit -m "Lesson 0: teach the trait formula and how to simulate a spatial field

Ported from simulateOccJSDMData.Rmd, which can then be deleted when the
lessons are published. sigma_bs and sigma_ts are now described as having
no effect in the simulator.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

Write the report to the path the controller gives, listing what was added where, the varPart numbers you read, and any sentence you could not support from code or output.

---

### Task 2: The simulator's help page

**Files:**
- Modify: `R/simulateData.R` (the roxygen block above `simulateOccJSDMData <- function`, lines 6 to 30 at fbe3ed6; no code); regenerate `man/simulateOccJSDMData.Rd`.

**Interfaces:**
- Consumes: Lesson 0's new appendix heading from Task 1, which the help page may mention in words (help pages do not link to unpublished lessons; name `vignette("simulateOccJSDMData")` for the worked example while it ships).

- [ ] **Step 1: Baseline.** `Rscript -e 'devtools::test(filter = "simulation")'` passes; record the summary.

- [ ] **Step 2: Expand the roxygen.** Keep the existing title and the read-count `@details` facts that are already there (read `R/simulateData.R` for them first). Use `\describe{\item{...}{...}}` lists as `R/runOccJSDM.R` does for its list arguments. Document, each checked against the code in `R/simulateData.R` and `simulateData()` in `R/jsdmfun.R`:
  - `list_datasettings`: every element (`n`, `S`, `g`, `M`, `P`, `K`, `ncov_psi`, `ncov_theta`) with its meaning and shape (`M` per site, `K` per sample and primer, as the code reads them).
  - `list_params`: every element, including the optional read-intensity ones already described.
  - `list_jsdmParams`: every element, including `tau` and `useSpatField`, which the current text omits. Say that `sigma_bs` and `sigma_ts` are accepted but have no effect on the simulated data. For the spatial field: squared-exponential kernel `sigma_s * exp(-dist^2 / (2 * l_s^2))`, so `sigma_s` is the field's marginal variance; `l_s` in coordinate units, coordinates drawn uniformly on the unit square and drawn whether or not the field is on; `ds = 0` gives one shared field, `ds > 0` one field per species correlated at rank `ds`; the field is added to the occupancy linear predictor only when `useSpatField = TRUE`.
  - A `@details` paragraph on how species' environmental coefficients are built: `B = t(Tr %*% G + A %*% C + Bt)`, with `G` and the free entries of `C` drawn from -1, 0 and 1, `A` standard normal, `Bt` normal with standard deviation `sigma_b`.
  - `@return`: the top-level structure, `data_list` (the input for `runOccJSDM()`) and `true_params`, naming the truth elements a user compares with (`z_true`, `w_true`, `jsdmParams_true$B`, `$eta`, `$varPart`); check the names against a real call.
  - A sentence pointing to `vignette("simulateOccJSDMData", package = "occJSDM")` for a worked example.
  No em-dashes; use `--` only where a dash is needed.

- [ ] **Step 3: Regenerate.** `Rscript -e 'roxygen2::roxygenise()'`. Then `git status --short`: only `R/simulateData.R` and `man/simulateOccJSDMData.Rd` may be modified. If any other Rd or `NAMESPACE` changed, revert it with `git checkout -- <path>` and say so in the report. Run `Rscript -e 'tools::checkRd("man/simulateOccJSDMData.Rd")'` (no output expected), view it with `Rscript -e 'tools::Rd2txt("man/simulateOccJSDMData.Rd")'` and read it once as a user would.

- [ ] **Step 4: Tests.** `Rscript -e 'devtools::test(filter = "simulation")'` passes; then run the Global Constraints link and em-dash checks.

- [ ] **Step 5: Commit.**

```bash
git add R/simulateData.R man/simulateOccJSDMData.Rd
git commit -m "Document every simulator setting, the spatial field and the trait formula

Help-page text only; no code changes. Ported from simulateOccJSDMData.Rmd.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: `plotCovariateEffect()` in Lesson 3, with its exporter and verifier

**Files:**
- Create: `dev/simstudy/vignette-lesson/covariate-effect-examples.Rmd`, `covariate-effect-export.R`, `covariate-effect-verify.R`; `vignettes/teaching-data/covariate-effect-data.rds`, `vignettes/teaching-data/covariate-effect-1.png` (written by the exporter).
- Modify: `vignettes/occJSDM-lesson-3.Rmd`, `.md`; `dev/simstudy/vignette-lesson/README.md` (the list of exporters and verifiers, if it has one).

**Interfaces:**
- Consumes: the full default fit in the archive (`default-fit.rds`, element `fit`), `vignettes/teaching-data/nonspatial-lesson.rds`, `output-lesson.rds` and `remaining-plots-data.rds` (read only), `helpers.R` and `score_lesson.R`.
- Produces: `covariate-effect-data.rds`, a list with `truth` (a tibble of the generating curve on the plotted grid: `Species`, `x` in original units, `truth`), `figures` (file, width, height, md5), `checks` (per species: the share of grid points where the truth lies inside the 95% band, and the largest absolute difference between median and truth), and `provenance` (as in `remaining-plots-export.R`: source hashes, input, lesson and output md5s, fit manifest, snippet, exporter and student-code md5s, API names, package files' md5s, session).

**Model.** Copy the structure of `remaining-plots-export.R` and `remaining-plots-verify.R` (read both in full first): the same archive and identity checks, the same "displayed chunks are the executed chunks" mechanism (`eval=FALSE, purl=TRUE` chunks purled and sourced), the same provenance. The only plotted call is `plotCovariateEffect()`.

- [ ] **Step 1: The snippet.** `covariate-effect-examples.Rmd` holds, in this order: a `covariate-effect-read` chunk (`purl=FALSE`, run in the lesson) that reads `teaching-data/covariate-effect-data.rds`; a `covariate-effect-setup` chunk (`eval=FALSE, purl=TRUE`) defining the truth table and `ggtern::theme_bw(base_size = 12)` as `remaining-plots-setup` does; and a `covariate-effect-gradient-1` chunk (`eval=FALSE, purl=TRUE`) that calls `plotCovariateEffect(fitmodel, covNames = "X_psi.EnvCov.1", idx_species = 1:10)`, takes the element by name (`[["X_psi.EnvCov.1"]]`, because the function returns a named list of plots), adds the theme, the generating curve as a black dashed line and labels in the style of `remaining-plots-gradient-1`, and assigns it to `covariate_effect_1`; then an `covariate-effect-gradient-1-image` chunk (`echo=FALSE, purl=FALSE`) including `teaching-data/covariate-effect-1.png`. Only gradient 1 is shown; the prose says the same call takes the second covariate's name.

- [ ] **Step 2: The exporter.** `covariate-effect-export.R FULL_FIT_DIRECTORY`: identity checks as in `remaining-plots-export.R`; API source check for `plotCovariateEffect`, `returnCovariateEffect`, `returnCovariateEffect_base`, `covariate_response_grid`, `plot_covariate_response`, `create_covariates_matrix`; the truth on the function's own grid: for each species, `plogis(B0 + B[1, ] * (x - mean_1) / sd_1 + B[2, ] * median(standardized X_psi column 2))`, with `B0`, `B` from `input$sim$true_params$jsdmParams_true`, `mean_1` and `sd_1` from `fitmodel$infos$list_X_psi_mat`, hidden factors at zero. Checked while planning: the simulator applies `B` to the standardized covariates, the fit's `X0_psi` is already standardized, and the function's grid runs from the minimum to the maximum observed value in original units. Then purl and source the snippet, save the PNG at width 10, height 8, dpi 150, white background, record the `checks`, and save the bundle with `compress = "xz"`. Print one closing line naming what was exported.

- [ ] **Step 3: The verifier.** `covariate-effect-verify.R FULL_FIT_DIRECTORY`: provenance and identity checks as in `remaining-plots-verify.R`; the purled snippet's md5 equals the recorded student-code md5; if Lesson 3 contains `covariate-effect-*` exportable chunks, they are identical to the snippet's (print "Lesson 3 displays the exact exported covariate-effect plotting body."); an independent recomputation of the plotted medians and 95% limits from the fit's draws at every grid point for every species (pool iterations within chains; standardize `x` with `mean_1` and `sd_1`; hold covariate 2 at the median of its standardized column), agreeing with the plot's data to 1e-10; the truth recomputed independently and agreeing to 1e-10; the PNG's md5 equals the recorded one; the recorded `checks` recomputed. Close with "Covariate-effect figure, truth, independent medians and exact displayed code verified."

- [ ] **Step 4: Run them.** From the worktree root:

```sh
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-PORT-PLAN/run-step.zsh p31-covariate-effect-export Rscript dev/simstudy/vignette-lesson/covariate-effect-export.R /Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-PORT-PLAN/run-step.zsh p32-covariate-effect-verify Rscript dev/simstudy/vignette-lesson/covariate-effect-verify.R /Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
```

Look at the PNG with the Read tool and describe it in the report.

- [ ] **Step 5: The lesson.** In Lesson 3, after the `remaining-plots-gradient-2` image and its paragraph (in "What does an effect mean for a species' distribution?" or wherever the occupancy-gradient curves now sit), add a short subsection, for example `### The same curves in original units`, containing the snippet's chunks verbatim and prose that:
  - says what differs from `plotOccupancyGradient()`: the same conditions (other numeric covariates at their medians, hidden factors and spatial terms at zero) but the horizontal axis in the covariate's original units and running over its full observed range, which is what a report or a reader outside modelling needs; and that the function returns a named list with one plot per covariate, so take the element by name;
  - restores the old file's point (`ORIG_occJSDM.Rmd` lines 362 to 371): the coefficient plot says whether an effect is credibly nonzero on the log-odds scale, this curve shows how large it is as a probability, and the same coefficient moves occupancy a lot or a little depending on the species' baseline;
  - reads the figure with the bundle's `checks` computed inline (for how many species the generating curve lies inside the band everywhere, and the largest departure), and says the departures are the fit's, as the gradient curves above showed, not the function's.
  Then update the two places that name the function without showing it: the sentence near the patchwork paragraph ("`plotCovariateEffect()` returns a named list ...", keep it accurate and point to the new subsection) and the "Find the function for your question" entry's truth check. In "Reproduce the teaching figures", add the new exporter and verifier commands beside the existing ones, in the same form.

- [ ] **Step 6: Render and check.** Render Lesson 3 with both commands; apply the PNG policy; rerun `covariate-effect-verify.R` through `run-step.zsh` as `p36-covariate-effect-verify` and confirm the "Lesson 3 displays ..." line; run `remaining-plots-verify.R` through `run-step.zsh` as `p37-remaining-plots-verify` to confirm it still passes; run the Global Constraints checks with `<files>` = `vignettes/occJSDM-lesson-3.Rmd`.

- [ ] **Step 7: Commit.**

```bash
git add dev/simstudy/vignette-lesson/covariate-effect-examples.Rmd dev/simstudy/vignette-lesson/covariate-effect-export.R dev/simstudy/vignette-lesson/covariate-effect-verify.R vignettes/teaching-data/covariate-effect-data.rds vignettes/teaching-data/covariate-effect-1.png vignettes/occJSDM-lesson-3.Rmd vignettes/occJSDM-lesson-3.md dev/simstudy/vignette-lesson/README.md
git commit -m "Lesson 3: show plotCovariateEffect() in original units

New exporter and verifier pair; no existing bundle changes. Ported from
ORIG_occJSDM.Rmd.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

(Drop `README.md` from `git add` if it needed no change.)

---

### Task 4: Three short additions in Lessons 2, 3 and 4

**Files:**
- Modify: `vignettes/occJSDM-lesson-2.Rmd`, `-3.Rmd`, `-4.Rmd` and their rendered `.md`.

- [ ] **Step 1: Lesson 2, threshold advice (keeper 3).** One sentence, where Lesson 2 explains `threshold = 1` (the paragraph beginning "`collCovariates = \"X_theta\"` names the collection covariate", or the `threshold` bullet in "The main settings are"; pick one and do not repeat it in the other). Content, from `ORIG_occJSDM.Rmd` lines 229 to 233 and the help page (`R/runOccJSDM.R` lines 295 to 299): keep the threshold at one read unless there is a specific reason, because the model already separates false positives through its false-positive rates, and raising the threshold discards weak true detections; the threshold must be at least one, and `threshold = 0` (modelling read counts directly) is not supported and stops with an error. Check that Lesson 2 does not already say this (lines 574 to 576 discuss threshold one); if it does, add only what is missing.

- [ ] **Step 2: Lesson 3, reading low laboratory detection rates (minor gap 2).** Beside "Laboratory true-positive and false-positive rates by primer", after the two primer images, two or three sentences from `ORIG_occJSDM.Rmd` lines 529 to 533: these rates are estimated for each primer and species separately, so they compare primers species by species; a species with low detection under one primer but not the other suggests a primer mismatch, and low detection under both suggests a species that leaves little DNA or a low collection rate upstream; and what to do on one's own survey (for example, keep both primers, or treat a species detected weakly by every primer with caution). Read the two primer figures (Read tool on `vignettes/teaching-data/native-plot-primer-1.png` and `-2.png`) and, if they show such a contrast, point to it; do not claim one that is not there.

- [ ] **Step 3: Lesson 4, `useEnvCov` (minor gap 4).** In the paragraph after the `prediction-native-call` chunk (the chunk is frozen; the paragraph is not), after the `useBiotic` sentence, add half a sentence that `useEnvCov` works the same way, so the call above includes the environmental covariates by default, and `useEnvCov = FALSE` predicts from the hidden factors (and spatial field) alone, a diagnostic rather than a prediction. Check the defaults against `R/output.R` lines 1555 to 1565 and the `resolveTerm()` lines below them.

- [ ] **Step 4: Render and check.** Render Lessons 2, 3 and 4 with both commands, apply the PNG policy (Lesson 2's `probability-and-state-1.png` and Lesson 3's two noise figures), run the Global Constraints checks with the three Rmds as `<files>`, and run `prediction-verify.R` through `run-step.zsh` as `p44-prediction-verify` with arguments `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples /Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples/prediction`.

- [ ] **Step 5: Commit.**

```bash
git add vignettes/occJSDM-lesson-2.Rmd vignettes/occJSDM-lesson-2.md vignettes/occJSDM-lesson-3.Rmd vignettes/occJSDM-lesson-3.md vignettes/occJSDM-lesson-4.Rmd vignettes/occJSDM-lesson-4.md
git commit -m "Port threshold advice, primer-rate reading and useEnvCov to Lessons 2 to 4

Ported from ORIG_occJSDM.Rmd.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Delete ORIG, update the living documents, verify everything, draft the pull request

**Files:**
- Delete: `vignettes/ORIG_occJSDM.Rmd`.
- Modify: `.Rbuildignore` (drop `^vignettes/ORIG_occJSDM\.Rmd$`), `vignettes/LESSON-PLAN.md`, `TODO.md`, `AGENTS.md`.
- Create (not committed): `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-PORT-PLAN/pr-body.md`.

- [ ] **Step 1: Delete and update.**
  - `git rm vignettes/ORIG_occJSDM.Rmd`; remove its `.Rbuildignore` line.
  - `vignettes/LESSON-PLAN.md`: in "Other teaching work still needed", replace the "Simulator reference" bullet with one saying its content is now in Lesson 0 and the simulator's help page, and that the file is deleted when Lesson 0 is published (`TODO.md`). In "Migration checklist for the original function walkthrough", replace the sentence that links `ORIG_occJSDM.Rmd` with one saying the original walkthrough was deleted on today's date after its remaining material was ported (git history before this pull request holds it). Append a dated decision-log entry (date from `date '+%-d %B %Y'`): "**\<date\>, port of the old vignettes:** pull request (link added when opened). Everything useful left in `ORIG_occJSDM.Rmd` and `simulateOccJSDMData.Rmd` moved to the lessons and the simulator's help page: Lesson 0 teaches the trait formula and, in its appendix, how to simulate a spatial field; Lesson 2 gives the threshold advice; Lesson 3 shows `plotCovariateEffect()` and how to read low primer detection rates; Lesson 4 names `useEnvCov`. `ORIG_occJSDM.Rmd` is deleted; `simulateOccJSDMData.Rmd` ships through the beta and is deleted when Lesson 0 is published." Delete any clause the tasks did not carry out.
  - `TODO.md`: append to the "Publish the teaching lessons after Doug's review" item: "When Lesson 0 is published, delete `vignettes/simulateOccJSDMData.Rmd` (its content is in Lesson 0 and `?simulateOccJSDMData` since \<date\>) and point the quickstart's "Where to go next" and the README's Getting started section to Lesson 0 instead of `vignette(\"simulateOccJSDMData\")`."
  - `AGENTS.md`: if its file list describes `ORIG_occJSDM.Rmd`, remove that; leave historical session notes unchanged.
  Commit with message "Delete ORIG_occJSDM.Rmd now its content is in the lessons" and the trailer.

- [ ] **Step 2: Every verifier.** Through `run-step.zsh` with step names `p90` to `p101`, the same list as Phase C's Task 9 Step 2 (`LESSON-STYLE-PLAN-C.md` lines 1008 to 1027), plus `covariate-effect-verify.R`, plus `Rscript -e 'devtools::test()'` (all package tests; report the summary), with `check-frozen-chunks.R` from this plan's directory. Expected lines as listed there, with "All 36 frozen chunks match main at fbe3ed6."

- [ ] **Step 3: Counts and scope.**

```sh
git diff --stat fbe3ed6 -- src DESCRIPTION NAMESPACE dev/simstudy/spatial-design-sweep dev/simstudy/jsdm-package-comparison
git diff --stat fbe3ed6 -- R man
git diff --stat fbe3ed6 -- vignettes/teaching-data
git log fbe3ed6..HEAD --format=%B | /usr/bin/grep -c "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git log --oneline fbe3ed6..HEAD | wc -l
```

Expected: the first is empty; the second shows only `R/simulateData.R` and `man/simulateOccJSDMData.Rd`; the third only the two new `covariate-effect-*` files; the trailer count equals the commit count.

- [ ] **Step 4: HTML for Doug's read.** Render Lessons 0, 2, 3 and 4 with the HTML command, apply the PNG policy, confirm `git status --short` is clean.

- [ ] **Step 5: Draft the pull request body** at the path above (markdown rules apply): Summary (one paragraph: what was ported and where, `ORIG_occJSDM.Rmd` deleted, `simulateOccJSDMData.Rmd` kept through the beta and why); a list mapping each audited item to its new home (keepers 1, 3 and 4; minor gaps 1 to 4 and 6) and the items deliberately not ported (variation partitioning, `plotBiplot()` arguments, and the obsolete material: `gt` in `listParams`, the planned count model, the old `returnOccupancyRates()` array shape, old WAIC numbers, the old latent-presence walkthrough and image, dead links, the array-location table, the residual-correlation clustering aside); the new exporter pair; Needs Alex: whether to drop `sigma_bs` and `sigma_ts` from the simulator or document them as unused, and confirmation of the help page's statements about the spatial field and the coefficient formula; checks run (Steps 2 and 3 output lines); the four HTML paths; the last line `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.

- [ ] **Step 6: Hand back.** Report the commit list, the check outputs and the PR body path.

Controller steps: push `codex/lesson-port`; `gh pr create --base main --title "Port the old vignettes into the lessons and the simulator help page" --body-file /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-PORT-PLAN/pr-body.md`; fill the decision-log link in a follow-up commit with the trailer.

---

## Self-review notes

- **Decision coverage.** Keeper 1: Task 1 Step 4 and Task 2 Step 2. Keeper 2: omitted by decision, named as not ported in Task 5 Step 5. Keeper 3: Task 4 Step 1. Keeper 4 with minor gap 3: Task 1 Step 3 and Task 2 Step 2. Minor gap 1: Task 3. Minor gap 2: Task 4 Step 2. Minor gap 4: Task 4 Step 3. Minor gap 5: skipped by decision, named in Task 5 Step 5. Minor gap 6: Task 1 Step 2, Task 2 Step 2 and the needs-Alex list. Partly covered items from the audit: coordinates drawn either way (Task 1 Step 2), the printed model-inference message (already covered; nothing to do). Deletion: Task 5 Step 1; the simulator guide's deletion is a `TODO.md` step by decision.
- **Placeholders.** None: every new text's content is specified by its source lines and the facts checked while planning; exact wording is the implementer's, under the writing rules.
- **Names.** Chunk prefix `covariate-effect-`, bundle `covariate-effect-data.rds`, figure `covariate-effect-1.png`, object `covariate_effect_1`, used consistently in Task 3 and the scope check of Task 5.

## Changes after planning (4 October 2026)

- **The help page moved to its own pull request.** Task 2's roxygen edit to `R/simulateData.R` made every fit-reading lesson verifier fail. Their source fingerprint covers every file in `R/`. Doug chose a separate pull request (`codex/simulator-help`), independent of PR #14, so this branch restores `R/` and `man/` to fbe3ed6 and keeps only the corrections to `simulateOccJSDMData.Rmd`.
- **Task 4b, added by Doug.** Lesson 3 gained an invented species-by-primer interaction figure, drawn by the real `plotDetectionRates()` from invented draws. Advice on keeping several primers or choosing one follows it.
- **Task 6, added by Doug.** Lesson 3's coloured figures use the colour-blind-safe Okabe-Ito palette. With Doug's approval this changed four chunks in two frozen snippets (scale lines only) and re-exported `native-plots.rds`, `remaining-plots-data.rds` and `covariate-effect-data.rds` and four exporter PNGs; no exporter or verifier changed.
- **A separate colours pull request.** The package's own plot defaults move to the same Okabe-Ito order in `codex/plot-colours`, so Lesson 3's figures keep their colours when it merges.

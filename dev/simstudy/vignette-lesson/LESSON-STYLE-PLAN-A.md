# Lesson Rewrite Phase A Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Rewrite Lessons 1, 0 and 2 in the teaching style, and write the new intuition lesson, so that each teaches rather than describes, with Doug reviewing one pull request per lesson.

**Architecture:** Each lesson is an R Markdown file under `vignettes/` that renders from a committed teaching bundle without fitting anything. A pass edits prose, moves a few sections, adds small run chunks that show data or read output, converts pipe tables to bullets, and gathers technical material into an appendix; it does not change bundles, verifiers or package code. The intuition lesson is a new file that computes everything by hand from the Lesson 0 bundle and small inline data. Every lesson is rendered twice (HTML for checking, unwrapped markdown for the site), verified by the existing scripts, and opened as a pull request.

**Tech Stack:** R Markdown (rmarkdown, knitr), dplyr, tidyr, tibble, ggplot2, the committed bundles in `vignettes/teaching-data/`, the verifiers in `dev/simstudy/vignette-lesson/` and `dev/simstudy/spatial-design-sweep/`, Node for the site test, `gh` for pull requests.

**Spec:** `dev/simstudy/vignette-lesson/LESSON-STYLE-SPEC.md` (committed 516d072 and later). The per-lesson diagnostic reports it argues from are in `dev/simstudy/vignette-lesson/style-diagnostic/`. The writing criteria are in `AGENTS.md`, section "Writing the vignettes and lessons".

## Global Constraints

- Work on a branch in its own worktree under `~/src/occJSDM-worktrees/`, never in `~/src/occJSDM` (Doug's checkout) and never on main. Branch names: `codex/lesson-1-teaching`, `codex/lesson-intuition`, `codex/lesson-0-teaching`, `codex/lesson-2-teaching`. Each task branches from main after the previous task's pull request has merged; if Doug's review is pending, branch from the previous task's branch and say so in the pull request.
- Markdown rules for every file touched: one line per paragraph, no hard wrapping, no em-dashes anywhere (use `--`), no pipe tables (bullet lists instead), escape a literal `>` `<` `~` `|` in prose as `\>` `\<` `\~` `\|`, never inside code spans or fenced blocks. Re-read a file immediately before and after editing it.
- Commit messages end with `Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>`. Pull request bodies end with `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.
- Two commits per lesson pass, in this order: first the teaching commit (additions, moves, new chunks, appendix, table conversions, factual corrections), then the copyedit commit (sentence splitting, "and/or", stops, repeated asides, wording), so Doug sees both passes separately. The intuition lesson has one writing commit and one copyedit commit.
- No changes to package code (`R/`, `src/`, `DESCRIPTION`, `NAMESPACE`), to any teaching bundle under `vignettes/teaching-data/`, to any verifier or exporter, or to the frozen protocol `dev/simstudy/spatial-design-sweep/PLAN.md`.
- Hard-coded result numbers in prose become inline R from the bundle wherever the bundle holds the number; a number the bundle lacks is stated with its source.
- The intended reader: an ecologist who runs eDNA or presence/absence surveys, knows basic R and the pipe, has not fitted an occupancy model or a JSDM, and needs reminding of the ideas rather than introducing. A quantitative reader finds the evidence in the appendix.
- R glosses: cut explanations of the pipe, `select()`, `filter()`, `mutate()`, `summarise()`, `group_by()`, `readRDS()`, lists and `$`; keep a gloss only where the operation is unusual (a join on species identity, a pivot, an array slice) and the ecology reads through the code.
- Cross-references keep the current lesson numbers; the Phase B renumbering changes them all at once. Refer to a lesson by what it teaches where a number would otherwise be the only identifier.
- Render commands, run from the worktree root (replace the stem):

```sh
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-1.Rmd", output_format="rmarkdown::html_vignette")'
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-1.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
```

- PNG policy: a render rewrites every figure file. Keep a changed PNG only when its content changed (a moved figure keeps its content; a new chunk adds a new file); revert the others with `git checkout -- vignettes/<stem>_files/` for the untouched ones, so the diff shows only real figure changes. Check with `git status --short vignettes/` before committing.
- Checks before every pull request, all from the worktree root, all must pass:

```sh
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
node --test dev/simstudy/lesson-site/test_lessons.js
git grep -n "—" -- vignettes ':!*.png' ':!*.rds' | wc -l
git grep -n "^|" -- 'vignettes/*.Rmd' | wc -l
```

  The two counts must be 0. The HTML output is gitignored and is not committed; the `.md` render and any new or changed PNGs are.

- Pull request body: what the pass changed (teaching additions by section, moves, new chunks, appendix, tables, factual corrections), what it deliberately left (deferred items from the report), a "needs Alex" list of any package fact the lesson now states, the checks run with their output lines, and a note that Doug reviews line by line.

---

## File structure

- `vignettes/occJSDM-lesson-1.Rmd` and `.md`: Lesson 1 pass (Tasks 1 and 2).
- `vignettes/occJSDM-lesson-intuition.Rmd` and `.md` (provisional stem; Phase B renames it): the new lesson (Task 3), plus registration in `vignettes/lesson-links.R`, `dev/simstudy/vignette-lesson/test_lesson_links.R`, `_config.yml`, `assets/js/lessons.js` and `.Rbuildignore`.
- `vignettes/occJSDM-lesson-0.Rmd` and `.md`: Lesson 0 pass (Task 4).
- `vignettes/occJSDM-lesson-2.Rmd` and `.md`: Lesson 2 pass (Task 5).
- `vignettes/LESSON-PLAN.md`: one decisions-log line per merged pass (each task).
- Figures under `vignettes/<stem>_files/figure-gfm/` and, for Lesson 2, `vignettes/teaching-data/lesson-2-*.png`.

## Shared text: the prerequisites block

Every lesson gets this block, as its own paragraph plus a bullet list, directly after the setup chunk and before the first teaching section, with the verb list trimmed to what that lesson actually uses (check with `git grep -o -h -E "\b(select|filter|mutate|summarise|group_by|arrange|left_join|inner_join|anti_join|pivot_longer|pivot_wider|across|case_when|slice_head|count)\(" vignettes/<stem>.Rmd | sort | uniq -c`):

```markdown
**What this lesson assumes you know.** The code uses base R and the tidyverse: the pipe `|>`, and from dplyr and tidyr the verbs listed below. If any are new, the two chapters of R for Data Science on [data transformation](https://r4ds.hadley.nz/data-transform) and [data tidying](https://r4ds.hadley.nz/data-tidy) teach everything used here in an afternoon. Operations that are unusual, such as joining two tables on species identity, are explained where they appear.

- `select()`, `filter()`, `mutate()`, `summarise()` and `group_by()` to choose, keep, add, summarise and group rows and columns.
- `left_join()` and `inner_join()` to combine tables on shared identifiers.
- `pivot_longer()` and `pivot_wider()` to move between one row per measurement and one column per measurement.
```

---

### Task 1: Lesson 1, teaching commit

**Files:**
- Modify: `vignettes/occJSDM-lesson-1.Rmd` (964 lines)
- Read: `dev/simstudy/vignette-lesson/style-diagnostic/lesson-1.md`, the Lesson 1 section of the spec, `AGENTS.md` lines 113 to 125

**Interfaces:**
- Produces: the section order, the appendix heading `## Appendix: evidence and reproduction` and the vocabulary (two-stage process, false negative, false positive, collection, detection, naive occupancy, credible interval, chains, Rhat, effective sample size) that Task 3 introduces and Task 4 leads into.

- [ ] **Step 1: Create the worktree and branch**

```bash
cd ~/src/occJSDM && git pull --ff-only origin main
git worktree add ~/src/occJSDM-worktrees/lesson-1-teaching -b codex/lesson-1-teaching main
cd ~/src/occJSDM-worktrees/lesson-1-teaching
```

- [ ] **Step 2: Read the report and the spec section, then the lesson end to end**

Read `style-diagnostic/lesson-1.md` in full, the Lesson 1 bullets of `LESSON-STYLE-SPEC.md`, and then `vignettes/occJSDM-lesson-1.Rmd` once from top to bottom as the reader would, before editing anything. Keep the report open while editing: every bullet in its section 2 with a line number is an edit to make or to record as deferred in the pull request.

- [ ] **Step 3: Record the before state**

Render the markdown once before editing and keep the file for comparison:

```sh
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-1.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
cp vignettes/occJSDM-lesson-1.md /tmp/lesson-1-before.md
git checkout -- vignettes/
```

- [ ] **Step 4: Add the prerequisites block**

Insert the shared prerequisites block after the paragraph that begins "We use saved results so that reading or knitting the lesson" (line 31 region), with the verb list trimmed to Lesson 1's verbs. Remove the sentence "`|>` passes a result into the next function. We explain the table operations as they appear." from the paragraph above it, since the block replaces it.

- [ ] **Step 5: Show the data the reader will handle**

Directly after the `load-teaching-data` chunk and its explanatory paragraph (line 48 region), add this chunk and one sentence introducing it ("These are the three tables the model receives, and the one that the lesson keeps aside as truth."):

```{r show-survey-data}
str(survey_data, max.level = 1)
head(survey_data$info)
head(survey_data$OTU[, 1:5])
survey_data$traits
names(known_truth)
```

Follow it with a paragraph that reads the output in the reader's terms: `info` has one row per PCR reaction with its site, sample and primer and the covariates (two occupancy covariates, two coordinates that this lesson does not use, one collection covariate); `OTU` has that reaction's read count for each species; `traits` has one row per species; and `known_truth` holds what the simulation knew and the model never sees, named here so the reader can recognise it later.

- [ ] **Step 6: Move convergence before the results**

Move the whole section `## Are the calculations stable enough to interpret?` (lines 765 to 800) to sit immediately after `## Now give occJSDM only the PCR observations` and its chunks, and before `## Calculate the errors ourselves`. Add one sentence at its start saying why it comes first: "Before reading any estimate, check that the chains agree, as the quickstart does; an estimate from chains that disagree is not an estimate." Update the one forward reference the original position had (search for "stable enough" and "convergence" in the prose and fix each sentence so it still reads correctly in the new order).

- [ ] **Step 7: Place the fitting reference and the recovery figure**

Keep `### A short fitting reference for your own data` (lines 202 to 239) as its own `##` section titled `## Fitting your own data: what the call needs`, placed after the convergence section moved in Step 6. Move the `occupancy-recovery` chunk (line 241 region) and the two paragraphs that read it out of that section into `## Calculate the errors ourselves`, as its opening figure. Convert the two pipe tables in the fitting reference (the data-structure table at line 206 and the settings table at line 216) to bullet lists, one bullet per row, in the form "`setting`: what you supply or choose". Add one sentence after the data-structure bullets saying that this lesson's survey has repeated sites and repeated sample identifiers, which is why the two-stage model was inferred, and that the message `runOccJSDM()` prints says which model it chose.

- [ ] **Step 8: Move the all-positive-samples table to the reveal section**

Move the paragraph beginning "To put the four examples in perspective" and the `all-positive-samples` chunk and its following two paragraphs (lines 700 to 761 region, after the alternative-priors figure) into `## Reveal the truth and compare it with the fit`, after the four cases have been revealed. Then write the reading the table lacks, as two paragraphs with these claims and no others, computing each number inline from `positive_sample_summary`:
  - Laboratory false positives are caught by PCR replication within the sample: `r` the count of samples in that category and their mean fitted DNA probability.
  - Field-stage false positives are a limit of field replication, not the model believing contamination: such a sample genuinely contains the species' DNA, so at the sample level it is indistinguishable from a true positive and the model rightly gives it a high DNA-presence probability; the only evidence against site presence is the other field sample at that site, and with two samples per site one positive and one negative is what a truly occupied site often produces, so the mean fitted site probability (`r` from the table) is the model weighing a contamination prior near 5% against a collection probability of 35 to 75%. State that the simulated contamination rates (field-stage 2 to 8%, laboratory 2.5 to 6.5%, from Lesson 0's settings) sit inside the default priors' assumption, so this is the default-prior fit on a good-practice survey.
  - The practical advice: more field samples per site, not more PCRs, guards site occupancy against contamination at collection; PCR replication guards against laboratory false positives.
  In the alternative-priors section that the table left, add one sentence making its separate point: loosening the contamination priors does not resolve the field-stage case, because the ambiguity is in the data, not the prior.

- [ ] **Step 9: State the pull towards the middle where its output appears**

After the `occupancy-recovery` figure (now in the errors section) and again after the band table (line 308 region), add the reading: estimates are pulled towards the middle, with the low band's true mean and estimated mean given inline from `error_summary` (the report's rendered values are 6.5% true against 20.7% estimated for the default fit; compute, do not type), and say in one sentence why a model with uncertainty does this (a site with few detections cannot be told apart from a site with low probability, so the posterior hedges).

- [ ] **Step 10: Work through the report's remaining teaching gaps**

For each bullet in section 2 of `style-diagnostic/lesson-1.md` not covered by Steps 4 to 9, make the edit it describes at the line it names, in the reader's terms, following the criterion it cites. In particular:
  - Line 106 region: gloss "collection failure" and "field-stage contamination" in a clause each where they first appear, and name the two-stage process once in full (field collection can miss a present species or pick up DNA from elsewhere; PCR can miss DNA in the sample or report DNA that is not there), since Lesson 0 and the intuition lesson will point here.
  - The `## Three questions, three different truths` table (line 102): convert to three bullets, each "question: model quantity, and what the simulation knows".
  - Line 184 region (`n_factors`, `n_lattrait`): say plainly there is no rule yet; start small and check that conclusions do not change when the count changes.
  - Line 95 region: explain the design (two samples, two primers, six PCRs) in one sentence as a choice that shows every stage with modest cost, and that the quickstart's example differs because the structure, not the counts, is what the model reads.
  - Lines 466 to 478: collapse the threshold-conversion explanation to one sentence and a link to Lesson 0's section on the simulator's read counts.
  - Lines 352 to 420: keep the maps, shorten their prose to what the figures show, and open with the question they answer.
  - Lines 640 to 653: in the case table, show "not applicable" for the true collection probability in the field-stage case, with the reason in the caption.
  - Line 432 region: do not add contamination-prior advice; list "a basis for changing the contamination priors, such as rates seen in blanks and negative controls" under "needs Alex" in the pull request body instead.
  - Line 964 region: replace the closing with three or four takeaways (estimates pulled towards the middle; laboratory false positives are well handled; field-stage false positives are a limit of field replication; field replicates matter for site occupancy), one line each, and the hook "Real surveys often pre-filter detections before modelling; a later addition contrasts that with modelling detection." Remove the "Paper2Agent" and "beta-release error target" disclaimers entirely (Doug's decision).

- [ ] **Step 11: Shorten the unequal-replication section and build the appendix**

Reduce `## Fit a survey with unequal replication` (lines 801 to 933) to about 40 lines: what was removed, why a lost sample is absent rows, what the fit needs, that it runs, and one sentence that unequal numbers of field samples per site, primers per sample and PCRs per sample and primer are all accepted by the fitter (it counts the rows in each block; the package's collection-alignment test fits such a survey). Move its full truth comparison and its second diagnostic system into a new final section `## Appendix: evidence and reproduction`, together with `## Reproduce the lesson and inspect its evidence` (lines 935 to 962) minus the removed disclaimers. The appendix opens with one sentence saying who it is for and that the lesson's conclusions do not depend on reading it. The closing section of the lesson (Step 10's takeaways) ends with one sentence saying what the appendix holds.

- [ ] **Step 12: Correct the facts the report doubted**

Line 33 region: the lesson has six not-run chunks (`reproduce-perfect-fit`, `reproduce-two-stage-fit`, `inspect-fitted-design`, `alternative-priors`, `unbalanced-fit-optional`, `extract-your-fit`), four of them fits; say so. Line 239 region: recompute the memory figure from the dimensions stated in the sentence and show the arithmetic in the sentence; if the sentence's dimensions are not in the bundle, state them as the recorded fit's settings.

- [ ] **Step 13: Convert the remaining pipe tables and cut the routine glosses**

Convert every remaining pipe table in the Rmd (the report lists lines 102, 206, 216 and 426; check with `git grep -n "^|" -- vignettes/occJSDM-lesson-1.Rmd`) to bullet lists. Cut the routine tidyverse explanations at lines 56, 260, 306, 356 and 488 and any other gloss of a verb in the prerequisites list; keep the explanations of joins on species identity and of the paired-row removal.

- [ ] **Step 14: Render, inspect, and commit the teaching pass**

```sh
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-1.Rmd", output_format="rmarkdown::html_vignette")'
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-1.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
git status --short vignettes/
```

Both renders must succeed with no new warnings. Read the rendered `.md` end to end once. Apply the PNG policy. Then:

```bash
git add vignettes/occJSDM-lesson-1.Rmd vignettes/occJSDM-lesson-1.md vignettes/occJSDM-lesson-1_files
git commit -m "Lesson 1 teaching pass: show the data, read every output, convergence first, appendix

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

---

### Task 2: Lesson 1, copyedit commit and pull request

**Files:**
- Modify: `vignettes/occJSDM-lesson-1.Rmd`, `vignettes/occJSDM-lesson-1.md`, `vignettes/LESSON-PLAN.md`

**Interfaces:**
- Consumes: Task 1's commit on `codex/lesson-1-teaching`.
- Produces: the merged Lesson 1 that Task 3's vocabulary must match.

- [ ] **Step 1: Copyedit the whole lesson, keeping every teaching addition**

Read the Rmd from top to bottom and, sentence by sentence: split sentences over about 30 words; replace every "and/or"; fix missing stops; remove any aside that repeats its sentence; check each factual claim in added prose against the rendered output and correct it; make sure every number that the bundle holds is computed inline. Do not remove or shorten a teaching addition from Task 1; if one is wrong, fix it. Section 5 of the report lists the sentences it found.

- [ ] **Step 2: Render and run every check**

```sh
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-1.Rmd", output_format="rmarkdown::html_vignette")'
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-1.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
Rscript dev/simstudy/vignette-lesson/test_lesson.R
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
node --test dev/simstudy/lesson-site/test_lessons.js
git grep -n "—" -- vignettes ':!*.png' ':!*.rds' | wc -l
git grep -n "^|" -- 'vignettes/*.Rmd' | wc -l
```

Expected: both renders clean; `test_lesson.R` passes (it checks the lesson helpers, not the prose); the link test prints "Shared lesson-link regression checks passed."; the site test reports 0 failures; both counts 0 for Lesson 1 (other lessons still have tables until their passes). The bundle verifier `verify_lesson.R` needs the raw build directory, which is not on this machine; the bundle is unchanged, so it is not run, and the pull request says so.

- [ ] **Step 3: Record the pass in the lesson plan and commit**

Append one line to the decisions log in `vignettes/LESSON-PLAN.md`: the date, "Lesson 1 teaching pass", the pull request, and one sentence on what changed. Then:

```bash
git add vignettes/occJSDM-lesson-1.Rmd vignettes/occJSDM-lesson-1.md vignettes/occJSDM-lesson-1_files vignettes/LESSON-PLAN.md
git commit -m "Lesson 1 copyedit pass

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
git push -u origin codex/lesson-1-teaching
```

- [ ] **Step 4: Open the pull request**

`gh pr create --base main --title "Lesson 1 in the teaching style"` with a body following the Global Constraints template: the changes by section, the deferred items from the report with the reason, the "needs Alex" list (at least the contamination-prior basis), the check output lines, and the closing attribution line. Report the URL.

---

### Task 3: The intuition lesson

**Files:**
- Create: `vignettes/occJSDM-lesson-intuition.Rmd` and its `.md` render and figures
- Modify: `vignettes/lesson-links.R` (add the stem to `lessons`), `dev/simstudy/vignette-lesson/test_lesson_links.R` (add the stem to its list), `_config.yml` (add `vignettes/occJSDM-lesson-intuition.md` and its figure directory to `exclude`, in the same form as the Lesson 0 entries), `assets/js/lessons.js` (add an entry with `published: false`, exactly as the Lesson 0 entry has, placed between the quickstart and Lesson 0; the Node test checks that `published` agrees with the `_config.yml` exclusions, figures included), `.Rbuildignore` (add `^vignettes/occJSDM-lesson-intuition\.(Rmd|md)$` and `^vignettes/occJSDM-lesson-intuition_files$`), `vignettes/LESSON-PLAN.md`
- Read: the spec section "The intuition lesson", the merged Lesson 1, Lesson 4 lines 219 to 285, and Doug's draft at `https://github.com/dougwyu/s-jSDM/blob/v0.2.1/Code/sjSDM_mojo_tutorial.Rmd` lines 40 to 180 (fetch with `gh api "repos/dougwyu/s-jSDM/contents/Code/sjSDM_mojo_tutorial.Rmd?ref=v0.2.1" --jq .content | base64 -d`)

**Interfaces:**
- Consumes: `vignettes/teaching-data/nonspatial-lesson.rds` elements `samples` (columns arm, species, Site, Sample, z, w, collection_probability, occupancy_truth), `input$sim$data_list$info` (Site, Sample, Primer, X_psi.EnvCov.1, X_psi.EnvCov.2, X_theta per PCR row) and `input$sim$true_params` (`z_true` 100 x 10, `w_true` 200 x 10).
- Produces: the lesson's section anchors that Phase B will link from Lessons 1, 3 and 4: `#what-an-occupancy-model-does`, `#the-sites-where-the-species-was-missed`, `#what-a-joint-model-does`, `#what-a-joint-model-adds-to-a-factorisation`.

- [ ] **Step 1: Worktree, branch and registration**

```bash
cd ~/src/occJSDM && git pull --ff-only origin main
git worktree add ~/src/occJSDM-worktrees/lesson-intuition -b codex/lesson-intuition main
cd ~/src/occJSDM-worktrees/lesson-intuition
```

Make the five registration edits listed under Files. Run `node --test dev/simstudy/lesson-site/test_lessons.js` and fix the entry until the test's withheld-lessons check passes (it compares `_config.yml` with the unpublished markers in `lessons.js`).

- [ ] **Step 2: Create the file skeleton**

```markdown
---
title: "What occupancy models and joint species distribution models do"
output:
  rmarkdown::html_vignette:
    toc: true
    css: teaching.css
vignette: >
  %\VignetteIndexEntry{What occupancy models and joint species distribution models do}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, collapse = FALSE, comment = "#>", fig.width = 8, fig.height = 5, dpi = 150, dev = "png")
source("lesson-links.R")
```

## Where this lesson fits

[one paragraph: this lesson explains the two ideas the rest of the lessons use, by hand and without fitting a model; read it before Lesson 1 if either idea is new or rusty; Lesson 0 (optional) explains where the survey comes from]

[the shared prerequisites block, verbs trimmed to this lesson]

```{r load, message=FALSE}
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

lesson <- readRDS("teaching-data/nonspatial-lesson.rds")
truth_samples <- as_tibble(lesson$samples) |> filter(arm == "default")
survey_info <- lesson$input$sim$data_list$info
```

## What an occupancy model does
## The sites where the species was missed
## Covariates say which sites those are
## How many replicates are enough
## What a joint model does
## What a joint model adds to a factorisation
## What you can now recognise in the other lessons
```

- [ ] **Step 3: Write the occupancy half, with these chunks**

The teaching sequence and the code for each section. The prose around each chunk follows the spec's "Occupancy half" bullet and the criteria; every printed object is read in the next paragraph.

`## What an occupancy model does`: say what a single-season occupancy model assumes (a site is occupied or not; each replicate sample of an occupied site detects the species with some probability; an unoccupied site never produces a detection, which is the no-false-positive assumption this half makes and Lesson 1 drops). Build the one-stage table from the simulation's true sample states so that assumption holds by construction, and say so:

```{r one-stage-table}
one_stage <- truth_samples |>
  filter(species == "OTU_5") |>
  mutate(detected = as.integer(z == 1 & w == 1)) |>
  select(Site, Sample, detected, z, collection_probability, occupancy_truth)

head(one_stage, 6)

site_history <- one_stage |>
  group_by(Site) |>
  summarise(
    samples = n(),
    detections = sum(detected),
    occupied = first(z),
    .groups = "drop"
  )

count(site_history, detections, occupied)
```

Read the count table: how many sites had two, one or no detections, and (because we are allowed to peek) how many of the no-detection sites were occupied. Then naive occupancy and the two-parameter maximum-likelihood fit:

```{r one-stage-mle}
naive_occupancy <- mean(site_history$detections > 0)

one_stage_loglik <- function(par, history) {
  psi <- plogis(par[1])
  p <- plogis(par[2])
  with(history, {
    detected_sites <- detections > 0
    ll_detected <- log(psi) + detections * log(p) + (samples - detections) * log(1 - p)
    ll_missed <- log(psi * (1 - p)^samples + (1 - psi))
    sum(ifelse(detected_sites, ll_detected, ll_missed))
  })
}

fit <- optim(c(0, 0), one_stage_loglik, history = site_history,
             control = list(fnscale = -1))
estimates <- c(psi = plogis(fit$par[1]), p = plogis(fit$par[2]))
round(estimates, 3)

c(naive = naive_occupancy,
  true_fraction_occupied = mean(site_history$occupied),
  true_mean_detectability = mean(one_stage$collection_probability[one_stage$z == 1]))
```

Read these: naive occupancy undercounts because missed sites are counted as absent; the model's psi is close to the realised fraction occupied; its p is the per-sample detectability, which in this simulation is the collection probability (a sample with DNA is almost always detected by twelve PCRs, so the one-stage p is close to the true mean collection probability for occupied sites; say this is a feature of the simulation, not a law).

`## The sites where the species was missed`:

```{r missed-sites}
M <- 2
psi_hat <- estimates[["psi"]]
p_hat <- estimates[["p"]]

conditional_occupancy <- psi_hat * (1 - p_hat)^M /
  (psi_hat * (1 - p_hat)^M + (1 - psi_hat))

missed <- site_history |> filter(detections == 0)

c(no_detection_sites = nrow(missed),
  conditional_probability = round(conditional_occupancy, 3),
  expected_missed_sites = round(nrow(missed) * conditional_occupancy, 1),
  actually_occupied = sum(missed$occupied))
```

Give the conditional probability its own paragraph: it is the probability that a site with no detections is nonetheless occupied, the ratio of "occupied and missed twice" to "occupied and missed twice, or empty"; multiply by the number of all-negative sites to get how many the survey missed; compare with the number actually occupied, which the simulation lets us see.

`## Covariates say which sites those are`: fit the same model with an occupancy covariate and a collection covariate, by hand:

```{r covariate-mle}
site_covariates <- survey_info |>
  distinct(Site, environment = X_psi.EnvCov.1)
sample_covariates <- survey_info |>
  distinct(Site, Sample, effort = X_theta)

history_cov <- one_stage |>
  left_join(site_covariates, by = "Site") |>
  left_join(sample_covariates, by = c("Site", "Sample")) |>
  mutate(environment = as.numeric(scale(environment)),
         effort = as.numeric(scale(effort)))

covariate_loglik <- function(par, d) {
  psi <- plogis(par[1] + par[2] * d$environment)
  p <- plogis(par[3] + par[4] * d$effort)
  per_sample <- ifelse(d$detected == 1, log(p), log(1 - p))
  by_site <- d |>
    mutate(per_sample = per_sample, psi = psi) |>
    group_by(Site) |>
    summarise(
      psi = first(psi),
      detected_any = any(detected == 1),
      ll_if_occupied = sum(per_sample),
      .groups = "drop"
    )
  with(by_site, sum(ifelse(detected_any,
    log(psi) + ll_if_occupied,
    log(psi * exp(ll_if_occupied) + (1 - psi)))))
}

fit_cov <- optim(c(0, 0, 0, 0), covariate_loglik, d = history_cov,
                 control = list(fnscale = -1), method = "BFGS")
round(setNames(fit_cov$par, c("psi_intercept", "psi_environment", "p_intercept", "p_effort")), 3)

missed_by_site <- history_cov |>
  group_by(Site) |>
  summarise(detected_any = any(detected == 1), environment = first(environment),
            occupied = first(z), .groups = "drop") |>
  filter(!detected_any) |>
  mutate(psi_site = plogis(fit_cov$par[1] + fit_cov$par[2] * environment),
         conditional = psi_site * (1 - plogis(fit_cov$par[3]))^M /
           (psi_site * (1 - plogis(fit_cov$par[3]))^M + (1 - psi_site))) |>
  arrange(desc(conditional))

head(missed_by_site, 8)
```

Read it: the environment slope says where the species is likelier; among the all-negative sites, the ones in suitable habitat now carry a higher conditional probability, and the simulation's `occupied` column shows those are indeed where the misses are; the effort slope says which samples were likelier to miss. Say in one sentence that these are occJSDM's `occCovariates` and `collCovariates`, named psi and theta there.

`## How many replicates are enough`:

```{r replicate-design, fig.height=4}
tibble(M = 1:8) |>
  mutate(chance_of_missing = (1 - p_hat)^M) |>
  ggplot(aes(M, chance_of_missing)) +
  geom_line() + geom_point() +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Field samples per site", y = "Chance of missing a present species",
       title = "With this detectability, how many samples does a site need?")
```

Read the figure with the formula 1 - (1 - p)^M in the prose, and the practical advice: the number of replicates is a design decision that the detectability sets; the same formula returns in Lesson 1 for contamination.

Close the half with a not-run chunk showing that the package fits this one-stage model when `info` has repeated rows per site and no `Sample` column, and that its detection probability is the per-sample one found here:

```{r package-one-stage, eval=FALSE}
one_stage_info <- survey_info |>
  distinct(Site, Sample, X_psi.EnvCov.1, X_psi.EnvCov.2, X_theta) |>
  select(-Sample)
one_stage_otu <- truth_samples |>
  mutate(detected = as.integer(z == 1 & w == 1)) |>
  select(species, Site, Sample, detected) |>
  pivot_wider(names_from = species, values_from = detected) |>
  arrange(Site, Sample) |>
  select(-Site, -Sample) |>
  as.matrix()

fit_one_stage <- runOccJSDM(
  data = list(info = one_stage_info, OTU = one_stage_otu),
  occCovariates = "X_psi.EnvCov.1",
  collCovariates = "X_theta",
  MCMCparams = list(nchain = 2, nburn = 2000, niter = 2000, nthin = 1)
)
# occJSDM prints "occJSDM has inferred occupancy data" for this shape.
```

Say that false positives are excluded here by construction and arrive in Lesson 1 with the second stage.

- [ ] **Step 4: Write the JSDM half, with these chunks**

`## What a joint model does`: open with the stacked-versus-joint contrast in the reader's terms (one model per species cannot say how species co-occur beyond the environment; a joint model adds the residual correlation among species). Then the toy factorisation, taken from Doug's draft with the nouns swapped to sites and species and the numbers kept:

```{r toy-factorisation}
species_scores <- rbind(                 # species x hidden gradients
  species_A = c(wet = 1.0, shade = 0.0),
  species_B = c(0.9, 0.1),
  species_C = c(0.6, 0.6),
  species_D = c(0.1, 0.9),
  species_E = c(0.0, 1.0)
)
site_scores <- rbind(                    # sites x hidden gradients
  site_1 = c(wet = 1.5, shade = 0.0),
  site_2 = c(1.2, 0.2),
  site_3 = c(1.0, 1.0),
  site_4 = c(0.2, 1.3),
  site_5 = c(0.0, 1.5),
  site_6 = c(0.3, 0.3)
)
suitability <- site_scores %*% t(species_scores)
round(suitability, 2)
round(svd(suitability)$d, 3)
```

Read it: thirty numbers, two directions; a site's row is high where its scores line up with a species' scores; the singular values show the table is rank two. Then the rotation check and the species-by-species cross-product from the draft (steps 3 of Doug's tutorial), with the lesson's point: the individual loadings are not identified, the species-by-species matrix is, and that matrix is what the package reports as residual correlations.

```{r toy-rotation}
theta <- pi / 6
rotation <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
max(abs((site_scores %*% rotation) %*% t(species_scores %*% rotation) - suitability))
round(species_scores %*% t(species_scores), 2)
```

Then step 4 of the draft (a measured covariate enters first, the hidden part explains what is left), as prose with the one-line formula, and step 5 (the same picture in ecology) as bullets mapping users, films, ratings, attributes, tastes and the similarity matrix to sites, species, occurrences, environmental covariates, site scores and the association matrix.

`## What a joint model adds to a factorisation`: move the three "additions" paragraphs from Lesson 4 lines 231 to 285 here verbatim apart from the lead-in (it factors the residual, not the data; it is a probability model of presence or absence; it treats the site scores as random and averages over them), followed by the paragraph "Put together, these define the joint probability distribution". Lesson 4 keeps its copy until Phase B trims it; say nothing about that in the lesson. End with one paragraph naming the package's terms for each piece: `n_factors` hidden site factors, species loadings, residual correlations from `returnResidualCorrelationMatrix()`, so Lesson 3's ordination and association sections can point back.

`## What you can now recognise in the other lessons`: five or six bullets, each naming an idea from this lesson and where it reappears (Lesson 1's two stages and its field-replication point; Lesson 1's priors as the good-practice assumption; Lesson 3's credible intervals, hidden factors and residual correlations; Lesson 4's conditional-prediction exercise), with links by current lesson number.

- [ ] **Step 5: Render, check, commit the writing pass**

```sh
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-intuition.Rmd", output_format="rmarkdown::html_vignette")'
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-intuition.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
node --test dev/simstudy/lesson-site/test_lessons.js
git grep -n "—" -- vignettes ':!*.png' ':!*.rds' | wc -l
```

Expected: both renders clean; every chunk ran (no `eval=FALSE` except `package-one-stage`); the printed estimates are sensible (psi within 0.1 of the true fraction occupied, p within 0.1 of the true mean collection probability; if not, the one-stage likelihood has a bug: check the missed-site term); the link test passes with the new stem; the site test passes; 0 em-dashes. Commit:

```bash
git add vignettes/occJSDM-lesson-intuition.Rmd vignettes/occJSDM-lesson-intuition.md vignettes/occJSDM-lesson-intuition_files vignettes/lesson-links.R dev/simstudy/vignette-lesson/test_lesson_links.R _config.yml assets/js/lessons.js .Rbuildignore
git commit -m "Add the intuition lesson on occupancy models and joint species distribution models

Co-Authored-By: Claude Fable 5.1 <noreply@anthropic.com>"
```

- [ ] **Step 6: Copyedit, record, push, pull request**

Copyedit as in Task 2 Step 1; re-render; re-run the checks; append a decisions-log line to `vignettes/LESSON-PLAN.md`; commit "Intuition lesson copyedit pass" with the trailer; push; `gh pr create --base main --title "Add the intuition lesson on occupancy models and JSDMs"` with the standard body, noting the provisional filename and that Phase B numbers it and trims Lesson 4 section 3.

---

### Task 4: Lesson 0 pass

**Files:**
- Modify: `vignettes/occJSDM-lesson-0.Rmd` (463 lines) and its `.md` and figures, `vignettes/LESSON-PLAN.md`
- Read: `style-diagnostic/lesson-0.md`, the spec's Lesson 0 bullets, the merged intuition lesson

**Interfaces:**
- Consumes: the intuition lesson's anchors and Lesson 1's vocabulary.
- Produces: the "No detection" label split (`True negative`, `Missed detection`) in Lesson 0's own tables only; Lesson 1 builds its own labels.

- [ ] **Step 1: Worktree and before-state**

```bash
cd ~/src/occJSDM && git pull --ff-only origin main
git worktree add ~/src/occJSDM-worktrees/lesson-0-teaching -b codex/lesson-0-teaching main
cd ~/src/occJSDM-worktrees/lesson-0-teaching
```

Read the report, the spec bullets and the lesson end to end. Render the before-state markdown as in Task 1 Step 3.

- [ ] **Step 2: Teaching commit**

Make these edits, then every other bullet of the report's section 2 at its line:
  - Prerequisites block after the setup chunk, verbs trimmed; cut the routine glosses the report lists (`readRDS()`, lists, `$`, `|>`, `seq()`, `rep()` and the dplyr verbs); keep the `anti_join()` and paired-row explanation.
  - Line 26 and line 461: replace the roadmap prose with one forward link each ("[Lesson 1](occJSDM-lesson-1.md) fits the model"; "the spatial lesson changes the ecological simulation"); the Phase B renumbering rewrites them.
  - Before `## Decide how collection and PCR can fail` (line 90): one passage laying out the two-stage process in the reader's terms (field collection can miss a present species, a false negative, or pick up DNA that is not from the site, a false positive; PCR can miss DNA in the sample or report DNA that is not there), then introduce `theta`, `theta0`, `p` and `q` as the four probabilities of those events before the settings chunk, and link to the intuition lesson for the one-stage version and to Lesson 1 for the full explanation.
  - After the settings table (now bullets): a paragraph on the good-practice assumption (Doug's bullet): the model assumes the survey was run so that contamination is infrequent, which is why the default priors on `q` and `theta0` are Beta(1, 20) with mean about 4.8%, and why this simulation's rates (2 to 8% and 2.5 to 6.5%) are chosen inside that assumption; Lesson 1 stress-tests it.
  - Convert the two pipe tables (lines 115 to 122 and 204 to 208) to bullets.
  - Keep "specify, then inspect"; add a `str()` or `head()` look at each object where it is first created (`observation_settings`, the ecological settings, the three input objects), each read in one sentence.
  - Read every printed output: the Site 3 and Sample 6 table (lines 263 to 266), a survey-wide `count(source)` of each kind of detection outcome added after line 266, a presence-by-detection cross-tabulation added after line 391 (`count(z, detected)` at the site level), and the OTU_1 against OTU_10 contrast in the maps (lines 333 to 393): say what the two species' different detectabilities do to their maps.
  - Line 256: split "No detection" into `True negative` (z == 0, no positive) and `Missed detection` (z == 1, no positive) in the `case_when`, and update the legend and the prose that reads the map.
  - Lines 430 and 434: rename the reassigned `samples_per_site` and `pcrs_per_primer` to `samples_per_site_unbalanced` and `pcrs_per_primer_unbalanced`; check every later use.
  - Lines 416 to 420 and 453 to 456: say that non-consecutive sample identifiers are accepted (Lesson 1's unbalanced fit kept the original identifiers) and that unequal numbers of field samples per site, primers per sample and PCRs per sample and primer are all accepted by the fitter, which counts the rows in each block and is tested on such a survey; point at the fit, not at `plotCumulativeSpeciesDetections()`, which takes single M and K values.
  - Appendix: the reproduction and provenance material at the end of the lesson moves under `## Appendix: evidence and reproduction`, with the closing section saying what it holds.
  Render, apply the PNG policy, commit "Lesson 0 teaching pass: the two-stage process, the good-practice assumption, every output read" with the trailer.

- [ ] **Step 3: Copyedit commit, checks, pull request**

As Task 2: copyedit; both renders; `test_lesson_links.R`; the site test; both counts 0 for Lesson 0; decisions-log line; commit "Lesson 0 copyedit pass"; push; `gh pr create --base main --title "Lesson 0 in the teaching style"` with the standard body (the "needs Alex" list is empty unless something new came up).

---

### Task 5: Lesson 2 pass

**Files:**
- Modify: `vignettes/occJSDM-lesson-2.Rmd` (445 lines), its `.md`, `vignettes/teaching-data/lesson-2-*.png`, `vignettes/LESSON-PLAN.md`
- Read: `style-diagnostic/lesson-2.md`, the spec's Lesson 2 bullets, `dev/simstudy/spatial-design-sweep/README.md` (for the numbers and the arrangement definitions)

**Interfaces:**
- Consumes: the bundle `vignettes/teaching-data/spatial-lesson.rds` (elements landscape, arrangements, statistics, oracle, fits, aggregate, reading, range_reading, paired, field_maps, lattice_maps, oracle_lattice, field_convergence, selected_fits, audit, provenance), unchanged.
- Produces: the lesson's closing advice paragraph, which the quickstart's spatial note can later link to.

- [ ] **Step 1: Worktree and before-state**

```bash
cd ~/src/occJSDM && git pull --ff-only origin main
git worktree add ~/src/occJSDM-worktrees/lesson-2-teaching -b codex/lesson-2-teaching main
cd ~/src/occJSDM-worktrees/lesson-2-teaching
```

Read the report, the spec bullets, the sweep README and the lesson end to end. Render the before-state markdown as in Task 1 Step 3 (this lesson's figures live under `teaching-data/`, so the PNG policy applies there).

- [ ] **Step 2: Teaching commit**

Make these edits, then every other bullet of the report's section 2 at its line:
  - Prerequisites block after the setup chunk; the lesson's wrangling code stays visible (spec decision).
  - Line 33: replace the bold central question with the one the sweep answers: how should 100 sites be placed if the spatial field is to be learned at all, and what does the answer cost in coverage.
  - Line 55: replace the PR 8 reference with the reasoning for using every surveyed location as a support point (the field is represented at the knots; with 100 sites, 100 knots lose nothing and keep the basis exact at the data).
  - Lines 71 to 77: keep the close-pairs advice with one sentence that the sweep found no measurable benefit from 20 close pairs at this budget and range.
  - Line 83: cut the sentence about an earlier computational validation.
  - Move the data-loading chunk (lines 93 to 116) to the start of 2A, after a sentence saying what the bundle holds and that it ships with the package's vignettes once the lessons are published, so the exercises run from the installed package.
  - Before the sweep's first result (2A): name, in the reader's terms, the species prevalence groups (two species at 5%, three at 25%, three at 75% of lattice cells), how each arrangement was laid out (spread at random; the first 80 spread sites plus 20 partners 0.01 away; ten clusters of ten within a radius of 0.02; a 10 by 10 grid), the three arms (an oracle given the true states and parameters; occJSDM fitted to the true states, which is its JSDM-only mode; occJSDM fitted to a two-stage eDNA survey with two samples, two primers and six PCRs), and the error measure (RMSE of the centred field at the surveyed sites against the true field, relative to assuming a flat field).
  - Move the reading-rule definitions (line 235) before their first use (line 223), and the reading-label table with its paragraphs (lines 225 to 237) from 2B into 2C, after the fit-recovery figure.
  - Move the convergence material (lines 334 to 345) to a new `## Appendix: evidence and reproduction`, together with the reproduction record; the lesson's closing section says what the appendix holds.
  - Check every conclusion typed beside a computed number (lines 237, 270, 345, 420, 426) against the rendered output and compute it inline where the bundle has it.
  - Close the lesson, before the exercises, with explicit advice matching the quickstart's spatial note: leave `spatCovariates` out in the beta unless sites are clustered within about one range of each other and the species of interest are common; and show, in one short chunk, how to put a plausible range for the reader's own system on the model's standardised scale (range divided by the standard deviation of the site coordinates on each axis, which must lie between 0.01 and 0.30 to be on the grid), using `sweep$statistics` columns `standardised_range_x` and `standardised_range_y` as the worked example.
  Render both outputs, open each changed `teaching-data/lesson-2-*.png` with the Read tool and confirm only moved or changed figures differ, commit "Lesson 2 teaching pass: name the design before the result, rules before labels, closing advice" with the trailer.

- [ ] **Step 3: Copyedit commit, checks, pull request**

As Task 2, plus the lesson's own verifier, which needs no raw archive:

```sh
Rscript dev/simstudy/spatial-design-sweep/verify-lesson.R .
```

Expected: "Lesson 2 bundle verified against the committed results." Decisions-log line; commit "Lesson 2 copyedit pass"; push; `gh pr create --base main --title "Lesson 2 in the teaching style"` with the standard body. The "needs Alex" list is empty; the pull request notes that the lesson will move to the end of the sequence in Phase B and keeps its current number until then.

---

## Self-review notes

- Spec coverage: Phase A's four items (Lesson 1, intuition lesson, Lesson 0, Lesson 2) are Tasks 1 to 5; the two-commit rule, the prerequisites block and tutorial link, the R-gloss decision, the appendix decision, the pipe-table and inline-number rules, the cross-reference rule and the per-lesson bullets of the spec each appear in a step. The Lesson 4 sjSDM item and all Lesson 3 and 4 work are Phase B and C, not in this plan.
- The one-stage likelihood in Task 3 handles each site's two samples with a closed form for the no-detection case and a product of Bernoulli terms otherwise; the covariate version aggregates per-sample log terms within site before mixing with psi, which is the standard single-season occupancy likelihood. The expected check in Step 5 (psi and p within 0.1 of their simulation values) catches a sign or mixing error.
- Type consistency: `truth_samples`, `survey_info`, `one_stage`, `site_history`, `estimates`, `p_hat`, `psi_hat`, `M` and `history_cov` are defined before use in the order the chunks appear; `scales::percent` needs the scales package, which ggplot2 imports.
- Registration of the new stem touches five files; the Node test and the link test both enumerate stems, so Step 1 of Task 3 runs the Node test before any writing.

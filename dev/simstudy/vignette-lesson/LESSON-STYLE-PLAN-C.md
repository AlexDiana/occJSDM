# Lesson Rewrite Phase C Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give Lessons 3 (outputs), 4 (prediction), 5 (four-JSDM comparison) and 6 (ten communities and traits) a teaching pass and a separate copyedit pass each, give Lesson 7 (spatial) its edges only, and clear the small cross-lesson items Phase B parked, in one branch and one pull request.

**Architecture:** Every lesson is an R Markdown file under `vignettes/` that renders from committed teaching bundles without fitting anything. A teaching pass edits prose, adds small run chunks that show data or read an output, adds not-run chunks that show a call the reader will make, moves or cuts repeated passages, and computes from the bundles every number the prose states; a copyedit pass then splits long sentences, fixes stops and wording, removes repeated asides and checks every claim against the rendered output. Chunks that a verifier compares or executes stay byte-identical, so no bundle, exporter or verifier changes. Each lesson is rendered to unwrapped markdown (committed) and HTML (gitignored, for Doug's read), checked by its verifiers and the shared link and site tests.

**Tech Stack:** R 4.5 with rmarkdown and knitr, dplyr, tidyr, tibble, ggplot2 and ggtern; the committed bundles in `vignettes/teaching-data/`; the occJSDM library in the full-fit archive for the verifiers that read fits; the verifiers in `dev/simstudy/vignette-lesson/`, `dev/simstudy/jsdm-package-comparison/` and `dev/simstudy/spatial-design-sweep/`; Node for the site test; perl and awk for the prose checks; git.

**Spec:** `dev/simstudy/vignette-lesson/LESSON-STYLE-SPEC.md`: "Sequence" (Phase C), "The reader", "The criteria, applied", "Decisions taken on 2 October 2026", "Lesson 3" and "Lesson 4" (the old Lesson 4 bullets now apply to Lessons 5 and 6), "What does not change" and "Verification for every lesson pull request". The owner's Phase C decisions of 3 October 2026 in `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/notes.md` override the spec where they differ: Lesson 7 gets edges only; Lesson 5 cites sjSDM fork release v0.1.0 now, drops the `SJSDM_MOJO_BACKEND` line, states the code-diff argument as an argument and adds a TODO item for the v0.1.0 rerun (no environment setup, no downloads); Lesson 3 has no line target. The diagnostic reports are `dev/simstudy/vignette-lesson/style-diagnostic/lesson-3.md` (now Lessons 3 and 4) and `lesson-4.md` (now Lessons 5 and 6); their line numbers refer to the lessons before Phase B, so every task below maps each item to the current lesson and section by quoted text. Format models: `LESSON-STYLE-PLAN-A.md` and `LESSON-STYLE-PLAN-B.md` in the same directory.

## Decisions this plan takes

Each is recorded so the reviewer and the owner can see it; the final report lists the ones the owner may want to revisit.

- **Task shape.** Lesson 4 (362 lines) is one task with two commits, teaching then copyedit. Lessons 3, 5 and 6 each have a teaching task and a copyedit task. Lesson 7's edges and the parked cross-lesson items are one task. A final task runs every check and drafts the pull request.
- **Shared tools.** Phase B's `run-step.zsh` and `revert-unchanged-png.zsh` change directory to the removed Phase B worktree, so Task 1 copies them into `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/` with the worktree path changed. `check-lesson-links.R` has no path in it and is reused from Phase B's directory unchanged. Two new tools are added beside the copies: `check-frozen-chunks.R` and `prose-checks.zsh`, both given in full below.
- **Frozen chunks.** Thirty-six chunks are compared or executed by verifiers (listed in `check-frozen-chunks.R`). They stay byte-identical, header line included, and may only move as whole blocks. Where a report item asks to change one (show one of the two near-identical native trait chunks; write the MCMC settings literally in `prediction-fit-one-factor`), the prose does the job instead and the pull request lists the item as deferred for that reason.
- **The ggtern theme.** Checked while planning, on 3 October 2026 with ggplot2 4.0.3 and ggtern 4.0.0: after `loadNamespace("occJSDM")` (which loads ggtern, an import), `theme_set(ggplot2::theme_bw())` makes the next plot fail with the error `The tern.axis.ticks.length.major theme element must be a <rel> object.`, while adding `theme_bw()` to a single plot and `theme_set(ggtern::theme_bw())` both work. Lesson 3 states this symptom in one sentence; why ggtern's elements fail ggplot2's check, and whether the package can avoid it, goes on the "needs Alex" list and into a `TODO.md` item. Lessons 4, 5 and 6 point to Lesson 3's sentence.
- **Hmsc schedule.** `dev/simstudy/jsdm-package-comparison/extension/fit-job.R` lines 27 and 28 and `extension/PLAN.md` line 27 confirm the schedule: four chains, 2,000 warm-up and 4,000 retained draws, then at most one longer attempt with 8,000 and 16,000; a fit still failing after it stays flagged. Lesson 6 states this, and the statement goes on the "needs Alex" list as the task instructions require.
- **New-site covariate spread.** The report doubted "standard deviation 10" for the new sites' environment (Lesson 4). `dev/simstudy/vignette-lesson/prediction-build.R` line 31 draws them with `rnorm(n_new * 2L, sd = 10)`, so the claim stands and only gains its source.
- **Link label.** "Quickstart and lesson guide" becomes "Quickstart" in Lessons 5, 6 and 7 now, because the quickstart says the lessons will follow after the beta and lists none.
- **Decision-log placeholders.** The five dated entries in `vignettes/LESSON-PLAN.md` that still say "(link added when opened)" get their pull request links (#19 for the four Phase A entries, #20 for Phase B), as the Phase A and B controller steps intended; no other dated text changes.
- **PR #14.** Lesson 4's WAIC subsection and Lesson 3's WAIC paragraph are edited (an undefined term and a sentence of history), which adds to the hand re-application PR #14 already needs. The pull request's collision paragraph names exactly the changed paragraphs.
- **Visible code stays visible.** The report suggested `echo=FALSE` for Lesson 3's truth-join chunk; it stays visible because every lesson shows all teaching code (Phase A decision for Lesson 7, applied to all).
- **Lesson 3 ordination order.** The unverified `extract-ordination` chunk duplicates three lines of the frozen `ordination-standard` chunk; it is deleted and `ordination-standard`, its paragraph and the ordinary-biplot image move to the head of the section, so the reader sees the calls they will run before the truth comparison.
- **Stale extension dev copies.** `extension/lesson-section.Rmd` and `extension/README.md` are hash-stamped sources of `vignettes/teaching-data/lesson-4-extension.rds` (its `sources` element), so they stay unchanged, as Phase B ruled; the pull request lists them as deferred.
- **No line target for Lesson 3.** The pass cuts repetition and seams and adds teaching; the pull request reports the before and after line counts without a target.

## Global Constraints

- Work only in the worktree `/Users/douglasyu/src/occJSDM-worktrees/lesson-phase-c`, on branch `codex/lesson-phase-c` (created from main at 7662f9a). Never work in `/Users/douglasyu/src/occJSDM`, never commit to main, never push (the controller pushes), never force anything (no `push --force`, no `reset --hard` of committed work, no `-f` on `git mv`). Reading the extension archive under `/Users/douglasyu/src/occJSDM/dev/simstudy/results/` for its verifier is read-only use.
- The full-fit archive is `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples`. Verifiers that read fits run from the worktree root through the Phase C step runner below, so each run is logged in the archive's `BUILD-LOG.md`. Agent shells keep no variables between calls, so commands spell paths out in full.
- Commit messages end with a blank line and `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`. The pull request body (drafted in Task 9, opened by the controller) ends with `🤖 Generated with [Claude Code](https://claude.com/claude-code)` and carries a "needs Alex" list containing only package facts the lessons newly state that the spec marks unconfirmed.
- Markdown rules for every `.md` and `.Rmd` touched, from `CLAUDE.md` and `AGENTS.md` lines 89 to 124: one line per paragraph, no hard wrapping; no em-dashes anywhere, including commit messages (use `--`); no pipe tables (bullet lists instead); escape a literal `>` `<` `~` `|` in prose as `\>` `\<` `\~` `\|` and write `-\>` for an arrow, never inside code spans or fenced blocks; keep inline code spans short; never add an `editor_options` block; no double blank lines. Re-read a file immediately before and after editing it (it may be open in RStudio); if an edit looks corrupted or duplicated, rewrite the whole file with Write.
- Writing rules, from the spec's "The reader" and "The criteria, applied" and `AGENTS.md`: the reader is an ecologist who runs eDNA or presence/absence surveys, knows basic R and the pipe, and has not fitted an occupancy model or a JSDM; explanations remind rather than introduce. Each section must name the concepts the reader needs in the reader's terms before using them, answer the question the reader is about to ask, show every part of the data or output the reader will handle, read every printed output in the reader's terms, spell out the reasoning behind every choice, and give practical advice for the reader's own survey. Cut glosses of the pipe, `select()`, `filter()`, `mutate()`, `summarise()`, `group_by()`, `readRDS()`, lists and `$`; keep a gloss only for an unusual operation (a join on species identity, a pivot, an array slice). Refer to a lesson by what it teaches where a number would otherwise be the only identifier.
- Numbers: every number the prose states that a bundle holds is computed inline with `` `r ...` `` from the bundle, so the prose cannot go stale. A number no bundle holds is stated with its source (a file and line, or a build-log step), or not stated. Exporters and bundles are out of scope, so nothing is added to a bundle.
- Out of scope: package code (`R/`, `src/`, `DESCRIPTION`, `NAMESPACE`); every teaching bundle under `vignettes/teaching-data/` (`*.rds`, `*.csv`, `lesson-4-calibration-tables.R`) and every exporter-written PNG there (`native-plot-*`, `native-traits-*`, `remaining-plots-*`, `ordination-*`); every exporter and verifier (a prose pass changes none); the four frozen snippet files `native-plot-examples.Rmd`, `native-traits-examples.Rmd`, `ordination-examples.Rmd` and `remaining-plots-examples.Rmd` in `dev/simstudy/vignette-lesson/`, and the two unbalanced mirrors beside them; everything under `dev/simstudy/spatial-design-sweep/`; the hash-stamped sources under `dev/simstudy/jsdm-package-comparison/` (scripts, `extension/lesson-section.Rmd`, `extension/calibration-section.Rmd`, `extension/README.md`, `extension/PLAN.md`); historical documents (completed plans and reports, `style-diagnostic/`, the dated decision log of `vignettes/LESSON-PLAN.md` apart from the five link placeholders, struck-through `TODO.md` items).
- The 36 frozen chunks (listed in `check-frozen-chunks.R` below) stay byte-identical, header line included. They may move only as whole blocks. Explain or qualify them in the prose around them.
- Headings that other files link to keep their exact text, so their anchors survive. Lesson 3: "What does an effect mean for a species' distribution?", "Put the true coefficient beside its estimate", "Residual species associations: did we recover what was put in?", "Ordination: compare the combined effect before naming the axes", "Check computation as well as ecological recovery", "When chains settle on two different explanations", "References and further reading", and the appendix headings "A real cancellation inside this simulated community", "Read native ordination plots after aligning their axes", "The mirror-labelling study", "Read trace draws from the fit object", "Reproduce the teaching figures". Lesson 4: "Predict occupancy at genuinely new sites", "Which true probability should a new-site prediction recover?", "Compare models using what actually occurred", "Check the additional fit and understand the WAIC limitation", "Appendix: evidence and reproduction". Lesson 5: "Exercise: predict one species given another", "The two sjSDM optima: evidence", "Reproduction record". Lesson 6: "Which fits passed, remained flagged or failed?", "Appendix: evidence and reproduction". A heading not on this list may change only if `check-lesson-links.R` still passes on every lesson and living document afterwards.
- Line numbers in this plan are from main at 7662f9a; the reports' line numbers are from before Phase B. Earlier tasks shift lines. Executors find the quoted heading, chunk label or sentence, and never trust a line number alone.
- Shell is zsh and `grep` and `ls` are aliased: use `/usr/bin/grep` and `/bin/ls` in commands.
- Render commands, run from the worktree root (replace `<stem>`). The four `library()` calls keep R 4.5's "built under" warnings out of the rendered `.md`:

```sh
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/<stem>.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/<stem>.Rmd", output_format="rmarkdown::html_vignette")'
```

Both must finish without a new warning in the console or in the `.md` (`/usr/bin/grep -c "built under" vignettes/<stem>.md` is 0, and any `Warning` line in the `.md` existed at 7662f9a or is removed by the pass). The HTML output is gitignored and never committed.

- PNG policy. A render rewrites every figure it draws. Keep a changed PNG only when its chunk changed on purpose; revert the others with the Phase C revert script, which compares pixels with the index copy. Four figures vary between renders for reasons unrelated to content (legend order or anti-aliasing): `vignettes/occJSDM-lesson-0_files/figure-gfm/environmental-maps-1.png`, `vignettes/occJSDM-lesson-2_files/figure-gfm/probability-and-state-1.png`, `vignettes/occJSDM-lesson-3_files/figure-gfm/mirror-labelling-chains-1.png` and `vignettes/occJSDM-lesson-3_files/figure-gfm/trait-cancellation-1.png`. If one of them is reported "kept, content changed" and its chunk did not change, restore it with `git checkout -- <path>`. Never point the revert script at exporter-written PNGs (`vignettes/teaching-data/native-plot-*`, `native-traits-*`, `remaining-plots-*`, `ordination-*`); renders do not write them. Lessons 5, 6 and 7 write their markdown figures to `vignettes/teaching-data/lesson-5-*`, `lesson-6-*` and `lesson-7-*`; Lessons 3 and 4 to `vignettes/occJSDM-lesson-3_files/` and `occJSDM-lesson-4_files/`. Check with `git status --short vignettes/` before every commit.
- Checks that must pass at the end of every task (`<files>` are the Rmd and md files the task changed):

```sh
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
node --test dev/simstudy/lesson-site/test_lessons.js
Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/check-frozen-chunks.R
Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/check-lesson-links.R vignettes/occJSDM-lesson-[0-7].Rmd vignettes/occJSDM-lesson-[0-7].md vignettes/LESSON-PLAN.md TODO.md
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/prose-checks.zsh <files>
git diff --name-only 7662f9a -- ':!*.png' ':!*.rds' | xargs /usr/bin/grep -n $'\u2014' | wc -l
git log 7662f9a..HEAD --format=%B | /usr/bin/grep -c $'\u2014'
```

The `$'\u2014'` argument is zsh's spelling of the em-dash character, so this plan never contains it. Expected: the link test ends "Shared lesson-link regression checks passed."; the site test reports `fail 0`; "All 36 frozen chunks match main at 7662f9a."; "Every lesson link in ... files names an existing lesson and heading."; every "must be 0" count in the prose checks is 0 and every listed long sentence has been read and either split or kept for a stated reason (a list, a quoted title or a table caption); the last two counts are 0.

- Stop and report to the controller, without working around it, if: a baseline verifier fails before this plan has changed anything; a verifier fails after an edit for any reason other than an accidental change to a frozen chunk (revert the chunk; never change a verifier); a render needs a bundle change; an out-of-scope file would need a change; the spec, `notes.md` and this plan conflict; or a reviewer finding conflicts with the spec, `notes.md` or this plan. Do not stop between lessons otherwise (owner's rule for Phase C).

---

## Shared tools (created in Task 1 Step 1, kept outside the repository)

All four Phase C tools live in `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/`, the controller's working directory. They are not committed.

`run-step.zsh` and `revert-unchanged-png.zsh` are Phase B's scripts (printed in full in `LESSON-STYLE-PLAN-B.md`, "Shared tools") with every `lesson-phase-b` replaced by `lesson-phase-c`. Usage is unchanged: `zsh run-step.zsh STEPNAME command args...` logs to the archive's `logs/STEPNAME.log` and appends a record to its `BUILD-LOG.md`; `zsh revert-unchanged-png.zsh PATHSPEC` reverts every modified tracked PNG under the pathspec whose pixels equal the index copy. Phase C step names start with `c` (for example `c11-L3-native-plots`), so they never collide with Phase B's numbered steps.

`check-lesson-links.R` is used from `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/` unchanged.

`check-frozen-chunks.R` (tested while planning: "All 36 frozen chunks match main at 7662f9a."):

```r
# Usage, from the worktree root: Rscript check-frozen-chunks.R
# Fails if a chunk that a verifier compares or executes differs from main at 7662f9a.
frozen <- list(
  "vignettes/occJSDM-lesson-0.Rmd" = c("unbalanced-select-samples", "unbalanced-paired-rows"),
  "vignettes/occJSDM-lesson-2.Rmd" = c("unbalanced-load", "unbalanced-fit-optional"),
  "vignettes/occJSDM-lesson-3.Rmd" = c(
    "native-plot-setup", "native-environment-example", "native-environment-2-example",
    "native-environment-perfect-1-example", "native-environment-perfect-2-example",
    "native-collection-example", "native-baseline-occupancy-example",
    "native-baseline-collection-example", "native-primer-1-example", "native-primer-2-example",
    "native-correlation-example", "native-effort-k-example", "native-effort-m-example",
    "native-traits-setup", "native-traits-gradient-1", "native-traits-gradient-2",
    "native-traits-perfect-gradient-1", "native-traits-perfect-gradient-2",
    "ordination-standard", "ordination-align", "ordination-sites", "ordination-loadings",
    "ordination-biplot",
    "remaining-plots-setup", "remaining-plots-gradient-1", "remaining-plots-gradient-2",
    "remaining-plots-stage1-fp", "remaining-plots-stage2-fp", "remaining-plots-detection"),
  "vignettes/occJSDM-lesson-4.Rmd" = c("prediction-fit-one-factor", "prediction-native-call"),
  "vignettes/occJSDM-lesson-5.Rmd" = c("reproduce-community")
)
chunk <- function(lines, label) {
  start <- grep(paste0("^```\\{r ", label, "[,}]"), lines)
  if (length(start) != 1L) return(NULL)
  end <- start + which(lines[(start + 1L):length(lines)] == "```")[1]
  lines[start:end]
}
bad <- character()
for (path in names(frozen)) {
  now <- readLines(path, warn = FALSE)
  base <- system2("git", c("show", paste0("7662f9a:", path)), stdout = TRUE)
  for (label in frozen[[path]]) {
    a <- chunk(base, label)
    b <- chunk(now, label)
    if (is.null(a) || is.null(b) || !identical(a, b)) bad <- c(bad, paste(path, label))
  }
}
if (length(bad)) {
  writeLines(c("Changed or missing frozen chunks:", bad))
  quit(status = 1)
}
cat("All", sum(lengths(frozen)), "frozen chunks match main at 7662f9a.\n")
```

`prose-checks.zsh` (tested while planning on Lesson 4 and `LESSON-PLAN.md`):

```zsh
#!/bin/zsh
# Usage, from the worktree root: zsh prose-checks.zsh FILE...
# Prints every rule breach in the given Rmd or md files; the counts under "must be 0" must be 0.
print "em-dashes (must be 0): $(/usr/bin/grep -c $'\u2014' "$@" | awk -F: '{s+=$NF} END {print s+0}')"
print "pipe-table lines (must be 0): $(/usr/bin/grep -c '^|' "$@" | awk -F: '{s+=$NF} END {print s+0}')"
print "editor_options blocks (must be 0): $(/usr/bin/grep -c 'editor_options' "$@" | awk -F: '{s+=$NF} END {print s+0}')"
print "double blank lines (must be 0):"
awk 'FNR==1{prev="x"} prev=="" && $0=="" {print "  " FILENAME ":" FNR} {prev=$0}' "$@"
print "unescaped -> or ~ in prose outside code (must be 0):"
perl -ne 'if(/^```/){$in=!$in} elsif(!$in && !/<!--/){my $l=$_; $l=~s/`[^`]*`//g; print "  $ARGV:$.: $_" if $l=~/(?<!\\)->|(?<![\\~])~(?!~)/} if(eof){close ARGV; $in=0}' "$@"
print "and/or (must be 0): $(/usr/bin/grep -c 'and/or' "$@" | awk -F: '{s+=$NF} END {print s+0}')"
print "sentences over 30 words in prose (review each; split unless a list or a quoted title):"
perl -ne 'if(/^```/){$in=!$in} elsif(!$in && /\S/){my $l=$_; $l=~s/`r [^`]*`/N/g; $l=~s/`[^`]*`/C/g; for my $s (split /(?<=[.!?])\s+(?=[A-Z*\[])/, $l){my @w=split /\s+/, $s; printf "  %s:%d: %d words: %.80s\n", $ARGV, $., scalar(@w), $s if @w > 30}} if(eof){close ARGV; $in=0}' "$@"
```

Run the prose checks on Rmd files, not on rendered `.md` files: the rendered markdown escapes and wraps differently, and its counts are not meaningful.

## Shared text

**The prerequisites block.** Lessons 3, 4, 5 and 6 get Plan A's block as its own paragraph plus a bullet list, directly after the paragraph that says how to run or knit the lesson and before the load chunk, with the verbs trimmed to what that lesson's code uses. List the verbs a lesson uses with:

```sh
/usr/bin/grep -o -E "(^|[^a-z_.:])(select|filter|mutate|transmute|summarise|group_by|ungroup|arrange|distinct|pull|left_join|inner_join|anti_join|bind_rows|pivot_longer|pivot_wider|across|case_when|if_else|coalesce|recode|count|rowwise|expand_grid|tribble|as_tibble|enframe)\(" vignettes/<stem>.Rmd | sed -E 's/^[^a-z_]//' | sort | uniq -c | sort -rn
```

The block's paragraph, verbatim; then one bullet per group of verbs the lesson uses, in the form of Lesson 2's list (`vignettes/occJSDM-lesson-2.Rmd`, the bullets after "What this lesson assumes you know"):

```markdown
**What this lesson assumes you know.** The code uses base R and the tidyverse: the pipe `|>`, and from dplyr and tidyr the verbs listed below. If any are new, the two chapters of R for Data Science on [data transformation](https://r4ds.hadley.nz/data-transform) and [data tidying](https://r4ds.hadley.nz/data-tidy) teach everything used here in an afternoon. Operations that are unusual, such as joining two tables on species identity, are explained where they appear.
```

If a lesson's unusual operations are specific (Lesson 3's array slices and `rowwise()`; Lesson 6's `tribble()` and nested pivots), name them in the paragraph's last sentence instead of the generic example, as Lesson 7's block does.

**The ggtern sentence** (Lesson 3, after its load chunk; Lessons 4, 5 and 6 point to it):

```markdown
The last line of the chunk sets the plotting theme from the ggtern package, which installing occJSDM also installs. occJSDM uses ggtern for its ternary plots, and once ggtern is loaded, `theme_set(ggplot2::theme_bw())` makes the next plot fail ggplot2's check of a ternary theme element (`tern.axis.ticks.length.major`); ggtern's own `theme_bw()` draws the same theme with those elements in place, so these lessons use it. If you plot with ggplot2 after loading occJSDM, do the same, or add `theme_bw()` to each plot instead of calling `theme_set()`.
```

The pointer for Lessons 4, 5 and 6, placed after their load chunk (adjust "last line" if the theme line is not last):

```markdown
The last line sets ggtern's version of `theme_bw()`; [Lesson 3](occJSDM-lesson-3.md#what-this-lesson-answers) explains why the lessons use it.
```

**The copyedit procedure** (Tasks 2, 3, 5 and 7). Read the Rmd from top to bottom and, sentence by sentence: split sentences over about 30 words (the prose checks list them); replace every "and/or" and every slash pair other than "presence/absence" with words; fix missing stops; remove any aside that repeats its sentence; replace curly quotes with straight ones; check each factual claim in prose against the rendered output and correct it; make sure every number the bundles hold is computed inline. Do not remove or shorten a teaching addition from the teaching commit; if one is wrong, fix it. Section 5 of the lesson's diagnostic report lists the sentences it found; map each to the current text.

---

## File structure

- `vignettes/occJSDM-lesson-3.Rmd`, `.md`, `occJSDM-lesson-3_files/`: Tasks 1 and 2.
- `vignettes/occJSDM-lesson-4.Rmd`, `.md`, `occJSDM-lesson-4_files/`: Task 3.
- `vignettes/occJSDM-lesson-5.Rmd`, `.md`, `vignettes/teaching-data/lesson-5-*.png`: Tasks 4 and 5.
- `vignettes/occJSDM-lesson-6.Rmd`, `.md`, `vignettes/teaching-data/lesson-6-*.png`: Tasks 6 and 7.
- `vignettes/occJSDM-lesson-7.Rmd` and `.md`, `vignettes/occJSDM-lesson-1.Rmd` and `.md`: Task 8.
- `TODO.md`: Task 1 (ggtern item), Task 4 (sjSDM v0.1.0 rerun item), Task 8 (one double parenthetical).
- `vignettes/LESSON-PLAN.md`: Task 8 (last-updated line, line 29, the five link placeholders), Task 9 (the Phase C decision-log line).
- Not committed: the four shared tools, `lesson-N-before.md` copies and `pr-body.md` in `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/`.

## Known collision with PR #14 (open, `codex/site-waic`)

PR #14 edits the WAIC subsection that is now in `vignettes/occJSDM-lesson-4.Rmd` ("Check the additional fit and understand the WAIC limitation"), the WAIC paragraph at the end of Lesson 3's diagnostics section, `prediction-export.R` line 57, `prediction-verify.R` line 46, `vignettes/LESSON-PLAN.md` line 64 and `TODO.md`. Phase C edits the first two (Task 1 Step 7, Task 3 Step 8) and `LESSON-PLAN.md` and `TODO.md` elsewhere. Whichever merges second re-applies the other's hunks by hand, re-renders Lessons 3 and 4 and re-runs `prediction-verify.R`. Task 9's pull request body names the edited paragraphs.

---

### Task 1: Lesson 3, teaching commit

**Files:**
- Modify: `vignettes/occJSDM-lesson-3.Rmd` (1,822 lines), `vignettes/occJSDM-lesson-3.md`, `TODO.md`; `vignettes/occJSDM-lesson-3_files/figure-gfm/` only if a chunk's figure changes on purpose (none is expected)
- Create (not committed): the four shared tools; `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/lesson-3-before.md`
- Read: `style-diagnostic/lesson-3.md` sections 2 to 5 (items up to old line 785 and from old line 1086 on; old lines 786 to 1085 are Lesson 4's, Task 3), the spec's Lesson 3 bullets, `notes.md`, `AGENTS.md` lines 89 to 124, Lesson 1 (`vignettes/occJSDM-lesson-1.Rmd`) for its vocabulary

**Verifiers that read this lesson:** `verify_native_plots.R`, `native-traits-verify.R`, `ordination-verify.R`, `remaining-plots-verify.R` (each compares the lesson's frozen chunks with its snippet file and prints "Lesson 3 displays the exact exported ...").

**Interfaces:**
- Consumes: `teaching-data/nonspatial-lesson.rds` (`lesson`: `input`, `observations`, `cases`, `cells`, `groups`, `rates`, `samples`, `diagnostics` with `perfect`, `default`, `alternative`, `manifests` with each fit's `mcmc` list `nchain`, `nburn`, `niter`, `nthin` and `priors`), `teaching-data/output-lesson.rds` (`outputs`: `coefficients` with columns `covariate`, `term`, `block`, `truth`, `estimate`, `lower`, `upper`, `excludes_zero`, `rhat`, `ess_bulk`, `ess_tail`, `arm`; `gradients`, `correlations` with `species1`, `species2`, `truth`, `estimate`, `lower`, `upper`, `arm`; `residual`, `baseline`, `collection`, `detection_effort` with `M`, `K`, `estimate`, `lower`, `upper`, `truth`; `trait_components`, `trait_hidden_correlation`), `teaching-data/diagnostics-lesson.rds` (`diagnostic_examples$traces` with `collection`, `detection`, `field_contamination`), `teaching-data/native-plots.rds` (`native_examples$truth$survey_counts`).
- Produces: the anchor `#what-this-lesson-answers` carrying the ggtern sentence (Tasks 3, 4, 6 link to it); the kept headings listed in the Global Constraints; the term definitions (posterior draw, chain, credible interval, log-odds, hidden site factors) that Lessons 4 to 6 point to instead of repeating.

- [ ] **Step 1: Create the shared tools and run the baseline**

```sh
mkdir -p /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C
sed 's#lesson-phase-b#lesson-phase-c#g; s#Phase B worktree#Phase C worktree#' /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/run-step.zsh > /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/run-step.zsh
sed 's#lesson-phase-b#lesson-phase-c#g' /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/revert-unchanged-png.zsh > /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/revert-unchanged-png.zsh
/usr/bin/grep -c "lesson-phase-c" /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/run-step.zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/revert-unchanged-png.zsh
/usr/bin/grep -c "lesson-phase-b" /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/run-step.zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/revert-unchanged-png.zsh
```

Expected: `run-step.zsh:2` and `revert-unchanged-png.zsh:1` for `lesson-phase-c`, and 0 for both files for `lesson-phase-b`. Write `check-frozen-chunks.R` and `prose-checks.zsh` into the same directory with the Write tool, exactly as printed under "Shared tools". Then run the baseline from the worktree root:

```sh
S=/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/run-step.zsh; A=/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/check-frozen-chunks.R
zsh $S c01-base-native-plots Rscript dev/simstudy/vignette-lesson/verify_native_plots.R $A
zsh $S c02-base-native-traits Rscript dev/simstudy/vignette-lesson/native-traits-verify.R $A
zsh $S c03-base-ordination Rscript dev/simstudy/vignette-lesson/ordination-verify.R $A
zsh $S c04-base-remaining-plots Rscript dev/simstudy/vignette-lesson/remaining-plots-verify.R $A
zsh $S c05-base-prediction Rscript dev/simstudy/vignette-lesson/prediction-verify.R $A $A/prediction
zsh $S c06-base-unbalanced Rscript dev/simstudy/vignette-lesson/unbalanced-verify.R $A/unbalanced
zsh $S c07-base-verify-teaching Rscript dev/simstudy/jsdm-package-comparison/verify-teaching.R archive dev/simstudy/jsdm-package-comparison
Rscript dev/simstudy/jsdm-package-comparison/extension/verify.R /Users/douglasyu/src/occJSDM/dev/simstudy/results/lesson-4-extension-20260923 .
Rscript dev/simstudy/spatial-design-sweep/verify-lesson.R .
```

Expected: "All 36 frozen chunks match main at 7662f9a."; every runner line ends `EXIT 0`, the four Lesson 3 logs contain "Lesson 3 displays the exact exported ...", and the last lines are "All prediction verification checks passed.", the unbalanced "... and displayed code." line, "Teaching numerical verification passed for every check the committed archive supports.", "... 626 scored results; 640 of 640 fits attempted." and "Lesson 2 bundle verified against the committed results." (the sweep verifier keeps its old wording; the sweep is out of scope). Any failure here is a stop.

- [ ] **Step 2: Read, then keep the before state**

Read `style-diagnostic/lesson-3.md` in full, the spec's Lesson 3 bullets and `notes.md`, then `vignettes/occJSDM-lesson-3.Rmd` once from top to bottom as the reader would, before editing anything. Keep the report open: every section 2 item for old lines 27 to 785 and 1086 to 2085 is a step below, done or recorded as deferred. Keep the before state:

```sh
git show 7662f9a:vignettes/occJSDM-lesson-3.md > /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/lesson-3-before.md
wc -l vignettes/occJSDM-lesson-3.Rmd
```

- [ ] **Step 3: The opening ("What this lesson answers")**

Map and make these edits, in the reader's terms:
- Report old 29 (purpose stated abstractly): after the first paragraph's three contrasts, add one sentence: with real data the truth is never known, so a simulation whose truth we know is the only place to see what each output can and cannot tell you.
- Report old 31 ("It replaces the output tour in the original occJSDM vignette"): replace that sentence with one saying how this survey differs from the quickstart's `sampledata` (three field samples, two primers and six PCRs per primer here; three samples, three primers and two PCRs there), and that the model reads the structure, not the counts, as Lesson 2's fitting reference explains. Delete "You do not need [the spatial lesson](occJSDM-lesson-7.md) first." (Phase B deferral: the spatial lesson is now last).
- Report old 33 (why the perfect-observation fit exists): add that it shows what the occupancy part recovers once detection error is removed, so the gap between the two fits is the cost of imperfect detection.
- Report old 33 (maintainer sentences): cut "This is not a before/after comparison of software versions." and "The fits use the verified model source recorded in Lesson 2; that source matches the code on main when this lesson was prepared."; the appendix's reproduction section already points to the records.
- Spec "hidden site factors where they are first mentioned": at "changing only the number of hidden site factors from two to one", gloss hidden site factors in a clause (unmeasured conditions at a site that several species respond to) with a link to [Lesson 1](occJSDM-lesson-1.md#what-a-joint-model-does).
- Report old 35 (`eval=FALSE`, knit, compact summaries): in the "All teaching code is visible" paragraph, say in one sentence that some chunks are shown but not run while knitting, because they need a full fit that takes too long to build, and that knitting means rendering this file.
- Phase B deferral (`fitmodel_perfect` used before its load chunk): in the paragraph that introduces `fitmodel`, add that the perfect-observation fit is called `fitmodel_perfect` in the optional examples and that the [reproduction instructions](#reproduce-the-teaching-figures) load both.
- Add a link to the function finder: "The [function finder](#find-the-function-for-your-question) near the end maps each question to the functions that answer it."

- [ ] **Step 4: Prerequisites block, the teaching objects, the terms and the theme**

Insert the prerequisites block (Shared text) after the "All teaching code is visible" paragraph and before the `load-results` chunk; Lesson 3's verbs are `filter()`, `mutate()`, `transmute()`, `select()`, `arrange()`, `distinct()`, `summarise()`, `group_by()`, `ungroup()`, `left_join()`, `bind_rows()`, `as_tibble()` and `enframe()`; name `rowwise()`, `expand_grid()` and array slices such as `draws[, "X_psi.EnvCov.1", "OTU_1"]` as the unusual operations explained where they appear. Re-run the verb command and adjust if this pass adds a verb.

Directly after the `load-results` chunk, add the ggtern sentence (Shared text), then this chunk (report old 39 to 64; spec "show the three teaching data objects after loading"):

```{r show-teaching-objects}
names(lesson)
names(outputs)
names(diagnostic_examples$traces)
head(outputs$coefficients)
```

Follow it with a bullet list that reads the output: `lesson` is Lesson 2's bundle (the simulated survey and its truth in `input`, one row per PCR observation with where each positive came from in `observations`, the four selected cases, per-site summaries, the saved diagnostic tables of the three fits and each fit's settings in `manifests`); `outputs` holds this lesson's summaries of the full fits, one table per output (name the tables the lesson uses); `diagnostic_examples$traces` holds the three trace excerpts of the diagnostics section; `known_truth` is the simulation's generating parameters, which no fit saw. Say they are teaching summaries exported from the full fits, not what `runOccJSDM()` returns, and name `outputs$coefficients`' columns in the reader's terms (report old 74 to 86).

Then, before "In the figures, **black crosses or lines show truth**", add a short bullet list headed "**Terms used throughout.**" (spec: "Explain posterior draws, chains, credible intervals and log-odds at first use"), one sentence each, in Lesson 1's and Lesson 2's vocabulary:
- a posterior draw: one plausible set of parameter values given the data and priors; a fit keeps thousands, and summaries are taken over them;
- a chain: one independent run of the sampler; the fits here ran `r lesson$manifests$default$mcmc$nchain` chains of `r format(lesson$manifests$default$mcmc$niter, big.mark = ",")` kept draws, `r format(with(lesson$manifests$default$mcmc, nchain * niter), big.mark = ",")` in all;
- a 95% credible interval: the range holding the middle 95% of the draws, the model's statement of uncertainty given its assumptions, not a measured distance from truth;
- log-odds: the scale the model adds effects on, `log(p / (1 - p))`; zero is 50%, and `plogis()` converts back to a probability;
- hidden site factors: as glossed in the opening, with the Lesson 1 link.

- [ ] **Step 5: The diagnostics section ("Check computation as well as ecological recovery")**

- Report old 1548 (overview after details): move the `fitting-diagnostics` chunk and the paragraph after it ("These are the saved parameter diagnostics, not a certification ...") from "### Summaries and what to do next" to directly after the section's second paragraph ("Think of each MCMC chain ..."), under a new heading `### Start with an overview of each fit`, and read the table there: name each fit's largest Rhat and smallest ESS inline from `diagnostic_summary`.
- Report old 1273 (thresholds before their justification): move the sentences that justify the 1.01 and 400 screens (from "These are screening rules, not sharp boundaries" to the end of that paragraph) to directly before the `flag-default-diagnostics` chunk, and add the practical sentence the report asks for at old 1300: use the table-returning function's `rhat` and `ess` for the screen, and the newer `posterior` diagnostics when an interval endpoint matters (tail ESS).
- Report old 1287 ("weaker low-contamination assumptions"): name the priors inline: the alternative fit used Beta(`r lesson$manifests$alternative$priors$a_q`, `r lesson$manifests$alternative$priors$b_q`) on both `q` and `theta0` (check that `a_theta0` and `b_theta0` match before writing "both"), against the default Beta(1, 20), with a link to [Lesson 2](occJSDM-lesson-2.md#what-changes-if-we-are-less-confident-about-low-contamination).
- Report old 1240 and 1402 (`/path/to/full-fits/...`): before `extract-full-fit-diagnostics`, say that the path is where the archived teaching fits live and that with your own data you skip these two lines and use your own `fitmodel`. Same sentence, short form, before `extract-field-trace` in the appendix.
- Report old 1324: cut "These are excerpts of the full fit, not newly fitted models."
- Report old 1440 aside: cut "and not a claim that only that many iterations were run" from the 48,000-draws paragraph, and compute 48,000 inline with `with(lesson$manifests$alternative$mcmc, nchain * niter)`.
- Report old 1485 (`B0` unexplained) and Phase B Task 3 Minor 2 (per-chain guidance stranded in the appendix): in the `extract-chain-means` paragraph and the warning-sign paragraph after it, say that `B0_output` holds the occupancy intercepts; move the appendix sentence "Compare chains on those four, not on average occupancy or `B0` alone." into the end of the warning-sign paragraph; move the appendix paragraph that begins "`theta0_output` and `jsdm_output$B0_output` are species by iteration by chain" to directly after the `extract-chain-means` chunk. The chunk itself is not frozen, but keep its code unchanged.
- Report old 1487 (dense array caveats): move "a fit with one species or one chain loses array dimensions when a row is selected" out of the prose into a one-line code comment in `extract-chain-means` only if the chunk can take it without other change; otherwise shorten the sentence to its practical instruction (use at least two chains and two species).
- Report old 1569 (repetition): the paragraph after `fitting-diagnostics` repeats the caveats of the screen paragraph; keep one statement that diagnostics check computation, not ecological accuracy, and cut the repeats.
- Report old 1571: in "The public diagnostics table does not include the trait matrix `G`", replace "the exporter also calculates diagnostics directly from ..." with the reader-level point: trait coefficients are not in the package's table, so this lesson computed their Rhat and ESS from the saved draws; the appendix shows how (`#where-to-find-other-parameter-draws`).
- Report old 1594 (WAIC repetition) and spec "WAIC before it is used": reduce the WAIC paragraph at the end of the section to two sentences: one defining WAIC (the widely applicable information criterion, meant to estimate how well a model would predict new data, smaller being better) and one saying that the package's current value cannot rank these two fits, with the existing link to Lesson 4's WAIC subsection. This paragraph is in the PR #14 collision.
- Phase B Task 3 Minor 7 seam: "Rhat compares chains, so it flags this only when chains actually land in different explanations." appears in "When chains settle on two different explanations" and again in the appendix's mirror study; keep the main-text one and cut the appendix repeat.

- [ ] **Step 6: Coefficients ("How many effects does each fit resolve?" and "Put the true coefficient beside its estimate")**

- Report old 68 (heading presupposes "undetected"): done in Phase B (heading now "How many effects does each fit resolve?"); confirm.
- Report old 70: replace "after accounting for the other model components" with what they are (the other gradient and the hidden site factors), with a forward link to the ordination section.
- Report old 85: after the `summarize-effect-evidence` table, read it: give each fit's count of environmental intervals excluding zero inline from `effect_counts`, and say in one sentence that an interval excluding zero is read as a resolved direction because 95% of the posterior lies on one side.
- Report old 88: give the reason that perfect observation still leaves the probability uncertain: each site gives one presence or absence, not its probability.
- Report old 92 and 1598: replace "The public function" with "The package function" everywhere in prose (`/usr/bin/grep -n "public" vignettes/occJSDM-lesson-3.Rmd`), and cut the reviewer language "archived PCR fit", "provenance" and "not a replacement fit" from the paragraph before `native-read-examples`, leaving: these are the package's own plots with the simulated truth added; knitting shows saved images of them; the commands run on your own `fitmodel`.
- Report old 106: replace the two exporter sentences ("The exporter joins each summary ..." and "The environmental predictor scale was checked ...") with one reader-level sentence: truth is put on the same standardised scale as the estimates, so a cross and a bar can be compared directly.
- Report old 246 ("native"): at the first "native 95% posterior intervals", say once that "native" in these captions means drawn by the package's own plotting function; elsewhere in prose write "the package's".
- Phase B deferral (463 and 465): merge "the figure below keeps one species order for both fits" and "This figure is custom because no package plot draws both fits in one panel, which is the comparison this section's headline question asks for." into one sentence.
- Report old 132: after "These coefficients are changes in log-odds for a one-standard-deviation change", add that `runOccJSDM()` standardises the occupancy covariates itself (the fit stores the standardised design matrix), so you supply raw values and read effects per standard deviation; the response-curve section shows the conversion back to raw units.
- Report old 1669 seam: the paragraph after the custom figure repeats "If the interval crosses zero ..." from the paragraph after the package plots; keep the fuller one (after the custom figure, three-step reading) and cut the repeat.

- [ ] **Step 7: Response profiles and baseline ("What does an effect mean for a species' distribution?")**

- Report old 138: the hidden site-factor contribution is defined in the terms list; link to it.
- Report old 140: trim the list of three things a profile is not to one sentence with forward links (fitted probabilities at surveyed sites: "Distinguish fitted probabilities ..." section; new-site averages: [Lesson 4](occJSDM-lesson-4.md#which-true-probability-should-a-new-site-prediction-recover)).
- Report old 150 to 177 and 166: the custom curves and `outputs$gradients` were removed in Phase B and the raw-unit conversion now sits under the native curves; record both as done.
- Report old 181: compute "24,000 rows" inline from `lesson$manifests$default$mcmc` and say "four chains of 6,000 kept draws" from the same values; "posterior-draw-by-species matrix" now relies on the terms list.
- Report old 195: after "`colMeans()` averages probabilities after transforming every draw", say which to report: the mean of the transformed draws, as `colMeans()` gives, because it is the posterior mean of the probability itself.
- Report old 1673 seam: "This is not the average occupancy across the landscape" and "It is not average occupancy across sites." both appear; keep the first and cut the second sentence's repeat.
- Read the `baseline-probabilities` table in prose: name the species whose interval misses its truth, computed inline from `outputs$baseline` (`filter(arm == "default", truth < lower | truth > upper)`), or say none does.

- [ ] **Step 8: Traits ("Traits ask a harder, different question")**

- Report old 207: add an ecological example: a measured trait such as drought tolerance that makes some species respond less steeply to a drying gradient.
- Report old 218: say that trait coefficients are per standard deviation of the trait across the species surveyed, so a report should give the trait's standard deviation beside the coefficient.
- Report old 242: turn "Only one of the four intervals excludes zero" into advice: for a survey of ten or so species expect weak trait results; more species, not more PCRs or sites, is what helps; compute the count inline from `outputs$coefficients` (`block == "Trait"`, `excludes_zero`).
- Report old 248: cut "Knitting displays the exported figures from that unchanged fit, including all 24,000 retained draws."
- Report old 250 to 305 (show one native trait chunk): deferred; the five trait chunks are frozen and `native-traits-verify.R` requires all of them. Instead, before `native-traits-gradient-2`, add "The second gradient uses the same call with its name." and before the perfect-observation pair keep the existing one-line lead-in.
- Report old 255 (theme with the angle reset): after `native-traits-setup`, add one sentence: the setup takes ggtern's theme for the reason given at the start and resets the axis text angle that the package's trait plot otherwise tilts.
- Report old 336 and 1029 (unmeasured species traits, `n_lattrait`): in the appendix's "A real cancellation inside this simulated community", gloss the unmeasured species trait (a species characteristic the survey did not measure that also shapes environmental responses) and say it is set by `n_lattrait` in `listParams`.
- Report old 340 ("the contributions add"): add the reason in one clause: a least-squares regression is linear in its response, so regressing each component on the same traits gives slopes that sum to the slope of the total.
- Report old 343 to 355: say that `known_truth$jsdmParams_true$B` holds the true environmental coefficients with one row per covariate and one column per species, so `B[1, ]` is every species' response to gradient 1.

- [ ] **Step 9: Associations ("Residual species associations: did we recover what was put in?")**

- Report old 393: open with the reader's definition: a residual correlation measures whether two species occur together, or apart, more often than their measured environmental responses predict. Keep the existing caution that it is not a biological interaction.
- Report old 393 (the factor count): one sentence: with `n_factors` hidden factors each species is described by that many loadings, so the correlations can express at most that many independent patterns; too few factors force unrelated pairs to share one.
- Report old 395 to 400: say that the returned array is quantile by species by species, with the lower limit, median and upper limit in its first dimension.
- Report old 1759 (now the long paragraph before `native-correlation-example`): split it into a reading paragraph and a caution paragraph; compute "All 45 pairs" inline from `outputs$correlations` (the default arm's pairs with `species1 < species2` whose `lower` is below 0 and `upper` above 0, out of all such pairs); cut "Unlike the crosses in the coefficient plots, these Xs are uncertainty markers; the numbers supply the truth." which repeats the two sentences before it.
- Phase B Task 3 Minor 7 seam: OTU_4's `NA` is explained twice (the long paragraph and the paragraph after the image); keep the one after the image, which gives the reason, and cut "`NA` means undefined: `OTU_4` has zero true residual variance. This is not a true correlation of zero." from the long paragraph.
- Report old 428: after "Compare the pattern and magnitude", add what to do with real data, where there is no truth: read the X markers and the intervals, and report only pairs whose interval excludes zero as resolved in direction.

- [ ] **Step 10: Ordination ("Ordination: compare the combined effect before naming the axes")**

- Report old 430: open with a definition: an ordination places sites and species on a few axes that summarise the variation the measured covariates leave unexplained.
- Report old 432: split the long formula paragraph into four: what the hidden factors are in words; the formula, with each symbol named in words before it appears; what the product carries (residual co-occurrence); and rotation invariance with its consequence for the reader: axis signs and order can flip between fits, so never interpret an axis by its number or sign alone.
- Report old 470 to 478 (ordinary calls first) and the decision above: delete the `extract-ordination` chunk; move `ordination-standard`, its lead-in sentence, the paragraph "The ordinary biplot uses the package's stored orientation ...", the `ordination-ordinary-biplot-image` chunk and the paragraph about short arrows to directly after the split explanation, then the combined-contribution comparison.
- Report old 456: after "These are fitted-site results: the observations at a site helped estimate its hidden scores.", add the consequence: fitted sites look better recovered than new sites would.
- Report old section-wide (practical advice): add one paragraph on what an ordination of your own survey is for: map the site scores to look for an unmeasured gradient, and check which species load together before naming a cause.
- Phase B Task 3 Minor 7 seams: "A loading is a species' response to a unit change in a hidden site score ..." and "A site lying farther in a species-arrow direction ..." appear in the main section and again in the appendix; keep the main-section copies and cut the appendix repeats, keeping the appendix's sentences that are specific to the aligned plots.
- Report old 525 (every tenth site): in the appendix paragraph before `ordination-score-table`, say the ten sites were chosen by their order, before any estimate was inspected, so the choice cannot favour well-recovered sites.
- Report old 581 (medians cluster near zero): in the appendix after the site-score figure, add why: each site gives little information about its hidden scores, so the estimates shrink towards zero; on a real ordination, expect site scores to look less spread out than the conditions they stand for. Shorten the circle-radius sentences to one: the circles are rough size guides built from the marginal intervals, not credible regions.
- Report old 651: the biplot caution paragraph in the appendix is split in Task 2.
- Variation partitioning (report old 653 to 686): the section was removed in Phase B by Doug's decision, with a `TODO.md` item to restore it; record all its items as not applicable.

- [ ] **Step 11: Collection effects and detection effort**

- Report old 689: give concrete collection covariates (water volume filtered, turbidity, time since sampling) so the reader can map them to their fieldwork.
- Report old 699: after the collection plot, say that it shows the slope; the intercept is the first row of the array returned by `returnCollectionCovariates()` below.
- Report old 721: split the paragraph that begins "Two related helpers answer different questions about collection." and cut its last sentence ("Its truth accounts for whether simulated reads actually pass the fitted threshold.").
- Report old 1717 (the read-threshold adjustment explained five times): keep the full explanation once, in "Laboratory true-positive and false-positive rates by primer"; elsewhere (the primer trace paragraph, "Separate false-positive and detection-rate plots", the detection-effort paragraph, the probability table section) replace it with a short reference to that subsection.
- Report old 739: after the paragraph on orange and native bars, say which to use when planning a survey: the expected count (the formula) to compare designs, and the survey-outcome interval to see how much one survey can vary around it.
- Report old 758: add the take-home with numbers computed inline from `outputs$detection_effort`: the true expected number of species detected with one field sample and six PCRs per primer (`M == 1, K == 6`) against two field samples with one PCR each (`M == 2, K == 1`), and the sentence that, in this community, a second field sample buys more than more PCRs.
- Slash pair "collection/PCR" in the native-bars paragraph: Task 2.

- [ ] **Step 12: Fitted probabilities, states and the latent-presence table**

- Report old 762: lead "Distinguish fitted probabilities, occupancy states and new-site predictions" with what `computePredictiveOccupancyProbs()` is for (the fitted habitat-based probability at surveyed sites, which Ji et al. reported), then the naming caution.
- Report old 766 to 774: add one line to `extract-conditional-occupancy` (not frozen): `predictive_occupancy <- occJSDM::computePredictiveOccupancyProbs(fitmodel)`, with a comment that it has the same site-by-species shape.
- Report old 1107: cut "No estimates are recalculated or replaced with truth."; replace the slash pair "sample/primer" with "sample and primer".
- Report old 1140 to 1159 (the site-build branch): split `display-native-state-table` into a visible chunk that builds and prints the gt table and runs only for HTML, and a hidden chunk that prints the plain table for Markdown; do the same for `display-native-probability-table`. The first becomes:

```{r display-native-state-table, eval=knitr::pandoc_to() %in% c("html", "html4", "html5")}
state_columns <- c(
  "Site", "Sample", "Primer", "PCR", "OTU", "Source",
  "TrueSite", "CondOccProb", "TrueSample", "CondSampleProb"
)

occJSDM::plotLatentPresences(
  state_table,
  species_name = "OTU_1",
  title = "Actual states beside the model's probabilities",
  columns = state_columns,
  container_height = 450
) |>
  gt::fmt_number(columns = c(PCR, TrueSite, TrueSample), decimals = 0) |>
  gt::fmt_percent(columns = c(CondOccProb, CondSampleProb), decimals = 1)
```

```{r display-native-state-table-markdown, echo=FALSE, eval=!(knitr::pandoc_to() %in% c("html", "html4", "html5"))}
state_columns <- c(
  "Site", "Sample", "Primer", "PCR", "OTU", "Source",
  "TrueSite", "CondOccProb", "TrueSample", "CondSampleProb"
)

state_table |>
  select(all_of(state_columns)) |>
  knitr::kable(digits = 3)
```

and the probability pair follows the same pattern with `probability_columns`, `probability_table`, `container_height = 350` and its `fmt_percent` line unchanged. Replace the paragraph after the first table with: `plotLatentPresences()` returns a gt table, which displays in HTML; the Markdown version of this page shows the same records as a plain table with decimals. Check after rendering that the `.md` tables are identical to the before-state tables.
- Report old 1161: after "At site ..., sample ... has no OTU_1 DNA", add that sample identifiers run across the whole survey, which is why a site's samples are not numbered 1 to 3.
- Report old 1167 to 1201 (`echo=FALSE` for the truth join): deferred by the decision "Visible code stays visible"; shorten the prose before `join-native-probability-truth` to what the join attaches.
- Report old 1130 and 1088: kept as they are (the report praises them).

- [ ] **Step 13: Finder, closing, references and appendix**

- Report old 2057: cut the last sentence of the appendix ("Student-facing figures use tidy tables ...").
- Report old 1529: cut "also used by our offline verification scripts" from "Where to find other parameter draws".
- Report old 1923 to 1925 (PR #13 history): in the appendix by spec; shorten "The corrected response-curve helper" to its practical content: `plotCovariateEffect()` works in original units and its median column is named `median`; the history sentence stays as one sentence with the PR link.
- Report old 1446 and 1448 (mirror study length): the main text has the short warning and check (Step 5); in the appendix, split the two long paragraphs in Task 2.
- Phase B Task 3 Minor 3 (finder texts): fix "Truth check: [Lesson 4](occJSDM-lesson-4.md): 300 ..." to a single colon ("Truth check: 300 independent non-spatial sites in [Lesson 4](occJSDM-lesson-4.md) ..."), and in the "How does each species respond" line say "the coefficient plots above, for both fits, and the response curves above, for the PCR fit".
- "Where to go next": keep its content; check it still names every appendix subsection after this pass's moves.
- Factual check, report old 780 (Ji et al. Supplementary Information 12): keep the sentence; Task 9 lists it for Doug, a coauthor, to confirm.

- [ ] **Step 14: The ggtern TODO item**

Add this bullet to `TODO.md` under "## Review and maintenance", directly after the "**Per-chain draws accessor ...**" bullet (re-read the file first):

```markdown
- **The ggtern theme workaround in the lessons (added 3 October 2026, lesson rewrite Phase C; needs Alex):** after `loadNamespace("occJSDM")`, which loads ggtern, `theme_set(ggplot2::theme_bw())` makes the next ggplot fail with the error `The tern.axis.ticks.length.major theme element must be a <rel> object.` (ggplot2 4.0.3, ggtern 4.0.0); adding `theme_bw()` to one plot and `theme_set(ggtern::theme_bw())` both work, so Lessons 3 to 6 use ggtern's theme and Lesson 3 explains why. Confirm the cause and whether occJSDM can avoid it, for example by calling ggtern only inside the ternary plotting function, so that users can keep `ggplot2::theme_bw()`.
```

- [ ] **Step 15: Check every spec bullet and report item for this lesson**

Go through the spec's Lesson 3 bullets and section 2 of the report (old lines 27 to 785 and 1086 to 2085) once more against the edited file. For each, confirm the edit or note why it is deferred; keep the list for the commit message body and Task 9. Spec bullets already done in Phase B and only to be confirmed: restructure, diagnostics first, package plots with truth, appendix moves, traces through the package route, uncited references removed, Lesson 2 statements fixed, `sampleresults` paragraph cut, pipe tables converted. The 1,150-line target is replaced by `notes.md` (no target).

- [ ] **Step 16: Render, verify, check, commit**

```sh
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/occJSDM-lesson-3.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/occJSDM-lesson-3.Rmd", output_format="rmarkdown::html_vignette")'
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/revert-unchanged-png.zsh vignettes/occJSDM-lesson-3_files
git status --short vignettes/
S=/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/run-step.zsh; A=/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
zsh $S c11-L3-teach-native-plots Rscript dev/simstudy/vignette-lesson/verify_native_plots.R $A
zsh $S c12-L3-teach-native-traits Rscript dev/simstudy/vignette-lesson/native-traits-verify.R $A
zsh $S c13-L3-teach-ordination Rscript dev/simstudy/vignette-lesson/ordination-verify.R $A
zsh $S c14-L3-teach-remaining-plots Rscript dev/simstudy/vignette-lesson/remaining-plots-verify.R $A
```

Expected: both renders clean; the revert script reverts every PNG (no chunk that draws a figure changed; `mirror-labelling-chains-1` and `trait-cancellation-1` may need `git checkout --`); each verifier ends `EXIT 0` and its log contains "Lesson 3 displays the exact exported ...". Read the rendered `.md` end to end once, and diff it against `lesson-3-before.md` to confirm that every computed table is unchanged except the new `show-teaching-objects` output. Run the Global Constraints checks with `<files>` = `vignettes/occJSDM-lesson-3.Rmd TODO.md`. Then:

```bash
git add vignettes/occJSDM-lesson-3.Rmd vignettes/occJSDM-lesson-3.md vignettes/occJSDM-lesson-3_files TODO.md
git commit -m "Lesson 3 teaching pass: name each idea before use, read every output, add practical advice

Terms list, teaching objects shown, ggtern theme explained, diagnostics overview
first, numbers computed from the bundles, repeated passages cut, latent tables
shown without the site-build branch. Deferred: one native trait chunk (frozen
for its verifier) and hiding the truth join (all teaching code stays visible).

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Lesson 3, copyedit commit

**Files:**
- Modify: `vignettes/occJSDM-lesson-3.Rmd`, `vignettes/occJSDM-lesson-3.md`
- Read: `style-diagnostic/lesson-3.md` section 5

**Verifiers that read this lesson:** the four of Task 1.

**Interfaces:**
- Consumes: Task 1's commit.
- Produces: the finished Lesson 3; its headings and terms list are what Tasks 3, 4 and 6 point to.

- [ ] **Step 1: Copyedit the whole lesson**

Apply the copyedit procedure (Shared text). The report's section 5 items, mapped to the current text: long paragraphs, split each into two or more (the mirror study's two appendix paragraphs that begin "One of our simulation studies" and "Rhat compares chains"; the paragraph that begins "In the study, the counts differed"; the ordination formula paragraph if Task 1 left any part over 600 characters; the residual-correlation reading paragraph; the Ji et al. paragraph "Which one should a study report?"; the collection-helpers paragraph; the `coda` and `posterior` paragraph; the appendix biplot caution paragraph that begins "A site lying farther"); long sentences ("That product is what carries residual co-occurrence", "`computeAverageCollectionProbs()` instead returns", "The two are most useful together", "For one species in that 300-site community", "The study ran 8 chains"); slash pairs ("before/after" if any remains, "collection/PCR", "simulation/fitting" in "Reproduce the teaching figures", and any other the prose checks or a read finds); the asides listed there that Task 1 has not already cut. Keep "presence/absence".

- [ ] **Step 2: Render, verify, check, commit**

Render both outputs as in Task 1 Step 16, apply the PNG policy, and run the four verifiers with step names `c15-L3-copy-native-plots`, `c16-L3-copy-native-traits`, `c17-L3-copy-ordination`, `c18-L3-copy-remaining-plots`. Expected as in Task 1 Step 16. Run the Global Constraints checks with `<files>` = `vignettes/occJSDM-lesson-3.Rmd`; every sentence over 30 words that remains has a stated reason. Report the line count (`wc -l vignettes/occJSDM-lesson-3.Rmd`) before Task 1 and now. Then:

```bash
git add vignettes/occJSDM-lesson-3.Rmd vignettes/occJSDM-lesson-3.md vignettes/occJSDM-lesson-3_files
git commit -m "Lesson 3 copyedit pass

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Lesson 4, teaching commit and copyedit commit

**Files:**
- Modify: `vignettes/occJSDM-lesson-4.Rmd` (362 lines), `vignettes/occJSDM-lesson-4.md`; `vignettes/occJSDM-lesson-4_files/` only for a figure whose chunk changes on purpose (none expected)
- Create (not committed): `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/lesson-4-before.md`
- Read: `style-diagnostic/lesson-3.md` items for old lines 786 to 1085 and the matching section 5 items; the finished Lesson 3 (terms list, diagnostics section, ggtern sentence)

**Verifiers that read this lesson:** `prediction-verify.R` (it parses and executes the frozen chunks `prediction-fit-one-factor` and `prediction-native-call` with stubbed API calls and compares the captured arguments and seeds).

**Interfaces:**
- Consumes: `teaching-data/nonspatial-lesson.rds` and `teaching-data/prediction-lesson.rds` (`prediction_examples`: `input` with `raw_covariates`, `outside_training_range`, `environmental_eta`, `selected_sites`, `public_prediction_seed`, `fitting_seed`; `manifests$two_factors` and `$one_factor` with `mcmc`, `priors`, `waic`, `factors`; `truth`, `cells`, `scores`, `paired_site_differences`, `public`, `diagnostics`, `probability_diagnostics`); Lesson 3's anchors `#what-this-lesson-answers` and `#check-computation-as-well-as-ecological-recovery`.
- Produces: the kept headings "Predict occupancy at genuinely new sites", "Which true probability should a new-site prediction recover?", "Compare models using what actually occurred" and "Check the additional fit and understand the WAIC limitation", which Lesson 5 (Task 4) links.

- [ ] **Step 1: Read, then keep the before state**

Read the report items for old lines 786 to 1085, then the lesson end to end. Run:

```sh
git show 7662f9a:vignettes/occJSDM-lesson-4.md > /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/lesson-4-before.md
```

- [ ] **Step 2: Opening, prerequisites, theme**

- In "What this lesson answers", insert the prerequisites block after the "All teaching code is visible" paragraph; Lesson 4's verbs are `filter()`, `mutate()`, `select()`, `group_by()`, `summarise()`, `left_join()`, `if_else()` and `pivot_wider()`, with `pivot_wider()` as the unusual operation, explained where it pairs the two models' scores.
- Add the ggtern pointer (Shared text) after the `load-results` chunk.
- Point to Lesson 3's terms list for posterior draws, chains and credible intervals in one clause of the "In the figures" paragraph.

- [ ] **Step 3: "Predict occupancy at genuinely new sites"**

- Report old 790 (choosing `n_factors`): after "The generating community has two factors, but that does not guarantee ...", add the practical advice: there is no rule yet for your own data; start small, compare candidate counts by held-out scores as this lesson does, and check that conclusions do not change, with a link to [Lesson 2's fitting reference](occJSDM-lesson-2.md#fitting-your-own-data-what-the-call-needs). The generating number need not win.
- Report old 806: after "The fitted model standardizes them using the **training** means and standard deviations.", say that `predictNewSites()` does this itself when you supply raw values with the training column names, so you never standardise new sites yourself. Add the source of "standard deviation 10" in one clause (the new sites were drawn as in `dev/simstudy/vignette-lesson/prediction-build.R`, `rnorm(..., sd = 10)`), or compute it inline as the rounded standard deviation of `prediction_examples$input$raw_covariates` if that reads better; keep "mean zero".
- Read the `prediction-load` output in prose: what the six rows show (raw values on the simulation's scale, some far from zero), and the count of sites beyond the training range, which is already inline.

- [ ] **Step 4: "Which true probability should a new-site prediction recover?"**

- Report old 829: before the `prediction-one-site-meaning` chunk, say that `lesson$input$jsdm$sigma_h` is the standard deviation of each hidden site factor in the simulation (`r lesson$input$jsdm$sigma_h` here), and that multiplying by the length of the species' loadings gives the spread of its hidden contribution.
- Report section 5 (hard-coded 93% and 82%): compute both inline from `site_example$conditional_truth` and `average_over_hidden_conditions` with `scales::percent(..., accuracy = 1)`.
- Report section 5 (curly quotes): replace the curly quotes around "convert the average log-odds" and "average the converted probabilities" with straight ones.

- [ ] **Step 5: "Use the package's new-site prediction function"**

The chunk `prediction-native-call` is frozen; the edits are prose.
- Report old 871 and 891 (seed before its reason): move "For each retained parameter draw, `predictNewSites()` also draws new hidden conditions." to the paragraph before the chunk and add: so set a seed first if you want the same quantiles each time.
- Report old 876: after the chunk, explain the arguments in one sentence each: `useSpatial = FALSE` leaves out the spatial field, which this non-spatial fit does not have; `confidence = 0.95` sets the interval's coverage; `verbose = FALSE` silences progress messages. Add that `useBiotic`, left at its default, includes the hidden-factor term for a fit with factors, drawing new hidden conditions as described.
- Read the `prediction-native-intervals` figure: say how many of the twenty intervals contain their cross, computed inline from `native_prediction_examples` (`conditional_truth` between `lower` and `upper`), and keep the existing caution that twenty examples cannot establish coverage.

- [ ] **Step 6: "Check point predictions against the appropriate truth"**

- Report old 921: replace "The exporter performs ..." and "The development helper is separate from the package's public `predictNewSites()` function." with what the reader can and cannot do: the package currently returns only quantiles (`summarised = FALSE` stops with "Only summarised version for now", `R/output.R` line 1667), so with your own data report the median and interval; the posterior mean scored here came from a development script that averages every draw's probability over hidden conditions, which the package cannot yet return. Compute 24,000 inline from `prediction_examples$manifests$two_factors$mcmc`.
- Read the `prediction-marginal-recovery` figure: say where the points fall relative to the diagonal and whether the crosses (sites beyond the training range) behave differently, computed inline from `prediction_cells` (mean absolute error for `outside_training_range` TRUE and FALSE).
- Report section 5 aside: cut "that last sentence is only an illustration of the unit, not a claim that all errors equal ten points" and keep the illustration.

- [ ] **Step 7: "Compare models using what actually occurred"**

- After the Brier and log-score definitions, add a reference point computed inline in a small run chunk: the Brier score of the true marginal probabilities themselves against the same outcomes (`mean((truth - z)^2)` over one model's rows of `prediction_cells`), so the reader sees how far from zero even perfect probabilities score.
- Report old 993 to 1012: reduce the two paragraphs after `prediction-paired-scores` to the teaching point in two or three sentences: both models predict the same sites, so compare them in pairs averaged within sites; the difference (inline) is tiny next to the numerical and survey-to-survey uncertainty, so it cannot rank the models. Move the rest (the site-based standard error's caveats, the `prediction-verify.R --score-mcse` calculation, the 0.00026 Monte Carlo standard error and the 1.3 ratio) to a new paragraph in the appendix headed by one sentence saying it records the numerical check; 0.00026 is not in the bundle, so keep its source (the build log of `prediction-verify.R --score-mcse`). Keep the ratio computed inline from the bundle's difference and that stated value.
- Read the `prediction-observed-scores` table in prose with both scores inline from `observed_scores`.

- [ ] **Step 8: "Check the additional fit and understand the WAIC limitation"**

This subsection is in the PR #14 collision; keep edits to the sentences named here.
- Report old 1029: after `prediction-fit-one-factor` (frozen), add one sentence: `n_lattrait = 1` keeps the one unmeasured species trait the baseline fit used (Lesson 3's appendix explains it), so only `n_factors` differs.
- Report old 1030 to 1031 (write the MCMC settings literally): deferred, because `prediction-verify.R` executes the chunk and compares its arguments; the next paragraph already states the settings. Compute them inline there from `prediction_examples$manifests$two_factors$mcmc` (four chains, 3,000 burn-in, 6,000 kept, no thinning).
- Report old 1036 to 1064: the screen sentence already points to Lesson 3's diagnostics; confirm and keep.
- Read the `prediction-diagnostics` table: give the largest prediction Rhat and smallest prediction ESS inline.
- Report old 1066 and Phase B Task 2 M2: replace "The old walkthrough extracted WAIC to compare model specifications. Here is the current extraction syntax, followed by the values from these same-data fits:" with a definition and lead-in: WAIC, the widely applicable information criterion, is meant to estimate how well a model would predict new observations, with smaller values better; here is the call and the values for these two fits.
- Report old 1083: before the bold sentence's explanation, add the plain version: the current WAIC rewards fitting the training survey's own hidden states, not predicting new sites. Keep the bold advice.

- [ ] **Step 9: Closing and appendix**

- "Where to go next": keep; make sure it says what the appendix now holds (the rebuild commands and the numerical check of the score difference).
- Appendix: the reproduction paragraph stays; add the moved paragraph from Step 7.

- [ ] **Step 10: Render, verify, check, commit the teaching pass**

```sh
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/occJSDM-lesson-4.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/occJSDM-lesson-4.Rmd", output_format="rmarkdown::html_vignette")'
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/revert-unchanged-png.zsh vignettes/occJSDM-lesson-4_files
git status --short vignettes/
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/run-step.zsh c21-L4-teach-prediction Rscript dev/simstudy/vignette-lesson/prediction-verify.R /Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples /Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples/prediction
```

Expected: renders clean; no PNG kept; the verifier ends "PASS: displayed fit/prediction chunk arguments and seeds captured without fitting" and "All prediction verification checks passed." with `EXIT 0`. Diff the `.md` against `lesson-4-before.md`: computed tables unchanged apart from the new reference-score output. Run the Global Constraints checks with `<files>` = `vignettes/occJSDM-lesson-4.Rmd`. Then:

```bash
git add vignettes/occJSDM-lesson-4.Rmd vignettes/occJSDM-lesson-4.md vignettes/occJSDM-lesson-4_files
git commit -m "Lesson 4 teaching pass: explain each argument, read each result, define WAIC

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 11: Copyedit, render, verify, check, commit**

Apply the copyedit procedure. Report section 5 items for this lesson: the long sentences the prose checks list (the opening's "This lesson supplies 300 ..." and "It continues Lesson 3 ...", the bold WAIC paragraph), any slash pair other than "presence/absence". Render both outputs, apply the PNG policy, rerun `prediction-verify.R` with step name `c22-L4-copy-prediction`, run the Global Constraints checks, and commit:

```bash
git add vignettes/occJSDM-lesson-4.Rmd vignettes/occJSDM-lesson-4.md vignettes/occJSDM-lesson-4_files
git commit -m "Lesson 4 copyedit pass

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Lesson 5, teaching commit

**Files:**
- Modify: `vignettes/occJSDM-lesson-5.Rmd` (848 lines), `vignettes/occJSDM-lesson-5.md`, `TODO.md`; `vignettes/teaching-data/lesson-5-*.png` only for a figure whose chunk changes on purpose (none expected)
- Create (not committed): `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/lesson-5-before.md`
- Read: `style-diagnostic/lesson-4.md` items for old lines 34 to 855 (opening to the reproduction record), the spec's Lesson 4 bullets, `notes.md` (sjSDM decision), `dev/simstudy/jsdm-package-comparison/SJSDM-STABILITY-REPORT.md` ("Two genuine local maxima" and "Twelve independent starts"), the finished Lessons 3 and 4

**Verifiers that read this lesson:** `dev/simstudy/jsdm-package-comparison/verify-teaching.R` (it evaluates the frozen `reproduce-community` chunk and checks the bundle).

**Interfaces:**
- Consumes: `teaching-data/jsdm-comparison.rds` (`comparison`: `training` with `x` and `y`, `test_x`, `truth` with `parameters`, `conditional_probability`, `scaled_coefficients`, `loadings`; `predictions` with `package`, `fit`, `target`, `site`, `species`, `estimate`, `truth`, `observed`, `band` whose levels are "Below 20%", "20% to below 80%" and "80% or above"; `curves`, `point_parameters$gllvm` and `$sjSDM` with `beta` (intercept and two slopes by species) and `loading` (two factors by species), `attempts` with `package`, `fit`, `seconds`, `diagnostics`, `selected`, `sjsdm_revision`); Lesson 4's anchors from Task 3; Lesson 3's anchors `#what-this-lesson-answers`, `#check-computation-as-well-as-ecological-recovery` and `#when-chains-settle-on-two-different-explanations`.
- Produces: the TODO item for the v0.1.0 rerun; the kept headings listed in the Global Constraints.

- [ ] **Step 1: Read, then keep the before state**

Read the report items, the spec's Lesson 4 bullets and the stability report sections named above, then the lesson end to end. Confirm from the stability report that it does not establish that the two maxima explain the unpenalised starts' disagreement (its "Two genuine local maxima" section reports the maxima under the penalty only); the Phase B hedge stays. Run:

```sh
git show 7662f9a:vignettes/occJSDM-lesson-5.md > /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/lesson-5-before.md
```

- [ ] **Step 2: Opening ("What are we comparing?")**

- Report old 38 (why compare with packages that have no detection model): add that only occJSDM models detection error; this lesson checks that its ecological core, the part it shares with established JSDMs, does as well as theirs on perfectly observed data.
- Report old 38 ("pure JSDM portion"): tie it to the quickstart's note that `info` with one row per site makes `runOccJSDM()` skip the detection stages and fit a JSDM to the observed presence/absence; that is the mode used here.
- Report old 38 (packages never introduced): one clause each: gllvm, generalised linear latent variable models fitted by approximate likelihood; sjSDM, a fast JSDM fitted by optimisation in PyTorch; Hmsc, a Bayesian hierarchical JSDM widely used in community ecology.
- Report old 40 and Phase B hedge: shorten the sjSDM revision history in the opening to one sentence: section 5 explains why the sjSDM fit needed several starts and what its two solutions mean. The hedged account stays in section 5.
- Phase B deferral: delete "You do not need the spatial lesson first."
- Phase B Task 5 M4: objective 6 becomes "Use what makes the model joint to predict one species from another."
- Insert the prerequisites block after the "All teaching code is visible" paragraph; Lesson 5's verbs are `filter()`, `select()`, `mutate()`, `transmute()`, `group_by()`, `summarise()`, `distinct()`, `pull()`, `left_join()`, `inner_join()`, `as_tibble()` and `pivot_longer()`; name base R's `sweep()` and `integrate()` as the unusual operations explained where they appear. Say that knitting needs only the packages loaded here; gllvm, sjSDM and Hmsc are needed only for the optional fits in the appendix.
- Report old 54 to 89: after `load-comparison`, add a run chunk `str(comparison, max.level = 1)` and read it in a bullet list naming each element the lesson uses.
- Report old 88: add the ggtern pointer (Shared text) after the chunk.

- [ ] **Step 3: Section 1 ("What the models receive, and what we keep secret")**

- Report old 104 to 106: before the tile figure, name `known_truth`'s parts the lesson uses: `conditional_probability`, `parameters`, `scaled_coefficients` and `loadings`.
- Report old 138: in "Optional: reproduce this community", say what the logit scale is (log-odds, as in Lesson 3's terms list) and that a slope of 1 multiplies the odds by about 2.7 per unit.
- Spec "covariate standardisation" and report old 202: after "Test sites use the **same** transformation", add the practical point for the reader's own data: `runOccJSDM()` standardises the occupancy covariates itself and `predictNewSites()` applies the training transformation to raw new-site values (Lesson 4), so with occJSDM you supply raw values; this lesson standardised beforehand only so that all four packages receive identical inputs.

- [ ] **Step 4: Section 2 ("How similar are the four models?") and the sjSDM version**

- Report old 206 (why the factor count is fixed): say it removes one source of difference between packages; on a real survey the count is unknown, and Lesson 4 shows how to compare candidate counts by held-out scores (link `occJSDM-lesson-4.md#compare-models-using-what-actually-occurred` if that heading survived Task 3, else the lesson). Fix "All four receive the same information, linear environmental predictors and two hidden factors" (report section 5): the models receive the predictors and estimate two hidden factors.
- Report old 208 to 213: before the bullets, one sentence splitting the packages into two families: Bayesian sampling (occJSDM, Hmsc) returns a posterior with uncertainty; optimisation (gllvm, sjSDM) returns one best estimate. Leave the bullets' technical terms for the appendix, glossing "variational approximation" in a clause.
- Report old 215: after the probit paragraph, point to Lesson 6's coefficient-coverage section for what the probit link means for coefficient checks.
- Spec sjSDM bullet and `notes.md`: replace the paragraph "The sjSDM R version is 1.0.7 inside Doug's fork release **v0.2.1**. ..." with:

```markdown
To fit sjSDM yourself, install it from Doug's fork, release [v0.1.0](https://github.com/dougwyu/s-jSDM/releases/tag/v0.1.0), which runs on Apple Silicon and uses PyTorch only. The fits in this lesson ran under the fork's later release v0.2.1 (sjSDM 1.0.7) with its optional Mojo backend switched off. We did not refit under v0.1.0. Our reason is an argument from reading the two releases' code, not a rerun: with the backend switch at zero, v0.2.1's fitting code reaches the same PyTorch loss as v0.1.0 line for line, apart from integer conversions of tensor sizes. A rerun under v0.1.0 is planned. We did not test the CRAN release of sjSDM.
```

- [ ] **Step 5: Sections 3 and 4 (the reminder and the two probability questions)**

- Phase B Task 5 M3: in section 3's reminder, link the build-up to [Lesson 1's joint-model section](occJSDM-lesson-1.md#what-a-joint-model-does) and keep the link to `#what-a-joint-model-adds-to-a-factorisation` for the additions; add that the exercise at the end of section 4 shows the first consequence in numbers.
- Spec "Cross-reference Lesson 3's marginal-versus-conditional and WAIC teaching" (now Lesson 4) and report old 296: at the start of section 4, one sentence reminding the reader that [Lesson 4](occJSDM-lesson-4.md#which-true-probability-should-a-new-site-prediction-recover) met these two targets as the conditional and the marginal probability; use "conditional" for the sampled-site question in the first bullet.
- Report old 302 to 324: after `averaging-illustration`, read its output with both numbers inline from `zero_factor_probability` and `marginal_probability`: averaging nearly doubles the probability here, because at a low score favourable hidden conditions raise the probability more than unfavourable ones lower it; the hidden spread of 2 was chosen large to make the effect visible.
- Report old 326: replace "gllvm's level-zero prediction ..." and "The inspected sjSDM environmental prediction does too." with the practical point: each package's ordinary prediction call sets the hidden scores to zero, so neither gives the marginal probability directly; this lesson computed it, as below.
- Report old 328: reduce the three integration methods to one sentence (each package's prediction was converted to the same marginal question) and move the rest of that paragraph to the appendix, before "What did these runs cost?", under the sentence "How each package's new-site prediction was made marginal:".
- Report old 366: after `marginal-extraction-example`, read the result with both numbers inline from `fitted_probability` and `true_probability` (estimated against true, for species_01 at the mean environment).
- Spec (why gllvm, not occJSDM) and report old 235 and 332 to 364: after the extraction paragraph, explain: occJSDM's `predictNewSites()` answers the marginal new-site question by drawing new hidden conditions (its `useBiotic` term, on by default for a fit with factors; Lesson 4 shows how to read its interval), but it has no mode that conditions on another species' record at the new site, so the exercise below uses gllvm's point parameters. Add this not-run chunk:

```{r occjsdm-new-site-route, eval=FALSE}
# fit_occJSDM is the fit from the appendix's occJSDM call.
# predictNewSites() draws new hidden conditions, so set a seed for repeatable quantiles.
set.seed(26092212)

test_quantiles <- occJSDM::predictNewSites(
  fit_occJSDM,
  X_psi = comparison$test_x,
  useSpatial = FALSE,
  useBiotic = TRUE,
  confidence = 0.95,
  verbose = FALSE
)

# Quantile by site by species: lower limit, median, upper limit. Species 1 is species_01.
test_quantiles[, 1:5, 1]
```

- Report old 258: in the exercise, say once that `beta` has one row for the intercept and one per slope, and `loading` one row per factor, one column per species each.
- Report old 283 and report section 5 ("somewhat smaller"): after `conditional-prediction`, read the table with numbers inline from `conditional_table`: the fitted gap and the true gap at the mean environment and at +1 SD, said as a fraction; and say why the baseline (species_10 alone) is also off: the fitted intercept and slopes differ from the generating ones, as section 8's curves show.
- Report old 285: add the practical advice: occJSDM has no conditional-prediction call, so for your own fitted occJSDM model this calculation, applied within each posterior draw, is the way to get it; Lesson 3's appendix shows where the draws are.

- [ ] **Step 6: Section 5 ("Check the computation before interpreting the ecology")**

- Report old 370: compute the chain settings inline from the occJSDM call they came from or state them from the appendix's `fit-occJSDM` chunk; say why they differ from the quickstart's (a longer run so that every monitored quantity passes the stricter screen), and that warm-up is the quickstart's burn-in, set by `nburn`.
- Report old 389: name the 1.01 and 400 screens as the published recommendations Lesson 3's diagnostics section explains, with its link; replace "bulk/tail ESS" with "bulk and tail ESS"; compute "1,135 monitored quantities" inline from `nrow(comparison$diagnostics)` after checking that it is the count the sentence means (all rows for occJSDM and Hmsc), otherwise from the `diagnostic_summary` sum.
- Report old 391: gloss EVA (gllvm's extended variational approximation) and "inconsistent objectives" (starts that reached different fitted log-likelihoods), and add the practical check for the reader's own gllvm fit: run several starts and compare their log-likelihoods.
- Report old 393: give the reason for the 0.1 stability check in one clause (a difference smaller than this in log-likelihood is too small to change predictions noticeably, which is why it was declared before fitting) only if `dev/simstudy/jsdm-package-comparison/DESIGN.md` or `FITTING-REPORT.md` records that reason; otherwise say it was declared before fitting as the stability criterion, with the file that records it.
- Report old 397: lead the two-optima subsection with the hills image, then the technical words.
- Report old 414: say why the selection rule was declared before truth was read: it prevents choosing the answer you like.
- Report old 450: add the transferable lesson: MCMC chains can also settle in different modes, which is why several chains are run and compared, with a link to [Lesson 3](occJSDM-lesson-3.md#when-chains-settle-on-two-different-explanations).
- Spec line 765 bullet: confirmed in Step 1 that the stability archive does not establish the causal link; keep the hedge in section 5 and the appendix exactly as Phase B wrote it.

- [ ] **Step 7: Sections 6 to 9 (results)**

- Report old 458 to 494: after the two probability figures, say that the four panels in each look nearly identical, so the figures mainly test the shared model and the data, not the packages.
- Report old 496 and 540 (repetition): the explanation of why sampled-site scatter is larger appears after the sampled-site figure and again in section 7; keep the fuller one in section 7 and reduce the first to one sentence pointing to it.
- Report section 5 (hard-coded 6.1, 7.1, 7.1, 6.3 and 12.1, 11.9, 13.3, 12.2): compute them inline from `overall_errors`, in package order; compute "about +1 point" and the "6.1 points" example in the next paragraph inline too.
- Report old 538: cut "These are measured results, not notional examples, and they are not the older occJSDM-only sample-size experiment."
- Report old 546: give the reason for the 20% and 80% band edges: they separate rare-at-this-site, middling and near-certain occurrences, the same bands Lesson 2 used for its pull towards the middle.
- Spec "Read the printed results in prose: ... the shrinkage towards the middle in the band tables" and report old 575: after the band tables, read them with numbers inline from `band_errors`: at sampled sites every package overestimates the "Below 20%" band and underestimates the "80% or above" band (for example occJSDM's signed errors `filter(band_errors, package == "occJSDM", question == "Reconstruct sampled sites")`), which is shrinkage towards the middle, explained as in [Lesson 2](occJSDM-lesson-2.md#why-rare-and-common-species-are-pulled-towards-the-middle); at new sites say which packages err in which direction in the low and high bands, from the same table, and that this one community cannot say why they differ.
- Report old 603 to 613: after `inspect-one-species`, say whether species_01's signed errors are negative while the overall averages are positive, inline from `species_errors`.
- Report section 5 ("a distance from matching simulated truth"): reword to "a distance from the simulated truth".
- Report old 620: after "no common uncertainty interval was calculated across all four fitting methods", point to Lesson 6's interval-coverage section for why.
- Spec "the four packages' identical response curves" and report old 666: after the gradient figures, say that in almost every panel the four curves lie on top of each other, so the misses (name the species whose curves visibly miss, checking the rendered figures, as the report's species_04 and species_09 on gradient 1 and species_03 on gradient 2) are shared by all four packages: they come from what 100 sites can reveal, not from one package. Cut the aside "they illustrate how to read a curve, rather than determining which results we include".
- Report old 686 to 689: after `held-out-outcome-scores`, add this run chunk and read it:

```{r reference-scores}
training_prevalence <- colMeans(training$y)

new_site_predictions |>
  filter(package == "occJSDM") |>
  mutate(prevalence = training_prevalence[as.character(species)]) |>
  summarise(
    `True probabilities` = mean((truth - observed)^2),
    `Each species' training prevalence` = mean((prevalence - observed)^2)
  ) |>
  knitr::kable(digits = 4, caption = "Brier scores of two references on the same 3,000 outcomes")
```

The first is the best a model could expect to score here; the second is what ignoring the environment scores; place the four packages between them, inline.
- Report old 689: add how to get such scores on your own survey: hold out some sites before fitting, predict them with `predictNewSites()` as Lesson 4 does, and score the predictions against what was observed there.

- [ ] **Step 8: Section 10 (closing)**

- Report old 837 to 839: add a take-home for the eDNA user: on perfectly observed presence/absence data, occJSDM's ecological model performed like the established packages in this community, so the choice among them can rest on whether you need the detection model, which only occJSDM has.
- Leave the "Quickstart and lesson guide" label for Task 8.

- [ ] **Step 9: Appendix (fitting calls, costs, evidence, record)**

- Report old 693: say which versions to install, pointing to the bullets of section 2 and the reproduction record.
- Report old 697: in the occJSDM subsection, say that `OTU = training$y` holds 0 and 1, so no read threshold applies, and gloss `n_lattrait = 0` (no unmeasured species trait).
- Spec ("what the section 10 arguments do, one sentence each, pointing to the dev scripts for the rest") and report old 714 to 730: after `fit-gllvm`, one sentence each for the residual-based start (`starting.val = "res"`), `jitter.var`, the `maxit` pair and `sd.errors = FALSE` (standard errors were skipped for the pilot; Lesson 6 replays selected fits to get them).
- `notes.md` sjSDM: delete the line `Sys.setenv(SJSDM_MOJO_BACKEND = "0")` and the blank line after it from `fit-sjSDM`. Replace the paragraph before the chunk ("Select the Python environment appropriate for your machine ... This switch explicitly disables Mojo in the recorded fork. ...") with: install sjSDM from fork release v0.1.0 and follow its installation instructions for the Python environment; with the later release v0.2.1, which the recorded fits used, run `Sys.setenv(SJSDM_MOJO_BACKEND = "0")` before loading sjSDM so that it uses PyTorch; keep the `df = 2` and seed sentences. No environment setup or download is part of this task.
- Report old 742 to 762: after the chunk, one sentence each for `sampling` (Monte Carlo draws per site for the likelihood), `step_size` (sites per optimisation batch), `lambda = 0` (no penalty on the environment or covariance) and `weight_decay` (the optimiser's own small penalty), checking each against sjSDM's documentation already quoted in `dev/simstudy/jsdm-package-comparison/DESIGN.md` or `FITTING-REPORT.md`; where neither records it, describe the argument only by what this lesson's text already says.
- Report old 765: the continuation is already stated plainly after the chunk; confirm.
- Report old 806 to 831: after `fitting-times`, read the runtime table with numbers inline from `runtime_summary`: seconds per selected fit for each package and the sjSDM total across attempts, as the practical answer to how long each takes.
- Report old 833: cut the parser-error history and the integration-accuracy sentence from the paragraph after the runtime table; `FITTING-REPORT.md` line 82 records the parser error, so link it in one clause.
- Reproduction record: keep the snapshot `d2ca508`, adding "(fork release v0.2.1, Mojo backend off)" after it; replace "input/source hashes" with "input and source hashes" in Task 5.

- [ ] **Step 10: The sjSDM TODO item**

Add this bullet to `TODO.md` under "## Review and maintenance", after the ggtern bullet Task 1 added (re-read the file first):

```markdown
- **Rerun the selected sjSDM start under fork release v0.1.0 (added 3 October 2026, lesson rewrite Phase C):** Lesson 5 now tells readers to install sjSDM from Doug's fork release v0.1.0 (`https://github.com/dougwyu/s-jSDM/releases/tag/v0.1.0`), while the recorded fits ran under v0.2.1 with `SJSDM_MOJO_BACKEND=0`. The lesson says the two give the same fit on an argument from reading the code: with the switch at zero, v0.2.1's dispatcher reaches the same PyTorch loss as v0.1.0 apart from integer casts on tensor sizes. Once the Python environment the stability work used is reachable, rerun the selected weak-penalty start (seed 26092811, the 3,000-epoch fit and its 1,000-epoch continuation in `dev/simstudy/jsdm-package-comparison/sjsdm-multistart.R`) under v0.1.0, compare its coefficients and loadings with the saved ones, and record the comparison in Lesson 5's appendix.
```

- [ ] **Step 11: Render, verify, check, commit**

```sh
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/occJSDM-lesson-5.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/occJSDM-lesson-5.Rmd", output_format="rmarkdown::html_vignette")'
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/revert-unchanged-png.zsh "vignettes/teaching-data/lesson-5-*"
git status --short vignettes/
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/run-step.zsh c31-L5-teach-verify-teaching Rscript dev/simstudy/jsdm-package-comparison/verify-teaching.R archive dev/simstudy/jsdm-package-comparison
```

Expected: renders clean; no PNG kept (no figure chunk changed); "Displayed simulation exactly reproduces the saved community and inputs." and "Teaching numerical verification passed for every check the committed archive supports." with `EXIT 0`. Diff the `.md` against `lesson-5-before.md`: computed tables unchanged apart from the new `str()` and reference-score outputs. Run the Global Constraints checks with `<files>` = `vignettes/occJSDM-lesson-5.Rmd TODO.md`. Then:

```bash
git add vignettes/occJSDM-lesson-5.Rmd vignettes/occJSDM-lesson-5.md vignettes/teaching-data TODO.md
git commit -m "Lesson 5 teaching pass: read every result, explain the choices, cite sjSDM fork v0.1.0

The SJSDM_MOJO_BACKEND line leaves the fitting chunk; the code-diff argument is
stated as an argument and TODO.md records the v0.1.0 rerun.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Lesson 5, copyedit commit

**Files:**
- Modify: `vignettes/occJSDM-lesson-5.Rmd`, `vignettes/occJSDM-lesson-5.md`

**Verifiers that read this lesson:** `verify-teaching.R`.

**Interfaces:**
- Consumes: Task 4's commit.
- Produces: the finished Lesson 5, linked from Lesson 6.

- [ ] **Step 1: Copyedit the whole lesson**

Apply the copyedit procedure. Report section 5 items for this lesson, mapped: long sentences (the exercise's opening, the probit paragraph, any the prose checks list); slash compounds "bulk/tail ESS" (if Task 4 left one), "input/source hashes", "gllvm/sjSDM"; the hypothetical-example tic ("For an explicitly hypothetical example" in section 4 and "As an explicitly hypothetical example" in section 7: reword one); "inspected" (if any remains); curly quotes around "converged" and "models predict better when deprived of observations".

- [ ] **Step 2: Render, verify, check, commit**

Render both outputs, apply the PNG policy to `"vignettes/teaching-data/lesson-5-*"`, rerun `verify-teaching.R` with step name `c32-L5-copy-verify-teaching`, run the Global Constraints checks with `<files>` = `vignettes/occJSDM-lesson-5.Rmd`, then:

```bash
git add vignettes/occJSDM-lesson-5.Rmd vignettes/occJSDM-lesson-5.md vignettes/teaching-data
git commit -m "Lesson 5 copyedit pass

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Lesson 6, teaching commit

**Files:**
- Modify: `vignettes/occJSDM-lesson-6.Rmd` (724 lines), `vignettes/occJSDM-lesson-6.md`, `vignettes/teaching-data/lesson-6-calibration-native-coefficient-intervals-1.png` (its chunk changes on purpose in Step 7; keep it only if its pixels change); other `lesson-6-*` PNGs only if their chunk changes
- Create (not committed): `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/lesson-6-before.md`
- Read: `style-diagnostic/lesson-4.md` items for old lines 857 to 1532 (sections 12 to 15 and the full results), the spec's Lesson 4 bullets on rare species, coverage and Hmsc, `dev/simstudy/jsdm-package-comparison/extension/PLAN.md` line 27 and `extension/fit-job.R` lines 25 to 35 (read-only), the finished Lessons 3 and 5

**Verifiers that read this lesson:** `dev/simstudy/jsdm-package-comparison/extension/verify.R` (bundle checks; it reads no lesson chunk).

**Interfaces:**
- Consumes: `teaching-data/lesson-4-extension.rds` (`extension`: `manifest` with `job`, `scenario`, `n_sites`, `n_species`, `response`, `use_traits`, `package`, `replicate`, `completed`, `fit_ok`, `diagnostic_pass`, `scored`; `overall` with `mae_pp`, `bias_pp`, `brier`; `species` with `species`, `mean_true_probability`, `training_presences`; `traits` with `trait`, `environment`, `estimate`, `lower`, `upper`, `truth`, `link`; `generating_curves`, `fitted_curves`, `all_attempted`), `teaching-data/lesson-4-calibration.rds` and `lesson-4-calibration-tables.R` (functions `calibration_report_table`, `calibration_comparison_table`, `calibration_summary_plot`, `calibration_checked_summary`).
- Produces: the kept headings "Which fits passed, remained flagged or failed?" and "Appendix: evidence and reproduction".

- [ ] **Step 1: Read, then keep the before state**

```sh
git show 7662f9a:vignettes/occJSDM-lesson-6.md > /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/lesson-6-before.md
```

- [ ] **Step 2: Opening ("What this lesson adds")**

- Phase B Task 4 M1: the objective "assess interval coverage alongside interval width" becomes "read interval coverage beside interval width, which the appendix tabulates".
- Insert the prerequisites block after the paragraph that says how sections 1 to 4 show their code; Lesson 6's verbs are `filter()`, `mutate()`, `transmute()`, `select()`, `arrange()`, `distinct()`, `group_by()`, `summarise()`, `across()`, `if_else()`, `coalesce()`, `recode()`, `left_join()`, `bind_rows()`, `as_tibble()`, `pivot_longer()` and `pivot_wider()`; name `tribble()` and `expand_grid()` as the unusual operations explained where they appear.
- Add the ggtern pointer (Shared text) after `load-extension-setup`.

- [ ] **Step 3: Section 1 ("Does more data help when ecology becomes harder?")**

- Phase B deferral: "The extension below repeats the experiment" and "The extension still assumes" become "This lesson repeats ..." and "This lesson still assumes ...".
- Report old 859: after "Within each community, the 100 training sites are part of the 300 training sites", add: so each pair of fits is paired, which is why lines connect the two site counts in the figure below.
- Report old 865 to 866: quantify the scenarios with a small run chunk before the scenario bullets: the mean true probability of each species in the baseline and rare scenarios, from `extension$species` for one package (`package == "occJSDM"`), 100 sites and ten species, averaged over communities; then say which three species become rare and how rare, inline, and that the baseline is the one-community design of Lesson 5 repeated.
- Report old 872 and section 5 slash: write out "package/dataset/model combinations" as "combinations of package, dataset and model", and add the breakdown: five ecological settings (baseline, rare, correlated, curved with a straight fit, curved with a squared term) by two site counts by ten communities by four packages is 400; two site counts by two species counts by traits supplied or omitted by ten communities by three packages is 240. Check both products against `extension$manifest` (`count(scenario, response, n_sites, use_traits, package)`) before writing them.
- Report old 901: replace the inline `if` sentence ("**Run status:** ...") with a fixed sentence (all 640 combinations were attempted), and add `stopifnot(extension$all_attempted)` as the first line after `readRDS()` in `extension-data` so the sentence cannot go stale.
- Report old 903 to 909: keep the passed, flagged and failed definitions and "failed is not crashed"; move the detailed acceptance rules (the gllvm and sjSDM agreement thresholds, the gradient rule) to a new appendix subsection "### Acceptance rules for each fit {.unlisted .unnumbered}" placed first in the appendix.
- Spec Hmsc bullet and report old 905: in the occJSDM and Hmsc paragraph, state the schedule: four chains, 2,000 warm-up and 4,000 retained draws; a fit that failed the checks got one longer attempt with 8,000 and 16,000; a fit still failing after that stays flagged and was not extended further (source: `extension/PLAN.md` line 27 and `extension/fit-job.R` lines 27 and 28). Add the flagged Hmsc count by scenario, inline from a small run chunk on `extension$manifest` (`package == "Hmsc", fit_ok, !diagnostic_pass`, counted by scenario), and say that whether still longer runs would have resolved them was not tested. This statement goes on Task 9's "needs Alex" list.
- Report old 907: add which scenarios the 14 failed gllvm fits come from, inline from `extension$manifest` (`package == "gllvm", completed, !fit_ok`, counted by scenario).
- Report old 923: the generating-curve axis label says "original simulation units"; add the clause that these are the simulator's raw values, before the training standardisation used for fitting.
- Report old 933 to 941: fix `extension-curved-design-example` so its comment matches its code, either by deleting "then use the training transformation" from the comment or by adding the transformation lines (centre and scale `test_predictors` with `training_predictors`' means and standard deviations); prefer adding the lines, as the reader asked how. Then say how the squared column reaches an occJSDM fit: add it as a column of `info` and name it in `occCovariates`.
- Spec "the sample-size and trait results" and report old 964 to 1017: after `extension-paired-summary`, read the table with numbers inline from `paired_summary`: in which scenarios more sites reduced error, by how much on average, and with how many paired communities; give `knitr::kable()` readable `col.names` (report old 1012).
- Report old 1019 to 1036: after `extension-fitted-curves`, read the figure (the straight fits cannot bend to the curve whatever the site count; the quadratic fits follow it), and make the caption say it is the curved scenario.

- [ ] **Step 4: Section 2 ("Do traits explain species responses ...?")**

- Report old 1051: gloss "species random slopes" in gllvm (each species gets its own response to each gradient, drawn from a shared distribution).
- Report old 1040 to 1125 (no trait-fitting code): before the trait figures, add the occJSDM call used in this experiment as a not-run chunk, taken from `extension/fit-job.R` lines 30 to 32, with one sentence that its inputs live in the fitting archive, not in this lesson's bundle, and that the trait-free arm omits `traits` from the list:

```{r occjsdm-trait-call, eval=FALSE}
# input$x: the standardised training environment; input$y: the 0/1 matrix;
# input$traits: one row per species with drought_tolerance and irrelevant_trait.
fit_with_traits <- occJSDM::runOccJSDM(
  list(info = input$x, OTU = input$y, traits = input$traits),
  occCovariates = names(input$x),
  listParams = list(n_factors = 2, n_lattrait = 0),
  MCMCparams = list(nchain = 4, nburn = 2000, niter = 4000, nthin = 1)
)
```

- Spec and report old 1055 to 1081: after `extension-trait-prediction-errors`, read it with numbers inline from `extension$overall` (`scenario == "traits"`): whether supplying traits lowered the average error, by how much, and whether more with 30 species than with 10, per package.
- Report old 1085 to 1121: after `extension-trait-coefficient-recovery`, read it inline from `extension$traits` (`link == "logit"`): how many drought-trait intervals contain the truth, and how many irrelevant-trait intervals exclude zero, by package and design.
- Report old 1123: fix "Its coefficient table records the true direction ..." so its referent is clear: the appendix's trait table ("Trait relationships in the same format") shows each relationship's truth beside the estimates and marks Hmsc's different link; or cut the sentence if that table does not show the true direction.

- [ ] **Step 5: Section 3 ("Across communities, are estimates biased?")**

- Report old 1171: replace "These are the September results" and "The older occJSDM-only validation study is a separate experiment and is not pooled here." with a description of the study (ten communities, each fitted at 100 and 300 training sites), dropping the dates.
- Report old 1176: define coverage in a clause where the summary bullets first use it, and fix "Coefficient coverage is more variable for occJSDM, gllvm and sjSDM" to say more variable than what (than probability coverage, across scenarios), checking the claim against the coverage tables.
- Report old 1181 and Phase B Task 4 M5: say why this section uses RMSE where section 1 used mean absolute error (RMSE gives larger errors more weight and is the usual measure across studies), so the reader does not compare section 1's numbers with these; restore the RMSE mechanics (square, average, root) and the two lost clauses: the prediction target is the probability averaged over unknown hidden conditions, and each community receives equal weight.
- Spec rare species and report old 1204: after "The rare-species scenario produces more overestimation, especially for occJSDM with 100 sites.", add what is known: the cause is not established; candidates are the priors, shrinkage towards the middle and sample size; link [Lesson 2](occJSDM-lesson-2.md#why-rare-and-common-species-are-pulled-towards-the-middle) and [Lesson 7](occJSDM-lesson-7.md)'s finding that the rarest species were unrecoverable in every arrangement, and name the bias recheck in `TODO.md`. Practical advice: treat rare-species estimates with caution and survey more sites where rare species matter.
- Report old 1214: add the clause that "original coefficient units" are the simulator's raw units, after undoing the standardisation.

- [ ] **Step 6: Section 4 ("Do 95% intervals contain the truth?")**

- Spec 82% coverage and report old 1240: after the occJSDM bullet in "How to read the baseline", say what the number means (about one interval in five misses, not one in twenty, so those intervals are too narrow) and that the cause is not established, with the same candidates and links as Step 5.
- Report old 1242: define "native intervals" where first used here: the intervals each package's own functions return.
- Report old 1256: after "More data cannot supply a curve that the fitted model leaves out.", add how to detect a missing curve in real data: compare fits with and without a squared term by their held-out scores, as Lesson 5 section 9 does.
- Phase B Task 4 M5: in the paragraph that begins "**occJSDM's baseline is not the whole story.**", restore the caveat that lower coverage in harder scenarios does not by itself identify a software defect.
- Report old 1276 and "What this lesson establishes": add a take-home for the occJSDM user: point predictions held up across communities; interval coverage did not always, so read occJSDM's intervals for environmental effects at 100 sites as too narrow, and add sites before trusting them.

- [ ] **Step 7: Appendix**

- Phase B deferral (`calibration_targets` renamed in `calibration-combined-targets`): in `calibration-combined-targets`, rename `calibration_targets` to `appendix_targets` and delete the redefinition of `baseline_calibration` (it equals the hidden setup's); in `calibration-combined-100`, `calibration-combined-300` and `calibration-combined-sensitivity`, use `appendix_targets` in place of `calibration_targets`. The rendered tables must be byte-identical to `lesson-6-before.md`.
- Phase B Task 4 note (`Warning in geom_text`): in `calibration-native-coefficient-intervals`, give the "Different link" labels their own data and map the package: replace the `geom_text(data = expand_grid(...), aes(x = "Hmsc", y = -Inf, label = "Different link"), ...)` layer with

```r
  geom_text(data = expand_grid(
      condition = c("baseline: linear", "correlated: linear", "curved: quadratic", "rare: linear"),
      measure = c("Coverage (%)", "Mean coefficient width")) |>
      mutate(package = "Hmsc"),
    aes(x = package, y = -Inf, label = "Different link"), inherit.aes = FALSE, vjust = -.7, size = 2.5) +
```

and confirm that the `.md` no longer prints the warning; keep the re-rendered PNG only if the revert script reports it changed.
- Phase B Task 4 M3: "These are original simulation units" becomes "The five settings are in original simulation units, before training-data standardisation."
- Report old 1532: cut the repeated probit and older-study caveats from the appendix's last paragraph, keeping the replay and refit facts.

- [ ] **Step 8: Render, verify, check, commit**

```sh
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/occJSDM-lesson-6.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
Rscript -e 'suppressMessages({library(dplyr);library(tidyr);library(tibble);library(ggplot2)}); rmarkdown::render("vignettes/occJSDM-lesson-6.Rmd", output_format="rmarkdown::html_vignette")'
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/revert-unchanged-png.zsh "vignettes/teaching-data/lesson-6-*"
git status --short vignettes/
Rscript dev/simstudy/jsdm-package-comparison/extension/verify.R /Users/douglasyu/src/occJSDM/dev/simstudy/results/lesson-4-extension-20260923 .
```

Expected: renders clean, with no `Warning in geom_text` in the `.md`; at most the coefficient-interval PNG kept; "Verified input identity, held-out states, independent probability truth and error summaries for 626 scored results; 640 of 640 fits attempted." Diff the `.md` against `lesson-6-before.md`: every computed table unchanged apart from the new run-chunk outputs. Run the Global Constraints checks with `<files>` = `vignettes/occJSDM-lesson-6.Rmd`. Then:

```bash
git add vignettes/occJSDM-lesson-6.Rmd vignettes/occJSDM-lesson-6.md vignettes/teaching-data
git commit -m "Lesson 6 teaching pass: read every result, state the schedules, restore the lost caveats

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Lesson 6, copyedit commit

**Files:**
- Modify: `vignettes/occJSDM-lesson-6.Rmd`, `vignettes/occJSDM-lesson-6.md`

**Verifiers that read this lesson:** `extension/verify.R`.

**Interfaces:**
- Consumes: Task 6's commit.
- Produces: the finished Lesson 6.

- [ ] **Step 1: Copyedit the whole lesson**

Apply the copyedit procedure. Items mapped from the report's section 5 and Phase B: the two coverage definitions at the start of section 4 (merge into one, Task 4 M2); "five fixed environmental settings" said twice in "All four packages together" (drop the first clause, Task 4 M2); "comparable trait-effect intervals" (say "trait-effect intervals on the same scale"); the long acceptance-rule sentences if any remain in section 1; slash compounds ("Rhat/ESS", "bulk/tail ESS", "gllvm/sjSDM", "grid/coefficient"); "The full results also show ..." and "The full tables below" where the appendix is meant (say "the appendix").

- [ ] **Step 2: Render, verify, check, commit**

Render both outputs, apply the PNG policy to `"vignettes/teaching-data/lesson-6-*"`, rerun `extension/verify.R` as in Task 6 Step 8, run the Global Constraints checks with `<files>` = `vignettes/occJSDM-lesson-6.Rmd`, then:

```bash
git add vignettes/occJSDM-lesson-6.Rmd vignettes/occJSDM-lesson-6.md vignettes/teaching-data
git commit -m "Lesson 6 copyedit pass

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Lesson 7 edges and the parked cross-lesson items

**Files:**
- Modify: `vignettes/occJSDM-lesson-7.Rmd` and `.md`, `vignettes/occJSDM-lesson-1.Rmd` and `.md`, `vignettes/occJSDM-lesson-5.Rmd` and `.md`, `vignettes/occJSDM-lesson-6.Rmd` and `.md`, `vignettes/LESSON-PLAN.md`, `TODO.md`

**Verifiers that read these files:** `dev/simstudy/spatial-design-sweep/verify-lesson.R` (bundle checks for Lesson 7); `verify-teaching.R` (Lesson 5); `extension/verify.R` (Lesson 6). Lesson 1 has no verifier.

**Interfaces:**
- Consumes: Tasks 5 and 7's finished Lessons 5 and 6.
- Produces: a reading path 0 to 7 in which every lesson ends with a navigation line.

- [ ] **Step 1: Lesson 7 opening and closing (edges only, `notes.md`)**

- At the end of the first paragraph of "What this lesson adds", add: This is the last lesson. It assumes [Lesson 2](occJSDM-lesson-2.md), whose two-stage model the eDNA-survey fits use, and reads the spatial fits' outputs in the way [Lesson 3](occJSDM-lesson-3.md) reads non-spatial ones: credible intervals, convergence checks and truth beside each estimate.
- Replace the closing navigation sentence after the exercises ("Return to the [Quickstart and lesson guide](occJSDM.md), [Lesson 2](occJSDM-lesson-2.md) or [Lesson 3](occJSDM-lesson-3.md).") with: "This is the last lesson. Return to the [Quickstart](occJSDM.md), to [Lesson 2](occJSDM-lesson-2.md) for fitting, or to [Lesson 3](occJSDM-lesson-3.md) for reading outputs."
- Check every cross-reference in Lesson 7 for the new order (`/usr/bin/grep -n -E "Lesson [0-9]|lesson-[0-9]" vignettes/occJSDM-lesson-7.Rmd`): each must name the lesson that now teaches the thing referred to. No other prose changes (owner's decision: no second prose pass).
- Rename the chunk labels `lesson-2-load` to `lesson-7-load` and `lesson-2-contents` to `lesson-7-contents`. Neither chunk draws a figure, and no verifier reads them (the only other mention is the sweep's historical `IMPLEMENTATION.md`, out of scope).

- [ ] **Step 2: The parked cross-lesson items**

- Lesson 1's closing: append, after the last bullet of "What you can now recognise in the other lessons", the paragraph "Continue to [Lesson 2: Fit the model and compare its answers with truth](occJSDM-lesson-2.md)."
- The link label: in Lessons 5 and 6, replace `[Quickstart and lesson guide](occJSDM.md)` with `[Quickstart](occJSDM.md)`; Lesson 7's was replaced in Step 1. Confirm with `/usr/bin/grep -n "lesson guide" vignettes/*.Rmd` (expected: no output).
- `vignettes/LESSON-PLAN.md`: set line 3 to "Last updated: " followed by today's date in the form `3 October 2026` (from `date '+%-d %B %Y'`); in line 29, delete "The spatial lesson is not a prerequisite for Lesson 3." and replace "Lesson N was renamed Lesson 4 on 23 September 2026 once its content was settled." with "The four-JSDM lesson, now Lesson 5, was first numbered Lesson 4 on 23 September 2026, once its content was settled."; in the dated log, replace "(link added when opened)" with "([PR #19](https://github.com/AlexDiana/occJSDM/pull/19))" in the four Phase A entries (Lesson 1 teaching pass, intuition lesson, Lesson 0 teaching pass, Lesson 2 teaching pass) and with "([PR #20](https://github.com/AlexDiana/occJSDM/pull/20))" in the Phase B entry. Confirm with `/usr/bin/grep -c "link added when opened" vignettes/LESSON-PLAN.md` (expected 0).
- `TODO.md`: in the bullet that begins "**Spatial amplitude and field structure after the site-arrangement sweep (Lesson 7) (added 2 October 2026):**", merge the double parenthetical into "(Lesson 7; added 2 October 2026)".

- [ ] **Step 3: Render, verify, check, commit**

Render Lessons 1, 5, 6 and 7 to markdown and HTML with the Global Constraints commands. Apply the PNG policy to `vignettes/occJSDM-lesson-1_files`, `"vignettes/teaching-data/lesson-5-*"`, `"vignettes/teaching-data/lesson-6-*"` and `"vignettes/teaching-data/lesson-7-*"` (no PNG should be kept). Confirm that `git diff vignettes/occJSDM-lesson-7.md` shows only the two prose changes, and the Lesson 1, 5 and 6 `.md` diffs only the navigation lines. Then:

```sh
Rscript dev/simstudy/spatial-design-sweep/verify-lesson.R .
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/run-step.zsh c41-L5-edges-verify-teaching Rscript dev/simstudy/jsdm-package-comparison/verify-teaching.R archive dev/simstudy/jsdm-package-comparison
Rscript dev/simstudy/jsdm-package-comparison/extension/verify.R /Users/douglasyu/src/occJSDM/dev/simstudy/results/lesson-4-extension-20260923 .
```

Expected: "Lesson 2 bundle verified against the committed results." (the sweep's old wording); the two comparison verifiers' lines as in Tasks 4 and 6. Run the Global Constraints checks with `<files>` = `vignettes/occJSDM-lesson-1.Rmd vignettes/occJSDM-lesson-5.Rmd vignettes/occJSDM-lesson-6.Rmd vignettes/occJSDM-lesson-7.Rmd vignettes/LESSON-PLAN.md TODO.md` (the long sentences in `LESSON-PLAN.md` are dated or status text outside this pass and are reported, not split). Then:

```bash
git add vignettes/occJSDM-lesson-1.Rmd vignettes/occJSDM-lesson-1.md vignettes/occJSDM-lesson-5.Rmd vignettes/occJSDM-lesson-5.md vignettes/occJSDM-lesson-6.Rmd vignettes/occJSDM-lesson-6.md vignettes/occJSDM-lesson-7.Rmd vignettes/occJSDM-lesson-7.md vignettes/LESSON-PLAN.md TODO.md
git commit -m "Lesson 7 edges and the parked cross-lesson fixes

Lesson 7 says it is last and what it assumes, ends with a navigation line and
renames its lesson-2 chunk labels; Lesson 1 gains its Continue line; the
Quickstart link label; LESSON-PLAN.md's date, reading-order sentence and pull
request links.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 9: Final verification, HTML renders and the pull request draft

**Files:**
- Modify: `vignettes/LESSON-PLAN.md` (the Phase C decision-log line)
- Create (not committed): `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/pr-body.md`
- No other repository change unless a check fails, in which case the fix is committed with a message naming the check.

**Interfaces:**
- Consumes: the branch after Task 8.
- Produces: a branch ready to push, every check's output, HTML renders of all eight lessons in `vignettes/` (gitignored) for Doug's read, and the drafted pull request body.

- [ ] **Step 1: Record the pass in the lesson plan and commit**

Append this entry to the end of the decision log in `vignettes/LESSON-PLAN.md` (after the Phase B entry, separated by one blank line), with today's date from `date '+%-d %B %Y'`, deleting any clause about an item that the pull request body lists as deferred:

```markdown
- **3 October 2026, Phase C prose passes:** combined Phase C pull request (link added when opened). Lessons 3 to 6 each had a teaching pass and a separate copyedit pass. Each lesson now lists the R it assumes, shows the data objects it loads, defines its terms before using them, reads every printed result with numbers computed from the bundles, explains each choice and argument, and gives advice for the reader's own survey. Lesson 3 explains the ggtern theme, starts its diagnostics with an overview and cuts its repeated passages; Lesson 4 defines WAIC and says what the package's new-site function can and cannot return; Lesson 5 cites sjSDM fork release v0.1.0 as the version to install, with the rerun under that release recorded in `TODO.md`, and shows occJSDM's own new-site route beside the gllvm exercise; Lesson 6 states the Bayesian fitting schedule and the failed and flagged fits by scenario, and restores the caveats Phase B lost. Lesson 7 had its edges only (Doug's decision): it says it is the last lesson and what it assumes, and ends with a navigation line. Lesson 1 gained its Continue line and the lessons' quickstart link is labelled "Quickstart".
```

```bash
git add vignettes/LESSON-PLAN.md
git commit -m "Record the Phase C prose passes in the lesson plan

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

- [ ] **Step 2: Every verifier**

```sh
S=/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/run-step.zsh; A=/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples
zsh $S c90-test_lesson Rscript dev/simstudy/vignette-lesson/test_lesson.R
zsh $S c91-verify_lesson Rscript dev/simstudy/vignette-lesson/verify_lesson.R $A
zsh $S c92-verify_outputs Rscript dev/simstudy/vignette-lesson/verify_outputs.R $A
zsh $S c93-verify_diagnostics Rscript dev/simstudy/vignette-lesson/verify_diagnostics.R $A
zsh $S c94-verify_latent_tables Rscript dev/simstudy/vignette-lesson/verify_latent_tables.R $A
zsh $S c95-verify_native_plots Rscript dev/simstudy/vignette-lesson/verify_native_plots.R $A
zsh $S c96-native-traits-verify Rscript dev/simstudy/vignette-lesson/native-traits-verify.R $A
zsh $S c97-ordination-verify Rscript dev/simstudy/vignette-lesson/ordination-verify.R $A
zsh $S c98-remaining-plots-verify Rscript dev/simstudy/vignette-lesson/remaining-plots-verify.R $A
zsh $S c99-prediction-verify Rscript dev/simstudy/vignette-lesson/prediction-verify.R $A $A/prediction
zsh $S c100-unbalanced-verify Rscript dev/simstudy/vignette-lesson/unbalanced-verify.R $A/unbalanced
zsh $S c101-verify-teaching Rscript dev/simstudy/jsdm-package-comparison/verify-teaching.R archive dev/simstudy/jsdm-package-comparison
Rscript dev/simstudy/jsdm-package-comparison/extension/verify.R /Users/douglasyu/src/occJSDM/dev/simstudy/results/lesson-4-extension-20260923 .
Rscript dev/simstudy/spatial-design-sweep/verify-lesson.R .
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
node --test dev/simstudy/lesson-site/test_lessons.js
Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/check-frozen-chunks.R
```

Expected last lines, in order: `test_lesson.R` passes; "All teaching-evidence checks passed."; "Trait decomposition, sampling-effort truth, source hashes and fit hashes verified."; "Diagnostic selection, fit hashes, source hashes and plotting identities verified."; "Truth joins survive row permutations and reject duplicate identities; native table builds."; "Twelve PNGs, compact export, full-fit hash, source hashes and exact teaching-code provenance verified."; the traits PASS line; the ordination PASS line; "Five PNGs, compact bundle, full-fit/library/source hashes, and exact displayed code verified."; "All prediction verification checks passed."; the unbalanced "... and displayed code." line; "Teaching numerical verification passed for every check the committed archive supports."; "... 626 scored results; 640 of 640 fits attempted."; "Lesson 2 bundle verified against the committed results."; "Shared lesson-link regression checks passed."; `fail 0`; "All 36 frozen chunks match main at 7662f9a." The four Lesson 3 logs each contain "Lesson 3 displays the exact exported ...".

- [ ] **Step 3: Counts and greps**

```sh
git diff --name-only 7662f9a -- ':!*.png' ':!*.rds' | xargs /usr/bin/grep -n $'\u2014' | wc -l
git log 7662f9a..HEAD --format=%B | /usr/bin/grep -c $'\u2014'
git log 7662f9a..HEAD --format=%B | /usr/bin/grep -c "Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git log --oneline 7662f9a..HEAD | wc -l
zsh /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/prose-checks.zsh vignettes/occJSDM-lesson-[0-7].Rmd vignettes/LESSON-PLAN.md TODO.md
Rscript /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-B/check-lesson-links.R vignettes/occJSDM-lesson-[0-7].Rmd vignettes/occJSDM-lesson-[0-7].md vignettes/occJSDM.Rmd vignettes/LESSON-PLAN.md TODO.md AGENTS.md dev/simstudy/vignette-lesson/README.md dev/simstudy/jsdm-package-comparison/README.md
git diff --stat 7662f9a -- R src DESCRIPTION NAMESPACE dev/simstudy/spatial-design-sweep vignettes/teaching-data ':!vignettes/teaching-data/*.png' dev/simstudy/vignette-lesson/*-examples.Rmd dev/simstudy/vignette-lesson/*verify*.R dev/simstudy/vignette-lesson/*export*.R dev/simstudy/jsdm-package-comparison
/usr/bin/grep -n "lesson guide\|link added when opened" vignettes/occJSDM-lesson-[0-7].Rmd vignettes/LESSON-PLAN.md
wc -l vignettes/occJSDM-lesson-[0-7].Rmd
```

Expected: 0; 0; the trailer count equals the commit count (eleven commits if no check needed a fix: this plan's commit and the ten of Tasks 1 to 9); the prose checks' "must be 0" counts are 0 for the lesson Rmds (long sentences remaining in Lessons 3 to 6 each have a reason recorded in the task reports; the counts for Lessons 0 to 2 and the living documents are reported, not edited); the link checker's success line; empty `git diff --stat` (nothing out of scope changed); one line from the grep, the new Phase C log entry's "(link added when opened)", which the controller fills; the line counts, compared with 586, 422, 1166, 1822, 362, 848, 724 and 609 at 7662f9a.

- [ ] **Step 4: HTML for Doug's read**

Render each of the eight lessons with the HTML command (Lessons 0 and 2 included, unchanged, so the set is complete and current). Apply the PNG policy afterwards to `vignettes/` (HTML renders of Lessons 0 to 4 rewrite nothing tracked; Lessons 5 to 7 write HTML figures to their gitignored `_files/figure-html/` folders). Confirm `git status --short` is clean and list the eight `vignettes/occJSDM-lesson-*.html` paths.

- [ ] **Step 5: Draft the pull request body**

Write `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/pr-body.md` (markdown rules apply; one line per paragraph; no em-dashes) with these sections:
- Summary: one paragraph on Phase C (teaching then copyedit passes on Lessons 3 to 6, two commits each; Lesson 7 edges only by Doug's decision; the parked cross-lesson items), and that no bundle, exporter, verifier or package file changed and every verifier passes.
- Per lesson: what the teaching pass added (by section, three to six bullets), what the copyedit pass changed, and the line count before and after.
- Doug's Phase C decisions and how each was carried out: Lesson 7 edges; sjSDM v0.1.0 citation, the dropped `SJSDM_MOJO_BACKEND` line, the code-diff argument stated as an argument and the `TODO.md` rerun item; no line target for Lesson 3.
- Decisions this plan took (the "Decisions this plan takes" bullets that the owner may want to revisit: the frozen chunks, the label "Quickstart", the PR links in the log, the visible truth join, the ordination order, the Hmsc schedule source).
- Needs Alex: (1) the cause of the ggtern failure and whether occJSDM can avoid it (Lesson 3's ggtern sentence, pointed to from Lessons 4 to 6; `TODO.md` item); (2) the Hmsc schedule statement in Lesson 6 (four chains of 2,000 warm-up and 4,000 retained draws, one longer attempt of 8,000 and 16,000, flagged fits not extended further; source `extension/PLAN.md` line 27 and `extension/fit-job.R` lines 27 and 28), and any Hmsc schedule sentence Lesson 5 now has. Add any other package fact a task added that the spec marks unconfirmed.
- For Doug to confirm: Lesson 3's statement that Ji et al. (2025) Supplementary Information 12 works through cases where the conditional and predictive probabilities disagree (the report could not check it from the repository).
- Deferred, with reasons: one of the near-identical native trait chunks and the literal MCMC settings in `prediction-fit-one-factor` (frozen for their verifiers); hiding the truth join (all teaching code stays visible); the stale extension dev copies and their README note (hash-stamped sources); verifier polish from Phase B (unasserted closing messages, the ordination arrow-start check: verifier changes are out of scope); the old lesson numbers in the header comments of `export-calibration.R` and `simulate-pilot.R` (hash-stamped); `FITTING-REPORT.md` line 3's displayed "Lesson 4" (historical prose, link target already fixed); the sweep verifier's "Lesson 2" message (sweep frozen); the `####` appendix headings that drop out of Lesson 5's table of contents (appendix navigation only; `toc_depth: 4` would list them if Doug wants); the message of 795c34a, which omits three amended changes (merged history, not rewritable); `git log --follow` on reused stems (use `git diff -B -M`); and any report item a task recorded as deferred.
- Collision with PR #14: the paragraph "Known collision with PR #14" of this plan, naming the edited paragraphs (Lesson 3's WAIC paragraph at the end of its diagnostics section; Lesson 4's "Check the additional fit and understand the WAIC limitation" lead-in sentence, its definition and the plain-language sentence before the bold advice; the score-difference paragraphs moved to Lesson 4's appendix).
- Checks run: the output lines of Steps 2 and 3.
- A note that Doug reviews line by line, with the eight HTML paths.
- The last line: `🤖 Generated with [Claude Code](https://claude.com/claude-code)`.

- [ ] **Step 6: Hand back to the controller**

Report: the commit list (`git log --oneline 7662f9a..HEAD`), the check outputs, the PR body path, and any check that needed a fix. The branch is ready to push.

Controller steps, not for an implementer:
- Push `codex/lesson-phase-c` and open the pull request with `gh pr create --base main --title "Lesson rewrite Phase C: teaching and copyedit passes on Lessons 3 to 6, Lesson 7 edges" --body-file /Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/pr-body.md`; then replace "(link added when opened)" in the Phase C decision-log entry with the pull request link in a follow-up commit with the trailer.

---

## Self-review notes

- **Spec coverage, Lesson 3 bullets.** Restructure (diagnostics first, package plots with truth, appendix guidance folded, prediction split, simulation-only material, mirror study and PR #13 history to the appendix, reproduction separated from the finder): done in Phase B, confirmed in Task 1 Step 15; the 1,150-line target is replaced by `notes.md` (no target), and Tasks 1 and 2 cut repetition and seams instead. Posterior draws, chains, credible intervals and log-odds at first use: Task 1 Step 4 (terms list); hidden site factors at first mention: Step 3; ordination before use: Step 10; variation partitioning: section removed in Phase B by Doug's decision, not applicable; WAIC before use: Step 5 (Lesson 3) and Task 3 Step 8 (Lesson 4); the three teaching objects shown: Step 4. Traces through the package route with internal slots in the appendix and the accessor sentence: done in Phase B, confirmed in Step 15. The `ggtern::theme_bw()` sentence: Step 4 with the symptom confirmed while planning, the cause on the "needs Alex" list and the `TODO.md` item in Step 14. Three uncited references removed, three Lesson 2 statements fixed, the `sampleresults` paragraph cut, three pipe tables converted: done in Phase B, confirmed in Step 15.
- **Spec coverage, old Lesson 4 bullets (now Lessons 5 and 6).** Restructure, the section 3 reminder, the exercise after section 4, section 10 and the record to the appendix, the full results reduced, the two-optima point with its evidence in the appendix, the pipe tables, the echo leak: done in Phase B. Read the printed results (identical curves, band shrinkage, sample-size and trait results, averaging illustration): Task 4 Steps 5 and 7, Task 6 Steps 3 and 4. Explain the choices (gllvm rather than occJSDM, with `useBiotic` shown alongside; the fixed factor count; the section 10 arguments): Task 4 Steps 4, 5 and 9. Rare-species overestimation and the 82% coverage: Task 6 Steps 5 and 6. Hmsc's 86 flagged fits: Task 6 Step 3, with the schedule confirmed from `extension/PLAN.md` and `fit-job.R`; the confirmed schedule includes one longer attempt, so the lesson says the fits were not extended beyond it rather than "not extended". Line 765 causal claim: Phase B hedged it; Task 4 Step 1 confirms the stability report does not establish the link, and Step 6 keeps the hedge. Cross-reference instead of re-teaching marginal versus conditional and WAIC: Task 4 Steps 4 and 5 (now Lesson 4's). Covariate standardisation: Task 4 Step 3 (and Lesson 3 in Task 1 Step 6, Lesson 4 in Task 3 Step 3). The section 3 move: done in Phase B. sjSDM version: Task 4 Steps 4, 9 and 10, as `notes.md` overrides (cite v0.1.0 now, drop the switch line, state the argument as an argument, TODO item for the rerun; the rerun itself is not done).
- **Report section 2 coverage.** Lesson 3 report, old lines 27 to 785 and 1086 to 2085: Task 1 Steps 3 to 13, each item named by its old line; old lines 786 to 1085: Task 3 Steps 3 to 8. Lesson 4 report, old lines 34 to 855: Task 4 Steps 2 to 9 (section 3's items 221 and 226 moved with that material to Lesson 1 in Phases A and B and are outside Phase C); old lines 857 to 1532: Task 6 Steps 3 to 7. Section 5 copyedit items: Tasks 2, 3 (Step 11), 5 and 7.
- **Phase B deferred items.** Plan B "Deferred to Phase C": prose passes (Tasks 1 to 7), prerequisites block (Tasks 1, 3, 4, 6; Lesson 7 already has one), R-gloss cuts (the copyedit procedure and Global Constraints), the ggtern explanation (Task 1), the sjSDM v0.1.0 citation and switch line (Task 4), the line 765 claim (Task 4 Steps 1 and 6), the `calibration_targets` rename (Task 6 Step 7), the Lesson 3 trim (replaced by `notes.md`). PR #20 "Deferred and Phase C": Lesson 1's Continue line, the "Quickstart and lesson guide" label and `LESSON-PLAN.md`'s date and line 29 (Task 8 Step 2); the 795c34a message (not rewritable, listed in the PR body); Lesson 3's "the figure below" seam, stranded per-chain guidance, `fitmodel_perfect` before its load and "You do not need the spatial lesson first" (Task 1 Steps 3, 5 and 6; Task 4 Step 2 for Lesson 5); Lesson 6's "The extension below", coverage defined twice, the lost clauses and the `geom_text` warning (Task 6 Steps 3, 5, 6 and 7, Task 7); Lesson 7's chunk labels (Task 8 Step 1); the stale extension dev copies (deferred, hash-stamped); the render-noise PNGs (PNG policy); `git log --follow` (PR body). Final review CAN WAIT items: Task 1 M1, M2 and M5 and Task 7 M2, M5 and the sweep note are verifier, process, hash-stamped or historical items listed as deferred in the PR body; Task 1 M3 is the PNG policy; Task 1 M4, Task 3's items and seams are in Task 1; Task 2 M1 and Task 4 M6 were resolved in Phase B; Task 2 M2 is Task 3 Step 8; Task 4 M2, M3, M5, M7 and the note are Tasks 6 and 7 or deferred; Task 5 M3 and M4 are Task 4 Steps 5 and 2; Task 5 M5 is deferred with its option; Task 6 M1 is the PR body and M2 Task 8; Task 7 M3 and M4 are Task 8. The Phase-C-labelled minors in `task-3-review.md` (Minors 2, 3, 5, 7), `task-4-review.md` (Minors 1 to 5, 7 and the note) and `task-5-review.md` (Minors 1, 3, 4, 5) are each placed above; Task 3 Minor 4 (`species_order` unused) no longer applies, because the restored custom figure uses it.
- **Placeholder scan.** No "TBD" or bracketed fill-ins. The `<stem>`, `<files>` and `<path>` in the Global Constraints commands are substitutions defined where they appear; the dates come from `date`. Where a step says "inline", it names the bundle object and the filter; where the wording of a sentence matters (the ggtern sentence, the sjSDM paragraph, the two TODO items, the decision-log entry, the navigation lines), it is quoted. Two steps leave a reasoned choice to the executor with both branches written out (the stability check's rationale in Task 4 Step 6, the sjSDM argument glosses in Step 9), because the recorded documents decide which branch is true.
- **Name consistency.** The tool paths `/Users/douglasyu/src/occJSDM-worktrees/.sdd/LESSON-STYLE-PLAN-C/{run-step.zsh,revert-unchanged-png.zsh,check-frozen-chunks.R,prose-checks.zsh}` and `.sdd/LESSON-STYLE-PLAN-B/check-lesson-links.R` are used identically in every task. New chunk labels (`show-teaching-objects`, `display-native-state-table-markdown` and its probability twin, `occjsdm-new-site-route`, `reference-scores`, `occjsdm-trait-call`, `lesson-7-load`, `lesson-7-contents`) are unique within their lessons and match none of the verifiers' patterns (`native-`, `ordination-`, `remaining-plots-` with `purl=TRUE`; `prediction-fit-one-factor`, `prediction-native-call`, `reproduce-community`). The anchor `#what-this-lesson-answers` is produced by Lesson 3's existing heading and consumed by the pointer text in Tasks 3, 4 and 6. Lesson 4's heading "Compare models using what actually occurred" is linked from Task 4 Step 4 only if it survived Task 3, which keeps it. Run-step names (`c01` to `c07`, `c11` to `c18`, `c21`, `c22`, `c31`, `c32`, `c41`, `c90` to `c101`) are unique.
- **Choices a reviewer may question.** Splitting the latent-table chunks into an HTML chunk and a hidden Markdown chunk changes visible code; it follows the report (the site-build branch is not for the reader) and the step requires the rendered tables to stay identical. Adding transformation lines to `extension-curved-design-example` changes a visible chunk, chosen over deleting the comment because the reader asked how. The Lesson 6 appendix rename touches four chunks, and the step requires byte-identical tables.

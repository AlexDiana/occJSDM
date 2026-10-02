# Lesson rewrite in the teaching style: design

Written 2 October 2026 from the five diagnostic reads in `style-diagnostic/` and Doug's decisions in discussion the same day. This is the spec the implementation plan argues from. Doug reviews each lesson line by line before it is published, so every change here is one he will see in a diff.

## Goal

Rewrite Lessons 0 to 4 so that every section teaches rather than describes, following the criteria in `AGENTS.md` ("Writing the vignettes and lessons"), with the quickstart as the exemplar. Keep the verified results, bundles and figures wherever the criteria do not require a change.

## The reader

An ecologist who runs eDNA or presence/absence surveys, knows basic R and the tidyverse pipe, and has not fitted an occupancy model or a JSDM before. A quantitative ecologist should still find the technical evidence, below the point where the first reader can stop.

## The criteria, applied

Each lesson, section by section, must: name the concepts the reader needs in the reader's terms before using them; answer the question the reader is about to ask; show every part of the data or output the reader will handle; spell out the reasoning behind every choice (numbers of samples and PCRs, factor counts, priors, which function or package demonstrates what); and give practical advice with complete examples, including what the reader should do differently on their own survey. Then a separate copyedit pass: split long sentences, replace "and/or", fix stops, remove asides that repeat their sentence, check every factual claim against the output. The two passes are two commits per lesson so both are visible in review.

## Decisions taken on 2 October 2026

- **Pass depth.** Lessons 0, 1 and 2: a prose pass with the small reorderings their reports name. Lessons 3 and 4: restructure first, then the prose pass, because a third of Lesson 3 and the appendices of Lesson 4 serve a reviewer or a methods reader rather than the ecologist, and teaching prose on top would make them longer and harder.
- **Splits and renumbering, once, after the beta.** New-site prediction and model comparison leaves Lesson 3 as its own lesson. The ten-community and trait extension (sections 12 to 15) leaves Lesson 4 as its own lesson. Lesson 2 moves to the end, whole, as already decided in `TODO.md`. All three happen in one renumbering step after the beta tag and Alex's review of PR #14, so every cross-link, the site list, the lesson plan and the memory notes change once. Proposed numbering for Doug to confirm: 0 simulate (optional); 1 fit the model and the detection stages; 2 understand the outputs; 3 predict at new sites and compare models; 4 compare four JSDMs on one community; 5 ten communities and traits; 6 spatial landscapes and survey design. If seven lessons is too many, 5 can stay inside 4 as a clearly separated second part.
- **Figures in Lesson 3.** The main figure of each section becomes the package's own plot (for example `plotOccupancyCovariates()`) with the simulated truth overlaid, so the reader sees the call they will use and the explanation in one place; the custom ggplot2 duplicate built from teaching summaries goes, and the appendix's unique guidance folds into the sections. A custom figure stays wherever the teaching point is a comparison the native plot cannot draw (for example perfect-observation against PCR-only fits in one panel), and the prose says why that figure is custom.
- **Technical material becomes an appendix in each lesson.** Reproduction records, hashes, exporter and verifier notes, the Lesson 3 mirror study, the sjSDM two-optima evidence, Lesson 4's full-results tables and the "fit each package yourself" arguments move below each lesson's closing section, under an "Appendix" heading, rather than out of the file. `dev/simstudy` is not on the Pages site and the repository may be inaccessible to ordinary ecologists; the appendix keeps the evidence reachable for quantitative readers and out of the reading path for the rest. The lesson's closing section says what the appendix holds.
- **R glosses.** Cut explanations of the pipe, `select()`, `filter()`, `mutate()` and the other routine verbs, `readRDS()`, lists and `$`. Keep a gloss where the operation is unusual (a join on species identity, a pivot, an array slice) and the ecology reads through the code.
- **New run chunks are in scope** when they show data the reader will handle or read an output the prose discusses: a `head()` of each loaded object, a count of each kind of detection error, a cross-tabulation of presence by detection. They are small, run at knit, and use only the bundle.
- **Pipe tables become bullet lists**, per the repository rule for files edited in Visual mode (Lesson 0 has 13 table lines, Lesson 1 27, Lesson 3 32, Lesson 4 16).
- **Hard-coded result numbers in prose** are replaced by inline R from the bundle wherever the number comes from the bundle, so the prose cannot go stale; a number that is not in the bundle is either added to the bundle by its exporter or stated with its source.
- **Cross-references** are written with the current numbers until the renumbering step, which updates them all; lesson text should refer to a lesson by what it teaches ("the outputs lesson") where a number would otherwise be the only identifier.

## Sequence

- **Phase A, now (the lessons are withheld from the beta, so nothing here touches its path).** Prose passes in reading order: Lesson 1 first, since every reader takes it; then Lesson 0; then Lesson 2. One branch and one pull request per lesson, two commits each (teaching, copyedit), Doug reviewing line by line. Each pass includes that lesson's small reorderings, its appendix, its pipe-table conversion and its factual corrections.
- **Phase B, after the beta tag and PR #14.** The combined restructure and renumbering: Lesson 2 to the end; prediction out of Lesson 3; the extension out of Lesson 4; appendices created in 3 and 4; Lesson 3's figures switched to package plots with truth overlaid; diagnostics moved to the front of Lesson 3; the Lesson 4 echo leak fixed. Every cross-link, title, `VignetteIndexEntry`, PNG prefix, the site list and its test, `LESSON-PLAN.md`, the Quickstart and README, and the memory notes updated once. One pull request, verified by `test_lesson_links.R`, the site test and every lesson verifier.
- **Phase C.** Prose passes on the restructured lessons in their new order (2 outputs, 3 prediction, 4 comparison, 5 extension, 6 spatial), one pull request each, Doug reviewing.

## Per-lesson work, with the recommended answers to each report's questions

Each report's "Questions for Doug" section is answered here with a recommendation; Doug overrides any of them in his review of this spec. Items marked "needs Alex" are facts about the package that the lesson should state only once confirmed.

### Lesson 0 (report: `style-diagnostic/lesson-0.md`)

- Add one passage laying out the two-stage process in the reader's terms before the detection settings: field collection can miss a species that is present (false negative) or pick up DNA that is not from the site (false positive); PCR can miss DNA that is in the sample or report DNA that is not there; `theta`, `theta0`, `p` and `q` are then introduced as the four probabilities of those events. The full explanation lives in Lesson 1, which every reader takes; Lesson 0 carries a short version and a link.
- Keep the order "specify, then inspect", which the title promises, and add a brief look at each object where it is first created.
- Replace the roadmap prose about Lesson 2's status (line 26) and the deferred dispersal work (line 461) with a single forward link each; the renumbering step rewrites them.
- Read every printed output: the Site 3 and Sample 6 table, a survey-wide count of each kind of error, and the OTU_1 against OTU_10 contrast that the settings set up but the maps never draw out.
- Split the "No detection" map label into true negatives and missed detections; Lesson 1 builds its own labels from its own bundle, so nothing downstream breaks.
- Rename the reused variables `samples_per_site` and `pcrs_per_primer` where they are reassigned; a small code change the prose then need not warn about.
- Sample identifiers: Lesson 1's unbalanced fit ran with the same three samples removed and their original identifiers kept, so non-consecutive sample IDs are accepted and the prose can say so; unequal numbers of PCRs per sample are not demonstrated anywhere, so the prose says that is untested (needs Alex to confirm whether it is supported).
- Cut the routine R glosses; keep the `anti_join()` and paired-row explanation, which carries the ecological point that a lost sample is absent rows.

### Lesson 1 (report: `style-diagnostic/lesson-1.md`)

- Show the data: `head()` of `survey_data$info`, `$OTU` and `$traits`, and a one-line overview of the `lesson` object, immediately after loading.
- State the two main findings in prose where their outputs appear: occupancy estimates pulled towards the middle (low band 6.5% true against 20.7% estimated), and field-stage false positives mostly believed (mean site probability 67.6% at sites that were all unoccupied) while laboratory false positives are handled well.
- Move convergence ("Are the calculations stable enough to interpret?") before the results, as the quickstart does.
- Keep the fitting reference for the reader's own data as its own section after the two fits, since it is practical advice; move the recovery figure out from under it to the results.
- Collapse the threshold-conversion explanation to one sentence and a link to Lesson 0.
- Keep the error maps but shorten them; they answer where errors occur, which readers ask.
- Keep the unequal-replication section as a demonstration of about 40 lines (what to remove, what the fit needs, that it runs), and move its full truth comparison and second diagnostic system to the appendix.
- Move provenance details and the "Paper2Agent" and "beta-release error target" disclaimers to the appendix or remove them; no reader of the lesson will recognise them.
- Prior advice (line 432): recommend a basis for changing the contamination priors, namely contamination rates seen in blanks and negative controls, phrased as a suggestion (needs Alex's endorsement before it is stated as package advice).
- Choosing `n_factors` and `n_lattrait` on real data: say plainly that there is no rule yet, start small, and check that conclusions do not change when the count changes.
- Explain why the lesson uses two samples, two primers and six PCRs while the quickstart's `sampledata` uses three, three and two, if the lesson plan records a reason; otherwise say that the structure matters and the counts are a choice.
- The case table: show "not applicable" for the true collection probability in the field-stage case.
- Leave a one-sentence hook for the pond case study in `TODO.md` ("real surveys often pre-filter detections; a later addition contrasts that with modelling them"); the case study itself is a later addition.
- Correct the facts the report doubted: there are six not-run chunks, four of them fits; recompute the memory figure at line 239 from the stated dimensions and state the arithmetic.
- Cut routine tidyverse explanations (lines 56, 260, 306, 356, 488), keep the non-obvious ones.

### Lesson 2 (report: `style-diagnostic/lesson-2.md`)

- Push the sweep half towards the ecologist: name the species prevalence groups, how each arrangement was laid out, the three arms (including that the true-state fit is the JSDM-only mode), the error measure and the survey's detection settings before they are used; define the informative, intermediate and uninformative thresholds before their first use; move the data-loading chunk into 2A and the reading-label table from 2B to 2C; move the convergence material to the appendix.
- Close with explicit advice matching the quickstart's line 64: leave `spatCovariates` out in the beta unless sites are clustered within about one range of each other and the species of interest are common; and show how to put a plausible range for the reader's own system on the model's standardised scale.
- Keep the plotting and wrangling code visible, for consistency with the other lessons' "all teaching code is shown", and because the exercises use the bundle.
- State that the teaching bundle ships with the package's vignettes once the lessons are published, so the exercises run from the installed package.
- Cut the line 83 sentence about an earlier computational validation the reader has not met; replace the PR #8 reference at line 55 with the reasoning for using every location as a support point; keep the close-pairs design advice (lines 71 to 77) with a note that the sweep found no measurable benefit at this budget and range.
- Replace the bold central question at line 33 with the question the sweep answers: how should 100 sites be placed if the spatial field is to be learned at all.
- Check the conclusions typed beside computed numbers (lines 237, 270, 345, 420, 426) against the rendered output, and compute them inline where the bundle has them.

### Lesson 3 (report: `style-diagnostic/lesson-3.md`), Phase B then C

- Restructure: diagnostics first; package plots with truth overlaid as main figures, custom figures only for comparisons; the appendix's unique guidance folded into the sections; new-site prediction and model comparison (lines 786 to 1085) split into their own lesson; the simulation-only material (trait cancellation 334 to 390, Procrustes alignment 460 to 649), the mirror study (1444 to 1514, reduced in the text to a short warning plus the per-chain check) and package history (PR #13) moved to the appendix; the reproduction section separated from the function finder. Target about 1,150 lines.
- Explain posterior draws, chains, credible intervals and log-odds at first use (line 66), hidden site factors where they are first mentioned (line 33), and ordination, variation partitioning and WAIC before they are used; show the three teaching data objects after loading.
- Traces: show the public route, `returnConvergenceDiagnostics()` and `plotTraceplot()`, in the text; reading internal slots goes to the appendix with a sentence that a per-chain accessor is planned (`TODO.md`).
- The `ggtern::theme_bw()` workaround: keep it with one sentence explaining that occJSDM registers ternary theme elements that make the plain theme fail validation (needs Alex to confirm the cause; if it can be fixed in the package, add a TODO item).
- Remove the three uncited references (Cai, Leibold, Pichler) until the internal-structure exercise exists; it is a post-beta item.
- Fix the three inconsistent statements about Lesson 2 (lines 31, 784 and 2074) to the current state now; the renumbering step rewrites the numbers.
- Cut the opening paragraph about the old `sampleresults` and "undetected" effects unless it answers a report that can be cited.
- Replace the three pipe tables (1251, 1520, 2059) with bullet lists.

### Lesson 4 (report: `style-diagnostic/lesson-4.md`), Phase B then C

- Restructure: the section 3 exercise after section 4, which it depends on; section 10 and the reproduction record to the appendix; sections 12 to 15 split into their own lesson; the full-results appendix reduced to what sections 14 and 15 do not already show; the two-optima teaching point kept in the text with its evidence in the appendix; the three pipe tables converted.
- Fix the echo leak: the `knitr::opts_chunk$set(echo = FALSE)` inside the chunk at line 1132 hides the code of every later chunk, including the "For example" chunk at 1525; scope it to the chunks that need it.
- Read the printed results in prose: the four packages' identical response curves (666), the shrinkage towards the middle in the band tables (575), the sample-size and trait results (964 to 1121), the averaging illustration (302 to 326).
- Explain the choices: why gllvm rather than occJSDM demonstrates conditional prediction (occJSDM's `predictNewSites()` has no conditional-on-another-species mode; its marginal route with the latent-factor term is `useBiotic = TRUE`, which the lesson should show alongside); why the factor count is fixed; what the section 10 arguments do (one sentence each, pointing to the dev scripts for the rest).
- State what is known about rare-species overestimation and occJSDM's 82% coefficient coverage at 100 sites: the cause is not established; the candidates are the priors, shrinkage and sample size; link to the Lesson 2 sweep's rare-species finding and the `TODO.md` item.
- Hmsc's 86 flagged fits: state the schedule used and that it was not extended (needs the dev README to confirm the schedule).
- Line 765: do not claim the weak penalty explains the unpenalised starts' disagreement unless the stability archive established it; otherwise say both explanations are open.
- Cross-reference Lesson 3's marginal-versus-conditional and WAIC teaching with a one-sentence reminder rather than re-teaching it.
- Covariate standardisation: `runOccJSDM()` standardises the occupancy covariates itself (the fit stores the standardised design matrix), so the prose at 202, 285 and 332 says so rather than asking the reader to scale first.

## What does not change

- Bundles, verifiers, committed results and the figures of Lessons 0, 1, 2 and 4, except where a figure moves or a chunk is added; Lesson 3's figures change by decision.
- The package code and defaults.
- The frozen protocol of the Lesson 2 sweep.

## Verification for every lesson pull request

- The lesson renders to HTML and to markdown with the unwrapped `github_document` call, with no new warnings; its verifier passes; `test_lesson_links.R` passes; the site test passes.
- `git grep` finds no em-dash in the changed files, and no pipe table in the Rmd.
- A reviewer reads the diff against the lesson's diagnostic report and confirms each listed gap is addressed or recorded as deferred, every factual claim in added prose is supported by the rendered output, and no verified number changed.

## Side note for Doug

The quickstart carries an `editor_options: markdown: wrap: none` block (its lines 11 to 13), added on 21 September 2026, which the project rules say never to add. The lessons do not have it and the passes will not copy it. Whether to remove it from the quickstart is Doug's call.

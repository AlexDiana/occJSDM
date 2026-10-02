# Lesson style diagnostic, 2 October 2026: summary

Five read-only diagnostics of the lessons against the writing criteria in `AGENTS.md` ("Writing the vignettes and lessons"), with the quickstart as the exemplar. One report per lesson in this folder: `lesson-0.md` to `lesson-4.md`, each with section-by-section findings and line numbers, a jargon list, structure notes, copyedit notes, a size estimate and questions for Doug. This file is the synthesis.

## Verdicts

- Lesson 0 (463 lines): mixed, large pass. Explains R mechanics more than the ecology; never lays out the two-stage process (field collection and PCR detection, with false negatives and false positives at each stage), so the detection parameters arrive as a symbol table; most outputs are printed but not read. Small structural change: one new passage, two pipe tables to bullets.
- Lesson 1 (964 lines): mixed, large pass. The middle teaches as well as the quickstart (probability versus state, priors in the reader's terms, the four detection cases). The input data are never shown with `head()`, and the lesson's two main findings (estimates pulled towards the middle; field-stage false positives mostly believed) appear in output but never in prose. Moderate structural change: four pipe tables to bullets, the fitting-reference subsection moved, the all-positive-samples table moved to the reveal section, convergence moved before results. Two facts doubted: "three optional fitting chunks" (four) and "about 800 MB" (arithmetic gives 1.6 GB).
- Lesson 2 (445 lines): mixed, large pass. The conceptual half teaches; the sweep half reads as a study report, never says what the result means for the reader's own survey, and uses the reading thresholds before defining them. Small within-lesson reordering: data loading into 2A, rule definitions before first use, reading-label table from 2B to 2C.
- Lesson 3 (2,085 lines): mixed, large pass, and the length is a structure problem. About 350 lines serve a reviewer (exporter, hashes, verifier scripts), a methods reader (axis alignment, trait cancellation, the mirror study on a different community) or package history; the main figures are custom plots from teaching summaries the reader will never have, while the package calls they will use sit in not-run chunks and a 424-line appendix that repeats the main text. Posterior draws, chains, credible intervals and log-odds are used from line 66 but explained late or never; diagnostics come last. Recommended restructure before the prose pass: diagnostics first, appendix guidance folded into the sections, simulation-only material and the mirror study moved out, reproduction split from the function finder; target about 1,330 lines, or about 1,150 if new-site prediction becomes its own lesson.
- Lesson 4 (1,532 lines, of which about 900 are code): mixed, large pass. Results are printed but not read (identical response curves, shrinkage in the band tables, the sample-size and trait results); choices lack reasoning (gllvm rather than occJSDM for the demonstrations, fixed factor count, unexplained fitting arguments); one real bug (hidden code at 1132 makes the chunk at 1525 render as output only). Required reordering: the section 3 exercise after section 4, section 10 and the reproduction record to the end. Recommended: sections 12 to 15 split into a separate lesson or article, the full-results appendix moved or cut, three pipe tables to bullets.

## Gaps common to all five

- Outputs are printed and not read: the prose rarely says what the table or figure shows and what it means.
- Concepts are used before they are named in the reader's terms: the two-stage process, posterior draws, chains, credible intervals, log-odds, hidden factors, ordination, variation partitioning, WAIC.
- The data objects the reader will handle are loaded but not shown, the step Doug added to the quickstart.
- Choices are made without their reasoning: numbers of samples and PCRs, factor counts, priors, which package demonstrates what.
- Little practical advice: no lesson ends with what the reader should do differently on their own survey, and Lesson 2 does not connect to the quickstart's advice to leave the spatial term out.
- Audience mixing: material for a reviewer or verifier (hashes, exporters, mirror studies, full-results appendices) sits inside teaching sections, mostly in Lessons 3 and 4, and R-mechanics glosses pitched below the stated reader sit in Lesson 0 and parts of Lesson 1.
- Pipe tables in the Rmd files (0: two, 1: four, 3: three, 4: three), which the repo rules say to avoid in Visual-mode files.
- Hard-coded result numbers in prose that go stale if the data are regenerated (Lessons 2, 3 and 4), and several factual claims the readers doubt (listed per report).

## Sizes

- Paragraphs needing a teaching addition: Lesson 0 about 22; Lesson 1 about 40; Lesson 2 about 26 of 45; Lesson 3 about 65 to 75 of 176, plus 35 to 45 to cut or move; Lesson 4 about 55.
- Sentences needing copyedit: 20, 35, 35, 60 to 80, 40.

## Decisions that shape the plan

- Whether Lessons 3 and 4 are restructured before their prose pass (both readers say the pass makes them longer and harder otherwise), and whether the restructure includes splitting new-site prediction (Lesson 3) and the ten-community extension (Lesson 4) into their own lessons.
- Where reviewer and verifier material goes: a separate reference article, an appendix, or the dev folder READMEs with a pointer.
- Whether the main figures should show the package's own plotting calls with a truth overlay instead of custom ggplot2 duplicates (Lesson 3 question 1).
- Sequencing against the post-beta renumbering in `TODO.md`: both touch every cross-link.
- Audience for the R-mechanics explanations in Lesson 0 and Lesson 1: keep as asides, or cut.
- Whether new run code chunks are in scope for a prose pass (several reports propose small ones, such as a `head()` of each data object or a count of each kind of detection error).
- The per-lesson questions in each report's final section (8, 14, 10, 10 and 11 questions).

## Side note

The quickstart carries an `editor_options: markdown: wrap: none` block (lines 11 to 13), which the project rules say never to add. The prose pass should not copy it into the lessons; whether to remove it from the quickstart is Doug's call.

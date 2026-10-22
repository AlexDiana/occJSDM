# Style diagnostic: Lesson 2, Spatial landscapes and survey design

Read against the criteria in `criteria.md` and the approved quickstart (`vignettes/occJSDM.Rmd`). Line numbers refer to `vignettes/occJSDM-lesson-2.Rmd` as of commit d15d198 (445 lines). Nothing was run or rendered; claims about the rendered numbers are judged from the inline R only.

Verdict in one line: the conceptual first half (lines 35-91) mostly teaches and needs light additions; the sweep half (lines 118-434) mostly describes, reads as a study report, and needs the larger pass.

## 1. What the lesson does well against the criteria

- Line 41: the two-valleys example names the spatial field in the reader's own terms (a valley that holds more of a species than its habitat predicts) before any notation, and flags its percentages as illustrative.
- Line 43: answers the eDNA reader's likely worry that the field is just a smoothed map of raw PCR detections.
- Line 45: defines range and strength plainly, and pre-empts the misreading of range as the species' geographical extent.
- Line 47: warns that the transformed distances are not kilometres and that river connectivity and downstream DNA transport are not modelled, which is exactly the question a stream-eDNA user would ask next.
- Line 53: practical caution that a smooth map or narrow interval is not evidence that extrapolation is reliable.
- Line 55: support points are named as computational anchors, not observations, with the argument that sets them and a sensitivity check.
- Lines 59-65: the prediction / association / causation distinction and the protected-valley example are well pitched; spatial confounding is named, explained by example and cited.
- Line 69: "a detection at a bridge does not automatically place an animal beside that bridge" is the kind of reader-terms concept the criteria ask for.
- Lines 71-77: concrete, actionable design advice (spread coverage, add close pairs, pilot at several separations) with the reasoning given for each.
- Line 81: explains what each replication level estimates and why none substitutes for another, plus the closure-period caveat.
- Lines 85, 89: "map pixel size is not ecological resolution" and "keep all samples and PCRs from a held-out site in the same fold" are strong practical advice.
- Line 120: gives the reasoning behind a design choice (the simulation matches occJSDM's assumptions, so it tests information, not model mismatch).
- Line 146: tells the reader what to expect before they over-read the maps (structure appears even without a spatial process; the submodel's job is only the patches).
- Line 177: "One trap is worth naming" teaches a real misconception (shrinking the study area changes nothing because axes are standardised) and links it to the arrangements.
- Line 181: the oracle is a good teaching device and is introduced with its purpose (an upper bound from the same states).
- Line 237: "a species at 5% occupancy has about five occupied sites among 100" turns a result into an intuition the reader can reuse.
- Line 241: the two gaps (estimation cost, detection cost) are named in one sentence the reader can carry through the rest of the lesson.
- Line 283: "Correlation ignores scale and the error does not: the fit shrinks the field towards zero" explains why two metrics disagree.
- Lines 418, 420: line 418 explains why the error map mirrors the field and limits the claim to one species; line 420 explicitly ties the result back to the conceptual section's coverage-versus-neighbours trade.
- Lines 430-434: exercises that use the saved tables, with exercise 1 asking the reader to apply the design statistic to their own system.

## 2. Section-by-section findings

### What this lesson adds (lines 29-33)

Describes. It announces structure and deferrals rather than telling the reader what they will be able to do or decide after the lesson.

- Line 31 (reader's next question): "the protocol and audited results are in the repository's development folder" assumes the reader has the source repository; suggest saying what the sweep is in one plain sentence (100 sites, four arrangements, simulated truth) and leaving the folder for the reproduction record.
- Line 31 (concept not in reader's terms): "same-scale environmental confounding" is used before spatial confounding is explained at line 61; suggest deferring the phrase or glossing it as an unmeasured habitat variable that varies at the same scale as the field.
- Line 33 (reader's next question): the bold central question ("why might a species still be absent") is never answered as such later; suggest either answering it in the closing section or replacing it with the question the lesson does answer (where should I put my sites, and will the spatial term help).
- Line 33 (practical advice): the lesson never tells the reader up front what it concludes for their own fits; suggest one sentence previewing the takeaway and connecting it to the quickstart's advice to leave `spatCovariates` out for now.

### Spatial effects, inference and sampling design (lines 35-91)

Mostly teaches. The gaps are links to the package (how do I do this in occJSDM?) and a few opaque sentences.

- Line 37 (choice without reasoning): says the section "can be read on its own" but does not say who should read it and who can skip to the sweep; suggest one sentence on its audience (anyone planning a survey) versus the sweep's (anyone deciding whether to fit a spatial term).
- Line 43 (concept not named): "log-odds scale" is used without a gloss; suggest a short clause saying adjustments are added before converting to a probability, so a fixed adjustment moves a 50% site more than a 95% site.
- Line 43 (concept not named): "Nearby sites are encouraged to have similar adjustments" does not say by what; suggest naming the Gaussian process (the term appears in `?runOccJSDM`) and the squared-exponential correlation that exercise 1 (line 432) relies on.
- Line 45 (concept not named): "strength" here becomes "field standard deviation" at line 120 and "amplitude" from line 181 on; suggest stating at line 45 that strength is measured as the field's standard deviation, also called its amplitude, and then using one term.
- Line 47 (practical advice): "Check that this distance model and its range grid can represent the spatial scales" gives no method; suggest telling the reader that each axis is divided by its standard deviation (confirmed in `transformCovariatesMatrix()`), so a range in their own units divided by the SD of their coordinates should land between 0.01 and 0.30.
- Line 47 (practical advice): "Enable spatial fitting by supplying two coordinate columns" contradicts nothing but omits the quickstart's recommendation (quickstart line 64) to leave `spatCovariates` out of real fits in the beta; suggest repeating it here with a forward pointer to the sweep as the evidence.
- Line 47 (choice without reasoning): the ten-value range grid is stated without saying why a grid rather than a continuous estimate, or what happens if the true range is outside it; suggest one sentence on each.
- Line 51 (incomplete example): `predictNewSites()` with coordinates is described but not shown, and line 37 points to Lesson 3 for only the non-spatial call; suggest either a short spatial call here or a sentence saying which argument carries the coordinates.
- Line 53 (concept not named): "interval coverage" and "beta-validation limitation" are package-development terms; suggest saying in plain words that it has not yet been checked whether 95% intervals contain the truth 95% of the time.
- Line 53 (concept order): "the spatial approximation" is mentioned before support points are introduced at line 55; suggest swapping the order of the two paragraphs or dropping the phrase.
- Line 55 (reader's next question): PR #8 and its merge date are development history; the reader wants to know whether to use every location as a support point, and why the sweep did; suggest replacing the PR reference with that reasoning (short-range fields need dense support points) and the cost.
- Line 55 (incomplete example): "Check sensitivity by increasing this number" has no code; suggest a one-line `listParams = list(n_supportpoints = ...)` example and naming what to compare.
- Line 65 (concept not named): "variance-partitioning plot" and "residual species correlations" assume the reader knows which functions produce them; suggest naming the functions or the lesson where they were introduced (this matters more once the lesson moves to the end).
- Line 73 (citation): Chipeta et al. has no year, unlike every other citation; suggest adding it.
- Line 75 (choice without reasoning): the 80 plus 20 split is given without saying why 20 rather than 10 or 40; suggest one sentence on the reasoning, and whether the sweep's "spread plus close pairs" arrangement is this design (lines 102-104 never say).
- Line 83 (concept not named): "calibration information or informative assumptions" is abstract; suggest one concrete example each (a positive control of known concentration; an informative prior on the false-positive rate via `listPriors`).
- Line 83 (reader's next question): the last sentence ("Computational validation with many independent binary observations sharing coordinates...") refers to a validation the reader has not met; suggest either explaining which study it means or cutting it.
- Lines 89-91 (practical advice): blocked cross-validation is recommended but the reader is not told whether occJSDM has any helper for it; suggest one sentence saying it must be done by hand (refit with `info` rows for held-out sites removed, then `predictNewSites()`), or pointing to where it is shown.
- Line 89 (concept not named): "fold" is cross-validation jargon; suggest a gloss on first use.

### 2A. One landscape, four surveys (lines 118-177)

Mixed. The framing paragraph and the trap at line 177 teach; the description of the simulated community and of the arrangements is missing.

- Line 93 (structure): the data-loading chunk sits at the end of the conceptual section, under the validation heading; suggest moving it to the start of 2A with a sentence saying what `sweep` is and whether a reader without the source repository can load it.
- Line 100 (part of the data not shown): `sweep` is a list whose parts are used throughout but never shown; following the quickstart's `str()` example, suggest a `names(sweep)` or `str(sweep, max.level = 1)` call and one line per part the lesson uses.
- Line 120 (part of the data not shown): the prose never says how many species each community has, or that species fall into 5%, 25% and 75% prevalence groups (these appear only in code at line 107 and first in prose as "common species" at line 223); suggest a sentence describing the community before the maps.
- Line 120 (choice without reasoning): the 3% range, 100 sites and three communities are stated without why; suggest one clause each (for example, short enough that the arrangements differ, the realistic budget, enough replication to read labels but not to give intervals).
- Line 120 (reader's next question): what 3% of the side would mean in a real study area is left to the reader; suggest translating it (for a 30 km region, patches about 1 km across).
- Line 146 (concept not named): "site-to-site noise" is not explained; suggest saying what it is in the simulation (an unstructured site effect, or the Bernoulli draw).
- Lines 148-158 (part of the data not shown): the four arrangements are shown as maps but the prose never says how each was built (cluster radius, number and spacing of pairs); suggest one bullet per arrangement, which would also answer line 75's question.
- Line 160 (concept not named): "correlated above 0.5" relies on a correlation function the reader has not been shown; suggest linking to the squared-exponential gloss proposed for line 43, and reconciling "about one range" here with the 1.18 ranges in exercise 1 (line 432).
- Line 172 (part of the output not shown): the column "Ratio of axis spreads" is never explained; suggest a sentence saying what it measures and why it matters (unequal axis standardisation, line 47).

### 2B. What the survey data contain (lines 179-237)

Mixed. The oracle idea teaches; the metric is not explained, the thresholds are used before they are defined, and the second half of the section reports occJSDM fits that belong in 2C.

- Line 181 (concept not named): "flat field" and the error being compared are not explained; suggest saying the error is the root-mean-square difference between estimated and true field at the sites, after centring, and that a flat field means predicting zero adjustment everywhere.
- Line 184 (choice without reasoning): the error is "centred" in the code but the prose never says why; suggest one sentence (the intercept absorbs any constant shift, so only the shape of the field is scored).
- Line 181 (concept not named): "intercept, slope, range and amplitude" uses amplitude before it is tied to strength; see line 45.
- Line 223 (concept order): "both informative thresholds" is used here but defined only at line 235; suggest moving the definition of the reading rules before this paragraph.
- Line 223 (reader's next question): why clustering raises the ceiling is implied by 2A but not restated; suggest one clause linking it to the "correlated neighbour" column of the design table.
- Lines 225-237 (structure): the reading-labels table and its two paragraphs report occJSDM's true-state fits, while the section is about the oracle and line 223 ends by saying "the next section asks whether occJSDM extracts it"; suggest moving lines 225-237 into 2C.
- Line 235 (practical advice): "Read the table rather than the prose for the result" tells the reader not to trust the prose; suggest removing it once the conditional prose at line 237 is replaced by fixed wording checked against the rendered numbers.
- Line 237 (concept not named): "the September studies" are unknown to the reader; suggest removing or replacing with what was expected and why.
- Line 237 (concept not named): "cells" (12 of them) is first used here without saying a cell is one arrangement by one species group; suggest a gloss.

### 2C. What occJSDM delivers (lines 239-345)

Describes for most of its length. The shrinkage explanation at line 283 teaches; the rest reads as a results section, and the reader is not told what any of it means for their own fit.

- Line 241 (concept not named): the three "arms" are never introduced as such in prose; suggest one bullet each, and say that the true-state fit is occJSDM fitted to the true presence/absence states, which is the JSDM-only mode the quickstart describes at its line 45.
- Line 241 (part of the data not shown): the survey arm's detection settings (collection and PCR detection probabilities, false-positive rates) are not given; suggest stating them so the reader can judge whether the detection cost is realistic for their own assays.
- Line 241 (copyedit and clarity): "two field samples, two primers and six PCRs per primer per sample, 12 PCRs per sample" reads as a list collision; see section 5.
- Line 267 (concept not named): "RMSE" first appears in a table caption; suggest glossing it at line 181 instead.
- Line 270 (reader's next question): the costs are given as RMSE increments with no sense of whether 0.0x is large; suggest comparing them to the flat-field error or to the oracle's reduction.
- Line 283 (reader's next question): the fit shrinks the field but the lesson does not say why (weak data plus a prior that pulls the amplitude towards zero); suggest one sentence, and name the prior argument in `listPriors` if one governs it.
- Line 283 (practical advice): the key practical consequence is unstated: in a reader's own fit, a small fitted amplitude does not show that spatial structure is absent; suggest saying so.
- Line 303 (concept not named): "posterior mass" and "one grid step" (the spacing of the ten-value grid, about 0.032) are not explained; suggest a gloss and the step size.
- Line 303 (part of the data not shown): "the 24 fits" is the first time the reader can infer that each fit covers one community's species jointly (4 arrangements by 2 arms by 3 communities); suggest stating the fit count and what one fit contains when the arms are introduced at line 241.
- Line 329 (concept not named): "points" means percentage points; suggest glossing on first use.
- Line 332 (reader's next question): the pull towards the middle is described but not explained, and the stronger pull in the survey arm is not linked to detection; suggest one sentence on each cause and one on what it means for a reader's occupancy maps (rare and ubiquitous species look more middling than they are).
- Lines 334-345 (choice without reasoning): "the prespecified rule for a longer run" is never stated; the table shows raw `key` and `phase` values; suggest either stating the rule or moving this material to the reproduction record.
- Line 345 (practical advice): the useful lesson for a reader (the spatial amplitude mixes slowly, so run longer chains and check its effective sample size in your own fits) is buried; suggest leading with it and saying whether `returnConvergenceDiagnostics()` reports the amplitude, since the quickstart (its line 86) says the table covers occupancy and detection coefficients.
- Line 345 (part of the output not shown): `results/field-convergence.csv` is a development-folder path the reader cannot open; suggest moving the reference to the reproduction record.

### 2D. Predicting unsurveyed locations (lines 347-420)

Mixed. The opening sentence, line 418 and line 420 teach; line 389 is a dense numerical report.

- Line 349 (part of the data not shown): "the lattice of 1,600 unsurveyed locations" is introduced without saying it is a 40 by 40 grid over the study area with known truth; suggest one clause.
- Line 389 (reader's next question): the eDNA-survey fits' positive bias is reported but not explained; suggest one sentence on its likely cause (false positives or the pull towards the middle at line 332) and its consequence for the reader's maps.
- Line 389 (structure within paragraph): the paragraph gives three results (fit versus environment-only, detection cost, oracle ceiling) in one block of nine inline numbers; suggest splitting it into three paragraphs, each led by its conclusion, so the conclusion "the environment term carries the prediction" is not buried mid-paragraph.
- Line 420 (reader's next question): the paragraph says the coverage-versus-neighbours trade does not show at this budget, but does not tell the reader whether to still follow the design advice at lines 71-77; suggest one sentence saying the advice stands for reasons the sweep did not test (a longer range, a stronger field, a larger budget) and how the reader could check their own case with `simulateOccJSDMData`.

### What this establishes, and what it does not (lines 422-434)

Describes. It summarises the sweep accurately but does not translate it into advice.

- Line 424 (concept not named): "three communities support the reading labels above, not confidence intervals" is opaque; suggest saying plainly that three replicates are enough to see whether a pattern is consistent but not to put an uncertainty interval on it.
- Line 426 (practical advice): the closing summary has no "what this means for your survey" sentence; suggest ending with the practical recommendation (leave the spatial term out in the beta unless sites are clustered within one range and species are common, and do not read a weak fitted field as no spatial structure).
- Line 428 (concept not named): "the separate contrasts still to come, with same-scale environmental confounding the first candidate" is roadmap language; suggest one sentence the reader can use, or cut.
- Lines 430-434 (practical advice): the exercises give no hints or expected answers, and the reader cannot run them without the source repository's `teaching-data` folder; suggest saying where the bundle is available and adding a one-line hint per exercise.
- Line 432 (practical advice): "Use a range that is plausible for your system" needs the conversion to the standardised scale (see line 47).

### Reproduction record (lines 436-445)

Describes, appropriately for its purpose; no teaching additions needed.

- Line 438 (part of the data not shown): `STUDY` and "the raw archive" are not available to a package user; suggest one clause saying this section is for maintainers and reviewers.
- Line 445 (structure): links name Lesson 1 and Lesson 3 and "Quickstart and lesson guide", which will be wrong once the lesson moves to the end and is renumbered; see section 7.

## 3. Jargon list

Terms used before they are explained to this reader, with the line of first use.

- same-scale environmental confounding, line 31 (spatial confounding is explained only at line 61)
- protocol and audited results, line 31
- spatial component, spatial submodel, lines 37 and 41 (explained by context; acceptable)
- log-odds scale, line 43
- spatial approximation, line 53 (support points explained at line 55)
- interval coverage, line 53
- beta-validation, line 53
- variance-partitioning plot, line 65
- residual species correlations, line 65
- calibration information, line 83
- computational validation with independent binary observations sharing coordinates, line 83
- fold, line 89
- spatial blocking, line 91
- community (a simulated set of species generated together), line 120
- field standard deviation, line 120 (not tied to strength)
- model mismatch, line 120
- site-to-site noise, line 146
- correlated above 0.5 (implies a correlation function never named), line 160
- standardised range, line 173 (explained at line 177)
- ratio of axis spreads, line 173 (never explained)
- oracle sampler, line 181 (purpose explained; "sampler" assumes MCMC)
- amplitude, line 181
- flat field, line 181
- centred RMSE, line 184 in code; RMSE first in prose at line 267
- informative, uninformative, intermediate, line 223 (defined at line 235)
- common species, rare species (prevalence groups), line 223 (groups defined only in code at line 107)
- September studies, line 237
- cells, line 237
- arm, survey arm, line 241
- true-state fit, line 241 (as "true-state")
- cost of estimation, cost of detection, line 241 (explained in the same sentence; acceptable)
- posterior median, line 283 (credible interval is in the quickstart; acceptable once glossed)
- posterior mass, grid step, line 303
- points (percentage points), line 329
- prespecified rule for a longer run, line 339
- bulk effective sample size, line 345 (ESS is in the quickstart; "bulk" is new)
- lattice, lattice predictions, line 345 in prose (line 123 in code)
- squared-exponential field, line 432

## 4. Length and structure

- Repetition: line 45's last clause ("it does not estimate a separate spatial range for each species") repeats the clause before it.
- Repetition: line 85's last sentence ("Increasing the number of support points only increases computational flexibility") repeats line 55's "More support points cannot replace missing field observations".
- Repetition: line 426 restates lines 223, 237, 283, 303 and 389 almost result by result; it could be cut to three sentences (the headline, the reason, the practical consequence) without losing teaching.
- Repetition: lines 235 and 237 restate the table they follow cell by cell; once the conditional prose is fixed, line 237 can keep only the two explanatory sentences (the 75% clustered cell, rarity as a limit).
- Could be shorter: lines 334-345 (convergence) serve a reviewer more than a learner; one or two sentences of practical advice would replace them in the lesson, with the detail moved to the reproduction record.
- Could be shorter: line 389 can lose about a third of its inline numbers if split as suggested, keeping one number per conclusion.
- Code volume: the second half shows about 200 lines of `dplyr`/`ggplot2` code that reads a bespoke bundle and calls no occJSDM function; the reader learns nothing from it about using the package. Hiding most of it (`echo = FALSE`) would shorten the visible lesson considerably; this is a decision for Doug (section 7).
- Structural changes needed for the pass to work (small, within the lesson):
  - Move the load chunk (lines 93-116) from the end of the conceptual section to the start of 2A.
  - Move the definition of the reading rules (line 235) before its first use at line 223.
  - Move the reading-labels table and its paragraphs (lines 225-237) from 2B to the start of 2C, where the true-state fits are introduced.
  - Optionally move the convergence material (lines 334-345) to the reproduction record.
- No merging or splitting of `##` sections is needed, and the lesson stays whole as decided.

## 5. Copyedit observations

- Long sentences (over about 40 words, or with three or more clauses): line 31 (second sentence), line 47 (third and fourth), line 55 (fifth), line 83 (fourth), line 120 (first and third), line 223 (second, fourth and last), line 241 (first), line 283 (fifth), line 303 (third), line 345 (third and fourth), line 389 (last sentence, about 70 words), line 420 (fourth), line 426 (third and fourth), line 428, line 432 (exercise 1), line 438 (second and fourth).
- Semicolon chains that would read better as two sentences: lines 31, 33, 43, 51, 223, 303, 345, 424, 426.
- "and/or": none found.
- Em-dashes: none found.
- Missing stops: none found; caption strings are consistent.
- Asides that repeat their sentence: line 45 (last clause), line 85 (last sentence); see section 4.
- Awkward list: line 241, "two field samples, two primers and six PCRs per primer per sample, 12 PCRs per sample".
- Inconsistent terms: strength (line 45), field standard deviation (line 120), amplitude (line 181 on); true-state fit (prose) versus "occJSDM: true states" (labels at line 106); eDNA-survey fit versus "occJSDM: eDNA survey".
- Citation format: Chipeta et al. (line 73) has no year.
- Factual claims to check:
  - Line 160 versus line 432: "within about one range, where field values are correlated above 0.5" versus 0.5 correlation at 1.18 ranges. With the kernel in `src/jsdm.cpp` (exp of minus d squared over 2 l squared), correlation at one range is about 0.61 and falls to 0.5 at 1.18 ranges, so line 432 is right and line 160 is loose; the two should agree in wording.
  - Line 237: hard-coded prose ("The exception is the clustered design for the 25% species, which is intermediate"; "its error reductions are all below 10%") sits beside two inline conditionals. If the second conditional took its first branch ("changes that for the common species"), the paragraph would contradict itself. Worth confirming which branch renders and replacing the conditionals with fixed text.
  - Line 270: "The clustered design, which has the highest ceiling, carries the largest of both" is hard-coded beside computed numbers; check it against the costs table.
  - Line 345: "every interval lies far below the true value" is hard-coded; line 283 computes the maximum upper bound, so check that it is indeed far below 1.
  - Line 420: "Adding close pairs to a spread design raised the oracle ceiling for the common species" sits uneasily with line 223, which groups pairs with spread and grid as low-ceiling arrangements; check against the oracle figure and give the number if kept.
  - Line 420: "about 3% of the error or less" is hard-coded; check against the rendered table.
  - Line 426: "the oracle clears the informative thresholds for the 25% species in every community" is hard-coded, while line 223 computes `passes(clustered_25)`; check it renders as 3 of 3.
  - Line 181: "Nothing can do better from the same states" is true on average for the posterior mean under squared error when the model is correct; worth softening to "on average" or "in expectation".
  - Line 55: the default of 20% of unique locations is confirmed by `getDefaultSupportPoints()` in `R/jsdmfun.R` (floor of 20%, capped at n minus 1); "approximately" is fine.
  - Line 47: the ten-value grid from 0.01 to 0.30 is confirmed in `R/runOccJSDM.R`; standardisation by mean and SD is confirmed in `transformCovariatesMatrix()`.

## 6. Estimate

- Paragraphs needing a teaching addition: about 26 of the roughly 45 prose paragraphs. About 9 are in the conceptual half (mostly one-sentence additions linking to package calls), about 17 in sections 2A to the closing section (several need two or three sentences).
- Sentences needing copyedit: about 35, concentrated in lines 223, 237, 283, 303, 345, 389 and 426.
- Structural moves: three small ones within the lesson (section 4).
- Size: large for this lesson. The conceptual half alone would be a small-to-medium pass; the sweep half needs new explanatory prose for the community, arrangements, arms, metric and practical consequences, plus replacement of the conditional prose.

## 7. Questions for Doug

- Who is the sweep half for? As written it serves a reviewer auditing the study (convergence rules, revision hashes, prespecified labels) as much as an ecologist learning survey design. A teaching pass would push it towards the ecologist; is that what you want, with the audit detail moved to the reproduction record and the dev folder?
- Should the plotting and wrangling code in 2A to 2D be shown? None of it calls an occJSDM function, and the reader cannot run it without the source repository. Hiding it would shorten the visible lesson by about 200 lines; showing it supports the exercises.
- Is `teaching-data/spatial-lesson.rds` going to be available to readers (installed with the package, or downloadable)? The exercises (lines 430-434) and the load chunk (line 100) depend on the answer.
- Should the lesson end with an explicit recommendation (for example, leave `spatCovariates` out in the beta unless sites are clustered within one range and the species are common), matching the quickstart's line 64? The sweep supports it, but it is your call how prescriptive to be.
- Keep or cut the line 83 sentence about computational validation with many binary observations sharing coordinates? It refers to an earlier study the reader has not met.
- Keep the PR #8 reference and date at line 55, or replace them with the reasoning for using every location as a support point?
- Should the design advice at lines 71-77 be softened, or kept with a note, now that the sweep shows no measurable benefit of close pairs at this budget and range (line 420)?
- After the move to the end of the sequence, what will this lesson be numbered, and should the title, `VignetteIndexEntry`, the Lesson 3 link at line 37 and the links at line 445 be updated in the prose pass or in a separate renumbering commit?
- The bold central question at line 33 ("why might a species still be absent") is not answered by the sweep; keep it as framing for the conceptual half, or replace it with the question the sweep answers?
- Side note, outside this lesson: the quickstart carries an `editor_options: markdown: wrap: none` block (its lines 11-13) that the project rules say never to add; lesson 2 has none, and the prose pass should not copy it across.

# Style diagnostic: Lesson 0 (Create and explore a simulated survey)

Lesson: `vignettes/occJSDM-lesson-0.Rmd`, 463 lines, read in full against the criteria in `criteria.md` and the quickstart exemplar. Read-only; nothing in the lesson was changed. Line numbers refer to the lesson as of commit d15d198.

Verdict: mixed. The lesson explains its R code carefully and gives several good cautions, but it explains the R mechanics far more than the ecology. The two-stage observation process that every setting encodes is never laid out in the reader's terms, and most outputs are printed without the prose telling the reader what to see in them.

## 1. What the lesson does well against the criteria

- Line 28 gives practical advice up front: who can skip the lesson and why someone would read it.
- Line 30 names a concept the reader will need, in plain terms: drawing a map is not the same as fitting a spatial model. Line 329 comes back to it with the reasoning (the environmental values were drawn independently of the coordinates, so points are shown rather than an interpolated surface).
- Line 73 answers the reader's next question before it is asked: six PCR replicates means six per primer per field sample, not six split between the primers.
- Line 92 (last sentence) and line 113 say plainly that the detection settings are choices made for the simulation, not findings about real primers.
- Line 124 names a subtle concept the reader would otherwise miss: an event in the simulation can produce zero reads, so the fitted `p` and `q` must be compared with the probability of a positive result at the threshold.
- Lines 179 and 455 give sound practical advice: after changing a setting, refit before comparing with truth; replacing a covariate after simulation breaks the data; a lost sample is absent rows, not zero reads and not `NA`.
- Line 212 gives practical advice that matters for real data: rows of `info` and `OTU` must stay paired, and sorting either one alone corrupts the data.
- Lines 182-202 show all three parts of the input (`info`, `OTU`, `traits`), matching the exemplar's "show every part of the data".
- Line 269 explains why the order of the `case_when()` conditions matters, with the reasoning for each case.
- Line 271 separates labels known from the simulation from anything the model produces, and says the labelled table is never passed to the fitter.
- The line 293 caption explains how the example sites were chosen (from truth and detection patterns, before fitted values were seen), which shows good practice without lecturing.
- Line 393 is the best teaching paragraph in the lesson: it tells the reader how to read the figure (left to right, row by row) and what each mismatch means.
- Line 397 explains why the removed samples were picked by a seed and not by looking at the data.
- Line 459 states the principle of the hand-off clearly: the fitter receives the survey, never the truth.

## 2. Section-by-section findings

### Where this lesson fits (lines 24-30)

Mostly describes: it is a roadmap and a status report more than an orientation.

- Line 26, missing reasoning: "a community whose true distribution is known" is the reason for the whole lesson, but the paragraph does not say why that helps. Suggest one sentence saying that known truth lets the reader check each model answer against the answer that generated the data, which a real survey never allows.
- Line 26, concept not named: the lesson never says, here or anywhere, what occJSDM does or what problem the two-stage model solves. Suggest a short sentence modelled on quickstart line 22 (false negatives from field collection or PCR, false positives from field or lab contamination), or a link to the quickstart.
- Line 26, describes instead of teaching: the long second sentence lists Lesson 2's contents and future work (oracle ceiling, spatial variation partitioning, held-out validation, dispersal contrasts). Suggest cutting it to one sentence on what Lesson 2 offers this reader, and leaving the status details to Lesson 2 or `TODO.md`.
- Line 30, the reader's next question unanswered: "the existing non-spatial simulation" reads as project history; the reader does not know what other simulation there might be. Suggest saying simply that this survey has no spatial structure and why that is a good place to start.

### Load the teaching data and the R tools (lines 32-51)

Teaches the R mechanics; leaves a practical gap and the main object unexplained.

- Line 34, practical advice missing: the teaching data are excluded from the built package (`.Rbuildignore` lists `vignettes/teaching-data`), so a reader who installed occJSDM cannot run this. Suggest a sentence telling the reader to clone or download the repository first, and where from.
- Line 42 and line 49, data not shown: `lesson` is loaded and dug into (`lesson$input$sim`, and later `lesson$cases` at line 294 and `lesson$input$source_commit` at line 179), but the reader never sees what it contains. Suggest a `str(..., max.level = 2)` call or a bullet list naming its parts.
- Line 46 and line 49, data not shown: `known_truth` is described only as "the answers"; its parts (`z_true`, `w_true`, `jsdmParams_true$eta`) first appear at lines 240-349 with no introduction. Suggest naming each one here in reader terms: which sites are truly occupied, which samples truly hold DNA, and each site's true occupancy score.
- Line 49, concept not named: "the two-stage model" is used before it is explained. Suggest a forward pointer to the new explanation proposed under line 92, or a brief gloss here.
- Line 51, concept not named in the reader's terms: "where occJSDM expects a rectangular block of numbers" is vague. Suggest saying that the fitter takes the read counts as a matrix with species as columns, which the lesson reshapes into a long table at line 218.
- Lines 49 and 51, pitched below the stated reader: the reader knows basic R and the pipe, yet the lesson explains `readRDS()`, lists, `$` and `|>`. See question 1.

### Specify the survey, one decision at a time (lines 53-88)

Partly teaches: the names are meaningful and line 73 is good, but no choice is explained.

- Line 55, missing reasoning: 100 sites, 10 species, 2 samples per site, 2 primers and 6 PCRs per primer are given without saying why. Suggest one sentence on whether this design is typical of a real eDNA survey or chosen to make the model's job feasible, and what the reader should change to match their own survey.
- Line 65 and line 73, concept not named in the reader's terms: "2,400 PCR observations per species" and the table label "PCR rows per species" suggest separate PCRs for each species. In metabarcoding, one PCR reaction yields reads for every species. Suggest saying there are 2,400 PCR reactions, each recording a read count for all ten species.
- Line 75, the reader's next question unanswered: `M` and `K` are vectors, one value per site and one per sample-primer combination, but the lesson does not say why. That is what allows unequal replication, which lines 395-455 rely on. Suggest one sentence linking the two.
- Line 75, missing reasoning: the order of `K` (sample within primer, or the reverse) is not stated. Suggest saying whether order matters here (it does not while every value is 6) and pointing to `?simulateOccJSDMData` for the case where it does.
- Lines 81-86, missing reasoning: two traits, two occupancy covariates and one collection covariate are set without saying what a collection covariate is in field terms (for example, filtered water volume or time since sampling). Suggest a short example of each kind of covariate.

### Decide how collection and PCR can fail (lines 90-124)

Mostly describes: it defines parameters by symbol before giving the reader a picture of the process they belong to.

- Line 92, concept not named in the reader's terms: this is where the reader needs the two-stage process laid out: a species occupies a site or not; a field sample captures its DNA or not (`theta`); DNA can also reach a sample from an unoccupied site (`theta0`, field contamination); each PCR detects DNA that is present (`p`) or gives a positive without it (`q`, lab contamination). Suggest a short paragraph, or a four-item list, before the code, naming each false-negative and false-positive route the way quickstart line 22 does.
- Line 92, missing reasoning: the settings make OTU_1 the hardest species to detect (`p` 0.35, lowest `theta_baseline`) and OTU_10 the easiest. The lesson never says so, although lines 333-393 map exactly these two species. Suggest one sentence here pointing ahead.
- Line 96, missing reasoning: `pmin(0.95, ...)` caps a probability that never exceeds 0.95 with these values. Suggest saying it is a guard for readers who change the first line.
- Lines 100-103, missing reasoning: `q` runs in opposite directions for the two primers, and `theta0` rises across species. Suggest one sentence saying why (for example, so the two primers' false positives do not fall on the same species).
- Lines 106-109, the reader's next question unanswered: `mu1 = 5` and `mu0 = 1.5` mean nothing until they are turned into reads. Suggest giving the typical counts (about 150 reads for a true detection, about 3 for a false positive) and the practical point that false positives tend to be low-read, which is why the threshold choice matters.
- Lines 115-122, format: a pipe table, which the repo rules forbid (Visual mode rewrites it). Suggest converting it to a bullet list as the quickstart does. The `theta_baseline` row (line 120) should also say where the covariate effect itself is set (the simulator draws it), since the reader will look for it.
- Line 124, the reader's next question unanswered: the warning about events versus results at the threshold is correct but abstract. Suggest a one-number illustration: with these settings a true event almost never gives zero reads, while roughly one false-positive event in seven does (my calculation, Doug to verify). That is why the fitted `q` should come out lower than the `q` set here.

### Specify the ecological model (lines 126-148)

Mostly describes: eleven settings with comment glosses, and no picture of the model they feed.

- Line 128, concept not named in the reader's terms: the lesson never says how occupancy probability is built. Suggest one or two sentences saying that each species' occupancy depends on its responses to the environment (shaped by its traits) plus a few hidden site factors that several species share, all combined on the logit scale.
- Line 128, concept not connected: "hidden community factor" is the latent factor that the quickstart's `n_factors` sets, and it is what creates associations between species that the covariates do not explain. Suggest naming that link, since associations between species are the point of a JSDM and the word "association" never appears in this lesson.
- Line 132, jargon: "one unmeasured trait dimension" (`gt`) is not explained. Suggest saying what it represents and why one is enough here.
- Lines 135 and 138, concept not named: `sigma_b` and `sigma_h` are glossed as "variation" and "scale" without saying what a larger value would do to the community. Suggest one clause each (for example, larger `sigma_b` makes species respond more differently to the same environment).
- Lines 136-140, incomplete example: only `sigma_bs` is marked as spatial; `sigma_ts`, `sigma_s` and `l_s` are spatial too but carry no comment. Suggest marking all four, or grouping them under one comment.
- Line 146, a choice whose reasoning is given: good. It could be one sentence shorter (see section 4).
- Line 148, missing reasoning: "We do not choose the seed to obtain an especially successful fit" is the right practice but needs its reason: choosing a seed that makes the model look good would mislead the reader about how well it works.

### Recreate the dataset if you want to (lines 150-179)

Partly teaches: the practical advice at line 179 is good, but the code's last block is unexplained.

- Line 152, jargon: "public simulator call" means an exported function, and "verified simulation" does not say verified how. Suggest plain wording, or one clause on what was checked.
- Line 163, the reader's next question unanswered: `model = "two_stage"` implies other options. Suggest one sentence naming the alternatives and pointing to the quickstart's presence/absence case.
- Lines 169-176, missing reasoning: the `dimnames()` labelling has only the comment "as in the saved teaching dataset". Suggest saying why it matters: the simulator returns unlabelled matrices, and the later joins match on site, sample and species IDs.
- Line 179, practical advice incomplete: the commit hash is shown, but the reader is not told what to do with it. Suggest one sentence on installing that revision from GitHub, and which parts of the software environment matter (R version, package versions).

### Understand the three input objects (lines 181-212)

Partly teaches: it shows every part, but the code comes without a lead-in, and two things the reader needs for their own data are missing.

- Line 183, structure: the section opens straight into code. Suggest one sentence saying what the three objects are and that the reader's own data must take this form.
- Lines 204-208, format: a second pipe table. Suggest converting it to bullets.
- Line 206, concept not named: the column-name prefixes (`X_psi.` for occupancy covariates, `X_theta.` for collection covariates, `Xs.` for coordinates) are never explained, although lines 281-284 use them. Quickstart line 41 explains them; suggest doing the same here.
- Line 208, jargon: OTU is never expanded or tied to the reader's own bioinformatic output. Suggest saying that OTU columns correspond to the reader's OTU table, and that raw read counts can be supplied as they are (the quickstart's line 42 point).
- Line 208, the reader's next question unanswered: the row names of `traits` must match the column names of `OTU`. Suggest saying so.
- Line 210, describes: the reading guide for the first six rows is good, but the closing remark on `slice_head()` is an R aside that interrupts it (see section 5).

### Follow the observations back to their known source (lines 214-271)

Teaches the mechanics well; does not say why the exercise matters or read its own output.

- Line 216, missing reasoning: "To explore individual detections" does not motivate the section. Suggest saying the aim: to see, with the answers known, the false positives and false negatives that the model will have to tell apart without them.
- Line 231 and line 235, the reader's next question unanswered: `reads >= 1` is the threshold choice, made silently. Suggest saying it matches `threshold = 1` in the fit, and linking it to line 124.
- Line 235, the reader's next question unanswered: "A missing read count remains missing" leaves the reader asking whether this dataset has any, and when a PCR result should be `NA` rather than zero. Suggest one sentence, and a forward link to the distinction at line 455.
- Line 256, concept not named: the "No detection" label lumps true negatives together with false negatives (missed detections), so the false negatives, half of the lesson's subject, never get a label. Suggest either splitting the label or saying explicitly that missed detections sit inside it.
- Line 258, concept not named in the reader's terms: "Field-stage false positive" is field contamination; "Laboratory false positive" is lab contamination. Suggest using the quickstart's words beside the labels.
- Lines 263-266, output not read: the table for OTU_1, Site 3, Sample 6 is printed but the prose never says what it shows. Line 271 calls it a "weak true-detection example" without defining weak; Lesson 1 line 657 says the sample has only two positives. Suggest one or two sentences reading the table: how many of the 12 PCRs are positive, the read counts, and why that makes the case weak.
- Data not shown: there is no survey-wide count of the source labels. Suggest a `count(source)` table so the reader sees how common each kind of error is in this survey.

### Put the sites and environmental conditions on a map (lines 273-329)

Mostly teaches; small gaps.

- Line 294, data not shown: `lesson$cases` appears with no explanation. Suggest one clause saying it holds the Lesson 1 example cases.
- Line 310 caption and line 329, the reader's next question unanswered: the scale of the covariate values (for example, standard normal) is not given. Suggest stating it once.
- Line 329: good reasoning. The forward reference to Lesson 2 will need its number changed by the post-beta renumbering (TODO line 116).

### Map probability, actual occurrence and detection separately (lines 331-393)

Teaches at the end (line 393); the middle is R mechanics, and the figure's key contrast is not drawn out.

- Line 333, missing reasoning: the two species are said to have been chosen "to connect the lessons", but not that they span the hardest and easiest to detect (see line 92). Suggest saying what the reader should expect to differ between their rows.
- Line 333, concept not named: "This deliberately simple detection rule ... is not an occupancy model" leaves the reader asking what is wrong with it. Suggest one sentence saying that "any positive" both misses sites where DNA went undetected and counts contaminated sites, so it overstates or understates occupancy depending on the species.
- Line 348 and line 379, jargon: "the generating model's scores" and "`eta`" are not named in the reader's terms. Suggest calling it the linear predictor on the logit scale, which `plogis()` (the inverse logit) turns into a probability.
- Line 379, length: a seven-sentence paragraph of function explanations, two of which repeat lines 237 and 269 (see section 4).
- Line 393, output not quantified: the paragraph says each kind of mismatch can happen but not how often. Suggest a small cross-tabulation of presence against any-positive for each of the two species, so the reader can see the false negatives and false positives in counts.
- Line 393, wording: "because of contamination" covers only the field route. A positive at an unoccupied site can also come from lab false positives (`q`); suggest naming both.

### Unequal numbers of field samples (lines 395-455)

Mostly teaches: the absent-rows point at line 455 is strong; the checks are not explained and two names are reused.

- Line 397, missing reasoning: "declared before fitting" has no stated purpose. Suggest one clause saying it prevents choosing which samples to drop after seeing results.
- Line 416, incomplete example: the removed sample and site numbers are hard-coded in prose although the chunk computes them. Suggest inline R, so the prose cannot drift from the output.
- Lines 430 and 434, the reader's next question unanswered: `samples_per_site` and `pcrs_per_primer` were numbers at lines 60 and 62 and are now overwritten with tables. A reader who reruns an earlier chunk will be confused or get errors. Suggest new names (a code change for Doug to approve), or a sentence warning about it.
- Lines 437-447, practical advice not stated: the `stopifnot()` block is the practical lesson (check the structure after every subset), but the prose never says what it does or why the reader should copy the habit. Suggest one sentence before the chunk.
- Line 418 and line 455, the reader's next question unanswered: will `runOccJSDM()` accept unequal samples per site, and do sample IDs need to be consecutive after removal? Line 416 keeps the original IDs without saying whether the fitter needs that. Suggest one sentence answering both (Doug to confirm the facts).

### Take the right objects into Lesson 1 (lines 457-463)

Line 459 teaches; line 461 describes the roadmap again.

- Line 459, jargon: "perfect-observation control" is used without definition. Suggest a clause saying it is a fit given the true presence/absence, which shows the best the model could do with perfect detection.
- Line 461, describes: future spatial work and deferred dispersal contrasts repeat lines 26, 30 and 329 and give this reader nothing to act on. Suggest cutting it to one sentence pointing to the spatial lesson, or removing it (see question 3).

## 3. Jargon list

Each term is listed at the line where it is first used before being explained to this reader.

- occJSDM (what it does): line 26, never explained in this lesson.
- oracle ceiling: line 26.
- spatial variation partitioning: line 26.
- held-out validation: line 26.
- dispersal (as a modelled process): line 26, partly explained at line 461.
- non-spatial simulation: line 30.
- knitting: line 34.
- two-stage model: line 49, never explained.
- tidyverse: line 51.
- rectangular block of numbers (the matrix the fitter expects): line 51.
- field sample: line 64.
- PCR observation: line 65.
- PCR replicate: line 73.
- simulator (`simulateOccJSDMData()`, named only at line 159): line 75.
- trait (measured): line 81.
- environmental covariate: line 85.
- collection covariate: line 86.
- true-detection event probability: line 92.
- laboratory detection event and laboratory false-positive event: lines 117-118.
- collection probability: line 120.
- log-read scale: line 121.
- threshold-positive result: line 124.
- ecological model: line 126.
- hidden community factor: line 128 (explained, but not linked to latent factors or species associations).
- unmeasured trait dimension (`gt`): line 132.
- spatial factor (`ds`): line 134.
- spatial field: line 142.
- seed: line 148.
- occupancy probability: line 148.
- public simulator call: line 152.
- package revision: line 179.
- OTU: line 208, never expanded.
- sample state: line 269.
- field-stage false positive: line 258, explained at line 269 but not in the reader's words.
- generating model's score, `eta`: line 348, explained only as "score" at line 379.
- realization: line 393.
- unequal replication: line 455.
- perfect-observation control: line 459.

## 4. Length and structure

- Repetition, roadmap: the future spatial lesson is described at lines 26, 30, 329 and 461. Lines 26 and 461 could each shrink to one sentence with no loss of teaching.
- Repetition, R mechanics: `rownames_to_column()` is explained at lines 237 and 379, and `left_join()` at lines 269 and 379. Line 379 could drop both repeats.
- Pitch of R explanations: for a reader who knows basic R and the pipe, the glosses of `readRDS()`, lists, `$`, `|>` (49-51), `seq()` (92), `rep()` (75), `rbind()` (113), `slice_head()` (210), `bind_cols()` (216), `distinct()` and `transmute()` (275), `group_by()` and `factor()` (379) take about a dozen sentences. Cutting most of them would make room for the ecological teaching without lengthening the lesson; see question 1.
- Line 146 could lose its middle sentence: the code comments already say the spatial settings are unused.
- Line 455 ends with "This is an example of preparing unequal replication, not a replicated experiment measuring the effects of sample loss." It is a defensive aside that could go.
- Structural change needed for the prose pass: small. Add one explanatory passage on the two-stage process, at the head of "Decide how collection and PCR can fail" (line 92) or as a short new section before "Specify the survey" (line 53). Convert the two pipe tables (lines 115-122 and 204-208) to bullet lists. No reordering or merging is required.
- Optional reordering (question 4): the reader meets eleven ecological settings and eight observation settings before seeing any data they produce. Moving "Understand the three input objects" ahead of the settings would let the settings explain data already seen. The current order matches the title ("Create and explore"), so this is a choice, not a requirement.
- Forward references to Lesson 2 at lines 26, 329 and 461 will all change under the post-beta renumbering in TODO line 116. If the prose pass happens first, keep those references short so the rename touches less.

## 5. Copyedit observations

- Line 26: the second sentence runs about 60 words, with a parenthetical list and a semicolon joining two unrelated clauses. Split it or cut it (see section 2).
- Line 26: three consecutive bolded sentences; bold that heavy dilutes emphasis.
- Line 30: "Mapping data and fitting a spatial model are different things." repeats the two sentences before it; keep it as the summary and shorten those two, or drop it.
- Line 92: "true-detection event probabilities" is a stacked noun phrase; "the probability that a PCR detects DNA that is present" reads more easily.
- Line 124: "positive/negative counts" is a slash construction of the "and/or" kind, and "counts" is the wrong word (they are results). Suggest "positive or negative results".
- Line 146: two long sentences joined by a semicolon; split them.
- Line 148: "We do not choose the seed to obtain an especially successful fit" can be misread as "we did not choose the seed". Reword so the meaning is clear: the seed was not picked to make the fit look good.
- Line 152: "displayed but not run when knitting" is bolded; the quickstart states the same point plainly at line 69, which is the model to follow.
- Line 179: four sentences on three topics (provenance, reproduction, refitting after a change, replacing covariates); split them into two short paragraphs or a bullet list.
- Line 210: the last sentence, on `slice_head()`, is an aside in a paragraph about reading the rows; move it or drop it.
- Line 269: "a negative PCR has no positive source to classify" is awkward; suggest "a negative PCR needs no source label".
- Line 310 caption: three sentences, the last an aside ("neither represents a real temperature or moisture measurement") that repeats line 275's point about arbitrary units.
- Line 379: seven sentences that each explain one function; break the paragraph up or cut it (section 4).
- Line 397: "The choice uses seed 3947, declared before fitting." is unclear about whose fitting and why it matters.
- Line 459: "It does not receive ..." has an ambiguous "It" after a sentence naming two fits; suggest "Neither fit receives ...".
- Line 461: "keeps the observation process recognisable" is vague.
- No "and/or" found. No missing full stops found in the prose; the code comments at lines 81-86 and 131-141 mix full sentences and fragments, which is acceptable.
- Factual doubt, line 65 and line 73: "PCR rows per species" and "2,400 PCR observations per species" are numerically right but suggest each species has its own PCRs. There are 2,400 reactions in total, each scored for all species.
- Factual doubt, line 393: "because of contamination" covers only the field route. Positives at unoccupied sites can also come from lab false positives (`q`), which line 257 classes separately.
- Factual doubt, line 416: the hard-coded removed samples (24 at Site 12, 61 at Site 31, 103 at Site 52) agree with two samples per site numbered consecutively, but they depend on the seed and the R version's sampler and are not tied to the output. Not wrong as far as I can tell, but fragile.
- Factual doubt, line 26: the description of Lesson 2's contents is a project status that will go stale, and Lesson 2's number will change under TODO line 116.
- Checked and consistent: line 124's conversion and zero-floor match `R/simulateData.R` lines 173-174; line 124's statement that the fitter uses a threshold matches `runOccJSDM()`, which requires `threshold >= 1`; line 141 (`tau` unused for binary) matches `R/jsdmfun.R` line 837, where `tau` enters only the continuous model; the 24,000 rows at line 235 and the "24 PCRs" in the line 381 caption agree with the design.
- Not verified: the glosses on `sigma_b` (line 135) and `sigma_h` (line 138) agree in broad terms with `R/jsdmfun.R` lines 751 and 774, but whether `sigma_b` is the spread around a trait-predicted mean, not the raw spread of coefficients, is for Doug to confirm before the prose explains it.

## 6. Estimate

- Paragraphs (including figure captions and the two tables) needing a teaching addition: about 22, plus one new passage on the two-stage process.
- Sentences needing a copyedit: about 20, plus about a dozen R-mechanics sentences to cut if question 1 is answered that way.
- Pass size: large for this lesson. The prose is fluent and the cautions are good, but the core ecological concept has to be added, most outputs need a reading, and the R glosses need re-pitching. Little existing text is wrong; most of the work is additions and trims.

## 7. Questions for Doug

1. Audience for the R explanations: the stated reader knows basic R and the pipe, but the lesson explains `readRDS()`, lists, `$`, `|>`, `seq()`, `rep()` and most dplyr verbs. Keep them, perhaps as a collapsible aside or a footnote for readers newer to R, or cut them so the prose can teach the ecology?
2. Lesson 0 is optional (line 28). Should the two-stage explanation live here, where the settings need it, or in Lesson 1, which every reader takes? If it goes in Lesson 1, Lesson 0 needs at least a short version and a link.
3. Roadmap prose: keep the status details of Lesson 2 at line 26 and the deferred dispersal work at line 461, or cut both to a single forward link, given the planned renumbering?
4. Order: keep "specify, then inspect", which matches the title, or show the three input objects first so the settings explain data already seen?
5. Additions with code: the report suggests a survey-wide `count(source)` table (after line 266) and a presence-by-detection cross-tabulation (after line 391). Both add chunks that are run, not just prose. Are new chunks in scope for the pass?
6. Line 256: split "No detection" into true negatives and missed detections, which changes a label that Lesson 1 may rely on, or explain the lumping in prose only?
7. Lines 430 and 434: may the pass rename the reused variables `samples_per_site` and `pcrs_per_primer`, which is a code change, or only warn about them in prose?
8. Lines 418 and 455: does `runOccJSDM()` require consecutive sample IDs after a sample is removed, and does it accept unequal PCRs per sample? The prose should answer the reader, but only you or Alex can confirm the facts.

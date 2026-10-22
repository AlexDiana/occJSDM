# Style diagnostic: Lesson 1

Lesson under diagnosis: `vignettes/occJSDM-lesson-1.Rmd` (964 lines), read in full against the five teaching criteria in `criteria.md` and the approved quickstart `vignettes/occJSDM.Rmd`. Line numbers refer to the Rmd as of commit d15d198. Rendered numbers quoted below come from reading `vignettes/occJSDM-lesson-1.md`; nothing was run or rendered.

Verdict: mixed. The conceptual core (lines 100-110, 422-434, 482-611) teaches very well, often better than the quickstart. Around it, the fitting reference, the diagnostics, the unequal-replication extension and the provenance section describe rather than teach, the reader never sees the input data, and two of the lesson's main results (the pull of occupancy estimates towards the middle, and the model believing most field-stage false positives) are shown in output but never said in prose.

## 1. What the lesson does well against the criteria

- Line 84 names the three ways an eDNA survey goes wrong in plain field terms (present but not sampled, sampled but not amplified, contamination), which is the same move Doug made in the quickstart's first paragraph.
- Line 86 states the lesson's two questions up front, in bold, so every later section can be read against them.
- Line 88 explains why simulated data are used, where truth is and is not given to the model, and that this is one dataset, not a general performance claim.
- Lines 100-110 are the best teaching in the lesson: the table separates probability, site state and sample state, line 108 says PCRs observe the sample, not the site, and line 110 answers the reader's likely confusion ("20% present") with a clear frequentist reading.
- Line 112's figure shows only truth and its caption says so, so the reader learns the probability/state distinction before any model output appears.
- Line 143 gives the reasoning for the perfect-observation control: it removes detection uncertainty but keeps the estimation problem.
- Line 145 spells out why the input is built the way it is (one row per site, select by site ID, omit `Sample` and `Primer`), which is the criterion-4 reasoning the quickstart gives at its line 45.
- Line 164 and line 182 explain individual arguments (`drop = FALSE`, `threshold = 1`, `spatCovariates = NULL`) at the point of use.
- Line 212 gives practical advice the reader needs: check the model `runOccJSDM()` prints against the survey you intended, and do not collapse replicates.
- Line 237 gives a practical warning in the reader's terms: a failed PCR is not a negative.
- Line 306 ends with practical advice: investigate missing fitted values rather than dropping them.
- Line 354 gives the reason OTU_1 and OTU_10 were chosen and disclaims selection on fit quality.
- Line 420 explains why the maps look patchy and warns against smoothing them, answering the question the maps raise.
- Lines 424, 432 and 434 explain priors in the reader's terms (good practice gives a reason to expect contamination to be rare), say what the priors do not mean, and refuse the folk rule about counting positives.
- Line 438 tells the reader to read the axis labels because the panels use different scales.
- Line 484 makes a point every eDNA user needs: at threshold one, 1 read and 1,000 reads count the same, so evidence means recurrence across replicates.
- Line 486 states the case-selection rule before any result is shown, and line 607 closes the loop by saying the field-stage case was not chosen for the model's mistake.
- Lines 601-607 interpret each case concretely with the fitted numbers, including where the model is wrong.
- Line 609 is the lesson's best practical takeaway: PCR replicates establish DNA in a tube, field replicates inform site occupancy, and neither guarantees the answer.
- Line 611 separates the conditional site-presence probability from the underlying occupancy probability with a worked number.
- Line 653 restates the collection probability as a plain question in bold.
- Lines 661 and 663 give the reasoning for the prior stress test and name its confound.
- Line 714 reports the uncomfortable result honestly (the permissive priors worsen recovery) and says the priors were not tuned.
- Line 763 points to the mirror-labelling check, which is the practical advice the TODO item asked for.
- Line 769 defines Rhat and ESS in the reader's terms and says what they do not check.
- Line 821 makes clear that a lost sample is removed rows, not zeros or `NA`.
- Line 913 is honest that one deletion is not a study of sample loss.
- Line 937 warns that "predictive" in the helper's name does not mean held-out prediction.
- Every fitting call is complete and shown (lines 166-177, 186-199, 665-679, 825-837), which meets criterion 5, and all of them already leave the spatial field out, consistent with the quickstart's advice.

## 2. Section-by-section findings

### Before you start (lines 27-80)

Describes. It sets up the files and helpers competently but shows the reader none of the data they will handle.

- Line 29 (practical advice): the long description of Lesson 2's contents does not help a Lesson 1 reader and will be wrong after the planned renumbering; the changed prose would say in one sentence what Lesson 0 offers and that this lesson is non-spatial.
- Line 33 (part of the output not shown): "posterior summaries" and "full MCMC draws" are used before either is explained; the changed prose would say, in a sentence, that a fit produces thousands of draws per quantity and that this file keeps only their averages and intervals.
- Line 43 (part of the data not shown): `survey_data` is loaded but never displayed; the added prose and chunk would show `head()` of its `info`, `OTU` and `traits`, exactly as Doug added to the quickstart, and say which `info` columns are occupancy, collection and coordinate columns.
- Line 43 (reader's next question): the reader who has just done the quickstart will ask why this lesson does not use `sampledata`; the added sentence would say this survey was simulated with its truth kept, so the fit can be checked (line 960 says this only at the very end).
- Line 41 (part of the data not shown): the `lesson` object has `input`, `cells`, `rates`, `cases`, `observations`, `samples`, `diagnostics` and `manifests`, each used later without introduction; a one-line `str(lesson, max.level = 1)` with one sentence per component would orient the reader.
- Line 48 (concept not named): "generating occupancy probability", "posterior mean" and "95% credible interval" appear here first; the changed prose would say "the true probability the simulation used", and give a one-clause meaning for posterior mean and credible interval.
- Line 56 (practical advice): explaining `select()` and `slice_head()` is below this reader's level; the prose could drop it and keep only the reason for the formatting helpers.

### What are we trying to learn? (lines 82-98)

Teaches. Strong framing; two gaps.

- Line 84 (concept not named): the three failure routes are described but not named as false negatives and false positives, the terms the quickstart now uses; the changed prose would attach those two names so the reader can carry them across lessons.
- Line 86 (concept not named): "occJSDM connects these stages" is vague; the changed prose would say the model gives each stage its own probability (occupancy, collection, PCR detection, and the two false-positive rates) and estimates them together, naming where each is met later.
- Line 93 (concept not named): "two hidden community factors describing additional variation among sites" is opaque to someone who has not fitted a JSDM; the changed prose would say these let species co-occur more or less than the measured environment predicts, as the quickstart puts it.
- Line 95 (reasoning behind a choice): two samples by two primers by six PCRs differs from the quickstart's three samples by three primers by two PCRs, with no reason; one sentence would say why this design was chosen, or that it is arbitrary.

### Three questions, three different truths (lines 100-139)

Teaches. The lesson's clearest section.

- Line 106 (concept not named): "collection failure and field-stage contamination" appear before either is explained; the changed prose would gloss both in a clause each, since the reader meets them again only at lines 430 and 613.
- Lines 102-106 (format): this is a pipe table, which the repo's markdown rules forbid; the pass would rewrite it as a three-item bullet list (see section 4).

### First give the JSDM perfect observations (lines 141-178)

Mostly teaches.

- Line 141 (concept not named): "JSDM" is used in the heading and line 143 but never expanded in this lesson; the changed prose would expand it and say in a clause what it models.
- Line 143 (reader's next question): the reader will ask why anyone would fit to perfect observations; the added sentence would say it is a benchmark that shows how much error comes from the ecological model alone, so the next fit's extra error can be attributed to detection.
- Line 174 (reasoning behind a choice): `n_lattrait = 1` is set with no explanation until line 184 and line 224; the changed prose would say in a sentence what an unmeasured trait dimension is (species similarity in response not captured by the measured traits).
- Line 176 (reasoning behind a choice): the MCMC settings are first used here but explained only in the next section; the explanation at line 184 should move to the first call, with a reason for four chains and these lengths rather than the defaults used in the quickstart.
- Lines 171-177 (reader's next question): the call omits `collCovariates` and `threshold`; one sentence would say why (with one row per site there are no detection stages, so neither applies).

### Now give occJSDM only the PCR observations (lines 180-200)

Mixed: settings are explained, but as a list of facts rather than reasoning.

- Line 182 (copyedit and reasoning): seven short sentences alternate between what the model must infer and individual argument notes; the changed prose would first say what is different from the perfect fit and why it is harder, then explain the arguments.
- Line 184 (practical advice): "The two community factors and one unmeasured trait dimension match the simulation" leaves the reader asking how to choose `n_factors` and `n_lattrait` for real data, where the truth is unknown; the added sentence would give the practical rule, or point to where it is covered.
- Line 184 (concept not named): "chains", "burn-in", "iterations" and thinning are named by their settings but a chain is never explained; one sentence would say that each chain is an independent run of the sampler and that agreement between chains is checked later.
- Line 194 (part of the data not shown): `collCovariates = "X_theta"` names a column the reader has never seen; this is fixed by showing `survey_data$info` at line 43.

### A short fitting reference for your own data (lines 202-260)

Describes. Useful reference material, but it is a lookup table interrupting the narrative, and the recovery figure that answers the lesson's first question sits inside it.

- Lines 206-210 (concept not named): "one-stage occupancy model" is introduced only as a table row; one sentence would say it is for surveys with repeated field visits or samples but a single detection step, with no PCR replicates.
- Line 212 (copyedit and reader's next question): the paragraph mixes seven topics; the last sentence leaves the reader with a one-row-per-site OTU table asking what to do, and the changed prose would say plainly that they must convert counts to 0/1 themselves, because the threshold is applied only when there are replicates (the code stops with "Counts model not supported yet", `R/runOccJSDM.R` line 230).
- Line 223 (concept not named): "residual species associations" needs the plain reading suggested for line 93.
- Line 226 (concept not named): "Prior settings" appears before priors are explained at line 424; either move the explanation earlier or add a forward pointer.
- Line 229 (concept not named): "design-matrix column names" and "category contrasts" are statistics jargon; the changed prose would say that a categorical covariate is split into one column per level and these new names are what the coefficient functions expect.
- Line 232 (part of the output not shown): the `colnames(fit$X_psi)` chunk is not run, so the reader never sees what the names look like; one example line of output in the prose would do.
- Line 237 (reasoning behind a choice): `NA` is supported only in the two-stage model, without saying what a one-stage user should do instead; one sentence would give the workaround or say there is none.
- Line 239 (copyedit and reader's next question): the memory paragraph is long and refers to "the original walkthrough" the reader has never seen; the changed prose would say when a reader needs draws (to get intervals for site occupancy, which line 958 also needs) and then give the memory cost.
- Line 241 (structure): the recovery figure, the lesson's first result, appears under the reference subheading; see section 4.
- Line 260 (part of the output not shown): the prose explains the plotting code but not what the figure shows; the added prose would name the pattern the reader should see, that low probabilities are overestimated and high ones underestimated in both fits, more so with PCR data, and that this pull towards the middle is the lesson's first answer.

### Calculate the errors ourselves (lines 262-350)

Mixed: the arithmetic is well explained, the result is not interpreted.

- Line 264 (structure): the code runs before absolute error is defined; the definition at line 310 should come before the chunk.
- Line 270 (reasoning behind a choice): the 20% and 80% band edges are not justified; one clause would say why these cut points (for example, they separate rare and near-certain occupancy, where the pull to the middle shows most).
- Line 306 (practical advice): explaining `abs()`, `group_by()` and `bind_rows()` is below this reader's level; the paragraph could keep only "each species-site pair is weighted equally" and the missing-values advice.
- Line 308 (reader's next question): "does not explain all of it" leaves the reader asking where the remaining 11 points come from even with perfect data; the added sentence would say that 100 binary outcomes carry little information about each probability and the prior pulls estimates towards the middle (the reason given only at line 480).
- Line 330 (part of the output not shown): the rendered numbers (low band: truth 6.5%, estimate 20.7%; high band: 91.2% and 69.5%) are quoted but "this pattern" is never named; the changed prose would name it as shrinkage towards the middle, give its size in points, and say what it means for a reader estimating a rare species' occupancy.
- Line 332 (reader's next question): the interval figure is shown without saying how often the reader should expect 95% intervals to contain the truth, or how often they do here; one or two sentences would say what to look for and point to the README's interval-coverage limitation.
- Line 350 (concept not named): "fitted-site estimates" and "inferred using the observations at these sites" need one plain sentence: these are estimates for sites the model saw, which is easier than predicting a new site.

### Where do the errors occur on the map? (lines 352-420)

Mixed: good reasoning for choices, no interpretation of what the maps show, and no stated reason to map at all.

- Line 354 (reasoning behind a choice): the section never says why an ecologist would map errors; the added sentence would say that clustered errors are the sign of missing spatial structure, which is what to look for in one's own fit.
- Line 356 (practical advice): the advice about site ID types is useful; the `left_join()` explanation is below the reader's level.
- Line 383 (practical advice): the `pivot_longer()` explanation is fine; its reason (one shared colour scale) is the teaching point and should lead.
- Line 398 (part of the output not shown): neither map is interpreted; one or two sentences would say what the reader sees (no geographical clustering of errors, consistent with line 420) so that line 420 answers an observation rather than predicting one.
- Line 420 (practical advice): the last sentence is roadmap material about Lesson 2 and dispersal; it would become one pointer after the renumbering.

### How does good practice enter the model? (lines 422-480)

Teaches in the opening, describes in the second half.

- Line 426 (concept not named): the table lists `p`, `q` and `theta0` but not `theta`, the collection probability, which completes the chain of stages and appears unannounced at line 613; a fourth item would close the gap.
- Line 428 (concept not named): "Beta(5, 1)" notation is never explained; one sentence would say a Beta distribution describes a probability between 0 and 1 and that the two numbers act like prior successes and failures, so Beta(1, 20) behaves like having seen one contamination in 21 tries.
- Line 432 (reader's next question): the reader will ask when to change these defaults for their own lab; the added prose would say what kind of evidence (for example, contamination seen in extraction blanks or PCR negatives) would justify a different prior, and point to `listPriors` and the stress test at line 659.
- Line 436 (copyedit and structure): a sentence ending in a colon is followed by another paragraph, not the figure; merge it into line 438.
- Line 440 (part of the output not shown): the rate-recovery figure is never interpreted; the added prose would say which rates are recovered well and which are not (for example, whether the small false-positive rates sit within their intervals).
- Lines 466-478 (reasoning behind a choice): the threshold-conversion detail explains why the black crosses differ from the simulator's nominal rates, which matters only to someone comparing with the simulator settings; see question 2.
- Line 480 (concept not named and structure): "logit scale" and "MCMC mix worse" are jargon, and the paragraph's subject (the occupancy-baseline prior) is not about good practice or detection; it would read better as its own short subsection that says in reader terms what the default prior does to rare and common species, which is also the explanation line 308 needs.

### Four examples: inspect the observations first (lines 482-541)

Teaches.

- Line 484 (part of the data not shown): the reader has not yet seen what `lesson$observations` holds; one sentence (or a `head()`) would say it is one row per species, sample, primer and PCR with reads and the 0/1 result.
- Line 486 (concept not named): "focal sample" and "genuine positives" are introduced in passing; the changed prose would define focal sample as the sample whose PCR pattern placed the case in its category, and say "true detections".
- Line 486 (reasoning behind a choice): "before looking at fitted probabilities" is stated but its reason (so the examples cannot be picked to flatter the model) appears only at line 607; move the reason here.
- Line 512 (reader's next question): the heading invites inspection but never asks the reader to judge; one sentence would ask the reader to decide which of the four patterns they would believe before revealing the truth, which is the section's teaching device.
- Line 534 (part of the output not shown): the case table's "Eligible samples" column (11, 541, 260, 44 in the rendered output) is never explained; one sentence would say it counts the samples that met each category's rule, which also tells the reader how common each situation is.

### Reveal the truth and compare it with the fit (lines 543-657)

Teaches. The case paragraphs are the model for the rest of the lesson; the gaps are about carrying the lesson to the reader's own data.

- Line 543 (structure): the heading is followed directly by a chunk with no sentence; one lead-in sentence would say what the colours now add.
- Line 582 (practical advice): the reader is never told where these site and sample probabilities come from in their own fit; one sentence would name the fitted object's `z_output` and `w_output` (posterior means by default) as the source.
- Line 594 (concept not named): the table says "Estimated chance site was occupied", and line 611 later calls it a "conditional site-presence probability"; introduce the term with the table, since line 611 depends on it.
- Line 609 (reader's next question): "could be harder still if that dependence is not represented by the model" leaves the reader unsure whether it is; the changed prose would say plainly that the model treats contamination events as independent, so shared contamination across samples or batches is not modelled.
- Line 613 (concept not named and structure): the collection probability is used here before line 653 defines it; move the bold question from line 653 to the start of this paragraph.
- Lines 640-650 (part of the output not shown): for the field-stage case the table shows a "True collection probability" of 75.0% (rendered output), which line 653 then says does not apply because the species was absent; the prose or the table would need to make that row's meaning clear before the reader misreads it.
- Line 655 (practical advice): the `plogis()` explanation is good; it would help to say once that the same logistic conversion links every covariate score to a probability in the model.
- Section length (structure): at 115 lines this section would read better split at line 613 into a subsection on collection conditions.

### What changes if we are less confident about low contamination? (lines 659-763)

Mixed: the stress test is well reasoned; the dataset-wide table that follows is unexplained and misplaced.

- Line 661 (reasoning behind a choice): Beta(1, 4) is not justified beyond "less confident"; one clause would say why mean 20% was chosen as the stress level.
- Line 677 (concept not named): the reader is not told that `a_q`, `b_q`, `a_theta0` and `b_theta0` are the two numbers of the Beta priors from line 428; one sentence would map them.
- Line 678 (reasoning behind a choice): the doubled chain lengths are explained only at line 767; a forward pointer here would answer the obvious question.
- Line 714 (practical advice): the result is reported, but the reader is not told what to do with real data; the added sentence would recommend a refit under a more permissive prior as a routine check of whether key conclusions depend on the contamination assumption.
- Line 716 (structure): the all-positive-samples table uses only the default fit, so it belongs with the case analysis, not with the prior comparison; see section 4.
- Line 761 (part of the output not shown): the table is never interpreted; the rendered output shows laboratory false positives are well handled (541 samples, mean fitted DNA probability 2.2%) but field-stage false positives are mostly believed (44 samples, mean fitted site probability 67.6% at sites that were all unoccupied); saying this is the dataset-wide answer to the lesson's second question, and is the most important missing sentence in the lesson.

### Are the calculations stable enough to interpret? (lines 765-799)

Mixed: Rhat and ESS are well defined; the placement and the link to the reader's own fit are not.

- Line 765 (structure): convergence is checked after every result has been interpreted, the opposite of the quickstart's "Before reading any estimate, check that the chains have converged"; see section 4.
- Line 767 (practical advice): the extended schedule is the practical remedy, but it is told as history; the changed prose would say explicitly that a slow-mixing parameter is fixed by a longer run, and that "Both versions are retained in the build archive" is for maintainers.
- Line 771 (practical advice): the reader is not told which function gives these diagnostics for their own fit; one sentence would name `returnConvergenceDiagnostics()` as in the quickstart.
- Line 777 (reasoning behind a choice): the 1.01 and 400 screens are called common but not justified; one clause each would do.
- Line 797 (reader's next question): the reader will ask why latent-factor coordinates are not checked; one sentence would say individual factors can swap or flip between chains without changing the fit, so their combined effect is checked instead.

### Fit a survey with unequal replication (lines 801-933)

Describes. It demonstrates correct input formatting, then repeats the full truth-comparison and diagnostic workflow at length for a result it says cannot be generalised.

- Line 803 (reader's next question): the section does not say why it is here; one opening sentence would say that real surveys lose samples and this shows how to supply them.
- Line 835 (copyedit): `2L`, `3000L` and `listPriors = list()` differ in style from the three earlier calls with no reason; make the calls consistent.
- Line 913 (practical advice): the result is stated, but the practical message (drop the rows, keep the IDs, nothing else changes) is spread over lines 237 and 821; the closing would restate it.
- Line 915 (concept not named): "public diagnostic table" and "R warning conditions" are internal language; the changed prose would say "the table from `returnConvergenceDiagnostics()`" and drop the warning remark.
- Line 933 (reader's next question): two Rhat implementations (`coda` and `posterior`) are reported without telling the reader which to use; one sentence would say the reader should use the package's function, and that the second check is the lesson's own extra.

### Reproduce the lesson and inspect its evidence (lines 935-964)

Describes. The extraction example is useful; the rest is provenance, and the lesson ends without saying what the reader has learned.

- Line 937 (practical advice): the chunk is good; it would be stronger with one sentence on how to get intervals for one's own fit, linking line 958's draws requirement back to `summarisedLatentPresences = FALSE` at line 239.
- Line 958 (concept not named): "applying a probability conversion to an average coefficient" is precise but opaque; one sentence would say why the order matters (the probability of the average is not the average of the probabilities).
- Lines 960-962 (practical advice): hashes, commit, seed and session details serve a reproducibility reader, not the ecologist; see question 5.
- Line 964 (practical advice): the closing lists what the lesson does not do and repeats line 29's Lesson 2 summary, but never states what the reader learned; the changed prose would give three or four takeaways (estimates are pulled towards the middle, laboratory false positives are well handled, field-stage false positives are often believed, field replicates matter for site occupancy).

## 3. Jargon list

Terms used before they are explained to this reader, with the line of first use.

- "oracle ceiling", "site-arrangement sweep", "variation partitioning", "held-out validation": line 29 (never explained in this lesson).
- "posterior summaries", "MCMC draws": line 33.
- "generating occupancy probability" (and "generating" throughout): line 48.
- "posterior mean", "95% credible interval": line 48.
- "default priors" (in the fit labels): line 61; priors are explained at line 424.
- "hidden community factors": line 93; never given an ecological reading.
- "realization": line 110.
- "JSDM": line 141; never expanded in this lesson.
- "binary data": line 143.
- `n_lattrait`, "unmeasured trait dimension": lines 174 and 184; defined only by its table row at line 224.
- "chains", "burn-in", "iterations", thinning: line 176 (code) and line 184; a chain is never explained.
- "one-stage occupancy model", "two-stage occupancy model", "pure JSDM": lines 206-210.
- "residual species associations": line 223.
- "design-matrix", "category contrasts": line 229.
- "latent quantity", "mixing": line 239.
- "fitted-site estimates": line 350.
- "Beta(5, 1)" and Beta notation: line 428.
- "threshold-adjusted generating rates": line 438.
- "nominal event probability", "log-read value": lines 466 and 468.
- "logit scale", "MCMC mix": line 480.
- "focal sample", "genuine positives": line 486.
- "conditional site-presence probabilities": line 611.
- "collection probability": line 613; defined at line 653.
- "screen" (Rhat 1.01, ESS 400): lines 777 and 799.
- "latent-factor coordinate", "derived quantity": line 797.
- "Monte Carlo uncertainty": line 913.
- "public diagnostic table", "R warning conditions": line 915.
- "classical `coda` Rhat", "`posterior` diagnostics", "reconstructed probability draws": line 933.
- "source hashes", "build archive": lines 767 and 960.
- "Paper2Agent interface", "beta-release error target": line 964.

## 4. Length and structure

- Repetition: line 29 and line 964 describe Lesson 2 in nearly the same words; line 486 and line 512 both say the case labels come from truth, not the model; line 260, line 350 and the caption at line 241 all say truth was not supplied to the fits; line 892 restates the signed and absolute error definitions from line 310.
- Repetition: the unequal-replication section (lines 840-933) repeats the truth scatter, interval, error-table and diagnostic workflow already taught at lines 241-350 and 765-799.
- Could be shorter: the dplyr and ggplot verb explanations at lines 56, 260, 306, 356 and 488 are below a reader who knows the tidyverse; the non-obvious ones (`pivot_longer()` at 383, the join key at 488, ID types at 356) are worth keeping.
- Could be shorter: lines 466-478 (threshold conversion), lines 801-933 (unequal replication, which could be cut to data preparation, the fit call, one figure and one paragraph), and lines 960-962 (provenance).
- Structural change needed (repo rule): four pipe tables in prose, at lines 102-106, 206-210, 216-227 and 426-430, break the repo's "no pipe tables" rule; the separator rows differ in style (lines 103 and 427 use `|----|`, lines 207 and 217 use `|---|`), which suggests two have already been canonicalised by Visual mode and two have not, so the latter will churn. The pass must convert all four to bullet lists.
- Structural change needed: the recovery figure at line 241 sits under the subheading "A short fitting reference for your own data" (line 202), so the answer to the lesson's first question is filed under a reference table. Either move lines 202-239 later (after the results, or to the end as a reference section) or give the figure its own heading.
- Structural change needed: the all-positive-samples table (lines 716-761) uses only the default fit and belongs at the end of the reveal section (after line 657), where it answers the second question across the dataset; the prior-comparison section would then hold only the stress test.
- Structural change recommended: convergence (lines 765-799) comes after all interpretation, contradicting the quickstart's teaching order; move it before "Calculate the errors ourselves" (line 262), or at least add a forward pointer there. Moving it is cleaner because line 767 also explains the alternative fit's longer run, which line 678 needs.
- Structural change recommended: absolute error is defined at line 310, after the chunk that computes it (line 264) and the sentence that reports it (line 308); reorder within the section.
- Structural change recommended: the occupancy-baseline prior paragraph (line 480) would become its own short subsection, and ideally sit near the error results, since it explains the pull towards the middle.
- Dependency: lines 29, 212, 239, 420, 763 and 964 cross-link Lessons 2 and 3, and the planned renumbering (TODO line 116) will change every one; the prose pass and the renumbering touch the same lines.

## 5. Copyedit observations

- Line 29: one sentence of about 55 words with a parenthetical list and two semicolons.
- Line 33: "The three optional fitting chunks" is wrong; there are four fitting chunks not run during knitting (lines 166, 186, 665, 825), plus two other unevaluated chunks (lines 231 and 939).
- Line 33: "The saved file contains the complete simulation" is incomplete; the unequal-replication section loads a second file (line 806).
- Line 182: an aside ("This fitting chunk is also optional.") interrupts a run of argument notes.
- Line 212: seven sentences on seven topics; split into model choice, sample IDs, and unsupported inputs.
- Line 212: "A one-row-per-site matrix of integer abundances greater than one is not currently a supported count-data JSDM" is awkward; a matrix is not a JSDM.
- Line 239: eight sentences, one of them a parenthetical correction of the help page; split by question (what is kept, what it costs, what thinning does).
- Line 239: "about 800 MB" is doubtful. 500 sites times 100 species times 4,000 iterations times 8 bytes is about 1.6 GB for `psi_output` alone, and `z_output` is kept at the same size; unless the walkthrough's figure assumed 2,000 iterations, the number is half what the arithmetic gives.
- Line 239 is correct against the code (only `z_output` and `psi_output` keep draws, `R/runOccJSDM.R` lines 954-962 and 1448-1455), but `?runOccJSDM` (`R/runOccJSDM.R` lines 315-322) says all four arrays keep draws when the flag is `FALSE`; the help page, not the lesson, is wrong. Outside this lesson, but worth a TODO.
- Line 310: "This is an arithmetic illustration; the reported averages come from the simulation." is an aside that repeats its own example; remove it.
- Line 385: the caption's "0–100%" uses an en dash in a string literal; it is data, so it may stay, though "0 to 100%" would match the rest of the prose.
- Line 420: "contrasting dispersal processes are deferred to later work" is roadmap, not teaching.
- Line 436: a sentence ending in a colon followed by a paragraph instead of the figure.
- Line 480: about ten sentences; "pure JSDM fits to binary data" appears three times; "clearly in spatial fits and only slightly in non-spatial ones" dangles after its parenthetical.
- Line 488: "`semi_join()` would keep matching rows without adding columns" explains a function that is not used; remove it.
- Lines 601, 603 and 607: "gives site presence X probability" and "assigns site presence X probability" are awkward; "a X probability that the site was occupied" reads better.
- Line 609: "could be harder still" is vague; see the section 2 finding.
- Line 653: "the site was absent" should be "the species was absent from the site".
- Line 767: "Both versions are retained in the build archive" is a maintainer aside.
- Line 797: "detection/error rates" uses a slash for "and"; "detection and false-positive rates" is clearer.
- Lines 210 and 212: "PCR/primer" uses a slash for a relation; "PCR replicates within primers" is clearer.
- Line 835: `2L`, `3000L` and `listPriors = list()` are inconsistent with the other three calls.
- Line 913: "about 17 percentage points" is hard-coded while every other number is inline R; it matches the rendered 17.1 now but could drift on a rebuild.
- Line 915: "The saved call raised no R warning conditions; the diagnostic flag is a separate issue." is an aside the reader does not need.
- Line 962: one long sentence with a parenthetical "(and likewise for the other fits)".
- Line 964: the link text "Quickstart and lesson guide" is stale, since the quickstart was rewritten and no longer guides the lessons; "Paper2Agent interface" and "beta-release error target" are internal terms.
- No "and/or" occurs in the lesson. No em-dash occurs in prose.
- Checked and correct: the prior means at lines 428-430 (5/6, 1/21); the baseline ranges at line 480 (the logistic of plus or minus 1.96 and 3.92); the counts at line 98 (200 samples, 2,400 PCRs) and line 803 (197 samples, 2,364 PCRs); the threshold conversion at lines 466-468 against the simulator (`R/simulateData.R` lines 164-175: false-positive log reads Normal(1.5, 1), counts rounded from the exponential minus one); the model-selection rules at lines 206-210 against `R/runOccJSDM.R` lines 150-217; `NA` support at line 237 against `R/runOccJSDM.R` line 585; "average covariate values" at line 480, given the occupancy covariates are standardised.

## 6. Estimate

- Paragraphs needing a teaching addition or change: about 40, of which about 10 are substantive (show the input data and the `lesson` object at line 43; name the pull to the middle at lines 260, 308 and 330; interpret the dataset-wide positives table at line 761; interpret the rate figure at line 440; explain the hidden factors and latent traits at lines 93 and 174; add practical advice on choosing priors and factors at lines 184 and 432; add closing takeaways at line 964), and the rest are one sentence or one clause each.
- Sentences needing copyedit: about 35, plus converting four pipe tables to bullet lists.
- Size: large. Lesson 1 is nine times the quickstart's length, needs real additions in most sections, and needs at least two structural moves (the reference subsection and the positives table) plus the table conversions; a convergence move and trims to three sections would follow from Doug's answers below.

## 7. Questions for Doug

- Audience of lines 202-239: is the fitting reference for the first-time reader of this lesson, or reference material that belongs at the end of the lesson, in the quickstart, or in `?runOccJSDM`?
- Lines 466-478: keep the threshold-conversion explanation, collapse it to one sentence, or move it to Lesson 0 where the simulator is described?
- Lines 352-420: the maps do not directly serve either of the lesson's two questions; keep, shorten, or drop?
- Lines 801-933: keep the unequal-replication section with its full truth comparison and two diagnostic systems, or cut it to a formatting demonstration of about 20 lines?
- Lines 960-962: keep the provenance details for reproducibility readers, or move them to the build README with a one-line pointer?
- Convergence order: move "Are the calculations stable enough to interpret?" before the results to match the quickstart, or keep it last as a separate lesson about numerical checks?
- Tidyverse explanations: the reader is assumed to know the pipe and dplyr; cut the routine verb explanations (lines 56, 260, 306, 356, 488) and keep only the non-obvious ones?
- Practical prior advice at line 432: can the lesson recommend a basis for changing the contamination priors (for example, contamination rates seen in blanks and negatives)? This is a substantive recommendation Alex may want to endorse.
- How to choose `n_factors` and `n_lattrait` on real data (line 184): is there advice the package authors want to give, or should the lesson say plainly that there is no rule yet?
- The pond case study (TODO line 120) proposes Lesson 1 as the place to contrast pre-filtering with modelling detection; should this pass leave room for it, or is it a later addition?
- Sequencing: the prose pass and the Lesson 2 renumbering (TODO line 116, after the beta tag and PR #14) touch the same cross-links at lines 29, 212, 239, 420, 763 and 964; which comes first?
- Data design: should the lesson explain why it uses two samples, two primers and six PCRs when the quickstart's `sampledata` uses three, three and two (line 95), or leave it?
- Lines 640-653: should the case table stop showing a "True collection probability" for the field-stage case, where line 653 says it does not apply, or keep it with a clearer note?
- Line 964: remove the "Paper2Agent" and "beta-release error target" disclaimers, which no reader of the lesson will recognise?

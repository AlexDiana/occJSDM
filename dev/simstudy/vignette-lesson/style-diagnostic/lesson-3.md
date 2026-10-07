# Style diagnostic: Lesson 3 (outputs compared with truth)

Lesson under diagnosis: `vignettes/occJSDM-lesson-3.Rmd`, 2,085 lines, about 176 prose paragraphs and 10,000 words of prose. Criteria: `.superpowers/style-diagnostic/criteria.md`. Exemplar: `vignettes/occJSDM.Rmd`. I read the lesson in full and in order, in six passes (lines 1-240, 240-540, 540-800, 800-1090, 1090-1400, 1400-1720 and 1720-2085), after the criteria and the quickstart. Nothing was run or rendered. I checked four claims against the package source (cited in section 5).

Reader assumed throughout: an ecologist who runs eDNA or presence/absence surveys, knows basic R and the pipe, and has fitted neither an occupancy model nor a JSDM.

## 1. What the lesson does well against the criteria

- Line 66 sets one visual convention for the whole lesson (black is truth, blue is estimate) and immediately names the concept the reader most needs: an interval is not a measure of how far the estimate is from truth.
- Line 132 answers the question every reader asks of a coefficient plot: what positive, negative and zero mean, and that the scale is log-odds per standard deviation, not percentage points.
- Line 134 gives a three-step reading recipe for each species and the two standard misreadings (an interval crossing zero does not mean no effect; excluding zero does not mean an accurate size). This is the exemplar's "practical advice" criterion at its best.
- Line 179 spells out the reasoning behind why one coefficient gives different probability changes for different species, in plain terms.
- Line 183 says what the baseline probability is for and what it is not (not landscape-average occupancy).
- Line 209 explains why trait effects are hard in the reader's terms: 100 sites per species response, but only ten species to compare, and more PCRs do not add species.
- Lines 309 and 332 teach a reusable tidyverse recipe for summarising any posterior array, and then explain each verb (`expand_grid()`, `rowwise()`). This is the most transferable teaching in the lesson.
- Line 426 explains a display choice (grey, not white, for undefined correlations) and why the other choice would mislead.
- Line 585 defines a loading in the reader's terms and blocks the causal misreading.
- Line 657 warns that the `Biotic` label does not mean biotic interactions were identified.
- Lines 725-737 and 758 are a model of the style: a restricted, reader-shaped question, the arithmetic as four lines of R the reader can run, and the ecological reason the one-sample curve levels off.
- Line 764 interprets a matrix entry concretely ("0.8 means an 80% posterior probability that the species was present").
- Line 780 answers "which one should a study report?" with a published precedent and a reason.
- Lines 788 and 812-817 frame new-site prediction as the reader's own situation, then name the two target probabilities (conditional and marginal) in plain terms. Lines 819-860 then work one site through with `plogis()`, `dnorm()` and `integrate()`, each explained.
- Lines 952 and 974 define signed error, absolute error, Brier score and log score in a sentence each.
- Line 1088 frames the latent-presence table as a practical question, and line 1161 reads two specific rows aloud, which is exactly "show the reader what they will see".
- Line 1130 explains why truth is joined by identifiers, not position, and what `relationship = "one-to-one"` protects against.
- Line 1233 explains MCMC chains in plain language and gives a four-step practical sequence.
- Line 1237 explains why the package's console thresholds are loose and when to use the stricter table.
- Lines 1285, 1347-1349 and 1590-1592 give practical advice of the kind Doug asked for: missing diagnostics are flagged, what a bad trace looks like, do not tune the sampler to hit truth, keep `nthin = 1`.
- Lines 1512-1514 give an honest, hedged procedure for the mirror case and say plainly what is interpretation rather than tested procedure.
- Lines 2004 and 2018 teach `patchwork` operators and the parenthesis trap.
- Throughout, every displayed figure built from saved summaries is paired with an `eval=FALSE` chunk showing the real call on a full fit (for example lines 94, 142, 185, 211, 395, 659, 708, 766, 868). Examples are complete in the exemplar's sense.

## 2. Section-by-section findings

### What this lesson answers (lines 27-67)

Describes. It tells the reader what was fitted and how the files are organised, but not why a known-truth simulation is the right way to learn the outputs, and it does not show the objects the reader will handle.

- Line 29 (concept not named): the purpose is stated as three abstract contrasts; add one sentence saying that with real data you never know the truth, so a simulation is the only place to see what each output can and cannot tell you.
- Line 31 (next question unanswered): "It replaces the output tour in the original occJSDM vignette" is history the reader cannot use; replace with a sentence saying how this community differs from `sampledata` in the quickstart (two samples, two primers, six PCRs here, against three samples, three primers, two PCRs there).
- Line 33 (choice without reasoning): the perfect-observation fit is introduced without saying why it exists; add that it shows what the occupancy part can recover once detection error is removed, so the gap between the two fits is the cost of imperfect detection.
- Line 33 (wrong audience): "verified model source recorded in Lesson 1; that source matches the code on main" and "This is not a before/after comparison of software versions" answer a maintainer's question; suggest cutting both or moving them to the reproduction section.
- Line 35 (concept not named): `eval=FALSE`, "knit" and "compact saved summaries" are used without explanation; say in one sentence that some chunks are shown but not run, because they need a full fit that takes too long to build.
- Lines 39-64 (data not shown): three teaching objects (`lesson`, `outputs`, `diagnostic_examples`) and `known_truth` are loaded but never shown, unlike the quickstart's `str()` and `head()`; add a short bullet list of what each holds and that they are teaching summaries, not package output.
- Line 62-63 (choice without reasoning): `ggtern::theme_bw()` with a comment about ternary-plot dependencies will puzzle the reader; either explain in one sentence or use `ggplot2::theme_bw()` if that works (question 7).

### Do environmental effects really go undetected? (lines 68-135)

Mixed. The coefficient-reading guidance at 132-134 teaches well; the opening answers a question this reader has not asked.

- Line 68 (concept not named): the heading presupposes an earlier claim that effects go undetected; reframe it as the reader's question, such as how well the model recovers each species' environmental response.
- Line 70 (concept not named): "after accounting for the other model components" leaves the reader to guess; name them (the other covariate and the hidden site factors, with a forward pointer).
- Line 72 (wrong audience): the paragraph about the older `sampleresults` object, with four inline counts, is a response to an issue; suggest cutting it or reducing it to one sentence linking back to the quickstart's `sampleresults` and why its truth is unknown.
- Lines 74-86 (data not shown): `outputs$coefficients` is used without showing its columns; add `head()` or name the columns (`arm`, `block`, `term`, `estimate`, `lower`, `upper`, `truth`, `excludes_zero`).
- Line 85 (concept not named): "95% intervals" and "excluding zero" are the lesson's main evidence criterion but are never defined; add one sentence on what a 95% credible interval is and why excluding zero is read as a resolved direction.
- Line 88 (next question unanswered): "It does not reveal the underlying occurrence probabilities" needs its reason: each site gives one presence or absence, so even perfect observation leaves the probability uncertain.
- Line 92 (concept not named): "The public function" is developer framing; say "the package function".
- Line 106 (wrong audience): "The exporter joins ..." and "The environmental predictor scale was checked by reconstructing the simulation's full ecological predictor from the fitted design matrix" serve package verification; replace with one reader-level sentence that truth was put on the same standardized scale as the estimates.
- Line 132 (next question unanswered): the per-SD scale implies `runOccJSDM()` standardizes covariates; say so, and that the reader therefore does not need to standardize their own.

### What does an effect mean for a species' distribution? (lines 136-204)

Mostly teaches. The gaps are jargon used ahead of its explanation and one piece of practical advice stranded in the appendix.

- Line 138 (concept not named): "hidden site-factor contribution" is first explained at line 432, 300 lines later; give a one-clause plain gloss here (unmeasured site conditions shared by several species) with a forward link.
- Line 140 (next question unanswered): it lists three things a profile is not (maps, fitted probabilities, new-site averages) before any of them has been introduced; either add forward links to the sections that cover each or trim to one sentence.
- Lines 150-177 (data not shown): `outputs$gradients` columns (`x`, `med`, `low`, `high`, `truth`) are not shown or named.
- Line 166 (practical advice missing): the x axis is in standard deviations; the conversion to raw units (`raw_value = mean + sd * standardized_value`) is at line 1855 in the appendix and should be stated here.
- Line 181 (concept not named): "posterior-draw-by-species matrix" and "pooling retained iterations across chains" are used before draws, iterations or chains are explained (line 1233); add one sentence on what a posterior draw is at its first use, or in the opening.
- Line 181 (next question unanswered): "24,000 rows" appears without its arithmetic; say four chains of 6,000 kept iterations (stated only at 464 and 1036).
- Line 195 (practical advice missing): having said the two ways of averaging differ, say which to report (the mean of the transformed draws, as `colMeans()` gives).

### Traits ask a harder, different question (lines 205-390)

Mixed. The opening (207-209) and the array recipe (309-332) teach; the native-plot subsection is repetitive; the cancellation subsection (334-390) is a research diagnosis that the reader cannot repeat on real data.

- Line 207 (concept not named): "the trait acts on the environmental response" is abstract; add an ecological example (for instance, a species trait such as body size that makes some species respond more strongly to a gradient).
- Line 218 (next question unanswered): it explains that fitting standardizes traits but not what this means for the reader's own reporting; add that trait coefficients are per standard deviation of the trait across the species surveyed.
- Line 242 (practical advice missing): the result (one of four intervals excludes zero) is stated but not turned into advice; add what this implies for a survey of ten or so species (expect weak trait results; more species, not more PCRs or sites, is what helps).
- Line 246 (concept not named): "native" first appears here and is used about 40 times; say once that it means the package's own plotting function, or replace it with "the package's".
- Line 248 (wrong audience): "Knitting displays the exported figures from that unchanged fit, including all 24,000 retained draws" reassures a reviewer, not a learner; cut.
- Lines 250-305 (incomplete or redundant example): the two native chunks differ only in `covName` and title; show one and say that the second is the same call with the other name.
- Line 255 (choice without reasoning): the custom theme with `axis.text` angle reset is unexplained.
- Line 336 (concept not named): "unmeasured species traits" (latent traits, set by `n_lattrait`, which appears unexplained at line 1029) are introduced without saying what they represent or how the reader controls them.
- Line 340 (reasoning compressed): "Because these regressions use the same predictors, the contributions add" relies on a property of least squares this reader will not know; add one sentence or cut the claim.
- Lines 343-355 (data not shown): `known_truth$jsdmParams_true$B` is used without saying it is the species-by-covariate matrix of true coefficients.
- Line 389 (scope): the subsection ends by saying it proves nothing general; together with line 387 this is a strong case for moving 334-390 to an appendix (question 2).

### Residual species associations: did we recover what was put in? (lines 391-429)

Describes, and thin: three prose paragraphs. The most useful reading guidance for this output is in the appendix at line 1759.

- Line 393 (concept not named): "residual correlation" is defined in model terms ("the model's shared hidden site component"); say in the reader's terms that it measures whether two species occur together, or apart, more often than their measured environmental responses predict.
- Line 393 (next question unanswered): the reader will ask how the number of factors chosen in `n_factors` limits these correlations; add one sentence.
- Lines 395-400 (output not shown): only the first array dimension is described; say the array is quantile by species by species.
- Line 428 (practical advice missing): "Compare the pattern and magnitude" works only with truth; say what to do with real data, which is the uncertainty display currently at line 1759.
- Line 1759 (misplaced): the finding that all 45 pairs are unresolved, and how to read the X markers, belong here, not in the appendix.

### Ordination: compare the combined effect before naming the axes (lines 430-652)

Mostly describes, and pitched above this reader. It is the longest interpretive section (223 lines) and most of its length is truth-assisted axis alignment that the reader cannot do with real data (line 462 says so).

- Line 430 (concept not named): "ordination" is never explained as a term; add a sentence saying it places sites and species on a few axes summarising the shared unexplained variation.
- Line 432 (concept not named, length): a 987-character paragraph that introduces the lesson's only formula, the notation `logit(psi_ij)`, scores, loadings, residual co-occurrence, two functions and rotation invariance; split into four paragraphs and say in words what each symbol is before the formula.
- Line 432 (next question unanswered): rotation invariance is stated but its consequence for the reader is not; add that axis signs and order can flip between fits, so do not interpret an axis by its number or sign alone.
- Line 456 (next question unanswered): "the observations at a site helped estimate its hidden scores" needs its "so what": this is why fitted sites will look better recovered than new sites.
- Lines 460-464 (wrong audience): Procrustes alignment, matrix algebra and "The exporter also verifies this for every draw" are simulation-verification material.
- Lines 470-478 (ordering): the ordinary calls, which are what the reader will run, come after the alignment explanation; move them to the start of the section.
- Lines 480-534 (wrong audience): a 55-line alignment loop that edits internal slots (`results_output$jsdm_output$U_output`); suggest an appendix or the dev scripts, leaving one sentence on what alignment does.
- Line 525 (choice without reasoning): every tenth site is shown; the comment gives the reason but the prose does not.
- Line 581 (next question unanswered): the key result, that site medians cluster near zero despite varied true scores, is not explained; add a sentence on why (each site gives little information on its hidden scores, so estimates shrink toward zero) and what it means for reading a real ordination.
- Line 581 (wrong level): the circle radius formula is package implementation detail; one sentence ("rough size guides, not credible regions") is enough.
- Line 651 (length): 685 characters mixing reading guidance with four caveats; split.
- Section-wide (practical advice missing): the reader is never told what an ordination of their own survey is useful for, for example mapping site scores to look for an unmeasured gradient.

### Variation partitioning: an allocation within the model (lines 653-686)

Mixed. The cautions are good; the concept itself is not introduced.

- Line 655 (concept not named): say what question variation partitioning answers for an ecologist (how much of each species' distribution is attributed to measured environment versus shared unexplained variation).
- Line 655 (next question unanswered): "In this non-spatial example the spatial fraction is zero"; say in one clause what changes for a spatial fit.
- Line 660 (output not shown, incomplete example): the returned table's structure is not shown, and `plotVariancePartitioning()`, listed in the function finder at line 2067, is never demonstrated; add the plot call.
- Line 681 (concept not named): the definition ("changes in the standard deviation of probabilities when model components are combined, with negative increments truncated and contributions normalized") needs a plain-language version first.
- Line 685 (ambiguity): "the preceding figure" means the ordination contribution plot at line 443, two sections back; name it.

### Collection effects and detection probabilities (lines 687-759)

Mostly teaches. The detection-effort subsection is the best example of the target style in the lesson.

- Line 689 (concept not named): give a concrete collection covariate (water volume filtered, turbidity) so the reader can map it to their own fieldwork.
- Line 699 (output not shown): the plot shows only the slope, though line 704 discusses the intercept; say so or show both.
- Line 721 (length, wrong audience): 752 characters joining two helpers, a Lesson 1 pointer and a threshold aside ("Its truth accounts for whether simulated reads actually pass the fitted threshold") that means nothing to the reader here; split and cut the aside.
- Line 739 (practical advice missing): it distinguishes the expectation interval from the survey-outcome interval but not which to use when planning a survey; add one sentence.
- Line 758 (practical advice missing): state the take-home with numbers from the figure (for this community, a second field sample buys more than more PCRs).

### Distinguish fitted probabilities, occupancy states and new-site predictions (lines 760-785)

Teaches, as a bridge section. One incomplete example and one factual inconsistency.

- Line 762 (next question unanswered): it opens by saying what is wrong with the function's name before saying what the function is for; lead with its use (line 780 has it).
- Lines 766-774 (incomplete example): only `computeConditionalOccupancyProbs()` is shown; add the one-line `computePredictiveOccupancyProbs()` call beside it, since line 780 recommends reporting it.
- Line 784 (factual inconsistency): "Spatial prediction remains part of the planned Lesson 2" contradicts line 31 and line 2074, which treat Lesson 2 as existing (see section 5).

### Predict occupancy at genuinely new sites (lines 786-1085)

Mixed. Lines 786-917 teach very well; lines 919-1085 shift to statistical verification of the comparison (Monte Carlo error, paired standard errors, prediction-fit diagnostics, WAIC internals).

- Line 790 (practical advice missing): the one- versus two-factor comparison is set up but the reader is never told how to choose `n_factors` for their own data; add a sentence (the comparison method at 970-1014 is the tool, and the generating number need not win).
- Line 806 (next question unanswered): it says the model standardizes with training constants; make explicit that `predictNewSites()` does this itself when given raw values (implied only at 866).
- Line 829 (data not shown): `lesson$input$jsdm$sigma_h` is used without saying it is the standard deviation of the hidden site factors.
- Lines 871 and 891 (reasoning out of order): `set.seed()` comes before the reason (the function draws random hidden conditions), which appears at 891; move the reason up.
- Line 876 (choice without reasoning): `useSpatial`, `confidence` and `verbose` are not explained.
- Line 921 (next question unanswered, wrong audience): the point estimate used for scoring comes from "the exporter" and "the development helper", not a package function, so the reader cannot reproduce it; say plainly whether the reader should use the `predictNewSites()` median or how to compute the mean.
- Lines 993-1012 (wrong level): paired site standard errors and an MCMC standard error from "the prediction verifier" are careful but heavy; the teaching point (the difference is too small to rank the models) fits in two sentences.
- Line 1012 (wrong audience): "documented in the prediction verifier" points to a dev script.
- Line 1029 (choice without reasoning): `n_lattrait = 1` is set without comment.
- Lines 1030-1031 (incomplete example): the MCMC settings and priors are passed as `prediction_examples$manifests$...`, so the reader cannot see the call's values; write them literally (they are stated in prose at line 1036).
- Lines 1036-1064 (ordering): diagnostic thresholds (`rhat > 1.01`, `ess < 400`) and Rhat, ESS and MCSE are used here before the diagnostics section introduces them at 1233-1300.
- Line 1066 (concept not named): WAIC is never defined, and "The old walkthrough extracted WAIC" is history; add one sentence on what WAIC is meant to measure.
- Line 1083 (wrong level): the limitation is right and the bold advice is excellent practical guidance, but "likelihood terms for the sampled, unobserved site and collection states" needs a plain version (it rewards fitting the training survey, not predicting new sites).

### Put the observations, inferred states and truth in one table (lines 1086-1228)

Teaches. The weight is in the truth-joining code, which the reader will not need.

- Line 1107 (wrong audience): "No estimates are recalculated or replaced with truth" is reviewer reassurance; cut.
- Lines 1140-1156 and 1159 (choice without reasoning, wrong audience): the `knitr::pandoc_to()` branch and the sentence about the Markdown rendering exist for the site build; hide the fallback and say only that `plotLatentPresences()` returns a `gt` table.
- Line 1161 (next question unanswered): "sample 3" at site 2 and "sample 6" at site 3 will confuse a reader expecting samples 1 and 2 within each site; say sample identifiers run across the whole survey.
- Lines 1167-1201 (wrong audience): 33 lines join simulator coefficients (`beta_theta_true`) to build true probabilities; consider `echo=FALSE` for the join, keeping the explanation at 1167 and the table.

### Check computation as well as ecological recovery (lines 1229-1595)

Mixed. Lines 1231-1442 teach well. The mirror subsection (1444-1515) is a research report on a different community. The closing subsections (1516-1595) repeat earlier caveats.

- Line 1229 (structure): the quickstart, which Doug approved, checks chains before reading any estimate (its lines 75-86), and line 1231 itself says to check diagnostics before interpreting uncertainty, yet this section comes after every interpretation; see section 4.
- Lines 1251-1258 (format): a pipe table, which the repo rules forbid in Visual-mode files; convert to bullets.
- Lines 1273-1274 (reasoning out of order): the `1.01` and `400` thresholds are applied before their justification at line 1300.
- Line 1287 (concept not named): "weaker low-contamination assumptions" should name the prior values, or say in words what was changed.
- Line 1300 (next question unanswered): it explains coda versus `posterior` diagnostics but not which the reader should use; add one practical sentence.
- Line 1304 (next question unanswered): extracting traces needs internal slots (`fitmodel$results_output$beta_theta_output`, `fitmodel$X_theta`); say whether there is an accessor, or that this is the supported way (question 6).
- Line 1324 (wrong audience): "These are excerpts of the full fit, not newly fitted models" is reassurance; cut.
- Lines 1402 and 1240 (incomplete example): `/path/to/full-fits/...` paths exist only for someone reproducing the archive; say what the reader substitutes (their own fit).
- Lines 1446 and 1448 (wrong audience, length): two paragraphs of 1,397 and 1,340 characters describing a separate 300-site study, its chain counts, seeds, burn-in and thinning; the reader needs the symptom, the check and the action (1485-1514), not the study design.
- Line 1485 (concept not named): `B0` appears without saying it is the occupancy intercept.
- Line 1487 (wrong level): array-dimension caveats for single-species and single-chain fits are correct but dense; move into code comments.
- Lines 1520-1527 (format, duplication): a second pipe table, partly overlapping the first; convert to bullets and consider moving to a reference appendix.
- Line 1529 (wrong audience): "also used by our offline verification scripts" is an aside for maintainers; cut.
- Line 1548 (ordering): the overview table ("Summaries and what to do next") comes after the detailed walk-through; an overview first would read better.
- Line 1569 (repetition): repeats the caveats of lines 1260 and 1285.
- Line 1571 (wrong audience): "the exporter also calculates diagnostics directly"; the reader-level point is only that trait coefficients are not in the public table.
- Line 1594 (repetition): repeats the WAIC limitation from line 1083 almost in full.

### Appendix: use the package's plotting functions (lines 1596-2019)

Describes. It is a gallery of the package's own plots with truth overlays, valuable as reference, but most of it repeats concepts already taught with custom plots, and several of its pieces of guidance belong in the main sections.

- Line 1598 (wrong audience): "archived PCR fit", "provenance" and "not a replacement fit" are reviewer language.
- Lines 1604-1613 (choice without reasoning): another `ggtern` theme, as at line 63.
- Line 1669 (repetition): repeats line 134.
- Line 1673 (repetition, wrong level): repeats line 183; the collection truth formula is verification detail.
- Line 1717 (repetition): the read-threshold adjustment is explained here and again at 1395, 1791 and 1954.
- Lines 1723 and 1742, 1796 and 1826, 1871 and 1898 (redundant examples): three pairs of near-identical chunks; show one of each and describe the other.
- Line 1759 (misplaced, length): 937 characters of essential reading guidance for the residual correlation heat map; move to the associations section and split.
- Lines 1785-1849 (duplication): the cumulative-detection plot repeats the concept of 723-758; merge with that section, where the expectation-versus-outcome contrast can be taught side by side.
- Lines 1851-1921 (duplication): repeats 136-177 with the package function; keep only the raw-scale conversion (1855) and move it to 166.
- Lines 1923-1925 (wrong audience): PR #13 bug history, including impossible probabilities of 1.772 to 1.815, and a column rename are release-note material; move to NEWS and keep, at most, one sentence that `plotCovariateEffect()` works in original units.
- Line 1952 (good practical tip) and line 1929 (good definition) should survive any cut.

### Reproduce the extraction or find a function (lines 2020-2076)

Two audiences in one section: reproduction commands for maintainers and a function finder for readers.

- Line 2022 (next question unanswered): `dev/simstudy/vignette-lesson/README.md` is excluded from the built package by `.Rbuildignore`, so the reader needs a GitHub link, not a relative path.
- Lines 2024-2046 (wrong audience): sixteen `Rscript` export and verify commands; the dev README already holds these, so a link would do.
- Line 2057 (wrong audience): "Student-facing figures use tidy tables ... the export and verification scripts document and check the array calculations" is a note to reviewers.
- Lines 2059-2074 (format, placement): the function finder is the most useful reference for the reader but is a pipe table, which the repo rules forbid; convert to a bullet list and consider moving it near the top as a map of the lesson.
- Line 2062 (wrong audience): "corrected conditional response curves are also available" refers to the bug history.

### References and further reading (lines 2077-2085)

- Lines 2079, 2083 and 2085: Cai et al. (2025), Leibold et al. (2021) and Pichler et al. (2025) are listed but cited nowhere in the lesson; only Ji et al. (2025) is cited (line 780). TODO.md line 146 plans a future exercise that would cite the first two; until then, remove or cite them (question 8).

## 3. Jargon list

Each term is listed with the line of first use and, where it is explained, the line of explanation.

- `eval=FALSE`, knit, knitting: line 35; never explained.
- perfect observation: line 33; explained there, but not why it is fitted.
- hidden site factors: line 33; explained at line 432.
- posterior, posterior uncertainty: line 66; never defined for this reader.
- public function: line 92; developer term, never explained.
- posterior draws, array: line 92; never defined.
- the exporter: line 106; never explained (the dev scripts at 2022 imply it).
- design matrix: line 106; never explained.
- 95% interval, credible interval: line 85 (prose "credible" first at line 581); never defined.
- log-odds: line 126 (axis label), line 132 (prose); never defined (the logit formula appears at 432).
- one-standard-deviation change, standardized predictor: lines 132 and 181; partly explained at 1855.
- response profile: line 140; explained there.
- hidden site-factor contribution: line 138; explained at 432.
- chains, retained iterations, pooling across chains: line 181; chains explained at 1233; burn-in first at 1036, never defined.
- probability scale, inverse-logit, `plogis()`: line 181 ("probability scale"); `plogis()` explained at 819; "inverse-logit" first at 862.
- trait standardization: line 218; explained there.
- native (function, plot): line 244; never defined.
- `purl`: chunk options from line 250; never explained (code, but visible).
- unmeasured species traits (latent traits, `n_lattrait`): line 336; not explained; `n_lattrait` unexplained at 1029.
- oracle diagnostic: line 340; explained there.
- residual correlation: line 393; defined in model terms only.
- loading: line 426; defined at 585.
- ordination: line 430 (heading); never defined.
- logit and the notation `psi_ij`, `beta_0j`, `U_i`, `L_j`: line 432; symbols explained there, logit not.
- rotation, reflection, orthogonal Procrustes: lines 432 and 462; explained briefly at 462.
- biplot: line 438; explained at 622.
- marginal 95% intervals, joint credible regions: lines 571 and 581; partly explained.
- variation partitioning: line 653; never defined in the reader's terms.
- conditional on presence: line 689; explained in context.
- latent 0/1 occupancy state: line 764; explained there.
- MCMC: line 790; explained at 1233.
- marginal probability, conditional probability: lines 571 (caption) and 817; explained at 817.
- Monte Carlo error, MCSE: line 921; partly explained at 1012 and 1064.
- RMSE: line 960; never defined (signed and absolute error are).
- Brier score, negative log score: line 974; explained there.
- Rhat, effective sample size (ESS): line 1036; Rhat explained at 1588 (and 1448), ESS at 1588; the console thresholds at 1237.
- WAIC: line 1066; never defined.
- `theta0`, `p`, `q`: `theta0` first in prose at line 1258 (table); `p` at 727; `q` at 1257; all explained in the table.
- bulk and tail ESS, rank-normalized split-chain Rhat: line 1300; pointed to the Stan guide.
- `B0`: line 1485; never explained.
- mirror (labelling): line 1446; explained there.

## 4. Length and structure

Verdict on Doug's concern: yes, the length is partly a structure problem, and a teaching pass on the current structure would make it worse. Section 2 lists roughly 100 gaps. Filling each with a sentence or two would add 150-200 sentences to a lesson that already has about 10,000 words of prose. The same lesson has at least 350 lines that serve package verification, maintainers or a separate study rather than this reader. Cutting or moving those first makes room for the teaching additions while still leaving a shorter lesson.

Three audiences are mixed in the lesson:

- The ecologist learning to read outputs (the intended reader): most of 27-450, 653-917, 1086-1165, 1229-1442 and 1590-1592.
- A reviewer verifying that the teaching figures are faithful to the fits: lines 33, 106, 248, 464, 921, 1012, 1107, 1324, 1529, 1571, 1598, 2022-2057, and the truth-joining code at 1167-1201.
- A methods reader interested in identifiability: Procrustes alignment (460-534), the trait cancellation (334-390), the paired-score and MCSE analysis (993-1012) and the mirror study (1444-1514).

The lesson repeats itself in these places:

- The read-threshold adjustment to true detection rates: lines 721, 1167, 1395, 1717, 1791 and 1954.
- The baseline is not average occupancy: lines 183 and 1673.
- An interval crossing zero is not an absent effect: lines 134, 719, 1669 and 1759.
- The WAIC limitation: lines 1083 and 1594.
- Diagnostics do not certify accuracy: lines 1260, 1285, 1569 and 1588.
- Response curves with site factors at zero: lines 138-177 and 1851-1921.
- Expected detections against survey outcomes: lines 723-758 and 1785-1849.
- Undefined correlations for OTU_4: lines 426 and 1759.
- "24,000 retained draws": lines 181, 248, 464, 704 and 921.

The figures are built twice. The main sections use custom `ggplot2` code on teaching summaries (`outputs$...`) that the reader will never have for their own fit. The package's own plot calls, which the reader will actually use, sit in `eval=FALSE` chunks and in a 424-line appendix. That is why so much of the appendix repeats the main text. One option is to make the package plot plus a truth overlay the main figure in each section and drop the custom duplicate (question 1). That alone would remove much of the appendix.

Structural changes needed for the prose pass to work:

- Move the diagnostics section (1229-1595) to straight after the opening, matching the quickstart's order and the section's own instruction at line 1231.
- Move the prediction-fit diagnostics (1036-1064) into the diagnostics section, so that thresholds are introduced before they are used.
- Fold the appendix's unique teaching into the main sections: the heat-map reading at 1759 into 391-429, the raw-scale conversion at 1855 into 166, the cumulative-detection plot at 1785-1849 into 723-758 and the rate plots at 1927-2018 into 687-759. Then either delete the remaining duplicate native-plot gallery or keep it as a short reference appendix.
- Move the Procrustes alignment code (480-534), the alignment-based native ordination plots (543-649) and the trait cancellation (334-390) to a "simulation checks" appendix, or to a separate article (question 2).
- Compress the mirror subsection (1444-1514) to the symptom, the `chain_summary()` check, the action and a link to the study report. Alternatively, move it whole to a separate article, with a pointer from the diagnostics section (question 3).
- Split the reproduction section. Replace the commands with a link to the dev README on GitHub, and turn the function finder into a bullet list placed near the top or as the final section.
- Move the PR #13 history (1923-1925) to NEWS.
- Consider splitting the lesson in two (question 4). The new-site prediction section (786-1085) is self-contained, with its own data file and its own extra fit, and amounts to 300 lines. It could be its own lesson: "Predict new sites and compare models". The planned post-beta rename (TODO.md line 116) renumbers every lesson anyway, so that would be the moment to do it.

Line budget if the cuts above are made (current line count first, then the target):

- Setup and opening (1-67): 67 to 60, with the teaching objects shown.
- Diagnostics, moved forward (1229-1595): 367 to 160, with the mirror case at about 25 and the parameter-array reference moved to the appendix.
- Environmental effects (68-135): 68 to 60.
- Response profiles and baseline (136-204): 69 to 65.
- Traits (205-390): 186 to 100, with the cancellation moved and one native chunk.
- Residual associations (391-429): 39 to 60, gaining the heat-map guidance.
- Ordination (430-652): 223 to 100, with alignment moved.
- Variation partitioning (653-686): 34 to 40, gaining the plot call.
- Collection and detection effort (687-759, plus appendix 1785-1849 and 1927-2018): 73 to 120 after merging.
- Fitted against conditional against predicted (760-785): 26 to 30.
- New-site prediction (786-1085): 300 to 180, whether it stays or becomes its own lesson.
- Latent-presence table (1086-1228): 143 to 100.
- Reference appendix (function finder, parameter arrays, optional simulation checks): about 150.
- References: 9 to 6.

That totals about 1,330 lines including the teaching additions, against 2,085 now. If the prediction section becomes its own lesson, this one is about 1,150 lines.

## 5. Copyedit observations

Long paragraphs that should be split (character counts):

- Line 1446 (1,397): the mirror study, with one sentence carrying a colon and a five-part list.
- Line 1448 (1,340): the chain-count sentence beginning "Counting the two chains of an earlier, shorter fit" is very hard to follow.
- Line 1514 (1,115).
- Line 432 (987).
- Line 1759 (937).
- Line 780 (831).
- Line 721 (752).
- Line 72 (730).
- Line 1300 (692).
- Line 651 (685).
- Line 1487 (675).
- Line 1083 (656).

Long sentences:

- Line 432, "That product is what carries residual co-occurrence: ...".
- Line 721, "`computeAverageCollectionProbs()` instead returns ...".
- Line 780, "The two are most useful together: ...".
- Line 1446, "For one species in that 300-site community ...".
- Line 1514, "The study ran 8 chains ...".

Asides that repeat their sentence, the pattern Doug removed in the quickstart:

- Line 966: "that last sentence is only an illustration of the unit, not a claim that all errors equal ten points".
- Line 1107: "No estimates are recalculated or replaced with truth".
- Line 1324: "These are excerpts of the full fit, not newly fitted models".
- Line 1440: "not a claim that only that many iterations were run".
- Line 1759: "Unlike the crosses in the coefficient plots, these Xs are uncertainty markers; the numbers supply the truth", which follows two sentences that already say so; the same paragraph's "This is not a true correlation of zero" repeats line 426.

Slashes and "and/or":

- No "and/or" in the lesson.
- Slash pairs to replace with words: line 33 "before/after", line 1107 "sample/primer", line 1789 "collection/PCR", line 2022 "simulation/fitting", line 2071 "site/sample", and in the tables "Collection intercept/effect" (1525) and "PCR detection/false-positive rate" (1526).
- "presence/absence" (six times) matches the quickstart and can stay.

Missing full stops: none found in prose.

Formatting against the repo's Visual-mode rules:

- Pipe tables at lines 1251, 1520 and 2059.
- Curly quotes at line 862, where the rest of the file uses straight quotes.
- Double blank lines at 1084-1085 and 2075-2076.
- The `editor_options: markdown: wrap: none` block at lines 13-15. The quickstart has the same block, so this is probably intended, but the CLAUDE.md wording ("Never add an `editor_options: markdown: wrap:` block") reads as forbidding it.
- No unescaped arrows or tildes in prose.

Hardcoded results in prose that will drift if the bundles are regenerated, while nearby numbers use inline R (lines 72, 338, 387, 1285, 1440, 1442):

- Line 860: 93% and 82%.
- Line 966: 3.8 and 10 percentage points.
- Line 1010: -0.00014 and 0.000067.
- Line 1012: 0.00025.
- Line 1064: 0.39 percentage points.
- Line 1759: all 45 pairs.
- Lines 181, 248, 704 and 921: 24,000.

Factual claims I doubt or could not verify:

- Line 784 against lines 31 and 2074. Line 784 calls Lesson 2 "planned", while line 31 treats it as an existing spatial lesson and line 2074 describes it as "a worked site-arrangement sweep". These cannot all be current. TODO.md line 116 also plans renumbering after the beta, which will change all three.
- Line 780, Ji et al. (2025). The claims that the paper reported predictive rather than conditional probabilities, and that its Supplementary Information 12 works through such cases, cannot be checked from the repo. Doug is a coauthor and should confirm the SI number.
- Line 806, "Normal distribution with mean zero and standard deviation 10" for the environmental covariates. Not verified. This is worth checking against the simulation settings, because an SD of 10 matters for the claim at 808 about values outside the training range.
- Line 1300, "classical `coda` Rhat and effective sample size". ESS is confirmed (`coda::effectiveSize`, R/diagnostics.R line 55); I did not trace the Rhat call.
- Lines 1446, 1480 and 1512, checked and correct. The default priors on `q` and `theta0` are both Beta(1, 20) (R/runOccJSDM.R lines 852-857). The prior mean is 1/21, about 0.048, and the share of Beta(1, 20) above 0.25 is 0.75 to the 20th power, about 0.3%.
- Line 1237, checked and correct. The console warnings fire at Rhat above 1.1 and ESS below 50 (R/diagnostics.R lines 478 and 487).
- Lines 181, 464, 1036 and 1440 are internally consistent: 4 chains of 6,000 make 24,000, and 4 chains of 12,000 make 48,000. Line 1227 (eight rows) and line 1161 (twelve PCRs, global sample numbers) are also consistent with the design at line 31.

## 6. Estimate

- Paragraphs needing a teaching addition or change: about 65-75 of the roughly 176 prose paragraphs, concentrated in the opening, ordination, variation partitioning, the second half of new-site prediction and the diagnostics section.
- Paragraphs or passages to cut or move for audience or repetition: about 35-45, plus the code blocks named in section 4.
- Sentences needing copyedit: about 60-80, mostly splitting the twelve long paragraphs listed in section 5, removing reviewer asides and replacing hardcoded numbers.
- Pass size: large. It is large because of the structure, not because the prose is poor; where the lesson teaches, it already teaches in the target style. A teaching-only pass without the structural changes would be medium in effort, but would leave a longer and harder lesson.

## 7. Questions for Doug

1. Main figures: should each section show the package's own plot with a truth overlay, so the reader sees the call they will use, and drop the custom `ggplot2` duplicate built from teaching summaries? This decides whether most of the 424-line appendix survives.
2. Simulation-only material (trait cancellation 334-390, Procrustes alignment 460-649): keep it in the lesson, move it to an appendix, or move it to a separate methods article? The reader cannot do any of it with real data.
3. The mirror case (1444-1514): keep a short warning plus the per-chain check here and move the study account to its own article, or keep it whole? It describes a different 300-site community.
4. Should new-site prediction and model comparison (786-1085) become its own lesson? If so, this could be timed with the post-beta renumbering in TODO.md line 116.
5. Should diagnostics move to the front, as in the quickstart? Line 1231 already argues for that order.
6. Is reading traces from internal slots (`fitmodel$results_output$...`, lines 1304-1415 and 1487-1510) the supported route for users, or should the lesson wait for an accessor (TODO.md line 109 mentions an unwritten per-chain helper)?
7. Can the `ggtern::theme_bw()` workaround (lines 63, 255, 1612, 1868) be replaced by `ggplot2::theme_bw()`, or does loading occJSDM really break the plain theme? If it is needed, the reader needs one sentence explaining why.
8. Should the three uncited references (Cai et al., Leibold et al., Pichler et al.) stay ahead of the planned internal-structure exercise, or come out until it exists?
9. Which statement about Lesson 2 is current (lines 31, 784 and 2074)? Should the lesson's cross-references be written now for the post-rename numbering, or left until the rename?
10. Is the opening paragraph about the old `sampleresults` and "undetected" effects (lines 68-72) answering a specific user report that should be cited, or can it go?

# Style diagnostic: Lesson 4 (four-package JSDM comparison)

Read-only diagnostic of `vignettes/occJSDM-lesson-4.Rmd` (1,532 lines) against the 2 October 2026 teaching criteria in `criteria.md`, with the quickstart (`vignettes/occJSDM.Rmd`) as the exemplar. The lesson was read in full, in order, in six passes (lines 1-250, 250-520, 520-800, 800-1080, 1080-1320, 1320-1532). Rendered outputs were checked against the committed `occJSDM-lesson-4.md` and two figure PNGs in `teaching-data/`; nothing was run or rendered. Line numbers refer to the Rmd. Of the 1,532 lines, about 899 are code inside chunks and about 177 are prose paragraphs or bullets, so the line count overstates the reading length; from line 1132 onward every chunk is hidden (`echo = FALSE`).

## 1. What the lesson does well against the criteria

- **Names concepts in the reader's terms.** The hidden factor as "an unmeasured soil property" (95); "perfect observation does not mean perfect knowledge of probability" with the 20% example (108); the recommender analogy for latent factors (221); the hypothetical 40% versus 70% site (298); warm sites that are nearly always dry for collinearity (943); drought-tolerant species for a trait effect (1042); the "hills" image for local optima (397).
- **Answers the next question.** "Why are the new-site errors lower here, even though sampled sites supply more observations?" is asked and answered (540, also 496); "failed does not necessarily mean the software crashed" (907); the absence of sjSDM from the trait part is not a failed result (1053); missing table entries are not zero coverage (1245).
- **Shows the data.** Both input tables are printed (99-106); the probability-versus-observation tile plot puts truth and data in one figure (110-134); the full generating parameter table is shown (140-143); a downloadable species error table is linked (601).
- **Spells out reasoning.** Why probabilities are compared instead of raw coefficients (215); why test sites reuse the training standardisation (202); why the covariance is monitored rather than factor axes (389); why the gradient is set to +1 SD in the exercise (235); why the curved scenario is fitted twice, separating too little information from an unsuitable shape (931); why the dashed and dotted reference lines differ for a nonzero trait effect (1083).
- **Practical advice and honest scoring.** Selection by training criterion, never by truth (414, 452); outcome scores for real data where truth is unknown (672); the replication unit is the community, not the species-by-site cell (1038, 1202); coverage must be read with width (1221, 1325); flagged fits stay in (909).
- **Structure that teaches.** Learning objectives up front (42-50); section 14 leads with its conclusions (1169-1177) then evidence; exercises in section 11 (841-847) use only saved tables, and exercise 5 extends the conditional-prediction example.

## 2. Section-by-section findings

### Opening: What are we comparing? (34-89)

Teaches the framing well but never tells an eDNA user why a comparison of four JSDM packages matters to them.

- 38, reader's next question: the reader will ask why an eDNA occupancy package is being compared with packages that have no detection model; say that only occJSDM handles detection error, and this lesson checks that its ecological core, the part shared with established JSDMs, is as good as theirs.
- 38, concept not in reader's terms: "the pure JSDM portion of occJSDM" should be tied to the quickstart's presence/absence mode (quickstart line 45), which is exactly what is fitted here.
- 38, concept not named: gllvm, sjSDM and Hmsc are never introduced; one clause each on what they are and who uses them would orient a reader who knows none of them.
- 40, choice without reasoning / audience: the sjSDM revision history is for returning readers; for a new reader it should become a one-sentence promise that section 5 explains why one fit needed several starts.
- 54-89, part of the data not shown: `comparison` is unpacked into `training`, `truth` and `predictions` without its structure being shown; a `str(comparison, max.level = 1)` beside the code, as in the quickstart, would show the reader every part they will handle.
- 88, missing practical advice: `ggtern::theme_bw()` silently requires ggtern, which the reader is never told to install, and the comment explains a cross-vignette build issue; say what is needed, or use the ggplot2 theme.

### 1. What the models receive, and what we keep secret (91-202)

Teaches.

- 104-106, part of the data not shown: `known_truth` is used (114, 141, 260) without its components being listed; one sentence naming `conditional_probability`, `parameters`, `scaled_coefficients` and `loadings` would prepare the reader.
- 138, concept not explained: "logit scale" is named and the reader is told a slope of 1 is not one percentage point, but not what it is; one sentence on log-odds, and that a slope of 1 multiplies the odds by about 2.7, would complete it.
- 202, missing practical advice: the reader will want to know whether `runOccJSDM()` standardises covariates itself or expects them pre-standardised, and how to apply the training transformation to their own prediction sites.

### 2. How similar are the four models? (204-217)

Describes: a specification table for readers who already know the four packages.

- 206, choice without reasoning: the factor count is fixed at the true value and "this lesson does not choose the count using WAIC", but the reader is not told why (it removes one source of difference between packages) or what to do in a real survey, where the count is unknown; cross-reference Lesson 3's WAIC section.
- 208-213, concepts not in reader's terms: "variational approximation", "PyTorch CPU optimisation", "weight decay", "low-step continuation", "covariance factor dimension", "unconstrained factors", "independent site level" and "latent traits" all appear in one table; a sentence before it should split the packages into two families (Bayesian sampling gives a posterior with uncertainty; optimisation gives a single best estimate) and leave the rest to section 10.
- 215, reader's next question: having learned that Hmsc uses probit while the truth is logit, the reader will ask whether that handicaps Hmsc; say briefly that probabilities are comparable and coefficients are not, pointing to 1337.
- 217, audience: "Doug's fork release v0.2.1" and Mojo mean nothing to this reader; they will instead ask whether the CRAN sjSDM gives the same results; move the fork detail to section 10 or the reproduction record and answer that question.

### 3. What a joint model adds to factoring a matrix (219-285)

Teaches; the strongest conceptual section in the lesson.

- 221, concept overstated: "the dimensions turn out to be genres that nobody labelled" invites the reader to name factors, which line 95 warns against; say the dimensions sometimes look like genres but need not correspond to anything nameable.
- 226, concept not in reader's terms: "Bernoulli outcome", "linear predictor" and "likelihood" arrive in one sentence; give each a gloss (a yes/no outcome, the species' score at the site, how probable the observed data are under given parameter values).
- 235, choice without reasoning: the demonstration uses gllvm's parameters, not occJSDM's, in an occJSDM lesson; the reader will ask why, and whether occJSDM's output can do this.
- 235, factual and structural: "the same gllvm parameters as the marginal example above" refers to the extraction at 332-364, which is in section 4, below; either move the exercise after section 4 or reword.
- 258, part of the data not shown: `point_parameters$gllvm` is used without showing that `beta` has rows for the intercept and two slopes and `loading` has one column per species; that explanation arrives only at 366.
- 283, output not read: the rendered table shows fitted gllvm putting species_10 alone at 89.2% against a true 83.9% at the mean environment (and 72.4% against 64.2% at +1 SD), and the prose comments only on the conditional gap; the reader will ask why the baseline is also off.
- 285, missing practical advice: tell the reader whether occJSDM offers conditional prediction for their own surveys, and if not, that this calculation is the way to get it.

### 4. Two probability questions that must not be mixed (287-366)

Teaches, though 326-328 shifts to a methods audience.

- 296, reader already met this: Lesson 3 (its line 817) introduced marginal and conditional probability; a cross-reference, and the word "conditional" for the sampled-site question, would join the two lessons.
- 302-324, output not read: the illustration prints 11.9% (zero hidden contribution) against 22.5% (averaged), and the prose never says that averaging nearly doubles the probability, or why (at low scores, favourable hidden conditions raise probability more than unfavourable ones lower it); also say that the hidden SD of 2 was chosen large to make the effect visible.
- 326, concept not in reader's terms: "level-zero prediction" and "the inspected sjSDM environmental prediction" are package internals; the practical point for a reader using those packages is which call gives a marginal prediction, or that none does.
- 328, audience: three integration methods in one paragraph, plus "a checked study calculation" the reader cannot see; this belongs in the reproduction notes, leaving one sentence that each package's prediction was converted to the same marginal question.
- 332-364, choice without reasoning: the worked marginal extraction again uses gllvm; an occJSDM reader will want the same calculation from an occJSDM fit, or a pointer to Lesson 3's new-site prediction section.
- 366, output not read: the result (14.6% estimated against 16.1% true for species_01) is printed and never commented on.

### 5. Check the computation before interpreting the ecology (368-452)

Mixed: 370-391 teaches; the two-optima subsection (395-452) mostly describes an optimisation investigation for verifiers, framed by one excellent teaching point at 450.

- 370, choice without reasoning: four chains of 2,000 warm-up and 4,000 kept iterations differ from the quickstart's two chains of 5,000 and 5,000; the reader will ask why, and "warm-up" should be tied to the quickstart's "burn-in" and the `nburn` argument.
- 389, choice without reasoning: the Rhat 1.01 and ESS 400 thresholds are not justified; name them as the standard published recommendations and cross-reference Lesson 3's diagnostics section, which also covers bulk and tail ESS.
- 391, concept not explained and missing practical advice: "EVA" and "inconsistent objectives" are undefined, and the reader is not told how to detect the same trap in their own gllvm fit (run several starts and compare their log-likelihoods).
- 393, choice without reasoning: the "declared 0.1 stability check" on the log-likelihood has no rationale.
- 397, concept order: "exact, deterministic optimiser", "stationary solutions", "curvature check" and "flat direction" come before the hills image that makes them intelligible; lead with the hills.
- 399, jargon: "native sjSDM starts" ("native" recurs at 1242, 1333, 1423).
- 414, reasoning worth teaching: the pre-declared selection rule is good practice; say in a sentence why declaring it before reading truth matters (it prevents choosing the answer you like).
- 450, missing practical advice: the transferable lesson for an occJSDM user is that MCMC can also land in different modes, which is why several chains are run and compared; say so.

### 6. Put the estimated probabilities beside truth (454-496)

Teaches.

- 458-494, output not read: the four panels in each figure look nearly identical, which the prose never says; the reader should be told that at this level the packages are hard to tell apart, so the figures mainly test the shared model and the data.
- 496, repetition: the explanation of why the sampled-site plot is more scattered is repeated at 540; keep one.

### 7. How far wrong, and in which direction? (498-614)

Teaches the definitions; describes the band results.

- 538, aside for another audience: "these are measured results, not notional examples, and they are not the older occJSDM-only sample-size experiment" refers to a study the reader has never met; drop it or say what that study was.
- 546, choice without reasoning: the 20% and 80% band edges are not explained.
- 575, output not read: the prose says how to read the band tables but not what they show; the sampled-site rows have every package overestimating low probabilities and underestimating high ones (occJSDM +8.8 and -10.4 points), which is shrinkage toward the middle and should be named and explained; at new sites gllvm and sjSDM err in the opposite direction from occJSDM and Hmsc in the low and high bands, which the reader will ask about.
- 603-613, output not read: species_01's signed errors are negative while the overall averages are positive; one sentence would use the example the reader just ran.

### 8. What environmental response does each model recover? (616-668)

Teaches. The species examples at 666 were checked against the rendered figures and are accurate.

- 620, reader's next question: "no common uncertainty interval was calculated across all four fitting methods" invites "why not?"; a forward pointer to section 15 (1245) answers it.
- 666, biggest unread result in the pilot: in almost every panel the four packages' curves lie on top of each other, so the misses (species_04 and species_09 on gradient 1, species_03 on gradient 2) are shared by all four; say that these errors come from what 100 sites can reveal, not from any one package.

### 9. If we did not know truth, how would we score predictions? (670-689)

Teaches; the one section a reader can apply directly to a real survey, and it is the shortest.

- 686-689, reader's next question: the scores (Brier about 0.182) have no reference point; the reader will ask whether that is good; because the truth is known here, scoring the true probabilities against the same outcomes would show the best achievable score, and a constant-prevalence model would show the worst sensible one.
- 689, missing practical advice: tell the reader to hold out some of their own sites, and how, so that they can compute these scores; cross-reference Lesson 3's prediction section if it covers this.

### 10. Optional: fit the same configurations yourself (691-833)

Describes: argument lists for readers who will rerun the study, plus audit notes.

- 693, missing practical advice: "with its recorded dependencies available" does not say which versions to install or where they are recorded.
- 697, concept and next question: the reader will notice that `OTU = training$y` holds 0/1, not read counts, and ask what `threshold` does here; also "latent species-trait structure" (`n_lattrait`) is not explained in reader's terms.
- 714-730, choice without reasoning: "residual-based initialisation", `jitter.var`, `maxit` and `sd.errors = FALSE` are unexplained; the last matters because section 15 later relies on gllvm standard errors (1423 says the fit was replayed to get them).
- 735, incomplete example: "Select the Python environment appropriate for your machine" gives no command; the Mojo switch is irrelevant to this reader.
- 742-762, choice without reasoning: `sampling`, `step_size`, `lambda = 0` alongside `weight_decay` are unexplained.
- 765, incomplete example: the selected sjSDM fit includes a 1,000-epoch continuation that is not in the chunk and lives in a `dev/` script excluded from the package, so the reader running the chunk does not reproduce the selected fit; either show it or say plainly that the chunk reproduces only the first stage.
- 806-831, output not read: the runtime table (rendered: occJSDM 18.6 s, Hmsc 24.0 s, gllvm 0.6 s per selected fit; sjSDM 664 s selected and about 8,000 s across 18 attempts) is the practical answer to "how long will this take", and the prose never reads it.
- 833, audience: the parser-error and driver-snapshot history, and the unrelated integration-accuracy sentence, are audit notes for the dev README.

### 11. What this lesson establishes, and what remains open (835-847)

Teaches (summary and exercises), but sits mid-lesson.

- 837-839, missing practical advice: there is no take-home for the eDNA user (for example: on presence/absence data occJSDM's ecological model performs like the established packages, so the choice can rest on detection modelling); add one.
- 839, structure: the paragraph signposts sections 12-15 from what reads as a conclusion; see section 4 of this report.

### Reproduction record (849-855)

Describes; serves verifiers. Commit snapshots, hashes and the excluded `dev/` directory; the "Return to" navigation at 855 is a lesson ending, in the middle of the lesson.

### 12. Does more data help when ecology becomes harder? (857-1038)

Mixed: the scenario design teaches; the run-status block is verification; the results are never stated.

- 859, choice without reasoning: the 100 training sites being a subset of the 300 makes the comparison paired; say so, since that is why lines connect points in the figure at 964.
- 865-866, concept not quantified: "Neither of the complications below" for the baseline, and "much less common" for rare species, need numbers (is the baseline the pilot's design? what prevalence do the rare species have?).
- 872, reasoning not shown: the 640 combinations are not broken down; one sentence (five ecological settings by two site counts by ten communities by four packages is 400; two site counts by two species counts by with/without traits by ten communities by three packages is 240) lets the reader reconstruct the design.
- 901, audience: the inline "still running" branch is study-management scaffolding; once the study is complete it can be a fixed sentence.
- 903-909, audience and length: seven sentences of acceptance rules (thresholds, gradients, objective agreement) serve verifiers; keep the passed/flagged/failed definitions and "failed is not crashed", move the rest.
- 905, reader's next question: 86 of 180 Hmsc fits are flagged and none of occJSDM's; the reader will ask why, and whether the Hmsc runs were simply too short.
- 907, reader's next question: which scenarios the 14 failed gllvm fits come from is not said.
- 923, concept: "original simulation units" here versus standardised units in the pilot needs one clause.
- 933-941, incomplete example: the comment says "then use the training transformation" but the code does not; and the reader is never shown how the squared column reaches a fit (in occJSDM, add the column to `info` and name it in `occCovariates`).
- 964-1017, output not read: neither the figure nor the paired table is described; the prose never says whether more sites reduced error, in which scenarios, or by how much; that answer only arrives in section 14, in a different metric.
- 1012, part of the output not shown well: the paired table is printed with raw column names, unlike every earlier table.
- 1019-1036, output not read and not identified: the fitted-curves figure has no prose, and its caption does not say it is the curved scenario.

### 13. Do traits explain species responses, and do they improve prediction? (1040-1125)

Mixed: the set-up teaches well; the results are not stated.

- 1051, concept not in reader's terms: "species random slopes" in gllvm is unexplained, and the reasoning is for specialists.
- 1055-1081, output not read: the prediction figure is followed by a general statement (1081) but not by what happened: did supplying traits lower error, and more with 30 species than with 10?
- 1085-1121, output not read: the coefficient-recovery figure is not described: did the intervals for the drought trait include the truth, and did those for the irrelevant trait cross zero?
- 1123, unclear and possibly wrong: "Its coefficient table records the true direction or absence of each relationship and explicitly identifies the different link scale" follows a sentence about gllvm, but the link scale belongs to Hmsc, and no such table appears in the lesson.
- 1040-1125, incomplete example: no trait-model fitting code is shown for any package; an occJSDM reader needs at least the occJSDM call, with traits supplied as in the quickstart (its lines 43 and 63).

### 14. Across communities, are estimates biased? (1167-1217)

Teaches; written as a report with conclusions first.

- 1171, concept not in reader's terms: "These are the September results" and "the older occJSDM-only validation study" are project history; describe the study, not its date.
- 1176, order: the summary uses "coverage" before it is defined at 1221; add a clause.
- 1181, choice without reasoning: the lesson switches from mean absolute error (sections 7 and 12) to RMSE without saying why; the reader will compare the pilot's 6.1-7.1 points with RMSE 7.1-7.4 and be confused.
- 1204, reader's next question and practical advice: rare species overestimated "especially for occJSDM with 100 sites" is the result most relevant to eDNA surveys, where many species are rare; say what drives it if known, and what a user should do (treat rare-species estimates with caution, survey more sites).
- 1214, concept: "original coefficient units" versus the standardised units used elsewhere needs a clause.

### 15. Do 95% intervals contain the truth? (1219-1276)

Teaches.

- 1221, best example is elsewhere: the worked coverage example at 1325 (an interval from 0.2 to 0.5 covers 0.4 but misses 0.7) is clearer than the definition here; move it up.
- 1225, reasoning arrives later: why five fixed settings instead of the 300 test sites is explained only at 1341 ("to keep the interval calculation inspectable"); bring the reason here.
- 1240, output not interpreted: occJSDM's 82% coverage of the first environmental effect at 100 sites is reported as "covered less often"; say what it means (about one interval in five misses, not one in twenty, so the intervals are too narrow) and whether the cause is known.
- 1242, jargon: "native intervals" is undefined here.
- 1256, missing practical advice: having shown that a missing curve ruins coverage, tell the reader how to detect it in real data (compare fits with and without a squared term using the held-out scores of section 9).
- 1274, ambiguous: "comparable trait-effect intervals" could mean similar or able to be compared.
- 1276, missing practical advice: there is no closing take-home for the occJSDM user from sections 12-15.

### Full results, interval widths and methods (1278-1532)

Describes; a reference appendix for verifiers that repeats much of sections 14-15.

- 1285, leftover framing: "We now use the same saved results to distinguish bias, typical error and interval coverage" reads as an introduction to material that sections 14-15 have already presented.
- 1287, 1299, 1325, 1353, repetition: the definitions of bias, RMSE, MCSE, coverage and the missing-entry labels repeat 1181, 1202, 1221 and 1245.
- 1289, repetition: this table repeats the 1183 table with MCSE added.
- 1327, misplaced teaching: the distinction between Bayesian credible and frequentist confidence intervals is useful to this reader and belongs in section 15.
- 1356-1365, code: `calibration_targets` and `baseline_calibration` are redefined, overwriting the hidden definitions at 1144-1151.
- 1415, 1417, 1484, repetition: these repeat 1256, 1254 and 1266, with the same numbers.
- 1421-1446, 1488, audience: covariance transformation, fit replay, Monte Carlo SE streams and the sensitivity filter are method verification.
- 1523-1530, incomplete example and a rendering bug: "so these tables can be filtered ... For example:" is followed by a chunk whose code is hidden by the `echo = FALSE` set at 1132; the committed `.md` (its lines 1637-1651) shows only the output tibble, so the reader never sees the filter they are invited to copy.
- 1532, repetition: the probit-arm and older-study caveats repeat 839, 1171 and 1337.

## 3. Jargon list

Terms used before they are explained to this reader, with the line of first use. Where Lessons 1 or 3 already teach a term, that is noted; those need a cross-reference rather than a definition.

- gllvm, sjSDM, Hmsc (as packages): 38, never introduced
- pure JSDM portion: 38
- local optima: 40, explained at 397
- interval coverage: 40, defined at 1221
- bias: 50, defined at 1181 (signed error at 502)
- standardised: 100 (code comment), explained at 202
- logit scale: 138, never explained here (Lesson 3 line 432 gives log-odds)
- linear environmental predictors: 206
- WAIC: 206, taught in Lesson 3
- latent traits: 210; latent species-trait structure: 697
- variational approximation (VA): 211, never explained
- unconstrained factors: 211
- PyTorch CPU optimisation: 212
- weight decay, penalty: 212, partly explained at 397
- low-step continuation: 212, partly explained at 765
- covariance factor dimension: 212
- probit: 213, partly explained at 215
- site level, random level: 213, partly at 769
- fork, Mojo: 217, never explained
- low-rank: 223
- Bernoulli outcome, linear predictor, likelihood: 226
- loadings, loading rows: 227 (loadings taught in Lesson 3)
- standard normal distribution, dot product: 227
- Monte Carlo error: 283
- level-zero prediction: 326
- integrate over: 326
- global parameters: 328
- condition on: 328
- normal-probit: 366
- warm-up: 370 (the quickstart says burn-in)
- bulk and tail ESS: 389 (Lesson 3 covers)
- EVA: 391, never explained
- objectives, objective scores: 391, 907
- log likelihood: 393
- deterministic optimiser, stationary solutions, curvature check, flat direction: 397
- native (starts, intervals, standard errors): 399
- fixed-grid predictions: 414
- percentage points: first prose use 414, defined at 505
- binary JSDM mode: 697
- residual-based initialisation: 714
- epochs, learning rate: 765
- elapsed wall time: 806
- parser error, driver file, input hashes, immutable driver snapshots: 833
- probit-generating scenarios, sensitivity to priors: 839
- snapshot (commit): 851
- maximum absolute gradient: 907
- species random slopes: 1051
- Monte Carlo standard error (MCSE): 1202 in prose (defined there)
- covariance matrix, covariance-matrix warning: 1270
- equal-tailed intervals: 1332
- normal approximation, 1.96 standard errors: 1333, 1421
- conditional logit coefficients: 1351

## 4. Length and structure

The length is mostly code, not prose: about 899 of 1,532 lines are inside chunks, and about 177 lines are prose. A teaching pass that adds a sentence or two to some 50 paragraphs would add perhaps 60-80 lines, roughly 5%. That alone does not make the lesson unreadable. The structure problem is different: the file holds two lessons and a reference appendix for three audiences, and the teaching pass would add prose to all three unless the audiences are separated first.

**Where it repeats itself**

- "No general ranking of packages" or equivalent: 40, 839, 909, 1125, 1177, 1254, 1274, 1321, 1417, 1488. Once in section 11 and once at the end of section 15 would do.
- The older occJSDM-only study not being pooled: 538, 1171, 1532.
- Species-by-site cells are not independent replicates: 1038, 1202, 1299, 1353.
- Bias and RMSE defined at 1181 and again at 1287; MCSE at 1202 and 1299; coverage at 1221 and 1325; missing-entry labels at 1245 and 1353.
- The same results stated twice: curved-scenario coverage (1256 and 1415), rare and correlated coverage (1254 and 1417), the wide gllvm rare-species interval (1266 and 1484), sjSDM's conditional SEs (1272 and 1423), the probit arm still needed (839, 1337, 1532).
- The same tables twice: baseline prediction bias and RMSE (1183 and 1289); baseline coverage (1227 and 1370-1382).
- Why sampled-site errors are larger: 496 and 540.
- Marginal versus conditional probability is taught again (296-298) after Lesson 3 taught it (Lesson 3 line 817); this repetition is defensible, but it should cite Lesson 3.

**Sections that serve a different audience (package verification or optimisation diagnostics, not ecology)**

- 217 (sjSDM fork and Mojo)
- 326-328 (how each package's prediction was made marginal)
- 391-452 (gllvm EVA rejection and the two sjSDM optima: the teaching point at 428 and 450 serves ecologists, the evidence at 397-426 and 440-448 serves verifiers)
- 804-833 (runtime accounting and the parser-error history)
- 849-855 (Reproduction record)
- 901-909 (run status and acceptance rules, partly)
- 1278-1532 (Full results), especially 1419-1502

**Structural changes needed for the prose pass to work**

- Required, small: move the exercise "predict one species given another" (231-285) to after section 4, since it uses the gllvm parameter extraction that section 4 introduces and refers to it as "above" (235).
- Required, small: move the Reproduction record (849-855) and section 10 with its runtime subsection (691-833) to the end of the lesson as appendices, so that section 11 can be the bridge between the pilot and the extension. This resolves the existing TODO item on section order.
- Recommended, the main decision: split the lesson in two. Part A (opening to section 9, then section 11 as its conclusion) is a self-contained teaching lesson on one community, about 800 lines including code. Part B (sections 12-15) is a simulation-study report in a different register, with hidden code and project-history references. Kept together, a reader finishing section 11 has met a conclusion, a reproduction record and navigation links, and then faces a second lesson of similar length. Split, each part can be taught at its own depth.
- Recommended: move "Full results, interval widths and methods" (1278-1532) out of the lesson into a separate reference article, or reduce it to a short appendix. Move the useful teaching in it (the coverage example at 1325, the credible-versus-confidence paragraph at 1327, the five-settings reasoning at 1341) into section 15 first.
- Convention: three pipe tables (208-213, 289-292, 863-868) break the repo's no-pipe-tables rule and must become bullet lists in the pass.

**Line budget if the lesson stays as one file** (current lines, then target; code that is moved is counted at its destination)

- Opening (34-89): 56, keep at about 56; line 40 shrinks, the motivation sentence is added.
- 1 (91-202): 112, keep.
- 2 (204-217): 14, keep at about 14 once the table becomes bullets and 217 moves.
- 3 (219-285): 67, keep, exercise moved after section 4.
- 4 (287-366): 80, about 72; 326-328 shrink to one sentence.
- 5 (368-452): 85, about 40; keep 370-391 and a short version of the two-optima lesson (the "hills" image, 428, 450), move 397-426 and 440-448 to the appendix.
- 6 (454-496): 43, keep.
- 7 (498-614): 117, keep; drop the aside at 538.
- 8 (616-668): 53, keep.
- 9 (670-689): 20, grows to about 26 with the reference scores.
- 11 (835-847): 13, keep, as the bridge.
- 12 (857-1038): 182, about 150; 903-909 shrink to three or four lines.
- 13 (1040-1125): 86, keep, plus result sentences.
- 14 and 15 (1167-1276): 110, keep.
- Appendix: fit it yourself and runtimes (691-833): 143, about 115 with 833 moved to the dev README.
- Appendix: reproduction record (849-855): 7, keep.
- Full results (1278-1532): 255, about 100 if kept, with the repetitions listed above removed; or 0 if moved to a reference article.

Net target as one file: about 1,250 lines before the teaching additions, about 1,320 after.

## 5. Copyedit observations

**Long sentences** (counted outside code spans)

- 221 (48 words): the recommender description.
- 227 (36 words): "The optimised packages here average over ... at new sites."
- 283 (40 words): the fitted-gap sentence with a semicolon.
- 905 (38 words) and 907 (38 words): the acceptance rules.
- 806, 1337 and 1341 are long paragraphs of medium-length sentences that would read better split.

**Slash compounds, the "and/or" family**

- "presence/absence" (97) is a term of art and can stay; "package/dataset/model" (872), "Rhat/ESS" (905, 1335), "bulk/tail ESS" (389, 905), "input/source hashes" (851), "grid/coefficient" (1523) and "gllvm/sjSDM" (1532) should be written out.

**Asides that repeat or distract**

- 538: "These are measured results, not notional examples, and they are not the older occJSDM-only sample-size experiment."
- 666: "they illustrate how to read a curve, rather than determining which results we include."
- 833: the final sentence on integration accuracy has nothing to do with the rest of the paragraph.
- 298 and 507: "For an explicitly hypothetical example" and "As an explicitly hypothetical example", the same verbal tic twice.

**Other wording**

- 206: "All four receive the same information, linear environmental predictors and two hidden factors" needs a colon, and the models do not receive hidden factors (95 says they never do); they estimate two.
- 326: "The inspected sjSDM environmental prediction does too" ("inspected" is vague).
- 599: "a distance from matching simulated truth" ("matching" is ambiguous).
- 1123: the referent of "Its" is unclear (see section 2).
- 1274: "comparable" is ambiguous.
- 391 and 496 use curly quotes, while the rest uses straight quotes or none.
- 87-88: the code comment explains a build-session issue to the reader.
- 1415: four long inline `subset()` expressions where the `report_pct()` helper used at 1256 would do.
- Double blank lines at 614-615 and 1281-1282.

**Factual claims I doubt**

- 52: "The code is visible in the worked lesson; sections 14-15 present the saved results as a report." The `echo = FALSE` at 1132 also hides every chunk in the Full results appendix, including the "For example" chunk at 1525, as the committed `.md` confirms.
- 206: says the models receive two hidden factors, which contradicts 95.
- 221: "The dimensions turn out to be genres that nobody labelled" overstates how interpretable latent dimensions are, and conflicts with 95.
- 226: "A recommender regresses 1s and 0s as if they were continuous scores" is true of the simplest factorisations but not of logistic or implicit-feedback recommenders; qualify it as "a basic recommender".
- 226: "sensible behaviour for rare species" is in tension with 1204 and 1254, where rare species give the worst bias and coverage for occJSDM.
- 227: "the dot product of their two loading rows", whereas the code (241, 245, 348) stores one column per species; harmless if the prose means a species-by-factor layout, but it should match the code.
- 235: "the marginal example above" is below, at 332.
- 283: "The fitted gap is somewhat smaller than the true gap" understates it; in the rendered table the gap is 3.5 points fitted against 6.3 true at the mean, and 7.2 against 10.6 at +1 SD, roughly half to two thirds.
- 765: "the weak penalty is what made the two local maxima identifiable" is a strong causal claim; 397 shows the two maxima exist under the penalty, but the lesson does not show that the unpenalised disagreement (393) came from unbounded coefficients rather than from the same two basins. Doug should confirm.
- 1123: "Its coefficient table records the true direction ... and explicitly identifies the different link scale" refers to a table that does not appear in the lesson.
- 1176: "Coefficient coverage is more variable for occJSDM, gllvm and sjSDM", but more variable than what? Hmsc is excluded from coefficient coverage, and gllvm's baseline is 94% at 100 sites.

## 6. Estimate

- Paragraphs needing a teaching addition: about 55 (opening 5, section 1 3, section 2 4, section 3 6, section 4 5, section 5 7, section 6 1, section 7 3, section 8 2, section 9 2, section 10 7, section 11 1, section 12 9, section 13 4, section 14 4, section 15 5; the Full results appendix adds more if kept). Many are "the output is printed and never read", which needs one or two sentences each.
- Sentences needing copyedit: about 40, including 11 factual points above.
- Size: **large**. The single-community pilot alone is a medium pass (the concepts are there; the gaps are unread outputs and unexplained choices). The extension and appendix make it large, because their results are mostly not stated in prose, and the pass cannot be done well until Doug decides the structure in section 4.

## 7. Questions for Doug

- Split or keep whole: should sections 12-15 become a separate lesson or article (Part B), with sections 1-11 as Lesson 4 proper? Or should the lesson stay one file with the appendices moved to the end?
- Full results appendix (1278-1532): keep in the lesson, reduce to a short appendix, or move to a separate reference article for verifiers?
- The two sjSDM optima (395-452): is this ecology teaching (a converged optimiser is not a unique answer) or a package-verification record? The pass would keep the teaching point either way; the question is whether the evidence stays.
- Who is section 10 for: readers who will rerun the study, or readers learning to fit each package? The answer decides whether its arguments need explaining or only pointing to the dev scripts.
- Why gllvm, not occJSDM, for the conditional-prediction exercise (237-281) and the marginal extraction (332-364)? Is there an occJSDM route the lesson should show instead or alongside?
- Does `runOccJSDM()` standardise covariates itself, and does it offer conditional or marginal prediction at new sites? The prose additions at 202, 285 and 332 depend on the answer.
- Rare-species overestimation and the 82% coefficient coverage for occJSDM at 100 sites (1204, 1240, 1254): is the cause known (prior, shrinkage, sample size), and how much should the lesson say about it?
- Hmsc's 86 flagged fits (905): were the runs too short for that scenario, and should the lesson say so?
- 765: is it established that the weak penalty, rather than the same two basins, explains the unpenalised starts' disagreement?
- Lesson 3 already teaches marginal versus conditional probability and WAIC: should Lesson 4 re-teach them (296-298) or cross-reference them?
- The planned renumbering (TODO line 116) makes this Lesson 3; should the pass wait for the rename, or go ahead, given that both touch every cross-reference?

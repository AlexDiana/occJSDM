# Teaching lessons: status, plans and notes

Last updated: 22 September 2026.

This is a working document for Doug and Alex. Edit it as the lessons develop, and add ideas or corrections in the notes section below. It is kept in Git but explicitly excluded from the source package by `.Rbuildignore`. It is not a vignette and does not need knitting.

The aim is to help an empirical ecologist understand what occJSDM does, read and learn from the R code, and judge its answers by comparing them with known simulated truth. A completed teaching example illustrates one dataset; it does not establish that the model is unbiased or ready for release.

## At a glance

| Lesson | Status | What it teaches | Source |
|---|---|---|---|
| Unnumbered Quickstart | Built; merged in PR #12 | Short fitting example, lesson navigation and where to start. Replaces the old standalone output tour. | [occJSDM.Rmd](occJSDM.Rmd) |
| Lesson 0: Create and explore a simulated survey | Built on main, including unbalanced-data extension | Optional introduction to sites, samples, primers, PCR replicates, covariates, traits and maps; removing whole samples while preserving paired rows. | [occJSDM-lesson-0.Rmd](occJSDM-lesson-0.Rmd) |
| Lesson 1: Fit the model and compare its answers with truth | Built on main, including fitting reference and unbalanced fit | Perfect-observation and PCR fits, occupancy errors, false positives, prior sensitivity and unequal field replication. | [occJSDM-lesson-1.Rmd](occJSDM-lesson-1.Rmd) |
| Lesson 2: Spatial landscapes and dispersal | Outline already on main; no worked spatial results | Smooth environmental gradients, additional spatial structure, contrasting dispersal and prediction away from sampled sites. | [occJSDM-lesson-2.Rmd](occJSDM-lesson-2.Rmd) |
| Lesson 3: Understand the model's outputs | Built, including native plots and the independent-site prediction comparison | Environmental and trait effects, response curves, species associations, ordination, variation partitioning, collection effects and detection effort, all with matching truth. | [occJSDM-lesson-3.Rmd](occJSDM-lesson-3.Rmd) |
| Lesson N: Compare four JSDMs | Built from the four-package pilot; sjSDM stability remains provisional | Compare the pure JSDM in occJSDM with gllvm, sjSDM and Hmsc using perfectly observed presence/absence. | [occJSDM-lesson-N.Rmd](occJSDM-lesson-N.Rmd) |

Read Lesson 0 if the data structure is unfamiliar, then Lesson 1 and Lesson 3. The unfinished spatial lesson is not a prerequisite for Lesson 3. Lesson numbering remains provisional.

## What has been built

The completed lessons use a shared non-spatial community with 100 sites, 10 species, two measured environmental covariates, two measured traits, two field samples per site, two primers and six PCR replicates per primer. The simulation also contains hidden site conditions and an unmeasured species trait. These are different sources of variation. The environmental values and hidden site conditions in this example were generated independently of coordinates: its maps show the sampled points, not smooth spatial surfaces.

Lesson 0 shows how the simulation is constructed and what the model receives. Lesson 1 compares fits supplied with actual presence/absence or PCR observations, explains both the direction and the average size of occupancy errors, and walks through identified weak true detections, strong true detections, laboratory false positives and field-stage false positives. Its alternative-prior example changes the field and laboratory contamination priors together; it cannot attribute the difference to just one of them.

Lesson 3 replaces the old output tour with comparisons against truth, practical diagnostics and supporting tables. The first migration added the native latent-presence table and nine native plots with matching truth: environmental and collection coefficients, baseline occupancy and collection rates, both primers' true/false-positive rates, residual correlations with uncertainty, and cumulative survey counts under different sampling effort. Estimates are paired with the correct generating values, including the scale conversions needed for trait and collection coefficients. It retains examples of poor recovery. In particular, it explains how a real measured-trait effect can be obscured by an opposing unmeasured-trait contribution in a community of only ten species. This is an explanation of the current simulation, not a diagnosis of every older fit.

The follow-up adds three truth-aligned native ordination figures, two native environmental gradients, three separate detection/false-positive rate plots, both native trait-effect plots, and the original bibliography. It also predicts 300 independent new sites and compares one-factor and two-factor fits to identical PCR training observations. The average absolute error against true marginal probabilities is about ten percentage points for both models. The tiny score difference does not establish a preferred factor count. Native new-site prediction intervals and marginal point predictions are compared with their respective correct probability targets.

The teaching R code is visible in HTML and Markdown. Long fitting commands are displayed but do not run during knitting. Compact, checked results are stored in `teaching-data/`; complete fits are archived separately. The original `sampledata` and `sampleresults` objects remain unchanged because they do not contain the full matching truth needed for these lessons.

The canonical lesson sources are the `.Rmd` files. Rendered `.md` files and figures are tracked for reading on GitHub; HTML is generated locally. Reproduction commands, provenance and numerical checks are in the [teaching build README](../dev/simstudy/vignette-lesson/README.md). The [original teaching design](../dev/simstudy/vignette-lesson/DESIGN.md) records the rationale and Paper2Agent assessment, but its proposal-era status statements are historical; use this document for the current lesson status.

## What is planned

### Lesson 2: space, habitat and dispersal

Implement the existing outline after the reviewed spatial correction is available and its source revision can be recorded. The outline was written while PR #8 was awaiting review; check its current status before starting. A merge by itself does not validate a new simulation or lesson.

1. **2A: Smooth environmental gradients.** Build environmental surfaces first, then generate species distributions from them. Show that geography can matter through measured habitat even without an extra spatial process.
2. **2B: Additional spatial structure.** Add a known spatial contribution representing unmeasured conditions. Compare models with and without spatial effects using the same observations. Show true, fitted and difference maps, plus matching variation partitioning.
3. **2C: Contrasting dispersal.** Use a separate, explicit movement or colonisation simulation for species with different dispersal abilities. Keep habitat and observation conditions comparable. This requires simulator work: changing a spatial Gaussian-process range is not the same as simulating dispersal, and the inspected spatial model uses a shared range rather than estimating a dispersal rate for each species.
4. **2D: Prediction beyond sampled sites.** Reserve sites and spatial blocks before fitting. Show interpolation, extrapolation and changes with distance from observations. Compare predictions with the appropriate known probabilities, keeping held-out observations out of fitting and preprocessing.

The detailed questions and constraints are already in the [Lesson 2 outline](occJSDM-lesson-2.Rmd). If the dispersal simulator cannot calculate probabilities directly, repeated independent simulations will be needed to estimate them, with their own Monte Carlo uncertainty shown separately.

### Lesson N: occJSDM, gllvm, sjSDM and Hmsc

The accepted approach is to begin with a small non-spatial comparison, then decide whether to expand it. The [detailed pilot design](../dev/simstudy/jsdm-package-comparison/DESIGN.md) was merged in PR #12, including the decision to retain the existing Mojo toolchain. The four-package pilot has now been fitted and checked on the agreed shared dataset. See the [fitting report](../dev/simstudy/jsdm-package-comparison/FITTING-REPORT.md), including saved truth comparisons, multiple-start failures, and the remaining sjSDM optimisation check. The [student-facing Lesson N](occJSDM-lesson-N.md) now teaches the shared simulation, fitting calls, prediction targets, actual signed and absolute errors, and true-versus-fitted environmental curves for all ten species. It renders saved results without refitting. Doug redirected work to teaching on 22 September; further sjSDM optimisation is parked, and the original selected sjSDM result remains explicitly provisional.

The completed pilot gave all four packages the same perfectly observed presence/absence data: 100 training sites, 300 independent test sites and 10 species. Use an independent simulator with two environmental covariates and two hidden site factors. Exclude traits, phylogeny, space and observation error from this first comparison so that the packages receive the same information.

Distinguish two questions: recovering probabilities at surveyed sites, where the observations help infer hidden conditions; and predicting probabilities at new sites, where those conditions are unknown. Verify that each package's predictions answer the same question before scoring them. Compare probability-scale response curves, signed errors, absolute errors, diagnostics and runtime. Explain the logit/probit difference rather than treating raw coefficients as interchangeable. One pilot is not enough to rank the packages generally.

Use Doug's sjSDM fork v0.2.1 with the PyTorch CPU backend explicitly selected for the main pilot. That release contains both PyTorch and optional Mojo paths. The detailed design pins the source commit and records the environment requirements. Doug decided not to upgrade Mojo; an optional backend comparison would use the existing Mojo 1.0.0 toolchain. Recheck installed package versions when running the experiment rather than silently updating them.

After the pilot works, possible extensions are repeated communities, 100/300/1,000 training sites, both logit- and probit-generating scenarios, and later comparisons of traits, space or community outputs. Set the computing budget and replication plan before launching those extensions.

### Other teaching work still needed

- **Spatial and cross-package prediction.** Lesson 3 now supplies the non-spatial, genuinely independent-site example. Spatial prediction remains in Lesson 2; the four-package pilot is now taught in Lesson N.
- **Choose the number of site-level latent variables using WAIC: unfinished migration work.** The original walkthrough used WAIC to choose the factor count. Lesson 3 currently shows WAIC values for one-factor and two-factor fits, but does not establish a valid selection. First validate or correct the WAIC calculation, including its treatment of unobserved occupancy and collection states. Then compare candidate factor counts using identical observations and otherwise matching settings, and show the criterion, its uncertainty and the selected count beside the known generating count of two. Do not assume that the selected count must equal the generating count. Keep the existing new-site prediction exercise as additional teaching; it does not replace this requested WAIC example.
- **How much species information helps trait recovery.** Lesson 3 shows a problem in one ten-species community. A repeated experiment varying the number of species is a possible follow-up, not an agreed or completed study.
- **Simulator reference.** `simulateOccJSDMData.Rmd` remains a separate reference. Decide later whether to retain it as a technical reference or consolidate overlapping material into Lesson 0 and Lesson 2.

### Variation partitioning among sites

Doug requested a source-code assessment of sjSDM's approach. The [implementation assessment](../dev/simstudy/site-variation-partitioning.md) and [reproducible algebra checks](../dev/simstudy/site-variation-partitioning-audit.R) are complete. They do not add a public function or change the existing species-level output.

The proposed new question is how much environment, space and residual associations help explain the whole observed community at each site. Start with perfectly observed presence/absence, then integrate occupancy and collection uncertainty when scoring PCR observations. Retain species intercepts in every model, use a consistent joint site likelihood, and divide shared improvement transparently. The inspected sjSDM implementation has allocation and edge-case problems, so adapt the idea rather than copying its calculation. Decide between a cheap within-fit diagnostic and a refitted-model comparison before implementation; these answer different questions.

### Paper2Agent companion

A feasibility assessment exists in the [original teaching design](../dev/simstudy/vignette-lesson/DESIGN.md); no occJSDM Paper2Agent integration has been installed or tested. The proposed first pilot would expose a small verified non-spatial workflow through the original R code: inspect a saved fit, check diagnostics and produce comparisons with truth, then demonstrate a small fresh fit within a declared runtime budget.

Assess faithful computation, whether changed inputs really produce new calculations, and whether explanations preserve the distinctions taught in the lessons. Successful execution alone would not establish unbiased ecological inference. Broader tools, hosting and analysis of user datasets remain later decisions.

## Teaching decisions and preferences to preserve

- Explain things for an empirical ecologist. Prefer a longer clear explanation to compressed statistical terminology.
- Use readable tidyverse code, descriptive names and blank lines between steps. Show the teaching code in the knitted HTML.
- Use **variation partitioning** in prose and figures; retain existing function names where the code requires them.
- Pair every ecological estimate with the corresponding simulated truth. Diagnostics and WAIC do not have a biological true value to invent or overlay.
- Distinguish occupancy probability, actual presence/absence, DNA collection and PCR detection. Match units, transformations, sites and species before comparing values.
- Explain how the assumption of good field and laboratory practice enters through priors. Show actual weak true, strong true and false-positive cases; do not promise that replication or high read counts make classification certain.
- Keep at most **six PCR replicates per primer per field sample**.
- For spatial teaching, include smooth environmental gradients and contrasting dispersal, introduced in separate sublessons where needed.
- Report signed error and absolute error separately. Keep poorly recovered examples and avoid choosing seeds because they give attractive results.
- Keep fitting separate from rendering, retain complete simulation truth, and record source versions, seeds, settings and diagnostics.

## Migration checklist for the original function walkthrough

The audit compares `vignettes/occJSDM.Rmd` at revision `8654ff1` with the teaching lessons. Lesson 3 is a rewrite and redistribution, not a complete transfer of every worked example. Listing a function in its index is not the same as teaching its use. Preserve useful workflows while replacing stale argument names, unmatched examples and unsupported claims.

The original walkthrough is preserved unchanged from that revision as [ORIG_occJSDM.Rmd](ORIG_occJSDM.Rmd). It is an archival reference, excluded from package builds. The current `occJSDM.Rmd` remains the Quickstart and lesson guide.

| Original content | Current home and status | Remaining work |
|---|---|---|
| Input structure and fitting arguments | Lessons 0 and 1 plus the Quickstart show current fitting code. Lesson 1 now includes a concise reference for model inference, identifiers, covariates, traits, priors, threshold, missingness and retained latent output. | Keep it aligned with the current API, including `n_lattrait` for fitting and the limited scope of `summarisedLatentPresences = FALSE`. |
| Unbalanced study design created by dropping whole field samples | Lesson 0 removes one whole sample at each of three declared sites, retaining 197 samples and 2,364 PCR rows. Lesson 1 presents a separate matching fit and its diagnostics against unchanged truth. | Keep the example limited to handling unequal effort. One removal pattern and fit cannot estimate the general consequences of sample loss. |
| MCMC settings | Lesson 1 explains chains, burn-in, retained iterations and thinning. | Keep aligned with the actual fitting API; do not recommend thinning as a repair for poor mixing. |
| Parameter-level fitting diagnostics | Restored in this follow-up to Lesson 3: public extraction, parameter-label guide, flagged rows including missing diagnostics, and comparison with newer coefficient diagnostics. | Verify future changes to diagnostic definitions; do not describe the public table as covering every model quantity. |
| Collection-covariate and primer traceplots | Restored in this follow-up to Lesson 3, with actual complete chain excerpts, code for a full fit, and correctly scaled truth lines. Also includes a flagged field-contamination example. | Retain these practical examples when reorganising the lesson. |
| Environmental and collection coefficients | Lesson 3 includes truth comparisons and worked native `plotOccupancyCovariates()` and `plotCollectionCovariates()` examples. | Retain the standardization explanation and truth crosses. |
| Baseline occupancy and environmental response curves | Lesson 3 includes native baseline plots and both native `plotOccupancyGradient()` curves with correctly matched truth. | `plotCovariateEffect()` has a reproduced link/scale defect and cannot yet supply a valid example. Its correction is awaiting review in [PR #13](https://github.com/AlexDiana/occJSDM/pull/13); revisit the example after that fix is merged. Keep zero-factor profiles distinct from marginal probabilities. |
| Conditional and fitted ecological occupancy probabilities | Worked cases, tables and maps in Lesson 1; definitions in Lesson 3. | No additional migration needed for the main distinction. |
| Prediction at new sites | Lesson 3 predicts 300 independent sites; shows the public call, actual-condition truth for native intervals, and marginal truth for posterior mean point predictions. | Spatial prediction remains in Lesson 2; the four-package comparison remains in Lesson N. Do not claim this validates joint community predictions. |
| Ordination scores, loadings and biplot | Lesson 3 now includes all five native return/plot functions, draw-wise truth-assisted rotation/reflection, matched truth overlays and the invariant combined contribution. | No additional migration needed. Preserve the simulation-only alignment caveat and the explanation that native circles are not joint credible regions. |
| Collection, detection and false-positive rate plots | Lesson 3 now includes collection, combined primer-specific rates, and all separate field-FP, laboratory-FP and detection plots with truth. | No additional migration needed. Preserve species/primer matching and read-threshold-adjusted laboratory truth. |
| Trait-by-environment effects | Expanded truth comparisons and cancellation explanation in Lesson 3, plus both native trait-coefficient plots with standardized truth. | A larger species-count investigation is optional future work, not a missing migration step. |
| Residual-correlation display | Lesson 3 now includes the native heat map, its interval-crossing-zero markers and true values matched to all 45 species pairs. Undefined true correlations remain explicit. | No additional migration needed for the native uncertainty example. |
| Variation partitioning, including ternary display | Lesson 3 compares environmental and residual shares in a non-spatial fit. | Add all three components and the corresponding display with a validated spatial example in Lesson 2. |
| Cumulative detections with different PCR and field replication | Lesson 3 now shows the native plot in both directions beside exact generating survey-count quantiles, plus the existing analytic-expectation example. | Retain the distinction between uncertainty in an expectation and variation in random survey counts. Both use at most six PCRs per primer. |
| WAIC selection of the number of site-level latent variables | **Unfinished.** Lesson 3 shows current WAIC values and their limitation. Its independent-site prediction comparison is an additional exercise, not a replacement for factor-count selection. | Validate or correct the criterion, then add the matched factor-count comparison with the known generating count and uncertainty. Replace the original hard-coded ranking with checked results. |
| `returnLatentPresences()` and its coloured table | Restored in Lesson 3 using the current fit. The native table is joined to actual states and generating probabilities by species/site/sample/primer/PCR identity. | Preserve the distinction between actual binary states, conditional state probabilities and generating probabilities. |

The migration now includes the native table, remaining valid native plotting examples, independent-site prediction, comparison of predictive performance, and original bibliography. **Migration remains unfinished:** the WAIC-based choice of factor count still needs validation and a worked example, and Lesson 2 still needs the spatial three-component/ternary example. The `plotCovariateEffect()` correction is awaiting review in [PR #13](https://github.com/AlexDiana/occJSDM/pull/13) before its original workflow can be taught as valid; Lesson 3 currently uses checked occupancy-gradient plots. Keep the existing prediction results, figures and explanations: their scientific content does not need rewriting to correct this completion status. Lesson N now provides the separate four-package pilot. A replicated benchmark, site-level variation partitioning implementation and the Paper2Agent pilot remain future projects.

## Doug and Alex's notes

Add free-form notes below. These can be questions, suggested examples, wording changes or decisions. A note is not automatically an agreed implementation task. When an idea is adopted, move it into the relevant plan above and record the decision; leave unresolved questions here.

### New notes

<!-- Add notes here. An optional format is: date; lesson/topic; note; next step or decision needed. -->

_Add notes here._

### Decisions made after this update

<!-- Record the date, decision and affected lesson here, then update its status above. -->

- **21 September 2026:** Doug requested pull request review of Lessons 3, N and 2. The review combines the completed Lesson 3, Quickstart and planning document with the Lesson N design. The Lesson 2 outline was already on main; the proposed changes repair its navigation. The planned spatial and four-package experiments remain unfinished.
- **21 September 2026, after PR #12 merged:** Doug asked whether the original vignette's content had all been retained. The audit identified missing practical instructions. He approved restoring hands-on fitting diagnostics first and recording every remaining migration gap above.

- **22 September 2026:** Doug approved starting the migration gaps and auditing sjSDM site-level variation partitioning. The first batch adds native tables/plots, fitting guidance and the unbalanced-data workflow. Further prediction, spatial and model-selection experiments remain planned.

- **22 September 2026, follow-up:** Doug approved completing the remaining non-spatial migration. One additional fit changes only the factor count, and 300 new sites were generated before inspecting predictions. Native plot and prediction exports retain all posterior draws. The spatial example remains deferred; newly reproduced output defects are documented, not silently repaired within the lesson branch.

- **22 September 2026, clarification:** Doug confirmed that the requested WAIC comparison is for choosing the number of site-level latent variables. Keep the existing prediction lesson unchanged, mark WAIC-based factor-count selection as unfinished migration work, and add the validated example later. The new-site experiment does not fulfil that requirement.

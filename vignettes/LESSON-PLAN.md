# Teaching lessons: status, plans and notes

Last updated: 21 September 2026.

This is a working document for Doug and Alex. Edit it as the lessons develop, and add ideas or corrections in the notes section below. It is kept in Git but explicitly excluded from the source package by `.Rbuildignore`. It is not a vignette and does not need knitting.

The aim is to help an empirical ecologist understand what occJSDM does, read and learn from the R code, and judge its answers by comparing them with known simulated truth. A completed teaching example illustrates one dataset; it does not establish that the model is unbiased or ready for release.

## At a glance

| Lesson | Status | What it teaches | Source |
|---|---|---|---|
| Unnumbered Quickstart | Built; merged into main | Short fitting example, lesson navigation and where to start. Replaces the old standalone output tour. | [occJSDM.Rmd](occJSDM.Rmd) |
| Lesson 0: Create and explore a simulated survey | Built; already on main | Optional introduction to sites, samples, primers, PCR replicates, environmental covariates, traits and maps of the simulated world. | [occJSDM-lesson-0.Rmd](occJSDM-lesson-0.Rmd) |
| Lesson 1: Fit the model and compare its answers with truth | Built; already on main | Perfect-observation and PCR fits, occupancy errors, false positives, the good-practice assumptions and sensitivity to contamination priors. | [occJSDM-first-lesson.Rmd](occJSDM-first-lesson.Rmd) |
| Lesson 2: Spatial landscapes and dispersal | Outline already on main; no worked spatial results | Smooth environmental gradients, additional spatial structure, contrasting dispersal and prediction away from sampled sites. | [occJSDM-lesson-2.Rmd](occJSDM-lesson-2.Rmd) |
| Lesson 3: Understand the model's outputs | Built; merged into main | Environmental and trait effects, response curves, species associations, ordination, variation partitioning, collection effects and detection effort, all with matching truth. | [occJSDM-lesson-3.Rmd](occJSDM-lesson-3.Rmd) |
| Lesson N: Compare four JSDMs | Design merged into main; no comparison fits or lesson yet | Compare the pure JSDM in occJSDM with gllvm, sjSDM and Hmsc using perfectly observed presence/absence. | [Pilot design](../dev/simstudy/jsdm-package-comparison/DESIGN.md) |

Read Lesson 0 if the data structure is unfamiliar, then Lesson 1 and Lesson 3. The unfinished spatial lesson is not a prerequisite for Lesson 3. Lesson numbering remains provisional.

## What has been built

The completed lessons use a shared non-spatial community with 100 sites, 10 species, two measured environmental covariates, two measured traits, two field samples per site, two primers and six PCR replicates per primer. The simulation also contains hidden site conditions and an unmeasured species trait. These are different sources of variation. The environmental values and hidden site conditions in this example were generated independently of coordinates: its maps show the sampled points, not smooth spatial surfaces.

Lesson 0 shows how the simulation is constructed and what the model receives. Lesson 1 compares fits supplied with actual presence/absence or PCR observations, explains both the direction and the average size of occupancy errors, and walks through identified weak true detections, strong true detections, laboratory false positives and field-stage false positives. Its alternative-prior example changes the field and laboratory contamination priors together; it cannot attribute the difference to just one of them.

Lesson 3 replaces the old output tour with ten figures and supporting tables. Estimates are paired with the correct generating values, including the scale conversions needed for trait and collection coefficients. It retains examples of poor recovery. In particular, it explains how a real measured-trait effect can be obscured by an opposing unmeasured-trait contribution in a community of only ten species. This is an explanation of the current simulation, not a diagnosis of every older fit.

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

The accepted approach is to begin with a small non-spatial comparison, then decide whether to expand it. The [detailed pilot design](../dev/simstudy/jsdm-package-comparison/DESIGN.md) is now merged into main, including the decision to retain the existing Mojo toolchain. Implementation and comparison fits remain to be done.

The pilot will give all four packages the same perfectly observed presence/absence data: 100 training sites, 300 independent test sites and 10 species. Use an independent simulator with two environmental covariates and two hidden site factors. Exclude traits, phylogeny, space and observation error from this first comparison so that the packages receive the same information.

Distinguish two questions: recovering probabilities at surveyed sites, where the observations help infer hidden conditions; and predicting probabilities at new sites, where those conditions are unknown. Verify that each package's predictions answer the same question before scoring them. Compare probability-scale response curves, signed errors, absolute errors, diagnostics and runtime. Explain the logit/probit difference rather than treating raw coefficients as interchangeable. One pilot is not enough to rank the packages generally.

Use Doug's sjSDM fork v0.2.1 with the PyTorch CPU backend explicitly selected for the main pilot. That release contains both PyTorch and optional Mojo paths. The detailed design pins the source commit and records the environment requirements. Doug decided not to upgrade Mojo; an optional backend comparison would use the existing Mojo 1.0.0 toolchain. Recheck installed package versions when running the experiment rather than silently updating them.

After the pilot works, possible extensions are repeated communities, 100/300/1,000 training sites, both logit- and probit-generating scenarios, and later comparisons of traits, space or community outputs. Set the computing budget and replication plan before launching those extensions.

### Other teaching work still needed

- **Genuine held-out prediction.** This is planned within Lesson 2 and Lesson N; the current fitted-site figures do not replace it. Decide whether a separate non-spatial prediction lesson would be useful.
- **Model selection.** A worked comparison needs models fitted to the same response data, with the generating model known and independent predictions assessed. The current perfect-observation and PCR fits cannot be compared by WAIC as if they fitted the same data.
- **How much species information helps trait recovery.** Lesson 3 shows a problem in one ten-species community. A repeated experiment varying the number of species is a possible follow-up, not an agreed or completed study.
- **Simulator reference.** `simulateOccJSDMData.Rmd` remains a separate reference. Decide later whether to retain it as a technical reference or consolidate overlapping material into Lesson 0 and Lesson 2.

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

## Doug and Alex's notes

Add free-form notes below. These can be questions, suggested examples, wording changes or decisions. A note is not automatically an agreed implementation task. When an idea is adopted, move it into the relevant plan above and record the decision; leave unresolved questions here.

### New notes

<!-- Add notes here. An optional format is: date; lesson/topic; note; next step or decision needed. -->

_Add notes here._

### Decisions made after this update

<!-- Record the date, decision and affected lesson here, then update its status above. -->

- **21 September 2026:** Doug requested merging Lessons 3, N and 2. The completed Lesson 3, Quickstart and planning document were merged into main, along with the Lesson N design. The Lesson 2 outline was already on main. This integration does not mark the planned spatial or four-package experiments as completed.

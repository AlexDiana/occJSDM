# Build the truth-based teaching lessons

The teaching sequence has eight lessons in `vignettes/`, numbered 0 to 7: `occJSDM-lesson-0.Rmd` (optional simulation and data orientation); `occJSDM-lesson-1.Rmd` (what occupancy models and joint species distribution models do, worked by hand without a fit); `occJSDM-lesson-2.Rmd` (fitting, truth comparisons and detection examples); `occJSDM-lesson-3.Rmd` (the model's outputs compared with truth); `occJSDM-lesson-4.Rmd` (prediction at new sites and model comparison); `occJSDM-lesson-5.Rmd` and `occJSDM-lesson-6.Rmd` (the four-package JSDM comparison, on one community and then across ten communities with traits); and `occJSDM-lesson-7.Rmd` (spatial landscapes and survey design). Lessons 0 and 2 show their teaching code and use the same original simulation and baseline fits. The unequal-replication extension described below adds one fit to a reduced version of that survey. Their maps display the existing sampled coordinates; the environmental values and hidden site factors are independent of geography in this dataset. No interpolation or spatial fit is implied.

The spatial lesson, `vignettes/occJSDM-lesson-7.Rmd`, is no longer an outline. PR #8 merged on 27 September 2026, and its worked sections now report the site-arrangement sweep, which is built and checked under `dev/simstudy/spatial-design-sweep/`. The spatial model has a shared range, not mechanistic or species-specific dispersal parameters; dispersal, species-specific ranges and same-scale environmental confounding remain deferred contrasts.

The fitting lesson is `vignettes/occJSDM-lesson-2.Rmd`. Its figures and tables come from `vignettes/teaching-data/nonspatial-lesson.rds`, a compact bundle containing the complete simulation, generating settings, identity mappings, selected cases, posterior summaries, diagnostics and full-fit hashes. Full MCMC objects remain in the build directory; they are not shipped with the package. The existing `sampledata` and `sampleresults` objects are unchanged.

The experiment uses main revision b53048a's non-spatial model, including the reviewed collection alignment, RNG and factor-reparameterisation fixes and the read-threshold correction. It does not test the spatial changes in PR #8. The documentation branch changes no code in `R/` or `src/`.

Rebuilt on 2 October 2026 with three field samples per site (previously two): 100 sites, 10 species, three field samples per site, two primers and six PCRs per primer and sample, so 300 samples and 3,600 PCR rows per species. Every seed, prior, MCMC setting, case-selection rule and fit argument is unchanged, and the site-level truth is identical to the two-sample survey. The bundles for Lessons 0, 2 and 3 were all re-exported from one fresh archive, `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples`, built on branch `codex/lesson-bundle-3-samples` from main at 5b058ed; the recorded source commit is eeb1675 (the unbalanced bundle records 093a16d, which has the same `helpers.R`). The package source now includes PR #13's covariate-effect correction. Fit times: perfect 26.5 s, default 78.6 s, alternative 76.5 s, longer alternative 150.5 s, one-factor prediction 72.7 s, unbalanced 76.0 s. Every command and its output is logged in `BUILD-LOG.md` in that archive. The dated validation records below describe the build checked on that date; the survey descriptions and recorded results have been updated to the three-sample build.

## Current rebuild after PR #14 (7 October 2026)

The non-spatial three-sample survey and its dependent bundles were rebuilt under merged revision `a56f554`. The new archive is `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-post-pr14-20261007`; the October 2 archive remains intact. The new archive has a dedicated `library`, a `BUILD-LOG.md`, per-step logs, exact commands in `run-rebuild.py`, the saved fitting session and `fit-comparison.csv`. All six posterior-results objects are exactly identical to those in the preceding archive: perfect observation, default PCR, both alternative-prior runs, one-factor prediction and unequal sampling. Seeds, priors, data and MCMC schedules are unchanged. Thread limits were set to one. Provenance records now refer to the merged source and the fresh fit files.

All existing independent verifiers pass for the eleven rebuilt bundles: nonspatial lesson, outputs, diagnostics, latent-presence table, native plots, aligned ordination, remaining plots, native traits, covariate effects, prediction and unequal sampling. This refresh covers the non-spatial survey used by Lessons 0-4. The four-package studies, spatial study and shipped `sampleresults` remain separate records and release tasks.

The corrected site-WAIC assessment has its own compact bundle and numerical records. Reproduce its two draw checks from the repository root using the new archive's library:

```sh
R_LIBS=/path/to/new-archive/library Rscript dev/simstudy/vignette-lesson/validate-waic.R /path/to/new-archive /path/to/new-output-1000 250
R_LIBS=/path/to/new-archive/library Rscript dev/simstudy/vignette-lesson/validate-waic.R /path/to/new-archive /path/to/new-output-4000 1000
```

Use `new-archive/waic-1000` and `new-archive/waic-4000` as the output directories when exporting the compact lesson bundle. From the repository root, then run `Rscript dev/simstudy/vignette-lesson/export-waic.R ARCHIVE` followed by `Rscript dev/simstudy/vignette-lesson/verify-waic.R ARCHIVE`. Each output directory must be new. Draws are evenly spaced within all four chains; the arguments specify draws per chain. Both runs independently check stricter quadrature and compare the final arithmetic with `loo::waic()`. The validator requires the matching archive library and exact fit-time source hashes. Historical `dev/simstudy/validate_site_waic.R` is unchanged and still reproduces the September layout. See [the site-WAIC report](../site-waic-validation.md) for results and their reliability limits.

## Generate the simulation and fits

Run from the repository root. Use an empty build directory and install this checkout into a dedicated R library so that the loaded package matches the source being recorded. These commands assume the package's dependencies are already installed, along with `posterior`, `testthat`, `rmarkdown` and `knitr`. The ordinary source install compiles the C++ code. Do not silently reuse a different installed occJSDM revision.

```sh
mkdir -p /tmp/occjsdm-lesson/library
R CMD INSTALL --preclean --library=/tmp/occjsdm-lesson/library .
export R_LIBS=/tmp/occjsdm-lesson/library
Rscript dev/simstudy/vignette-lesson/test_lesson.R
Rscript dev/simstudy/vignette-lesson/build_lesson.R /tmp/occjsdm-lesson prepare
Rscript dev/simstudy/vignette-lesson/build_lesson.R /tmp/occjsdm-lesson perfect
Rscript dev/simstudy/vignette-lesson/build_lesson.R /tmp/occjsdm-lesson default
Rscript dev/simstudy/vignette-lesson/build_lesson.R /tmp/occjsdm-lesson alternative
Rscript dev/simstudy/vignette-lesson/build_lesson.R /tmp/occjsdm-lesson alternative long
Rscript dev/simstudy/vignette-lesson/summarise_lesson.R /tmp/occjsdm-lesson
Rscript dev/simstudy/vignette-lesson/verify_lesson.R /tmp/occjsdm-lesson
```

The initial alternative-prior fit is retained, but the published example uses its longer fit because an individual field-contamination parameter (OTU_6's `theta0`) had ESS about 186 initially (Rhat 1.025 and ESS about 230 in the two-sample survey). `summarise_lesson.R` uses a saved `*-long-fit.rds` when one is present. The table in the lesson reports the selected fits' actual diagnostics; no scientific error target is imposed.

`prepare` refuses to overwrite an input bundle. The fitting step refuses to overwrite a fit, and checks the source hashes against the input. The verification step checks generator hashes, fit hashes and input identity. To change the simulation or fit source, use a new build directory and regenerate the evidence rather than mix old fits with new truth. Selected cases are saved to `cases-selected-before-fitting.csv` before any fit is run. MCMC seeds are independent of the simulation seed and are fixed per fit.

**Source checks since 4 October 2026.** The 11 verifiers that read saved fits now call `lesson_source_check()` in `source-check.R`. So do the six figure exporters (`export_native_plots.R`, `native-traits-export.R`, `ordination-export.R`, `remaining-plots-export.R`, `covariate-effect-export.R`, `prediction-export.R`). Each previously required `identical(..., lesson_source_hashes())`. The new check ignores comments and compares each changed `R/` file function by function, against the version in git history whose md5 the fit recorded. Only the functions in `lesson_plot_only_changes` may differ: `plotFPTPStage2Rates()`, `plotDetectionRates()`, `plotStage2FPRates()`, `plotTraceplot()` and `okabe_ito()`, added for PR #23's colour scales. It also fails if any other code in `R/` refers to one of them, so they cannot be called while fitting. `DESCRIPTION`, `NAMESPACE` and `src/` must still match byte for byte. A new entry in `lesson_plot_only_changes` needs a dated justification beside it. The fitting, build and summarise scripts keep the strict check on purpose: `build_lesson.R`, `summarise_lesson.R`, `summarise_outputs.R`, `summarise_diagnostics.R`, `summarise_latent_tables.R`, `unbalanced-build.R`, `prediction-build.R` and `remaining-plots-diagnose-covariate.R`. A fit or summary should record exactly the code it ran with. Several bundles also record these scripts' md5 (in `generator_hashes`, `summary_source_hashes` and similar fields), so editing the scripts would break the verifiers. Consequently, after any change to `R/`, plot-only changes included, these scripts stop on the existing archive, and re-running one needs a new build directory and new fits. In the archive, `library` was reinstalled from current main on 4 October 2026, because the figure verifiers compare the installed functions with `R/`; the fit-time build is kept as `library-fbe3ed6`. Re-exported bundles keep the fit-time `source_hashes` and add `export_source_hashes`, the tree they were exported from, which is for information only.

The source categories describe this simulator: a positive with `w=0` is a laboratory false positive; a positive with `z=0,w=1` is a field-stage false positive; a positive with `z=1,w=1` is a true detection. Missing observations are separate. The first sample in each declared category is selected in lexicographic species-name, numeric-site and numeric-sample order. Consequently `OTU_10` sorts before `OTU_2`; this is intentional and does not use fitted results. Category totals and the all-positive-sample table prevent the selected cases from being presented as a general error-rate estimate.

The prior stress test changes both `q` and `theta0` from Beta(1,20) to Beta(1,4), leaving `p` unchanged. It cannot distinguish which of the two altered priors is responsible for a change. It illustrates sensitivity to the joint low-contamination assumption and is not a recommendation to change the defaults.

## Render without rerunning MCMC

The small RDS is tracked so a fresh checkout can render immediately. These commands run from the repository root:

```sh
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-2.Rmd", output_format="rmarkdown::html_vignette")'
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-2.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
```

Use the same commands with `occJSDM-lesson-0.Rmd` and `occJSDM-lesson-7.Rmd` to render the optional introduction and the spatial lesson; every other lesson renders the same way with its own file name. Render each in a fresh R environment when checking standalone execution. Teaching chunks are visible by default; only document-formatting setup is hidden. Simulation and MCMC examples marked `eval=FALSE` are displayed but do not run while knitting. The Lesson 0 simulator call reproduces the existing settings without calling the build helper; the helper remains unchanged to preserve its recorded provenance.

The GitHub Markdown and its PNGs are tracked as readable review artifacts. The HTML is generated locally. The R Markdown source remains canonical. Neither render loads full MCMC fits or starts fitting. To reproduce the original numerical results exactly, use the recorded source revision and R/package versions; different platforms can still introduce small numerical differences. If `library()` warns that a package was built under a newer version of R, knitr writes that warning into the lesson. R only checks this when a package is first attached, so attaching dplyr, tidyr, tibble and ggplot2 quietly before calling `rmarkdown::render()` in the same R session keeps the warnings out.

### Editing links in RStudio

Write links between lessons as ordinary Markdown, for example `[Lesson 2](occJSDM-lesson-2.md)` or `[Diagnostics](occJSDM-lesson-3.md#check-computation-as-well-as-ecological-recovery)`. Do not put inline R expressions inside a link destination: RStudio's Visual editor can encode those expressions as URL text, preventing knitr from evaluating them.

Each of the nine teaching documents (the quickstart and the eight lessons) sources `vignettes/lesson-links.R` in its hidden setup chunk. This shared document hook changes recognised lesson destinations from `.md` to `.html` only in HTML output. It preserves section anchors, external links, other documents and code examples. It does not edit the `.Rmd` source. When adding a lesson, add its filename stem to the helper's explicit list and source the helper in that lesson's setup, and add it to the lesson list at the top of `assets/js/lessons.js`, which gives the GitHub Pages site its links between lessons. Each entry needs the lesson's file stem, a short label, its exact YAML title and a `published` flag, which must agree with the lesson's lines in `_config.yml`. Two checks also list the lessons and will fail until they are updated: the expected labels and position text in `dev/simstudy/lesson-site/test_lessons.js` (run with `node --test`), and the page count in `build_preview.R`. The helper is tracked and included in the source package.

The GitHub Pages site needs one more thing. Pages builds the lessons from the tracked `.md` files and converts their lesson links to `.html` with the jekyll-relative-links plugin, which only recognises a link whose text sits on a single line. pandoc wraps Markdown at 72 columns by default, and a link whose text crossed a line break was left pointing at the `.md` file, which the site serves as raw Markdown: 15 lesson links were broken this way on the live site on 1 October 2026. So always render the GitHub Markdown with `pandoc_args = "--wrap=none"`, as the commands in this README do, which writes each paragraph on one line. On GitHub itself the `.md` links work either way.

Run the focused regression check from the repository root:

```sh
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
```

The check renders a small fixture to HTML, Markdown and HTML again in one R session. It verifies the destinations and anchors, preserves literal examples and unrelated links, and rejects dynamic or encoded lesson destinations in the nine lesson sources. It also reads the nine committed `.md` files and fails if the text of any lesson link wraps onto a second line. No MCMC is run. When changing the link convention, also edit and save disposable copies in RStudio Visual mode and render those copies to both formats; a normal command-line render alone does not exercise the editor's Markdown rewriting.

Validation on 21 September 2026: all five documents rendered to HTML and Markdown. Copies of all five were then opened, edited and saved in the installed RStudio Visual editor and rendered again to both formats. All 19 lesson links retained the correct destinations, and the diagnostic section anchor was present. The focused regression check passed. A source-package build included the helper, excluded the archived original and successfully rendered the Quickstart after unpacking. The existing generated Markdown and figures were unchanged; no model code, simulated data or fits changed.

## What is checked

- Helper tests distinguish source categories despite reordered observations/species, retain missingness and verify the read-rounding correction.
- The evidence check confirms that the binary fit uses the same actual z and that the two-stage fits share their simulation and source hashes.
- Occupancy summaries are reconstructed from posterior model components; the independent check uses a separate component-wise calculation for all sites and draws of OTU_1, including interval endpoints. For both two-stage fits, all 1,000 cell means are also compared with the fitter's saved psi means.
- Reported sample/site probabilities are compared with their saved fitted arrays; all band errors are recalculated.
- Numerical diagnostics are exported for reported probabilities and named parameter blocks. Neither converged sampling nor passing these checks proves unbiased ecological inference.

## Remaining tutorial work

### Teaching-code revision checked on 20 September 2026

Both completed lessons and the Lesson 2 outline render independently to HTML and GitHub Markdown. The HTML embeds `vignettes/teaching.css` for larger, more widely spaced code; format-specific links connect the corresponding HTML or Markdown lessons. Only formatting/link setup is hidden.

The visible Lesson 0 simulation reproduces the saved observations and truth. Its clearer variable names change the incidental row labels of the input `p` matrix, not the values or primer order. The visible tidyverse transformations reproduce all 24,000 observation/source rows, every signed/absolute error summary and the existing perfect-observation input. The optional public extraction example matches all 1,000 default-fit occupancy means. Map joins, both-sample case counts, the 24 helper assertions and the independent full-fit evidence checks pass. The selected cases, model source, simulation bundle, full fits and build-helper hashes are unchanged; no new MCMC was run.

Independent review identified missing truth-matrix labels in the optional fresh simulation call. Explicit site, sample and species labels were added and verified. The new maps and revised figures were visually inspected. Lesson 2 remains an outline: its smooth environmental and explicit dispersal simulations still need to be implemented and validated after spatial review.

### Lesson 3: output interpretation, added 21 September 2026

`vignettes/occJSDM-lesson-3.Rmd` extends the agreed design to environmental and trait effects, response profiles, baseline probabilities, residual correlations, combined ordination contributions, variation partitioning, collection effects and expected detection with additional sampling. `vignettes/occJSDM.Rmd` is now the short quickstart/lesson guide; the old unmatched output tour is replaced by Lesson 3 and the existing Lesson 1. The original example data and fits remain unchanged. Since Phase B of the lesson rewrite (3 October 2026), Lesson 3 no longer covers variation partitioning, which was removed until a working spatial model gives it something to show, or new-site prediction, which is now Lesson 4; the Lesson 1 named here is now Lesson 2.

The new lesson reuses the complete `perfect-fit.rds` and `default-fit.rds` from the same archive. Their model source hashes match current main revision 8654ff1: the changes since the recorded b53048a source are documentation changes. No new MCMC or model changes are needed. Run the following from the repository root, with a matching occJSDM library:

```sh
Rscript dev/simstudy/vignette-lesson/summarise_outputs.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/verify_outputs.R /path/to/full-fits
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-3.Rmd", output_format="rmarkdown::html_vignette")'
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-3.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE, pandoc_args="--wrap=none"))'
```

The exporter writes `vignettes/teaching-data/output-lesson.rds`, including the original fit manifests and hashes of its source, the original compact bundle and the legacy saved example used for the historical interval-count check. It checks actual training observations and reconstructs the full true predictor to establish the environmental scale. True trait effects are converted using trait standard deviations; the simulator standardizes its environmental predictors but not its traits. Changing exporter code requires re-exporting and re-verifying the compact bundle.

Checks independently reconstruct raw coefficient summaries, response-profile endpoints, selected posterior products of scores and loadings, correlation intervals and the non-spatial partition truth. They also verify the trait-decomposition algebra, collection standardization and detection-effort expectations. The effort check uses an independent sum over the possible number of collected field samples. Fits remain in the external archive, not in the package.

Validation on 21 September 2026: the independent numerical checks passed for both fits, including direct chain diagnostics for the environmental and trait coefficients. Lesson 3 and the quickstart rendered to HTML and GitHub Markdown. All ten figures were visually inspected; local links, embedded teaching CSS and visible code blocks were checked. Scientific review found no significant issues. The existing lessons, simulation bundle and model source are unchanged; this documentation change did not rerun the package-wide test suite or MCMC.

Later integration validation on 21 September 2026: Lesson 3 and the Lesson 4 design were combined for pull request review; the Lesson 2 outline was already present on main. The package test suite passed 641 expectations with zero test failures or warnings; the long coverage study remained explicitly skipped. All five teaching documents rendered to HTML and Markdown, and their local links were checked. Encoded inline-R lesson links introduced by RStudio formatting were repaired, and Lesson 1 now directs readers to Lesson 3. A source-package archive retained the lesson sources while excluding `vignettes/LESSON-PLAN.md` and the development plans. Production model code and existing simulation results were unchanged.

Interpretation choices are deliberate:

- Response curves use hidden site factors set to zero, matching the public gradient function. They are not new-site probabilities marginalized over unmeasured conditions.
- Trait effects are distinguished from the realized regression of true species slopes on traits. In this ten-species community, the measured Trait_1 and an unmeasured species trait correlate by chance, partially cancelling their contributions to the first environmental response. This illustrates a difficulty, not a general diagnosis of older fits without saved truth.
- True correlations involving OTU_4's zero-variance factor contribution remain undefined, not zero. They appear grey.
- Ordination is checked using the posterior mean of the draw-wise score/loading product, not products of posterior means or unaligned axes.
- Variation partitioning uses the same probability-standard-deviation allocation as the simulator and fitter, including the environmental intercept. Its residual component does not establish biotic interactions.
- Detection-effort curves describe expected true detections conditional on all ten species occupying the site, mean collection conditions, independent field samples and two primers. They exclude false positives and respect six PCRs per primer. Their intervals describe uncertainty in the expectation, unlike the public cumulative-detection plot's simulated survey-outcome intervals. They do not show the benefit of refitting with more data.

Remaining teaching work: spatial examples after the reviewed correction, genuine held-out prediction, a matched model-selection experiment and the separate package-comparison lesson. Paper2Agent remains a feasibility proposal and has not been installed or run.

### Lesson 3: practical fitting diagnostics restored, 21 September 2026

The diagnostics section now demonstrates public diagnostic extraction, parameter labels, filtering individual warnings (including missing diagnostics), array indexing and traceplots. It distinguishes the public function's classical `coda` Rhat/ESS from the newer coefficient diagnostics calculated with `posterior`, and distinguishes numerical convergence from recovery of ecological truth. `vignettes/LESSON-PLAN.md` records the remaining migration gaps against the original walkthrough at revision `8654ff1`; that planning document remains excluded from package builds.

Run these additional commands from the repository root with the matching occJSDM library before rendering the updated Lesson 3:

```sh
Rscript dev/simstudy/vignette-lesson/summarise_diagnostics.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/verify_diagnostics.R /path/to/full-fits
```

The exporter writes `vignettes/teaching-data/diagnostics-lesson.rds` (about 1 MB). It contains every retained draw in all four chains for the collection slope and primer-1 detection rate of OTU_1 and OTU_6 from `default-fit.rds`, plus OTU_6's field-contamination rate from `alternative-long-fit.rds`. The last example is selected by the lowest public ESS among field-contamination parameters, not by its error against truth. No additional thinning or model fitting occurs. The original simulation and output-summary bundles are unchanged.

The exporter checks the archived fits' source, input and file hashes and reproduces their saved public diagnostic tables. Its independent verifier compares every exported draw and chain with the original arrays, checks every native traceplot's iteration and chain identities, independently recalculates the standardized collection truth and threshold-adjusted PCR truth, and verifies the selection of the flagged field-contamination example. The new bundle records the exporter hash, original teaching-bundle hash and fit manifests; re-export and re-verify it if the exporter changes.

Lesson 3 uses `ggtern::theme_bw()` so both ordinary plots and native occJSDM traceplots continue to render after loading ggtern, including successive renders in one R process. With ggplot2 4.0.3 and ggtern 4.0.0, resetting the ordinary ggplot2 theme after ggtern has registered its theme elements otherwise produces a theme-validation error. This is a lesson-level compatibility choice, with no changes to the plotting package or model code.

Validation on 21 September 2026: every exported draw, chain, iteration and truth line passed the independent archive check. All five optional extraction/diagnostic examples ran against the full saved fits, and all three native traceplots built successfully. Lesson 3 rendered to HTML and Markdown; the three new figures were visually inspected and the original ten figures were unchanged. Scientific review found no significant issues. R reported that dplyr, tibble and ggplot2 were built under R 4.5.2; these startup warnings did not prevent the examples or renders from completing. No new MCMC or package-wide test run was needed for this documentation-only change.

Added 30 September 2026: `vignettes/teaching-data/mirror-labelling-chains.csv` (about 3 KB) holds, for species 6 of `design-qfar_K6-sites300-05`, each chain's posterior mean and SD of `theta0` and of mean occupancy over the 100 scored sites, with the generating values, for the original four-chain fit and the 16 fresh chains of the [convergence-flag diagnosis](../convergence-flag-diagnosis/REPORT.md); Lesson 3's mirror-labelling subsection draws it. Its values are copied as text, with no fit or exporter, from that study's committed results: anatomy/chain-summary.csv (md5 dcd6db73c8bec84fba440a0c9da7ce73), extended/chain-summary.csv (md5 bb83a00e826c4e40cbb3f1123089e8e0), extended/mode-mass.csv (md5 f8ae61ef5a41ad91ad02730907fe7769) and modes/anchor-validation.csv (md5 e034955775dbe08803b49ff052d07f86), as recorded in the CSV's comment header.

### Migration additions: native tables and plots, 22 September 2026

Lesson 3 now shows `returnLatentPresences()` and `plotLatentPresences()` using the same complete default PCR fit. The compact `latent-presence-lesson.rds` retains all 36,000 native records; Lesson 3 selects OTU_1 at two declared sites, joins by species/site/sample/primer/PCR and displays known states beside conditional probabilities. A separate table compares generating probabilities with fitted ecological, collection and PCR probabilities. HTML uses the native coloured `gt` table; GitHub Markdown uses a readable ordinary table. Colours distinguish groups of rows, not confidence or correct classification.

The native plotting appendix uses `native-plots.rds` and nine `native-plot-*.png` files. It demonstrates occupancy and collection coefficients, baseline probabilities, both primers' true/false-positive rates, residual-correlation uncertainty and cumulative detection counts with PCR or field replication on the horizontal axis. Every ecological plot adds matching truth. The cumulative plots compare survey-outcome intervals with exact generating count quantiles, not with intervals around an expected count. Raw ordination axes remain deferred until score/loading alignment is explicitly validated.

From the repository root, with the archive's library (`library`, current main; the fit-time build is `library-fbe3ed6`) first in `R_LIBS`:

```sh
Rscript dev/simstudy/vignette-lesson/summarise_latent_tables.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/verify_latent_tables.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/export_native_plots.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/verify_native_plots.R /path/to/full-fits
```

No new fitting occurs in these four commands. The latent-table check independently matches exported values to the full fit and original reads, shuffles rows before truth joins and rejects duplicate identity keys. The plotting exporter evaluates the exact displayed plotting bodies in `native-plot-examples.Rmd`; Lesson 3 contains the same bodies, not a runtime child-document dependency on the excluded `dev/` directory. The verifier checks that equivalence, plotted intervals and axes, all 45 correlation uncertainty markers and truth labels, and all 24 survey-count distributions by exhaustive enumeration of the ten binary species outcomes. Exported source, fit and image hashes guard against mixing versions.

To change a native example, edit its canonical snippet in `dev/simstudy/vignette-lesson/native-plot-examples.Rmd`, copy the revised section into Lesson 3, then re-export, verify and render. The source package includes the compact bundles and PNGs; the much larger full fit stays in the external archive. Lesson 1 also has a concise fitting reference, including model inference, missing observations, current latent-trait arguments and the actual scope of retained latent draws.

The separate [site-partitioning audit](../site-variation-partitioning.md) records how sjSDM attributes fit to components for individual sites, where its implementation needs correction, and a proposed occJSDM approach. Its lightweight algebra script runs without fitting models or importing sjSDM's Python backend:

```sh
Rscript dev/simstudy/site-variation-partitioning-audit.R /path/to/s-jSDM
```

This is an implementation assessment. No site-partitioning public function or change to occJSDM's current partition calculation is included.

### Unequal field replication, 22 September 2026

Lesson 0 now removes one whole sample at each of three distinct sites, chosen with seed 3947 before fitting: Site 12/Sample 36, Site 31/Sample 91 and Site 52/Sample 154 (Samples 24, 61 and 103 in the two-sample survey). One `(Site, Sample)` key selects matching rows in both metadata and observations. Exactly 36 PCR rows are removed; 297 samples and 3,564 PCR rows remain at all 100 sites, with 97 sites keeping three samples and three keeping two. All ten species, traits, retained values and the complete original latent truth remain unchanged. Each retained sample still has two primers and six PCRs per primer. Missing samples are not zero-filled.

Lesson 1 reconstructs the input from those saved keys, displays the matching fitting call and compares the reduced-survey estimates with truth and the original default fit. It includes all-site and selected-site figures, signed/absolute errors and diagnostics. These are a single preparation/fitting demonstration, not a replicated study of the consequences of sample loss.

The archive is `unbalanced/` inside `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples`; it holds complete `input.rds`, `unbalanced-fit.rds`, the summary and the removal keys. The two-sample build's archive was `work/unbalanced-lesson-20260922` in an earlier task's external workspace. Reproduce from the repository root with a matching occJSDM library and a fresh directory:

```sh
mkdir -p /path/to/new-unbalanced-archive
Rscript dev/simstudy/vignette-lesson/unbalanced-build.R /path/to/new-unbalanced-archive prepare
Rscript dev/simstudy/vignette-lesson/unbalanced-build.R /path/to/new-unbalanced-archive fit
Rscript dev/simstudy/vignette-lesson/unbalanced-build.R /path/to/new-unbalanced-archive summarise
Rscript dev/simstudy/vignette-lesson/unbalanced-verify.R /path/to/new-unbalanced-archive
```

The fit uses seed 20260924, default priors, two hidden site factors, one latent trait dimension and four chains with 3,000 burn-in and 6,000 retained iterations per chain. It ran once in about 76 seconds. The fit raised no R warning conditions, and the public diagnostic table has no flag: the maximum classical Rhat is 1.004197 and the minimum ESS 823.7 (the two-sample fit had one flag, OTU_3's first environmental slope at Rhat 1.015488). The 1,000 reconstructed occupancy probabilities have maximum rank-normalized Rhat 1.008026 and minimum ESS 1,220.854. Full-community occupancy MAE is 15.4039 percentage points, against 15.4352 for the default fit to the complete survey. These diagnostics and ecological errors answer different questions; neither result was hidden by selecting another seed or rerunning the fit.

The package receives only the compact `unbalanced-lesson.rds` and rendered teaching figures. The two `unbalanced-lesson-*.Rmd` files under `dev/` are source snippets mirrored in Lessons 0 and 1, not child documents required for knitting an installed vignette. Keep displayed code and these snippets in step when editing the example.

Validation for this migration batch, 22 September 2026: the independent native-table, native-plot, unbalanced-fit and site-partition algebra checks passed. All five teaching documents rendered to HTML and Markdown; the shared link regression test passed. New static figures were inspected, and the native HTML tables were checked in the generated markup. The source-package build succeeded with vignette building disabled, included all new compact data and images, and excluded the archived original, lesson plan and `dev/`. Lessons 1 and 3 then rendered successfully from the unpacked source package, without the excluded development snippets or external full fits. No package-wide test run was needed for these teaching additions; production R/C++ code is unchanged. The original archived vignette remains byte-identical to revision 8654ff1.

## Remaining non-spatial migration: ordination, native plots and prediction

The follow-up on `codex/lesson-migration-completion` adds the remaining valid non-spatial plotting examples, the original bibliography, and a matched independent-site prediction comparison to Lesson 3 (moved to Lesson 4 in Phase B of the lesson rewrite). It changes no production R/C++ code. The first migration batch was already merged as e0e813d.

### Ordination and remaining plotting functions

Run from the repository root with the archive's library (`library`, current main; the fit-time build is `library-fbe3ed6`) first in `R_LIBS`:

```sh
Rscript dev/simstudy/vignette-lesson/ordination-export.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/ordination-verify.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/remaining-plots-export.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/remaining-plots-verify.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/covariate-effect-export.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/covariate-effect-verify.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/native-traits-export.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/native-traits-verify.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/remaining-plots-diagnose-covariate.R /path/to/full-fits
```

`remaining-plots-diagnose-covariate.R` reproduces the pre-PR #13 covariate-response bug against the old archive and package only; it cannot pass on code after PR #13, which removed `plotCovariateEffect_base()` and corrected the calculation, so it is not run in a rebuild from current main.

The canonical snippets `ordination-examples.Rmd` and `remaining-plots-examples.Rmd` are copied into Lesson 3, not sourced as child documents. Their verifiers check that the displayed bodies exactly match the exported calculations. Keep them synchronized when editing. The compact RDS files and ten native PNGs are included in the source package. Neither rendering nor these plot exports reruns MCMC.

Ordination uses draw-wise orthogonal Procrustes alignment against the known generating loadings, transforming scores and loadings jointly. It is a simulation-only aid, not a claim that an axis has been identified independently of truth. All 24,000 draws are retained. The verifier uses an independent analytic two-dimensional rotation/reflection calculation and checks unchanged score/loading products and Gram matrices, native marginal quantiles, identities, plotted coordinates, circle widths, arrow scaling and image hashes. The lesson explains that native circles summarize two marginal interval widths, not joint credible regions. All 100 sites appear in the biplot; ten sites chosen by original order have separate common-axis panels, and all ten species have loading panels.

The `native-traits-examples.Rmd` snippet is likewise mirrored in Lesson 3. It restores both `plotTraitsCoefficients()` examples and matches generating coefficients to the fitted trait scale by multiplying each raw generating effect by the corresponding sample trait standard deviation. Its exporter and verifier check all 24,000 G draws, four native intervals and matching plotted truths.

The remaining plots cover both environmental gradients and separate field false-positive, laboratory false-positive and laboratory true-positive rates. Their direct all-draw verification checks 800 gradient summaries and 50 rate intervals, correct species/primer identities and plotted dodge coordinates. Probability curves hold other standardized predictors at their observed medians and site factors at zero. Laboratory truth is adjusted for the observed read threshold. See [the detailed plotting notes](remaining-plots-README.md).

The `covariate-effect-examples.Rmd` snippet is likewise mirrored in Lesson 3, in its own exporter and verifier pair so that no earlier bundle changes. It shows `plotCovariateEffect()` for environmental gradient 1 in original units, with the generating curve on the function's own 200-point grid (gradient 2 at its standardized median, site factors at zero). The verifier recomputes all 2000 medians and 95% limits from the 24,000 draws and the truth from the raw survey table. It also checks the plotted ribbon, line and truth layers by panel, recomputes the recorded coverage checks, and requires the lesson's displayed chunks to equal the snippet's.

The migration reproduced a separate defect in `returnCovariateEffect()` / `plotCovariateEffect()`: the numeric helper adds a log-odds intercept to an already inverse-logit-transformed partial effect, and its supposed raw grid is already standardized and then standardized again. The diagnostic demonstrated out-of-range median probabilities on the two-sample default fit. [PR #13](https://github.com/AlexDiana/occJSDM/pull/13) has since merged: it corrected `returnCovariateEffect()` and `plotCovariateEffect()` and removed `plotCovariateEffect_base()`, so the diagnostic reproduces the bug only against the old archive and package. Lesson 3 keeps the verified `plotOccupancyGradient()` examples for its standardized response curves and, since 4 October 2026, also shows the corrected `plotCovariateEffect()` in original units (see above). Two md5-stamped dev snippets kept stale sentences after PR #13 and the three-sample rebuild, because editing a snippet fails its verifier's `snippet_md5` check until its bundle is re-exported. `native-plot-examples.Rmd` (line 274) said the fitted data contained "two field samples per site", and `remaining-plots-examples.Rmd` (line 73) said `plotCovariateEffect()` "cannot currently be restored". Both were corrected on 4 October 2026, and the native-plot, remaining-plots and covariate-effect bundles were re-exported; the last is included because it records the remaining-plots bundle's md5. No PNG changed.

### Matched models and genuinely new sites

The full-fit archive is `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples`. New inputs, one additional full fit and full prediction exports are in its `prediction/` directory. Both remain outside the package. The two-sample build used `work/vignette-lesson-20260919` and `work/prediction-lesson-20260922` in an earlier task's work directory. The tracked compact bundle is `vignettes/teaching-data/prediction-lesson.rds`.

Create a separate output directory, then run these commands using the original archive's library. `prepare` and `fit` refuse to overwrite their saved inputs/fits:

```sh
mkdir -p /path/to/new-site-check
Rscript dev/simstudy/vignette-lesson/prediction-build.R /path/to/full-fits /path/to/new-site-check prepare
Rscript dev/simstudy/vignette-lesson/prediction-build.R /path/to/full-fits /path/to/new-site-check fit
Rscript dev/simstudy/vignette-lesson/prediction-export.R /path/to/full-fits /path/to/new-site-check
Rscript dev/simstudy/vignette-lesson/prediction-verify.R /path/to/full-fits /path/to/new-site-check
# Optional full export reproduction and score Monte Carlo uncertainty:
Rscript dev/simstudy/vignette-lesson/prediction-verify.R /path/to/full-fits /path/to/new-site-check --score-mcse
```

Preparation seed 20260925 creates 300 independent sites numbered 101-400 from the original raw environmental distribution, retaining training standardization constants and the original ecological parameters. New hidden conditions and Bernoulli occupancy states are generated without changing the training survey. Thirteen new sites have a covariate beyond the observed training range; none is discarded. Sites 101-110 are selected for the native-call illustration before inspecting predictions. The simulator's new hidden conditions and occupancy states are used only for evaluation, never for fitting either model.

The additional fit uses seed 20260926 and changes `n_factors` from two to one. It retains the original PCR observations, two environmental predictors, collection covariate, traits with one latent trait, all priors, threshold one, four chains, 3,000 burn-in and 6,000 retained iterations per chain with no thinning. Exact training observations and fitted design matrices are checked against the original two-factor fit (seed 20260921). There is one fit per specification, not replicated communities or independent repeated MCMC runs.

The public `predictNewSites()` illustration uses seed 20260927, ten predeclared new sites and all retained draws. Its output slices are lower quantile, median and upper quantile. Because it draws unknown new-site conditions as well as parameters, the lesson compares its intervals with the generating probabilities conditional on the sites' actual simulated hidden conditions. The current implementation draws factors separately within species; its per-species marginal outputs must not be presented as coherent joint community realizations.

For point prediction, `prediction-math.R` deterministically integrates each draw's logistic probability over Gaussian hidden conditions, then averages over all 24,000 posterior draws. It uses the full residual standard deviation `sigma_h * sqrt(sum(L^2))` for that species and draw, with raw environmental inputs transformed by the frozen training constants. This is a development helper, not a new public API. It avoids confusing the native median or the zero-factor response profile with the marginal posterior mean. Quadrature is selected by an explicit error check over locations -30 to 30 and spreads up to the largest stored spread: 61 nodes for the two-factor fit and 31 for the one-factor fit, with grid discrepancies of 3.0e-7 and 4.2e-6 respectively. Generating marginal truth is computed separately from known parameters.

The independent-site scores are based on exactly the same 3,000 species-site outcomes for both models. Mean absolute error against marginal truth is 8.9543 percentage points for one factor and 8.9101 for two; signed errors are +0.5250 and +0.3606 points. Brier scores are 0.1871684 and 0.1868378, respectively, and negative log scores 0.5546200 and 0.5538165. The one-minus-two paired Brier difference is +0.0003306, slightly favouring two factors (in the two-sample survey it was -0.0001365, slightly favouring one), with standard error 0.0000664 computed across 300 site averages. The species within a site are not counted as independent replicates. This site-only uncertainty excludes MCMC error, repeated-training-sample uncertainty and variation among communities. It is not evidence for a generally superior model or a recovered true factor count.

All 200 public parameter diagnostic rows are below Rhat 1.01 and above ESS 400, with no missing values, and the new fit emitted no warnings. Across species' mean predicted probabilities, maximum Rhat is 1.006579, minimum ESS about 1,218, and largest Monte Carlo standard error 0.002591 (0.259 percentage points). These are computation checks, not proof that the tiny score difference is resolved or that ecological prediction errors are acceptable.

The lesson also displays the archived legacy stored scores, without ranking models by them. The exporter and verifier read `results_output$WAIC` explicitly; `extractWAIC(..., type = "legacy")` retrieves the same scalar. `runOccJSDM()` combines likelihood terms for sampled occupancy/collection states with observed PCR terms. This is not a new-site observed-data criterion integrating the unobserved states; fitting both models to identical observations does not resolve that target mismatch. PR #14 adds an observed-data site criterion, approved by Alex on 6 October 2026; reliability validation on the current three-sample fits remains separate work. Existing one-factor and two-factor new-site scores teach a valid same-data predictive comparison while leaving that limitation explicit.

The independent prediction verifier also estimates a first-order Monte Carlo standard error for each score by propagating the posterior-draw variation through the score gradient. With independent fits, the two fit-specific numerical variances add for their difference. The Brier difference has estimated MCSE 0.0002601 (its magnitude is 1.27 MCSE); the negative-log-score difference has MCSE 0.0007551 (about 1.06 MCSE). These numerical uncertainties are distinct from the new-site sampling SE. No longer MCMC run or seed selection was performed to resolve a model ranking.

A preliminary wrapper probe also found that reshaping an empty factor-loading array can fail in `predictNewSites()` before reaching C++, even with `useBiotic = FALSE`. This was tested by replacing the archived fit's loading array with a zero-factor array, not by fitting a complete zero-factor model. Reproduce the wrapper probe with `Rscript dev/simstudy/vignette-lesson/prediction-zero-factor-probe.R /path/to/full-fits`. A genuine zero-factor fit/regression test remains a separate follow-up; the teaching comparison uses one and two factors and does not claim to validate the zero-factor path.

Validation of the completed follow-up, 22 September 2026: ordination, remaining rate/gradient plots, native trait plots and the full prediction verifier passed against the archived matching library. The full prediction pass independently reproduced all 6,000 exported posterior means/intervals and corresponding probability diagnostics using every retained draw, checked native call arguments and seeds without refitting, and reproduced the first-order score MCSE. Independent adaptive integration agreed with selected posterior means within 1.3e-11. The response-curve defect and separate empty-array wrapper probe reproduced as documented.

All five lessons rendered to HTML and GitHub Markdown. Lesson 3 was rendered again after the final trait addition; all twelve new figures were visually inspected. Shared link regression checks passed, new image/data paths exist, and there are no duplicate chunk names. The final source-package build succeeded with vignette building disabled, included all sixteen new lesson data/image assets, excluded `dev/`, `LESSON-PLAN.md` and the archived original, and rendered Lesson 3 successfully after unpacking without external full fits. Production R/C++ files and the archived original remain unchanged. Before merging, the source test suite was run: 641 expectations passed, with no failures, errors or test warnings; the opt-in coverage study was skipped. The ordination, gradient/rate, trait and standard prediction archive verifiers, Lesson 3 HTML/Markdown renders and shared-link checks were rerun successfully. The fresh source-package build included all sixteen new teaching assets, excluded the development scripts, lesson plan and archived original, and rendered Lesson 3 after unpacking. No new model fitting or release-bias study was run during this merge check.

### Lesson 3: package plots for gradient 2, the perfect-observation fit and the ordinary biplot, added 3 October 2026

Phase B of the lesson rewrite makes the package's own plot with truth the main figure of each Lesson 3 section. Three exporters were extended and re-run on `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-3samples` without refitting: `export_native_plots.R` adds `plotOccupancyCovariates()` for gradient 2 and for both gradients of the perfect-observation fit; `native-traits-export.R` adds `plotTraitsCoefficients()` for both gradients of the perfect-observation fit; `ordination-export.R` saves the ordinary, unaligned biplot without a truth layer. The matching verifiers re-derive every new plotted interval or point from all 24,000 draws of the relevant fit, re-derive the truth independently, and still require the lesson's displayed code to equal the exported code. A `plotVariancePartitioning()` export with the true shares was added and then withdrawn the same day: the owner will include variation partitioning once a working spatial model exists, so `remaining-plots-examples.Rmd`, `remaining-plots-export.R` and `remaining-plots-verify.R` are back to their earlier content and `remaining-plots-data.rds` was re-exported from them. The runs are logged as steps 40 to 57 in that archive's `BUILD-LOG.md`.


### PR #14 integration and historical archive verification (7 October 2026)

The new namespace, C++ and R scoring code is outside `source-check.R`'s existing plot-only compatibility allowlist. In addition, `prediction-export.R` itself is recorded in the compact bundle's `source_md5`, so its explicit legacy-score extraction changes the recorded source identity even though the numerical value is unchanged. Run existing full-archive exporters and verifiers with their matching historical source and library. They are not a compatibility check against the merged package. Do not change recorded hashes or broaden the allowlist to bypass these gates. The 7 October rebuild above supplies a new non-spatial archive under the merged source. The earlier integration preserved all historical bundles; those archives still require their matching historical source. Rendering uses those compact bundles without executing archive verification or the new scoring example.

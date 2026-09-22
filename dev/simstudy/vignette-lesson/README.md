# Build the truth-based teaching lessons

The teaching sequence now starts with `vignettes/occJSDM-lesson-0.Rmd` (optional simulation and data orientation), followed by `vignettes/occJSDM-lesson-1.Rmd` (fitting, truth comparisons and detection examples). Both show their teaching code and use the same original simulation and baseline fits. The unequal-replication extension described below adds one fit to a reduced version of that survey. Their new maps display the existing sampled coordinates; the environmental values and hidden site factors are independent of geography in this dataset. No interpolation or spatial fit is implied.

`vignettes/occJSDM-lesson-2.Rmd` is explicitly a future outline. After PR #8 is reviewed, it will introduce smooth environmental gradients, residual spatial structure, a separate dispersal simulation with contrasting species, and held-out prediction. The currently reviewed spatial model has a shared range, not mechanistic or species-specific dispersal parameters. That distinction constrains the later simulation and interpretation.

The readable lesson is `vignettes/occJSDM-lesson-1.Rmd`. Its figures and tables come from `vignettes/teaching-data/nonspatial-lesson.rds`, a compact bundle containing the complete simulation, generating settings, identity mappings, selected cases, posterior summaries, diagnostics and full-fit hashes. Full MCMC objects remain in the build directory; they are not shipped with the package. The existing `sampledata` and `sampleresults` objects are unchanged.

The experiment uses main revision b53048a's non-spatial model, including the reviewed collection alignment, RNG and factor-reparameterisation fixes and the read-threshold correction. It does not test the spatial changes in PR #8. The documentation branch changes no code in `R/` or `src/`.

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

The initial alternative-prior fit is retained, but the published example uses its longer fit because an individual field-contamination parameter had Rhat 1.025 and ESS about 230 initially. `summarise_lesson.R` uses a saved `*-long-fit.rds` when one is present. The table in the lesson reports the selected fits' actual diagnostics; no scientific error target is imposed.

`prepare` refuses to overwrite an input bundle. The fitting step refuses to overwrite a fit, and checks the source hashes against the input. The verification step checks generator hashes, fit hashes and input identity. To change the simulation or fit source, use a new build directory and regenerate the evidence rather than mix old fits with new truth. Selected cases are saved to `cases-selected-before-fitting.csv` before any fit is run. MCMC seeds are independent of the simulation seed and are fixed per fit.

The source categories describe this simulator: a positive with `w=0` is a laboratory false positive; a positive with `z=0,w=1` is a field-stage false positive; a positive with `z=1,w=1` is a true detection. Missing observations are separate. The first sample in each declared category is selected in lexicographic species-name, numeric-site and numeric-sample order. Consequently `OTU_10` sorts before `OTU_2`; this is intentional and does not use fitted results. Category totals and the all-positive-sample table prevent the selected cases from being presented as a general error-rate estimate.

The prior stress test changes both `q` and `theta0` from Beta(1,20) to Beta(1,4), leaving `p` unchanged. It cannot distinguish which of the two altered priors is responsible for a change. It illustrates sensitivity to the joint low-contamination assumption and is not a recommendation to change the defaults.

## Render without rerunning MCMC

The small RDS is tracked so a fresh checkout can render immediately. These commands run from the repository root:

```sh
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-1.Rmd", output_format="rmarkdown::html_vignette")'
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-1.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE))'
```

Use the same commands with `occJSDM-lesson-0.Rmd` and `occJSDM-lesson-2.Rmd` to render the optional introduction and planned spatial outline. Render each in a fresh R environment when checking standalone execution. Teaching chunks are visible by default; only document-formatting setup is hidden. Simulation and MCMC examples marked `eval=FALSE` are displayed but do not run while knitting. The Lesson 0 simulator call reproduces the existing settings without calling the build helper; the helper remains unchanged to preserve its recorded provenance.

The GitHub Markdown and its PNGs are tracked as readable review artifacts. The HTML is generated locally. The R Markdown source remains canonical. Neither render loads full MCMC fits or starts fitting. To reproduce the original numerical results exactly, use the recorded source revision and R/package versions; different platforms can still introduce small numerical differences.

### Editing links in RStudio

Write links between lessons as ordinary Markdown, for example `[Lesson 1](occJSDM-lesson-1.md)` or `[Diagnostics](occJSDM-lesson-3.md#check-computation-as-well-as-ecological-recovery)`. Do not put inline R expressions inside a link destination: RStudio's Visual editor can encode those expressions as URL text, preventing knitr from evaluating them.

Each of the five teaching documents sources `vignettes/lesson-links.R` in its hidden setup chunk. This shared document hook changes recognised lesson destinations from `.md` to `.html` only in HTML output. It preserves section anchors, external links, other documents and code examples. It does not edit the `.Rmd` source. When adding a lesson, add its filename stem to the helper's explicit list and source the helper in that lesson's setup. The helper is tracked and included in the source package; `ORIG_occJSDM.Rmd` remains unchanged and excluded.

Run the focused regression check from the repository root:

```sh
Rscript dev/simstudy/vignette-lesson/test_lesson_links.R
```

The check renders a small fixture to HTML, Markdown and HTML again in one R session. It verifies the destinations and anchors, preserves literal examples and unrelated links, and rejects dynamic or encoded lesson destinations in the five lesson sources. No MCMC is run. When changing the link convention, also edit and save disposable copies in RStudio Visual mode and render those copies to both formats; a normal command-line render alone does not exercise the editor's Markdown rewriting.

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

`vignettes/occJSDM-lesson-3.Rmd` extends the agreed design to environmental and trait effects, response profiles, baseline probabilities, residual correlations, combined ordination contributions, variation partitioning, collection effects and expected detection with additional sampling. `vignettes/occJSDM.Rmd` is now the short quickstart/lesson guide; the old unmatched output tour is replaced by Lesson 3 and the existing Lesson 1. The original example data and fits remain unchanged.

The new lesson reuses the complete `perfect-fit.rds` and `default-fit.rds` from the same archive. Their model source hashes match current main revision 8654ff1: the changes since the recorded b53048a source are documentation changes. No new MCMC or model changes are needed. Run the following from the repository root, with a matching occJSDM library:

```sh
Rscript dev/simstudy/vignette-lesson/summarise_outputs.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/verify_outputs.R /path/to/full-fits
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-3.Rmd", output_format="rmarkdown::html_vignette")'
Rscript -e 'rmarkdown::render("vignettes/occJSDM-lesson-3.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE))'
```

The exporter writes `vignettes/teaching-data/output-lesson.rds`, including the original fit manifests and hashes of its source, the original compact bundle and the legacy saved example used for the historical interval-count check. It checks actual training observations and reconstructs the full true predictor to establish the environmental scale. True trait effects are converted using trait standard deviations; the simulator standardizes its environmental predictors but not its traits. Changing exporter code requires re-exporting and re-verifying the compact bundle.

Checks independently reconstruct raw coefficient summaries, response-profile endpoints, selected posterior products of scores and loadings, correlation intervals and the non-spatial partition truth. They also verify the trait-decomposition algebra, collection standardization and detection-effort expectations. The effort check uses an independent sum over the possible number of collected field samples. Fits remain in the external archive, not in the package.

Validation on 21 September 2026: the independent numerical checks passed for both fits, including direct chain diagnostics for the environmental and trait coefficients. Lesson 3 and the quickstart rendered to HTML and GitHub Markdown. All ten figures were visually inspected; local links, embedded teaching CSS and visible code blocks were checked. Scientific review found no significant issues. The existing lessons, simulation bundle and model source are unchanged; this documentation change did not rerun the package-wide test suite or MCMC.

Later integration validation on 21 September 2026: Lesson 3 and the Lesson N design were combined for pull request review; the Lesson 2 outline was already present on main. The package test suite passed 641 expectations with zero test failures or warnings; the long coverage study remained explicitly skipped. All five teaching documents rendered to HTML and Markdown, and their local links were checked. Encoded inline-R lesson links introduced by RStudio formatting were repaired, and Lesson 1 now directs readers to Lesson 3. A source-package archive retained the lesson sources while excluding `vignettes/LESSON-PLAN.md` and the development plans. Production model code and existing simulation results were unchanged.

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

### Migration additions: native tables and plots, 22 September 2026

Lesson 3 now shows `returnLatentPresences()` and `plotLatentPresences()` using the same complete default PCR fit. The compact `latent-presence-lesson.rds` retains all 24,000 native records; Lesson 3 selects OTU_1 at two declared sites, joins by species/site/sample/primer/PCR and displays known states beside conditional probabilities. A separate table compares generating probabilities with fitted ecological, collection and PCR probabilities. HTML uses the native coloured `gt` table; GitHub Markdown uses a readable ordinary table. Colours distinguish groups of rows, not confidence or correct classification.

The native plotting appendix uses `native-plots.rds` and nine `native-plot-*.png` files. It demonstrates occupancy and collection coefficients, baseline probabilities, both primers' true/false-positive rates, residual-correlation uncertainty and cumulative detection counts with PCR or field replication on the horizontal axis. Every ecological plot adds matching truth. The cumulative plots compare survey-outcome intervals with exact generating count quantiles, not with intervals around an expected count. Raw ordination axes remain deferred until score/loading alignment is explicitly validated.

From the repository root, with the archive's matching library first in `R_LIBS`:

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

Lesson 0 now removes one whole sample at each of three distinct sites, chosen with seed 3947 before fitting: Site 12/Sample 24, Site 31/Sample 61 and Site 52/Sample 103. One `(Site, Sample)` key selects matching rows in both metadata and observations. Exactly 36 PCR rows are removed; 197 samples and 2,364 PCR rows remain at all 100 sites. All ten species, traits, retained values and the complete original latent truth remain unchanged. Each retained sample still has two primers and six PCRs per primer. Missing samples are not zero-filled.

Lesson 1 reconstructs the input from those saved keys, displays the matching fitting call and compares the reduced-survey estimates with truth and the original default fit. It includes all-site and selected-site figures, signed/absolute errors and diagnostics. These are a single preparation/fitting demonstration, not a replicated study of the consequences of sample loss.

The new archive is `work/unbalanced-lesson-20260922` in this task's external workspace; it holds complete `input.rds`, `unbalanced-fit.rds`, the summary and logs. The original archive was not changed. Reproduce from the repository root with a matching occJSDM library and a fresh directory:

```sh
mkdir -p /path/to/new-unbalanced-archive
Rscript dev/simstudy/vignette-lesson/unbalanced-build.R /path/to/new-unbalanced-archive prepare
Rscript dev/simstudy/vignette-lesson/unbalanced-build.R /path/to/new-unbalanced-archive fit
Rscript dev/simstudy/vignette-lesson/unbalanced-build.R /path/to/new-unbalanced-archive summarise
Rscript dev/simstudy/vignette-lesson/unbalanced-verify.R /path/to/new-unbalanced-archive
```

The fit uses seed 20260924, default priors, two hidden site factors, one latent trait dimension and four chains with 3,000 burn-in and 6,000 retained iterations per chain. It ran once in about 67 seconds. The fit raised no R warning conditions, but the public diagnostic table has one flag: OTU_3's first environmental slope has classical Rhat 1.015488 and ESS 1,601. The 1,000 reconstructed occupancy probabilities have maximum rank-normalized Rhat 1.007584 and minimum ESS 940.211. Full-community occupancy MAE is 17.0439 percentage points. These diagnostics and ecological errors answer different questions; neither result was hidden by selecting another seed or rerunning the fit.

The package receives only the compact `unbalanced-lesson.rds` and rendered teaching figures. The two `unbalanced-lesson-*.Rmd` files under `dev/` are source snippets mirrored in Lessons 0 and 1, not child documents required for knitting an installed vignette. Keep displayed code and these snippets in step when editing the example.

Validation for this migration batch, 22 September 2026: the independent native-table, native-plot, unbalanced-fit and site-partition algebra checks passed. All five teaching documents rendered to HTML and Markdown; the shared link regression test passed. New static figures were inspected, and the native HTML tables were checked in the generated markup. The source-package build succeeded with vignette building disabled, included all new compact data and images, and excluded the archived original, lesson plan and `dev/`. Lessons 1 and 3 then rendered successfully from the unpacked source package, without the excluded development snippets or external full fits. No package-wide test run was needed for these teaching additions; production R/C++ code is unchanged. The original archived vignette remains byte-identical to revision 8654ff1.

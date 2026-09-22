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

## Remaining non-spatial migration: ordination, native plots and prediction

The follow-up on `codex/lesson-migration-completion` adds the remaining valid non-spatial plotting examples, the original bibliography, and a matched independent-site prediction comparison to Lesson 3. It changes no production R/C++ code. The first migration batch was already merged as e0e813d.

### Ordination and remaining plotting functions

Run from the repository root with the original archive's matching installed library first in `R_LIBS`:

```sh
Rscript dev/simstudy/vignette-lesson/ordination-export.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/ordination-verify.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/remaining-plots-export.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/remaining-plots-verify.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/native-traits-export.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/native-traits-verify.R /path/to/full-fits
Rscript dev/simstudy/vignette-lesson/remaining-plots-diagnose-covariate.R /path/to/full-fits
```

The canonical snippets `ordination-examples.Rmd` and `remaining-plots-examples.Rmd` are copied into Lesson 3, not sourced as child documents. Their verifiers check that the displayed bodies exactly match the exported calculations. Keep them synchronized when editing. The compact RDS files and ten native PNGs are included in the source package. Neither rendering nor these plot exports reruns MCMC.

Ordination uses draw-wise orthogonal Procrustes alignment against the known generating loadings, transforming scores and loadings jointly. It is a simulation-only aid, not a claim that an axis has been identified independently of truth. All 24,000 draws are retained. The verifier uses an independent analytic two-dimensional rotation/reflection calculation and checks unchanged score/loading products and Gram matrices, native marginal quantiles, identities, plotted coordinates, circle widths, arrow scaling and image hashes. The lesson explains that native circles summarize two marginal interval widths, not joint credible regions. All 100 sites appear in the biplot; ten sites chosen by original order have separate common-axis panels, and all ten species have loading panels.

The `native-traits-examples.Rmd` snippet is likewise mirrored in Lesson 3. It restores both `plotTraitsCoefficients()` examples and matches generating coefficients to the fitted trait scale by multiplying each raw generating effect by the corresponding sample trait standard deviation. Its exporter and verifier check all 24,000 G draws, four native intervals and matching plotted truths.

The remaining plots cover both environmental gradients and separate field false-positive, laboratory false-positive and laboratory true-positive rates. Their direct all-draw verification checks 800 gradient summaries and 50 rate intervals, correct species/primer identities and plotted dodge coordinates. Probability curves hold other standardized predictors at their observed medians and site factors at zero. Laboratory truth is adjusted for the observed read threshold. See [the detailed plotting notes](remaining-plots-README.md).

The migration reproduced a separate defect in `returnCovariateEffect()` / `plotCovariateEffect()`: the numeric helper adds a log-odds intercept to an already inverse-logit-transformed partial effect, and its supposed raw grid is already standardized and then standardized again. The diagnostic demonstrates out-of-range median probabilities on the unchanged default fit. TODO now records the correction implemented in the separate [PR #13](https://github.com/AlexDiana/occJSDM/pull/13), awaiting Alex's review; no production fix or replacement fit is included in this lesson branch. The lesson uses the verified `plotOccupancyGradient()` alternative.

### Matched models and genuinely new sites

The original full-fit archive is `work/vignette-lesson-20260919` under the current task's work directory. New inputs, one additional full fit and full prediction exports are in `work/prediction-lesson-20260922`. Both remain outside the package. The tracked compact bundle is `vignettes/teaching-data/prediction-lesson.rds`.

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

For point prediction, `prediction-math.R` deterministically integrates each draw's logistic probability over Gaussian hidden conditions, then averages over all 24,000 posterior draws. It uses the full residual standard deviation `sigma_h * sqrt(sum(L^2))` for that species and draw, with raw environmental inputs transformed by the frozen training constants. This is a development helper, not a new public API. It avoids confusing the native median or the zero-factor response profile with the marginal posterior mean. Quadrature is selected by an explicit error check over locations -30 to 30 and spreads up to the largest stored spread: 61 nodes for the two-factor fit and 31 for the one-factor fit, with grid discrepancies below 2.4e-7. Generating marginal truth is computed separately from known parameters.

The independent-site scores are based on exactly the same 3,000 species-site outcomes for both models. Mean absolute error against marginal truth is 9.9349 percentage points for one factor and 10.0174 for two; signed errors are -3.8053 and -3.8321 points. Brier scores are 0.1846826 and 0.1848191, respectively. The one-minus-two paired Brier difference is -0.0001365, with standard error 0.0000670 computed across 300 site averages. The species within a site are not counted as independent replicates. This site-only uncertainty excludes MCMC error, repeated-training-sample uncertainty and variation among communities. It is not evidence for a generally superior model or a recovered true factor count.

All 200 public parameter diagnostic rows are below Rhat 1.01 and above ESS 400, with no missing values, and the new fit emitted no warnings. Across species' mean predicted probabilities, maximum Rhat is 1.008944, minimum ESS about 746, and largest Monte Carlo standard error 0.003886 (0.389 percentage points). These are computation checks, not proof that the tiny score difference is resolved or that ecological prediction errors are acceptable.

The lesson also extracts actual current WAIC values, without ranking models by them. `runOccJSDM()` combines likelihood terms for sampled occupancy/collection states with observed PCR terms. This is not a new-site observed-data criterion integrating the unobserved states; fitting both models to identical observations does not resolve that target mismatch. A validated observed-data criterion or site-level cross-validation remains separate work. Existing one-factor and two-factor new-site scores teach a valid same-data predictive comparison while leaving that limitation explicit.

The independent prediction verifier also estimates a first-order Monte Carlo standard error for each score by propagating the posterior-draw variation through the score gradient. With independent fits, the two fit-specific numerical variances add for their difference. The Brier difference has estimated MCSE 0.0002539 (its magnitude is only 0.54 MCSE); the negative-log-score difference has MCSE 0.0006905 (about 1.02 MCSE). These numerical uncertainties are distinct from the new-site sampling SE. No longer MCMC run or seed selection was performed to resolve a model ranking.

A preliminary wrapper probe also found that reshaping an empty factor-loading array can fail in `predictNewSites()` before reaching C++, even with `useBiotic = FALSE`. This was tested by replacing the archived fit's loading array with a zero-factor array, not by fitting a complete zero-factor model. Reproduce the wrapper probe with `Rscript dev/simstudy/vignette-lesson/prediction-zero-factor-probe.R /path/to/full-fits`. A genuine zero-factor fit/regression test remains a separate follow-up; the teaching comparison uses one and two factors and does not claim to validate the zero-factor path.

Validation of the completed follow-up, 22 September 2026: ordination, remaining rate/gradient plots, native trait plots and the full prediction verifier passed against the archived matching library. The full prediction pass independently reproduced all 6,000 exported posterior means/intervals and corresponding probability diagnostics using every retained draw, checked native call arguments and seeds without refitting, and reproduced the first-order score MCSE. Independent adaptive integration agreed with selected posterior means within 1.3e-11. The response-curve defect and separate empty-array wrapper probe reproduced as documented.

All five lessons rendered to HTML and GitHub Markdown. Lesson 3 was rendered again after the final trait addition; all twelve new figures were visually inspected. Shared link regression checks passed, new image/data paths exist, and there are no duplicate chunk names. The final source-package build succeeded with vignette building disabled, included all sixteen new lesson data/image assets, excluded `dev/`, `LESSON-PLAN.md` and the archived original, and rendered Lesson 3 successfully after unpacking without external full fits. Production R/C++ files and the archived original remain unchanged. Before merging, the source test suite was run: 641 expectations passed, with no failures, errors or test warnings; the opt-in coverage study was skipped. The ordination, gradient/rate, trait and standard prediction archive verifiers, Lesson 3 HTML/Markdown renders and shared-link checks were rerun successfully. The fresh source-package build included all sixteen new teaching assets, excluded the development scripts, lesson plan and archived original, and rendered Lesson 3 after unpacking. No new model fitting or release-bias study was run during this merge check.

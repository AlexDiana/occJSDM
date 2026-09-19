# Build the first truth-based lesson

The readable lesson is `vignettes/occJSDM-first-lesson.Rmd`. Its figures and tables come from `vignettes/teaching-data/nonspatial-lesson.rds`, a compact bundle containing the complete simulation, generating settings, identity mappings, selected cases, posterior summaries, diagnostics and full-fit hashes. Full MCMC objects remain in the build directory; they are not shipped with the package. The existing `sampledata` and `sampleresults` objects are unchanged.

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
Rscript -e 'rmarkdown::render("vignettes/occJSDM-first-lesson.Rmd", output_format="rmarkdown::html_vignette")'
Rscript -e 'rmarkdown::render("vignettes/occJSDM-first-lesson.Rmd", output_format=rmarkdown::github_document(html_preview=FALSE))'
```

The GitHub Markdown and its PNGs are tracked as readable review artifacts. The HTML is generated locally. The R Markdown source remains canonical. Neither render loads full MCMC fits or starts fitting. To reproduce the original numerical results exactly, use the recorded source revision and R/package versions; different platforms can still introduce small numerical differences.

## What is checked

- Helper tests distinguish source categories despite reordered observations/species, retain missingness and verify the read-rounding correction.
- The evidence check confirms that the binary fit uses the same actual z and that the two-stage fits share their simulation and source hashes.
- Occupancy summaries are reconstructed from posterior model components; the independent check uses a separate component-wise calculation for all sites and draws of OTU_1, including interval endpoints. For both two-stage fits, all 1,000 cell means are also compared with the fitter's saved psi means.
- Reported sample/site probabilities are compared with their saved fitted arrays; all band errors are recalculated.
- Numerical diagnostics are exported for reported probabilities and named parameter blocks. Neither converged sampling nor passing these checks proves unbiased ecological inference.

## Remaining tutorial work

Extend the agreed design to environmental and trait effects, residual correlations, variation partitioning, spatial examples and genuine held-out prediction. Replace the older reference walkthrough's unmatched saved examples as those lessons are finished. Paper2Agent remains a feasibility proposal and has not been installed or run.

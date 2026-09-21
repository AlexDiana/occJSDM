occJSDM: quickstart and lesson guide
================

occJSDM estimates species occurrence while allowing for imperfect field
collection, PCR detection and false positives. It can also fit a JSDM to
directly observed presence/absence data. The teaching lessons show the R
code and compare the results with the truth used to simulate their data.

## Choose a lesson

- [Lesson 0: Create and explore a simulated
  survey](occJSDM-lesson-0.md). Optional preparation: sites, field
  samples, primers, PCR replicates, environmental covariates, traits and
  maps of the known truth.
- [Lesson 1: Fit the model and compare its answers with
  truth](occJSDM-first-lesson.md). The main starting point: fit the
  model, understand occupancy error, and inspect actual examples of weak
  true detections, strong true detections and false positives.
- [Lesson 2: Spatial landscapes and dispersal](occJSDM-lesson-2.md). A
  planned lesson on smooth environmental gradients, residual spatial
  structure and a separate dispersal simulation. It is not yet a worked
  spatial validation.
- [Lesson 3: Understand the model’s outputs by comparing them with
  truth](occJSDM-lesson-3.md). Environmental responses, trait effects,
  species associations, ordination, variation partitioning and detection
  effort. It includes a function index and explains why a real simulated
  effect can remain uncertain.

You can read Lesson 3 after Lesson 1; the spatial lesson is not a
prerequisite. The earlier output tour has been replaced by Lesson 3. The
old `sampledata` and `sampleresults` objects remain available, but their
generating truth is not stored with them, so the lessons use a complete
matching simulation and fit instead.

## A minimal fitting example

From the repository’s `vignettes` directory, the following optional code
fits the same non-spatial survey used in the lessons. It is shown
without running during knitting because MCMC takes time. Lessons 1 and 3
render saved, checked results and show their matching truth.

``` r
library(occJSDM)

lesson <- readRDS("teaching-data/nonspatial-lesson.rds")
survey_data <- lesson$input$sim$data_list
known_truth <- lesson$input$sim$true_params

set.seed(20260921)

fitmodel <- runOccJSDM(
  data = survey_data,
  listParams = list(n_factors = 2, n_lattrait = 1),
  threshold = 1,
  occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
  collCovariates = "X_theta",
  spatCovariates = NULL,
  MCMCparams = list(nchain = 4, nburn = 3000, niter = 6000, nthin = 1),
  summarisedLatentPresences = TRUE
)
```

`survey_data` is what the model receives. `known_truth` is kept
separately for evaluation. `n_factors` controls the hidden site factors,
whereas `n_lattrait` controls the unmeasured species-trait structure.
These are different parts of the model.

For reproduction with the original source and complete saved fits,
follow `dev/simstudy/vignette-lesson/README.md`. Using the same seed
with a different package revision does not guarantee the original
results. The lessons record source and fit hashes and keep rendering
separate from fitting.

## What to inspect next

Start with the diagnostics and probability comparisons in Lesson 1, then
the coefficient, trait and community-output figures in Lesson 3. An
interval crossing zero does not demonstrate that an ecological effect is
absent; an interval excluding zero does not demonstrate that its
magnitude is accurate.

For a particular function, use the index at the end of Lesson 3 or the
package’s R help. Genuine new-site prediction, spatial validation and
comparisons with other JSDM packages require their own matched examples
and are not established by the current non-spatial lessons.

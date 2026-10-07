Subject: occJSDM v0.1.0 beta: occupancy modelling for eDNA metabarcoding

We have released the beta of occJSDM, an R package combining a joint species distribution model with the two-stage eDNA occupancy model of Ji et al. (2025). It estimates occupancy while accounting for false negatives and false positives in field collection and laboratory detection. It returns species-specific true-positive and false-positive field collection probabilities and lab true-positive and false-positive detection probabilities, with the lab probabilities estimated separately for each primer. Studies using multiple primers can therefore be analysed as a single dataset.

The package supports environmental and collection covariates, species traits, ordination, residual species correlations, prediction at new sites, optional spline-based nonlinear environmental responses, and spatial effects. Spatial inference is still limited in the beta.

occJSDM uses site and sample identifiers to identify the study design. It fits one-stage occupancy models for studies with repeated field samples but only one PCR observation per field sample, or JSDM-only models for one observation per site. With repeated field samples per site, pooling PCR products before sequencing such that only one observation remains per field sample leads to a one-stage model.

The one- and two-stage eDNA models analyse binary detections. Read counts can be supplied directly and are converted to detections or non-detections at a user-specified threshold (default 1, meaning even one read counts as a potential detection). Users should not set low-read detections to zero solely because their read counts are low: occJSDM uses all detections to infer false positives. Standard sequence-quality and taxonomy checks still apply. JSDM-only models accept binary presence/absence or continuous observations, using a Gaussian response model for continuous data. Count-response JSDMs and continuous-intensity detection models are not supported in this beta.

occJSDM's false-positive inference assumes that practitioners have applied careful field and lab practice, so that contamination probabilities are expected to be truly low. If contamination is high in the dataset, check for disagreement between chains, which can settle on different explanations of the same data.

Validation is ongoing. Occupancy probabilities can be pulled towards the middle, interval coverage has not been established, and spatial fields can be poorly recovered when sites are far apart relative to the spatial range. Single-species occupancy fits currently encounter a dimension-handling error; its correction is deferred until after beta. Please read the [Known limitations](https://github.com/AlexDiana/occJSDM#known-limitations) before using estimates in an analysis.

A quickstart and a simulator guide are available. Vignette Lessons 0-7 are first drafts awaiting a line-by-line review.

An optional MCP pilot using [Paper2Agent](https://github.com/jmiao24/Paper2Agent) exposes four tools: `validate_data`, `fit_model`, `diagnostics` and `summarise_fit`. It supports non-spatial one-stage occupancy, two-stage eDNA and binary JSDM fits. Continuous JSDM and spatial fitting are available through the ordinary R interface. The pilot uses the verified scientific revision `b7b7001`, pre-dating the beta, and awaits review of scientific explanations and teaching material. Ordinary R installation does not need MCP or Python. Pilot source, a separate download and setup instructions are linked from the release.

Install the beta with its vignettes:

```r
remotes::install_github("AlexDiana/occJSDM@v0.1.0-beta", build_vignettes = TRUE)
vignette("occJSDM", package = "occJSDM")
vignette("simulateOccJSDMData", package = "occJSDM")
```

[Release notes and downloads](https://github.com/AlexDiana/occJSDM/releases/tag/v0.1.0-beta). Feedback and bug reports are welcome through [GitHub issues](https://github.com/AlexDiana/occJSDM/issues). This is a GitHub beta; CRAN submission and broader interval-calibration work follow later.

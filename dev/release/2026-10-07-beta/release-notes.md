# occJSDM 0.1.0 beta

This GitHub prerelease combines joint species distribution modelling with the two-stage eDNA occupancy model of Ji et al. (2025). It estimates occupancy while allowing for false negatives and false positives in field collection and laboratory detection.

The package supports environmental and collection covariates, species traits, nonlinear environmental responses, spatial effects, ordination, residual species correlations, variation partitioning and prediction at new sites. Simpler survey designs support classical occupancy and JSDM-only models.

Recent corrections address collection-covariate alignment, random-number reproducibility, residual correlations, spatial sampling, read thresholds and environmental-response outputs. The shipped example is refitted with the current implementation. A quickstart and a simulator guide are available. Lessons 0-7 are first drafts only and are awaiting a full line-by-line review.

The beta was refreshed on 7 October 2026 to include the merged teaching rebuild, site-WAIC assessment and optional MCP execution pilot. The R model code, priors and shipped example data are unchanged from the initial beta. The numbered lessons remain excluded from the installed package and lesson site while their full review is pending.

The optional MCP pilot exposes validation, non-spatial fitting, diagnostics and fit summaries through the existing R implementation. Its source is included under `mcp/paper2agent/` in the tagged repository, and `occJSDM-cloud-pilot.zip` is a separate release download with setup instructions and the required scientific source payload. Ordinary R installation does not need MCP or Python. The pilot deliberately uses the verified scientific revision `b7b7001`, predating the beta; changing that runtime requires new native references and verification. Cloud checks and unscored model rehearsals establish an execution checkpoint, while scientific and Indonesian review and workshop adoption remain pending. Students must learn from the vignettes and lessons.

`extractWAIC()` now computes observed-data site-level WAIC by default for non-spatial binary, occupancy and two-stage fits. `computeSiteWAIC()` exposes its diagnostics and `compareSiteWAIC()` compares fits of the same observations. Spatial and continuous fits require explicit legacy extraction, `extractWAIC(fit, type = "legacy")`, whose stored score is unsuitable for choosing models for independent new-site predictions. Older occupancy/two-stage fits require their actual fitting threshold to be supplied. Check the numerical and WAIC reliability warnings before interpreting a comparison; factor-count selection has not yet been validated.

Validation is continuing. Occupancy probabilities can be pulled towards the middle; wider baseline priors do not generally fix this, and stronger collection priors can hide real effects. Interval coverage has not been established, and spatial fields can be poorly recovered when sites are far apart relative to the spatial range. High contamination can make chains settle on different explanations of the same data. Read the [Known limitations](https://github.com/AlexDiana/occJSDM#known-limitations) before using estimates in an analysis.

Install this fixed beta revision with its vignettes:

```r
remotes::install_github("AlexDiana/occJSDM@v0.1.0-beta", build_vignettes = TRUE)
vignette("occJSDM", package = "occJSDM")
vignette("simulateOccJSDMData", package = "occJSDM")
```

Please report problems through [GitHub issues](https://github.com/AlexDiana/occJSDM/issues). This is a GitHub beta; CRAN submission and broader interval-calibration work follow later.

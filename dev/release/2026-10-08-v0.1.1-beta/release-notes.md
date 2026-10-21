# occJSDM 0.1.1 beta

The second GitHub prerelease of occJSDM. The model code, priors and shipped example data are unchanged from 0.1.0, so existing fits and results stand. This release adds a posterior-draws accessor, fixes single-species fits and a plotting side effect, and starts a changelog in `NEWS.md`.

## New features

-   `returnPosteriorDraws()` returns the posterior draws of one parameter as a labelled array with chains kept separate: occupancy intercepts and slopes, collection coefficients, PCR true- and false-positive rates, or field-contamination rates. Use it for your own convergence checks and plots.
-   `plotTraceplot()` now labels its panels with the array's own covariate and species names by default, so it works directly on the output of `returnPosteriorDraws()`.

## Bug fixes

-   Single-species occupancy and two-stage fits now run. Previously a one-column state matrix lost its matrix shape inside the sampler and the fit failed.
-   `plotVariancePartitioning()` no longer changes your active ggplot2 theme, and with ggtern 4.0.0 it no longer leaves ordinary ggplot2 themes failing after a ternary plot.

## Documentation

-   `?runOccJSDM` now gives the actual default of `n_lattrait` and the values `threshold` accepts.
-   The quickstart explains the field-contamination prior and what is known about its interval coverage.
-   Help pages now say "variation partitioning". The function names are unchanged.

The full list is in [NEWS.md](https://github.com/AlexDiana/occJSDM/blob/v0.1.1-beta/NEWS.md). The known limitations from 0.1.0 still apply, including undercoverage of some intervals and weak recovery of spatial fields; see the README's [Known limitations](https://github.com/AlexDiana/occJSDM#known-limitations). Lessons 0-7 remain first drafts under review and are not included in the installed package.

Install this release with its vignettes:

```r
remotes::install_github("AlexDiana/occJSDM@v0.1.1-beta", build_vignettes = TRUE)
vignette("occJSDM", package = "occJSDM")
```

This release was checked on R 4.6.1 (macOS, Apple silicon). Please report problems through [GitHub issues](https://github.com/AlexDiana/occJSDM/issues).

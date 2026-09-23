# occJSDM

An R package for fitting a combined occupancy and joint species distribution model (occJSDM), optionally accounting for environmental and detection covariates, species traits, spatial autocorrelation, and for eDNA-style data, a two-stage observation process (false-negative and false-positive detection errors in the field and in the lab).

#### **N.B. This is beta software, and we are still in bugfixing mode.** 

## Installation

``` r
# install.packages("remotes")
remotes::install_github("AlexDiana/occJSDM", build_vignettes = TRUE)
```

Note the `build_vignettes = TRUE` -- without it, `remotes::install_github()` skips building vignettes by default, and `vignette("occJSDM", package = "occJSDM")` will report that no vignette was found.

## Getting started

Start with the teaching lessons. They show readable R code alongside the simulated truth and matching results:

- [Lesson 0 (optional): Create and explore a simulated survey](vignettes/occJSDM-lesson-0.md).
- [Lesson 1: Fit the model and compare its answers with truth](vignettes/occJSDM-lesson-1.md).
- [Lesson 2: Spatial landscapes and dispersal](vignettes/occJSDM-lesson-2.md), a plain-language account of what the spatial field learns, how to choose sample spacing and extent, and how to validate spatial prediction. Its worked example is pending review of the spatial submodel in PR #8.
- [Lesson 3: Understand the model outputs](vignettes/occJSDM-lesson-3.md), with true and fitted environmental effects, trait effects, species associations, variation partitioning and detection-effort curves.
- [Lesson 4: Compare four JSDMs](vignettes/occJSDM-lesson-4.md), a worked pure-JSDM pilot with occJSDM, gllvm, sjSDM and Hmsc, comparing probabilities and environmental responses with simulation truth. sjSDM turns out to have two local optima, which the lesson explains.

After installing with vignettes, open the completed lessons in R:

``` r
vignette("occJSDM-lesson-0", package = "occJSDM")
vignette("occJSDM-lesson-1", package = "occJSDM")
vignette("occJSDM-lesson-3", package = "occJSDM")
vignette("occJSDM-lesson-4", package = "occJSDM")
```

The quickstart guide and simulator reference are also available:

``` r
vignette("occJSDM", package = "occJSDM")
vignette("simulateOccJSDMData", package = "occJSDM")
```

## How to cite

If you use `occJSDM`, please cite the methods paper describing the underlying two-stage occupancy model, and cite this repository for the specific software implementation/version used.

**Methods paper (primary citation):**

> Ji, Y., Diana, A., Li, X., Matechou, E., Griffin, J. E., Liu, S., Luo, M., Wu, C., Bai, R., Yao, C., Yin, T., Dong, F., Wu, F., Wang, K., Yu, Z., Chen, X., Jiang, X., Che, J., Yu, D. W., & Popescu, V. D. (2025). High Quality, Granular, Timely, Trustworthy and Efficient Vertebrate Species Distribution Data Across a 30,000 km<sup>2</sup> Protected Area Complex. *Ecology Letters*, *28*(12), e70302. <https://doi.org/10.1111/ele.70302>

**Software (secondary citation):**

> Diana, A., & Yu, D. W. (2026). occJSDM: Occupancy Joint Species Distribution Models (Version 0.1.0) [Computer software]. <https://github.com/AlexDiana/occJSDM>

A machine-readable citation is also available in [`CITATION.cff`](CITATION.cff) -- GitHub uses this to populate the "Cite this repository" button in the sidebar of the repo page.

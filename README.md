# occJSDM

An R package for fitting occupancy joint species distribution models
(JSDMs) to eDNA metabarcoding reads data, accounting for a two-stage
detection process (field collection + PCR amplification/lab detection),
species traits, and optional spatial autocorrelation.

## Installation

```r
# install.packages("remotes")
remotes::install_github("AlexDiana/occJSDM")
```

## Getting started

See the package vignettes for a walkthrough of fitting a model with
`runOccJSDM()` and simulating data with `simulateOccJSDMData()`:

```r
vignette("occJSDM", package = "occJSDM")
vignette("simulateOccJSDMData", package = "occJSDM")
```

## How to cite

If you use `occJSDM`, please cite the methods paper describing the
underlying two-stage occupancy model, and cite this repository for the
specific software implementation/version used.

**Methods paper (primary citation):**

> Ji, Y., Diana, A., Li, X., Matechou, E., Griffin, J. E., Liu, S., Luo,
> M., Wu, C., Bai, R., Yao, C., Yin, T., Dong, F., Wu, F., Wang, K., Yu,
> Z., Chen, X., Jiang, X., Che, J., Yu, D. W., & Popescu, V. D. (2025).
> High Quality, Granular, Timely, Trustworthy and Efficient Vertebrate
> Species Distribution Data Across a 30,000 km<sup>2</sup> Protected Area
> Complex. *Ecology Letters*, *28*(12), e70302.
> <https://doi.org/10.1111/ele.70302>

**Software (secondary citation):**

> Diana, A., & Yu, D. W. (2026). occJSDM: Occupancy Joint Species
> Distribution Models (Version 0.1.0) [Computer software].
> <https://github.com/AlexDiana/occJSDM>

A machine-readable citation is also available in
[`CITATION.cff`](CITATION.cff) -- GitHub uses this to populate the "Cite
this repository" button in the sidebar of the repo page.

# occJSDM

An R package for fitting a combined occupancy and joint species distribution model (occJSDM), optionally accounting for environmental and detection covariates, species traits, spatial autocorrelation, and for eDNA-style data, a two-stage observation process (false-negative and false-positive detection errors in the field and in the lab).

#### **N.B. This is beta software, and we are still in bugfixing mode.** 

Validation is still in progress. Read the [Known limitations](#known-limitations) before relying on the results.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("AlexDiana/occJSDM", build_vignettes = TRUE)
```

Note the `build_vignettes = TRUE` -- without it, `remotes::install_github()` skips building vignettes by default, and `vignette("occJSDM", package = "occJSDM")` will report that no vignette was found.

## Getting started

The [quickstart](vignettes/occJSDM.md) shows the input data, a fitting call and a first look at the output. The simulator guide shows how to simulate a survey whose truth you know. After installing with vignettes, open them in R:

``` r
vignette("occJSDM", package = "occJSDM")
vignette("simulateOccJSDMData", package = "occJSDM")
```

Teaching lessons that compare every output with simulated truth are being reviewed and will be published after the beta.

## Known limitations

occJSDM is beta software and validation is still in progress. These limitations come from simulation studies of particular designs, so treat the numbers as examples rather than as the errors to expect in your own survey.

- **Occupancy probabilities are pulled towards the middle.** Low probabilities are estimated too high and high ones too low. The bias is strongest with limited replication, and it appears even when the overall average is close to the truth. In non-spatial two-stage (eDNA) simulations at 100 sites (two field samples per site, six PCR replicates per primer), true probabilities below 20% were overestimated on average by about 17 to 20 percentage points and those above 80% underestimated by about 20 to 25, depending on contamination; more field samples or sites reduced this without removing it, and perfect presence/absence data still left about 10 points at both extremes. See the [current-code comparison](https://github.com/AlexDiana/occJSDM/blob/main/dev/simstudy/current-main-recheck/REPORT.md) and its [summary table](https://github.com/AlexDiana/occJSDM/blob/main/dev/simstudy/current-main-recheck/results/paired-summary.csv).
- **Widening the occupancy-baseline prior is not a general fix.** The experimental `listPriors$sigma_b0` option sets the prior standard deviation (SD) of each species' baseline occupancy on the logit scale; the default is 1 and larger values widen the prior. In simulations, SD 2, 3 and 5 reduced the overestimation of low probabilities, clearly in spatial binary presence/absence fits, slightly in non-spatial binary fits and by up to about 3 points in two-stage (eDNA) fits. In two-stage fits, however, they made the MCMC chains mix worse and increased the bias of probabilities between 20% and 80%, by under a point at SD 2 and by more at 3 and 5, so the default is unchanged. Occupancy-only and continuous fits were not tested. See the [study report](https://github.com/AlexDiana/occJSDM/blob/main/dev/simstudy/occupancy-intercept-prior/REPORT.md) and `?runOccJSDM`.
- **Stronger collection priors can hide real collection effects.** Tightening the prior on collection-covariate slopes below its default variance of 2 weakened real collection effects without fixing the occupancy bias. See the [archived bias recheck and prior comparisons](https://github.com/AlexDiana/occJSDM/blob/main/dev/simstudy/nonspatial-bias-recheck.md).
- **Interval coverage has not been established.** 95% posterior intervals for site-level occupancy probabilities can contain the truth less often than 95% of the time. With the default prior in the [occupancy-baseline prior study](https://github.com/AlexDiana/occJSDM/blob/main/dev/simstudy/occupancy-intercept-prior/REPORT.md), two-stage (eDNA) fits did so about 72% to 79% of the time for true probabilities below 20% or above 80% and 90% to 93% in between; in spatial binary fits the figures were about 45% to 47% and 68%. Earlier checks on older code recorded intervals that were mostly too wide for the field false-positive rate `theta0` and too narrow for `B0` in continuous-response models. A full interval-calibration study is planned after the beta; see [Interval calibration](https://github.com/AlexDiana/occJSDM/blob/main/TODO.md#interval-calibration) in the to-do list.
- **Spatial fields are poorly recovered when sites are far apart relative to the spatial range.** In nine simulated spatial binary communities at 100 sites, the fitted fields were nearly flat: each site gave one presence/absence record per species, at the shortest range most sites had no neighbour with correlation above 0.5, and even a check given the true amplitude, range and coefficients reduced field error only about 6% below that of a flat field. Rare species were hardest, and the experimental half-Cauchy amplitude prior (`sigma_bs_prior`) did not help. See the [spatial-amplitude comparison](https://github.com/AlexDiana/occJSDM/blob/main/dev/simstudy/spatial-amplitude-prior/README.md) and its [diagnosis](https://github.com/AlexDiana/occJSDM/blob/main/dev/simstudy/spatial-amplitude-prior/diagnosis/REPORT.md).

## How to cite

If you use `occJSDM`, please cite the methods paper describing the underlying two-stage occupancy model, and cite this repository for the specific software implementation/version used.

**Methods paper (primary citation):**

> Ji, Y., Diana, A., Li, X., Matechou, E., Griffin, J. E., Liu, S., Luo, M., Wu, C., Bai, R., Yao, C., Yin, T., Dong, F., Wu, F., Wang, K., Yu, Z., Chen, X., Jiang, X., Che, J., Yu, D. W., & Popescu, V. D. (2025). High Quality, Granular, Timely, Trustworthy and Efficient Vertebrate Species Distribution Data Across a 30,000 km<sup>2</sup> Protected Area Complex. *Ecology Letters*, *28*(12), e70302. <https://doi.org/10.1111/ele.70302>

**Software (secondary citation):**

> Diana, A., & Yu, D. W. (2026). occJSDM: Occupancy Joint Species Distribution Models (Version 0.1.0) [Computer software]. <https://github.com/AlexDiana/occJSDM>

A machine-readable citation is also available in [`CITATION.cff`](CITATION.cff) -- GitHub uses this to populate the "Cite this repository" button in the sidebar of the repo page.

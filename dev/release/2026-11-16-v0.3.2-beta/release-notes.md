# occJSDM 0.3.2 beta

The sixth GitHub prerelease of occJSDM. It removes a cost at the end of a fit that grew with the square of the number of coefficients. Estimates are unchanged from 0.3.1.

## Who is affected, and what to do

-   **Affected:** fits with two or more chains and thousands of coefficients in one block. In practice that is the occupancy coefficients of a fit with many covariates, many species, or both, since that block has one coefficient for each covariate and species. After the sampler had finished, the convergence summary that `runOccJSDM()` prints took minutes and gigabytes: with 2 chains of 5,000 draws, about 2 minutes and 3 GB for 3,000 coefficients and 8 minutes and 8 GB for 6,000, and for 12,500 it did not finish within 16 GB, so a fit of that size could run out of memory after sampling and be lost. In this release the summary for 12,500 coefficients takes about 4 seconds and 13 MB. **Update before running fits of this size.** Nothing needs refitting.
-   **Not affected:** results. The summary reports the same R-hat and effective sample size and raises the same warnings, and a fit stores exactly what it stored under 0.3.1.
-   **If you have not yet updated to 0.3.1:** this release includes its fix for fits with one row per site whose `Site` column was not in sorted order. The [0.3.1 notes](https://github.com/AlexDiana/occJSDM/releases/tag/v0.3.1-beta) say which models must be refitted.

The full list is in [NEWS.md](https://github.com/AlexDiana/occJSDM/blob/v0.3.2-beta/NEWS.md). The known limitations still apply; see the README's [Known limitations](https://github.com/AlexDiana/occJSDM#known-limitations).

Install this release with its vignettes:

```r
remotes::install_github("AlexDiana/occJSDM@v0.3.2-beta", build_vignettes = TRUE)
vignette("occJSDM-00-quickstart", package = "occJSDM")
```

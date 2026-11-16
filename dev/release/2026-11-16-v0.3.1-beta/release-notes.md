# occJSDM 0.3.1 beta

The fifth GitHub prerelease of occJSDM. It fixes a defect that gave wrong results in every earlier beta for one kind of data. Nothing else about fitting has changed.

## Who is affected, and what to do

-   **Affected:** fits with one row per site, that is presence/absence (binary) or continuous data with no repeated samples, whose `Site` column was not in sorted order. Site names in file order such as `FHT1`, `FHT2`, ..., `FHT10` are not in sorted order, because `FHT10` sorts before `FHT2`. In these fits every earlier beta (`v0.1.0-beta` to `v0.3.0-beta`) attached the occupancy covariates and the site coordinates to the wrong sites, without an error or a warning. Covariate effects, the spatial field, the hidden environmental covariates, predictions and site WAIC were all wrong. **Update to this release and refit these models.**
-   **Not affected:** one-row-per-site fits whose data had no `Site` column, or whose `Site` was already in sorted order (for example `1:n`, or zero-padded names such as `S001` in file order). Their estimates are unchanged, draw for draw.
-   **Never affected:** occupancy and two-stage (eDNA) fits, which have repeated rows per site.

## Other change

-   A fit with one row per site now stores its site names in `infos$siteNames`, as occupancy and two-stage fits already did, so ordination scores and plots label its sites by name rather than by row number. With no `Site` column it still stores `1:n`. `?runOccJSDM` has a new section, *Row order of the results*.

The full list is in [NEWS.md](https://github.com/AlexDiana/occJSDM/blob/v0.3.1-beta/NEWS.md). The known limitations still apply; see the README's [Known limitations](https://github.com/AlexDiana/occJSDM#known-limitations).

Install this release with its vignettes:

```r
remotes::install_github("AlexDiana/occJSDM@v0.3.1-beta", build_vignettes = TRUE)
vignette("occJSDM-00-quickstart", package = "occJSDM")
```

# occJSDM 0.2.0 beta

The third GitHub prerelease of occJSDM. It adds parallel chains and a faster spatial sampler. The model and priors are unchanged, but seeded fits give different numbers from 0.1.1, so refit anything you need to reproduce.

## New features

-   `runOccJSDM()` can run its chains on separate cores: set `ncores` in `MCMCparams`, for example `MCMCparams = list(nchain = 4, ncores = 4)`. `ncores` is the total core budget, so four chains with `ncores = 4` run on four single-threaded workers and `ncores = 8` gives each of the four workers two threads. The default, `ncores = 1`, runs the chains one after another. The returned `infos` records the sampling settings, the seeds of every chain and how the fit was executed.

## Behaviour changes

-   Every chain is now seeded on its own, so a fit is identical whatever `ncores` says (stated for a single-threaded BLAS) and no chain continues another's random stream. A fit under a given `set.seed()` value therefore differs from 0.1.1. The shipped `sampleresults` was refitted under the new seeding.
-   The spatial range step and the spatial coefficient draw form their weighted Gram matrices with the symmetric BLAS routine, which halves that arithmetic: 1.25 to 1.45 times faster per iteration on a 900-site spatial fit, depending on the BLAS. The result agrees to rounding error, not bit for bit. The shipped `sampleresults` is a spatial fit made before this change, so refitting it with the same seed now gives draws that differ at rounding level.

## Documentation

-   The README explains how to use parallel chains and how to link R to an optimised BLAS, which decides how long a spatial fit takes (about 25 times faster per iteration at 900 sites with Apple's Accelerate).

The full list is in [NEWS.md](https://github.com/AlexDiana/occJSDM/blob/v0.2.0-beta/NEWS.md). The known limitations from the earlier betas still apply, including undercoverage of some intervals, occupancy probabilities pulled towards the middle, and weak recovery of spatial fields; see the README's [Known limitations](https://github.com/AlexDiana/occJSDM#known-limitations). The experimental community mean for the occupancy intercepts is not in this release. Lessons 0-7 remain first drafts under review and are not included in the installed package.

Install this release with its vignettes:

```r
remotes::install_github("AlexDiana/occJSDM@v0.2.0-beta", build_vignettes = TRUE)
vignette("occJSDM", package = "occJSDM")
```

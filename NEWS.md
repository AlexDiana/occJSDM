# occJSDM (development version)

## New features

-   `runOccJSDM()` can run its chains on separate cores: set `ncores` in `MCMCparams`. `ncores` is the total core budget, so `ncores = 4` runs four chains on four single-threaded workers and `ncores = 8` on four workers of two RcppParallel threads each. A parallel fit prints the `tail -f` command for its progress file. The returned `infos` records the sampling settings, the seeds of every chain and how the fit was run.

-   Experimental, off by default: `listPriors = list(b0_prior = "community")` centres the prior on each species' occupancy intercept on an estimated community mean, `mu_b0`, instead of on 0, with a normal prior on that mean whose standard deviation, `mu_b0_sd`, defaults to 2. Its draws are returned in `jsdm_output$mu_b0_output` and reported by `returnConvergenceDiagnostics()`, and the prior is recorded in `infos$intercept_prior`. In occupancy and two-stage fits the sampler also proposes, each iteration, a common shift of every species' occupancy intercept and `mu_b0`, accepted with the latent occupancy states integrated out, so that the community mean mixes at a realistic number of species; `b0_shift`, `b0_shift_sd` and `b0_shift_proposals` in `MCMCparams` control it and `infos$b0_shift` records its acceptance rates. Default fits are unchanged. The option is under review and is not yet a recommended setting; see `?runOccJSDM`.

## Behaviour changes

-   Every chain is now seeded on its own, so a fit is identical whatever `ncores` says and no chain continues another's random stream. Because of this, a fit under a given `set.seed()` value gives different numbers from version 0.1.1; the shipped `sampleresults` was refitted under the new seeding with the same seed and settings. Identity across `ncores` is stated for a single-threaded BLAS; see `?runOccJSDM`.
-   The spatial range step and the spatial coefficient draw now form their weighted Gram matrices with the symmetric BLAS routine (the rows of the design scaled by the square root of the weights, then multiplied by their own transpose), which halves that arithmetic. The result is the same quantity to rounding error, not bit for bit, so a seeded spatial fit no longer reproduces the draws of the previous build exactly; measured on a 900-site binary spatial fit, an iteration is 1.25 to 1.45 times faster depending on the BLAS.

## Deprecations

-   `returnVariancePartitioning()` and `plotVariancePartitioning()` are renamed `returnVariationPartitioning()` and `plotVariationPartitioning()`, to match the documentation. The old names still work but give a deprecation warning and will be removed in a future release.

# occJSDM 0.1.1

Second beta, released on GitHub as `v0.1.1-beta`. The model code, priors and shipped example data are unchanged from 0.1.0.

## New features

-   `returnPosteriorDraws()` returns the posterior draws of one parameter as a labelled array with chains kept separate: occupancy intercepts and slopes, collection coefficients, PCR true- and false-positive rates, or field-contamination rates. Use it for your own convergence checks and plots.
-   `plotTraceplot()` now labels its panels with the array's own covariate and species names by default, so it works directly on the output of `returnPosteriorDraws()`.

## Bug fixes

-   Single-species occupancy and two-stage fits now run. A one-column state matrix lost its matrix shape when expanded from sites to field samples, which broke the field false-positive sampler.
-   `plotVariancePartitioning()` no longer changes your active ggplot2 theme. ggtern is now loaded only when the function is called, and with ggtern 4.0.0 the function corrects two invalid tick-length defaults that otherwise made ordinary ggplot2 themes fail after a ternary plot.

## Documentation

-   `?runOccJSDM` now gives the actual default of `n_lattrait` and the values `threshold` accepts.
-   The quickstart explains the field-contamination prior and what is known about its interval coverage.
-   Help pages and lessons now say "variation partitioning". The function names are unchanged.

## Internal changes

-   Removed unused internal R and C++ helpers. No exported function changed; the removed code is kept under `deprecated/` in the repository, outside the package.

# occJSDM 0.1.0

First public beta, released on GitHub as `v0.1.0-beta` on 7 October 2026. It combines joint species distribution modelling with the two-stage eDNA occupancy model of Ji et al. (2025), and ships a quickstart and a simulator guide. See the [release notes](https://github.com/AlexDiana/occJSDM/releases/tag/v0.1.0-beta) for its features and known limitations.

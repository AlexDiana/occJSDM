# occJSDM (development version)

In development, not yet released.

## Bug fixes

-   The Polya-Gamma draws that `runOccJSDM()` makes in binary, occupancy and two-stage fits had half the correct variance in every cell whose linear predictor was beyond 12 in absolute value, which is a fitted probability below about 6 in a million or above about 0.999994. Their mean was correct. Such cells are species that the model is almost certain are absent from a site, or present at it. In the Quickstart fit they are 0.05% of the occupancy cells drawn, with at least one in 44% of iterations, and 15 of the 37 million collection cells drawn. The draw in those cells is now an inverse Gaussian, which is what the Polya-Gamma distribution tends to for a large predictor, in place of a normal with half the variance. Every beta from `v0.1.0-beta` to `v0.3.2-beta` has the defect. In the Quickstart example its effect is too small to see: refitted with five seeds before the fix and five after, no posterior mean moves by more than 0.07 posterior standard deviations, which is less than two seeds differ by. Its effect has not been measured on data in which a species is nearly absent or nearly everywhere. A fit that contains such a cell now gives different numbers from the same seed, from the first such cell onwards; a fit with none is unchanged, draw for draw. The continuous model makes no Polya-Gamma draws and is unaffected.

# occJSDM 0.3.2

Sixth beta, released on GitHub as `v0.3.2-beta`. It removes a cost at the end of a fit that grew with the square of the number of coefficients, which made fits with many covariates or many species slow to finish and could exhaust memory. Estimates are unchanged from 0.3.1.

## Bug fixes

-   The convergence summary that `runOccJSDM()` prints when a fit with two or more chains finishes (`computeDiagnostics()`) now takes time and memory in proportion to the number of parameters. It used to grow with the square of the size of the largest parameter block, so a fit with a few thousand occupancy coefficients (many covariates, many species, or both) spent minutes and several gigabytes after the sampler had finished, and a larger one could exhaust memory and lose the fit. The R-hat and effective sample size it reports are the same statistics as before, and the printed summary and its warnings are unchanged. `returnConvergenceDiagnostics()` is unaffected.

# occJSDM 0.3.1

Fifth beta, released on GitHub as `v0.3.1-beta`. It fixes a defect that made fits with one row per site wrong whenever the `Site` column was not in sorted order. Every earlier beta has the defect.

## Bug fixes

-   `runOccJSDM()` attached the occupancy covariates and the site coordinates to the wrong sites in a fit with one row per site (the binary and continuous models) whenever `data$info$Site` was not already in sorted order, without an error or a warning. The covariates were sorted by `Site` while the response kept its input order. Site names in file order such as `FHT1`, `FHT2`, ..., `FHT10` are not in sorted order, and with them a planted slope of 1.5 was estimated as about 0. Fits were correct when `data$info` had no `Site` column, or when `Site` was already sorted (for example `1:n`, or zero-padded names such as `S001` in file order), and those fits' estimates are unchanged, draw for draw. Occupancy and two-stage fits, which have repeated rows per site, were never affected. **Refit any one-row-per-site model whose `Site` column was not in sorted order:** its covariate effects, spatial field, hidden environmental covariates, predictions and site WAIC were all computed with covariates on the wrong rows. To check an existing fit, compare `fit$X_psi` with your standardised covariates row by row.

## Behaviour changes

-   A fit with one row per site now stores its site names in `infos$siteNames`, as occupancy and two-stage fits already did, so `returnOrdinationScores()` and `plotOrdinationScores()` label its sites by name. It used to store the row numbers `1:n` whatever the `Site` column held. With no `Site` column it still stores `1:n`. No estimate changes; code that used `infos$siteNames` of such a fit as a row index should use `seq_len(fit$infos$n)`.

## Documentation

-   `?runOccJSDM` has a new section, *Row order of the results*. With one row per site, every per-site result is in the row order of `data$info`. With repeated rows per site, the sites are in sorted order by `Site` and named in `infos$siteNames`.

# occJSDM 0.3.0

Fourth beta, released on GitHub as `v0.3.0-beta`. The community-mean prior on the occupancy intercepts is now the default for multi-species data, so fits that leave `b0_prior` unset give different numbers from 0.2.0.

## Behaviour changes

-   `runOccJSDM()` now centres the prior on each species' occupancy intercept on an estimated community mean by default (`listPriors$b0_prior = "community"`, added in 0.2.0) for binary, occupancy and two-stage data with more than one species. Simulation studies showed that it removes most of the overestimation of rare species' occupancy at little cost in ordinary communities. Continuous data and single-species fits keep the old prior, `b0_prior = "fixed"`, and asking for `"community"` with one species is now an error. Fits that leave `b0_prior` unset therefore give different numbers from 0.2.0, and their output gains `mu_b0_output`; set `listPriors = list(b0_prior = "fixed")` to reproduce them. The shipped `sampleresults` and `quickstartresults` have been refitted at the new default, with their previous seeds and settings, so their numbers differ from 0.2.0 too.
-   `runOccJSDM()` now rejects a read `threshold` below 1. Values between 0 and 1 used to be accepted; for integer read counts they gave the same detections as the default of 1.
-   `simulateOccJSDMData()` now names the rows of `jsdmParams_true$eta` and `z_true` by site and of `w_true` by sample, matching the `Site` and `Sample` columns of `info`, and names their columns by species, matching `OTU`. Truth can be joined to observations by identifier without labelling the matrices first. Only names are added; the simulated values are unchanged.

## Bug fixes

-   `computeSiteWAIC()`'s default `quadrature` now continues to orders 101 and 161 after 65. A draw with a large residual standard deviation could fail to converge by 65 and stop the whole calculation, as happened for the one-factor fit in Lesson 4 after the change of default prior. Draws that converge earlier are unaffected.
-   `runOccJSDM()` now decides between the one- and two-stage models by looking for repeated (`Site`, `Sample`) pairs, not repeated `Sample` values. A one-stage survey whose sample IDs restart at each site (1, 2, 3 at every site) used to be fitted as two-stage. Sample IDs may restart at each site in either design.
-   The error for `NA` observations outside the two-stage model no longer ends in a stray "FALSE", and now suggests deleting a row that is wholly missing.

## Documentation

-   The quickstart vignette now tours every main result: occupancy at each site, environmental and trait-by-environment effects, collection, detection and false-positive rates, residual associations and ordination, variation partitioning, expected species richness, prediction at new sites and the value of extra samples or PCRs, each also saved as a table. It no longer explains the model; the lessons do.
-   New example data for it: `quickstartdata`, a simulated river eDNA survey modelled on Ji et al. (2025), with a categorical occupancy covariate, and `quickstartresults`, its fit without a spatial field. `sampledata` and `sampleresults` are unchanged.

# occJSDM 0.2.0

Third beta, released on GitHub as `v0.2.0-beta`. It adds parallel chains and a faster spatial sampler. Seeded fits give different numbers from 0.1.1.

## New features

-   `runOccJSDM()` can run its chains on separate cores: set `ncores` in `MCMCparams`. `ncores` is the total core budget, so `ncores = 4` runs four chains on four single-threaded workers and `ncores = 8` on four workers of two RcppParallel threads each. A parallel fit prints the `tail -f` command for its progress file. The returned `infos` records the sampling settings, the seeds of every chain and how the fit was run.

-   Experimental, off by default: `listPriors = list(b0_prior = "community")` centres the prior on each species' occupancy intercept on an estimated community mean, `mu_b0`, instead of on 0, with a normal prior on that mean whose standard deviation, `mu_b0_sd`, defaults to 2. Its draws are returned in `jsdm_output$mu_b0_output` and reported by `returnConvergenceDiagnostics()`, and the prior is recorded in `infos$intercept_prior`. In occupancy and two-stage fits the sampler also proposes, each iteration, a common shift of every species' occupancy intercept and `mu_b0`, accepted with the latent occupancy states integrated out, so that the community mean mixes at a realistic number of species; `b0_shift`, `b0_shift_sd` and `b0_shift_proposals` in `MCMCparams` control it and `infos$b0_shift` records its acceptance rates. Default fits are unchanged. The option is under review and is not yet a recommended setting; see `?runOccJSDM`.

## Behaviour changes

-   Every chain is now seeded on its own, so a fit is identical whatever `ncores` says and no chain continues another's random stream. Because of this, a fit under a given `set.seed()` value gives different numbers from version 0.1.1; the shipped `sampleresults` was refitted under the new seeding with the same seed and settings. Identity across `ncores` is stated for a single-threaded BLAS; see `?runOccJSDM`.
-   The spatial range step and the spatial coefficient draw now form their weighted Gram matrices with the symmetric BLAS routine (the rows of the design scaled by the square root of the weights, then multiplied by their own transpose), which halves that arithmetic. The result is the same quantity to rounding error, not bit for bit, so a seeded spatial fit no longer reproduces the draws of the previous build exactly; the shipped `sampleresults` is a spatial fit made before this change, so refitting it with the same seed now gives draws that differ at rounding level; measured on a 900-site binary spatial fit, an iteration is 1.25 to 1.45 times faster depending on the BLAS.
-   A partial `MCMCparams` list in `runOccJSDM()` now fills a missing `nburn` or `niter` with 5,000, as the help page documents, instead of 1,000. For example, `list(nchain = 1)` used to run 1,000 burn-in and 1,000 kept iterations and now runs 5,000 of each, so it takes about five times as long. Fits that give both `nburn` and `niter`, or leave `MCMCparams` out, are unchanged. To keep the old length, give `nburn = 1000` and `niter = 1000`.

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

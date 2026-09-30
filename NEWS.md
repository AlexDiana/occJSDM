# occJSDM 0.1.0

This is the first public beta of occJSDM. Validation is still in progress; see Known limitations below before relying on the results.

## Opt-in prior options

The default priors are unchanged. Each option below takes effect only when it is set in `listPriors` for `runOccJSDM()`; see `?runOccJSDM` for the details.

- `sigma_b0` (experimental) sets the standard deviation of the Normal prior on each species' occupancy baseline `B0` (default 1). In a simulation study, wider values reduced the overestimation of low occupancy probabilities in binary fits but made two-stage (eDNA) fits mix worse; keep the default for occupancy and two-stage data.
- `sigma_bs_prior = "half_cauchy"` (experimental) puts a half-Cauchy prior, with scale `sigma_bs_scale` (default 1), on the residual spatial-coefficient standard deviation in spatial fits. The default `"inverse_gamma"` keeps the existing prior. The option does not establish improved spatial recovery.
- `tau_prior = "half_cauchy"` puts a half-Cauchy prior, with scale `tau_scale` (default 1), on the noise standard deviation for continuous responses. The default `"inverse_gamma"` keeps the existing prior, with `a_tau` and `b_tau` (both default 5).
- `b_betatheta_slope_var` sets the prior variance of the collection-covariate slopes (default 2). It is provisional; treat other values as a diagnostic rather than a recommended setting.

## Notable fixes

- Collection covariates are aligned with the field samples they describe ([#5](https://github.com/AlexDiana/occJSDM/pull/5)).
- Random draws in the sampler run on R's main thread, so `set.seed()` reproduces a fit on the same platform whatever thread count is requested ([#6](https://github.com/AlexDiana/occJSDM/pull/6), [#9](https://github.com/AlexDiana/occJSDM/pull/9)).
- The residual species covariance is preserved when the latent factors are reparameterised ([#7](https://github.com/AlexDiana/occJSDM/pull/7)).
- Read thresholds greater than one no longer erase qualifying detections.
- Spatial fits use a consistent complete basis for fitting and prediction, and a valid update of the spatial range ([#8](https://github.com/AlexDiana/occJSDM/pull/8)).

## Known limitations

Occupancy probabilities are pulled towards the middle, interval coverage has not been established, and sparse site designs carry little information about spatial fields. See [Known limitations](https://github.com/AlexDiana/occJSDM#known-limitations) in the README for the evidence.

# Optional continuous-noise prior: review and validation

**Review after PR #8.** This branch, `codex/continuous-noise-prior`, is based on the spatial correction so GitHub shows only the additional noise-prior work. The current inverse-gamma default is unchanged. Alex can review whether the half-Cauchy option should be available separately from any later decision to change the default. No binary or detection prior changes.

The runtime changes are in `read_noise_prior()` and its public wiring in `R/runOccJSDM.R`, plus `sample_tau_half_cauchy()` and its dispatch in `R/jsdmfun.R`. Read `tests/testthat/test-noise-prior.R` alongside them: it checks direct integration of the conditional density, changes of response units, validation of prior inputs and actual public fitting behaviour. The help file documents the new settings. The CSVs and reference script support the evidence below.

The [spatial review guide](spatial-beta-recheck.md) contains the basis/range correction and binary evidence. The [non-spatial investigation](https://github.com/AlexDiana/occJSDM/blob/codex/nonspatial-bias-recheck/dev/simstudy/nonspatial-bias-recheck.md) is independent of this branch. This split reorganises completed work without rerunning simulations or changing saved results.

## What the option does

- Continuous-response fits now accept an explicit noise prior through `listPriors`, and save the applied choice in `infos$noise_prior`. `tau_prior = "half_cauchy"` places a half-Cauchy prior on each species' noise standard deviation; `tau_scale` defaults to 1 in response units. `tau_prior = "inverse_gamma"` retains the existing family on noise variance, with configurable `a_tau` and `b_tau`, both defaulting to 5. The current default remains inverse-gamma while the default choice is being decided. These settings do not alter the binary or detection models.

The flexible noise option allows small noise values that the previous prior strongly discouraged. It does not force the noise towards the simulated truth. Tests cover residual scales 0.1, 1 and 3, and changing measurement units together with the prior scale.

## Why the noise prior matters

In the original continuous checks, each of 100 locations had only one observation per species, with true noise SD 0.1 and spatial SD 1. The existing inverse-gamma prior produced mean noise estimates 0.787, 0.619 and 0.552 across the three spatial ranges. The estimated spatial fields were much too weak: their slopes against the true fields were 0.415, 0.654 and 0.725, where a slope of 1 would represent the correct magnitude.

An experimental half-Cauchy noise prior, keeping the same 99 support points, changed those field slopes to 0.842, 0.894 and 0.943. Noise estimates fell to 0.357, 0.278 and 0.116. Allowing all 100 locations as support points then gave field slopes 0.918, 0.950 and 0.945, field RMSE 0.169, 0.173 and 0.224, and noise estimates 0.222, 0.150 and 0.110. These are paired mechanism checks, not proof of unbiased estimation. Several individual coefficient and noise traces still had convergence warnings; one observation per location supplies limited information for distinguishing noise from a short-range field.

An independent calculation for the middle range fixes the regression means at their true values and integrates the range, spatial SD and noise SD in observation space. With 99 support points, its mean noise SD is 0.611 under the original prior and 0.277 under the half-Cauchy prior. Using the full generating GP covariance with the half-Cauchy prior reduces this to 0.159. Thus both the prior and the spatial approximation contribute; the remaining noise inflation is present in the specified posterior, rather than demonstrating a faulty noise sampler. Changing the spatial-variance prior alone has a much smaller effect in this comparison.

The original half-Cauchy experiments used an exact rejection sampler for the conditional noise distribution. The public option uses an auxiliary-variable Gibbs update with the same stationary distribution and no rejection loop. It has been checked by direct numerical integration and independent mathematical review. Let `v = tau^2` and `A = tau_scale`. The [half-Cauchy mixture identity](https://arxiv.org/abs/1508.03884) gives

```text
v | auxiliary ~ inverse-Gamma(1/2, 1/auxiliary)
auxiliary ~ inverse-Gamma(1/2, 1/A^2)
auxiliary | v ~ inverse-Gamma(1, 1/v + 1/A^2)
v | auxiliary, residuals ~ inverse-Gamma((n+1)/2, RSS/2 + 1/auxiliary).
```

Here `inverse-Gamma(a,b)` means that the reciprocal has a Gamma distribution with shape `a` and rate `b`. Each update draws the auxiliary variable from the current noise SD, then draws the new variance. The auxiliary is refreshed before every noise update. No detection prior is relaxed by this change.

## Informative continuous checks

The public half-Cauchy option was checked in three independent continuous datasets with 100 locations, ten measurements per location, 12 species, true noise SD 0.1 and 100 support points. These use two chains with 300 burn-in and 500 retained iterations. All three generating ranges were recovered. Noise estimates are 0.1006, 0.0993 and 0.1011, and spatial SD estimates are 0.997, 0.992 and 0.972. Total predicted-response RMSE is 0.0314 to 0.0319, with mean error between -0.0012 and +0.0022 in response units. No fitting warnings were recorded.

After centering each species' field around its own average, spatial recovery slopes are 0.994, 0.999 and 1.000. A slope of 1 means that the estimated spatial variation has the correct strength. Uncentered field RMSE is larger, 0.113 to 0.205, because the intercept and the field's average can offset one another; the total predicted level is estimated much more accurately. An [independent observation-space Gaussian integration](spatial-beta-continuous-range-reference.csv) at six dispersed variance states per dataset confirms the strongly concentrated range posteriors. The constant range traces are therefore supported by the posterior calculation, rather than taken as evidence of convergence on their own. [Compact continuous results](spatial-beta-continuous-results.csv) retain all three datasets' metrics.

## Reproduction and provenance

The runner saves data, truth, seeds, requested settings, the applied noise prior, source and loaded-library hashes, fits, warnings and metrics. It rejects an explicit noise-prior request if the chosen source ignores it. The existing source-default mode can reproduce older fits.

```sh
Rscript dev/simstudy/validate_spatial_continuous.R --source=/path/to/compiled-source --out=/path/to/new-continuous-results --cell=2 --repeats=10 --tau=.1 --knots=100 --tau-prior=half_cauchy --tau-scale=1 --burn=300 --iter=500 --chains=2
Rscript dev/simstudy/verify_spatial_continuous_range.R --input=/path/to/new-continuous-results --out=/path/to/new-range-reference
```

Use cells 1, 2 and 3 in separate fresh output directories. The original weak-information comparison uses one repeat and 99 or 100 knots, choosing either prior explicitly. The unchanged inverse-gamma settings are `--tau-prior=inverse_gamma --a-tau=5 --b-tau=5`. These are saved reproduction commands; they were not rerun during the PR split.

The long fits used immutable source snapshots containing the spatial correction and public half-Cauchy option; saved metadata confirms the option was applied. Later selector normalization and validation guards do not alter the sampling path used by those runs. The subsequently approved residual-correlation fix does not affect these zero-factor, zero-trait fits. Raw fits, data and provenance remain in the task's `work/spatial-noise-resume-20260913` directory.

## Historical verification and remaining decisions

Before the PR split, the full-support, grouped-arithmetic and noise-prior tests fail before their respective changes and pass afterward. Independent review checked the mathematics, public wiring and update order. It identified two cases in which a requested prior could be silently ignored; named prior selectors are now normalized, and the validation runner checks the applied prior. The source suite passed 498 assertions before the two added named-selector regressions, which subsequently passed with all 16 noise-prior assertions. After incorporating `main` on 13 September, 174 focused alignment, spatial and noise assertions passed. The 13 September freshly built installed-package check passed 464 assertions with zero test failures or test warnings; its eight skips cover opt-in coverage, CRAN-excluded recovery and source-only checks.

That 13 September package check reports zero errors, three existing warnings and three existing notes. The warnings concern the compiler's R-header warning option, the undocumented `verbose` argument in `predictNewSites.Rd`, and GNU Makevars syntax. The notes concern the worktree's hidden Git file, LICENSE metadata and existing undefined globals. No new package warning was introduced. This verification used arm64 macOS and R 4.5.0; it does not replace cross-platform checking.

On 14 September, after incorporating `main` at `80d449d` and the approved PR #7, the combined source suite passed 735 assertions with zero failures, errors or test warnings. The one skip is the opt-in coverage study. R also printed an environment warning that the installed `testthat` was built under R 4.5.2; the test run used R 4.5.0 and completed successfully. This was a fresh source-suite integration check, not a repeat of the installed-package check or long spatial simulations. The updated TODO removes accidentally committed merge-conflict text and records the three approved code fixes under *Fixed bugs* 49-51; the spatial review and separate beta checks remain open.

**Decision still needed:** Alex should review the optional prior and decide separately whether the continuous-model default should ever change. This PR keeps the inverse-gamma default and can be reviewed on that basis. Evidence supports the alternative in the tested small-noise settings, not every continuous dataset. Keep the earlier weak-information results and limitations visible. Rare species, support-point adequacy and other beta checks remain in the [consolidated TODO](https://github.com/AlexDiana/occJSDM/blob/codex/nonspatial-bias-recheck/TODO.md).

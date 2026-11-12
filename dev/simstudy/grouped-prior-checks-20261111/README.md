# Checks of the grouped covariate prior's default scales

12 November 2026. R 4.6.1, on the `claude/grouped-covariate-prior` branch at `3461213` (this check changes no package code). Task 8 of the part 1 plan for the [grouped covariate prior](../../superpowers/specs/2026-11-11-grouped-covariate-prior-design.md).

**The default axis scale fails the prior check.** The design asks for a prior median between 1 and 2 for the standard deviation across sites of the embedding's contribution to the linear predictor. With the default scales the median is **16.3 with 64 embedding columns** and **24.0 with 128**, both with 8 latent traits. No default was changed. The axis scale that brings the median to 1.5 is about 0.72 divided by the number of embedding columns, where the default is 1 divided by its square root.

## What was run

``` bash
Rscript dev/simstudy/grouped-prior-checks-20261111/prior_check.R
```

Nothing is fitted. Each of 2,000 draws is one draw from the prior of one species' coefficients on the embedding columns, with no measured traits:

- each axis scale from a half-Cauchy with scale `axis_scale`, and the embedding leftover scale from a half-Cauchy with scale 1;
- the community row of `G` with standard deviation equal to the axis scale, and 8 rows of `C` with standard deviation equal to the axis scale divided by the square root of 8;
- the species' 8 latent traits, standard Normal;
- a leftover with standard deviation equal to the leftover scale times the axis scale.

The coefficients are the community row, plus the latent traits times `C`, plus the leftover. They multiply a standardised matrix of 250 sites, built with the simulator's factor model for an embedding block, and the standard deviation across sites is recorded.

The script calls the package's own `read_coef_prior()` for the default scales and `coef_prior_sds()` for how the scales combine, so it checks the prior as the package defines it. Its header comment lists the source lines for every other part of the draw.

## Numbers

``` text
Prior check of the grouped covariate prior
R 4.6.1, 2000 prior draws, 250 sites, target median 1 to 2

K = 64 embedding columns, gt = 8
  defaults: axis_scale = 0.12500 (1 / sqrt(K)), sigma_emb_scale = 1, c_base = 0.35355
  default scales, scales drawn (the check)       median   16.280   5%    3.470   95%    325.585
  median inside 1 to 2: NO
  scales held at their prior scales              median    1.557   5%    0.886   95%      2.614
  independent columns, scales held               median    1.707   5%    1.395   95%      2.125
  independent columns, scales drawn              median   16.565   5%    4.350   95%    354.279
  axis_scale for a median of 1.5: 0.01152 = 0.737 / K = 0.092 / sqrt(K)
  axis_scale = 0.01152, fresh draws              median    1.557   5%    0.337   95%     33.449
  sigma_emb_scale = 0.01, axis_scale default     median   10.479   5%    2.936   95%    117.998

K = 128 embedding columns, gt = 8
  defaults: axis_scale = 0.08839 (1 / sqrt(K)), sigma_emb_scale = 1, c_base = 0.35355
  default scales, scales drawn (the check)       median   24.009   5%    5.261   95%    448.913
  median inside 1 to 2: NO
  scales held at their prior scales              median    1.542   5%    0.903   95%      2.657
  independent columns, scales held               median    1.706   5%    1.428   95%      2.052
  independent columns, scales drawn              median   23.191   5%    5.964   95%    423.155
  axis_scale for a median of 1.5: 0.00552 = 0.707 / K = 0.062 / sqrt(K)
  axis_scale = 0.00552, fresh draws              median    1.437   5%    0.305   95%     28.516
  sigma_emb_scale = 0.01, axis_scale default     median   14.377   5%    3.972   95%    166.916
```

The check itself, against the target of 1 to 2:

- 64 columns: median 16.3, 5% point 3.5, 95% point 326. Outside the target, about 11 times the middle of it.
- 128 columns: median 24.0, 5% point 5.3, 95% point 449. Outside the target, about 16 times the middle of it.

## What the numbers mean

**The design's arithmetic is right for scales held at their prior scales, and wrong once the scales are drawn.** With every axis scale held at `axis_scale` and the leftover scale held at 1, the median is 1.56 and 1.54 on the factor-model columns and 1.71 on independent columns, which is the 1.7 the design expects, and it does not change with the number of columns. With the scales drawn from their half-Cauchy priors it is 10 to 16 times larger.

**The cause is the tail of the half-Cauchy on the axis scales, not the correlation between columns.** Independent columns give almost the same medians (16.6 and 23.2). The contribution's variance is close to a sum of squared axis scales over the columns. A squared half-Cauchy draw has no mean, so that sum is carried by its few largest terms and grows with the square of the number of columns, not in proportion to it. Dividing the scale by the square root of the number of columns therefore leaves a contribution that still grows with that square root: doubling the columns from 64 to 128 raised the median by a factor of 1.47.

**So the default does not meet the design's aim** that the embedding's prior contribution should not grow with the number of embedding columns.

## The scale that would bring the median to 1.5

The contribution is exactly proportional to `axis_scale`, because every coefficient is its axis scale times something that does not depend on it. The value for a median of 1.5 follows by proportion and was confirmed with fresh draws.

- 64 columns: `axis_scale` 0.0115, where the default is 0.125 (10.9 times smaller). Fresh draws at that value: median 1.56.
- 128 columns: `axis_scale` 0.0055, where the default is 0.0884 (16.0 times smaller). Fresh draws at that value: median 1.44.
- Both are close to 0.72 divided by the number of columns (0.737 and 0.707). By proportion from the medians above, a default of 0.7 divided by the number of columns would give medians of 1.42 and 1.49. That rule was not itself run, and it rests on two column counts at 8 latent traits.

**The leftover scale is not a usable lever.** Setting `sigma_emb_scale` to 0.01, which all but removes the leftover, still leaves medians of 10.5 and 14.4, because the community row and the blends do not depend on it.

**A smaller axis scale fixes the median and not the spread.** At the values above the 5% and 95% points are 0.34 and 33 with 64 columns, and 0.31 and 29 with 128: one prior draw in twenty still puts the embedding's contribution above about 30 on the logit scale. That is the same half-Cauchy tail. The design names the regularised horseshoe, which caps large axis scales, as its fallback if the axis scales run away; whether to rescale `axis_scale`, adopt that fallback, or both is Doug's decision.

## Limits of the check

- It draws one species with no measured traits. A fit with traits has one more row of `G` per trait column, each with the axis scale as its prior standard deviation, so its prior contribution is larger still. That case was not simulated.
- One site matrix per setting and 2,000 draws. The fresh-draw medians of 1.56 and 1.44 at scales chosen to give 1.5 show a Monte Carlo error of about 5%.
- Only 8 latent traits were run. By design the base scale of `C` makes the blends' share of the prior variance the same for any number of them, so the result should carry over, but other numbers were not run.

## Fit tests added with this check

The tests are at the end of `tests/testthat/test-grouped-prior-fits.R` and all pass (run with `NOT_CRAN=true`; the whole file takes 67 seconds).

**Unused embedding columns are shrunk more than used ones.** The test simulates 150 sites and 30 species with 10 embedding columns, of which 3 carry an effect, and fits 2 chains of 1,500 + 1,500 iterations under the default scales (`axis_scale` 0.316). It asserts that the smallest posterior median axis scale among the 3 used columns is above the largest among the 7 unused ones. Three further fits of the same data, 49 seconds each, gave:

- used columns: 1.65, 0.15, 0.50; then 1.54, 0.14, 0.46; then 1.54, 0.14, 0.45;
- unused columns: between 0.006 and 0.083 in all three fits;
- smallest used over largest unused: 1.9, 2.2 and 1.7.

The test passed in every run made (the two runs of the test file and these three fits), but the margin on the weakest used column is about a factor of 2, and in one fit the two chains put that column's median at 0.08 and 0.19. The test fit is not seeded.

**Grouped fits run on the remaining model types and layouts.** Short fits return finite scales for the occupancy and two-stage model types, with 40 embedding columns on 30 sites, with no latent traits, on two parallel chains, and with a spatial field. To simulate the first two, `simulate_grouped_fixture()` in `tests/testthat/helper-grouped-prior.R` now fills in the detection parameters and replicate counts as `simulate_fixture()` does; its output for the binary model type is unchanged.

# Checks of the grouped covariate prior's default scales

Two runs of one prior check, both on 12 November 2026, for the [grouped covariate prior](../../superpowers/specs/2026-11-11-grouped-covariate-prior-design.md). The first run, under the default axis scale of the time (`1 / sqrt(K)` for K embedding columns), failed. Doug then changed the default to `0.7 / K`, keeping the half-Cauchy, and the check was run again. The second run is recorded first because it describes the package as it now is. The record of the first run follows it and is kept as written, apart from its headings, which are one level lower, and one marked correction. Later the same day the whole-branch review led to two changes to the sampler, to a chain's starting state and to a `chol()` call, and the fit tests were repeated; that record is the section *After the whole-branch review*, placed before the two runs because it supersedes their fit-test numbers. The prior check's own numbers are unaffected, since it fits nothing.

**Under the present default of `0.7 / K` the prior median is inside the target, and the upper tail is still wide.** The design asks for a prior median between 1 and 2 for the standard deviation across sites of the embedding's contribution to the linear predictor. The median is now **1.42 with 64 embedding columns** and **1.49 with 128**, both with 8 latent traits. The 95% points are 28.5 and 27.8.

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

The script calls the package's own `read_coef_prior()` for the default scales and `coef_prior_sds()` for how the scales combine, so it checks the prior as the package defines it. Its header comment names the function that defines every other part of the draw.

## After the whole-branch review: starting state, shrinkage test and the `chol()` warning

12 November 2026, after both runs below. R 4.6.1 with the Accelerate BLAS, on mains power with Low Power Mode off. The whole-branch review found that the sampler's starting state distorted the fits recorded under *Fit tests under the present default* below, and it settled what the `chol()` warning was. Both were fixed (commits `d2c19bb` and `885adf9`) and the fits were repeated. The numbers in the two runs below are kept as dated history.

### What was wrong with the starting state

A chain started `C`, the responses of the latent traits, at 1 on every covariate (on and above the diagonal of `C`, and 0 below it). The first update of the axis scales treats that starting `C` as a draw from its prior, whose standard deviation on an embedding column is the axis scale divided by the square root of the number of latent traits. A value of 1 is far too large for that prior, so the first update threw every axis scale to about the number of latent traits divided by the square root of the count of species, trait rows and latent traits together, whatever `axis_scale` was. By the review's measurement the chains had not come back after 2,000 sweeps, and because every chain was thrown the same way, R-hat could not show it.

The embedding columns of `C` now start at their prior scale: the starting axis scale divided by the square root of the number of latent traits. Measured columns, and every fit made without the grouped prior, start exactly as before.

Measured on a small fit of 60 sites and 12 species with 2 measured and 8 embedding columns, 4 latent traits and 2 chains, where the axis scales start at 0.0875 and the value they were thrown to is 0.918:

- before the change, the sixteen axis scales after the first sweep were between 0.80 and 1.11;
- after it, they are between 0.050 and 0.107.

### The shrinkage fit after the change

This is the fit of the shrinkage test: 150 sites, 30 species, 10 embedding columns of which 3 are used, seed 4, and 2 chains of 1,500 + 1,500 iterations. The posterior median axis scales are given for chain 1, for chain 2 and pooled over both. The same fit before the change is listed by chain under *Fit tests under the present default*.

- `X_emb.1`, used: 1.4331, 1.3566, pooled 1.3966;
- `X_emb.2`: 0.0387, 0.0283, pooled 0.0348;
- `X_emb.3`, used: 0.1155, 0.0889, pooled 0.1005;
- `X_emb.4`: 0.0196, 0.0428, pooled 0.0302;
- `X_emb.5`, used: 0.3856, 0.3837, pooled 0.3848;
- `X_emb.6`: 0.0062, 0.0323, pooled 0.0257;
- `X_emb.7`: 0.0404, 0.0448, pooled 0.0424;
- `X_emb.8`: 0.0251, 0.0284, pooled 0.0263;
- `X_emb.9`: 0.0326, 0.0185, pooled 0.0264;
- `X_emb.10`: 0.0302, 0.0377, pooled 0.0341.

What changed:

- The two strongly used columns are as clear as before. `X_emb.1` is 33 times the largest pooled median of an unused column (25 before) and `X_emb.5` is 9.1 times it (6.9 before). Within either chain alone, `X_emb.5` is at least 8.6 times the largest unused median of that chain.
- The weakly used column, `X_emb.3`, is now 2.4 times the largest unused pooled median, where it was 1.10. It is above every unused column in both chains: 0.1155 against a largest of 0.0404 in chain 1, and 0.0889 against 0.0448 in chain 2.
- The pooled medians of the unused columns are between 0.026 and 0.042, where they were between 0.030 and 0.058.
- The axis scales have still not converged at this chain length. R-hat on the log scale is 1.00, 1.14 and 1.00 for the three used columns and reaches 1.49 among the unused ones. Before the change the largest was 1.60.

Three other datasets of the same design, seeds 5, 6 and 7, were refitted too. Pooled medians after the change:

- seed 5: used 0.26, 1.54 and 0.47; unused between 0.015 and 0.214; weakest used over largest unused 1.22, where it was 1.31;
- seed 6: used 0.96, 0.44 and 1.82; unused between 0.026 and 0.052; ratio 8.4, where it was 4.66;
- seed 7: used 0.95, 1.08 and 0.43; unused between 0.010 and 0.091; ratio 4.7, where it was 4.28.

By chain these three fits still disagree about some unused columns: 0.0525 and 0.2792 for `X_emb.1` in seed 5 (R-hat 6.9), 0.1216 and 0.0121 for `X_emb.3` in seed 6, and 0.1760 and 0.0392 for `X_emb.4` in seed 7. In chain 1 of seed 5 the weakest used column (0.2167) is level with an unused one (0.2165). So the new starting state removed a transient that all chains shared. It did not make the axis scales of unused columns converge at 1,500 + 1,500 iterations, and a weakly used column is still not reliably told from an unused one at that length.

### What the fit tests now assert

The shrinkage test had been passing on the ratio of 1.10 for the weakly used column. It now asserts what the sampler delivers reliably, and is renamed "axes with a strong effect keep a much larger scale than axes with no effect": each of the two strongly used columns, `X_emb.1` and `X_emb.5`, must have a pooled median axis scale at least 3 times the largest pooled median of an unused column. The observed factors are 33 and 9.1, so the weaker of the two has a threefold margin. Nothing is asserted about `X_emb.3`. How well a weakly used axis is recovered is for the validation study to measure.

A new test, "a chain starts with C at its prior scale, so the first scale draw agrees with the start", makes the small fit described above and asserts that every axis scale after the first sweep is below 0.283, the geometric mean of the starting value and the value they used to be thrown to, and above the starting value divided by the same factor of 3.2. It failed before the change, with a largest value of 1.11.

Another new test checks the scale updates by simulation. With no data, a draw of the coefficients from their prior is alternated with the scale update and the draws of `G`, `C` and `A`, so that each scale must keep its half-Cauchy prior. Over 50,000 sweeps the share of draws below each prior quartile is within 0.018 of 0.25, 0.50 and 0.75 for all four scales. The tolerances are 0.045 to 0.09, five or more Monte Carlo standard deviations as measured over 24 seeds. Breaking the update on purpose (a wrong shape, a wrong rate, or the prior scale left out) made the test fail each time.

### The `chol()` warning, settled

The warning is described under *A warning from the sampler in the seed 5 fit* below. The review read Armadillo's source and established what it is:

- Before factorising, `chol()` runs a rough symmetry check. It compares two pairs of mirror-image elements, in the corners of the matrix, and warns when a pair differs by more than 10,000 times machine epsilon (about 2.2e-12), both absolutely and relative to the element.
- The matrix is formed as a general matrix product divided by the square of a leftover standard deviation. The two elements of a pair are the same sum computed in two orders, so they can differ by rounding. The recorded mismatch was 2.4e-11 of the element's value. It passed the relative tolerance because that sum happened to be close to zero, and the absolute tolerance because the division by a very small variance multiplied it by several thousand.
- The factorisation reads only the upper triangle of the matrix. One element of each pair never enters the draw, so the draw was never affected.

The fix (commit `885adf9`) wraps the argument of `chol()` in `sampleBuniv()` in `arma::symmatu()`, which copies the upper triangle onto the lower. The check then passes by construction for any finite matrix, and the factorisation is given the same upper triangle as before. Nothing else in the function reads the symmetrised matrix. What was checked:

- The seed 5 fit that printed the warning three times was repeated with the fix and with the old starting state put back in memory. It printed no warning and its `results_output` is identical to that of the fit that printed three.
- The four shrinkage fits (seeds 4 to 7) made with the new starting state printed the warning eight times between them before the fix. After it they print none, and all four `results_output` are identical before and after.
- 200 direct calls of `sampleBuniv()` on inputs of the kind recorded for seed 5 (30 species, a latent trait summing almost to zero, a leftover standard deviation of 0.0094): 154 printed the warning before the fix, none after, and the 200 draws are identical. Twenty of these are now a test.
- Fits made without the grouped prior are unchanged: the four reference fits of the first test in `tests/testthat/test-grouped-prior.R`, and eleven further cases compared with a fresh build of commit `4f2ee6c`, the last one before the grouped prior. Ten of the eleven are fits covering no traits, no latent traits, a factor and a spline, no covariates, two-stage and continuous data with a spatial field, parallel chains and a fixed intercept prior, and all ten `results_output` are identical. The eleventh, a one-species fit, stops with the same error in both builds, so it compares nothing.

The other `chol()` calls in the sampler were read. The joint coefficient draw, `sampleB_SoR()`, factorises a matrix that is made exactly symmetric before a diagonal prior precision is added, so it cannot trip the check. The remaining calls are in the draw of the site factor scores and on the detection side. They form their matrices without dividing by a small variance, the grouped prior does not touch them, and none printed the warning in any run made for this record, including the whole test suite, so they were left alone.

### Test results after these changes

With the code at commit `dcf7bf0`; the commit that follows it changes help text, comments and this record only.

- `tests/testthat/test-grouped-prior.R`: 172 expectations pass, none fail and none are skipped. Its first test, which compares four fits made without the grouped prior with the saved reference fits, passes 4 of 4.
- `tests/testthat/test-grouped-prior-fits.R`, run with `NOT_CRAN=true`: 186 expectations pass, none fail, with no warnings, in 67 seconds.
- The whole suite: 2,164 expectations pass, none fail, with no warnings, and 6 are skipped. Five of the skips need a declared single-threaded reference BLAS and one needs `OCCJSDM_SIMSTUDY=1`.

## Second run: the present default, `0.7 / K`

12 November 2026. R 4.6.1, on the `claude/grouped-covariate-prior` branch, at the commit that changes the default (Task 10 of the part 1 plan). The script is the same as in the first run except that the line reporting the defaults now prints the rule the package uses.

### Numbers

``` text
Prior check of the grouped covariate prior
R 4.6.1, 2000 prior draws, 250 sites, target median 1 to 2

K = 64 embedding columns, gt = 8
  defaults: axis_scale = 0.01094 (0.7 / K), sigma_emb_scale = 1, c_base = 0.35355
  default scales, scales drawn (the check)       median    1.425   5%    0.304   95%     28.489
  median inside 1 to 2: YES
  scales held at their prior scales              median    0.136   5%    0.078   95%      0.229
  independent columns, scales held               median    0.149   5%    0.122   95%      0.186
  independent columns, scales drawn              median    1.449   5%    0.381   95%     30.999
  axis_scale for a median of 1.5: 0.01152 = 0.737 / K = 0.092 / sqrt(K)
  axis_scale = 0.01152, fresh draws              median    1.557   5%    0.337   95%     33.449
  sigma_emb_scale = 0.01, axis_scale default     median    0.917   5%    0.257   95%     10.325

K = 128 embedding columns, gt = 8
  defaults: axis_scale = 0.00547 (0.7 / K), sigma_emb_scale = 1, c_base = 0.35355
  default scales, scales drawn (the check)       median    1.485   5%    0.326   95%     27.775
  median inside 1 to 2: YES
  scales held at their prior scales              median    0.095   5%    0.056   95%      0.164
  independent columns, scales held               median    0.106   5%    0.088   95%      0.127
  independent columns, scales drawn              median    1.435   5%    0.369   95%     26.181
  axis_scale for a median of 1.5: 0.00552 = 0.707 / K = 0.062 / sqrt(K)
  axis_scale = 0.00552, fresh draws              median    1.437   5%    0.305   95%     28.516
  sigma_emb_scale = 0.01, axis_scale default     median    0.890   5%    0.246   95%     10.327
```

The check itself, against the target of 1 to 2:

- 64 columns: median 1.42, 5% point 0.30, 95% point 28.5. The median is inside the target.
- 128 columns: median 1.49, 5% point 0.33, 95% point 27.8. The median is inside the target.

### What the numbers mean

**The median meets the target and no longer grows with the number of columns.** Doubling the columns from 64 to 128 moved the median from 1.42 to 1.49, a factor of 1.04, where under the old default it was a factor of 1.47.

**These are the first run's draws rescaled, not new draws.** The script uses the same seeds, and the contribution is exactly proportional to `axis_scale`, so every line that depends on the default is the first run's line multiplied by 0.7 divided by the square root of the number of columns (0.0875 at 64 columns, 0.0619 at 128). The medians are the 1.42 and 1.49 that the first run predicted by proportion. The run confirms that the package now applies the rule; it adds no information about Monte Carlo error, which the first run put at about 5%.

**The upper tail is as wide as before, relative to the median.** The 95% point is 20 times the median at 64 columns and 19 times at 128. One prior draw in twenty puts the standard deviation of the embedding's contribution above about 28 on the logit scale. This is the tail of the half-Cauchy on the axis scales, which the change of default does not touch. The validation study is to examine it.

**The diagnostic lines read differently now.** They were written to explain the failure of the first run, and each still runs:

- Scales held at their prior scales: 0.14 and 0.10 (0.15 and 0.11 on independent columns). The design's arithmetic, which holds the scales at their prior scales, no longer gives about 1.7. The median reaches the target only once the scales are drawn, which makes it 10 to 16 times larger, as in the first run.
- The scale that gives a median of 1.5: 0.737 and 0.707 divided by the number of columns, as in the first run. The default of 0.7 is 5% below the first and 1% below the second.
- A leftover scale of 0.01: medians of 0.92 and 0.89, which is about 60% of the medians under the defaults. The community row and the blends do not depend on the leftover scale, so it remains a weak lever.

The limits listed under the first run apply unchanged: one species with no measured traits, one site matrix per setting, 8 latent traits, and two column counts.

### Fit tests under the present default

Dated history: this subsection and the next record the fits as they were at commit `3df7253`, before the starting state was changed. The fits were repeated afterwards, and the current numbers and assertions are in the section *After the whole-branch review* above.

Run on mains power with Low Power Mode off. Final results, at commit `3df7253`:

- `tests/testthat/test-grouped-prior.R`: 124 expectations pass, none fail.
- `tests/testthat/test-grouped-prior-fits.R`, run with `NOT_CRAN=true`: 135 expectations pass, none fail.

**One test had the old default written into it and was updated.** The test "the starting scales leave measured columns at exactly 1" builds a prior with two embedding columns and expected their starting axis scales to be 1 divided by the square root of 2, which is 0.71. When the default was changed (commit `95ee99b`) that was the file's one failure, with 133 expectations passing, in 68 seconds: the axis scales now start at 0.7 divided by 2, which is 0.35, as they should, and the measured column still starts at exactly 1. Commit `3df7253` changed the expected value to `c(1, 0.7 / 2, 0.7 / 2)` and added one expectation that the starting axis scales equal the `axis_scale` of the prior the test builds. No fit-based test failed at either commit.

**Unused embedding columns are still shrunk more than used ones, by a thinner margin.** The test "axes with no effect are shrunk more than axes with an effect" passes. With 10 embedding columns its default `axis_scale` fell from 0.316 to 0.07. A fit that repeats the test's simulation and fit (49 seconds) gave these posterior median axis scales, pooled over the two chains:

- used columns: 1.45, 0.064 and 0.40;
- unused columns: between 0.030 and 0.058;
- smallest used over largest unused: 1.10, where the first run's three fits gave 1.7 to 2.2.

The same fit by chain, as the median in chain 1 and then in chain 2:

- `X_emb.1`, used: 1.4408 and 1.4588;
- `X_emb.2`: 0.0443 and 0.0269;
- `X_emb.3`, used: 0.0895 and 0.0317;
- `X_emb.4`: 0.0667 and 0.0517;
- `X_emb.5`, used: 0.3875 and 0.4159;
- `X_emb.6`: 0.0500 and 0.0238;
- `X_emb.7`: 0.0381 and 0.0572;
- `X_emb.8`: 0.0320 and 0.0424;
- `X_emb.9`: 0.0344 and 0.0265;
- `X_emb.10`: 0.0340 and 0.0503.

The two chains disagree about the weakly used column, `X_emb.3`: chain 1 puts it above every unused column (0.0895 against a largest of 0.0667) and chain 2 puts it below four of the seven (0.0317), so the axis scales had not converged at 1,500 + 1,500 iterations in this small test. The test passes on the pooled medians.

The fit is seeded through `simulate_grouped_fixture(seed = 4L)`, and a second run of it gave the same numbers to four decimal places, so the test's outcome is fixed on this machine. A margin of 10% is small, and a change to the order in which the sampler uses its random numbers could reverse it.

Three further fits used seeds 5, 6 and 7 in place of 4. The seed also sets the simulated data, so these are three other datasets of the same design, not repeats of one. On the medians pooled over both chains, all three put every used column above every unused one:

- seed 5: used 0.26, 1.48 and 0.54; unused between 0.029 and 0.196; ratio 1.31;
- seed 6: used 0.91, 0.42 and 1.77; unused between 0.023 and 0.090; ratio 4.66;
- seed 7: used 0.95, 1.03 and 0.41; unused between 0.019 and 0.095; ratio 4.28.

The seed 5 fit was also looked at by chain. In chain 1 one unused column (0.2747) is above the weakest used column (0.2660), and two unused columns differ between the chains by a factor of more than 10 (0.0045 and 0.0513; 0.0576 and 0.0032). Seeds 6 and 7 were not looked at by chain.

### A warning from the sampler in the seed 5 fit

The seed 5 fit printed `chol(): given matrix is not symmetric` three times and finished. The fit was repeated with the calls to the C++ coefficient draws wrapped so that the call that printed each warning could be recorded; the wrappers draw no random numbers and the repeat gave the same axis scales. What was observed:

- All three warnings come from `arma::chol` in `sampleBuniv()` in `src/jsdm.cpp` (line 528 at the time), called from `sample_GC()` in `R/jsdmfun.R`. That call draws one occupancy column's community response, trait effects and blend weights. The draw of the species' coefficients and the draw of the latent traits printed nothing.
- All three were in one sweep: chain 2, sweep 1,146 of 3,000, which is in the burn-in.
- They were for embedding columns 4, 7 and 8, the three with the smallest axis scales in that sweep: 0.032, 0.077 and 0.073. The embedding leftover scale was 0.294, so the leftover standard deviations of those columns, by whose squares the draw divides, were 0.0094, 0.0228 and 0.0215. The other seven embedding columns, with axis scales of 0.12 to 0.96, gave no warning in that sweep.
- In that sweep the second latent trait summed almost to zero across the 30 species, which made one off-diagonal element of the matrix small (9.8e-05 before the division). Recomputed in R, that element and its mirror image differ by about 2.4 parts in 100 billion of their value for all three columns.
- The draw was not affected. Armadillo's symmetry check compares pairs of mirror-image elements at a tolerance of 1e4 times machine epsilon, the factorisation reads only the upper triangle of the matrix, and the recorded relative mismatch of 2.4e-11 sits between an element that is read and one that is not. This was established by the whole-branch review, and the warning has since been silenced without changing any draw; the mechanism, the fix and the checks are under *The `chol()` warning, settled* above. This bullet first read that the effect on the draw had not been verified.
- The warning was not seen in the fits made without the grouped prior. Only four such fits were run in this work, the short reference fits in the first test of `tests/testthat/test-grouped-prior.R`. It was also not seen in the seed 4, 6 and 7 fits or in any run of the fits test file.
- Whether the old default gives the warning on the seed 5 data was not checked.

## First run: the old default, `1 / sqrt(K)`

12 November 2026. R 4.6.1, on the `claude/grouped-covariate-prior` branch at `3461213` (this check changes no package code). Task 8 of the part 1 plan for the [grouped covariate prior](../../superpowers/specs/2026-11-11-grouped-covariate-prior-design.md).

**The default axis scale fails the prior check.** The design asks for a prior median between 1 and 2 for the standard deviation across sites of the embedding's contribution to the linear predictor. With the default scales the median is **16.3 with 64 embedding columns** and **24.0 with 128**, both with 8 latent traits. No default was changed. The axis scale that brings the median to 1.5 is about 0.72 divided by the number of embedding columns, where the default is 1 divided by its square root.

### Numbers

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

### What the numbers mean

**The design's arithmetic is right for scales held at their prior scales, and wrong once the scales are drawn.** With every axis scale held at `axis_scale` and the leftover scale held at 1, the median is 1.56 and 1.54 on the factor-model columns and 1.71 on independent columns, which is the 1.7 the design expects, and it does not change with the number of columns. With the scales drawn from their half-Cauchy priors it is 10 to 16 times larger.

**The cause is the tail of the half-Cauchy on the axis scales, not the correlation between columns.** Independent columns give almost the same medians (16.6 and 23.2). The contribution's variance is close to a sum of squared axis scales over the columns. A squared half-Cauchy draw has no mean, so that sum is carried by its few largest terms and grows with the square of the number of columns, not in proportion to it. Dividing the scale by the square root of the number of columns therefore leaves a contribution that still grows with that square root: doubling the columns from 64 to 128 raised the median by a factor of 1.47.

**So the default does not meet the design's aim** that the embedding's prior contribution should not grow with the number of embedding columns.

### The scale that would bring the median to 1.5

The contribution is exactly proportional to `axis_scale`, because every coefficient is its axis scale times something that does not depend on it. The value for a median of 1.5 follows by proportion and was confirmed with fresh draws.

- 64 columns: `axis_scale` 0.0115, where the default is 0.125 (10.9 times smaller). Fresh draws at that value: median 1.56.
- 128 columns: `axis_scale` 0.0055, where the default is 0.0884 (16.0 times smaller). Fresh draws at that value: median 1.44.
- Both are close to 0.72 divided by the number of columns (0.737 and 0.707). By proportion from the medians above, a default of 0.7 divided by the number of columns would give medians of 1.42 and 1.49. That rule was not itself run, and it rests on two column counts at 8 latent traits.

**The leftover scale is not a usable lever.** Setting `sigma_emb_scale` to 0.01, which all but removes the leftover, still leaves medians of 10.5 and 14.4, because the community row and the blends do not depend on it.

**A smaller axis scale fixes the median and not the spread.** At the values above the 5% and 95% points are 0.34 and 33 with 64 columns, and 0.31 and 29 with 128: one prior draw in twenty still puts the embedding's contribution above about 30 on the logit scale. That is the same half-Cauchy tail. The design names the regularised horseshoe, which caps large axis scales, as its fallback if the axis scales run away; whether to rescale `axis_scale`, adopt that fallback, or both is Doug's decision.

### Limits of the check

- It draws one species with no measured traits. A fit with traits has one more row of `G` per trait column, each with the axis scale as its prior standard deviation, so its prior contribution is larger still. That case was not simulated.
- One site matrix per setting and 2,000 draws. The fresh-draw medians of 1.56 and 1.44 at scales chosen to give 1.5 show a Monte Carlo error of about 5%.
- Only 8 latent traits were run. By design the base scale of `C` makes the blends' share of the prior variance the same for any number of them, so the result should carry over, but other numbers were not run.

### Fit tests added with this check

The tests are at the end of `tests/testthat/test-grouped-prior-fits.R` and all pass (run with `NOT_CRAN=true`; the whole file takes 67 seconds).

**Unused embedding columns are shrunk more than used ones.** The test simulates 150 sites and 30 species with 10 embedding columns, of which 3 carry an effect, and fits 2 chains of 1,500 + 1,500 iterations under the default scales (`axis_scale` 0.316). It asserts that the smallest posterior median axis scale among the 3 used columns is above the largest among the 7 unused ones. Three further fits of the same data, 49 seconds each, gave:

- used columns: 1.65, 0.15, 0.50; then 1.54, 0.14, 0.46; then 1.54, 0.14, 0.45;
- unused columns: between 0.006 and 0.083 in all three fits;
- smallest used over largest unused: 1.9, 2.2 and 1.7.

The test passed in every run made (the two runs of the test file and these three fits), but the margin on the weakest used column is about a factor of 2, and in one fit the two chains put that column's median at 0.08 and 0.19.

Correction, made with the second run: this record first said that the test fit is not seeded. It is seeded indirectly. `simulate_grouped_fixture(seed = 4L)` sets R's seed, and the fit that follows draws its chain seeds from R's random numbers, so the test's own fit is the same every time. The three further fits above differ from one another, so they were not each made straight after that call.

**Grouped fits run on the remaining model types and layouts.** Short fits return finite scales for the occupancy and two-stage model types, with 40 embedding columns on 30 sites, with no latent traits, on two parallel chains, and with a spatial field. To simulate the first two, `simulate_grouped_fixture()` in `tests/testthat/helper-grouped-prior.R` now fills in the detection parameters and replicate counts as `simulate_fixture()` does; its output for the binary model type is unchanged.

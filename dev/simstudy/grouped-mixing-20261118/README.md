# Slow mixing of the grouped prior's scales

Started 18 November 2026. Nothing here has been run yet. The fitting-worker budget (five workers in total) was fully used by the AlphaEco simulation run when this was written.

## What was seen

The AlphaEco simulation run (package 32847ce, 2 chains of 5,000 + 5,000, 320 real sites, 48 simulated species, 64 embedding columns) found every grouped fit failing the design contract's convergence rule, in the first 15 fits. The report was from the harness session and its figures are interim and unverified here.

- 12 to 55 of the 64 log axis scales over R-hat 1.1 per fit (worst 1.3 to 5.0), minimum effective sample size 6 to 21 of 10,000.
- The embedding leftover scale over 1.1 in 11 of 15 fits, the measured-group scale in 7 of 15.
- 25 to 65 percent of the embedding B entries over 1.1 in the null and low-rank truths, 2 to 10 percent in the others.
- Fixed-prior arms converge (B R-hat about 1.01 to 1.03).

The package spec had already seen it: at 1,500 + 1,500 the axis scales of inactive columns did not converge (R-hat 1.49 to 7), and it warned that 5,000 + 5,000 might not be enough.

## What the code does

Each sweep draws all of B jointly given the current prior standard deviations, then draws sigma_b, sigma_emb and every axis scale lambda\_k from the leftovers of that B (`sample_coef_scales()` in `R/coef-prior.R`, called from `R/jsdmfun.R`), then G and C, then A. Each scale update is conjugate given the coefficients it scales.

## Hypotheses, to be tested before any fix

Nothing below is established. This is the systematic-debugging record: no fix until the root cause is shown.

- H1, the funnel. The coefficients of a column (its S leftovers, its G entries and its C entries) are drawn with prior standard deviation proportional to lambda\_k, and lambda\_k is then drawn from those same coefficients. This is the centred parameterisation. When the data say little about a column (the null truth, and the inactive columns of the low-rank truth) the coefficients simply follow lambda\_k and lambda\_k follows them, so log lambda\_k does a random walk with a step of about 1 over the square root of twice the number of coefficients tied to it (S plus the rows of G plus the rows of C, about 55). Predictions: the autocorrelation time of log lambda\_k for inactive columns is in the hundreds of sweeps; it grows with S; it is much shorter for columns that the data identify (active columns).
- H3, the ridge. The prior standard deviation of a leftover is sigma_emb times lambda\_k, so the two scales trade off, and sigma_emb is informed by all S times K leftovers at once, so its conditional step is tiny (about 1.3 percent at 48 by 64). G and C break the ridge only weakly. Predictions: the product sigma_emb times lambda\_k mixes much better than lambda\_k alone, and log sigma_emb is negatively correlated with the mean of log lambda\_k.
- H2, the tail. The half-Cauchy lets lambda\_k wander to very large values for unidentified columns, and B follows. Prediction: traces of large lambda\_k excursions that coincide with large B entries. A change in the prior (half-t) would help this and not H1 or H3.

## Plan

1. `repro.R` runs one synthetic fit (all data simulated here, no dataset content) and writes, per fit, the R-hat and autocorrelation time of log lambda\_k split by active and inactive columns, of the product, of log sigma_emb, the correlation of log sigma_emb with mean log lambda, and the share of embedding B entries over 1.1.
2. Runs, one worker each, at 1,500 + 1,500 first (the spec's small-test length, minutes per fit), then 5,000 + 5,000 for the one that decides: null and low-rank truths at S = 24, 48, 96 with n = 256 and K = 64. H1 predicts the null autocorrelation time rises with S and is much longer than the active columns'. H3 predicts the product's autocorrelation time is much shorter than lambda\_k's.
3. Only then candidate fixes, one at a time, each against the same runs: an interweaving or scale-expansion move (rescale a column's leftovers, G and C together with lambda\_k and accept on the likelihood ratio, which makes the coefficients follow the scale without changing the prior), a joint move of sigma_emb against all lambda\_k, and the half-t or regularised horseshoe fallbacks the spec names. If three fixes fail, stop and question the parameterisation with Doug.

## Rules

- Installed package only, never `devtools::load_all()`.
- One worker per run, and only while the five-worker budget has room. Doug approved synthetic reproductions in occJSDM on 18 November 2026.
- No dataset content. Results here are synthetic and may be committed as summaries.

## First results, 18 November 2026

Six synthetic fits, package 3b9d21a, 256 sites, 64 embedding columns, 5 measured covariates, 4 latent traits, 2 chains of 1,500 + 1,500, one seed each. Null and low-rank truths at 24, 48 and 96 species. One run per cell, so these are indications and not estimates. Tau is kept draws per effective draw, out of 3,000.

- The failure reproduces on simulated data, at every species count and in both truths: 47 to 60 of the 64 log axis scales over R-hat 1.1, 66 to 85 percent of the embedding B entries over 1.1, and the log of sigma_emb over 1.1 in all six (1.2 to 10). It is not specific to AlphaEco's sites.
- Each axis scale mixes slowly on its own. Median tau of an inactive column's log scale is about 130 to 240 sweeps in all six runs (a posterior sd of about 0.9 in log units, so a factor of 2.5), about 50 in the active columns of the low-rank runs (median scale 0.24 against 0.010 for the inactive ones). At 5,000 + 5,000 that is an effective sample size near 50 per column, which with two chains puts many columns over 1.1 by chance and a few far over it.
- The columns look independent. The mean pairwise correlation of log scales across columns is about 0, and the mean over all 64 columns has the same tau as a single column (160 to 300), so there is no common slow mode.
- Not supported: H3. The product sigma_emb times the axis scale does not mix better than the scale alone (median tau 150 to 220 against 130 to 240), and the correlation of log sigma_emb with the mean log scale changes sign between runs. sigma_emb is itself slower than the scales (tau 180 to 460, and 780 in one run).
- Not supported as stated: H1's step size prediction. The funnel argument says tau should grow with the number of coefficients tied to a scale, about 31 at 24 species and 103 at 96, which is a factor of about 3. The null runs give 158, 131 and 154 for 24, 48 and 96, and the low-rank runs 177, 130 and 220: flat to a factor of 1.3.
- Not supported: H2, the tail. The slowest columns are the ones with small scales (rank correlation of tau with the median scale from -0.2 to -0.46), not large ones. The slow part is the low end of the scale, where the data say nothing and the half-Cauchy prior is nearly flat in log scale.

What stands: independent slow random walks of log lambda\_k in columns the data do not identify, with a posterior that is wide in log scale. Why the walk is slow, and not faster with fewer coefficients tied to the scale, is not explained, so a fix chosen now would be a guess.

## Next

A minimal test of the centred coupling as the cause: a prototype scale-expansion move on each axis scale (rescale the column's leftovers, G and C together with the scale, accept on the likelihood ratio and the half-Cauchy prior; the scaled coefficients' prior and Jacobian cancel). If the coupling is the cause, tau should drop from about 150 sweeps to a few. If it does not, the coupling is not the cause and the half-Cauchy auxiliary-variable update itself is the next suspect. This is a change to `R/jsdmfun.R` and the likelihood evaluation, and is the first fix attempt.

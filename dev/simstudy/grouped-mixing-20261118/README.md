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

## Fix attempt 1: a scale-expansion move on each axis scale, 18 November 2026

Branch `claude/grouped-mixing`, cut from 3b9d21a. The prototype is commit eba6a20. Nothing was pushed. This is attempt 1 of at most 3. In short: the move fixes the axis scales and does not fix sigma_emb, so the decision rule below is met on its first condition and not on its second.

### Evidence added before the fix

Three measurements from the traces the six first runs saved (every fifth sweep, 8 columns per run). They change two of the readings under "First results".

- The funnel's step size is right. The change in a log axis scale over five sweeps has a standard deviation of 0.33, 0.27 and 0.21 at 24, 48 and 96 species in the null runs, and 0.33, 0.27 and 0.20 in the low-rank runs. H1 predicts the square root of 5 over the number of coefficients tied to a scale (31, 55 and 103), which is 0.40, 0.30 and 0.22. So the walk does take smaller steps with more species.
- The flat tau is most likely a ceiling of the estimator together with a narrower posterior, not evidence against H1. `coda::effectiveSize()` fits an autoregression within each chain, so it cannot see much more persistence than about a tenth of a 1,500-draw chain, which is where the values of 130 to 240 sit. In the same traces the lag-5 autocorrelation rises with species (0.905, 0.932 and 0.936 in the null runs, 0.869, 0.908 and 0.950 in the low-rank runs) and the within-chain sd of a log axis scale falls (0.78, 0.75 and 0.61 in the null runs). At 5,000 + 5,000 the same baseline gives a tau near 300 (below). The AlphaEco session's chain-length diagnostic, reported to this session on 18 November and not verified here, found that four times the draws only doubled the effective sample size of the axis scales. Read every baseline tau in this file as "at least".
- sigma_emb has the same funnel at a much larger size. Its step per sweep in log units is 0.026, 0.017 and 0.013 at 24, 48 and 96 species, which is 1 over the square root of species times columns (0.026, 0.018 and 0.013): it is drawn from every leftover of the embedding block at once, and those leftovers are then drawn given it. Its two chains sit up to 2.6 log units apart in the six runs, so crossing its range takes on the order of 10,000 sweeps. This also explains the H3 result: within one run sigma_emb hardly moves, so its product with an axis scale mixes like the axis scale.

### The move

After the scale draws and before the G, C and A draws, each embedding column gets one Metropolis proposal. A Normal step is added to the column's log axis scale, and the column's leftovers, its column of G and its column of C are multiplied by the same factor, so its row of B is too. The code is `axis_scale_move()` in `R/coef-prior.R`, called from `update_jSDMcoef()` in `R/jsdmfun.R`.

- The algebra in the brief was checked and holds. The three prior standard deviations from `coef_prior_sds()` are all proportional to the axis scale. With m coefficients scaled by a factor r, their prior density at the new scale is r to the power -m times the old one, and the Jacobian of the scaling is r to the power m, so the two cancel. What is left is the likelihood ratio times the ratio of the half-Cauchy density of the axis scale, each multiplied by the scale itself because the walk is in log scale.
- The likelihood is the one the coefficient draw uses, given the current Polya-Gamma weights, and not the logistic likelihood with the weights integrated out. The factor-score draw later in the same sweep still conditions on those weights, so a move that ignored them would not leave the joint distribution unchanged. Given the weights the log likelihood is quadratic in the scaling factor, so the ratio costs two small matrix products per column.
- The proposal step starts at 0.5 for every column and is adapted per column during burn-in only, toward an acceptance rate of 0.44, with the recursion `b0_shift_step()` uses. It is frozen for the kept sweeps. This was decided before any comparison fit was run: the prior-invariance run showed a fixed 0.5 accepted 88 percent of the time on a column the data do not inform, which is far too small a step. The brief allowed either a fixed step or burn-in adaptation.
- First version: grouped fits of model type "binary" with embedding columns only. Every other fit skips the move and draws the random numbers it drew before. A fit that ran the move records the steps and acceptance rates in `infos$axis_move`.

### Tests

`tests/testthat/test-axis-scale-move.R`, nine tests, written before the code and watched to fail.

- The acceptance decision against the target density written out term by term (Normal priors at the standard deviations `coef_prior_sds()` gives, the half-Cauchy prior in log scale, the likelihood given the weights) with the Jacobian of the scaling, over 20 states, replaying the move's random numbers. This is the check of the algebra.
- An accepted move multiplies the leftovers, G, C and the row of B by one factor; measured and rejected columns come back bit for bit.
- Prior invariance: with a constant likelihood, 50,000 sweeps of the scale draws, the move, and the G, C and A draws leave each scale half-Cauchy. Over 8 seeds the largest error of a quartile share was 0.016 for an axis scale. Leaving out the log-scale factor moves the shares by 0.47 to 0.49 and a prior scale wrong by a factor of 2 by 0.30 to 0.34 (`checks/invariance_mutations.R`).
- A spy inside a real fit: what the move returns for B and the axis scales is what the sweep returns, and the linear predictor it is given is the one implied by the coefficients it is given.
- The adaptation recursion, one step per column, the fit-level record, and the guard (off in an occupancy fit and under the fixed prior).

The full suite passes: 350 tests, 0 failures, 6 skips, none of them new. The fixed-prior reference test ran (the BLAS matches the fixture) and passes untouched. Between the two installed libraries, fits under the fixed prior for all four model types, grouped occupancy, two-stage and continuous fits, and a grouped binary fit with no embedding give identical draws; only the grouped binary fit with an embedding differs (`checks/identity_between_libraries.R`).

### The experiment

Baseline 3b9d21a and prototype eba6a20, each installed from a git export into its own library. Null and low-rank truths, 48 species, 256 sites, 64 embedding columns, 2 chains of 1,500 + 1,500, seeds 1 and 2, the same seeds in both arms. Then one seed of each truth at 5,000 + 5,000. `compare.R` makes the tables in `results/attempt1`.

The decision rule, fixed before the runs: the hypothesis is supported if the median tau of the inactive log axis scales falls at least 5 times, and log sigma_emb and the share of embedding B entries over R-hat 1.1 also improve, and the posterior summaries agree with the baseline. If tau does not fall 5 times the hypothesis is rejected and the move is not tuned to rescue it.

Results at 1,500 + 1,500, baseline then prototype.

- Null, seed 1. Tau of the inactive log axis scales 131 and 9.8. Log axis scales over R-hat 1.1: 57 and 1 of 64, largest 6.4 and 1.11. Log sigma_emb R-hat 9.8 and 1.86. Share of embedding B over 1.1: 0.87 and 0.20.
- Null, seed 2. Tau 175 and 10.3. Over 1.1: 50 and 3, largest 4.7 and 1.61. Log sigma_emb R-hat 1.45 and 7.3. B share 0.74 and 0.62.
- Low rank, seed 1. Tau of the inactive columns 130 and 12.8, of the active columns 83 and 21. Over 1.1: 59 and 19, largest 8.7 and 3.9. Log sigma_emb R-hat 10.1 and 1.12. B share 0.86 and 0.62.
- Low rank, seed 2. Tau of the inactive columns 166 and 17.8, of the active columns 67 and 33. Over 1.1: 41 and 24, largest 3.1 and 7.0. Log sigma_emb R-hat 1.70 and 19.3. B share 0.28 and 0.88.

Results at 5,000 + 5,000, seed 1, baseline then prototype.

- Null. Tau 292 and 11.5. Over 1.1: 60 and 0, largest 10.8 and 1.03. Log sigma_emb R-hat 3.5 and 1.14. B share 0.87 and 0.056.
- Low rank. Tau of the inactive columns 314 and 16.2, of the active columns 86 and 28. Over 1.1: 50 and 2, largest 8.5 and 2.2. Log sigma_emb R-hat 17.7 and 1.00. B share 0.78 and 0.35.

The move itself: acceptance after burn-in has a median of 0.44 to 0.45 in every fit. The adapted step has a median of 1.9 to 2.0 in log units in the null fits, four times the starting value, and in the low-rank fits goes down to 0.05 for the most strongly active columns. The tau of log sigma_emb is 350 to 830 in the prototype fits at both lengths, no better than the baseline's, and all of these are "at least".

### Does the answer stay the same

- In the fits above the two arms were compared block by block with the yardstick `compare.R` describes: the difference between the arm means against the differences between the two chains within each arm, where a ratio near 1 means the arms differ as much as their own chains do. In the null fits every block is at or below 1.03, except log sigma_emb at 5,000 + 5,000 (1.47, a single quantity). In the low-rank fits the blocks are between 0.66 and 1.30, except the embedding rows of B in seed 1, at 1.63 at the short length and 1.96 at the long one. Those are runs in which neither arm has converged and the chains sit in different configurations (below), which this yardstick cannot handle: when both chains of an arm happen to share a configuration it understates that arm's error. The embedding's contribution to the linear predictor, which is what predictions depend on, agrees: its mean spread across sites is 2.119 and 2.122 in the low-rank fit and 0.198 and 0.200 in the null fit at 5,000 + 5,000.
- A check with more power, on a problem small enough for both samplers to converge: a low-rank truth with 150 sites, 16 species and 8 embedding columns, 4 chains of 3,000 + 12,000, largest R-hat 1.10 or less in every run (`checks/agree_small.R`). The baseline was run twice on different chain seeds to see how much two runs of the same sampler differ. For the embedding rows of B the posterior means differ by at most 0.14 posterior sd between the two baseline runs and by at most 0.14 between baseline and prototype; for the measured rows 0.057 and 0.050; for the log axis scales 0.18 and 0.18. So the prototype sits inside the baseline's own run-to-run variation. On this problem the median effective sample size of a log axis scale rose from 330 to 981 of 48,000 draws.

### Cost

About 6 percent per sweep, so the move stays in R. Three alternating pairs of short fits (200 + 200) on a quiet machine, on mains power with Low Power Mode off: 31 s for the baseline and 33 s for the prototype in two rounds each (the third round of each arm took 56 to 57 s and is put down to scheduling). The same pairs under a load of four other workers gave medians of 80 s and 83 s. A direct timing of the move at this size gives about 3 ms per sweep, against about 39 ms for a whole sweep on the quiet machine.

Two wall times in the tables are not timings. The prototype fit of null seed 1 at 1,500 + 1,500 shows 4,620 s because the laptop slept on battery for about an hour during it. The 5,000 + 5,000 fits ran under different loads (the prototype pair with five workers busy, the baseline pair with two or three), so their seconds cannot be compared.

### Verdict against the rule

- Tau of the inactive log axis scales falls at least 5 times: met, in every pair. It fell 9 to 17 times at 1,500 + 1,500 and 19 to 25 times at 5,000 + 5,000. In the null truth the log axis scales now pass R-hat 1.1 at 5,000 + 5,000 (0 of 64 over it, largest 1.03).
- Log sigma_emb and the B share also improve: not met. Log sigma_emb R-hat improved in four of the six pairs and got worse in two (seed 2 of each truth), and its tau shows no consistent change (higher in four pairs, lower in two). With two chains on a walk this slow, R-hat mostly records where the chains happen to be. The B share follows it: where the two sigma_emb chains agree the share falls a long way (0.87 to 0.056 in the null fit at 5,000 + 5,000), and where they do not it stays high.
- Posterior summaries agree: met, with the limits given above.

So the hypothesis is supported for what it was about: the centred coupling between an axis scale and its coefficients is why the axis scales walk slowly, and the move removes it. It is not a cure for the fits, because sigma_emb still mixes slowly. The move was not tuned after the comparison.

### What remains

- sigma_emb, in both truths. It is the same centred walk with a step of 1 over the square root of species times columns.
- The low-rank truth has two configurations, and chains do not move between them in 10,000 sweeps. In one, log sigma_emb is about 1.7 to 2.0, the active columns' axis scales are small (mean log scale about -2.5 to -3.3) and the leftovers carry the signal; the fitted active rows of B are shrunk (rms 0.56 to 0.61 against a truth of 0.86 to 0.97). In the other, log sigma_emb is about -0.9 to -1.4, the active columns' axis scales are larger (about -0.9 to -1.9), the leftovers are small and the community response and blends carry the signal; the active rows are closer to the truth (0.64 to 0.83). Of the six low-rank baseline chains, four are in the first configuration; of the six prototype chains, one is. Whether the first is a real mode or a start-up trap that the slow sigma_emb walk cannot leave is not known. Chains that share data and seeds across the two run lengths are not independent, so these counts are small.
- Within the second configuration a few axis scales are still slow in the low-rank truth (2 of 64 over 1.1 at 5,000 + 5,000, largest 2.2, with the two sigma_emb chains in agreement). An untested reading is slow switching of which of several correlated columns carries a signal, which the AlphaEco session's diagnostic also suggests.

### Proposed next, not implemented

Fix attempt 2 would be the same kind of move for sigma_emb: propose a Normal step in log sigma_emb and multiply every leftover of the embedding block by the same factor. The leftovers' prior and the Jacobian cancel as before, and the acceptance is the likelihood ratio given the Polya-Gamma weights (quadratic again, one sites-by-species product per proposal) times the half-Cauchy ratio of sigma_emb in log scale. Where the leftovers are small, as in the null truth and in the second low-rank configuration, the likelihood is nearly flat in that direction and sigma_emb should move freely.

This is not the joint move of sigma_emb against all the axis scales that the brief names, and the reason is a reading of the evidence, not a test. The joint move travels along the ridge on which the leftover standard deviation stays constant. The slow direction measured here is the level of the leftover standard deviation itself: the product of sigma_emb with an axis scale mixes no better than its parts, and the step is set by the number of leftovers. Along the ridge the 64 half-Cauchy priors of the axis scales hold the position fairly tightly, so there is less to gain. If attempt 2 leaves the two low-rank configurations unconnected, the joint move is the candidate for attempt 3, and after that the parameterisation should be questioned with Doug as the plan says.

A change of prior (half-t or regularised horseshoe) is not indicated by these runs: the slow columns are the small-scale ones, and a move that leaves the prior alone removed the slowness of the axis scales.

### Not verified

- One seed at 5,000 + 5,000 and two at 1,500 + 1,500, at one size (48 species, 256 sites, 64 columns). No other species counts were rerun with the move.
- Which low-rank configuration holds the posterior mass.
- The move with a spatial field. Such a fit runs and returns finite scales (an existing test), and the linear predictor handed to the move includes the spatial term by construction, but the wiring test covers a non-spatial fit only.
- Occupancy, two-stage and continuous fits, where the move is off. It would be valid there given the latent states, but that was not built or tested.
- `R CMD check` was not run. The package documentation, NEWS and defaults are unchanged; the grouped prior's help page still says the axis scales mix slowly.
- Nothing was run on the AlphaEco data. The AlphaEco session's figures quoted above are that session's report.

A correction to the six first results: `repro.R` took the embedding rows of B to be every row after the fifth, which included one measured row (a simulated measured covariate is a factor with two contrast columns). Their two B columns counted 65 rows, one of them measured. The script now takes the rows from the fit's own labels, and every figure in this section uses the 64 embedding rows.

## Fix attempt 2: a scale-expansion move on sigma_emb, 18 November 2026

Doug approved building this on 18 November 2026, after the report on attempt 1. It is attempt 2 of at most 3. This part was written and committed before any code or fit of attempt 2 existed.

### Hypothesis and decision rule, fixed in advance

The hypothesis: sigma_emb mixes slowly for the same reason the axis scales did. It is drawn given every leftover of the embedding block, and those leftovers are then drawn given it, so its log moves by about 1 over the square root of species times columns per sweep.

The test: one Metropolis proposal per sweep, after the scale draws and before the axis scale moves. A Normal step is added to log sigma_emb and every leftover of the embedding block is multiplied by the same factor, so each embedding row of B moves by that factor times its leftover; G, C, A and the axis scales stay where they are. The leftovers' prior and the Jacobian of the scaling cancel, as in attempt 1, which leaves the likelihood ratio given the Polya-Gamma weights times the half-Cauchy ratio of sigma_emb in log scale. One proposal per sweep, with a step that starts at 0.5 and is adapted during burn-in only, exactly as in attempt 1. Same scope: grouped fits of model type "binary" with embedding columns. The move is added on top of attempt 1's, which stays as it is.

The comparison: attempt 1's prototype (eba6a20) is the baseline, so the only difference between the arms is the new move. The same four runs at 1,500 + 1,500 (null and low-rank truths, seeds 1 and 2, 48 species, 256 sites, 64 columns), then seed 1 of each truth at 5,000 + 5,000. The baseline arm reuses attempt 1's prototype fits, which were made with the same library, script and seeds.

The rule. The hypothesis is supported if all of these hold.

- The tau of log sigma_emb falls at least 5 times in every pair. The baseline's values (350 to 830) are lower bounds, so this is the easy direction for the test to fail in and the hard one to pass by accident.
- In the null-truth fits, log sigma_emb is under R-hat 1.1 at 5,000 + 5,000 and the share of embedding B entries over 1.1 is lower than the baseline's in every null pair.
- The axis scales are not made worse: the median tau of the inactive log axis scales is no more than 1.5 times the baseline's in every pair.
- The posterior agrees with the baseline's on the small problem of attempt 1, judged against the difference between two baseline runs.

If the tau of log sigma_emb does not fall 5 times, the hypothesis is rejected and the move is not tuned to rescue it. The low-rank truth is reported separately and is not part of the rule, because of its two configurations: what is recorded there is which configuration each chain ends in, and whether any chain crosses.

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

### What was built

Prototype commit 4437375, on top of attempt 1. `sigma_emb_move()` in `R/coef-prior.R`, called from `update_jSDMcoef()` just before `axis_scale_move()`, which is handed the B and linear predictor the new move leaves. The algebra was checked as for attempt 1 and holds: a leftover's prior standard deviation is sigma_emb times its column's axis scale, so scaling sigma_emb and every embedding leftover by one factor leaves the leftovers' prior times the Jacobian unchanged. The new move shares attempt 1's state, guard and burn-in adaptation, and `infos$axis_move` gains its step and acceptance rates.

Tests, written first and watched to fail: six in `tests/testthat/test-sigma-emb-move.R`. The acceptance decision against the written-out target density over 40 states; the joint scaling of every embedding leftover; prior invariance with both moves in the sweep (largest quartile error over 8 seeds 0.013 for sigma_emb and 0.030 for an axis scale; leaving out the log-scale factor moves the sigma_emb shares by 0.39 or makes the run fail, and a prior scale wrong by a factor of 2 moves them by 0.20, `checks/invariance_mutations_sigma_emb.R`); the wiring and order inside a real fit; the record. The full suite passes: 356 tests, 0 failures, 6 skips, none new. Fits under the fixed prior and grouped fits with the moves off are bit-identical to the baseline and to attempt 1.

One existing test was changed. "a chain starts with C at its prior scale" in `test-grouped-prior-fits.R` bounded the first stored axis scales within a factor of 3.2 of their start, on the reasoning that the first stored value is the first Gibbs draw. With a scale move in the sweep that is no longer so, and a first move of that size is not rare. The test passed under attempt 1 by its seed and failed under attempt 2 when the random stream shifted. It now reads the first Gibbs draw itself from `sample_coef_scales()`, with its bounds unchanged, and passes.

### Results

Attempt 1 then attempt 2, at 1,500 + 1,500.

- Null, seed 1. Tau of log sigma_emb 467 and 33, R-hat 1.86 and 1.05. Log axis scales over 1.1: 1 and 3 of 64, largest 1.11 and 1.20. Tau of the inactive log axis scales 9.8 and 9.8. Share of embedding B over 1.1: 0.20 and 0.38.
- Null, seed 2. Tau of log sigma_emb 411 and 32, R-hat 7.3 and 1.08. Over 1.1: 3 and 3, largest 1.61 and 1.45. Tau of the inactive log axis scales 10.3 and 10.7. B share 0.62 and 0.47.
- Low rank, seed 1. Tau of log sigma_emb 828 and 131, R-hat 1.12 and 1.07. Over 1.1: 19 and 6, largest 3.9 and 1.37. Tau of the log axis scales 12.8 and 12.7 for inactive columns, 21 and 19 for active ones. B share 0.62 and 0.54. All four chains are in the second configuration.
- Low rank, seed 2. Tau of log sigma_emb 348 and 370, R-hat 19.3 and 1.08. Over 1.1: 24 and 6, largest 7.0 and 1.48. Tau of the log axis scales 17.8 and 27.1 for inactive columns, 33 and 75 for active ones. B share 0.88 and 0.025. Attempt 1 has one chain in each configuration; attempt 2 has both in the first.

At 5,000 + 5,000, seed 1.

- Null. Tau of log sigma_emb 820 and 26, R-hat 1.14 and 1.02. Log axis scales over 1.1: 0 and 0, largest 1.03 and 1.05. Tau of the inactive log axis scales 11.5 and 11.5. B share 0.056 and 0.088.
- Low rank. Tau of log sigma_emb 648 and 180, R-hat 1.00 and 1.02. Over 1.1: 2 and 3, largest 2.2 and 4.7. Tau of the log axis scales 16.2 and 15.1 for inactive columns, 28 and 31 for active ones. B share 0.35 and 0.39. All four chains are in the second configuration.

The move itself: acceptance after burn-in is 0.30 to 0.51. The adapted step is 1.0 to 1.2 in log units in the null fits, 0.10 to 0.18 in the low-rank chains in the second configuration, and 0.03 in the two chains in the first.

### Verdict against the rule

The rule is not met. The move was not tuned.

The rule as written above is ambiguous on one point, and that is recorded here and not resolved in the move's favour. Its first and third conditions say "in every pair", and its last paragraph says the low-rank truth is not part of the rule. What was meant was that the configurations would be reported separately. Read with every pair, the first condition fails and the hypothesis as stated is rejected. Read with the low-rank truth left out, the first condition is met and the rule still fails, on the share of embedding B in the null pairs. The rest of this section takes the stricter reading.

- Tau of log sigma_emb falls at least 5 times in every pair: not met. It fell 13 to 14 times in the null fits at the short length and 31 times at the long one. In the low-rank truth it fell 6.3 times in seed 1 at the short length, 3.6 times at the long one, and not at all in seed 2.
- Null truth, log sigma_emb under 1.1 at 5,000 + 5,000: met (1.02). Share of embedding B lower in every null pair: not met (higher in two of the three pairs).
- Axis scales not made worse: met in five pairs, and missed by a hair in low-rank seed 2 (27.1 against a limit of 26.7), where the two arms have different mixes of configurations.
- Posterior agrees on the small problem: met. Against the first baseline run the prototype's posterior means differ by at most 0.12 posterior sd for the embedding rows of B, 0.04 for the measured rows, 0.22 for the log axis scales and 0.11 for log sigma_emb; the two baseline runs differ from each other by 0.14, 0.06, 0.18 and 0.13. Against the second baseline run every one of the prototype's differences is smaller than those.

What the results say is narrower than the hypothesis. Where the embedding leftovers carry no signal, sigma_emb was stuck for the reason proposed and the move frees it: the likelihood is nearly flat when the leftovers are scaled, so the step is about 1. Where the leftovers carry signal the likelihood pins them, the step that keeps a 0.44 acceptance rate falls to between 0.03 and 0.18, and the move can do little. In the first low-rank configuration that step is about the size of the Gibbs draw's own (0.018), and nothing improves. There the product of sigma_emb with an axis scale is fixed by the data, and what is left slow is the trade-off between the two along that ridge. That is the joint move the original brief named, and the reasoning given under attempt 1 for setting it aside was wrong for this regime.

In the large fits the yardstick of `compare.R` gives log sigma_emb a ratio of 3.9 and 7.3 in two null pairs. That is the baseline arm's doing: attempt 1's two chains start at the same value and had hardly left it (tau 470 to 820), so they agree with each other and not with the posterior. Attempt 2's six null chains have mean log sigma_emb between -1.05 and -1.56.

### The R-hat of the embedding coefficients

The share of embedding B entries over 1.1 does not come down with the scales, so the diagnostic was checked against two others on the same draws (`checks/rhat_check.R`, attempt 2's library, seed 1). The package uses coda's point estimate, which carries a degrees-of-freedom correction that grows when the draws are heavy-tailed, and an embedding coefficient's draws are (median kurtosis about 12, against 3 for a Normal). The other two are the plain ratio without that correction and the rank-normalised split R-hat, taking the larger of its bulk and folded versions.

- Null truth, 1,500 + 1,500. Share of embedding B over 1.1: coda 0.38, plain 0.12, rank-normalised 0.45.
- Null truth, 5,000 + 5,000. Coda 0.088, plain 0.004, rank-normalised 0.035 (largest 1.21). Measured B, the log axis scales, sigma_emb and the embedding's contribution per species have nothing over 1.1 on any of the three.
- Low-rank truth, 5,000 + 5,000, both chains in the second configuration. Coda 0.39, plain 0.084, rank-normalised 0.35. Measured B and the embedding's contribution have nothing over 1.1 on any of the three (largest 1.07). Two axis scales are over 1.1 on the rank-normalised version (largest 1.73), and their columns' coefficients are over 1.1 in 73 and 100 percent of entries; in the other 62 columns the share is 0.34.

So coda's correction inflates the share several times over, but the coefficients are not simply fine: the rank-normalised version, which also looks at spread, is nearly as high as coda's in the low-rank truth. A null fit at 5,000 + 5,000 is close to clean on it. In the low-rank truth about a third of the individual embedding coefficients still differ between chains while their combined effect on the linear predictor and the measured coefficients do not, which is what trade-offs among correlated columns would look like. That reading has not been tested.

### The two low-rank configurations

No chain has crossed. In the 18 low-rank chains run so far (the baseline and both attempts, both run lengths), log sigma_emb stays between 1.1 and 2.4 in seven chains and between -2.6 and 0.0 in eleven, over whole runs of up to 10,000 sweeps. By arm, the first configuration holds four of six baseline chains, one of six with attempt 1 and two of six with attempt 2. Chains that share data and seeds across the two run lengths are not independent.

### Cost

Both moves together add about 5 percent per sweep. Three alternating rounds of short fits (200 + 200), run alone on battery with Low Power Mode off: 31, 32 and 31 s for the baseline, 32, 33 and 32 s for attempt 1, and 33, 34 and 32 s for attempt 2. The wall times in attempt 2's tables are not timings: the fits started under Low Power Mode and the laptop slept several times.

### Where this leaves the problem

- The axis scales: fixed by attempt 1 where the data do not inform a column, and much improved elsewhere.
- sigma_emb in the null truth: fixed by attempt 2. With both moves a null fit at 5,000 + 5,000 has every scale, every measured coefficient and the embedding's contribution under R-hat 1.1.
- The low-rank truth: not fixed, and not a scale that walks slowly. It has two configurations that no chain crosses, a few axis scales that disagree between chains inside one configuration, and individual embedding coefficients that trade off.

Two of the three attempts are used. A third move of the same kind (the joint move of sigma_emb against the axis scales) would be expected to speed sigma_emb inside the first configuration. It is not expected to connect the configurations, because moving between them means the leftovers handing the signal to the community response and blends, or back, and no change of scale does that. So the recommendation is to take the parameterisation question to Doug now, with one attempt in hand, and not to spend the third attempt on that move first. The question is that the leftovers and the blends can each explain the same coefficients, and the prior does not choose firmly between them. Things to decide between, none of them tried: which configuration the model is meant to find and whether chains should start in it; whether the half-Cauchy prior on sigma_emb, which allows values near 7, is wider than intended; and whether the validation's convergence rule should be judged on the scales at all in a truth like this one.

### Not verified, attempt 2

- One seed at 5,000 + 5,000, two at 1,500 + 1,500, one size.
- Which low-rank configuration holds the posterior mass, and whether the first is reached only through the start.
- The R-hat comparison was run on one seed of each truth.
- A spatial field with the new move, other model types, and `R CMD check`, as for attempt 1. The package documentation, NEWS and defaults are unchanged.
- Nothing was run on the AlphaEco data.

## The prior scale of sigma_emb: a settings-only experiment, 18 November 2026

Not a fix attempt, and no package code is changed. Doug approved running it on 18 November 2026 after a discussion of the parameterisation. This part was written before the fits were run.

sigma_emb is one number for the whole embedding block: the ratio of the size of a column's leftovers to the size of its shared part (community average, traits and blends). The two low-rank configurations are two values of that ratio, about 0.3 (blends carry the signal) and about 6 (leftovers carry it, and nothing is pooled across species). Scored against the simulated truth, the 18 low-rank chains run so far recover the active embedding coefficients better in the blend configuration: the slope of the estimate on the truth is 0.69 to 0.75 there and 0.40 to 0.55 in the leftover configuration, with the measured coefficients about equally well recovered in both. Which configuration a chain reaches depends mostly on the data set: seed 1 goes to blends in 10 of 12 chains, seed 2 to leftovers in 5 of 6.

The experiment asks whether a smaller half-Cauchy scale on sigma_emb (0.3 in place of the default 1, set through `sigma_emb_scale`) puts every chain in the blend configuration, and what it costs. Four data sets: the low-rank truth with seeds 1 and 2, the null truth with seed 1, and an "unstructured" truth in which 10 embedding columns have coefficients that are independent across species, with no community average and no blends, which is the case a small ratio should fit worst. Each is fitted at both scales with attempt 2's library (4437375), 4 chains of 1,500 + 1,500. `prior_scale_experiment.R` writes the tables to `results/prior-scale`.

What will be read off: how many chains end in each configuration; the recovery of the active embedding coefficients, the measured coefficients and the embedding's contribution to the linear predictor; and the R-hat of the scales, the measured coefficients and the embedding's contribution. `sigma_emb_scale` is also the chain's starting value of sigma_emb, so a clean result would not separate a better prior from a better start. That would need a second, start-only run.

### Results

Doug settled the question of intent on 18 November 2026, while the fits were running: the purpose of the grouped prior is to learn habitat gradients when they exist. So in a low-rank truth the leftover configuration is an unwanted solution and not an alternative answer.

The first eight fits looked like support for the smaller scale. In low-rank seed 2 the default put 1 of 4 chains in the blend configuration and a scale of 0.3 put 4 of 4 there, with the slope of the active coefficients on the truth rising from 0.42 to 0.62; seed 1 was 4 of 4 at both scales; and nothing was lost in the null truth or in the unstructured truth, where sigma_emb went to about 8 under either scale and the active coefficients were recovered identically (slope 0.85 and 0.84). Three more runs were then added to separate the prior from the start, to add a third data set and to check the null truth at the longer length. `--start` in the script sets the starting value of sigma_emb apart from the prior scale, inside the script's R session only. They change the reading.

Low-rank chains in the blend configuration, of 4, by prior scale and starting value of sigma_emb.

- Seed 1. Scale 1, start 1: 4. Scale 0.3, start 0.3: 4.
- Seed 2. Scale 1, start 1: 1. Scale 1, start 0.3: 3. Scale 0.3, start 1: 0. Scale 0.3, start 0.3: 4.
- Seed 3. Scale 1, start 1: 0. Scale 0.3, start 0.3: 2.

What this says.

- The starting value decides, not the prior scale. In seed 2 the smaller scale with the default start puts no chain in the blend configuration, and the default scale with the smaller start puts three there. A half-Cauchy is too heavy-tailed for its scale to matter at values near 7: moving the scale from 1 to 0.3 lowers the prior density there by a factor of about 3.
- A smaller start is not reliable either. In seed 3 it rescues two chains of four, and those two do not agree with each other: log sigma_emb is about -3.6 in both, but the slope of the active coefficients is 0.57 in one and 0.83 in the other, the embedding's contribution is over R-hat 1.1 for 75 percent of species, and 12.5 percent of the measured coefficients are too (largest 1.34). So the blend configuration is itself more than one place.
- A low start does not trap chains in the blend configuration. In the unstructured truth chains started at 0.3 moved to about 8 within burn-in.
- The leftover configuration is not rare at the default: 3 of 4 chains in seed 2 and 4 of 4 in seed 3, with the active coefficients at half their true size (slope 0.38 to 0.50). In those fits the scales and coefficients agree across chains, so R-hat reports nothing wrong: in seed 3 at the default, log sigma_emb has a rank-normalised R-hat of 1.10 and only 3 axis scales are over 1.1.
- The smaller scale has a cost in the null truth that the first eight fits did not show. At 5,000 + 5,000 with 2 chains, 26 percent of the embedding coefficients are over 1.1 on the rank-normalised R-hat at a scale of 0.3, against 3.5 percent at the default. The scales, the measured coefficients and the embedding's contribution are under 1.1 at both (largest 1.04).

So this experiment does not support changing `sigma_emb_scale`, and the default is left alone. What it does establish is that the configuration is chosen during burn-in by where the chain starts and is then kept, which is a posterior with several separated regions that the sampler does not cross, and that the unwanted region is the usual destination from the default start in two of three data sets. Which region holds the posterior mass is still not known.

Not verified: three low-rank data sets, one size, 4 chains of 1,500 + 1,500; whether the recovered blends match the true ones (only coefficients and their combined effect were scored); held-out predictions.

## Fix attempt 3: a staged burn-in that holds sigma_emb small at first, 18 November 2026

Doug approved building this on 18 November 2026, after the prior-scale experiment. It is the third and last attempt before the parameterisation itself is questioned. This part was written and committed before any code or fit of attempt 3 existed.

### Hypothesis and decision rule, fixed in advance

The hypothesis: chains reach the leftover configuration because of what happens early in burn-in. From the default start the leftovers can take up the signal within a few sweeps, while the community response and blends need many sweeps to organise, and once the leftovers have it the chain does not come back. If the leftovers are kept small until the blends have organised, a low-rank truth will end in the blend configuration and stay there, and a truth with no shared structure will still move to the leftover configuration once it is free to.

The test: for the first half of burn-in (the whole number of sweeps at or below half), sigma_emb is held at one tenth of its prior scale, which is 0.1 at the default. The chain starts there, the draw of sigma_emb is skipped and so is its scale move; everything else runs as in attempt 2, including the draws of the axis scales, which use the held value. From the next sweep on the sampler is attempt 2's, unchanged. Burn-in is discarded, so the kept draws come from the same sampler as before and its stationary distribution is not touched. Same scope as the two moves: grouped fits of model type "binary" with embedding columns. The length of the hold and its value are fixed here and will not be adjusted after the fits.

The comparison: attempt 2's library (4437375) is the baseline, and its fits in `results/prior-scale` at the default scale and start are reused. With attempt 3: the low-rank truth with seeds 1, 2 and 3, the unstructured truth and the null truth, each with 4 chains of 1,500 + 1,500 at the default prior; then low-rank seeds 2 and 3 with 4 chains and the null truth with 2 chains at 5,000 + 5,000. The tables are made by `prior_scale_experiment.R` and go to `results/attempt3`.

The rule. The hypothesis is supported if all of these hold.

- Every one of the 12 low-rank chains at 1,500 + 1,500 and every one of the 8 at 5,000 + 5,000 is in the blend configuration (mean log sigma_emb below 0.5 over the kept draws), and none of them has a kept draw of log sigma_emb above 0.5.
- In low-rank seeds 2 and 3 the slope of the active coefficients on the truth, pooled over chains, is higher than the baseline's (0.42 and 0.50). In seed 1 it is within 0.02 of the baseline's 0.78.
- In the unstructured truth all 4 chains are in the leftover configuration over the kept draws, and the pooled slope is within 0.02 of the baseline's 0.85.
- In the null truth at 5,000 + 5,000 no log axis scale and not log sigma_emb is over 1.1 on the rank-normalised R-hat, and the error of the embedding's contribution is no more than 1.1 times the baseline's at 1,500 + 1,500 (0.094).

If any low-rank chain ends in the leftover configuration the hypothesis is rejected, and the hold is not lengthened or lowered to rescue it. Reported but not part of the rule: how well the chains inside the blend configuration agree with each other (the R-hat of the measured coefficients and of the embedding's contribution), since seed 3 showed that they may not.

## The rotation test, 19 November 2026

Not a fix attempt, and no package code is changed. Doug asked for it on 19 November 2026. This part was written before any of its fits were run. The results of attempt 3 are written up after it, below.

An embedding's columns have no meaning one by one. Multiplying them by an orthogonal matrix gives another embedding that encodes the same thing, and the model can represent the same truth exactly in either. The grouped prior, though, gives every column a heavy-tailed scale of its own, so it asks which individual columns matter. The test asks whether that matters: is the answer the same when the same data are fitted with the columns rotated?

The design. The three low-rank data sets (seeds 1, 2 and 3), each fitted with the embedding columns as simulated and with them rotated by two random orthogonal matrices, the same two for every data set. Sites, species and outcomes are identical; only the columns handed to the fit differ, and the fit standardises each column as it always does. Attempt 3's library (a3d5bb7), 4 chains of 5,000 + 5,000, the same chain seeds for every rotation of a data set. `rotation_test.R` writes the tables to `results/rotation`. Everything scored is something a rotation does not change: the embedding's contribution to the linear predictor at each site and species, against the truth; the slope of that contribution on the true one; the measured coefficients; and how well the four chains agree on these.

What counts as a material change, fixed in advance. For a data set, the rotated fits differ materially from the unrotated one if, in both rotations, the error of the embedding's contribution is more than 10 percent larger or the slope of the contribution on the truth is more than 0.05 lower. If that holds in at least two of the three data sets, the prior depends on how the embedding happens to be rotated in a way that matters, and an option that does not treat columns one by one moves to the front. If it holds in none, the per-column scales do no harm of this kind.

Two things to keep in mind when reading the result. The simulated truth is sparse in the columns as simulated (10 of 64 carry signal), which is the case the per-column scales are built for, so some loss after rotation is expected and by itself shows only that the prior uses the sparsity it is given. A real embedding has no such privileged set of columns, so the rotated fits are the more realistic of the two. The second reading is of the chains: if they agree better in one set of columns than the other, that says whether the choice among columns is what makes the blends hard to settle. That comparison is reported and is not part of the rule.

# Refit at the community-mean intercept prior: what changes in the lessons

4 November 2026. Since 1b4a2cd, `runOccJSDM()` defaults to `b0_prior = "community"` for binary, occupancy and two-stage data with more than one species. Every shipped and lesson fit had been made with the fixed Normal(0, 1) prior. This note records the refit at the new default and lists every number or claim in the vignettes and README that changes, for Doug's review of Lessons 2 to 4. No lesson or Quickstart prose was edited; the `.md` renders of Lessons 0, 2, 3 and 4 were regenerated.

## Summary

- **In the default-prior fits, almost nothing moves.** This community's true mean occupancy intercept is about -0.1 on the logit scale, where the fixed prior was already centred, so the community mean changes estimates by at most about half a percentage point. The pull towards the middle in Lesson 2 is unchanged (for example, PCR fit below 20%: 11.8 points too high before, 11.7 now).
- **mu_b0 converges in every default-prior fit at 3,000 + 6,000** (R-hat at most 1.0003, ESS at least 8,394), so the lesson burn-in was not lengthened.
- **The alternative-prior fit (Beta(1, 4) contamination priors) now has a stuck chain.** In both its shorter and its longer run, chain 4 settles in a mirror explanation for OTU_10: high occupancy (intercept 1.58, truth -0.68) with most detections blamed on field contamination (`theta0` 0.69, truth 0.08). The other three chains agree with each other (intercept -1.8, `theta0` 0.11). The longer fit's largest R-hat is 9.24 and 27 parameters are flagged (2 before); mu_b0's R-hat is 1.28. A burn-in of 15,000 does not free the chain. With the same seeds and `b0_prior = "fixed"`, no chain enters the mode (4 flags, largest R-hat 1.013). That is one seed, so it shows the mode is reachable under the new default; it does not show the community prior causes it. Lesson 2's prior-sensitivity section and Lesson 3's diagnostics section are built on this fit, and Lesson 3 currently says such a mirror explanation was found only in a separate study.
- **The Quickstart's featured flag is now a covariate effect, not `theta0`.** In the refitted `quickstartresults` (same seed and settings) no `theta0` is flagged (largest R-hat 1.004, 1.039 before), which stopped the render at chunk `featured-flag`. Doug chose option (c): ship the refit and feature the flagged parameter with the highest R-hat. That is `beta_psi` for `protectedprotected` in Mammalia_Carnivora_3 (R-hat 1.015). Only code was changed; the sentences around it that no longer fit are listed below with suggested replacements, and the Quickstart renders against the installed package.
- **`computeSiteWAIC()` fails with its default quadrature orders on the new one-factor fit.** One posterior draw has a residual SD of 4.07, and the default ladder (15, 25, 41, 65) does not converge there. Lesson 4's displayed call (`computeSiteWAIC(fit_one_factor, threshold = 1)`, all draws) would stop with "Joint quadrature did not converge". The lesson's WAIC bundle was computed with orders 101 and 161 added (and 181 and 241 on the stricter check), which converge and agree with each other.

## What was refitted

- Build archive: `/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-community-prior-20261104`, with its own `library` (HEAD 1b4a2cd source, R 4.6.1, Accelerate BLAS), `run-rebuild.py` (the exact steps), `logs/` and `build-steps.jsonl`. Recorded source commit 739a1f2 (same `R/` and `src/` as 1b4a2cd). AC power, Low Power Mode off, one thread per fit.
- Same seeds, data, priors, MCMC schedules and fit arguments as the 7 October build; `b0_prior` not passed. Fit times: perfect 35 s, default 119 s, alternative 125 s, longer alternative 216 s, one-factor prediction 79 s, unbalanced 84 s.
- Shipped data: `data/sampleresults.rda` (seed 20261007, 2 chains of 5,000 + 5,000, 1.1 minutes) and `data/quickstartresults.rda` (seed 20261101, 2,000 + 2,000, `nthin = 2`, 0.5 minutes), by `dev/release/2026-11-04-community-prior-refit/refit-shipped.R`, which also saved their diagnostics and provenance there.
- Teaching bundles regenerated in `vignettes/teaching-data/`: `nonspatial-lesson`, `output-lesson`, `diagnostics-lesson`, `latent-presence-lesson`, `native-plots`, `ordination-examples`, `remaining-plots-data`, `native-traits-data`, `covariate-effect-data`, `prediction-lesson`, `unbalanced-lesson`, `site-waic-lesson`, and the PNGs their exporters write. WAIC tables are in `community-prior-validation/`.
- Not refitted, because they come from separate studies: `jsdm-comparison.rds`, `lesson-4-calibration.rds`, `lesson-4-extension.rds` and `lesson-4-species-errors.csv` (Lessons 5 and 6), `spatial-lesson.rds` (Lesson 7) and `mirror-labelling-chains.csv` (interval-calibration study). All were fitted with the fixed prior; Lessons 5 to 7 still describe fixed-prior fits.

## Convergence of mu_b0 and of each fit

From `community-prior-validation/mu-b0-convergence.csv` (`check-mu-b0.R`) and the shipped-fit diagnostics. "Flagged" uses R-hat above 1.01 or ESS below 400 over `returnConvergenceDiagnostics()`.

- Perfect observation (binary, 3,000 + 6,000): mu_b0 R-hat 1.000, ESS 17,091; 0 of 31 flagged.
- Default PCR fit: mu_b0 R-hat 1.000, ESS 9,260; 0 of 101 flagged.
- Alternative priors, shorter (3,000 + 6,000): mu_b0 R-hat 1.241, ESS 3,522; 30 of 101 flagged, largest R-hat 9.30.
- Alternative priors, longer (6,000 + 12,000), the one the lessons use: mu_b0 R-hat 1.278, ESS 6,258; 27 of 101 flagged, largest R-hat 9.24.
- Diagnostic only, alternative priors at 15,000 + 6,000: mu_b0 R-hat 1.24; 28 flagged, largest R-hat 9.21 (chain 4 still stuck).
- Diagnostic only, alternative priors with `b0_prior = "fixed"` at 3,000 + 6,000: 4 of 100 flagged, largest R-hat 1.013.
- One-factor prediction fit: mu_b0 R-hat 1.000, ESS 8,775; 2 of 101 flagged (`beta_psi` for X_psi.EnvCov.1 of OTU_2 and OTU_3, R-hat 1.01).
- Unbalanced fit: mu_b0 R-hat 1.000, ESS 8,395; 0 of 101 flagged.
- `sampleresults`: mu_b0 R-hat 1.001, ESS 614; 34 of 131 flagged (31 of 130 before), largest R-hat 1.07 (1.10 before); the spatial blocks (`Bs`, `Gs`, `As`, `Cs`) still warn.
- `quickstartresults`: mu_b0 R-hat 1.000, ESS 1,832; 10 of 221 flagged (4 of 220 before): 6 for R-hat (largest 1.02, all `beta_psi` or `beta0_psi`) and 4 for ESS (smallest 240).
- Quickstart data at 10,000 + 10,000, `nthin = 1` (claim check, not shipped): 3 of 221 flagged (R-hat 1.01 to 1.02); mu_b0 R-hat 1.00, ESS 3,738; 1.4 minutes.

The collapsed shift step accepted 44% of proposals in every two-stage fit.

## Decisions for Doug

- **Quickstart traceplot section (decided 4 November: option (c)).** The featured parameter is now chosen by code; the prose around it still describes `theta0` and needs Doug's edits, listed under Quickstart below.
- **The alternative-prior example.** Keep the stuck chain and teach it (it is a real mirror explanation, in the lessons' own data, and Lesson 3 already has a section for it), or refit that arm with another seed, or with `b0_prior = "fixed"` to keep the old story. If the last, Lessons 2 and 3 should say why that arm differs from the default.
- **WAIC quadrature.** Raise `computeSiteWAIC()`'s default orders, or make it extend the ladder before it stops, so the call Lesson 4 displays works on these fits. Until then, Lesson 4's displayed code fails for a reader with the full fits.

## Changes by file

Source lines are in the `.Rmd`; "md" lines are in the regenerated `.md` render. Items marked **prose** are hard-coded text that is now wrong and needs an edit. The others are inline R or code output and update by themselves, but read differently.

### README.md

- Line 83, **prose**: "`sampleresults` was refitted under per-chain seeding and retains some diagnostic flags ... See the refit report (2026-10-13)". It is now the 4 November refit at the community prior (`dev/release/2026-11-04-community-prior-refit/`), with 34 flags.
- Line 74, **prose**: the numbers for pull towards the middle come from the fixed-prior `current-main-recheck` study, not from the new default. Simulation studies show that the community mean removes most of the overestimation of rare species' occupancy, so this bullet describes the old default.

### R/data.R (not edited)

- `quickstartresults` docs (lines 161 to 166): "The script that made it, with its record, is in `dev/release/2026-11-01-quickstart-data/`". It is now made by `dev/release/2026-11-04-community-prior-refit/refit-shipped.R`, and the documented call and "default priors" still hold. "A few parameters miss the quickstart's convergence rule" still holds (10 do). Not edited, because any change in `R/` changes the source hashes that every lesson fit records, and the build scripts would then refuse the archive.

### Quickstart (`occJSDM-00-quickstart.Rmd`)

**Code changed** (option (c), Doug's choice of 4 November). Chunk `featured-flag` (lines 159 to 192) now features the flagged parameter with the highest `rhat`, using the reader's rule (`rhat > 1.01`). It compares it with the same parameter (same covariate or primer) in the species where that parameter's `rhat` is lowest. The guard is kept, as `stopifnot(featured_rhat > 1.01)` on the picked parameter. The code works for any parameter type, so a later refit cannot silently feature an unflagged one. The visible chunk `traceplot` (lines 196 to 200) now plots `featured_draws` with `param_name = featured_param`, because it used to plot `theta0` directly. In the shipped fit the featured parameter is the effect of protection (`beta_psi`, `label1` `protectedprotected`) on Mammalia_Carnivora_3, R-hat 1.015 (shown as 1.02 in the table), ESS 604. The comparison species is Mammalia_Carnivora_2, R-hat 1.000. The chain means are 1.79 and 1.67 on the logit scale, and the chains overlap broadly.

**Sentences that no longer fit** (line numbers are in the edited `.Rmd`; quoted text shows inline R as rendered now):

- Line 143, inline: "The table has 11.05 rows for each of the 20 species". The diagnostics table gains one `mu_b0` row (221 rows, not 220). Suggested: "The table has `r (nrow(diagnostics) - 1) / n_species` rows for each of the `r n_species` species, and one for the whole community:". Then add a bullet after line 149: "- `mu_b0`: the community's average baseline occupancy, on the logit scale, which each species' baseline is drawn towards (one for the whole community)."
- Line 194. Old: "10 of the 221 parameters are flagged here: 6 for a high `rhat` and 4 for a low `ess`. Take the field false-positive rate (`theta0`) of Mammalia_Carnivora_3, whose `rhat` is 1.015. `plotTraceplot()` shows its two chains beside those of Mammalia_Carnivora_2, whose `theta0` has the lowest `rhat` (1.000):". (Before the refit: "4 of the 220 ...: 4 for a high `rhat` and 0 for a low `ess`. Take the field false-positive rate (`theta0`) of Mammalia_Artiodactyla_2, whose `rhat` is 1.039".) Suggested: "`r nrow(flagged)` of the `r nrow(diagnostics)` parameters are flagged here: `r n_rhat_flags` for a high `rhat` and `r n_ess_flags` for a low `ess`. Take the one with the highest `rhat`: the effect of protection (`beta_psi` for `protectedprotected`) on `r featured_species`, whose `rhat` is `r sprintf("%.3f", featured_rhat)`. `plotTraceplot()` shows its two chains beside those of `r healthy_species`, whose protection effect has the lowest `rhat` (`r sprintf("%.3f", healthy_rhat)`):". If the wording should survive a future refit, "the effect of protection (`beta_psi` for `protectedprotected`)" can be made inline from `featured$param` and `featured$label1`.
- Line 202. Old: "For Mammalia_Carnivora_2 the two chains overlap throughout. For Mammalia_Carnivora_3, chain 1 makes more frequent and higher excursions than the other, so its mean is 1.79 against 1.673. In short, the chains disagree about how high the rate can go, and that disagreement is what `rhat` measures." The figure does not show larger excursions in one chain, and "the rate" no longer applies. Suggested: "For `r healthy_species` the two chains overlap throughout. For `r featured_species` they overlap too, but chain `r high_chain` sits slightly higher than the other, with a mean of `r round(max(chain_means), 2)` against `r round(min(chain_means), 2)` on the logit scale. That small disagreement between chains is what `rhat` measures, and here it is just over the rule's limit." The rounding changes from 3 to 2 digits, because the current output mixes 1.79 with 1.673.
- Line 206. Old: "For instance, a flagged field false-positive rate (`theta0`) matters less if you are interested in covariate effects (`beta_psi` and `beta_theta`), but it matters more if you want to know where a species really occurs." The flagged parameter is now itself a covariate effect. Suggested: "For instance, this flag matters if your question is how Red Listed species respond to protection, because `r featured_species` is Red Listed, but much less if you are interested in the field or PCR error rates."
- Line 207. Old: "For this dataset and fit, we have found that a run of 10,000 burn-in and 10,000 kept iterations, without thinning, cleared every flag. More iterations do not always help, though." At the new default, 10,000 + 10,000 (seed 20261101, 1.4 minutes) leaves 3 of 221 flagged: `beta0_psi` of Mammalia_Rodentia_3, `beta_psi` for `protectedprotected` of Mammalia_Artiodactyla_3 and `theta0` of Amphibia_Anura_1, all with R-hat 1.01 to 1.02 and ESS above 1,000. The featured parameter is cleared. Suggested: "For this dataset and fit, we have found that a run of 10,000 burn-in and 10,000 kept iterations, without thinning, cleared the flag above and left three parameters just over the rule (`rhat` 1.01 to 1.02). More iterations help, but they do not always clear every flag."
- Line 208. Old, as rendered: "Mammalia_Carnivora_3 is estimated to occupy about 52 of the 100 sites. Its field false-positive rate can only be learned from the sites it does not occupy, and here there are few of them." The second sentence is about `theta0`, and its inline condition (occupied sites over half) now picks "few of them" for 48 unoccupied sites. Before, the sentence was about Mammalia_Artiodactyla_2 at about 80 of 100. Suggested: "`r featured_species` is estimated to occupy about `r featured_occupied` of the `r n_sites` sites. Its response to protection can only be learned from the contrast between its protected and unprotected sites, and the survey has 46 protected and 54 unprotected sites, so each side holds only a couple of dozen occupied sites." Alternatively, drop the species-specific sentence and keep only the general advice. The flagged-species list in the same bullet updates by itself (now six species).
- The four ESS flags (`beta_psi` for `protectedprotected` of Mammalia_Rodentia_1, Mammalia_Eulipotyphla_2 and Amphibia_Anura_3, ESS 240 to 399, and the collection intercept of Mammalia_Rodentia_2, ESS 388) are new; before, no parameter failed on ESS. No sentence is wrong because of them, but the section never shows an ESS flag. An optional sentence after line 202: "The other kind of flag, a low `ess`, means the chains agree but move slowly, so fewer independent draws stand behind the summary; here the slowest is the protection effect of Mammalia_Eulipotyphla_2, with an `ess` of 240."
- Line 87, "Fitting live takes just over a minute": the refit took 0.5 minutes on 2 cores (AC power, Low Power Mode off), so "under a minute" may now fit. Not changed by the prior as such; check on Doug's machine.

**Inline values that changed** (no edit needed): baseline range (line 235) 0.22 to 0.78, now 0.21 to 0.78; trait coefficients (line 343) `log10_body_mass` on `hmi` -0.55 (-0.91 to -0.22), now -0.54 (-0.87 to -0.24), and `red_listed` on protection 0.81 (0.25 to 1.38), now 0.80 (0.22 to 1.35); positively associated pairs (line 467) 38, now 37; typical interval width (line 493) 3.1, now 3; protection effect at new sites (line 593) 0.14 to 0.43, now 0.14 to 0.44.

**Records and tests:** `tests/testthat/test-quickstart-data.R` expects 2 to 10 flags in `quickstartresults`; the refit has exactly 10, so it passes at the upper edge. The test was left unchanged. `dev/release/2026-11-01-quickstart-data/RECORD.md` describes the fixed-prior fit (its flags, the 10,000 + 10,000 run that cleared every flag, the trait coefficients, 76 of 190 pairs); it is a dated record and was not edited.

### Lesson 0

- md "The simulation was made at occJSDM revision": old a56f554, R 4.5.0; new 739a1f2, R 4.6.1.
- md 1004: `identical(field_off$data_list, survey_data)` old FALSE, new TRUE. The bundle's simulation is now made by the same package and R version as the reader's, so the regenerated survey matches exactly.

### Lesson 2

- Line 330, inline plus **prose** (md 408): old "the alternative-prior fit has 2 parameters outside one screen or the other, with a largest Rhat of 1.013 and a smallest ESS of 397"; new 27, 9.239 and 395. The table that follows (md 417 to 445) now lists 27 rows, including mu_b0 (R-hat 1.278) and every OTU_10 parameter.
- Line 339, **prose**: "Treat small differences under those alternative priors cautiously; we do not claim that every parameter has fully converged." This understates it now: one chain has a different explanation for OTU_10.
- Diagnostics summary table (md 402 to 404): alternative fit largest R-hat old 1.013, new 9.239; largest R-hat of occupancy probabilities old 1.007, new 1.531; flagged old 1, new 27. Default fit largest R-hat old 1.009, new 1.006; smallest ESS old 835, new 907.
- Error table (md 588 to 597) and lines 473 to 481 inline: perfect fit 20% to 80% signed 0.5, now 0.4; PCR fit all-pairs MAE 15.4, now 15.5; PCR signed below 20% 11.8, now 11.7; PCR signed above 80% -14.7, now -15.0; PCR signed 20% to 80% 3.1, now 3.0. Perfect-fit bias numbers (6.2 and 8.8) unchanged.
- Coverage (md 632 to 634, inline): perfect 89.2%, now 89.5%; PCR 76.8%, now 77.0%; misses 23.2%, now 23.0%.
- Case table and case paragraphs (md 862 to 903, inline): weak true detection site 54.9%, now 54.4%, sample 17.7%, now 17.4%; laboratory false positive site 98.7%, now 98.3%, sample 1.7%, now 1.8%; strong true detection site 98.4%, now 98.3%; field-stage false positive site 0.7%, now 0.6%. Fitted underlying probability at the weak case's site 86.1%, now 86.3% (md 977).
- Field-stage summaries (md 1028 to 1121, inline): mean site probability of field-stage false positives 47.4%, now 47.3%; lone contaminated sample average 44.4%, now 44.3%, maximum 89.7%, now 89.4%; highest fitted occupancy 82.1%, now 81.7%; and other changes of 0.1 to 0.4 points in the same tables.
- Line 998, inline plus **prose** (md 1193): old "the field-stage false-positive case receives 0.6% probability of site presence, close to its 0.7% under the default priors. Although this case is rejected under both priors"; new 24.9% against 0.6%. The case is OTU_10 at site 1, the species whose chain 4 is stuck, so "close to" and "rejected under both priors" are now false. Mean site probability of field-stage false positives under the alternative priors: 41.4%, now 43.8%. Mean absolute error, default then alternative: 15.4 to 16.2 before, 15.5 to 16.9 now. "The alternative priors therefore worsen overall recovery" still holds, but part of that is now the stuck chain.
- Figure `prior-sensitivity-cases` (line 960): the grey triangle for the field-stage case moves from near 0 to about 25%.
- Lines 509 to 515, **prose**, section "Why low and high probabilities are pulled towards the middle": the second paragraph says "The prior for baseline occupancy, `B0` ... the default puts about 95% of its weight on baselines between 12% and 88%, so a species whose baseline lies outside that range is pulled towards it." This is no longer the default behaviour. By default, each species' `B0` now has a Normal(`mu_b0`, 1) prior around an estimated community mean, so species are pulled towards their own community's average, not towards 50%. The 12% to 88% range now describes how far a species can sit from the community mean, and the fixed centre applies only with `b0_prior = "fixed"`. The experimental `sigma_b0` blockquote (line 515) describes studies run with the fixed prior. Lesson 6's line 614 points to this paragraph for "how the default prior pulls rare and common species towards the middle", so it needs the same change.
- Line 1046, **prose** summary bullet: "Occupancy estimates are pulled towards the middle: rare species' probabilities come out too high". Still true of this lesson's fits, whose species are not rare as a community. With the community prior that is mainly the within-species pull described in line 511.
- Line 1042 (md 1244), inline: complete-survey MAE 15.4, now 15.5 (the reduced fit stays 15.4). The phrase "against 15.4 for the complete survey" now reads "15.4 ... against 15.5". Replication table (md 1356 to 1361) changes by 0.1 to 0.2 points; signed error for the three affected sites -0.6, now -0.8.
- Unbalanced diagnostics (md 1363 to 1390, inline): largest parameter Rhat 1.004, now 1.009; occupancy-probability Rhat 1.008, now 1.007; minimum ESS 1,221, now 1,279.
- Line 1201 (md 1423), inline: revision a56f554 and R 4.5.0, now 739a1f2 and R 4.6.1.
- Line 294 says all fits use four chains of 3,000 + 6,000; still true.

### Lesson 3

- Line 142 (md 171), inline: perfect fit's table 30 parameters, now 31; smallest ESS 2,097, now 2,092. PCR fit 100 parameters, now 101; largest Rhat 1.009, now 1.006; smallest ESS 835, now 907. Both gain the `mu_b0` row.
- Diagnostic summary table (md 168 to 169) and the slope and trait diagnostics (md 450 to 453): small changes, for example PCR environment ESS 736, now 817.
- Line 200 onwards, section on the longer alternative-prior fit: the flagged table (md 244 to 270) now has 27 rows, not 2, and includes OTU_10's `theta0` (R-hat 9.24, ESS 18,764).
- Line 213, **prose**: "`theta0` for OTU_6 has a low effective sample size even though its Rhat is below 1.01 ... The first environmental slope for OTU_6 is flagged for Rhat instead." Now OTU_6's `theta0` has R-hat 1.015 and ESS 395, so it is flagged on both, and OTU_6's first slope is no longer flagged.
- Line 292, **prose**: the selection rule (lowest `theta0` ESS) still picks OTU_6, so the trace figure `flagged-field-trace` still shows OTU_6. The parameter the table now points to, though, is OTU_10's `theta0`.
- Line 316 (md 383), inline: "The package's ESS for this parameter is about 397" now renders as "about 1.8764^{4}, 395", because two `theta0` rows are flagged and the inline code returns both. This is a rendering defect to fix with the prose.
- Line 322, **prose**: "One of our simulation studies found such a mirror explanation when laboratory contamination was set far above what the default priors assume." The lessons' own longer alternative fit now shows one (OTU_10, chain 4). Lines 1924 to 1939, the appendix, say "This is the only such case found in these studies".
- Line 748 (md 794), inline: the generating curve lies inside the 95% band for six of ten species, now five; largest gap 0.32, now 0.33 (OTU_5).
- Baseline table (md 844 to 853, inline): every species' baseline estimate moves by at most 2.4 points (OTU_6 lower bound 57.4%, now 55.0%).
- Latent-presence tables (md 1539 to 1610 and 1680 to 1691): site-presence probability for OTU_1 at site 1, 0.987, now 0.983, and similar third-decimal changes.
- Ordination scores (md 1881 to 1889): small changes. The TODO note on the ordination contribution (correlations 0.67 and 0.49, slopes 0.31 and 0.04, from `output-lesson.rds`) should be rechecked against the new bundle.
- Figures regenerated: `collection-trace`, `environmental-coefficients`, `flagged-field-trace`, `ordination-contribution`, `primer-trace`, `trait-effects`, and the native-plot, native-trait and covariate-effect PNGs in `teaching-data/`.

### Lesson 4

- Line 230 (md 232), inline: average absolute error within the training ranges 8.8, now 8.9 (10.9 outside, unchanged).
- Line 246 (md 253), inline: mean overestimation 0.52 and 0.36 points, now 0.33 and 0.20; average absolute error "about 8.9" unchanged.
- Scores (md 296 to 299, inline): negative log scores 0.5546 and 0.5538, now 0.5547 and 0.5541; Brier with two hidden covariates 0.1868, now 0.1870.
- Line 304 (md 322), inline: paired Brier difference 0.00033, now 0.00021, still "slightly favouring the fit with two". SE 6.6e-05, now 5.3e-05 (also md 479).
- Line 358, **prose** (md 381): "There are no flagged rows under these checks." The one-factor fit's parameter table now has 2 flags (R-hat 1.01, `beta_psi` of OTU_2 and OTU_3), though the probability checks the sentence refers to still pass (largest Rhat 1.004, smallest ESS 1,089, both inline). The diagnostics table (md 378 to 379) shows 101 parameters and 2 flags.
- Line 387, displayed code `computeSiteWAIC(fit_one_factor, threshold = 1)`: fails on the new full fit with the default quadrature (see Summary).
- Line 395, **prose**: "PR #14 merged on 7 October 2026. We rebuilt this three-sample survey with the same seeds, priors and MCMC schedules. All six fits' posterior results were identical to those in the previous archive." No longer true: the fits are now the 4 November community-prior refits, whose posteriors differ. "checked on 4,000 of each fit's 24,000 retained draws" still holds, but with the wider quadrature orders.
- Out-of-sample WAIC (md 402 to 403): 24,495.31 and 24,501.00, now 24,498.09 and 24,504.65. Site WAIC (md 439 to 440): 23,121.60 and 23,120.29, now 23,121.33 and 23,119.88; sites flagged for large penalties 99 and 99, now 100 and 99. Comparisons (md 455 to 456): difference at 1,000 draws 0.47, now 1.68; at 4,000 draws 1.31, now 1.45; paired SE 0.76, now 0.67.

## Verifier and script changes made for this rebuild

- `verify_diagnostics.R`: expected flag count in the longer alternative fit, 2 before, 27 now (dated comment).
- `prediction-verify.R`: tolerance of the single-draw quadrature check, 1e-8 before, 1e-7 now. At the draw with residual SD 4.07, the verifier's own 121-node rule differs from adaptive integration by 3.6e-8; exported means agree to 5e-10.
- `validate-waic.R`: optional fourth argument `wide`, which appends orders 101 and 161 (181 and 241 on the stricter check). Both WAIC runs used it.
- `export-waic.R`: records the export date and the archive's source commit instead of the hard-coded 7 October values, and takes the record folder as an optional argument.
- `check-mu-b0.R`: new; tabulates mu_b0 and overall convergence for every lesson fit.
- Pre-existing drift, not fixed: `verify_native_plots.R`, `ordination-verify.R`, `remaining-plots-verify.R`, `native-traits-verify.R` and `covariate-effect-verify.R` compare their example Rmd chunks with Lesson 3's, and `unbalanced-verify.R` compares with Lesson 2's. Doug's prose edits (13d9eb4, fdcfec3) changed a plot title and a theme line in Lesson 3, and Lesson 2's displayed fit now adds `ncores = 4L`. These comparisons fail regardless of the refit. In this rebuild, each verifier was run as committed. Where only that comparison failed, its data checks were rerun from a copy without it (`verify-relaxed.sh` and the `*-data-only` or `*-snippet-chunk` logs in the archive), and every data check passed. Updating the example Rmd files to match the lessons is left for after Doug's review.

## Reproduce

From the repository root, with the archive created as above:

```sh
A=/Users/douglasyu/src/occJSDM-worktrees/lesson-archive-community-prior-20261104
R CMD INSTALL --preclean --library=$A/library .
python3 $A/run-rebuild.py          # the lesson chain; waits for $A/logs/sampleresults.done
$A/run-shipped.sh                  # sampleresults, then quickstartresults
R_LIBS=$A/library Rscript dev/simstudy/vignette-lesson/check-mu-b0.R $A dev/simstudy/vignette-lesson/community-prior-validation/mu-b0-convergence.csv
```

The diagnostic alternative-arm fits are `$A/diagnose-alternative/fit.R` (logs and summaries beside it), and the Quickstart claim check is `$A/quickstart-long-check.R`.

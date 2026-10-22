# Task 4: the other flagged fits and the impact on earlier conclusions

`impact.R` (tests in `test-impact.R`) reads saved fits and archived tables and fits nothing. Numbers are in percentage points of occupancy probability unless stated, and come from the CSVs in `impact/`. The last section is interpretation and is marked as such.

## Answers

- **The other four flagged fits need no extended run.** None has a separated species. Their flags come from scored Rhats between 1.05 and 1.075, nearly all on original-site cell probabilities and all of a single species per fit, and from a package warning on `B_output` in two fits; the four-sample fit has only the warning (`impact/other-flagged-fits.csv`).
- **The extended run of community 5 is unflagged under the full pr11 flag rule.** Applied to its 16 chains bound into one fit, the rule finds a largest group Rhat of 1.0040, a largest element Rhat of 1.0059, no non-finite Rhat and no package warning; the same code reproduces the pr11 fit's recorded flag to within 5e-15 (`impact/flag-rule.csv`).
- **No practical conclusion of the current-code comparison reverses.** With community 5 of the high-contamination 300-site design scored from the near-truth chains, the mirror chains or the extended run, every affected result with a published observed-chain sensitivity range stays inside it, and every statement holds except one that fails only if community 5 is scored from the mirror chains alone; it holds with the published pooled estimate, which weights the two labellings equally, with the near-truth chains and with the extended run (`impact/affected-comparisons.csv`).
- **The intercept-prior phase B verdict stays FAIL for SD 2, 3 and 5.** With the extended run as the control for community 5, criterion 4 still fails at both contamination levels for every SD, even if each SD's own community 5 fit is also counted as unflagged; criteria 1 to 3 keep their pass or fail result for every SD, although one component (SD 3, overall MAE at high contamination) moves from pass to fail (`impact/intercept-prior-phase-b.csv`).

## The other four flagged fits (Step 1)

- `design-qnear_K6-sites300-02`: species 8 drifting (`theta0`, collection intercept and slope, both occupancy slopes; separation at most 0.66, all-chain Rhat at most 1.074). Flag: 35 scored diagnostics above 1.05, all of species 8 (31 original-site cell probabilities, 2 occupancy slopes, the collection intercept and `theta0`; at most 1.074), and one package warning on `B_output`.
- `design-qfar_K6-sites300-07`: species 10 drifting (second occupancy slope; separation 0.77, Rhat 1.053). Flag: 18 scored diagnostics of species 10 (17 cell probabilities, 1 slope; at most 1.055); no package warning.
- `design-qfar_K6-sites300-09`: species 4, 6, 7 and 9 drifting (separation at most 0.505, Rhat at most 1.050). Flag: 12 original-site cell probabilities of species 9 (at most 1.054); no package warning.
- `design-qfar_K6-field4-09`: species 3 and 6 drifting (first occupancy slope; separation at most 0.71, Rhat at most 1.046). Flag: one package warning on `B_output` only; no scored diagnostic exceeds 1.05.
- The largest separation in these four fits is 0.505 to 0.77. The unflagged reference fits `-02` and `-03` show comparable drift (largest separations 0.47 and 0.31; Task 1), all far below the separation threshold of 3.

## Is the extended run an unflagged control?

The pr11 rule (`current-main-recheck/select.R` on the diagnostics of its `run.R`) flags a fit with any fitting warning, any group or element Rhat above 1.05, or any non-finite or nonpositive Rhat. A single-chain fit cannot be judged by it, since the package computes no Rhat for one chain, so the 16 extended chains were bound along the chain dimension into one fit and judged as a fit that had run them together (`bind_chains`). The group and element Rhats are those of the intercept-prior study's `fit_flags` for phase B, and the fitting warnings those of the package's own `computeDiagnostics` on the bound output.

- Validation on the saved pr11 fit: the rule gives 7 fitting warnings, a largest group Rhat of 1.7091 and a largest element Rhat of 1.7375, flagged; these match the intercept-prior control record and the current-code diagnostic summary within 4.7e-15, and the recomputed warnings are the fit's 7 saved warnings.
- Extended run, 16 chains of 10,000 draws: 0 fitting warnings, largest group Rhat 1.0040, largest element Rhat 1.0059, no non-finite Rhat: unflagged.

## Community 5 scored four ways (Step 2)

Original 100 sites, 1,000 cells; bands by the generating probability (`impact/community5-scores.csv`). Coverage is the phase B 95% interval coverage, the only scorer here that defines one.

- pr11 fit, all four chains pooled (as published): MAE 21.17; signed error below 0.2 +17.00, above 0.8 -28.41; coverage 0.808; largest cell Rhat 1.74.
- pr11 chains 1 and 3 (near-truth): MAE 19.16; +12.46 and -25.52; coverage 0.803; largest cell Rhat over these two chains 1.04.
- pr11 chains 2 and 4 (mirror): MAE 23.71; +21.54 and -31.30; coverage 0.751; 1.01.
- Extended run, 16 chains pooled: MAE 19.17; +12.44 and -25.75; coverage 0.816; largest cell Rhat 1.006.

## Affected comparisons of the current-code report

Each item quotes the report, then gives the published value, the value with community 5 from the near-truth chains, the mirror chains and the extended run, and the published observed-chain sensitivity range. Only the 300-site fit of community 5 is replaced. The CSV quotes each sentence in full.

- "Using 300 sites gives 15.20 and 16.97 points." At high contamination: 16.97; 16.77, 17.22, 16.77; range 16.68 to 17.31. Inside the range.
- "Every community improves relative to the baseline in each of these comparisons." Baseline to 300 sites, high contamination: reduction 3.55 (95% interval 2.31 to 4.79, 10 of 10 improve); 3.75 (10 of 10), 3.30 (9 of 10), 3.75 (10 of 10); range 3.21 to 3.84. Inside; the reduction stays positive throughout. With the mirror chains alone community 5 no longer improves on its baseline, so the statement would fail in that case only. It holds with the published pooled estimate, which weights the two labellings equally, with the near-truth chains and with the extended run, whose 16 chains all stayed in the near-truth region (Task 3).
- "In the two-stage model, four field samples reduce error from 17.97 to 15.26 points under low contamination and from 20.52 to 17.23 under high contamination." 3.29 in every case, since it does not use this fit.
- "The additional advantage of 300 sites over four field samples is small: 0.06 points under low contamination and 0.26 under high contamination. Their paired intervals span zero (-0.98 to 1.11 and -0.90 to 1.42 points). These data do not identify a clear winner between those designs." At high contamination: 0.26 (-0.90 to 1.42); 0.46 (-0.49 to 1.42), 0.01 (-1.57 to 1.58), 0.46 (-0.49 to 1.41); range -0.12 to 0.61. Inside; the interval spans zero in every case, and 6 of 10 communities favour 300 sites in every case.
- "Low probabilities remain too high and high probabilities remain too low." At 300 sites and high contamination, below 0.2: +15.34; +14.89, +15.80, +14.89; range +14.72 to +15.89. Above 0.8: -18.82; -18.53, -19.11, -18.55; range -19.37 to -18.29. Inside; the signs hold.
- "Every paired 95% interval for these eleven changes includes zero." Current minus archived code, 300 sites, high contamination: 0.249 (-0.197 to 0.695); 0.047 (-0.140 to 0.235), 0.503 (-0.495 to 1.501), 0.048 (-0.138 to 0.235). The interval still includes zero; no observed-chain range was published for this change.
- "Across the ten-community summaries, overall mean absolute error changed by at most 0.025 percentage points in binary JSDM and 0.249 points in the two-stage comparisons." With the near-truth chains or the extended run the two-stage figure becomes 0.226, attained by the known-site-factor arm at high contamination; with the mirror chains, 0.503. A reported value, not a reversal.

## Intercept-prior phase B (Step 3)

The extended run's phase B scores replace the control's community 5 at high contamination, and that control is counted as unflagged, as the pr11 rule finds above. Descriptive only; the study's files are unchanged.

- Criterion 4, high contamination: control flagged fits 3 become 2, against 4, 9 and 10 for SD 2, 3 and 5, so it fails for every SD, by a wider margin. Low contamination is unchanged (control 1; SDs 3, 7 and 9) and fails for every SD on its own. The FAIL verdict therefore cannot change.
- Criterion 4 if each SD's own flagged community 5 fit at high contamination were also counted as unflagged: 3, 8 and 9 against the control's 2, which still fails for every SD, and low contamination fails independently. So FAIL holds even allowing for the caveat below.
- Criterion 1: the control's pooled low band goes from 14.88 to 14.65; SD 2, 3 and 5 stay better (13.21, 12.41, 12.11). SD 2 improves in 9 of 10 communities instead of 10 (7 needed). Passes for every SD, as before.
- Criterion 2, high contamination: the control's bands become 14.89 (low), 1.75 (middle) and 18.55 (high), and its MAE 16.77 instead of 16.97. SD 3's MAE is then 0.68 worse, over the 0.5 limit (was 0.48), so that component now fails; SD 3 already failed criterion 2 on the middle band. SD 2's MAE difference becomes +0.19 (was -0.01) and SD 5's +0.93 (was +0.73). The criterion-level results are unchanged: pass for SD 2, fail for SD 3 and 5.
- Criterion 3: the control's coverage changes by at most 0.0031 (high band 0.716 to 0.719); passes for every SD, as before.
- Caveat: the new arms' own community 5 fits at high contamination are flagged too (largest element Rhat 1.54, 1.28 and 1.26 for SD 2, 3 and 5) and were not replaced, so criteria 1 to 3 are not a like-for-like comparison of unflagged fits.

## Checks

33 checks pass (`impact/reproduction.csv`). The published community 5 rows of `community-scores.csv` are reproduced from the saved pr11 fit within 4.5e-12 (every metric and numeric column), its band errors within 4.4e-16, its per-chain scores within 9.4e-16, and its phase B control scores within 4.3e-14. The pr11 flag rule reproduces the fit's recorded flag within 4.7e-15, with its 7 saved warnings. Every affected published summary is recomputed from the published scores and from the reproduced scores within 5e-14 points, and the phase B gate detail within 6.0e-14 with identical pass values. Each substitution leaves every other row unchanged. The extended chains reproduce their stored posterior means within 1.2e-14.

## Interpretation (not a result)

- The published pooled estimate of community 5 averaged two chains in each mode. The unflagged extended run agrees with the near-truth chains to within 0.02 points of MAE, so the published value is about 2.0 points higher for community 5's error, and for the ten-community mean of the high-contamination 300-site design about 0.20 higher. Using the extended run slightly favours 300 sites and removes most of the apparent code-version change in that arm, but changes no conclusion.
- In phase B an unflagged control is slightly more accurate, which makes the wider priors look slightly worse; it strengthens rather than weakens the FAIL verdict.

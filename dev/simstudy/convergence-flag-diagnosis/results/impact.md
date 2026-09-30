# Task 4: the other flagged fits and the impact on earlier conclusions

`impact.R` (tests in `test-impact.R`) reads saved fits and archived tables and fits nothing. Numbers are in percentage points of occupancy probability unless stated, and come from the CSVs in `impact/`. The last section is interpretation and is marked as such.

## Answers

- **The other four flagged fits need no extended run.** None has a separated species. Their flags come from scored Rhats between 1.05 and 1.075, nearly all on original-site cell probabilities and all of a single species per fit, and from a package warning on `B_output` in two fits; the four-sample fit has only the warning (`impact/other-flagged-fits.csv`).
- **No practical conclusion of the current-code comparison reverses.** With community 5 of the high-contamination 300-site design scored from the near-truth chains, the mirror chains or the converged extended run, every affected result with a published observed-chain sensitivity range stays inside it, and every statement holds except one that fails only with the mirror chains alone (`impact/affected-comparisons.csv`).
- **The intercept-prior phase B verdict stays FAIL for SD 2, 3 and 5.** With a converged control for community 5, criterion 4 still fails at both contamination levels for every SD; criteria 1 to 3 keep their pass or fail result for every SD, although one component (SD 3, overall MAE at high contamination) moves from pass to fail (`impact/intercept-prior-phase-b.csv`).

## The other four flagged fits (Step 1)

- `design-qnear_K6-sites300-02`: species 8 drifting (`theta0`, collection intercept and slope, both occupancy slopes; separation at most 0.66, all-chain Rhat at most 1.074). Flag: 35 scored diagnostics above 1.05, all of species 8 (31 original-site cell probabilities, 2 occupancy slopes, the collection intercept and `theta0`; at most 1.074), and one package warning on `B_output`.
- `design-qfar_K6-sites300-07`: species 10 drifting (second occupancy slope; separation 0.77, Rhat 1.053). Flag: 18 scored diagnostics of species 10 (17 cell probabilities, 1 slope; at most 1.055); no package warning.
- `design-qfar_K6-sites300-09`: species 4, 6, 7 and 9 drifting (separation at most 0.50, Rhat at most 1.050). Flag: 12 original-site cell probabilities of species 9 (at most 1.054); no package warning.
- `design-qfar_K6-field4-09`: species 3 and 6 drifting (first occupancy slope; separation at most 0.71, Rhat at most 1.046). Flag: one package warning on `B_output` only; no scored diagnostic exceeds 1.05.
- Largest separation of any quantity in these fits: 0.77, against the separated threshold of 3. Drift of the same size appears in the unflagged reference fits `-02` and `-03` (Task 1).

## Community 5 scored four ways (Step 2)

Original 100 sites, 1,000 cells; bands by the generating probability (`impact/community5-scores.csv`). Coverage is the phase B 95% interval coverage, the only scorer here that defines one.

- pr11 fit, all four chains pooled (as published): MAE 21.17; signed error below 0.2 +17.00, above 0.8 -28.41; coverage 0.808; largest cell Rhat 1.74.
- pr11 chains 1 and 3 (near-truth): MAE 19.16; +12.46 and -25.52; coverage 0.803; largest cell Rhat over these two chains 1.04.
- pr11 chains 2 and 4 (mirror): MAE 23.71; +21.54 and -31.30; coverage 0.751; 1.01.
- Extended run, 16 chains pooled: MAE 19.17; +12.44 and -25.75; coverage 0.816; largest cell Rhat 1.006.

## Affected comparisons of the current-code report

Each item gives the published value, then the value with community 5 from the near-truth chains, the mirror chains and the extended run, then the published observed-chain sensitivity range. Only the 300-site fit of community 5 is replaced.

- Overall MAE, 300 sites, high contamination ("Using 300 sites gives ... 16.97"): 16.97; 16.77, 17.22, 16.77; range 16.68 to 17.31. Inside the range.
- Baseline to 300 sites ("every community improves"): reduction 3.55 (95% interval 2.31 to 4.79, 10 of 10 improve); 3.75 (10 of 10), 3.30 (9 of 10), 3.75 (10 of 10); range 3.21 to 3.84. Inside; the reduction stays positive throughout. With the mirror chains alone community 5 no longer improves on its baseline, so "every community improves" would fail for that estimate only.
- Baseline to four field samples: 3.29 in every case, since it does not use this fit.
- Four samples to 300 sites ("small advantage ... no clear winner"): 0.26 (-0.90 to 1.42); 0.46 (-0.49 to 1.42), 0.01 (-1.57 to 1.58), 0.46 (-0.49 to 1.41); range -0.12 to 0.61. Inside; the interval spans zero in every case, and 6 of 10 communities favour 300 sites in every case.
- Low-probability bias, 300 sites ("low probabilities remain too high"): +15.34; +14.89, +15.80, +14.89; range +14.72 to +15.89. Inside; still positive.
- High-probability bias, 300 sites ("high probabilities remain too low"): -18.82; -18.53, -19.11, -18.55; range -19.37 to -18.29. Inside; still negative.
- Current minus archived code, 300 sites, from the section "How much did probability accuracy change?" (every interval includes zero): 0.249 (-0.197 to 0.695); 0.047 (-0.140 to 0.235), 0.503 (-0.495 to 1.501), 0.048 (-0.138 to 0.235). The interval still includes zero; no observed-chain range was published for this change.

## Intercept-prior phase B (Step 3)

The extended run's phase B scores replace the control's community 5 at high contamination, and that control is counted as unflagged. Descriptive only; the study's files are unchanged.

- Criterion 4, high contamination: control flagged fits 3 become 2, against 4, 9 and 10 for SD 2, 3 and 5, so it fails for every SD, by a wider margin. Low contamination is unchanged (control 1; SDs 3, 7 and 9) and fails for every SD. The FAIL verdict therefore cannot change.
- Criterion 1: the control's pooled low band goes from 14.88 to 14.65; SD 2, 3 and 5 stay better (13.21, 12.41, 12.11). SD 2 improves in 9 of 10 communities instead of 10 (7 needed). Passes for every SD, as before.
- Criterion 2, high contamination: the control's bands become 14.89 (low), 1.75 (middle) and 18.55 (high), and its MAE 16.77 instead of 16.97. SD 3's MAE is then 0.68 worse, over the 0.5 limit (was 0.48), so that component now fails; SD 3 already failed criterion 2 on the middle band. SD 2's MAE difference becomes +0.19 (was -0.01) and SD 5's +0.93 (was +0.73). The criterion-level results are unchanged: pass for SD 2, fail for SD 3 and 5.
- Criterion 3: the control's coverage changes by at most 0.0031 (high band 0.716 to 0.719); passes for every SD, as before.
- The new arms' own community 5 fits at high contamination are flagged too (largest element Rhat 1.54, 1.28 and 1.26 for SD 2, 3 and 5) and were not replaced, so this is not a like-for-like converged comparison.

## Checks

28 checks pass (`impact/reproduction.csv`). The published community 5 rows of `community-scores.csv` are reproduced from the saved pr11 fit within 4.5e-12 (every metric and numeric column), its band errors within 4.4e-16, its per-chain scores within 9.4e-16, and its phase B control scores within 4.3e-14. Every affected published summary is recomputed from the published scores and from the reproduced scores within 5e-14 points, and the phase B gate detail within 6.0e-14 with identical pass values. Each substitution leaves every other row unchanged. The extended chains reproduce their stored posterior means within 1.2e-14.

## Interpretation (not a result)

- The published pooled estimate of community 5 averaged two chains in each mode. The converged extended run agrees with the near-truth chains to within 0.02 points of MAE, so the published value overstated community 5's error by about 2.0 points, and the ten-community mean of the high-contamination 300-site design by about 0.20. Correcting it slightly favours 300 sites and removes most of the apparent code-version change in that arm, but changes no conclusion.
- In phase B a converged control is slightly more accurate, which makes the wider priors look slightly worse; it strengthens rather than weakens the FAIL verdict.

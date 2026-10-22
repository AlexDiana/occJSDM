# Validation record

Completed 28 September 2026. These checks validate the implementation and saved calculations. They do not establish unbiased recovery or nominal posterior coverage.

## Frozen production and fits

Production source was archived from `d3d710e406c5b022df9cb84b7e29e783b4b225c4` and installed in the study's dedicated library. Every fit records production/library/fitting-script hashes, the input hash, starting RNG state, MCMC budget and warnings. The generator, runner and scorer were unchanged throughout the 81 initial and 26 longer fits. Each worker used one sampler thread; at most four fitting workers ran concurrently.

Seven short pilots checked dimensions, truth alignment, reconstruction and runtime. They are excluded from scientific recovery results. A pilot-only checker defect was corrected before the initial batch: absent stored binary probabilities now produce an explicit unavailable check instead of a vacuous maximum over an empty vector. Repeated pilot fits and warnings were bit-identical before and after that checker correction; the original pilots remain archived.

## Automated checks

| Check | Result |
| --- | --- |
| Frozen installed-package suite, with NOT_CRAN=true | 772 expectations passed; zero failures or test warnings; one explicit opt-in full coverage-grid skip |
| Spatial subset of package tests | 88 expectations passed |
| Study generator/scorer tests | 38 expectations passed |
| Selection/aggregation/sensitivity tests | 24 expectations passed |
| Observed-chain/range diagnostic tests | 24 expectations passed |
| Final research-test rerun, 28 September | All 86 expectations passed; zero failures, test warnings or skips |
| Independent saved-draw audit | All 81 selected fits passed; largest numerical discrepancy 3.524958e-13 |
| Observed-chain analysis | All 81 selected result hashes checked; pooled values reproduce final aggregate and paired estimates |

The independent audit does not source the study generator or scorer. It reconstructs all selected posterior probability draws through native spatial projection, independently checks interval endpoints, containment and group scores, and checks fit/result identities and actual chain dimensions. Separate scorer checks compare independently built bases at every range-grid value against native package matrices and compare selected predictor draws.

The research tests exercise reproducible rare-species truth, paired observations, probability-scale averaging, absent binary summaries, cancellation of signed errors, fixed-range stratification, exact MCMC protocols, community-level longer-run changes, observed-chain extrema and constant/differing range-chain cases. Installed dependencies emit a build-version notice for R 4.5.2; these notices are separate from test failures or test warnings.

## Reporting and independent review

The final HTML report renders from the compact result bundle. Direct parsing confirms 11 second-level sections, 11 tables and four nonempty embedded PNG figures; all sections are eligible for its depth-two floating table of contents. The report requires 81 selected fits, every prescribed longer phase and independent-audit hashes matching the selected results. Negative checks on the final report reject wrong audit hashes, missing hashes and incomplete longer phases. An export-path issue found during final inspection was fixed by preserving absolute figure paths during knitting; all four final HTML images are verified base64-embedded rather than unresolved relative paths.

All four final scientific PNGs were inspected directly. The spatial-range figure was adjusted to include the full range of generating truths. Browser policy blocked the local file URL, so browser-based visual inspection was unavailable; no alternate browser or local-server workaround was used. This limitation does not affect the completed HTML render or direct figure inspection.

An independent reviewer recomputed the aggregate results to within 4.44e-15, checked all 81 audit identities and confirmed the reported errors, paired comparisons, sensitivities, prior context and coverage caveats. Two minor wording/rounding corrections were applied. No important scientific reporting issue remained.

One native convergence warning remains for `range6-rep03-high-k050`. The warning and diagnostic disagreement are preserved, all communities remain included, and initial/longer and observed-chain sensitivity are reported. No additional fits were selected after seeing scientific outcomes. No release bias criterion is cleared by this validation record.

Logs and full provenance remain in the raw study archive documented in [results/README.md](results/README.md). The report's compact tables, plots, hash records and audit results are versioned with the reproduction scripts.

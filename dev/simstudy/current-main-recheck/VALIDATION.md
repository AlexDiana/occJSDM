# Validation record

The 27 September 2026 integration check uses PR #11 with production revision `2a75bf1` merged into it. Production files under `R/` and `src/`, bundled example data, and saved teaching data match that revision. The changes here concern simulation scoring, research evidence and documentation.

The simulation-truth regression file passes 27 expectations. Its checks cover read-event versus positive-read truth, numeric collection-covariate transformations, current numeric-trait transformations, absent parameter blocks and ambiguous metadata mappings. The new trait-scale assertions failed before the correction and passed afterward. Independent review also checked the trait contribution identity on a current fitted object.

Full `devtools::check()` completed with zero errors, three existing warnings and six notes, including construction and rebuilding of package vignettes. Installed-package tests passed 769 expectations with zero failures or test warnings and three skips. The skips cover opt-in coverage and source-only checks.

The warning text matches the already checked production baseline after normalizing elapsed timings and audit-log paths:

- The R header requests an unknown compiler warning option, `-Wfixed-enum-extension`.
- `predictNewSites.Rd` does not document `verbose`.
- `src/Makevars` contains GNU extensions.

Existing notes remain; this is not a claim of CRAN readiness or completed beta validation. The check and baseline comparison are saved in the execution archive as `package-check.log`, `package-check-result.rds`, `verify-package-check.R` and `verify-package-check.log`. The installed test log is `package-check/occJSDM.Rcheck/tests/testthat.Rout`.

Seven short setup pilots were independently reconstructed before the main comparison, checking 28 original-site probability-band scores. The largest posterior-mean discrepancy was `1.332e-15`. Their short-chain diagnostics are not scientific results. The final selected-fit verification, provenance, convergence flags and initial-versus-longer sensitivity are retained in `results/` and described in the [report](REPORT.md).

All 110 initial fits and 51 prespecified longer fits completed. Independent reconstruction of all 110 selected fits and 440 original-site probability-band scores passed, with a maximum posterior-mean discrepancy of `2.597922e-14`. All 20 production, library and fitting/scoring fingerprints were unchanged, and the historical result directories have no changes. The selected fits retain five diagnostic flags and nine fitting warning messages; numerical reconstruction is not a convergence claim.

The observed-chain sensitivity check reconstructs each chain in all five flagged fits, verifies its pooled result against the independently checked mean, and retains all ten communities in every comparison. It exports 100 chain/group rows, including the pooled reference, plus ten-community MAE and signed-bias ranges and design-comparison ranges. Independent review checked its reconstruction and extremum calculations; the report explicitly excludes a confidence-interval or convergence interpretation. All three final figures were visually inspected.

Independent review checked all 110 archived input mappings and the known-site-factor control, then reviewed the selection, verification and summary code. Findings about preserving infinite diagnostic values, actual job identities, MCMC schedules, missing probability bands and unequal rare-species denominators were corrected before exporting the final results. No archived fit or result was overwritten.

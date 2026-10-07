# Post-beta package maintenance, 7 October 2026

The user approved the three maintenance tasks in TODO.md: review the local simulation archive, retire remaining dead package functions and then declare genuine data-masked columns. This report covers the code work; the separate [archive review](../archive-review-20261007/REPORT.md) covers retention and provenance. The base revision is `0c5b1ff19c4de43be329fbd06730aee24548ca60`.

## Dead-code decisions

The package's R call graph was traced from all 42 public exports, `.onLoad()`, `thinOutput()` and the four previously retained experimental functions. Both function calls and variable references were followed, including functions supplied as callbacks. Repository searches covered tests, vignettes, development scripts and MCP code outside archived outputs and pinned runtime copies. Generated native wrappers were checked separately from internal C++ callers. These are repository-local caller checks; undocumented external use of internal functions cannot be ruled out.

Twenty-four distinct R helpers, represented by 25 bodies because `buildGrid()` had duplicate definitions, were moved unchanged into three dated files under `deprecated/R/`. They cover the obsolete grid and spline builders, superseded R samplers, unused likelihood/variance-partitioning helpers, `plotCovariateTrend()`, `sampleB_m()` and `computeMinESS()`. [retired-helpers.csv](retired-helpers.csv) records every function, original file/line range, archive location and body SHA256. The moved files are excluded from package builds by the existing `deprecated/` rule.

Thirteen C++ bodies were moved unchanged into two dated files under `deprecated/`: `sample_w_cpp`, `sample_betatheta_cpp_parallel_old`, the four `isPointInBand*` helpers, `findClosestPoint`, `dist_matrix`, `gpCovMatrix`, `K`, its private helper `k_cpp`, `convert_to_correlation` and `XsBs`. Their retired R callers were removed first. Shared helpers, including `sample_beta_nocov_cpp_TS`, remain in production source.

Four native reference bodies stay unchanged but lose unused R wrappers: `sample_w_cim_cipp`, `sample_pq_cpp`, `samplePGvariables` and `sampleB`. The last remains an internal C++ callee of `sample_U_cpp`. `Rcpp::compileAttributes()` regenerated both export files after changing annotations. In total 16 internal R wrappers/native registrations were retired; public exports, model defaults, priors, native seed handling and sampler thread behavior are unchanged.

The following apparently unreachable functions were deliberately preserved:

- `sample_BBsL()` and its native dependencies are exercised by name-based regression tests. `sample_BBsL_cpp`, `sample_BBsL_parallel`, `sampleB_SoR`, `XtOmegaX_SoR` and `sample_betatheta_cpp` retain their wrappers because production code, regression tests or native-algebra research checks use them.
- `loglik_spatialEffect()` is used by `dev/test_sample_l.R`; preserving this research reference avoids invalidating historical investigations.
- `computePredictiveProbs`, `partition_r2`, `returnSpatialEffectMean` and `plotSpatialEffect` remain as agreed in `e90a5b2`, along with their dependencies. `.onLoad()` is an R hook, and `thinOutput()` retains its separate review/retention decision.
- `K2`, `k2_cpp`, `spatEffectMeanCpp`, `KsBproduct`, `setOccJSDMSeed` and all live samplers remain. No existing regression test was removed.

Parsed expressions of every original non-wrapper R function were compared after the move; all match their live or archived counterparts. The 13 archived and four retained native reference bodies match the original SHA256 values. This verifies preservation, not that obsolete experimental code is correct.

## Data-masked columns

After cleanup, 41 distinct names were verified as columns evaluated in ggplot2 aesthetics, dplyr data masks or tidyselect. [data-masked-columns.csv](data-masked-columns.csv) maps each declaration to its source functions. `R/globals.R` declares these names with `utils::globalVariables()`, and DESCRIPTION declares the utils dependency, with an explicit roxygen import regenerated into NAMESPACE. No plotting or fitting expression changed.

Nine actual undefined variables remain visible: `fitModel`, `Ks_new`, `Xs_centers_new`, `Tr`, `X_ord`, `Bs_output`, `Ks`, `Xs_centers` and `list_Xs`. Missing helper calls in the intentionally retained experimental code also remain visible. They were excluded from declarations, so the R-code NOTE is reduced rather than concealed. The user's approval of this task supersedes the historical July decision to leave the NSE NOTE permanently unchanged.

## Validation

Validation is complete. The full source suite and a fresh installed-package suite each pass 1,075 expectations with zero failures or test warnings and one opt-in full-coverage-study skip. The baseline and corrected final `R CMD check --no-manual` runs both have zero errors, three warnings and three NOTEs. Their installed check suites each pass 1,039 expectations with eight expected skips under CRAN/source-availability rules. Both package vignettes build and rebuild successfully. Full baseline and changed-package checks use R 4.5.0 with the installed R 4.5 dependency library. An isolated local launcher/configuration ensures compilation links the R 4.5 library instead of the current R 4.6 framework; no system R installation was altered. The first baseline build crashed in `K2()` before package edits under the mixed-version configuration. The corrected baseline builds both package vignettes successfully.

The unchanged warning sections are the R-header compiler warning, undocumented `verbose` in `predictNewSites` help and GNU extensions in Makevars. The NOTEs cover LICENSE metadata, genuine retained R-code problems and a `.git` file included from the manual worktree; the last is a check-artifact limitation, not a new source-code problem. An initial changed-package check additionally flagged the utils dependency as unimported; the explicit roxygen import removed that new NOTE in the corrected check. No existing warning was hidden or changed.

The public API remains 42 exports; native registrations fall from 38 to 22, with both tested `sample_BBsL` entry points retaining 17 arguments. Rcpp wrapper regeneration is byte-idempotent. Namespace usage diagnostics fall from 195 baseline lines to 178 after cleanup and 54 when declarations are recognised; these line counts include unused-local diagnostics and are not counts of package NOTEs. [validation-summary.json](validation-summary.json) records the results. [tested-source.csv](tested-source.csv) records source hashes: all 21 R/C++/header/namespace/Makevars files match the checked source archive byte-for-byte, and all original DESCRIPTION fields survive build normalisation. The build-ignore file is recorded separately because it is not shipped.

Bulky logs, original function expressions, call-graph records, the source-test result object, build/check trees and the isolated runtime stay under ignored `dev/simstudy/results/post-beta-maintenance-20261007/`. The implementation plan is under `dev/superpowers/plans/2026-10-07-post-beta-maintenance.md`.

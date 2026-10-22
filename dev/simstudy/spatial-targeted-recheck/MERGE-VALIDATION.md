# Integration validation, 28 September 2026

The completed spatial study and rare-species diagnosis were integrated with main `43351bb`, which includes the approved PR #13 output correction. This branch changes research scripts, compact evidence, reports and planning records. All changed paths are excluded by `.Rbuildignore`; installed-package source, help, tests and vignettes are identical to current main. The full package check completed earlier today for that same package input with zero errors, three existing warnings and five notes, including vignette builds and 834 passing installed-package expectations. Its logs remain in `dev/simstudy/results/pr13-integration-20260928/`.

Fresh checks on the combined source pass 837 package expectations (one opt-in coverage-study skip) and all 86 research expectations, with zero failures or test warnings. The diagnostic checks pass explicit latent-state enumeration, independent R/C++ likelihood agreement and adaptive binomial integration. All 41 compact result SHA256 identities and both research/calculation source manifests match. The compact report renders with 11 second-level sections, 11 tables and four embedded PNG figures.

Independent merge review reproduces all 660 study aggregate means within 4.89e-15, all 81 audit identities, the 26 longer selections and single retained convergence flag, and all 28 diagnostic summaries from 504 conditional calculations within 6.39e-16. No critical or important review findings remain. Full audit reproduction still requires the original raw archive; the compact bundle is sufficient for inspecting and rendering the reported results.

The original study remains frozen at `d3d710e`; this integration does not change any saved fit, result, model, prior or sampler. The tightened TODO preserves the spatial-amplitude half-Cauchy comparison before beta and occupancy-intercept widening after beta. PR #14 remains pending Alex. These checks do not establish unbiased recovery, nominal interval coverage or completion of the beta release criteria.

Integration logs, test results and rendered HTML are retained locally in `dev/simstudy/results/spatial-merge-20260928/`.

# Posterior draws accessor, 8 October 2026

Implemented on `codex/posterior-draws`, based on `1953512`. The approved scope is an exported accessor for the six occupancy/detection parameter blocks already covered by `returnConvergenceDiagnostics()`: `beta0_psi`, `beta_psi`, `beta_theta`, `p`, `q` and `theta0`.

`returnPosteriorDraws()` returns fitted-scale values in a labelled four-dimensional array: coefficient or primer, species, saved iteration, chain. Species-only parameters use a length-one first axis. Species and covariates are selected by name; primers by identifier, accepting equivalent numeric and character forms. Selection retains all dimensions and requested ordering. Missing parameter draws and invalid selectors produce errors. `plotTraceplot()` uses array labels by default, preserving explicit label overrides. No sampler, prior or fitting behavior changed.

Lesson 3 now uses the public accessor for collection, primer and field-contamination traces, per-chain means and newer coefficient diagnostics. Help, exports, the pkgdown index and the completed-work record are updated. Existing teaching fits, compact bundles and figure contents are preserved.

## Verification

- Before implementation, the source suite passed 1,075 expectations. New tests failed because the accessor was missing and traceplots ignored array labels.
- The final full source suite passed 1,132 expectations: zero failures, errors or test warnings; one expected opt-in coverage-study skip. The installed-package diagnostics/accessor tests passed 128 expectations with no failures or test warnings.
- An isolated R 4.5.0 install succeeded. Its exported accessor reproduced every saved value and the iteration/chain extents for the default PCR, perfect-observation and longer alternative-prior teaching fits. All five revised optional Lesson 3 examples executed against those saved fits, including direct `posterior` diagnostics. The independent reviewer additionally checked all six parameter blocks in shipped `sampleresults`.
- Lesson 3 rendered to HTML and unwrapped GitHub Markdown. Shared lesson-link checks, generated-help validation, `pkgdown::check_pkgdown()` and `git diff --check` passed. One rerendered PNG differed only by colour-channel rounding of at most 1/255; its unchanged content was visually checked and the committed original restored.
- Independent code review found no functionality issue. Its minor missing roxygen Markdown directive was corrected, and the help links regenerated and verified.

R emitted startup notices that installed `testthat` and `tibble` were built under R 4.5.2; neither produced test warnings. A full `R CMD check` and the opt-in coverage study were not run for this interface change. Task-specific saved-fit verification code and preview outputs are retained in the worktree's ignored `dev/simstudy/results/posterior-draws-20261008/` directory. No teaching fits were regenerated.

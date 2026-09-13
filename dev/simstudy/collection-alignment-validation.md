# Collection-covariate alignment validation

10 September 2026. Branch: `codex/align-collection-covariates`. Implementation ready for Alex's review.

The fix builds collection covariates from canonical typed `(Site, Sample)` rows in the sampler's order. This removes lexical ordering errors for numeric identifiers and collisions between pasted character labels. The public modelling features and prior settings are retained.

**Paired recovery check.** Five simulated datasets in each of the existing `base` (spatial and traits) and `traits_isolated` (traits, no spatial field) scenarios, for 20 fits before/after. Each dataset has 100 sites, 10 species, two samples per site, two primers and three PCR replicates per primer. Fits use two chains, 1,000 burn-in iterations and 1,000 retained iterations per chain. The spatial arm uses 20 knots. Both versions use identical simulated data and fitting seeds, the unchanged C++ backend, default priors including collection-slope variance 2, and one TBB thread.

Generating slopes are multiplied by the sample-level covariate standard deviation before comparison, so truth and estimates are on the fitted scale. Sign-group summaries average species within each dataset and then give each dataset equal weight. The reported bias MCSE is the standard error across dataset-level mean errors, not a posterior interval or a coverage result.

**Spatial and traits.**

- Negative slopes: mean truth -0.993; estimated mean -0.033 before, -0.951 after; mean error after +0.042 (MCSE 0.153).
- Zero slopes: mean truth +0.000; estimated mean +0.077 before, -0.048 after; mean error after -0.048 (MCSE 0.134).
- Positive slopes: mean truth +0.993; estimated mean +0.116 before, +1.053 after; mean error after +0.061 (MCSE 0.107).

Across all species and datasets, mean absolute slope error fell from 0.698 to 0.326 (53% lower).

**Traits, no spatial field.**

- Negative slopes: mean truth -0.991; estimated mean -0.232 before, -1.029 after; mean error after -0.038 (MCSE 0.115).
- Zero slopes: mean truth +0.000; estimated mean -0.058 before, +0.046 after; mean error after +0.046 (MCSE 0.056).
- Positive slopes: mean truth +0.991; estimated mean +0.401 before, +1.065 after; mean error after +0.075 (MCSE 0.134).

Across all species and datasets, mean absolute slope error fell from 0.585 to 0.366 (37% lower).

**Interpretation and remaining gates.** The severe attenuation of nonzero collection slopes disappears in these scenarios. Five datasets per scenario provide a focused diagnostic, not evidence that all parameter biases have been resolved. Individual slope estimates still have appreciable absolute error; see the saved MAE and RMSE summaries. The separate RNG, residual-correlation, spatial, high-`q` and `B0` gates remain open.

- `base`: mean `B0` error changed from -0.302 to -0.219; mean occupancy-probability error from -0.042 to -0.029. These residual offsets still need the planned bias recheck.
- `traits_isolated`: mean `B0` error changed from -0.164 to -0.116; mean occupancy-probability error from -0.026 to -0.019. These residual offsets still need the planned bias recheck.

**Verification.** The full test suite passes 367 assertions, including 70 alignment assertions. Two existing tests are skipped: the opt-in broad coverage study and the known multi-thread reproducibility gap. The new regression tests fail on the original implementation and pass with this fix; they check actual sample/covariate/observation/index pairing, rather than only agreement between shuffled fits. Identifier columns requested as covariates have their own regression test.

The installed-package check finishes with 0 errors, 3 warnings and 5 notes, with its tests passing. Vignettes and the PDF manual were not rebuilt. Warnings concern unchanged compiler/header and unused-function diagnostics, an undocumented `verbose` argument in `predictNewSites`, and GNU extensions in `src/Makevars`. Notes concern existing R code checks, the LICENSE declaration, the worktree's `.git` file, inability to verify network time, and `xcrun_db` temporary files. No source, documentation or Makefile changes were made to those areas.

**Reproduction and provenance.** Run from the package root:

```sh
Rscript --vanilla dev/simstudy/validate_collection_alignment.R 5 dev/simstudy/results/collection-alignment 6c837bc
```

The before revision is `6c837bc2dc60dcafbcc8764cbe8978230db62c79`. The corrected `R/runOccJSDM.R` Git blob is `b2fe3edfed43f81513f67a5029dd2c9b38e63401`. The diagnostic uses `simstudy_seed_for()` unchanged and resets the fitting seed for each before/after pair. It refuses to compare revisions with different C++ source.

The runner writes `slopes.csv`, `slope-summary.csv`, `fit-metrics.csv`, `paired-results.rds` and `provenance.rds`. Provenance records exact before/after R source, the starting revision and diff, MCMC settings, thread count and R session. Sign-group reporting was refined after fitting to use equally weighted dataset means with matching MCSEs; the final summaries were recomputed from saved fits using the reviewed runner's reporting block and checked independently. This reporting change did not require new MCMC.

All 20 fits completed and were saved. Editing the running diagnostic script caused a trailing parse error afterward. A separate post-processing run completed successfully, verified that the saved fitting source matches the final code, and checked the corrected summaries. The committed diagnostic script passes a fresh parse check.

`TODO.md` retains this item as `ALEX TO REVIEW`; the historical Fixed bugs and Completed work records have not been altered.

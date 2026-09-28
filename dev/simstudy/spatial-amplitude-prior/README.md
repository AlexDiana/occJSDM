# Spatial-amplitude half-Cauchy experiment

The [prespecified protocol](PLAN.md), [compute amendment](COMPUTE-AMENDMENT.md) and [pre-outcome estimand amendment](ESTIMAND-AMENDMENT.md) compare the existing inverse-gamma spatial-amplitude prior with an opt-in half-Cauchy(scale 1) prior using nine saved binary communities at 100 spatial support points. Production implementation and fitting protocol were frozen at `4509629`; the old fitting reference is `d3d710e`, and the branch starts from main `74e33a5`. Default priors remain unchanged. This experiment is separate from continuous observation noise and from the post-beta occupancy-intercept-prior task.

**Use the amended median analysis.** In these full-rank binary fits, the unbounded half-Cauchy leaves the posterior amplitude mean infinite and spatial-field means non-integrable. Both priors therefore use pointwise field medians and rank/quantile diagnostics; bounded occupancy probabilities retain their posterior means. The original mean calculations remain archived, but their extension gate is mathematically inapplicable. No fitting or prior change was made in response to this discovery.

## Files and archives

- `metrics.R` reconstructs field draws in bounded blocks and checks amplitude and per-species field summaries.
- `run.R` reads the original immutable input archive, uses a separately installed experimental package, and saves complete fits before scoring. It fingerprints both installed R code and the compiled library.
- `analysis.R` and `summarise.R` preserve the original selection/calculation workflow. Their mean-based recovery outcomes and gate are superseded by the estimand amendment.
- `continue.R` uses available worker slots to run the finite list of required longer fits and the separate initialization check, then stops at the binary comparison. The original cap was four; Doug authorized a total cap of eight during execution. Each sampler remains single-threaded. It does not automatically start two-stage fits.
- `robust.R`, `rescore.R` and `robust-summarise.R` implement the separately versioned median/quantile analysis, retain every previously triggered longer run, and apply the amended diagnostics symmetrically. Immutable legacy results and full fits are fingerprinted.
- `verify-robust.R` audits every selected fit and the separate start check through the native spatial projection. It independently checks probabilities, medians, interval containment and diagnostic calculations. `extension-decision.csv` adds completed auditing with zero unresolved diagnostic threshold crossings to the amended statistical gate.
- `plot.R` and `report.Rmd` produce scientific figures and an HTML report from audited selected results.

Raw archive for this run: `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-amplitude-20260928`. Original input/control archive: `/Users/douglasyu/src/occJSDM/dev/simstudy/results/spatial-targeted-20260927`. Both are ignored by Git. Keep them for full reproduction; compact evidence will not substitute for all saved posterior draws.

## Reproduction

Use absolute paths for `REPO`, `STUDY` and `REFERENCE`. `STUDY/source` must contain the frozen production files and `STUDY/library` the separately installed experimental package. Keep the installed library unchanged throughout the run. Run from this directory when invoking the standalone research tests or the audit.

```sh
Rscript test-metrics.R
Rscript test-analysis.R
Rscript test-execution.R
Rscript test-robust.R
Rscript run.R --repo=REPO --study=STUDY --reference=REFERENCE --mode=baseline --prior=inverse_gamma --workers=4
Rscript run.R --repo=REPO --study=STUDY --reference=REFERENCE --mode=initial --workers=4
Rscript continue.R REPO STUDY REFERENCE 8
Rscript rescore-checked.R REPO STUDY 4
Rscript robust-summarise.R REPO STUDY select
Rscript robust-summarise.R REPO STUDY initial
Rscript robust-summarise.R REPO STUDY final
Rscript verify-robust.R STUDY STUDY/robust-v1/summary-binary-final 4
Rscript plot.R REPO STUDY/robust-v1/summary-binary-final
python3 export.py STUDY REPO
```

Run up to four rescoring/audit workers after fitting ends, or one alongside fits. Inspect the amended `select` output before final summarization: any additional initial diagnostic flag requires the same single paired longer schedule and then rescoring. Do not rerun or overwrite existing fits. Render `report.Rmd` with `params = list(summary = 'STUDY/robust-v1/summary-binary-final')` or the compact `results` directory, using an isolated evaluation environment. Check `extension-decision.csv` before undertaking the matching low/high-contamination comparison. A failed field-recovery, convergence or audit criterion stops that extension; no prior-default change is implied.

The checked rescoring wrapper verifies the frozen hashes in `scoring-dependencies.csv` before and after extraction and binds every cached result to those dependencies. It complements the core scorer's own hashes without changing previously frozen scoring code. Summary selection requires this dependency manifest as well as matching input, fit, legacy-result and scorer fingerprints.

## Numerical diagnostic audit

Two reconstructions can agree to machine precision yet produce slightly different folded-rank Rhat because the two observations surrounding an even-sample median have theoretically equal absolute deviations; floating-point arithmetic can assign different ranks to that pair. The independent pilot showed field/trace/probability differences below 8.1e-15, exact stored-diagnostic reproduction, and no Rhat threshold crossing. The audit therefore checks exact diagnostics from the saved traces separately from numerical native reconstruction. It records native Rhat differences, changed folded ranks and any crossing of 1.05. Review a crossing conservatively under the longer-sampling rule; do not alter or round frozen production traces after inspecting results. Any unresolved threshold crossing blocks the extension pending review.

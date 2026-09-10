---
title: "TODO"
output: html_document
---

```{=html}
<!-- DO NOT hard-wrap this file. Write each paragraph as a single line.

     There is deliberately no `editor_options: markdown: wrap:` block, and
     `occJSDM.Rproj` sets `MarkdownWrap: None`. Doug edits these docs in
     RStudio's Visual mode, which rewrites the file to canonical form on
     every save. While the wrap was set to 72 columns, anything written by
     hand or by an agent got reflowed on the next open, producing a diff
     nobody authored, over and over.

     Measured 30 July 2026: RStudio's canonicalisation cannot be
     reproduced outside RStudio. Feeding RStudio's own output back through
     pandoc changes 234 of 768 lines, across several flag combinations. So
     matching it by hand or by hook is not achievable, and the fix is to
     remove the line-breaking step entirely: with no wrapping there is no
     algorithm left to mismatch.

     Three related rules:
     - NO EM-DASHES. Use a plain double hyphen. Em-dashes were the single
       largest driver of reflow churn under the old setting, because
       pandoc writes them 2 characters wider and every later wrap point
       shifted.
     - NO PIPE TABLES in this file. Use bullet lists instead. RStudio
       recomputes a table's separator row (|----|----|) from the cell
       contents and the available width, and it switches between
       padding to content width and redistributing across the line at a
       width threshold. So any edit to any cell can rewrite the
       separators, and 873f419 is an example. Lists have no
       width-derived formatting and are stable. PLAN.md keeps its
       tables, and accepts the churn, because the data there earns them.
     - Keep inline `code spans` short. Under the old 72-column wrap, a
       span landing on the boundary could be mangled: in 2242b34 that
       turned `* 2` into `- 2` in the logDetKuu note, silently asserting
       something false. Less dangerous now, but still good practice.

     What this file is FOR: a list of to-dos with enough context to act
     on, addressed to Alex. It is not an investigation log. Evidence and
     superseded reasoning go to AGENTS.md; completed work goes to the
     Fixed bugs section below, which is the record. If an item is
     fixed, write the Fixed bugs entry and DELETE the item; do not leave
     a struck-through stub. -->
```

# **v0.1.0-beta Public release**

**Release criterion agreed 10 September 2026:** retain all advertised modelling features, fix incorrect or materially biased point estimates, and allow undercoverage or overcoverage to wait. This triage assumes a GitHub beta. The CRAN submission and paper work can follow later.

**Beta gates: review the collection-covariate fix, complete the other three code workstreams, then run the targeted bias and release checks.** The `B0` and high-`q` findings remain conditional release gates: reassess them after the code fixes before deciding whether any prior change is needed.

## Required code work (Alex, with Doug validating)

1. **ALEX TO REVIEW: align collection covariates with the samples used by the sampler.** Implemented on branch `codex/align-collection-covariates`. In `R/runOccJSDM.R`, one canonical table of typed `(Site, Sample)` pairs now supplies `M` and `X_theta`, in the same order as `P`, `K` and the latent-state indices. Collection covariates are read directly from those rows, preserving identifier columns when requested as covariates. The pasted `SiteSample` label no longer controls covariate grouping or ordering.

    **Review:** check the source change and `test-collection-alignment.R`. Its 70 assertions verify actual covariate/observation/index pairing for numeric sites `1`, `2`, `10`, global and within-site sample IDs, shuffled rows, unequal primer/PCR replication, colliding character labels, categorical covariates, identifier covariates and intercept-only models. The original implementation fails the pairing checks; the revised implementation passes.

    **Before closing:** review the paired recovery results in `dev/simstudy/collection-alignment-validation.md`, then record Alex's decision before moving this item to *Fixed bugs*. The reproducible runner is `dev/simstudy/validate_collection_alignment.R`. It compares negative, zero and positive slopes on the fitted covariate scale using unchanged priors and backend code. This alignment fix does not close the separate `B0`, high-`q`, RNG, correlation or spatial gates. Any material slope bias exposed by subsequent checks remains a beta blocker; interval coverage alone does not.

2. **Ensure random draws on the public fitting path are safe and independent.** `src/rng.h`, `src/functions.cpp`, `src/jsdm.cpp`, and `runOccJSDM()`. Combine the two existing RNG items: TBB workers obtain the same thread ID through `omp_get_thread_num()`, duplicating random streams; `sample_beta_nocov_cpp_TS()` also calls the R-RNG-based `sample_beta_cpp()` inside a worker, creating a race. The shared-stream problem also reaches the live JSDM `PG_Worker`, so changing the collection sampler alone is insufficient.

    **Minimum beta change:** run all RNG-bearing sampling steps serially on the main thread, using the existing serial implementations where available. This retains every modelling feature; full parallel sampling can follow later. Deterministic worker calculations can remain parallel. Enforce this in the package's execution path, including when the user requests multiple TBB threads. `options(mc.cores = 1)` does not control these workers, and a README instruction alone is insufficient.

    **Alternative:** implement independent streams keyed to the sampling step, chain, iteration and species/element, and remove all R RNG calls from worker threads. Do not simply reseed a generator from `(base_seed, species)` on every call, which would repeat its draws each iteration. Restoring the `_TS` leaf alone fixes only the race, not the duplicated streams.

    **Done when:** audit every live RNG call path; same-seed repeated fits reproduce, consecutive fits consume different draws, and the package cannot enter the unsafe path under a multi-thread environment setting. If retaining parallel random draws, also verify stream independence and equivalence of the sampled distributions. Reproducibility alone does not establish independence.

3. **Preserve residual species correlations during factor reparameterisation.** `reparamFactorModel()` in `R/jsdmfun.R`, its calls in `runOccJSDM()`, and the correlation outputs in `R/output.R`. The current transform preserves `U %*% L` but rescales factors unequally, so `cov2cor(crossprod(L))` changes. A direct check of current code changes one species pair from +0.316 to -0.316 without changing the linear predictor. This is an algebraic output defect; the withdrawn coverage argument is unnecessary, and even the signs are not a safe workaround.

    **Do:** use an orthogonal QR rotation with any needed sign convention, without the diagonal magnitude rescaling. Handle the one-factor case with a sign change rather than division by a loading's magnitude. Review both the `U`/`L` and `A`/`C` call sites. If retaining a non-orthogonal transform instead, carry its transformed latent covariance into every affected output and prediction; do not substitute an empirical covariance of fitted site scores without deriving the target.

    **Done when:** deterministic checks preserve both the linear predictor and the residual correlation matrix before/after transformation for one and several factors, including unequal factor scales and mixed-sign correlations. Ordination, trait products and prediction still work. Add a focused recovery check with non-degenerate true correlations; the old nearly full-width intervals cannot validate correlation magnitudes.

4. **Resolve the spatial length-scale boundary behavior on the full fitting path.** `update_jSDMcoef()`, `computePsiCoef()` and `precomputeSORmatrices()` in `R/jsdmfun.R`, together with the spatial coefficient update in `src/jsdm.cpp`. The recorded failure is that different generating ranges lead to the largest grid value. *Fixed bugs* 48 closed the isolated `sample_ls()` investigation; a successful test using a supplied GP draw does not validate the inputs produced during an actual fit.

    **Do:** reproduce the full-fit symptom with the current simulator and a nonzero spatial field. Trace the `SE` passed to `sample_ls()` after the coefficient updates, and verify that its coordinates, scale and covariance representation match `Ks_all`, `Lm1_grid` and `logDetKuu_grid`. Fix the demonstrated mismatch in field construction or scoring. Do not repeat the already ruled-out amplitude and log-determinant experiments without new evidence.

    **Done when:** several sufficiently informative simulated datasets with distinct interior grid ranges no longer all select the same upper boundary, and spatial-field and occupancy point estimates recover their generating pattern and level. Use the coordinate scale actually fitted and a fixed knot count. A moving chain alone is insufficient; exact range recovery in every replicate and nominal interval coverage are not beta requirements. Keep spatial fitting available.

## Required bias recheck and release preparation (Doug and Alex)

- **Recheck `B0` and the collection-prior trade after fixing sample alignment.** The existing variance comparison gives `B0` bias -0.160, -0.106 and -0.044 at slope-prior variances 2, 0.5 and 0.1. However, the saved nonzero collection slopes are more attenuated at the tighter priors, which pooled signed bias hid. Do not change the default to 0.1 merely to improve `B0`. Keep Alex's existing choice of 2 while diagnosing the data-pairing defect, then compare `B0`, sign-stratified collection slopes and occupancy probabilities on paired data. If material bias remains, resolve its cause or validate a prior specification that improves the affected estimates without transferring the bias elsewhere. *Fixed bugs* 46 remains the historical decision; it is not evidence of zero bias under today's release criterion.

- **Recheck the biased half of the high-`q` item; defer its coverage-only half.** Saved results with true `q` in 0.15-0.30 give mean bias -0.0365 at K = 3 and -0.0314 at K = 30. That is point-estimate bias, not just undercoverage. After the alignment/RNG fixes, rerun the existing paired `qnear_K3`, `qnear_K30`, `qfar_K3` and `qfar_K30` cells. If the offset persists, distinguish prior sensitivity from a remaining likelihood/latent-state defect and validate the remedy using `q`, `p` and occupancy estimates. Keep an identifiable true-/false-detection model; do not flatten the priors indiscriminately. Correct the old causal claim: the current `sample_pq_cpp_parallel()` computes Beta parameters and draws `p`/`q` via `R::rbeta()` on the main thread. They are not directly Polya-Gamma updates, so similar coverage patterns do not establish a common PG bug. The near-prior arm's undercoverage can wait if point estimates are acceptable.

- **Use a focused bias gate, not a new nominal-coverage study.** First validate row alignment, coefficient/truth scale and element mapping. Reuse paired seeds and existing scenarios for the changes above; preserve `simstudy_seed()`. Check mean errors separately by true slope sign and by probability/prevalence range, alongside absolute error and RMSE. Inspect occupancy probability levels as well as correlations, which cannot detect a level shift. Choose practically acceptable errors on the probability scales before evaluating the rerun, with Monte Carlo uncertainty assessed across replicate datasets. A pooled signed mean near zero is insufficient. If the targeted checks expose further material bias, it remains a beta blocker; broad paper-quality coverage experiments can wait.

- **Validate and refresh the version actually released.** Run the normal tests and an installed-package check after the code changes, fixing installation failures, crashes and new functional regressions. Refresh shipped `sampleresults` and any dependent vignette numbers/plots against the corrected code and matching data. Record the code revision, seeds, priors and thread setting used for the bias checks. Review the numerical changes before tagging the beta. Existing cosmetic check NOTEs and a larger coverage study are outside this beta gate.

- **Update public limitations and the announcement before release.** Put the retained coverage limitations and prior-sensitivity guidance in the README and fitting documentation as well as the announcement; internal TODO notes are excluded from the public documentation. Describe the actual tested beta behavior. Remove the old claims that correlation signs are trustworthy, that near-prior `q` intervals are calibrated, and that moderate K removes the high-`q` bias. `useSpatField` is a simulator setting; the fitting API selects the spatial field through `spatCovariates`. Retain all advertised modelling features once their gates above clear. No `NEWS.md` is required for this beta, per the existing decision.

## Announcement draft (Doug; finalise after the release checks)

> Subject: occJSDM beta: joint species distribution modelling with two-stage eDNA detection
>
> We are releasing the beta of occJSDM, an R package combining joint species distribution modelling with the two-stage eDNA occupancy model of Ji et al. (2025). It estimates false-negative and false-positive detection at field and lab stages, with primer-specific lab rates.
>
> Features include environmental and collection covariates, species traits, nonlinear environmental responses, spatial effects, ordination, residual species correlations, variance partitioning, and prediction at new sites. Simpler study designs support classical occupancy and JSDM-only models.
>
> This is beta software. Credible intervals can under- or overcover; nominal interval coverage has not been established across all supported designs. False-positive models require informative assumptions, and users should examine prior sensitivity, especially with weak detection or higher contamination rates. The README and vignettes describe the tested settings and remaining limitations.
>
> Installation and examples: <https://github.com/AlexDiana/occJSDM>. Feedback and bug reports are welcome.

# **Future work after beta**

Every outstanding item from the previous TODO is accounted for below or in the required work above. Deferred means still open; it does not mean fixed. The previous crashes/API section had no open items.

## Interval calibration

- **Collection-slope and near-prior `q` undercoverage:** defer the remaining interval-width investigation once the point-estimate checks above clear. Keep the two mechanisms separate until evidence connects them. The collection row-alignment defect itself stays in the beta list.
- **`B0` undercoverage in the continuous model:** defer the second clean-configuration experiment. Its recorded concern is interval width with negligible bias. This is separate from `B0` bias in occupancy/two-stage fits, which is rechecked above.
- **`theta0` overcoverage:** defer. The recorded point-estimate bias fell substantially; excess interval width is allowed for beta. If revisited, investigate the previously untested collection-prior mean change rather than tightening `theta0`'s own prior to force a coverage target.

## Review and maintenance

- **Review `thinOutput()`:** defer Alex's review and the keep/delete decision while it remains internal and unused by the public workflow. Its rewrite is still pending review, not newly added to *Fixed bugs*. It thins the second-to-last array axis and preserves posterior means, nested JSDM arrays and scalar WAIC; tests distinguish iteration counts from site/species dimensions. Refit the shipped example rather than relying on this helper for release preparation.
- **Move remaining dead functions to `deprecated/`:** defer. Recheck current callers before removal; the old list incorrectly includes the now-live `sample_BBsL_cpp()`. Keep the four previously retained functions unless Alex revisits that decision. Change C++ export annotations and regenerate wrappers rather than editing `RcppExports.R` manually. `.onLoad()` is called by R and is not dead code.
- **`globalVariables()` for data-masked columns:** defer until dead-code cleanup is done. Only declare genuine NSE column names; do not hide undefined variables in executable code.
- **Repair `sample_rnb()`:** defer with count-data support. Supply the current size vector explicitly, implement the intended prior and choose a usable proposal scale before wiring it into fitting. Count models are currently rejected.
- **Inert simulator `sigma_ts` and `sigma_bs` inputs:** defer API removal or implementation. Keep their limitations explicit in simulation documentation and do not treat the supplied `sigma_bs` as a generated truth for validating the fitted spatial variance. `sigma_ts` is unused; the simulated residual spatial coefficients are zero, although `sigma_bs` is live in fitting. This is not evidence that the fitting-side parameter is biased.

## Paper and broader validation

- **Categorical species traits:** defer implementation. State the currently supported trait encoding; retaining species-trait modelling does not require introducing a new encoding in this beta.
- **Reproduce the Ecology Letters analyses:** defer the full reproduction and decision about including it in the repository. The beta needs the targeted checks above.
- **Repeat the complete simulation grid after fixes:** defer the comprehensive paper run. For a deliberate production-grid run, specify `base,binary,d_overfit,d_underfit,low_information,occupancy,primers_3,spatial_isolated,species_20,traits_isolated`; a bare runner invocation also selects additional experimental cells.
- **Choose the paper's replicate count:** defer the R = 200-500 calibration study and any claim of nominal coverage. The existing R = 100 study remains a historical baseline.
- **Simulation-study presentation:** already decided: regenerate the pkgdown validation article from `validation-data.rds`. No new presentation decision is required; update its data before quoting new results.
- **Publish the pkgdown site:** defer as a beta dependency; the README and vignettes can serve beta users. The site scaffolding is built, but publication requires the manual workflow and Alex's Pages configuration change. Include current limitations wherever documentation is published.

## Performance and parallelisation

All speed work can wait once the unsafe RNG path is removed from beta. Preserve the live serial `sample_BBsL_cpp()` call unless a separately validated change replaces it.

- **Alternative Polya-Gamma sampler:** defer the performance experiment. Profile again after the correctness fixes; the historical hotspot was the collection-covariate update.
- **Parallel chains:** defer. When implemented, use portable PSOCK workers, independent process streams, one safe sampler per process, and correctly merge per-chain WAIC accumulators and posterior means.
- **Remove `.onLoad()`'s global `mc.cores` setting:** defer the session-state/CRAN cleanup. This setting provides no protection for the current TBB RNG problem; that protection must be implemented in the beta work above.
- **Repeated `computePsiCoef()` calls:** defer optimisation. The three current calls use changed coefficients or factors; there is no demonstrated redundant-call saving.
- **Precompute/fuse the `c_imk` update:** defer moving invariant work out of the chain loop and avoiding the repeated `w_all` gather.
- **Optional, cheaper WAIC:** defer the new option and likelihood optimisation; preserve the current estimator's semantics.
- **Vectorised initial values:** defer replacing nested `w`/`z` initialisation loops with grouped reductions.
- **Selective posterior storage:** defer a `keep` argument and allocation/thinning improvements. Some large arrays are genuinely filled and used; classify them before suppressing storage.
- **Avoid repeated matrix inversions:** defer wiring in and validating the precision/Cholesky implementation. The unused `_TS_opt` name does not make it worker-safe: its normal draw still uses `arma::randn()`.
- **Fully reproducible parallel random draws:** defer the full stream redesign if beta uses main-thread serial draws. Key streams or persistent state by all relevant sampling identifiers so iterations, chains and parameter blocks cannot repeat one another's stream. Re-enable parallel sampling only after resolving both the race and the identical-stream defect, including the JSDM PG path.

## Future modelling features

- **Improved model-selection criterion:** defer replacing the current criterion's tendency to overfit; avoid implying that a selected model is necessarily the true model.
- **Count-data models:** defer, including the `sample_rnb()` work above.
- **Source-sink inference scenario:** defer a dedicated simulation with opposing environmental and spatial effects.
- **Separate environmental, spatial and latent-factor contributions:** defer restricted/orthogonalised alternatives intended to keep environmental effects stable when additional components are added. This is a modelling extension, separate from the correlation correction required for beta.
- **Site-level variance partitioning:** defer the new output and the naming decision between variance and variation partitioning.
- **Spike-in-based abundance changes:** defer the eDNAPlus extension.
- **Regularised covariate selection:** defer shrinkage/selection for environmental and spatial covariates, including high-dimensional embeddings.
- **Nonlinear-response simulation:** defer exposing the simulator's spline switch and defining response-curve recovery statistics. The existing fitting option `listParams$splineVars` should be described accurately in the beta fitting documentation; a new simulation interface is not needed to retain that fitting feature.
- **Further speedups:** covered by the performance list above; no separate beta task.

## Historical record

The *Fixed bugs* and *Completed work* sections below are retained unchanged. Their old group letters, item numbers and status wording refer to the pre-triage document. Use the subject-based decisions above for current release status; deferred review items have not been promoted into the fixed record.


# **Fixed bugs**

This is the authoritative record of audit items that have been resolved. Items 1-10 were fixed in commit `b7b6aa2` (26 July 2026); items 11-20 followed on 27 July 2026. Line numbers point at the *post-fix* code.

Items 16 and 18 are marked **partially fixed**: the crash in each is gone, but part of the original defect remains, and the remainder is tracked as a live item in group B above. Nothing in this list has a test attached -- every entry was verified by reading the code, which is what the *MEE paper* testing item is meant to address.

**Commit b7b6aa2:**

1.  ~~**`data$info` is re-sorted but `data$OTU` is not.**~~ **FIXED.** The sort permutation is now applied to `y`/`OTU` along with `data_info`. `R/runOccJSDM.R:476-558`.

2.  ~~**`sample_pq_cpp()`'s catch-all `else` inflates the Stage 2 FP counts.**~~ **FIXED.** The C++ version now counts the four cases independently, gated on `primerIdx == l`, matching the R reference version. `src/functions.cpp:652-668`.

3.  ~~**`sampleBuniv()` drops the residual precision (operator precedence).**~~ **FIXED.** The expression now correctly evaluates to `1.0 / (sigma * sigma)`. `src/jsdm.cpp:579`.

4.  ~~**`Bs <- list_params$B` in `update_jSDMcoef()`.**~~ **FIXED.** Now correctly reads `list_params$Bs`. `R/jsdmfun.R:940`.

5.  ~~**`sigma_bs` is sampled but never used.**~~ **FIXED.** The spatial block of the prior covariance now correctly uses `sigma_bs^2`. `R/jsdmfun.R:711-713`.

6.  ~~**`useBiotic` is ignored in `computeNewOutputs()`.**~~ **FIXED.** Now tests `useBiotic` where it should (in addition to `useEnvCov`). `src/jsdm.cpp:501`.

7.  ~~**`simulateData()` adds the spatial field *after* the responses are drawn.**~~ **FIXED.** The spatial field is now added before drawing responses. `R/jsdmfun.R:429-461`.

8.  ~~**`tau_output` is returned entirely `NA`.**~~ **FIXED.** The `tau` parameter is now saved in the MCMC loop. `R/runOccJSDM.R:934-1236`.

9.  ~~**WAIC running means divide by the raw iteration counter.**~~ **FIXED.** A dedicated `currentWAICiter` counter is initialised to 1 (`R/runOccJSDM.R:831`) and incremented only inside the WAIC block (`:1200`), so the running means now divide by the number of WAIC accumulations rather than by the raw MCMC iteration index. A guard, `if (numIters != (currentWAICiter - 1)) stop("Current wAIC iter wrong")` (`:1251`), asserts the two agree. Confirmed introduced by `b7b6aa2` via `git log -S currentWAICiter`. *(This fix was not recorded when the other `b7b6aa2` fixes were logged; recovered on 27 July while reconciling cross-references.)*

10. ~~**Non-thread-safe RNG inside the OpenMP loops.**~~ **PARTIALLY FIXED, THEN REGRESSED 3 August -- see group B item 9.** `sample_beta_nocov_cpp_TS()` (`src/functions.cpp:308`) was switched to call the thread-safe `sample_beta_cpp_TS()` rather than the non-TS `sample_beta_cpp()`, which closed the race in `sample_betatheta_cpp_parallel()` (`src/functions.cpp:607`) -- the thread-safe path the audit described as "half-finished" was wired up on that side. **The `samplePGvariables()` path (`src/jsdm.cpp:398`) is still affected** via `randinvg()`'s `R::rnorm` call; that residual was tracked as a group B item and is now closed as *Fixed bugs* 28. *(Also not recorded at the time; recovered 27 July.)* **This entry no longer describes the code.** `43f2342` reverted `sample_beta_nocov_cpp_TS()` back to calling `sample_beta_cpp()`, reopening the exact race this entry records having closed. Verified empirically: two same-seed fits differ by up to 3.6 on `beta_theta_output` at 6 threads, bit-identical at 1. Do not read this entry as current status; read group B item 9.

**Post-audit fixes:**

11. ~~**Collection covariates are grouped by `Sample` alone.**~~ **FIXED.** Now groups by `c("Site", "Sample")` to match the site-then-sample ordering of `X_theta`. `R/runOccJSDM.R:682`.

12. ~~**The GP length-scale is never updated.**~~ **FIXED.** Enabled the `sample_ls()` call in `update_jSDMcoef()`. `R/jsdmfun.R:1063`.

13. ~~**`transformCovariatesMatrix()` never applies the stored factor levels.**~~ **FIXED.** Now correctly uses `df[[col]]` instead of `df[,col]` to re-level the column. `R/runOccJSDM.R:38`.

14. ~~**`listPriors$prior_beta_psi` / `prior_beta_psi_sd` have no effect.**~~ **FIXED.** Removed dead references and cleaned up the prior handling. `R/runOccJSDM.R:738-739`.

15. ~~**`summarisedLatentPresences = FALSE` errors on two-stage models.**~~ **FIXED.** Now correctly allocates and writes `w_output_chain` and `theta_output_chain`. `R/runOccJSDM.R:1152`.

16. ~~**`thinOutput()` cannot run at all.**~~ **PARTIALLY FIXED.** The function still exists and now runs: `niter` is read from `fitModel$results_output$jsdm_output$B0_output` instead of the long-gone `beta_ord_output`, and the `thin` argument is honoured rather than hard-coding `by = 5`. Two of the three original defects survived this fix, plus a third nobody had spotted. All were fixed by Claude on 4 August 2026 and are awaiting review as the `thinOutput()` item in group A; read that item, not this one, for current status.

17. ~~**`computeDiagnostics()` runs on `psi_output`.**~~ **FIXED.** Added `psi_output` to the skip list to avoid meaningless diagnostics. `R/diagnostics.R:403-404`.

18. ~~**`predictNewSites()` does not honour its documented no-op behaviour.**~~ **PARTIALLY FIXED -- this entry previously overstated what landed; corrected 27 July 2026.** The gating on the presence of covariates was done (`R/output.R:1457,1470`). The **`NULL` defaults were not**: `formals(predictNewSites)` shows `X_psi` and `X_s` still have no defaults at all, so `predictNewSites(fit, X_psi = X)` still fails on the missing `X_s` promise even when `useSpatial = FALSE`. **That residual is now closed: see *Fixed bugs* 34.**

    Third entry in this list found to overstate a fix (with 16 and 20), which is the argument for the test suite: every one of these was recorded from reading a diff rather than from running the code.

19. ~~**`set.seed(1)` inside `computeSpeciesDetected()`.**~~ **FIXED incidentally in b7b6aa2.** The RNG-resetting line was deleted. `R/output.R:2222`.

20. ~~**`computeSpeciesDetected()` crashes: `B` (bootstrap draw count) is referenced but never assigned.**~~ **PARTIALLY FIXED.** The crash is gone: `idxObs <- unique(floor(seq(1, min(500, niter), length.out = 500)))` and `B <- length(idxObs)` (`R/output.R:2221-2222`) both assign `B` and guard the fewer-than-500-draws case that the original hard-coded `B <- 200` did not. **The sampling is still a prefix, not a random draw** across chains, which was the substantive half of the original critique. **Now fully fixed** -- see Fixed bugs 23.

**Fix of 27 July 2026:**

21. ~~**`plotFPTPStage2Rates()` ignores `primerName`.**~~ **FIXED.** `p_output`/`q_output` are now subset to the requested primer *before* the quantiles are computed (`p_output[idx_primer, , , , drop = FALSE]`, `R/output.R:790-791`), instead of `apply(p_output, 2, ...)` pooling over primer, iteration and chain. `primerName` is matched against `fitModel$infos$primerNames` via `match(as.character(...))`, so both `primerName = 2` and `primerName = "2"` work (`primerNames` is stored as an *integer* vector, so the string form would otherwise silently fail); an unrecognised primer now errors with the list of available primers rather than quietly plotting the pooled result. The dead `idx_speciesprimer <- str_match(...)` line is deleted, and the plot title now reports which primer is shown.

    **Verified against `sampleresults`** (`P = 3`): the three primers now give three different plots (previously identical); each primer's intervals match hand-computed `quantile()` values on that primer's slice; no primer reproduces the old pooled result; the default still equals primer 1; `idx_species` subsetting and rendering are unaffected.

    **Note:** this was the *only* use of `stringr` anywhere in `R/`, so the package no longer needs it. Three `@import stringr` tags remain (`R/output.R:757`, `:857`, `:1008` -- the latter two were already spurious) along with `import(stringr)` in `NAMESPACE` and `stringr` in `DESCRIPTION`. Removing all three is a small, safe cleanup, best folded into the group C dead-code pass.

22. ~~**REGRESSION: `runOccJSDM()` crashes on every fit that has no spatial covariates.**~~ **FIXED** by Alex in `d6b70b1` ("Fixed bug on sample l_s"), the same day it was reported. Any call without `spatCovariates` had been dying in the first MCMC iteration with `missing value where TRUE/FALSE needed`, thrown from the acceptance test in `sample_ls()`, because `precomputeSORmatrices()` leaves `logDetKuu_grid`/`Lm1_grid` entirely `NA` when there is no spatial field but `update_jSDMcoef()` called `sample_ls()` anyway.

    **The actual fix was a one-word correction, and a better one than the fix suggested when this was filed.** The guard was `ps <- ncol(Bs)`; it is now `ps <- nrow(Bs)` (`R/jsdmfun.R`, in `update_jSDMcoef()`'s "read state variables" block). `Bs` is `[ps x S]` -- spatial basis dimensions by species -- so `ncol(Bs)` was the *species* count, which is never 0. `if (ps > 0)` was therefore always true regardless of whether a spatial field existed. `nrow(Bs)` is the spatial dimension and is correctly 0 with no spatial covariates, so `sample_ls()` is simply not reached. This addresses the root cause (a transposed dimension in the guard) rather than the workaround originally proposed here, which was to re-gate on `X_centers > 0`.

    **Verified after the fix** against the shipped `sampledata`: both `runOccJSDM(..., spatCovariates = c("Xs.1","Xs.2"))` and the same call *without* `spatCovariates` now complete. Tracing confirms the mechanism -- `precomputeSORmatrices()` still reports `X_centers = 0, grids all NA = TRUE` in the non-spatial case, but `sample_ls()` is no longer called, so the `NA`s are never read.

    **Still open from the original report:** the small-`n` `kmeans()` failure found while reproducing this (`number of cluster centres must lie between 1 and nrow(x)` at `n = 30` *with* spatial covariates) was not addressed by `d6b70b1` and has not been re-tested. It constrains the feasible grid for the simulation-study vignette, so confirm it before designing that.

23. ~~**`computeSpeciesDetected()` uses a prefix of the draws, not a sample across chains.**~~ **FIXED.** `R/output.R:2229` now reads `idxObs <- unique(floor(seq(1, niter, length.out = 500)))` -- the `min(500, niter)` cap is gone, so the 500 indices are spread evenly across the *entire* flattened draw set rather than landing in its first 500 rows. Since `beta_theta_output` is flattened iter-within-chain at `:2222`, spanning the full range means every chain now contributes. This resolves the substantive half of the original audit critique (see Fixed bugs 20, which fixed only the crash).

    Note the result is a *systematic* thin rather than the random draw from `seq_len(niter * nchain)` originally suggested. That is equivalent for this purpose and is standard practice for thinning an MCMC chain, so no further change is needed.

**Fixes of 28-29 July 2026 (Alex, `42198d9` / `e60e3ad`).** Each verified against the source rather than the commit message.

24. ~~**`sigma_h` is never sampled.**~~ **FIXED.** A new `sample_sigmah(U, a_sigmah, b_sigmah)` (`R/jsdmfun.R:1173`) is now called each iteration (`:1694`), so `sigmah_output` varies instead of sitting at its initialised 1. This was inert for fitted models -- `U` at the training sites is drawn under a hard-coded unit-variance prior regardless -- but it is the factor-score SD `predictNewSites()` uses when simulating `U` at *new* sites, so out-of-sample predictions no longer silently assume unit factor variance.

25. ~~**The prior mean for collection-covariate slopes is 1, not 0.**~~ **FIXED.** `b_betatheta <- rep(0, ncov_theta)` (`R/runOccJSDM.R:757`); previously `rep(1, ...)`, with only element 1 (the intercept) overwritten, so every actual covariate slope carried a prior centred on +1. The prior variance was also widened from `diag(1)` to `diag(2)`.

    This was the most clearly evidenced item in the whole audit: at R = 100 the intercept -- the one coefficient whose prior mean *was* overridden -- covered at 0.937, while the slopes covered at 0.496 with bias +0.113, against true slopes averaging +0.01. Across the full grid `beta_theta` undercovered in **every** cell (0.676-0.730).

    **Confirmed by the 29 July re-run**: `beta_theta` improved in every cell, to 0.709-0.771, with two-thirds of the bias removed. It is still well short of nominal, so a second cause remains -- tracked as the slope-overconfidence item in group B, and since narrowed to the collection-covariate slopes.

26. ~~**Non-thread-safe RNG in the hottest OpenMP loop.**~~ **FIXED.** `randinvg()` (`src/jsdm.cpp:86`) now draws from the `thread_local` `rnorm()` rather than `R::rnorm`; the old line is commented out beside it. That closes the last hole on the `samplePGvariables()` path, whose Polya-Gamma helpers were already converted in `53c38f1`.

    **Two consequences worth acting on.** A fixed seed should now reproduce on Linux and Windows, which (i) removes the reason `sampleresults.rda` had to be refit on macOS or an OpenMP-disabled build, and (ii) unblocks CRAN item 16. It also means the tier-1 constraint in `dev/simstudy/PLAN.md` §5.1 -- structural assertions only, never numeric equality -- can be revisited; it was adopted solely because of this race. Verify reproducibility on a multi-threaded platform before relying on any of that.

27. ~~**Stage 2 hyperparameters documented as settable but never read.**~~ **FIXED.** `a_p`, `b_p`, `a_q` and `b_q` now read from `listPriors` (`R/runOccJSDM.R:751-754`), matching the behaviour `@param listPriors` already claimed, and the documented defaults were corrected to match the code (5/1 and 1/20). A user running a low-detection study can now override the prior without editing the package.

    **Only the wiring is closed.** What the default *should be* remains open and is a design decision, not a defect. It was tracked as a group B item until Alex removed it in `093f2bb`; if that removal meant the decision is made, the chosen values should be recorded here, and if not the item needs restoring.

28. **`set.seed()` did not control any of the C++ samplers, so `runOccJSDM()` was not reproducible.** Found 29 July 2026 while writing the regression test for Fixed bugs 26; fixed the same day.

    Fixed bugs 26 replaced `randinvg()`'s use of R's global RNG inside an OpenMP loop with `thread_local` engines, which correctly closed the data race. But neither replacement engine ever read R's RNG state: `get_rng()` was seeded from the literal `12345 + omp_get_thread_num()`, and `mvrnormArmaQuick_TS()` from `std::random_device{}()`. Measured before the fix: two fits of the same fixture under the same `set.seed(4242)` differed by 5.09 on `B0_output`.

    Neither biased the posterior -- both are valid streams -- but users could not reproduce a fit, and the literal `12345` was per *process* rather than per run, so every simulation-study worker started its Polya-Gamma stream at the same position and replicates at the same position in different workers shared random numbers.

    **The fix** (`src/rng.h`, new). Both translation units now share one per-thread `mt19937`, seeded via `std::seed_seq` from a base seed that `runOccJSDM()` draws from R with `sample.int()` and installs through the new `setOccJSDMSeed()`. `set.seed()` therefore controls the whole fit, and consecutive fits stay independent because R's RNG advances between them. Three subtleties, each of which cost a debugging round:

    (a) A `thread_local` engine is constructed once per thread and would keep its original seed for the life of the session, so a generation counter triggers re-seeding when R installs a new base seed.
    (b) The generation must *not* be seed material. Mixing it in made the same base seed yield a different stream on the second fit of a session -- reproducible-looking in a fresh process, broken in a suite.
    (c) `std::normal_distribution` caches the second Box-Muller deviate, and that cache survives `rng.seed()`. Without an explicit `dist.reset()`, the first `rnorm()` of a fit could return a value left over from the previous fit, so reproducibility depended on the parity of normal draws earlier in the process. This is why the test passed standalone and failed inside the suite.

    `arma::randn()` elsewhere in the package was never affected: RcppArmadillo sets `ARMA_RNG_ALT` to route Armadillo's RNG to R's, so those draws already honoured `set.seed()`.

    **Scope.** Reproducibility holds for a given thread count. Threads derive separate streams, so if the package is ever built with OpenMP actually enabled (it is not on the macOS dev machine -- see group D), changing the thread count changes which stream produces which element. Inherent to per-thread streams, not a defect.

    **Tests.** `test-regression-bugs.R` covers same-seed equality, different-seed inequality, and independence of consecutive fits. Verified separately that the simulation harness gives identical output for a repeated replicate and different output across replicates.

    **Reviewed by Alex, 31 July 2026: "We don't care about reproducibility."** Taken as a decision not to invest further, not a request to revert, and the fix stays. Worth recording why: reproducibility here is **load-bearing internally even though it is not a user-facing priority**. The simulation study's paired design depends on it, and that pairing is what produced the strongest evidence in the whole study, that only 104 of 49,978 `resid_cor` coverage decisions flipped between the pre- and post-fix runs on identical truths. Remove the R-derived seeding and every future before/after comparison loses that power. The tier-1 test at `test-regression-bugs.R:243` guards it and should stay.

29. **The sparse-GP knot default no longer floors at 30 or crashes below 31 sites.** Filed 27 July 2026 (then group B item 3) after the simulation study hit it; fixed by Alex in `42198d9`, unlogged. Verified 29 July: `getDefaultSupportPoints()` (`R/jsdmfun.R:875`) is now `min(floor(n * 0.2), n - 1)`. The old `max(30, floor(n * 0.2))` fed `kmeans(X_s, centers = ps)` and so was a constant 30 for any dataset below 150 sites -- roughly one knot per site at n = 31, defeating the point of a sparse GP -- and errored outright below 31. The `n - 1` cap is what removes the crash.

30. **`ds = 0` no longer produces a null spatial field.** Filed 27 July 2026 (then group B item 4); fixed by Alex in `42198d9`, unlogged. The simulator's cross-species spatial covariance used to collapse to jitter at `ds = 0` -- measured `sd(spatField)` of 0.0019 against \~1.0 at `ds = 2` -- so any scenario built at `ds = 0` was silently a null-field test. Verified 29 July at seed 42: `sd(spatField)` is 0.598 at `ds = 0`, 0.467 at `ds = 1`, 0.678 at `ds = 2`. The study grid still uses `ds = 2`, now by choice rather than necessity.

    *Both of these were removed from group B by `42198d9` without a Fixed-bugs entry, which is why they are recorded here late. The check that caught it: `TODO.md`'s group B numbering had a gap at 5 and 6.*

31. **`plotCollectionRates()` errored on every input.** Reported by Doug 29 July 2026, fixed the same day. Failed with `object 'Min' not found` for any `fitModel`, with or without `idx_species`.

    `plotSpeciesRates()` (`R/output.R:819`) had been extracted as a shared helper and never wired up to its only caller. Three independent breakages in the same call path, which is why nothing had ever run it successfully:

    (a) the helper read columns `Min`/`Max`, while `plotCollectionRates()` passed the `2.5%`/`97.5%` that `quantile()` names;
    (b) it filtered on a `Species` column the caller never created -- the code that built one is still there, commented out, from before the extraction;
    (c) it referenced `speciesNames` as a free variable, present in neither its arguments nor the package namespace.

    Its roxygen documented a third state again: columns `Species`, `2.5%`, `97.5%`, and parameters `orderSpecies`/`subset` that the signature does not have.

    **Fixed** by giving the helper an explicit contract -- one row per species with `Min`/`Max`, plus `idx_species` and the *unsubset* `speciesNames` -- and doing the subsetting inside it, so the ordering is computed on exactly the rows plotted. Also added a clear error for models with no collection stage, which previously fell through to a subscript failure.

    **Note, since resolved:** this fixed the species-ordering defect for this function only. `plotOccupancyRates()`, `plotFPTPStage2Rates()` and `plotStage1FPRates()` still ordered on the full species set while indexing a filtered one (this was filed as group C at the time, not group B as originally written here). The test added here asserted label-to-value pairing rather than mere absence of error, and served as the template for fixing those three -- see Fixed bugs 32.

32. **`plotOccupancyRates()`, `plotFPTPStage2Rates()` and `plotStage1FPRates()` shared the species-ordering defect `plotCollectionRates()` had (Fixed bugs 31).** Filed as group C at the time; fixed by Claude 29 July 2026. Each computed `order()` on the *filtered* `idx_species` subset and then used the result to index the *unfiltered* `speciesNames`, so for any `idx_species` other than a prefix `1:k` the factor levels named the wrong species and bars silently vanished.

    `plotOccupancyRates()` and `plotStage1FPRates()` now delegate to the `plotSpeciesRates()` helper fixed in Fixed bugs 31, which subsets first and derives labels from the subset. `plotFPTPStage2Rates()` has a two-interval (`p`/`q`) layout that doesn't fit that helper, so it was fixed inline: order and labels are both now derived from the filtered `data_plot`. `plotStage2FPRates()` was already correct and untouched.

    Verified against a live fit with `idx_species = c(3, 1, 10)`: all four functions now plot exactly that subset with matching labels. Full test suite passes (119/119). `R/output.R`.

33. **`returnCovariateEffect()`/`plotCovariateEffect()` had no `idx_species` default, and fixing that exposed two further bugs in the code they call.** Filed as group C at the time; fixed by Claude 29 July 2026.

    Both functions declared `idx_species` with no default, so `returnCovariateEffect(fit, covName)` errored instead of defaulting to all species -- the same gap as `predictNewSites()` (Fixed bugs 34). Gave both a `NULL` default resolving to all species, matching every other return/plot function in `R/output.R`.

    That default immediately reached two pre-existing bugs that a narrow, prefix `idx_species` had been masking:

    (a) `plotCovariateEffect_base()` (`R/jsdmfun.R`) re-applied `apply(B_output, c(1,2), c)` to arrays its only caller, `plotCovariateEffect()`, had already collapsed. The second `apply()` collapsed the species margin again, leaving the array's third dimension sized by `ncov_psi` instead of `S` -- any species index beyond `ncov_psi` (2 in the test fit) errored with `subscript out of bounds`, exactly the region "all species" reaches immediately. Renamed the formals to `B0_output_vec`/ `B_output_vec` to match the already-collapsed inputs and dropped the redundant `apply()` calls.

    (b) `returnCovariateEffect_base()`'s per-species loop labelled each row `speciesNames[i]` (the loop counter) instead of `speciesNames[sp_idx]` -- the same class of defect as Fixed bugs 32 -- mislabelling every species whenever `idx_species` wasn't the prefix `1:k`.

    Verified against a live fit: the default now facets all 10 species correctly; `idx_species = c(3, 1)` plots and labels `OTU_3`/`OTU_1` correctly (previously mislabelled `OTU_1`/`OTU_2`, then errored once the default was exercised). Full test suite passes (119/119). `R/jsdmfun.R`, `R/output.R`.

34. ~~**`predictNewSites()` could not be called without supplying both `X_psi` and `X_s`.**~~ **FIXED 30 July 2026** (Claude; `R/output.R`, `src/jsdm.cpp`). Residual of the original audit's B.4, previously tracked in group C.

    **The filed defect.** `X_psi` and `X_s` had no defaults, so the `is.null()` guards the author had written could never fire: R raised `argument "X_psi" is missing, with no default` first.

    **Worse than filed.** Those guards used `&`, not `&&`. Since `&` evaluates both sides, the missing promises were forced even when the caller had asked for neither term, so `predictNewSites(fit, useEnvCov = FALSE, useSpatial = FALSE)` failed too. There was no way to call the function at all without both matrices.

    **The fix.** Both default to `NULL`. `useEnvCov` and `useSpatial` adopt the tri-state `useBiotic` already used: `NULL` uses the term if the fit estimated it, `TRUE` uses it and errors if it did not, `FALSE` skips it. That is what their roxygen always claimed, and was unreachable because `useSpatial` hard-stopped on a fit with no spatial field instead of ignoring it. Strictly more permissive: every call that worked before still works identically.

    **Three further defects, exposed by making those paths reachable, all pre-existing and all fixed here:**

    (a) `computeNewOutputs()` sliced `Ks_all` and `Bs_output` unconditionally, so `useSpatial = FALSE` aborted with `Cube::slice(): index out of bounds`. Moved inside the `useSpatial` guard.
    (b) It read the new-site count as `X.n_rows` unconditionally, so `useEnvCov = FALSE` silently returned a zero-row result. Now taken from whichever term is active; both being off is an explicit error, since nothing then determines how many sites to predict for.
    (c) A fit with no spatial field returns `Bs_output` with a zero-length first dimension, which collapses under `apply()` and broke the `aperm()` with `'perm' is of wrong length`. Only reshaped when the spatial term is in play.

    **Tested beyond absence of error** (`test-api-contracts.R`): each term measurably changes the prediction when toggled, probabilities stay in `[0,1]`, and quantiles stay ordered. Without the first of those the switches could have been cosmetic and the shape assertions would still have passed.

    **Noticed, not fixed:** `computeNewOutputs()` prints `Computing species i out of S` to stdout via `Rcout` on every call, unconditionally, and it cannot be silenced. Filed separately in group C.

35. **The vignette could not be built, so `R CMD check` never reached code inspection.** **FIXED 30 July 2026** by Doug regenerating `data/sampleresults.rda`.

    Two failures in sequence, each hidden behind the previous one. First `plotCollectionRates()` errored on every input (Fixed bugs 31). With that fixed, the build failed in `plotCovariateEffect()`, apparently for want of a `covNames` default; naming the covariate in the chunk did not fix it either, and it then failed with `'from' must be a finite number`.

    **The real cause was neither function.** Both read `fitModel$infos$X0_psi`, the raw occupancy covariates, to build a prediction grid. `X0_psi` was added by `e60e3ad` ("Added GAMs"), and the shipped `sampleresults.rda` had last been regenerated on 26 July, before that. `min(NULL, na.rm = TRUE)` is `Inf`, and `seq(Inf, -Inf, ...)` throws. `X0_psi` was the only name a current fit's `infos` carried that the shipped object lacked.

    **Confirmed after the refit:** `X0_psi` is present (100 x 2), `list_X_psi_mat` carries `bs_info`/`target_spline_vars`, both vignette chunks run, and `devtools::check()` completes with **0 errors** for the first time. Remaining status: 2 warnings, 3 notes, tracked separately.

    **The lesson worth keeping**: a stale shipped dataset presents as a bug in whatever function touches it first, and moves to the next function each time one is fixed. Two functions were investigated and one was needlessly suspected before the data was. When an example object fails and a fresh fit does not, suspect the object.

36. ~~**Both exported GAM functions failed for any user with a categorical occupancy covariate.**~~ **FIXED 30 July 2026** (Claude; `R/occJSDM-package.R`, `NAMESPACE`, commit `032dcff`). Logged 1 August: it was completed but never given a *Fixed bugs* entry, so it was invisible in this file.

    **The defect.** `tidyr::pivot_longer()` was called but never imported, and `tidyr` was in `DESCRIPTION` `Imports:` with nothing imported from it in `NAMESPACE`. It is reached in the categorical-covariate branch of `returnCovariateEffect()` and `plotCovariateEffect()`, so both would fail with "could not find function" from a properly installed namespace. `stats::setNames` and `stats::rnbinom` were missing the same way.

    **Why it survived the test suite.** `devtools::load_all()` resolves unimported symbols through the global environment, so a functional test passes whether or not the import exists. The regression test therefore asserts on the imports environment directly, which holds under both `load_all()` and `R CMD check`.

    **Not cosmetic, despite arriving as an `R CMD check` NOTE.** The "checking dependencies in R code" NOTE is now gone entirely, 3 notes to 2. `dnbinom` and `bs` remain on the undefined-globals list, correctly: they are reached only by dead code, and importing `bs` would add `splines` to `DESCRIPTION` to support code that should not ship.

    **Reviewed by Alex, 2 August 2026: "Agree with the fix."**

37. ~~**Ten dead functions moved to `deprecated/`.**~~ **DONE 30 July 2026** (Claude; commits `7bae018`, `f2e2701`). Logged 1 August so the deprecation has a stable anchor; it was previously cited only by section-A position, which has since been reused.

    `sample_z`, `sample_w`, `sample_cimk`, `sample_betatheta` from `R/jsdmfun.R`, plus six dead samplers and plots from `R/mcmcfun.R` and `R/output.R`. All were unreachable from any exported entry point. Roughly 26 more remain, tracked in group D, deliberately deferred until the sampler rewrite lands so the two do not collide.

    **Reviewed by Alex, 2 August 2026: "Ok to deprecate these functions."**

38. ~~**`main` would not load at all, and two fixes were needed to restore it.**~~ **FIXED 31 July 2026** (Claude; `src/Makevars`, `src/Makevars.win`, `DESCRIPTION`, `NAMESPACE`, `R/occJSDM-package.R`, `R/runOccJSDM.R`). Closed 2 August 2026. After `8f9f315` the test suite went from 167 passing to 25 passing, 3 failures and 33 errors, because the built `.so` had an undefined `RcppParallel::tbbParallelFor` and failed to load, taking every function in the package with it.

    **Change 1, the linking.** `RcppParallel` was in `LinkingTo`, which makes the headers visible but does not link the libraries. Added the documented `RcppParallel::RcppParallelLibs()` call to both `Makevars` files, **and** added `RcppParallel` to `Imports` with an `importFrom`. **Both halves are needed**: the Makevars line alone fixes the missing symbol but the load still fails, because on macOS `libtbb` is referenced through `@rpath` and `devtools` copies the `.so` to a temp directory. Importing the package makes its namespace load first and set up the TBB paths. The reasoning is in a comment in `src/Makevars` so the next person does not stop at the first fix and conclude it did not work.

    **Change 2, the debugging scaffold at `R/runOccJSDM.R:404`.** It had live assignments overwriting every argument with values referencing `occ_data_effort`, a dataset not in the package, plus a live `summarisedLatentPresences`. That block was fully commented before `8f9f315`. **Re-commented rather than deleted**, since Alex evidently uses it, with a note saying it must stay commented and what happened when it did not.

    **Two things left alone at the time, both since resolved by Alex:**

    - `verbose` in `computeNewOutputs()` defaults to `T`, so `predictNewSites()` still prints one line per species unless a caller opts out. It is suppressible now, which was the harder half, but the default behaviour, the test-output noise and the CRAN-reviewer exposure are all unchanged. **Alex: "We can leave verbose on, people won't think of turning it on."**
    - `sample_z_cpp_parallel()` was exported but not yet called at the time; `runOccJSDM` still used `sample_z_cpp`, so only `sample_w_cim_cipp_parallel()` was wired in and half the parallelisation work was unreachable. **Alex: "sample_z_cpp_parallel now used."** Confirmed: it is called at `R/runOccJSDM.R:1113` as of `41abe69`, which leaves `sample_z_cpp()` itself with no callers -- whether *that* should now be deprecated is a new open item in group A.

    The `verbose`-default choice noted above is the same decision closed separately as *Fixed bugs* 39.

39. ~~**`computeNewOutputs()` prints to stdout on every call and cannot be silenced.**~~ **FIXED** (Alex added a `verbose` argument to `predictNewSites()`, threaded through to the `Rcout` call in `src/jsdm.cpp`). Closed 2 August 2026.

    **The original defect.** `src/jsdm.cpp` ran `Rcout << "Computing species ..."` inside the species loop unconditionally, so every `predictNewSites()` call printed one line per species, and `suppressMessages()` did not catch it because `Rcout` is stdout rather than R's condition system.

    **Verified live, not just from the response.** `verbose` reaches the C++ `if(verbose)` gate around the `Rcout` call (`src/jsdm.cpp:482`) via `computeNewOutputs()` (`R/output.R:1683`). Ran both ways on a fitted model: `verbose = FALSE` suppresses all four per-species lines; `suppressMessages()` alone, with `verbose` left at its default, still lets all four through -- expected, not a gap, given the default chosen below.

    **The item's fix spec asked for a default of `FALSE`; Alex chose `T` instead, deliberately.** `predictNewSites()` still prints by default unless a caller opts out, which is not what was originally asked for -- but it is a design decision, not an unfinished fix. Same decision already recorded under *Fixed bugs* 38 (the RcppParallel-linking fix): **"We can leave verbose on, people won't think of turning it on."** **Closed by decision, not by elimination.**

40. ~~**`sample_z_cpp()` was exported but had no callers.**~~ **DE-EXPORTED 2 August 2026** (Claude; `src/functions.cpp`, `R/RcppExports.R`, `src/RcppExports.cpp`).

    `41abe69` wired `sample_z_cpp_parallel()` into `runOccJSDM()`, leaving the serial version reachable only as an unused wrapper. Tag removed and `Rcpp::compileAttributes()` re-run; the C++ body stays in `src/` as the serial reference implementation of the parallel one. Detail and decision provenance in `AGENTS.md`, "Detail behind Fixed bugs 40-42".

41. ~~**`BBSL_Worker` called the non-thread-safe `sampleB_SoR()`, racing on R's RNG from every TBB thread.**~~ **FIXED 2 August 2026** (Claude; `src/jsdm.cpp`). Introduced by `41abe69`; found the same day while running the `continuous` arm for the `B0` item in group B.

    `sampleB_SoR()`'s `arma::randn()` reached R's global RNG from every TBB worker thread, permitting duplicate and torn draws, so the chain could fail to target the intended posterior. It now draws from `rnorm()` in `src/rng.h` instead. The race is closed and results are statistically valid again.

    **Bit-reproducibility above one thread is not restored**, because TBB work-stealing varies which thread draws for which species even at a fixed thread count. Set `RCPP_PARALLEL_NUM_THREADS=1` whenever a run has to be reproducible; the suite is 248 passing single-threaded, and two tests in `test-regression-bugs.R` are pinned to one thread for this reason, with a companion skip naming the gap. Closing it is *MEE paper* Alex to-do 8. Measurements and full mechanism in `AGENTS.md`, "Detail behind Fixed bugs 40-42".

42. ~~**`sampleB_SoR_TS()` was exported, had no callers, and its name invited the exact error it was written to prevent.**~~ **DEPRECATED 2 August 2026** (Claude; `src/jsdm.cpp`, `deprecated/jsdm-sampleB_SoR_TS.cpp`, `R/RcppExports.R`, `src/RcppExports.cpp`). Closes the decision left open at the end of *Fixed bugs* 41.

    Written as the thread-safe variant of `sampleB_SoR()` for the race above, and left purposeless when that race was closed another way. Worse than merely dead: it seeds from OS entropy, so anything built on it would ignore `set.seed()` entirely, despite the `_TS` name reading as the endorsed choice. De-exported and moved to `deprecated/`, not deleted. Group D's dead-wrapper count drops to 11 of 35, re-measured rather than decremented; the suite is unchanged at 248 passing. Detail in `AGENTS.md`, "Detail behind Fixed bugs 40-42".

43. ~~**The `runOccJSDM()` debugging scaffold went live for the fourth time.**~~ **FIXED 2 August 2026** (Claude; `R/runOccJSDM.R`, `tests/testthat/test-regression-bugs.R`). First occurrence `8f9f315` (*Fixed bugs* 38); second `46d8804`, fixed in `e5d0105`; third some time before `11981a1`, which Alex fixed himself; fourth `522b89e`, minutes later, apparently as a side effect of editing nearby lines rather than a deliberate change.

    Re-commented as before. This time also added a static regression test, since four occurrences of the identical defect is not a coincidence to keep fixing by hand: it reads `R/runOccJSDM.R` and asserts the live `data = occ_data_effort` line is not present, catching the defect before any fit is attempted rather than after. Verified: tier 1 passes clean, 144 of 144.

44. ~~**Two of the "assorted smaller items" in group C.**~~ **FIXED** (Alex; `R/runOccJSDM.R`, `R/output.R`). Closed 2 August 2026, Alex marked both in `f4a59e4`.

    (a) `createDataIdx()`'s `maxP` for `model = "occupancy"`. Verified moot rather than fixed as originally described: `maxP` is used nowhere outside `if (model == "two_stage")` blocks, so the described crash path is not reachable as the code now stands. Confirmed live via the tier-1 occupancy smoke test, which passes.

    (b) `computeSpeciesDetected()`'s roxygen no longer documents any signature at all -- title, description and `@noRd` only -- so the stale Beta-approximation `@param` block is simply gone.

45. ~~**`listPriors$b_betatheta_slope_var`, a new prior hook.**~~ **APPROVED 2 August 2026** (Alex: "Happy with the new fix"). Closes group A item 1.

    Exposes `B_betatheta`'s previously hard-coded slope variance, default 2 unchanged. Built as a diagnostic for the `beta_theta` slope item, which is still open; approving the hook is not a decision on the value, which remains the `b_betatheta` variance item in group B. Alex's `522b89e` also removed the explanatory comment at the call site that had marked it as provisional, consistent with treating it as a permanent feature now.

46. ~~**`B0`'s bias doubled between the pre- and post-fix runs.**~~ **CLOSED BY DECISION 2 August 2026** (Alex: "we can close the point on B0 and assume that the bias is only due to the confounding with beta_theta"). Not closed by elimination -- the quantitative link was never established, only that `B0` is unbiased in both arms tested with `beta_theta` absent entirely (`binary` +0.0122 at 1.0 SE, `continuous` +0.0066 at 1.2 SE, `PLAN.md` 16.4), against -0.0633 at 3.3 SE with only the intercept present. Alex accepted that as sufficient. Full investigation, including the two hypotheses tested and ruled out along the way, in `AGENTS.md`.

    The coverage sub-finding this item also carried (`continuous` at 0.879 against nominal, `PLAN.md` 16.5) was not addressed by this decision and is carried forward as its own item in group B.

47. ~~**Three assorted smaller items (C2b, C2c, C2d), all verified fixed (Alex).**~~

48. ~~**Item 1 in group B (`sample_ls()` "scores the wrong density").**~~ **CLOSED BY EVIDENCE 6 September 2026** (Alex + Doug; `R/jsdmfun.R`). Three rounds of evidence converged on `sample_ls()` not being the defect. (1) Four candidate causes inside `sample_ls()` itself were ruled out in July (missing amplitude, wrong amplitude passed, `logDetKuu` factor, weak data; full detail in `AGENTS.md` under "*Evidence behind the group B items*" `/ sample_ls(): four candidate causes ruled out`). (2) Alex's reply on the joint TODO exchange: "I run a simple test in the file `R/test_sample_l.R`. It always go to the true value (or very close by)", and later "The code has been fixed. The code only checks that sample_ls works, which it does in the sense that it converges to the true value. If there is a problem, it is somewhere else." (3) Doug's hand-built-equivalent run: `l_s_grid`, `Lm1_grid`, `logDetKuu_grid` derived directly from `K2()` + `chol()`; `SE = t(LU) %*% Z` at `l_s_true = 0.12` (truth at grid idx 6). Starting `idx_ls` at 1, 20, and 6, all three chains settle on idx 6 within \~100 iterations and stay there (mode share \~1.0 post-burn). What the test does not establish is what `SE` `update_jSDMcoef()` actually constructs once `Bs` has been redrawn -- that and the `precomputeSORmatrices()`-vs-current-`Bs` alignment are now item 10 in group B (rails-at-top upstream of `sample_ls()`). Items fully removed from group B is 1; the finding lives on as item 10. Note that `simulateData()` errors with `invalid 'sizeargument` on the *current* `dev/test_sample_l.R` (moved there from `R/` on 8 September 2026, where its top-level code broke package installation), so Alex's run depends on a session-internal snapshot of the helpers and is not reproducible against `main` head without manual setup. **(b) `d > NULL` errors on single-species input.** `get_param(listParams, "n_factors")` returns 0 by default; the old cap `if (d > ncol(OTU))` errored when `OTU` was a single-species vector because `ncol()` returns `NULL`. Fixed by branching first on `ncol(OTU) > 1`; the `else` branch sets `d <- 0` with a message. `R/runOccJSDM.R:759-770`.

    **(c) `reparamFactorModel(A_output, C_output)` fails when `gt > ncov_psi`.** `C_output` is `[gt x ncov_psi]`; when `gt > ncov_psi`, `qr.Q()` returns a non-square `[gt x ncov_psi]` matrix and `Q %*% diag(diag(R), nrow = gt)` fails on conformability. Fixed by `gt_default <- floor(sqrt(min(S, ncov_psi)))`, which ensures the default `gt <= ncov_psi`. A caller who manually sets `listParams$n_lattrait > ncov_psi` could still trigger it, but that is a usage error. `R/runOccJSDM.R:773`.

    **(d) Spatial-covariate `is.numeric` check ran before the names-present check.** A mistyped `spatCovariates` name gave `undefined columns selected` rather than the intended "Covariate names provided not in data\$info". Fixed by reordering: the `%in% colnames(data$info)` check (`:570`) now precedes `sapply(data_info[,spatCovariates], is.numeric)` (`:574`). `R/runOccJSDM.R:570-576`.

# **Completed work**

Finished work, kept for context rather than as tasks. Bug fixes live under *Fixed bugs* above; this is everything else.

1.  ~~**Purge `data/traitdata_caiwang.rdata` from git history.**~~ **DONE.** Purged via `git filter-repo` and force-pushed the rewritten `main` to `origin`. All commit hashes from `11c5449` onward changed as a result, so **Alex needed to re-clone or hard-reset**; that has since happened.

2.  ~~**Port the model-diagnostics functions from the GLGS-eDNA repo.**~~ **DONE.** Alex's note: *"I ALREADY ADDED DIAGNOSTICS SO NOT NEEDED ANYMORE"*. `computeRhat()`, `summarisePosterior()`, `plotTraceplot()` and `returnConvergenceDiagnostics()` are in `R/diagnostics.R`; the first two were later de-exported (`@noRd`) as redundant with `returnConvergenceDiagnostics()` for user-facing purposes. A "Model diagnostics" vignette section was added to `vignettes/occJSDM.Rmd`.

    **Not ported, and still available if wanted:** `compare_to_true()`/`plot_estimated_vs_true()`, which need `true_params` from a simulation. The simulation suite (item 3 below) now covers that ground more systematically, so these are probably redundant rather than outstanding.

3.  ~~**Build the simulation test suite and run the coverage study.**~~ **DONE 28 July 2026**, except the three follow-ups still listed under *MEE paper / Doug to dos 2*. Full specification and results: `dev/simstudy/PLAN.md`. Provenance and lessons: `AGENTS.md`.

    - **Tier 1** (`6d9526d`): 89 structural assertions, \~7 s, **ships to CRAN**. `helper-fixtures.R`, `test-smoke-configs.R`, `test-regression-bugs.R`, `test-api-contracts.R`. `test-placeholder.R` deleted.

    - **Tier 2** (`f63eeeb`): recovery canary, \~30 s, `skip_on_cran()`, thresholds measured over eight seed sets.

    - **Tier 3 + runner** (`f63eeeb`, `6123036`, `e0718f1`, `32b56a2`): env-gated study plus `dev/simstudy/run_study.R`, with checkpointing, `--resume`, `--caffeinate` and progress reporting.

    - **The run itself:** re-run 29 July 2026 after Alex's fixes and the RNG seeding fix. 10 scenarios x 100 replicates, 1000 fits, 0 failures, 285 min on 5 cores, 155,578 interval checks. Results in `dev/simstudy/results/`, tabulated in `PLAN.md` §12. The 27-28 July pre-fix run is retained for comparison (`PLAN.md` §12.6).

    - **It is a paired comparison, and that is worth protecting.** `draw_truth()` seeds on (scenario, replicate), so the simulated data and true values are bit-identical between the pre- and post-fix runs -- verified, `max|truth difference| = 0`. Every difference is therefore attributable to the code rather than to sampling variation, which is what makes statements like "only 104 of 49,978 `resid_cor` decisions flipped" possible. **Do not change `simstudy_seed()`**, or future runs lose comparability with these two.

    It found six defects -- three while the tests were being written, three from the run -- none of which the static audit had caught. Four of the six are now fixed (Fixed bugs 24-27 and group B); the rest are the rails-at-top (`l_s` ranges, item 10), `reparamFactorModel()` and `beta_theta` slope items in group B.

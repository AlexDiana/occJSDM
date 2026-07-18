# occJSDM

R package for fitting occupancy models to eDNA metabarcoding reads data (joint species distribution model / JSDM framework), accounting for a two-stage detection process (field collection + PCR amplification/lab detection), species traits, and optional spatial autocorrelation.

## Package layout

- `R/runOccJSDM.R` -- main model-fitting entry point (`runOccJSDM()`), plus data-prep helpers (`process_covariates()`, `createDataIdx()`, `get_param()`).
- `R/jsdmfun.R` -- core JSDM machinery, including `simulateData()` (the lower-level data simulator) and coefficient/variance-partitioning helpers (`computeBtcoef()`, `computePsiCoef()`, `computeVariancePartitioning()`).
- `R/simulateData.R` -- exports `simulateOccJSDMData()` (the higher-level wrapper that simulates a full occJSDM-style dataset -- occupancy + two-stage detection + read counts -- from `simulateData()`) and `toRunOccJSDMFormat()` (converts `simulateOccJSDMData()` output into the `info`/`OTU` list format `runOccJSDM()` expects, i.e. the same shape as `sampledata`).
- `R/mcmcfun.R`, `R/output.R`, `R/diagnostics.R` -- MCMC sampling and post-processing/output/plotting functions used after `runOccJSDM()` (e.g. `returnOccupancyCovariates()`, `plotVariancePartitioning()`, `returnLatentPresences()`, `plotResidualCorrelationMatrix()`; see NAMESPACE for the full exported list).
- `src/jsdm.cpp`, `src/functions.cpp` (via `RcppExports.R`) -- Rcpp/Armadillo backend, including the spatial kernel `K2()` (squared-exponential GP kernel over site coordinates) and `sample_w_cpp()` (samples latent collection state `w` from continuous read intensities when `threshold = 0`).
- `vignettes/occJSDM.Rmd` -- walkthrough of fitting a model with `runOccJSDM()` and using the output/plotting functions.
- `vignettes/simulateOccJSDMData.Rmd` -- walkthrough of `simulateOccJSDMData()`: parameter lists, simulating and converting to `runOccJSDM()` format via `toRunOccJSDMFormat()`, simulating with/without spatial autocorrelation via `useSpatField`, checking variance partitioning, visualizing occupancy and the linear predictor spatially, and how trait covariates (`Tr`) shape species-specific occupancy coefficients via the `G` matrix.
- `TODO.Rmd` -- structured outstanding feature list, organized as **v0.1.0-beta Public release** (Alex to dos / Doug to dos), **MEE paper** (Doug to dos / Alex to dos), and **Future versions**. See "Current work status" below for the current item list.
- `analysis/analysis.R` -- ad hoc analysis script (not part of package build).
- `tests/testthat/` -- unit tests (testthat edition 3, set up via `usethis::use_testthat(3)`); currently covers `toRunOccJSDMFormat()` (dimensions, column names, intercept handling, mismatched-dimension errors).
- `data/` -- `sampledata.rda` / `sampledata_orig.rda` (example dataset used in `occJSDM.Rmd`), `sampleresults.rda` (precomputed fit used by the vignette, `nchain=2, nburn=5000, niter=5000, nthin=1`), `data_out_env_trait.rdata`, and `CaiWang_data/` (a data directory the user is actively reorganizing; `traitdata_caiwang.rdata` was recently removed in favor of this).
- `CITATION.cff` -- citation metadata; primary citation is the Ji et al. (2025) *Ecology Letters* methods paper, secondary is the software repo (`https://github.com/AlexDiana/occJSDM`). No CRAN/Zenodo listing exists, so this positioning follows standard practice for a research tool wrapping a published method.
- `README.md` -- minimal readme with install instructions, vignette pointers, and a "How to cite" section (since GitHub's citation button is easy to miss).
- `DESCRIPTION` -- `Authors@R` lists Alex Diana (`Alex.diana92@yahoo.it`, maintainer/`"cre"`) and Douglas W. Yu (`dougwyu@mac.com`, `"aut"`).

## Data structures

**Model input (`runOccJSDM()`)**: a list with `info` (data.frame, one row per PCR replicate, with `Site`, `Sample`, `Primer`, occupancy covariates `X_psi.*`, and collection covariates `X_theta.*`) and `OTU` (matrix of read counts, rows matching `info`, columns = species).

The `threshold` argument to `runOccJSDM()` controls how `OTU` is interpreted: - `threshold > 0` (default `1`): reads are truncated to binary presence/absence (`OTU >= threshold` -\> detection), and `w` is sampled via `sample_w_cim_cipp()`. - `threshold == 0`: reads are modeled continuously via `logy1 = log(OTU + 1)` and a two-component Normal mixture -- given a true detection, `logy1 ~ Normal(mu1, sigma1)`; given a false-positive/contamination detection, `logy1 ~ Normal(mu0, sigma0)` -- sampled via `sample_w_cpp()`. `mu1_output`/`sigma1_output`/ `mu0_output`/`sigma0_output` in `results_output` are only populated in this mode.

**Simulated data (`simulateOccJSDMData()`)** returns `list(true_params, data_list)`:

- `data_list`: `X_psi` (occupancy covariates), `X_theta` (collection covariates), `Xs` (site coordinates, `n` x 2 -- always simulated, only affects occupancy if `useSpatField = TRUE`), `Tr` (species traits, `S` x `g` -- always affects occupancy via the trait-response matrix `G`), `y` (simulated **read counts**, matching the `threshold = 0` mode of `runOccJSDM()` -- not binary detections; see below).
- `true_params`: true occupancy/detection states and coefficients, including `jsdmParams_true$varPart` (per-species Environmental/Spatial/Biotic variance partitioning), `jsdmParams_true$eta` (linear predictor, useful for inspecting the spatial field without the noise added by thresholding to binary occupancy), and `mu1_true`/`sigma1_true`/`mu0_true`/`sigma0_true` (the true read-intensity parameters used to generate `y`).

`simulateOccJSDMData()` takes `list_datasettings` (`n`, `S`, `g`, `M`, `P`, `K`, `ncov_psi`, `ncov_theta`), `list_params` (`p`, `q`, `theta0`, `theta_baseline`, and optionally `mu1`, `sigma1`, `mu0`, `sigma0` -- defaulting to `5`, `1`, `1.5`, `1` -- controlling simulated read intensity), `list_jsdmParams` (`gt`, `d`, `ds`, `sigma_b`, `sigma_bs`, `sigma_ts`, `sigma_h`, `l_s`), and `useSpatField` (default `FALSE`). Read counts are generated as `y = round(exp(logy1) - 1)` where `logy1` is drawn per-replicate from `Normal(mu1, sigma1)` (true detections), `Normal(mu0, sigma0)` (false-positive/contamination), or `0` (no reads). With the defaults, this produces a distribution roughly similar in shape to `sampledata$OTU` (\~two-thirds zeros; nonzero median around 120), though the upper tail runs somewhat heavier.

**`toRunOccJSDMFormat(sim, n, M, P, K, drop_theta_intercept = TRUE)`** expands `sim$data_list$X_psi`/`X_theta` (site/sample-level) to the PCR-replicate level (via the same indexing as `createDataIdx()`), builds `Site`/`Sample`/`Primer` id columns, and returns `list(info, OTU = sim$data_list$y)`. `n`/`M`/`P`/`K` must match the `list_datasettings` originally used to create `sim` -- they are not stored in `sim` itself, so passing inconsistent values will error or silently mis-index. `X_theta`'s intercept column (always column 1) is dropped by default since `runOccJSDM()` adds its own via `process_covariates()`.

## Known issues in existing (pre-`toRunOccJSDMFormat`) code

- `predictOccupancyProbs()` (`R/output.R`) is currently broken: it references an undefined `X_ord` object instead of its own `X_s` argument, and errors even earlier because `nchain`/`niter` are derived from `beta_psi_output`'s dimensions, which don't have the expected shape on at least some fitted models. Confirmed broken via live testing; not fixed as of this session.
- `returnLatentPresences()` (`R/output.R`) references an undefined `varPart_output` object (likely leftover from a copy-paste) and will error at runtime.
- Both bugs already have `@note` roxygen comments flagging them in the source.
- `computeAverageCollectionProbs()` and `computeConditionalSamplePresenceProbs()` were confirmed working (via live testing against fitted model objects) as long as the fitted model actually populated `results_output$theta_output` / `results_output$w_output` (always true for a fresh `runOccJSDM()` fit with default `summarisedLatentPresences = TRUE`).
- `computeMinESS()` (`R/diagnostics.R`, bug introduced in commit `af5ebe3`): an `ESS_beta0psi` matrix is allocated to hold ESS for the new `beta0_psi_output` (intercept) term, but the loop body still populates the pre-existing `ESS_betapsi` from `beta_psi_output`, leaving `ESS_beta0psi` unfilled (all NA) and excluded from the final `min()`. Flagged as TODO.Rmd item 1.5 (Alex to dos); not yet fixed.
- Trait-reading fragility in `runOccJSDM()`: it checks `data$traits` via partial name-matching against `data$traitsMatrix` (relying on `$`'s partial matching, since `traitsMatrix` is the only element of `data` starting with `"traits"`), not the unused `traitsMatrix` function argument. Flagged as TODO.Rmd item 1.4 (Alex to dos).
- Scalar `p`/`q` indexing bug in `simulateOccJSDMData()`: `list_params$p`/`list_params$q` are scalars but the function body indexes them as matrices (`p_true[idx_p_k[i], s]`); needs `matrix(p, nrow = P, ncol = S)` expansion before use. Flagged as TODO.Rmd item 1.3 (Alex to dos).
- Counts model (`data_type == "counts"`, auto-detected from `OTU` values) is unsupported downstream: `stop("Counts model not supported yet")`. No explicit user-facing `count=` argument exists (see TODO.Rmd item under MEE paper / Alex to dos, "ability to analyse count data").

## Model inference logic (`inferDataModel()`, `R/runOccJSDM.R:110-155`)

`runOccJSDM()` classifies the fitted model type based on **row-level duplication of `Site`/`Sample` in `data$info`**, not directly on M/K/P:
- No repeated values in `Site` → JSDM-only (M=K=P=1).
- `Site` repeats, `Sample` does not → classical occupancy model (M>1, K=P=1).
- `Sample` repeats (due to K>1 and/or P>1) → two-stage eDNA model, including the case M=1 with K>1 or P>1 (a single site sampled with multiple PCR replicates/primers).

## Output functions of note

- `returnVariancePartitioning(fitModel)` -- per-species table of (Env, Spatial, Biotic, StDev) variance fractions; feeds `plotVariancePartitioning()`.
- `returnResidualCorrelationMatrix(fitModel, confidence = .95)` -- returns a 3 × S × S array (quantile × species × species) of posterior credible intervals for pairwise residual correlations. The 50% slice gives the median correlation matrix; comparing bounds' signs gives a significance flag.
- `plotCumulativeSpeciesDetections(fitModel, K, primer = 0, alpha = .95)` -- credible interval for cumulative species detected as a function of PCR replicates (K), via bootstrapped per-species Beta-distributed detection probabilities. No M-based (sample/visit-level) equivalent exists yet (TODO.Rmd, MEE paper / Alex to dos).
- Ordination is incomplete: `returnOrdination()` and `plotOrdinationScores()` (both unexported, near-duplicates) return a `ggplot` of site factor scores with credible-interval error bars, despite `returnOrdination()`'s docstring claiming to return a plain 3 × sites × factors quantile array. Neither handles species loading scores. TODO.Rmd item 1.1 (Alex to dos) calls for both site and species scores as plain tables/matrices plus ggplot2-based plotting, still unaddressed.

## Git and build artifacts

- `src/*.o` and `src/*.so` files are **tracked in git** despite being ignored in `.gitignore` (added in lines 7-8). These are compiled object files and should not be committed. To fix: `git rm --cached src/*.o src/*.so` and commit. They were tracked before `.gitignore` rules were added and remain in git history until explicitly removed.

## ggtern and plotting

- `plotVariancePartitioning()` (via `plotVarPart()` in `R/jsdmfun.R`, line 1466) produces a warning: "Ignoring unknown labels: L, T, R" because it calls `labs(L = "Environment", T = "Biotic", R = "Spatial")`. The `ggplot2::labs()` function doesn't recognize ternary-axis parameters `L`, `T`, `R`. **Fix**: Remove the `labs()` call entirely — the axis labels are correctly set by the aesthetic names (`Env`, `Biotic`, `Spatial`) in the `aes()` mapping and don't require additional specification. The theme elements `tern.axis.title.T`, `tern.axis.title.L`, `tern.axis.title.R` already control styling.

## `man/*.Rd` tracking

`man/*.Rd` files are currently tracked and committed in git (not gitignored), despite an earlier commit (`dd158a9`, prior to this session) whose message claimed intent to "stop generating/tracking man/\*.Rd" -- that intent was never followed through in `.gitignore`. Regenerating docs via `devtools::document()` after adding/changing roxygen tags will produce `man/*.Rd` diffs that should be committed (or the gitignore situation resolved) deliberately.

## Current work status

- **Most recent commits** (as of July 2026, newest first): `03a0dd7`, `da32649`, `f555bb0`, `ae317bf` (all "Update/clean up TODO.Rmd"), `e957292` (comprehensive `occJSDM.Rmd` vignette overhaul -- filled placeholders, rewrote the M/K/P model-inference bullet list to describe the Site/Sample row-duplication mechanism accurately, added sections on cumulative species detections and residual correlation structure, renamed several section headings, added references to Cai et al. 2025, Leibold et al. 2021, Pichler et al. 2025), `4099ae8` (Update TODO.Rmd), `3d3cd6f` (marked items done, restructured MEE paper section into Doug/Alex sub-lists, flagged the `computeMinESS()` bug), `276e3ec` (added `CITATION.cff`/`README.md`, fixed `DESCRIPTION`'s `Authors@R`), `24df008` (merge from Alex's fork), `af5ebe3` (Alex's pulled commit -- fixed several real bugs: mislabeled "continuous" data-type message, `model == "occupancy"` assigning the fallback ID to `Site` instead of `Sample`, moved OTU-reading/NA-handling earlier, threaded `list_jsdmParams$tau` through to `simulateData()` instead of hardcoding `NULL`, guarded factor-model reparameterization with `ncov_psi > 0 & gt > 0`).
- **Working tree**: clean, `main` up to date with `origin/main` as of this session.
- **TODO.Rmd structure** (current, see file for full detail): organized as **v0.1.0-beta Public release** (1. Alex to dos: ordination table/plot, `computePredictiveOccupancyProbs()`, scalar p/q indexing bug, trait-reading fragility, `computeMinESS()` bug; 2. Doug to dos: purge `traitdata_caiwang.rdata` from git history, port model diagnostics functions from GLGS-eDNA repo), **MEE paper** (1. Doug to dos: reproduce Ecology Letters results as a package test; 2. Alex to dos: Overleaf math vignette, GAMs for JSDM, auto-calculate `gt`, count-data support, M-based `plotCumulativeSpeciesDetections()`, source-sink inference scenario, remove space's effect on env covariates, site-level variance/"variation" partitioning), and **Future versions** (spike-in abundance-change estimation, model selection via regularisation/shrinkage e.g. for geospatial foundation model embeddings).

## Notes

- When editing files that are also open in the RStudio editor, prefer re-reading the file immediately before and after edits -- a prior session saw a file get corrupted/duplicated after an `edit` call, likely from an editor/disk desync; rewriting the whole file with `write` fixed it.
- This repo's git history shows a pattern of committing by concern (e.g. testthat setup, generated docs, vignette content fixes kept separate from unrelated whitespace-only diffs) rather than one large commit -- worth continuing when making multiple unrelated changes.
- TODO.Rmd tracks near-term and future development priorities; keep synchronized with actual work being done.

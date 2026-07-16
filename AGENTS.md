# occPlus

R package for fitting occupancy models to eDNA metabarcoding reads data (joint species distribution model / JSDM framework), accounting for a two-stage detection process (field collection + PCR amplification/lab detection), species traits, and optional spatial autocorrelation.

## Package layout

- `R/runOccPlus.R` -- main model-fitting entry point (`runOccPlus()`), plus data-prep helpers (`process_covariates()`, `createDataIdx()`, `get_param()`).
- `R/jsdmfun.R` -- core JSDM machinery, including `simulateData()` (the lower-level data simulator) and coefficient/variance-partitioning helpers (`computeBtcoef()`, `computePsiCoef()`, `computeVariancePartitioning()`).
- `R/simulateData.R` -- exports `simulateOccPlusData()` (the higher-level wrapper that simulates a full occJSDM-style dataset -- occupancy + two-stage detection + read counts -- from `simulateData()`) and `toRunOccPlusFormat()` (converts `simulateOccPlusData()` output into the `info`/`OTU` list format `runOccPlus()` expects, i.e. the same shape as `sampledata`).
- `R/mcmcfun.R`, `R/output.R`, `R/diagnostics.R` -- MCMC sampling and post-processing/output/plotting functions used after `runOccPlus()` (e.g. `returnOccupancyCovariates()`, `plotVariancePartitioning()`, `returnLatentPresences()`, `plotResidualCorrelationMatrix()`; see NAMESPACE for the full exported list).
- `src/jsdm.cpp`, `src/functions.cpp` (via `RcppExports.R`) -- Rcpp/Armadillo backend, including the spatial kernel `K2()` (squared-exponential GP kernel over site coordinates) and `sample_w_cpp()` (samples latent collection state `w` from continuous read intensities when `threshold = 0`).
- `vignettes/occJSDM.Rmd` -- walkthrough of fitting a model with `runOccPlus()` and using the output/plotting functions.
- `vignettes/simulateOccPlusData.Rmd` -- walkthrough of `simulateOccPlusData()`: parameter lists, simulating and converting to `runOccPlus()` format via `toRunOccPlusFormat()`, simulating with/without spatial autocorrelation via `useSpatField`, checking variance partitioning, visualizing occupancy and the linear predictor spatially, and how trait covariates (`Tr`) shape species-specific occupancy coefficients via the `G` matrix.
- `TODO.Rmd` -- outstanding feature list (e.g. NA handling in `data$OTU`, model selection via `loo`, filling in `spatCovariates` and `traitsMatrix` in the vignette, better documenting the `threshold` parameter, MCMC testthat tests).
- `analysis/analysis.R` -- ad hoc analysis script (not part of package build).
- `tests/testthat/` -- unit tests (testthat edition 3, set up via `usethis::use_testthat(3)`); currently covers `toRunOccPlusFormat()` (dimensions, column names, intercept handling, mismatched-dimension errors).
- `data/` -- `sampledata.rda` / `sampledata_orig.rda` (example dataset used in `occJSDM.Rmd`), `sampleresults.rda`, `data_out_env_trait.rdata`, and `CaiWang_data/` (a data directory the user is actively reorganizing; `traitdata_caiwang.rdata` was recently removed in favor of this).

## Data structures

**Model input (`runOccPlus()`)**: a list with `info` (data.frame, one row per PCR replicate, with `Site`, `Sample`, `Primer`, occupancy covariates `X_psi.*`, and collection covariates `X_theta.*`) and `OTU` (matrix of read counts, rows matching `info`, columns = species).

The `threshold` argument to `runOccPlus()` controls how `OTU` is interpreted: - `threshold > 0` (default `1`): reads are truncated to binary presence/absence (`OTU >= threshold` -\> detection), and `w` is sampled via `sample_w_cim_cipp()`. - `threshold == 0`: reads are modeled continuously via `logy1 = log(OTU + 1)` and a two-component Normal mixture -- given a true detection, `logy1 ~ Normal(mu1, sigma1)`; given a false-positive/contamination detection, `logy1 ~ Normal(mu0, sigma0)` -- sampled via `sample_w_cpp()`. `mu1_output`/`sigma1_output`/ `mu0_output`/`sigma0_output` in `results_output` are only populated in this mode.

**Simulated data (`simulateOccPlusData()`)** returns `list(true_params, data_list)`:

- `data_list`: `X_psi` (occupancy covariates), `X_theta` (collection covariates), `Xs` (site coordinates, `n` x 2 -- always simulated, only affects occupancy if `useSpatField = TRUE`), `Tr` (species traits, `S` x `g` -- always affects occupancy via the trait-response matrix `G`), `y` (simulated **read counts**, matching the `threshold = 0` mode of `runOccPlus()` -- not binary detections; see below).
- `true_params`: true occupancy/detection states and coefficients, including `jsdmParams_true$varPart` (per-species Environmental/Spatial/Biotic variance partitioning), `jsdmParams_true$eta` (linear predictor, useful for inspecting the spatial field without the noise added by thresholding to binary occupancy), and `mu1_true`/`sigma1_true`/`mu0_true`/`sigma0_true` (the true read-intensity parameters used to generate `y`).

`simulateOccPlusData()` takes `list_datasettings` (`n`, `S`, `g`, `M`, `P`, `K`, `ncov_psi`, `ncov_theta`), `list_params` (`p`, `q`, `theta0`, `theta_baseline`, and optionally `mu1`, `sigma1`, `mu0`, `sigma0` -- defaulting to `5`, `1`, `1.5`, `1` -- controlling simulated read intensity), `list_jsdmParams` (`gt`, `d`, `ds`, `sigma_b`, `sigma_bs`, `sigma_ts`, `sigma_h`, `l_s`), and `useSpatField` (default `FALSE`). Read counts are generated as `y = round(exp(logy1) - 1)` where `logy1` is drawn per-replicate from `Normal(mu1, sigma1)` (true detections), `Normal(mu0, sigma0)` (false-positive/contamination), or `0` (no reads). With the defaults, this produces a distribution roughly similar in shape to `sampledata$OTU` (\~two-thirds zeros; nonzero median around 120), though the upper tail runs somewhat heavier.

**`toRunOccPlusFormat(sim, n, M, P, K, drop_theta_intercept = TRUE)`** expands `sim$data_list$X_psi`/`X_theta` (site/sample-level) to the PCR-replicate level (via the same indexing as `createDataIdx()`), builds `Site`/`Sample`/`Primer` id columns, and returns `list(info, OTU = sim$data_list$y)`. `n`/`M`/`P`/`K` must match the `list_datasettings` originally used to create `sim` -- they are not stored in `sim` itself, so passing inconsistent values will error or silently mis-index. `X_theta`'s intercept column (always column 1) is dropped by default since `runOccPlus()` adds its own via `process_covariates()`.

## Known issues in existing (pre-`toRunOccPlusFormat`) code

- `predictOccupancyProbs()` (`R/output.R`) is currently broken: it references an undefined `X_ord` object instead of its own `X_s` argument, and errors even earlier because `nchain`/`niter` are derived from `beta_psi_output`'s dimensions, which don't have the expected shape on at least some fitted models. Confirmed broken via live testing; not fixed as of this session.
- `returnLatentPresences()` (`R/output.R`) references an undefined `varPart_output` object (likely leftover from a copy-paste) and will error at runtime.
- Both bugs already have `@note` roxygen comments flagging them in the source.
- `computeAverageCollectionProbs()` and `computeConditionalSamplePresenceProbs()` were confirmed working (via live testing against fitted model objects) as long as the fitted model actually populated `results_output$theta_output` / `results_output$w_output` (always true for a fresh `runOccPlus()` fit with default `summarisedLatentPresences = TRUE`).

## `man/*.Rd` tracking

`man/*.Rd` files are currently tracked and committed in git (not gitignored), despite an earlier commit (`dd158a9`, prior to this session) whose message claimed intent to "stop generating/tracking man/\*.Rd" -- that intent was never followed through in `.gitignore`. Regenerating docs via `devtools::document()` after adding/changing roxygen tags will produce `man/*.Rd` diffs that should be committed (or the gitignore situation resolved) deliberately.

## Notes

- When editing files that are also open in the RStudio editor, prefer re-reading the file immediately before and after edits -- a prior session saw a file get corrupted/duplicated after an `edit` call, likely from an editor/disk desync; rewriting the whole file with `write` fixed it.
- This repo's git history shows a pattern of committing by concern (e.g. testthat setup, generated docs, vignette content fixes kept separate from unrelated whitespace-only diffs) rather than one large commit -- worth continuing when making multiple unrelated changes.

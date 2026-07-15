# occPlus

R package for fitting occupancy models to eDNA metabarcoding reads data
(joint species distribution model / JSDM framework), accounting for a
two-stage detection process (field collection + PCR amplification/lab
detection), species traits, and optional spatial autocorrelation.

## Package layout

- `R/runOccPlus.R` -- main model-fitting entry point (`runOccPlus()`),
  plus data-prep helpers (`process_covariates()`, `createDataIdx()`).
- `R/jsdmfun.R` -- core JSDM machinery, including `simulateData()` (the
  lower-level data simulator) and coefficient/variance-partitioning
  helpers (`computeBtcoef()`, `computePsiCoef()`,
  `computeVariancePartitioning()`).
- `R/simulateData.R` -- `simulateOccPlusData()`, the higher-level
  wrapper that simulates a full occJSDM-style dataset (occupancy +
  two-stage detection) from `simulateData()`. Not yet exported in
  NAMESPACE (only available via `devtools::load_all()`); run
  `devtools::document()` if you add `@export`.
- `R/mcmcfun.R`, `R/output.R`, `R/diagnostics.R` -- MCMC sampling and
  post-processing/output/plotting functions used after
  `runOccPlus()` (e.g. `returnOccupancyCovariates()`,
  `plotVariancePartitioning()`, `returnLatentPresences()`,
  `plotResidualCorrelationMatrix()`; see NAMESPACE for the full
  exported list).
- `src/jsdm.cpp` (via `RcppExports.R`) -- Rcpp/Armadillo backend,
  including the spatial kernel `K2()` (squared-exponential GP kernel
  over site coordinates).
- `vignettes/occJSDM.Rmd` -- walkthrough of fitting a model with
  `runOccPlus()` and using the output/plotting functions.
- `vignettes/simulateOccPlusData.Rmd` -- walkthrough of
  `simulateOccPlusData()`: parameter lists, simulating with/without
  spatial autocorrelation via `useSpatField`, checking variance
  partitioning, visualizing occupancy and the linear predictor
  spatially, and how trait covariates (`Tr`) shape species-specific
  occupancy coefficients via the `G` matrix.
- `TODO.Rmd` -- outstanding feature list (e.g. NA handling in
  `data$OTU`, model selection via `loo`, filling in `spatCovariates`
  and `traitsMatrix` in the vignette).
- `analysis/analysis.R` -- ad hoc analysis script (not part of package
  build).
- `data/` -- `sampledata.rda` / `sampledata_orig.rda` (example dataset
  used in `occJSDM.Rmd`), `sampleresults.rda`, plus
  `data_out_env_trait.rdata` and `data_todoug_20260714.csv`.

## Data structures

**Model input (`runOccPlus()`)**: a list with `info` (data.frame, one
row per PCR replicate, with `Site`, `Sample`, `Primer`, occupancy
covariates `X_psi.*`, and collection covariates `X_theta.*`) and `OTU`
(matrix of read counts, rows matching `info`, columns = species).

**Simulated data (`simulateOccPlusData()`)** returns
`list(true_params, data_list)`:

- `data_list`: `X_psi` (occupancy covariates), `X_theta` (collection
  covariates), `Xs` (site coordinates, `n` x 2 -- always simulated,
  only affects occupancy if `useSpatField = TRUE`), `Tr` (species
  traits, `S` x `g` -- always affects occupancy via the trait-response
  matrix `G`), `y` (simulated detections).
- `true_params`: true occupancy/detection states and coefficients,
  including `jsdmParams_true$varPart` (per-species
  Environmental/Spatial/Biotic variance partitioning) and
  `jsdmParams_true$eta` (linear predictor, useful for inspecting the
  spatial field without the noise added by thresholding to binary
  occupancy).

`simulateOccPlusData()` takes `list_datasettings` (`n`, `S`, `g`, `M`,
`P`, `K`, `ncov_psi`, `ncov_theta`), `list_params` (`p`, `q`, `theta0`,
`theta_baseline`), `list_jsdmParams` (`gt`, `d`, `ds`, `sigma_b`,
`sigma_bs`, `sigma_ts`, `sigma_h`, `l_s`), and `useSpatField` (default
`FALSE`).

## Notes

- When editing files that are also open in the RStudio editor, prefer
  re-reading the file immediately before and after edits -- a prior
  session saw a file get corrupted/duplicated after an `edit` call,
  likely from an editor/disk desync; rewriting the whole file with
  `write` fixed it.

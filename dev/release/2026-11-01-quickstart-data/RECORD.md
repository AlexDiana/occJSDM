# Quickstart example data and fit: record

1 November 2026. Design: `dev/superpowers/specs/2026-11-01-quickstart-expansion-design.md`.

## Files

- `simulate-quickstartdata.R` writes `data/quickstartdata.rda` and `truth.rds` (planted values; not shipped, not used by the Quickstart).
- `fit-quickstartresults.R NBURN NITER NTHIN OUT.rds` runs the Quickstart's live-fit call and saves the fit, its warnings and provenance.
- `slim-quickstartresults.R FIT.rds 1` saves the chosen fit as `data/quickstartresults.rda` and writes `diagnostics.csv`, `warnings.txt` and `provenance.rds`. Its thinning option was not needed: `runOccJSDM()` thins itself (`nthin`), keeping `niter` draws per chain, so the shipped object is exactly the live-fit call.

## Runs

All on the M4 MacBook Air, R 4.6.1, Accelerate BLAS, AC power, Low Power Mode off, package at `aca06cb` (clean `R/`, `src/`, `data/`), seed 20261101, 2 chains on 2 cores. Fits saved in `~/occjsdm-runs/quickstart/` (not in the repository). Flags use the Quickstart's rule, `rhat > 1.01` or `ess < 400`, over the 220 rows of `returnConvergenceDiagnostics()`.

- 1,000 burn-in + 1,000 kept, `nthin = 1`: 0.3 minutes, 76 flagged; worst `theta0` of `Mammalia_Rodentia_2` (R-hat 1.04, ESS 106). `runOccJSDM()` warned that `sigmah` has R-hat above 1.1.
- 10,000 + 10,000, `nthin = 1`: 4.1 minutes, 0 flagged, no warnings. Compresses to about 81 MB, too large to ship.
- 20,000 + 20,000, `nthin = 1`: 8 minutes, 1 flagged (`beta0_psi` of `Aves_Passeriformes_2`, R-hat 1.01, ESS 923), no warnings.
- 2,000 + 2,000, `nthin = 2`: 1.2 minutes, 9 flagged, no warnings.
- **2,000 + 2,000, `nthin = 3`: 1.4 minutes, 3 flagged, no warnings. Shipped.** 16.3 MB.

## The shipped fit's flags

All three involve the two species planted with field contamination (`theta0 = 0.12`, against 0.01 for the rest), which also have low true occupancy (0.15 and 0.28): with few occupied sites, a detection is hard to attribute to presence or to contamination. All have R-hat of about 1.00; ESS is what fails.

- `beta0_psi` (baseline occupancy) of `Aves_Passeriformes_2`, ESS 280.
- `theta0` (field false-positive rate) of `Mammalia_Rodentia_2`, ESS 384.
- `beta_theta` intercept (collection) of `Mammalia_Rodentia_2`, ESS 380.

The Quickstart features `theta0` of `Mammalia_Rodentia_2`. A longer run (10,000 + 10,000) cleared every flag.

## Names in the fit

The factor's coefficient is labelled `protectedprotected` (column name plus level, from `model.matrix()`); pass that as `covName` where a function needs it.

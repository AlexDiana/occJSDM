# Quickstart example data and fit: record

1 November 2026. Design: `dev/superpowers/specs/2026-11-01-quickstart-expansion-design.md`.

## Files

- `simulate-quickstartdata.R` writes `data/quickstartdata.rda` and `truth.rds` (planted values; not shipped, not used by the Quickstart).
- `fit-quickstartresults.R NBURN NITER NTHIN OUT.rds` runs the Quickstart's live-fit call and saves the fit, its warnings and provenance.
- `slim-quickstartresults.R FIT.rds 1` saves the chosen fit as `data/quickstartresults.rda` and writes `diagnostics.csv`, `warnings.txt` and `provenance.rds`. Its thinning option was not needed: `runOccJSDM()` thins itself (`nthin`), keeping `niter` draws per chain, so the shipped object is exactly the live-fit call.

## Runs

All on the M4 MacBook Air, R 4.6.1, Accelerate BLAS, AC power, Low Power Mode off, seed 20261101, 2 chains on 2 cores. Fits saved in `~/occjsdm-runs/quickstart/` (not in the repository). Flags use the Quickstart's rule, `rhat > 1.01` or `ess < 400`, over the 220 rows of `returnConvergenceDiagnostics()`.

First version of the data (body mass in kg, Red List effect on protection 1.4, hidden-covariate loadings SD 0.7):

- 1,000 burn-in + 1,000 kept, `nthin = 1`: 0.3 minutes, 76 flagged. `runOccJSDM()` warned that `sigmah` has R-hat above 1.1.
- 10,000 + 10,000, `nthin = 1`: 4.1 minutes, 0 flagged. Compresses to about 81 MB, too large to ship.
- 20,000 + 20,000, `nthin = 1`: 8 minutes, 1 flagged.
- 2,000 + 2,000, `nthin = 2`: 1.2 minutes, 9 flagged; `nthin = 3`: 1.4 minutes, 3 flagged.

Three changes to the simulation followed, each to make a figure tell its story; each was refitted at 2,000 + 2,000:

- Body mass stored as log10 (raw kg, from 0.01 to 225, let two ungulates dominate the standardised trait). `nthin = 3`: 7 flagged.
- Red List effect on protection raised from 1.4 to 2.0 (its trait-by-environment interval had crossed zero). `nthin = 3`: 6 flagged.
- Hidden-covariate loadings SD raised from 0.7 to 1.5 (no species pair was clearly associated). `nthin = 3`: 1 flagged; `nthin = 1`: 43 flagged; **`nthin = 2`: 1.2 minutes, 4 flagged, no warnings. Shipped**, 16.5 MB.

On the shipped data, 10,000 + 10,000 with `nthin = 1` took 3.4 minutes and flagged nothing, with no warnings; the Quickstart says so.

## The shipped fit's flags

All four fail on R-hat (1.01 to 1.04); every ESS is above 400.

- `theta0` (field false-positive rate) of `Mammalia_Artiodactyla_2`, R-hat 1.04. The Quickstart features it: one chain makes more and higher excursions than the other. The species is estimated to occupy about 80 of the 100 sites, leaving few unoccupied sites to inform its field false-positive rate.
- `beta0_psi` of `Mammalia_Carnivora_1` and of `Mammalia_Eulipotyphla_2`, and `beta_psi` for `hmi` of `Mammalia_Artiodactyla_2`, R-hat 1.01.

In the shipped fit, the trait-by-environment coefficients are -0.55 (95% interval -0.91 to -0.22) for log10 body mass on `hmi` and 0.81 (0.25 to 1.38) for Red List status on protection, and 76 of the 190 species pairs have residual correlations whose interval excludes zero.

## Names in the fit

The factor's coefficient is labelled `protectedprotected` (column name plus level, from `model.matrix()`); pass that as `covName` where a function needs it.

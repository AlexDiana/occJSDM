# Rare-dominated communities, 15 October 2026

The R = 300 paper simulation study found species-level occupancy compressed towards the middle in the two-stage cells: in `base`, species with true occupancy 0.1 to 0.2 were estimated at about 0.35, and species at 0.8 to 0.9 at about 0.62. The simulator draws true B0 from Normal(0, 1), the fitter's own fixed B0 prior, so fewer than 3% of the study's species have true occupancy under 0.2, and over a whole community the upward pull on rare species and the downward pull on common ones largely cancel. Real eDNA communities are mostly rare species. This run measures the bias when they are.

## Design

Three new arms, each a copy of a paper-study cell with the mean of the true B0 draw moved from 0 to -2.75:

-   `base_rare`, paired with `base` (two-stage).
-   `occupancy_rare`, paired with `occupancy` (one-stage).
-   `binary_rare`, paired with `binary` (occupancy observed directly).

The parents are replicates 1 to 100 of the R = 300 paper study, not rerun. The new arms use 100 replicates each with the study runner, 5 replicate processes and one thread per fit, at the runner's default 1000 burn-in and 1000 sampling iterations with 2 chains, and the fitter's default B0 prior, Normal(0, 1).

## Why -2.75

The brief suggested a mean of -1.5, expecting a median true occupancy of about plogis(-1.5) = 0.18. The realised occupancy is higher than the occupancy at the intercept, because the covariates, latent factors and spatial field spread each species' linear predictor across sites, and averaging plogis over that spread lifts a rare species. Measured over replicates 1 to 100 of `base` and `binary` without fitting, a mean of -1.5 gives a median true occupancy of 0.31 with about 21% of species under 0.2, and -2.75 gives a median of 0.175 to 0.19 with about 57% under 0.2. The arms use -2.75, the shift that gives the community the brief describes. In the same check, 216 of `base_rare`'s 1,000 species and 247 of `binary_rare`'s had true occupancy of 0.1 or less, against 6 in `base` and 2 in `binary`.

## The spatial field is kept

The brief allowed dropping the spatial field if the [spatial 2x2 report](../spatial-field-2x2-20261015/REPORT.md) supported it. It does not, for three reasons. The 2x2 found no measurable effect of the field on B0 in `base`, but it did not run `binary`, and the range test the same evening found B0 coverage in `binary` partly recovered when the simulated range was moved inside the fitted grid, so B0 there does depend on the spatial part. Dropping the field would break the pairing: the field's own draws shift the random stream before occupancy is simulated, so the datasets would no longer differ from the parents' only through B0. And it would save little: the 2x2's three arms, two of them fitting no spatial term, ran at 6.7 fits a minute, about the paper study's rate.

## Pairing

The simulator draws B0 inside the package with no argument to move it. The harness setting `b0_mean` builds a copy of the simulator whose B0 line adds the shift to the same Normal(0, 1) draw, without changing package code or the package namespace (`simstudy_shifted_simulator()` in `tests/testthat/helper-simstudy.R`). Each arm takes its parent's seeds, so at every replicate its true B0 is the parent's minus 2.75, element for element, and every other part of the truth is identical: covariates, traits, B, G, the latent factors, coordinates, the spatial field, the collection covariates, `beta_theta`, `theta0`, `p` and `q`. The occupancy states, and so the data, differ. `tests/testthat/test-simstudy-b0-shift.R` checks this at replicate 1 for all three pairs, checks that a zero shift reproduces the package simulator bit for bit, and checks that cells without `b0_mean` simulate exactly as before.

## Run

Launched detached from this directory with:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh > launch.log 2>&1 < /dev/null &
```

Monitor with `tail -f run.log`. The log, provenance and results are ignored while the run is in flight so the runner records a clean tree.

## Results

The run finished on 16 October 2026 at 00:38 with all 300 fits and no failures. In the two eDNA cells, species with true occupancy under 0.1 are estimated at about 0.34 whatever their truth, and the community's mean occupancy at 0.37 against a true 0.20; at a given true occupancy the bias matches the parent cells'. The findings, and a question for Alex on estimating the community mean of B0, are in [REPORT.md](REPORT.md).

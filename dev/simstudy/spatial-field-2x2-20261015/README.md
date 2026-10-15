# The spatial 2x2, 15 October 2026

The R = 300 paper simulation study ([its calibration report](../paper-study-20261015/REPORT.md)) found the B intervals too narrow in every cell, worst for the largest true effects. Doug's hypothesis is that the simulated spatial field, which the fitted spatial term cannot fully represent, attenuates B. This run separates that from the effect of fitting a spatial term at all, on the base cell's settings.

## Design

Four arms, crossing whether a spatial field is simulated with whether a spatial term is fitted:

-   `base`: field simulated, spatial term fitted. Replicates 1 to 100 of the R = 300 paper study, not rerun.
-   `field_unfitted`: field simulated, no spatial term fitted. New.
-   `nofield_fitted`: no field simulated, spatial term fitted. New.
-   `no_spatial`: neither. Defined in the harness since August, never run before.

The three new runs use 100 replicates each with the study runner, 5 replicate processes and one thread per fit, at the runner's default 1000 burn-in and 1000 sampling iterations with 2 chains. The harness gains a `fit_spatial` setting that decouples fitting from simulation; production cells leave it unset and fit exactly as before, which a refit of four production cells confirmed bit for bit.

## Pairing

All four arms draw `base`'s seeds. The simulator draws the whole JSDM layer before the spatial field, so every arm has the same true B0, B, G, A, C and L. `field_unfitted` has exactly `base`'s datasets, and `nofield_fitted` exactly `no_spatial`'s, so those two comparisons differ only in the fit. Between field and no-field arms the datasets differ, because the field's own draws shift the random stream before occupancy and detection are simulated. The harness tests in `tests/testthat/test-simstudy-spatial-2x2.R` check this.

## Run

Launched detached from this directory with:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh > launch.log 2>&1 < /dev/null &
```

Monitor with `tail -f run.log`. The log, provenance and results are ignored while the run is in flight so the runner records a clean tree.

## Results

Filled in after the run.

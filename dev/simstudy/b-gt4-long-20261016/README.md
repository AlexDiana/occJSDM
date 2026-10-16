# 4 latent traits with long chains, 16 October 2026

The non-spatial candidates run ([its report](../b-nonspatial-candidates-20261015/REPORT.md)) found two levers on B's 95% coverage in `no_spatial`, each tested alone on identical datasets: fitting 4 latent traits instead of 2 raised it from 0.906 to 0.936 (paired +0.030, SE 0.006), and 5,000 + 5,000 iterations per chain instead of 1,000 + 1,000 raised it to 0.923 (+0.017, SE 0.005). This run fits the two together, to see whether they bring B to nominal and whether their effects add.

## Design

One cell, `no_spatial_gt4_long`, at 100 replicates and seed label `base`: `no_spatial`'s datasets fitted with 4 latent traits and 5,000 burn-in and 5,000 sampling iterations per chain, 2 chains. It is defined in `tests/testthat/helper-simstudy.R` beside `no_spatial_gt4` and `no_spatial_long`, and its own MCMC settings replace the runner's. It pairs on identical data with three arms already run at revision `c40f385`, whose package code is unchanged here: `no_spatial` (2 latent traits, 1,000 + 1,000), `no_spatial_gt4` (4 latent traits, 1,000 + 1,000) and `no_spatial_long` (2 latent traits, 5,000 + 5,000). Those three are read from `../b-nonspatial-candidates-20261015/results/simstudy-20261016-004116.rds`, not rerun.

The harness test in `tests/testthat/test-simstudy-b-candidates.R` checks that the cell hands the fitter 4 latent traits and 2 factors, carries `no_spatial_long`'s MCMC settings, and fits exactly `no_spatial`'s data at replicate 1.

## Run

The run uses the study runner with `--pkg`, so the workers draw L'Ecuyer streams from the runner's seed and every replicate simulates the same dataset as in the earlier run. 5 replicate processes, one thread per fit. Launched detached from this directory with:

```sh
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' sh launch.sh > launch.log 2>&1 < /dev/null &
```

Monitor with `tail -f run.log`. The log, provenance and results are ignored while the run is in flight so the runner records a clean tree.

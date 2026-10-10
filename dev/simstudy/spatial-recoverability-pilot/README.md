# Spatial recoverability pilot

An oracle-only positive control for the spatial representation question. It asks in which sampling regimes spatial structure is recoverable from single-season binary occupancy data, for two kinds of spatial truth: a few common latent gradients shared by all species, and an independent field per species whose range is set by a dispersal trait. No package is fitted. The oracle knows every generating parameter except the latent fields, estimates the fields by elliptical slice sampling, and predicts a held-out 40 by 40 lattice by Gaussian-process conditioning. The specification is `dev/superpowers/specs/2026-10-10-spatial-recoverability-pilot-design.md` and the implementation plan is `dev/superpowers/plans/2026-10-10-spatial-recoverability-pilot.md`.

## Results

See `REPORT.md`. Compact results, the verification record and the recoverability map are in `results/`. Raw job files, with thinned site-field draws and convergence traces, are on the research volume at `/Volumes/occJSDM Research/occJSDM/spatial-recoverability-20261010/jobs`; the run index, inputs and logs are in the ignored `dev/simstudy/results/spatial-recoverability-20261010/`.

A single-season snapshot cannot attribute species-specific spatial structure to dispersal as a mechanism. Every result here is a statement about the statistical footprint, not the cause.

## Reproduction

All commands run from the repository root through the pinned R 4.5.0 launcher `run-r`, never system `Rscript`. The launcher points at the private R home built for the first spatial comparison and forces single-threaded numerics.

```sh
P=dev/simstudy/spatial-recoverability-pilot
for t in seeds kernel generator oracle score run; do $P/run-r -f $P/test-$t.R || exit 1; done
STUDY=$PWD/dev/simstudy/results/spatial-recoverability-20261010
JOBS="/Volumes/occJSDM Research/occJSDM/spatial-recoverability-20261010/jobs"
$P/run-r -f $P/run.R --args --mode=prepare --study=$STUDY --jobs="$JOBS" --workers=8
$P/run-r -f $P/run.R --args --mode=timing --study=$STUDY
nohup $P/run-r -f $P/run.R --args --mode=fit --study=$STUDY --truth=t2 --workers=8 > $STUDY/logs/fit-t2.log 2>&1 &
nohup $P/run-r -f $P/run.R --args --mode=fit --study=$STUDY --truth=t1 --workers=8 > $STUDY/logs/fit-t1.log 2>&1 &
$P/run-r -f $P/run.R --args --mode=score --study=$STUDY
$P/run-r -f $P/verify.R --args --study=$STUDY
$P/run-r -f $P/plot.R --args --study=$STUDY
```

`prepare` refuses to run twice on the same study directory, writes the thresholds and budget to `config.rds` before generating any community, and records an MD5 for every input. `fit` never refits a job whose file exists. Each `nohup` step must finish before the next starts; the logs end with a completion line.

## Files

- `run-r`: the pinned launcher.
- `seeds.R`: label-hashed seed namespace, distinct from the sweep's.
- `kernel.R`: squared-exponential kernel, stable factorisation with Cholesky or truncated-eigen route, lattice conditional, Gauss-Hermite quadrature.
- `generator.R`: species table, arrangements with scaled cluster separation, both truths, intercept solving.
- `joint-ellipse.cpp`: joint elliptical slice sampler over two shared fields and twelve independent parts; reduces exactly to the sweep's `ellipse.cpp` for one species and one field.
- `oracle.R`: truth 2 per-species oracle (the sweep's sampler, unmodified), truth 1 joint oracle, convergence screen.
- `predict.R`: lattice prediction, uninformed reference, layout geometry.
- `score.R`: per-species scores, aggregation, gate and map reading.
- `run.R`: driver with modes prepare, timing, fit and score.
- `verify.R`: independent recomputation of every gate decision and the map, plus input, seed and conditioning checks.
- `plot.R`: copies compact results and draws `results/recoverability-map.png`.
- `test-*.R`: the validation tests, one per module.

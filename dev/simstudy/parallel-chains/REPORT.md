# Parallel chains on the quickstart fit: timings, memory, identity and BLAS provenance

Measured on 13 October 2026. The question is what `MCMCparams$ncores` buys on the quickstart fit, what it costs in memory, and whether the fit and the random-number stream stay identical whatever the worker count.

## What was run

- The fit is the quickstart call on the shipped `sampledata`, with four chains of 5,000 burn-in and 5,000 kept iterations and no thinning, preceded by `set.seed(20261007)` each time. Only `ncores` changes between configurations, apart from the thread pin in `serial-1thread`.
- Five configurations, each run twice: `serial-ambient` (`ncores` 1, the in-process fit with the ambient thread count, which was 10), `serial-1thread` (`ncores` 1 with `RCPP_PARALLEL_NUM_THREADS` set to 1), `ncores-2`, `ncores-4` and `ncores-8`. With four chains, `ncores` 8 gives four workers of two threads each.
- Machine: Apple M4 laptop (Mac16,13) with ten cores (four performance and six efficiency) and 32 GB of memory, macOS 27.0.1, on mains power, R 4.6.1 (2026-06-24), system R at `/usr/local/bin/Rscript`.
- The package was frozen before measuring: `launch.sh` exports `HEAD` with `git archive`, installs the export into a dedicated library under the run directory, and the benchmark loads it from there. The main run's `revision.txt` names `dd5b1b7a8512fcc57100ae3ddaefa06cbc56427d`. That commit was later amended, for its trailer only, to `8780cb9` with the same tree. The main run started at about 21:41 and finished at 22:16:42 (36 minutes).
- The rerun of the two ambient serial fits (below) installed its own export; its `revision.txt` names `b556db5bdba221e268712bd3f4985a10c626a781`, the commit that added the configuration selector. `R/`, `src/`, `DESCRIPTION` and `NAMESPACE` are unchanged between the two revisions and since `f9fc489`. The two exports were compared: `R/`, the `src/` sources, `DESCRIPTION` and `NAMESPACE` are identical (only the compiled objects differ, because they embed the build path).
- Memory is the resident set size in KB sampled once a second by `sample-rss.sh`: the parent R process, and every process whose command line contains `workRSOCK`. `summarise.R` keeps only the worker PIDs that the fit's own execution record names.
- Identity: each fit, with its execution record dropped, is compared with `identical()` to the first fit of its own run (the in-process `serial-ambient` repeat 1). Identical RNG means the parent session's `.Random.seed` after the fit, and the next draws from `rnorm`, `runif` and `sample`, equal those after the reference fit.
- The quickstart call, as in `benchmark.R`:

```r
runOccJSDM(data = sampledata, listParams = list(n_factors = 2), threshold = 1,
           occCovariates = c("X_psi.EnvCov.1", "X_psi.EnvCov.2"),
           collCovariates = c("X_theta.1", "X_theta.2"),
           spatCovariates = c("Xs.1", "Xs.2"),
           MCMCparams = list(nchain = 4, nburn = 5000, niter = 5000, nthin = 1, ncores = ncores),
           summarisedLatentPresences = TRUE)
```

## BLAS evidence and declaration

The facts below come from `blas-evidence.rds` (collected by `blas-provenance.R` at 20:41 UTC on 13 October 2026, from the main run). The rerun collected the evidence again at 21:18 UTC and it is identical apart from the timestamp.

- The BLAS that R loads, from `extSoftVersion()`, is the file below. LAPACK is `libRlapack.dylib` in the same folder (LAPACK 3.12.1).

```
/Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib
```

- It is a regular file, not a symlink. It is 209,392 bytes, modified 2026-06-25 22:58:14 (BST), owned by root with group admin, a Mach-O arm64 shared library.
- Its md5 is `9dfbe56d4ceb13138673fdaf297b4b56`.
- Its linkage, from `otool -L`, is `libgfortran.5.dylib`, `libquadmath.0.dylib` and `/usr/lib/libSystem.B.dylib` (all but the last under the same `R.framework` folder). The first line `otool -L` prints is the library's own install name, `libRblas.dylib`, which `otool -D` confirms, so it is not a dependency. No library name mentions Accelerate, vecLib, OpenBLAS, MKL, BLIS or ATLAS.
- It defines `dgemm_` itself (`nm -gU` lists `T _dgemm_`).
- The installer receipt `org.R-project.R.fw.pkg` is present: version 4.6.1, volume `/`, location `Library/Frameworks`, installed at 2026-10-07 13:05:43 UTC. The receipt lists `libRblas.0.dylib`. The same receipt also lists a vecLib alternative, `libRblas.vecLib.dylib`, which is not the library loaded here.
- The compatibility symlink `libRblas.dylib` in the same folder points at `libRblas.0.dylib` (it was modified on 7 October 2026), so the unversioned name resolves to the reference library too.
- Declaration: `OCCJSDM_REFERENCE_BLAS` was set to 1 in the launching environment, with the basis text recorded in `blas-declaration-basis.txt`, quoted here.

```
R 4.6.1 installed by the CRAN installer package org.R-project.R.fw.pkg, whose receipt lists libRblas.0.dylib; the file is regular, not a symlink, links only to libgfortran, libquadmath and libSystem, defines dgemm_ itself, md5 9dfbe56d4ceb13138673fdaf297b4b56; declared by Doug Yu and Claude on 2026-10-13
```

- The declaration was made by a person on this evidence, after the evidence matched what had been established on the same machine earlier that day; the scripts make no verdict. Every execution record (all twelve fits) has `blas_declared` TRUE and names the same library.
- The identity claims below hold under this declaration, that is, with the reference BLAS and no threaded BLAS. They are not claims about any other BLAS.

## Results

Main run, from `summarise.R`. Speed-up is the mean `serial-ambient` wall time over the mean wall time of the configuration. Each worker's peak is listed for both repeats.

- **serial-ambient** (ncores 1, 1 worker, 10 threads each): wall 324 s and 318 s (speed-up 1.00x over serial-ambient), parent peak 1002 MB, no worker process (in-process fit), identical fit TRUE/TRUE, identical RNG TRUE/TRUE.
- **serial-1thread** (ncores 1, 1 worker, 1 thread each): wall 329 s and 302 s (speed-up 1.02x over serial-ambient), parent peak 1032 MB, no worker process (in-process fit), identical fit TRUE/TRUE, identical RNG TRUE/TRUE.
- **ncores-2** (ncores 2, 2 workers, 1 thread each): wall 162 s and 168 s (speed-up 1.95x over serial-ambient), parent peak 1087 MB, workers together 566 MB at most, each worker 285/281 MB in repeat 1 and 282/284 MB in repeat 2, identical fit TRUE/TRUE, identical RNG TRUE/TRUE.
- **ncores-4** (ncores 4, 4 workers, 1 thread each): wall 117 s and 114 s (speed-up 2.78x over serial-ambient), parent peak 1059 MB, workers together 1042 MB at most, each worker 260/259/261/262 MB in repeat 1 and 258/257/258/258 MB in repeat 2, identical fit TRUE/TRUE, identical RNG TRUE/TRUE.
- **ncores-8** (ncores 8, 4 workers, 2 threads each): wall 112 s and 115 s (speed-up 2.83x over serial-ambient), parent peak 1059 MB, workers together 1039 MB at most, each worker 258/259/259/262 MB in repeat 1 and 263/258/261/257 MB in repeat 2, identical fit TRUE/TRUE, identical RNG TRUE/TRUE.

Identity held in all ten fits. The parent's memory is about 1 GB whatever the worker count (1002 to 1087 MB), and each worker adds about 260 to 285 MB. The parent and worker peaks are reported separately and are not summed, because the parent peak and the workers' peaks need not fall in the same second.

The speed-ups above use the original serial timings, which the next section revisits. The full per-fit timings, including each chain's own elapsed time, are in `results.csv`.

## The rerun of the ambient serial fits

- Why: Doug closed several applications at the start of the first serial runs, so the two original `serial-ambient` timings may have been taken under more load than the later runs. The rerun measures the same fit again, in a fresh run directory and from a fresh export, under the same conditions as the parallel runs. It ran with `launch.sh "$OUT" serial-ambient`.
- Original: 323.9 s and 318.2 s (mean 321.1 s; chains 78 to 81 s each).
- Rerun: 309.4 s and 305.2 s (mean 307.3 s; chains 75 to 77 s each), 4 percent faster. Both reruns were identical to the rerun's first fit, with identical RNG state, and the rerun's reference fit equals the main run's reference fit (an `identical()` comparison of the two `reference-fit.rds` files, and of the two RNG records, made once after the runs; it is not part of the scripts).
- Speed-ups over the original serial mean (321.1 s): `serial-1thread` 1.02x, `ncores-2` 1.95x, `ncores-4` 2.78x, `ncores-8` 2.83x.
- Speed-ups over the rerun serial mean (307.3 s): `serial-1thread` 0.97x, `ncores-2` 1.87x, `ncores-4` 2.66x, `ncores-8` 2.71x.
- Either baseline gives the same picture: two workers about 1.9 times, four workers about 2.7 times, eight cores no better than four, and pinning the serial fit to one thread makes no measurable difference. The 4 percent gap between the two serial pairs is smaller than the gap between the two repeats of `serial-1thread` (329 s and 302 s, 9 percent), so with two repeats each this run cannot show that load caused the difference.

## The chains-outnumber-workers case

With `ncores` 2 and four chains, each of the two workers runs two chains one after the other.

- Wall times were 162 s and 168 s, 1.95 times faster than the original serial mean and 1.87 times faster than the rerun's.
- The parent peaked at 1087 MB. The two workers peaked at 285 and 281 MB in repeat 1 and 282 and 284 MB in repeat 2, 566 MB together at most. A worker that runs two chains peaks only about 10 percent higher than one that runs a single chain (257 to 263 MB at `ncores` 4 and 8).
- Each chain took 78 to 83 s, the same as a chain run alone, so two workers did not slow each other down.

## Interpretation

This paragraph is interpretation, not result. Four single-threaded workers gave about 2.7 to 2.8 times rather than four times, although `serial-1thread` shows that the threaded parts of the sampler are worth almost nothing, so the shortfall is not lost threading. Eight cores gave no more than four because four chains can occupy only four workers; the extra cores became a second thread per worker, which is worth almost nothing. The per-chain times show where the time goes. A chain takes about 75 to 81 s alone and about 78 s with two workers, but 108 to 111 s when four workers run at once, about 1.4 times longer, and four divided by 1.4 is roughly the 2.8 observed. The cost of starting the workers and collecting the chains is small, about 5 s (wall time minus the slowest chain, at `ncores` 4 and 8). A plausible cause is this laptop's M4, which has four performance and six efficiency cores, so that four concurrent workers, with the parent and the rest of the system, do not all run on performance cores. This is a hypothesis: nothing here measures core placement, and a server with uniform cores may behave differently.

## Limits

- Windows was not run.
- No threaded BLAS was available, so nothing was measured about BLAS threads competing with the workers; the identity claims hold for the declared reference BLAS only.
- Two repeats per configuration, on one machine. The repeats of one configuration differ by up to 9 percent (`serial-1thread`), which is as large as most differences between the serial timings.
- Workers were identified by their command line (`workRSOCK`) and kept only if the fit's recorded PIDs named them; memory is resident set size at one-second sampling, so a peak shorter than a second can be missed.
- The possible load difference that motivated the rerun: the machine is an ordinary working laptop (a desktop session and a backup agent were running when the rerun was launched), and load was not measured in either run.

# Parallel-chains example refit, 13 October 2026

This folder records the refit of the shipped example `data/sampleresults.rda` under per-chain seeding. It changes the shipped example and its evidence, not model code, priors or the input data. The refit is a workflow example whose convergence flags are disclosed. It is not a certified or fully converged fit, and successful scripts do not establish scientific convergence.

This is the second refit of 13 October 2026, redone because the build-fingerprint and in-process error-path fix (`efa40f9`, which this report first cited as `a02d49a`; its commit message trailer was amended afterwards, which changed the hash but not the tree) changed the hashed set, so the provenance had to be re-established on the fixed build. The first refit ran at revision `70d1166` and shipped the fit with md5 `93af272acbe1e2068b1f5168a2510920`. Neither change touches sampling, and the redone fit is identical to the first apart from its execution record: the same chain seeds, the same draws (checked once with `identical()` after dropping `infos$execution`, outside the scripts) and a byte-identical `diagnostics.csv`. Everything below describes the second refit.

## Why the example was refitted

Per-chain seeding changes the numbers of every seeded fit. Each chain now draws its own R and C++ seeds from the session's random number stream, so a fit made with the same seed before and after the change is not the same fit. The 7 October example was made under the old scheme, so it no longer matches what the package returns for the quickstart call. The old file is replaced rather than kept beside the new one.

The refit reuses the quickstart's call, the seed `20261007` and the unchanged `data/sampledata.rda` (md5 `f98e109b8c55be17b44c32e48a2d88cc`). No seed was chosen to obtain a better-looking result.

## Settings

- Model: two-stage detection, 100 sites, 10 species, two latent species factors, occupancy covariates `X_psi.EnvCov.1` and `X_psi.EnvCov.2`, collection covariates `X_theta.1` and `X_theta.2`, spatial covariates `Xs.1` and `Xs.2`, threshold 1, posterior means for latent presences.
- No prior override was supplied, so the default priors apply.
- MCMC: two chains, 5,000 burn-in iterations and 5,000 retained iterations per chain, thinning 1.
- Random number generator: Mersenne-Twister, Inversion and Rejection, seed `20261007`.
- Software: R 4.6.1 on aarch64-apple-darwin23, package version 0.1.1, installed from the frozen export at `dev/simstudy/results/parallel-chains-20261013-4WOVWS/library`. The full session is in [session-info.txt](session-info.txt).

## The two fits

The script fitted the example twice from the same seed, once with `ncores = 1` and once with `ncores = 2`, and stopped unless the two fits were identical apart from the execution record. They were: the gate passed (`identical()` on the fits with `infos$execution` removed, and on the session random number state after each fit), and the script's assertions on both execution records held: the core budget, the workers, the threads per worker, the thread setting, the BLAS declaration, the progress file, the requested BLAS thread variables and the worker process ids.

- One core: one worker in the parent process, thread setting "ambient", so RcppParallel used its default of 10 threads on this machine. Elapsed 159.06 seconds in all (the two chains took 78.5 and 77.8 seconds).
- Two cores: two workers of one thread each, thread setting "budget". Elapsed 86.64 seconds in all (83.5 and 83.3 seconds per chain, run concurrently).
- The two-core object is the one shipped. Its worker process ids were 22161 and 22162.
- The first refit took 159.95 and 91.58 seconds.

The identity holds with the reference BLAS declared below. It is not a claim about any other BLAS.

## Chain seeds

The seeds are recorded in `infos$chain_seeds` of the shipped fit and in the provenance. Row `r` is the R seed and row `cpp` the C++ seed of each chain. They are identical to the first refit's, as they must be: the seed and the order of the draws are unchanged.

- Chain 1: `r` 1972798908, `cpp` 2045380901.
- Chain 2: `r` 1724957942, `cpp` 252297060.

## BLAS evidence and declaration

The facts come from [blas-evidence.rds](blas-evidence.rds), collected afresh for this refit by the benchmark's `blas-provenance.R` at 22:40 UTC on 13 October 2026, after the fresh export was installed, and copied unchanged into the evidence directory and here: its md5 is `04cd9a698e8ee17ea236cdc5a591640a` in all three places. Apart from its timestamp it is identical to the benchmark's evidence file (collected at 20:41 UTC), which the first refit reused. Both execution records in the refit's provenance name the same BLAS and LAPACK libraries, and the BLAS file's md5 is still the recorded one at the time of writing.

- The BLAS R loads is `/Library/Frameworks/R.framework/Versions/4.6/Resources/lib/libRblas.0.dylib`, and LAPACK is `libRlapack.dylib` in the same folder (LAPACK 3.12.1).
- It is a regular file, not a symlink, 209,392 bytes, with md5 `9dfbe56d4ceb13138673fdaf297b4b56`.
- It links only to libgfortran, libquadmath and libSystem, names no other BLAS, and defines `dgemm_` itself.
- The installer receipt `org.R-project.R.fw.pkg` (version 4.6.1) lists `libRblas.0.dylib`.
- `OCCJSDM_REFERENCE_BLAS` was set to 1 in the launching environment. Both execution records have `blas_declared` TRUE. The basis text is recorded in [blas-declaration-basis.txt](blas-declaration-basis.txt) and quoted here.

```
R 4.6.1 installed by the CRAN installer package org.R-project.R.fw.pkg, whose receipt lists libRblas.0.dylib; the file is regular, not a symlink, links only to libgfortran, libquadmath and libSystem, defines dgemm_ itself, md5 9dfbe56d4ceb13138673fdaf297b4b56; declared by Doug Yu and Claude on 2026-10-13
```

- The declaration was made by Claude, as the approved implementation plan prescribes: the plan fixes the basis text, which names Doug Yu and Claude, and has the declaration made once every fact matches the facts recorded in the design spec. Every fact matched, so the declaration was made for this run on the same basis as the benchmark's. The scripts themselves make no verdict.
- The workers were asked to run their BLAS single-threaded: the two-core execution record lists `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `MKL_NUM_THREADS` and `VECLIB_MAXIMUM_THREADS` as 1 and `MKL_DOMAIN_NUM_THREADS` as `MKL_DOMAIN_ALL=1`. These are requests the package cannot verify. The one-core fit requests nothing of its own process. No threaded BLAS was available, so nothing was measured about BLAS threads.

## Warnings seen

Each fit raised the same seven fitting warnings, in the same order (provenance field `warnings`, one list per fit), as in the first refit:

- Convergence issue in `theta0_output`, Max R-hat above 1.1.
- Convergence issue in `Bs_output`, Max R-hat above 1.1.
- Low ESS in `Bs`, Min ESS below 50.
- Convergence issue in `Gs_output`, Max R-hat above 1.1.
- Convergence issue in `As_output`, Max R-hat above 1.1.
- Low ESS in `As`, Min ESS below 50.
- Convergence issue in `Cs_output`, Max R-hat above 1.1.

The `Bs`, `Gs`, `As` and `Cs` blocks are the spatial blocks, which the public occupancy and detection diagnostic table does not include. The `theta0_output` warning was not raised by the 7 October fit, whose largest `theta0` Rhat was 1.083, below the 1.1 threshold. The refit's own log, kept in the run's evidence directory, repeats these warnings with their block summaries but is not copied here.

## Convergence limits

This is a workflow example, not a fully converged scientific analysis. Of the 130 rows returned by `returnConvergenceDiagnostics()`, 31 exceed Rhat 1.01 or fall below ESS 400 (the 7 October fit had 35 of 130). The counts below were recomputed from this refit's file and equal the first refit's, as expected. The rows are in [diagnostics.csv](diagnostics.csv) and the flagged ones in [diagnostics-flagged.csv](diagnostics-flagged.csv).

- 23 rows have Rhat above 1.01 and 24 have ESS below 400; 16 have both.
- By parameter, the flagged rows are 17 of 20 `beta_psi`, 7 of 10 `beta0_psi`, 4 of 30 `beta_theta` and 3 of 10 `theta0`. None of the 30 `p` and 30 `q` rows is flagged: their Rhats are at most 1.0018 and their ESS values at least 2,126.
- The largest Rhat is 1.100938, for `OTU_3`'s `theta0`. The smallest ESS is 90.99, for the occupancy slope of `X_psi.EnvCov.2` for `OTU_7`.

The existing two-chain schedule is kept so the saved result matches the teaching call. The quickstart's disclosure that flagged estimates need longer chains and reassessment before scientific interpretation applies to this fit too. No convergence claim is inferred from the identity gate or from the checks below.

## Hashes and revision

- The 7 October fit, which the first refit replaced: `01ed6e3da2347ddd99c2643ffc285180`, 52,102,192 bytes.
- Previous fit (the first refit of 13 October): `93af272acbe1e2068b1f5168a2510920`, 51,924,360 bytes.
- New fit: `d78e9ab78cab9c1fc5473f1f1e533caf`, 51,924,368 bytes. It differs from the previous fit only in its execution record (worker process ids, progress file path and timings).
- The refit ran at repository revision `ee3c18f` (tree `c42f3825cf2efd90c8b7ab288c6b99ca273a7159`), the branch HEAD at the time. The provenance records that commit under its hash at the time, `cac97a78c3865b51330403630c492e5827064cdd`: the commit message trailers of four commits (`efa40f9`, `ee3c18f`, `e5e00a9` and `67e3acf`) were amended after the refit and the final check, which changed their hashes but not their trees, so the recorded tree is what identifies the source. The production files it hashes (`R/`, `src/`, `DESCRIPTION`, `NAMESPACE`, 25 files) last changed in `efa40f9` (`a02d49a` before the trailer amendment), the fingerprint and error-path fix; against the first refit's provenance only `R/chains.R` and `R/runOccJSDM.R` differ. The script confirmed the hashes were unchanged after the fits.

## Relation to the 7 October refit

The 7 October verify script (`dev/release/2026-10-07-beta/verify-release.R`) now fails by design. It checks the 7 October fit's hash and the 7 October production source hashes, and both the fit and 14 of the 21 hashed source files have changed since. Running it against this checkout stops at its fit-hash assertion. The 7 October folder is a historical record and is left untouched.

## Files and verification

- [refit-sampleresults.R](refit-sampleresults.R) reproduces the refit and [refit-README.md](refit-README.md) gives the launch recipe. The script replaces `data/sampleresults.rda` in the checkout it is given.
- [verify-release.R](verify-release.R) checks the shipped data against [provenance.rds](provenance.rds): the seed, the MCMC parameters, the identity gate, the worker counts, the requested BLAS thread variables, the fit and input hashes, the production source hashes, and the shipped object's model, size, chain seeds and dimensions. Run with the repository path alone, it prints that provenance, source hashes, seed, settings and data identity are verified, as it did for this commit.
- The check against an installed copy (the second argument) was not run against the frozen library, because that library holds the data from before the refit. It was run from a fresh export instead, and is recorded in the next bullet.
- Final check, 13 October 2026: revision `e5e00a9` (hash `26b373999e406ffa454ef66ee55d9ae85b17dc2f` before the trailer amendment; tree `7f460fcc1aba7febf06f2a4c5e5a319336429cec`, identical to the checked export) was exported with `git archive`, built, and checked with `R CMD check --as-cran` under `OCCJSDM_REFERENCE_BLAS=1`, `_R_CHECK_LIMIT_CORES_=TRUE` and `NOT_CRAN=false`, giving `Status: 2 WARNINGs, 4 NOTEs` and no errors: the warnings are the Rd `\usage` mismatch for `predictNewSites` and the GNU make extensions in `src/Makevars`, as on the main branch. The verify script, given the installed export as its second argument, printed its success line, and the pkgload suite had 0 failures, 1 skip and 1383 passes. Another session's simulation study started at 23:45 BST, after the refit had finished, and kept four cores busy throughout the check and the suite, so their timings are not comparable with the first final check's. The full log is [final-check.log](final-check.log).

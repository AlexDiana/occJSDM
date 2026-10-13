# Parallel-chains example refit, 13 October 2026

This folder records the refit of the shipped example `data/sampleresults.rda` under per-chain seeding. It changes the shipped example and its evidence, not model code, priors or the input data. The refit is a workflow example whose convergence flags are disclosed. It is not a certified or fully converged fit, and successful scripts do not establish scientific convergence.

## Why the example was refitted

Per-chain seeding changes the numbers of every seeded fit. Each chain now draws its own R and C++ seeds from the session's random number stream, so a fit made with the same seed before and after the change is not the same fit. The 7 October example was made under the old scheme, so it no longer matches what the package returns for the quickstart call. The old file is replaced rather than kept beside the new one.

The refit reuses the quickstart's call, the seed `20261007` and the unchanged `data/sampledata.rda` (md5 `f98e109b8c55be17b44c32e48a2d88cc`). No seed was chosen to obtain a better-looking result.

## Settings

- Model: two-stage detection, 100 sites, 10 species, two latent species factors, occupancy covariates `X_psi.EnvCov.1` and `X_psi.EnvCov.2`, collection covariates `X_theta.1` and `X_theta.2`, spatial covariates `Xs.1` and `Xs.2`, threshold 1, posterior means for latent presences.
- No prior override was supplied, so the default priors apply.
- MCMC: two chains, 5,000 burn-in iterations and 5,000 retained iterations per chain, thinning 1.
- Random number generator: Mersenne-Twister, Inversion and Rejection, seed `20261007`.
- Software: R 4.6.1 on aarch64-apple-darwin23, package version 0.1.1, installed from the frozen export at `dev/simstudy/results/parallel-chains-20261013-wL6Toy/library`. The full session is in [session-info.txt](session-info.txt).

## The two fits

The script fitted the example twice from the same seed, once with `ncores = 1` and once with `ncores = 2`, and stopped unless the two fits were identical apart from the execution record. They were: the gate passed (`identical()` on the fits with `infos$execution` removed, and on the session random number state after each fit), and the script's own assertions on the execution records held.

- One core: one worker in the parent process, thread setting "ambient", so RcppParallel used its default of 10 threads on this machine. Elapsed 159.95 seconds in all (the two chains took 79.2 and 78.0 seconds).
- Two cores: two workers of one thread each, thread setting "budget". Elapsed 91.58 seconds in all (88.2 seconds per chain, run concurrently).
- The two-core object is the one shipped. Its worker process ids were 5907 and 5908.

The identity holds with the reference BLAS declared below. It is not a claim about any other BLAS.

## Chain seeds

The seeds are recorded in `infos$chain_seeds` of the shipped fit and in the provenance. Row `r` is the R seed and row `cpp` the C++ seed of each chain.

- Chain 1: `r` 1972798908, `cpp` 2045380901.
- Chain 2: `r` 1724957942, `cpp` 252297060.

## BLAS evidence and declaration

The facts come from [blas-evidence.rds](blas-evidence.rds), which is the benchmark's evidence file (collected at 20:41 UTC on 13 October 2026) copied unchanged: its md5 is `cb726c3269a177729ad3889d5f30ec17` in both places. The refit did not collect the evidence again. Both execution records in the refit's provenance name the same BLAS and LAPACK libraries, and the BLAS file's md5 is still the recorded one at the time of writing.

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

Each fit raised the same seven fitting warnings, in the same order (provenance field `warnings`, one list per fit):

- Convergence issue in `theta0_output`, Max R-hat above 1.1.
- Convergence issue in `Bs_output`, Max R-hat above 1.1.
- Low ESS in `Bs`, Min ESS below 50.
- Convergence issue in `Gs_output`, Max R-hat above 1.1.
- Convergence issue in `As_output`, Max R-hat above 1.1.
- Low ESS in `As`, Min ESS below 50.
- Convergence issue in `Cs_output`, Max R-hat above 1.1.

The `Bs`, `Gs`, `As` and `Cs` blocks are the spatial blocks, which the public occupancy and detection diagnostic table does not include. The `theta0_output` warning was not raised by the 7 October fit, whose largest `theta0` Rhat was 1.083, below the 1.1 threshold. The refit's own log, kept in the run's evidence directory, repeats these warnings with their block summaries but is not copied here.

## Convergence limits

This is a workflow example, not a fully converged scientific analysis. Of the 130 rows returned by `returnConvergenceDiagnostics()`, 31 exceed Rhat 1.01 or fall below ESS 400 (the 7 October fit had 35 of 130). The rows are in [diagnostics.csv](diagnostics.csv) and the flagged ones in [diagnostics-flagged.csv](diagnostics-flagged.csv).

- 23 rows have Rhat above 1.01 and 24 have ESS below 400; 16 have both.
- By parameter, the flagged rows are 17 of 20 `beta_psi`, 7 of 10 `beta0_psi`, 4 of 30 `beta_theta` and 3 of 10 `theta0`. None of the 30 `p` and 30 `q` rows is flagged: their Rhats are at most 1.0018 and their ESS values at least 2,126.
- The largest Rhat is 1.100938, for `OTU_3`'s `theta0`. The smallest ESS is 90.99, for the occupancy slope of `X_psi.EnvCov.2` for `OTU_7`.

The existing two-chain schedule is kept so the saved result matches the teaching call. The quickstart's disclosure that flagged estimates need longer chains and reassessment before scientific interpretation applies to this fit too. No convergence claim is inferred from the identity gate or from the checks below.

## Hashes and revision

- Previous fit (7 October refit): `01ed6e3da2347ddd99c2643ffc285180`, 52,102,192 bytes.
- New fit: `93af272acbe1e2068b1f5168a2510920`, 51,924,360 bytes.
- The refit ran at repository revision `70d116603ff01112b442f49d49d007b27b1022f2`, which is the branch HEAD at the time and is recorded in the provenance. The production files it hashes (`R/`, `src/`, `DESCRIPTION`, `NAMESPACE`, 25 files) are unchanged since `f9fc489`, and the later commits touch only benchmark and release scripts and their records. The script confirmed the hashes were unchanged after the fits.

## Relation to the 7 October refit

The 7 October verify script (`dev/release/2026-10-07-beta/verify-release.R`) now fails by design. It checks the 7 October fit's hash and the 7 October production source hashes, and both the fit and 14 of the 21 hashed source files have changed since. Running it against this checkout stops at its fit-hash assertion. The 7 October folder is a historical record and is left untouched.

## Files and verification

- [refit-sampleresults.R](refit-sampleresults.R) reproduces the refit and [refit-README.md](refit-README.md) gives the launch recipe. The script replaces `data/sampleresults.rda` in the checkout it is given.
- [verify-release.R](verify-release.R) checks the shipped data against [provenance.rds](provenance.rds): the seed, the MCMC parameters, the identity gate, the worker counts, the requested BLAS thread variables, the fit and input hashes, the production source hashes, and the shipped object's model, size, chain seeds and dimensions. Run with the repository path alone, it prints that provenance, source hashes, seed, settings and data identity are verified, as it did for this commit.
- The check against an installed copy (the second argument) is not run here, because the frozen library holds the data from before the refit. It runs against a fresh export in Task 9.
- Final check, 13 October 2026: revision `7904cde634f0a8e97f5ecd8a003acb111ab89840` was exported with `git archive`, built, and checked with `R CMD check --as-cran` under `OCCJSDM_REFERENCE_BLAS=1`, `_R_CHECK_LIMIT_CORES_=TRUE` and `NOT_CRAN=false`, giving `Status: 2 WARNINGs, 4 NOTEs` and no errors: the warnings are the Rd `\usage` mismatch for `predictNewSites` and the GNU make extensions in `src/Makevars`, as on the main branch. The verify script, given the installed export as its second argument, printed its success line, and the pkgload suite had 0 failures, 1 skip and 1373 passes. The full log is [final-check.log](final-check.log).

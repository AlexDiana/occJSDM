# The end-of-fit convergence summary: cost now grows with the number of parameters, not its square

Measured on 16 November 2026 on branch `claude/linear-diagnostics`, against `main` at `48b2b3a` (the code released as 0.3.0). The question was why `runOccJSDM()` spends minutes and gigabytes after the sampler has finished when a fit has thousands of coefficients, and whether that can be removed without changing what the summary reports.

The one-sentence answer: `computeDiagnostics()` handed each parameter block whole to `coda::gelman.diag()`, which builds several matrices with one row and one column per parameter even when asked for univariate statistics only, and replacing that call with a calculation that treats each parameter on its own gives the same numbers in a few seconds and a few megabytes.

## What was wrong

`runOccJSDM()` ends by calling `computeDiagnostics()` whenever there are two or more chains. For each parameter block (`B`, `B0`, `beta_theta`, `p`, `q` and so on) it built a `coda::mcmc.list` holding every parameter of the block and called `coda::gelman.diag()`, with `autoburnin = FALSE` and `multivariate = FALSE`, for R-hat and `coda::effectiveSize()` for the effective sample size.

`multivariate = FALSE` only stops `gelman.diag()` computing the multivariate statistic. It still forms the full covariance matrix of the parameters within every chain (`var(x)` on a draws-by-parameters matrix), averages those matrices cell by cell with `apply(S2, c(1, 2), mean)`, and forms two more parameters-by-parameters matrices for the covariance terms of its degrees-of-freedom correction. It then reads only the diagonal of each. With P parameters, n draws and m chains that is m times n times P squared arithmetic for the covariances, and an R function call and a small allocation for each of the P squared cells of the average, which from reading the code is probably where most of the memory goes (the old code was measured as a whole, not profiled).

`coda::effectiveSize()` was never quadratic: it already loops over parameters. It was slow for a different reason, a call to `lm()` and a call to `ar()` for every parameter of every chain.

The block that matters is the occupancy coefficients `B`, which has covariates times species parameters. A fit with 250 embedding columns and 50 species has 12,500.

## What changed

`computeDiagnostics()` now takes both statistics from a new internal function, `blockRhatESS()`, in `R/diagnostics.R`. Nothing in it is shared between parameters.

- `chain_moments_ess_cpp()` in the new `src/diagnostics.cpp` reads the block of draws in place, with no copy, and returns the mean, the sample variance and the effective sample size of every chain of every parameter. The effective sample size is `coda::effectiveSize()`'s, reimplemented: the test for a chain that is a straight line (which is how coda gives a chain that never moved zero effective draws), then a Yule-Walker autoregression with its order chosen by AIC up to the default maximum of `stats::ar()`, then the chain's variance over the fitted spectral density at frequency zero.
- `gelmanRhat()` computes the R-hat point estimate from the chain means and variances. It follows the arithmetic of `gelman.diag()` line for line and keeps its variable names, but computes only the diagonal entries that `gelman.diag()` went on to use.
- For a chain whose autoregression is numerically singular (a non-finite or tied AIC), the C++ code declines to answer and `blockRhatESS()` passes that one parameter to `coda::effectiveSize()` itself. No chain in any test or benchmark triggered this, apart from draws of magnitude 1e160, where coda stops with an error and so, still, does the new code.

The printed summary, the two warnings and the rules for skipping a block or dropping a fixed parameter are untouched. So are `returnConvergenceDiagnostics()` and the helpers behind it, the sampler, and everything a fit stores.

## What was measured

- One synthetic block of P parameters with 2 chains of 5,000 draws, P = 300, 3,000, 6,000 and 12,500, as the only parameter set of a minimal `results_output`. Each parameter and chain is an AR(1) series with coefficient 0.5. The same seed gives the same draws to both versions, so their printed summaries can be compared.
- "Before" is `main` at `48b2b3a`. "After" is this branch at `d644362`. Each was exported with `git archive` and installed into its own library with `R CMD INSTALL`, so both are optimised, byte-compiled builds, as a user has.
- Each size ran in a fresh R process (`bench-diagnostics.R`, driven by `run-bench.sh`). Time is the elapsed time of the one call to `computeDiagnostics()`.
- Memory is reported two ways. "R heap above input" is the most R heap in use at any moment of the call, less what was in use just before it: what the function allocated on top of the draws it was given. "Process peak" is the peak memory footprint of the whole R process from `/usr/bin/time -l`, which includes R itself, the loaded packages (about 150 MB) and the block of draws (8 bytes per draw: 23, 229, 458 and 954 MB at the four sizes).
- Conditions: Apple M4 laptop (four performance and six efficiency cores, 32 GB), macOS 27.2.1, on mains power with Low Power Mode off (`pmset -g` reported `lowpowermode 0` at the start of every run), R 4.6.1 with the Accelerate BLAS, coda 0.19.4.1, Apple clang 21 at R's default `-O2`. The calculation is single-threaded in both versions and uses no BLAS. The laptop was not idle: see the note on variation under Results. The old code's run at 12,500 parameters also overlapped with a compile and two runs of the test suite, which does not affect its outcome.
- The old code could not finish 12,500 parameters on this machine. That run was given a limit of 16 GB of memory footprint and 40 minutes, and `run-bench.sh` stopped it when it passed the first. (`footprint` reports whole gigabytes at that size, so the stop came at the first reading of 17 GB.)

## Results

Old code, `main` at `48b2b3a`:

- 300 parameters: 2.6 s. R heap above input 226 MB. Process peak 0.48 GB.
- 3,000 parameters: 125 s. R heap above input 3.0 GB. Process peak 3.6 GB.
- 6,000 parameters: 479 s. R heap above input 8.4 GB. Process peak 10.0 GB.
- 12,500 parameters: did not finish. Stopped after 23 minutes (1,395 s) when its footprint passed the 16 GB limit; the last reading was 17 GB and still rising. The square law from 6,000 parameters puts the whole run at about 35 minutes and about 36 GB of R heap, more than the laptop has.

New code, this branch at `d644362`:

- 300 parameters: 0.11 s. R heap above input 3 MB. Process peak 0.15 GB.
- 3,000 parameters: 0.97 s. R heap above input 6 MB. Process peak 0.38 GB.
- 6,000 parameters: 2.0 s. R heap above input 8 MB. Process peak 0.63 GB.
- 12,500 parameters: 4.3 s. R heap above input 13 MB. Process peak 1.19 GB.

Reading them:

- Time. The old code takes 48 times longer for 10 times the parameters (300 to 3,000) and 3.8 times longer for twice (3,000 to 6,000): quadratic. The new code takes 8.8 and 2.1 times longer over the same steps, and 2.1 times longer again for the 2.1-fold step to 12,500: linear, at about 0.3 ms per parameter. It is 24 times faster at 300 parameters, 130 times at 3,000 and 240 times at 6,000.
- Memory the function adds. The old code needed 13 to 18 times the size of its input in extra R heap at 3,000 and 6,000 parameters, growing with the square of the block. The new code needs 13 MB at 12,500 parameters, about 3 MB fixed plus under 1 KB per parameter, for the per-chain matrices and the vectors of the R-hat arithmetic; its C++ work buffer is a fixed 1.3 MB.
- Memory of the whole process. For the new code the process peak equals the footprint the process already had before the call, to within 1 MB at every size (152, 380, 629 and 1,190 MB before; 152, 381, 629 and 1,190 MB at peak). The 1.19 GB at 12,500 parameters is the block of draws itself (954 MB), R with the package loaded (about 150 MB) and some garbage left from building the block. So the process total is above 1 GB only because the input alone is 0.95 GB; the diagnostics add nothing measurable to it.
- Resident set size understates the old code's cost. At 6,000 parameters `/usr/bin/time` reported a maximum resident set size of 5.2 GB beside the peak footprint of 10.0 GB, because macOS compresses idle pages. Footprint is what Activity Monitor shows and what runs a laptop out of memory, so it is the figure used here.

Run-to-run variation is large, and the final pass was the slow one. The laptop's load average was about 9 during it (another R session and a backup were running). A first pass on a quieter machine, with an earlier revision of the scripts, gave 1.3, 77 and 423 s for the old code at 300, 3,000 and 6,000 parameters, and 0.08, 0.83, 1.7 and 3.5 s for the new code at the four sizes (`evidence/first-pass.txt`). The R heap figures barely moved between passes (225, 2,719 and 8,571 MB for the old code; 3, 6, 8 and 13 MB for the new). Neither pass changes any conclusion.

For scale: on 16 November a reviewer measured a real fit with 133 covariates and 48 species (6,384 occupancy coefficients, 2 chains of 10 sweeps) at 265 s and about 6 GB against 16.5 s for one chain, and stopped a 261-covariate fit after 14 minutes at about 16 GB. Those figures are of the same order as the 6,000-parameter run here, although that fit had 10 draws per chain and this block 5,000, which suggests the part of the cost that does not depend on the number of draws (the cell-by-cell average over P squared cells) dominates. That was not profiled.

## The values are the same

- On the benchmark blocks the two versions print the same summary to every digit shown (R-hat to three decimals, ESS to one) at the three sizes the old code finished: see `evidence/before.txt` and `evidence/after.txt`.
- On a real fit (`compare-fit.R`: the quickstart two-stage spatial call on `sampledata`, 2 chains, 100 burn-in and 200 kept iterations, fixed seed) the two builds store an identical `results_output` and print an identical summary with the same 25 warnings (`evidence/compare-fit.txt`).
- `tests/testthat/test-diagnostics.R` compares `blockRhatESS()` with coda directly on small blocks: autocorrelation from -0.9 to 0.98, one parameter, chains of 10 to 5,000 draws, one chain stuck, every chain stuck at a different value (R-hat `Inf`, ESS 0), chains that are copies of each other (R-hat `NaN`), a chain that is a straight line, a constant parameter, a parameter with an infinite draw, and one chain. It also runs `computeDiagnostics()` beside the previous implementation, kept verbatim in `tests/testthat/helper-legacy-diagnostics.R`, and requires identical console output on a synthetic `results_output` that exercises every branch of the reporting loop, on chains of two and three draws, and on the two-stage and continuous fixture fits.
- A looser check outside the suite, 400 random blocks with 3 to 1,000 draws, 2 to 5 chains, spreads from 1e-4 to 1e4 and offsets up to 1e5, agreed with coda to a median relative difference of 2e-16 in R-hat and 2e-15 in ESS, and 2e-8 at worst in both. The worst cases are blocks whose draws sit at an offset 1e7 to 1e8 times their spread, where subtracting the mean costs that many digits in either implementation; no other block differed by more than 1e-10.
- The full suite passes against the installed build of `d644362`: 281 tests, 1,854 expectations, no failures, and the 6 skips that `main` has on this laptop (5 need a declared single-threaded reference BLAS, 1 is the coverage study behind `OCCJSDM_SIMSTUDY`).

Two places where "the same" needs a qualifier:

- The values agree to rounding, not bit for bit. The summary prints three decimals of R-hat and one of ESS, so a digit can differ only if a value falls within rounding error (typically 1e-15, relatively) of a rounding boundary. The same holds for the two thresholds behind the warnings (R-hat above 1.1, ESS below 50) and for the threshold that marks a parameter as fixed (pooled variance at most 1e-12).
- With a single kept draw per chain and two or more chains, the old code stopped with an error from `stats::ar()` at the end of the fit, and the new code does too, by design: this change does not alter behaviour, even where the behaviour looks like a defect. See the follow-ups.

## Follow-ups, not done here

- `returnConvergenceDiagnostics()` is linear in the number of parameters but slow: it calls coda once per parameter and builds a one-row tibble for each. `summarisePosterior()` took 2.7 s for 300 parameters with 2 chains of 5,000 draws, about 9 ms per parameter, so a 12,500-parameter block would take about 2 minutes and a 500,000-parameter one over an hour. It could use `blockRhatESS()`, but its R-hat follows a slightly different rule for chains that never moved (`NA` where `computeDiagnostics()` keeps coda's value), so that is a change of its own.
- A fit that keeps one draw per chain on two or more chains fails in `computeDiagnostics()` after the sampler has finished (`'order.max' must be >= 1`). Nobody runs such a fit on purpose, but it loses the fit.
- The new code is single-threaded. At the measured rate, about 0.3 ms per parameter for 2 chains of 5,000 draws, a block of 500,000 parameters (2,000 species by 250 covariates) would take about 2.5 to 3 minutes. It parallelises trivially over parameters if that ever matters; the draws themselves, 40 GB at that size, will matter first.
- Branch `claude/grouped-covariate-prior` (draft PR 41) and the study on `claude/alphaeco-geospatial-embeddings-452a98` were cut before this change. Both need it merged or cherry-picked (commits `0d81436` and `d644362`) before their wide fits are run.

## Reproducing

``` sh
git archive 48b2b3a | tar -x -C <src-before>
R CMD INSTALL --no-docs --library=<lib-before> <src-before>
git archive d644362 | tar -x -C <src-after>
R CMD INSTALL --no-docs --library=<lib-after> <src-after>

sh run-bench.sh before <lib-before> evidence/before.txt "300 3000 6000" 20 16
sh run-bench.sh before <lib-before> evidence/before.txt "12500" 40 16
sh run-bench.sh after <lib-after> evidence/after.txt "300 3000 6000 12500" 10 8

Rscript --vanilla compare-fit.R fit <lib-before> before.rds
Rscript --vanilla compare-fit.R fit <lib-after> after.rds
Rscript --vanilla compare-fit.R compare before.rds after.rds
```

Do not run the old code at 12,500 parameters without the limits: it needs more memory than the laptop has.

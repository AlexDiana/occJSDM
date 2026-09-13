# RNG safety validation, 11 September 2026

The beta implementation runs every RNG-bearing sampling body on R's main thread. It retains spatial fields, traits, latent factors, splines, collection covariates and multiple primers. Deterministic TBB workers remain available. This item is ready for Alex's review; it does not resolve the separate spatial, correlation or point-estimate bias gates.

The baseline is collection-alignment commit `90026d6`. Implementation commit `9d95045` on branch `codex/serial-sampling-rng` changes four RNG worker invocations and the C++ stream storage. The subsequent commit adds this report, reproducible diagnostics, public fitting documentation and removal of the obsolete known-gap test. Priors, conditional distributions, parameter updates and the serial traversal order are unchanged. The seed mapping remains `seed_seq{base_seed, 0}`. The generation counter triggers reseeding and normal-cache reset but contributes no seed material. Chains and iterations continue the stream instead of repeatedly reseeding it.

## Baseline and regression evidence

The unmodified baseline passed its existing source suite, which skipped the optional coverage study and the known RNG gap. The new `test-rng-safety.R` failed nine assertions before implementation. With four threads requested, a constant 256 by 32 PG predictor matrix produced only 2,447 distinct draws out of 8,192; repeated blocks visibly replayed the same stream. PG outputs changed across one/four thread requests, and collection coefficient draws differed from the existing serial implementation. Whole binary, occupancy and two-stage fits changed across thread requests and repeated four-thread fits. Continuous fits already reproduced in this check. This confirms that the live PG path matters independently of collection sampling.

After implementation, all 29 new assertions pass. They compare complete `results_output` objects and final R RNG states across requested one/four threads and repeated four-thread fits for all four data models. The fixtures include spatial fields, traits, latent factors, collection covariates and two primers where applicable. The collection wrapper also exactly matches the existing serial sampler, including the final R RNG state. A repeated PG call advances the C++ stream. Existing regressions confirm that a changed seed and consecutive fits without resetting the seed produce different draws. The package preserves the requested RcppParallel setting instead of changing it to one.

The final source suite passes 396 assertions with zero failures, errors or test warnings. Its only skip is the opt-in coverage study. `R CMD INSTALL --install-tests` succeeds, including loading from temporary and final library locations. The installed-package suite passes 360 assertions with zero failures, errors or test warnings and eight expected skips: the coverage study, five CRAN-gated checks, and two checks that require the source tree. Source tests exercise the latter seven checks. The environment emits a pre-existing warning that testthat was built under R 4.5.2, while this run uses R 4.5.0. Roxygen regeneration reports two existing empty `@details` tags in `R/output.R`; only the RNG and fitting help pages change.

## Complete live RNG path audit

- Spatial support-point preparation calls `stats::kmeans()` from `computeSpatialSummaries()` on the R thread; its random centre selection precedes sampling. `runOccJSDM()` obtains the C++ base seed with R's `sample.int()` on that same thread. The chain and iteration loops are ordinary sequential R loops. The model initialisation and posterior reparameterisation sections introduce no worker sampling.
- The occupancy/collection `z` update calls `sample_z_cpp_parallel()`. `ZProbWorker` computes probabilities using arithmetic and logarithms. After `parallelFor()` returns, a main-thread loop makes the `R::rbinom()` draws.
- The two-stage `w` update calls `sample_w_cim_cipp_parallel()`. `ProbWorker` computes conditional probabilities using arithmetic and logarithms. Its `R::rbinom()` draws also occur in the main-thread loop after `parallelFor()` returns.
- The collection coefficient update calls `sample_betatheta_cpp_parallel()`. `BetaThetaWorker` now executes directly as `worker(0, S)`. Its `sample_beta_nocov_cpp_TS()` leaf reaches the C++ PG generator and `sample_beta_cpp()`'s Armadillo normal generator. Both now execute serially. The historical `_TS` suffix confers no permission to call R RNG from a worker.
- The JSDM binary likelihood augmentation calls `samplePGvariables_parallel()`. `PG_Worker` now executes directly in the existing column-major order. Its PG rejection sampler uses the single advancing C++ uniform/normal stream in `rng.h`.
- The joint coefficient update remains the live serial `sample_BBsL_cpp()`. Its `sampleB_SoR()` normal draws use the serial C++ stream. The serial `sample_GC()` and `sample_A()` updates call `sampleBuniv()` for observed/latent trait coefficients. `sample_U_cpp()` iterates serially and draws through `sampleB()`. Spatial counterparts use the same paths. `sample_sigmab()`, `sample_sigmah()` and continuous-model `sample_tau()` call `rinvgamma_cpp()`, whose `R::rgamma()` remains on the R thread. `sample_ls()` uses ordinary R `runif()` calls on that thread.
- Lab detection rates call `sample_pq_cpp_parallel()`. `BetaParamWorker` only counts outcomes and computes Beta shape parameters. The subsequent `R::rbeta()` loop is serial. Collection false-positive rates use `sample_theta0()` and R's `rbeta()` in a sequential species loop.

There are exactly three executable `parallelFor()` sites left in package C++ source: `ZProbWorker`, `ProbWorker` and `BetaParamWorker`. None calls a random generator or another sampling function. No executable OpenMP sampling pragmas remain. A source review also checked compatibility/dormant paths: `SampleOmegaWorker` and `BBSL_Worker` now execute directly on the main thread, so neither can reintroduce unsafe scheduling if its helper is used. The old serial collection helpers and serial PG entry point stay serial. The unused `_TS_opt` leaf still uses R-backed Armadillo randomness and must remain on the main thread; it has no callers. No dormant helper is newly wired into the public fitter.

## Runtime thread audit and independent distribution checks

`audit_rng_threads.py` creates a separate minimal package copy and instruments all actual R/Armadillo RNG statements in both C++ files plus the shared C++ uniform/normal leaves. A guard captures the thread running `setOccJSDMSeed()` and rejects any instrumented RNG call on a different thread before drawing or updating the audit counters. The production checkout contains none of this instrumentation.

The instrumented package passed 18 short fits: requested one thread, four threads, and a repeated four-thread fit for each of binary, continuous, occupancy, two-stage, spline two-stage and intercept-only zero-factor two-stage configurations. Every pair of complete results and final R RNG states matched exactly. Each fit used two chains, five burn-in iterations and eight retained draws per chain. Convergence/ESS warnings are expected for these deliberately short path checks and are saved in `fit-warnings.txt`; these runs make no convergence or parameter-recovery claim.

Including the distribution checks below, the guard recorded 699,888 C++ uniform calls, 60,455 C++ normal calls, 56,636 Armadillo normal calls, 124,800 R binomial calls, 3,744 R Beta calls and 1,638 R gamma calls, all on the main thread. Ordinary R-level RNG calls are covered by the sequential call-path audit above.

Independent moment checks use 32,768 draws per distribution, with four threads requested. PG(1,0) has known mean 0.25 and variance 1/24. A zero-design Gaussian update contributes no likelihood precision, so the collection and C++ coefficient draws must each follow the supplied Normal prior with mean 0.7 and variance 1.6. This makes the reference independent of fitted posterior output and of the sampler's own summaries.

- PG(1,0): mean 0.249351, variance 0.0416207, adjacent-pair correlation 0.000158, zero duplicate draws. The mean error is -0.576 Monte Carlo standard errors.
- Collection/R normal: mean 0.699088, variance 1.612793, adjacent-pair correlation 0.001100, zero duplicate draws. The mean error is -0.130 Monte Carlo standard errors.
- C++ normal: mean 0.697552, variance 1.589726, adjacent-pair correlation -0.000773, zero duplicate draws. The mean error is -0.350 Monte Carlo standard errors.

The runner checks mean error below five analytic Monte Carlo standard errors, relative variance error below 6%, absolute adjacent correlation below 0.04 and no duplicated draws. These are modest distribution/stream diagnostics, not a proof of statistical independence for every finite pseudorandom sequence. Together with the single advancing stream and the absence of RNG worker scheduling, they directly address the duplicated-stream mechanism that reproducibility alone could miss.

## Reproduction and limits

Run from the package root, using a fresh destination for the audit copy:

``` sh
Rscript --vanilla -e 'devtools::test(filter = "rng-safety", stop_on_failure = TRUE)'
Rscript --vanilla -e 'devtools::test(stop_on_failure = TRUE)'
python3 dev/simstudy/audit_rng_threads.py . /absolute/output/thread-audit
Rscript --vanilla dev/simstudy/validate_rng_safety.R --source=/absolute/output/thread-audit --out=/absolute/output/audit-results
Rscript --vanilla dev/simstudy/validate_rng_safety.R --out=/absolute/output/plain-results
```

The diagnostic seed is 6182 for the C++ stream, with R seed 6123; the short fits use seed 12831. Raw logs and regenerable CSV/RDS output from this review are under `/Users/douglasyu/Documents/Codex/2026-09-10/fam/work/rng-validation`, outside the tracked checkout. `moments.csv`, `fit-checks.csv`, `main-thread-counts.csv` and `session-info.txt` are produced by the tracked scripts. The tested environment is R 4.5.0 on arm64 macOS 26.6.2 with RcppParallel/TBB and requested thread counts one and four.

Serial execution preserves the existing one-thread seed mapping and conditional sampling formulas while removing dependence on TBB task scheduling. It may reduce sampling throughput. Exact draws across operating systems, C++ standard libraries or R/package versions are not promised. Linux/Windows CI and final combined-branch release checks remain appropriate. This work does not estimate posterior point bias, establish interval coverage, redesign parallel streams, or change any prior. The spatial, residual-correlation, `B0` and high-`q` checks remain separate release work.

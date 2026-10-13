# Parallel MCMC chains for runOccJSDM()

Date: 13 October 2026. Status: design agreed in conversation with Doug on 13 October 2026, provisionally; written for Doug's review, including an external review he will commission, before an implementation plan is made. Nothing has been implemented or run. This is the `TODO.md` item "Run MCMC chains in parallel" under *Needed for the paper*. It will go to Alex as a pull request against `AlexDiana/occJSDM`, in the style of PR #14.

## Purpose

`runOccJSDM()` runs its MCMC chains one after another in one process. This design runs them on separate cores, one chain per worker process at a time, with reproducible per-chain seeds, an explicit core budget, and exact combination of the posterior arrays, the latent-state means and the WAIC accumulators. The result of a fit must not depend on how many workers ran it.

Why now, decided 13 October 2026: the paper's simulation study should run on the code that ships, so that its validation also covers the chain-combination code. The production simulation run waits until Alex merges this. Parallel chains will not speed that study, which already parallelises over replicates; the gain is for users, lesson builds and small studies.

## What the code does today

Verified on 13 October 2026 against `main` at `d72daaf`. Line numbers refer to `R/runOccJSDM.R`, which no commit since `880ddb0` touches.

- Line 1051 calls `setOccJSDMSeed(sample.int(.Machine$integer.max, 1L))` once per fit, seeding the C++ Mersenne-Twister stream in `src/rng.h` from R's RNG. Then `for (chain in 1:nchain)` at line 1053 runs the chains serially. Every chain resets to identical starting values (lines 1095 to 1199), so chains differ only by their stretch of the stream.
- A chain also consumes R's own RNG, not only the C++ stream: `sample_theta0()` draws with `rbeta()` (`R/mcmcfun.R` line 24), the inverse-gamma draws go through `R::rgamma()` under an `RNGScope` (`src/jsdm.cpp` line 286), and the length-scale index uses `sample.int()` (`R/jsdmfun.R` line 1211). Per-chain seeding therefore has to set both generators at the start of every chain.
- Posterior draws are stored per chain in `*_output_chain` arrays and copied into full arrays with trailing dimensions `[..., niter, nchain]` after each chain (lines 1367 to 1402).
- With `summarisedLatentPresences = TRUE`, the default, the latent-state means of `z`, `psi`, `theta` and `w` accumulate inside the loop as `(1 / (niter * nchain))` times the value, once per retained iteration (lines 1296 to 1314). `theta` and `w` are always returned as means; only `z` and `psi` have a full-array option.
- The legacy complete-data WAIC accumulates across all chains with `update_waic_summary()` (`R/mcmcfun.R` line 35), which keeps per observation a running mean of the log-likelihoods, a running sum of squared deviations (`M2`) and a running mean of the likelihoods, updated with a counter that runs from 1 over all chains (lines 1340 to 1360). `compute_waic()` (`R/mcmcfun.R` line 56) finalises it with `niter * nchain` draws (lines 1406 to 1420). `computeSiteWAIC()` in `R/site-waic.R` works on the stored arrays and indexes draws as iteration within chain, so it is unaffected by how the chains were run.
- `reparamFactorModel()` (`R/jsdmfun.R` line 6) rotates each draw separately, so it can stay in the parent after combination and cannot depend on the worker count.
- The roxygen `@details` (lines 449 to 455) say random draws run serially on the main thread whatever the RcppParallel thread count. That is correct. `tests/testthat/test-rng-safety.R` was run against the installed 0.1.1 build on 13 October 2026 and all 29 assertions pass, including fits bit-identical at 1 and 4 threads for the binary, continuous, occupancy and two-stage models. `AGENTS.md`'s August 2026 statement that fits are bit-reproducible only at one thread (the section "Detail behind Fixed bugs 40-42", and the Fixed bugs 41 advice to set `RCPP_PARALLEL_NUM_THREADS=1`) is stale; it was superseded by Fixed bugs 50 on 14 September 2026.
- A pre-existing oddity: with `summarisedLatentPresences = FALSE`, a per-chain `theta_output_chain` array is allocated (line 1063) but never written or read; `theta_output` is always the posterior mean.
- `.onLoad()` in `R/zzz.R` sets `options(mc.cores = 2)`. Nothing in the package reads it, and this design does not use it.

## Decisions made on 13 October 2026

- The worker count lives in `MCMCparams$ncores`, default 1, so fits are serial unless asked.
- `ncores` is the total core budget: workers = min(nchain, ncores), and each worker gets floor(ncores / workers) RcppParallel threads.
- A parallel fit shows a progress file plus a summary. The start message includes the `tail` command.
- The shipped `sampleresults` is refitted in this pull request under the new code.
- Architecture: one single-chain function, driven by `lapply()` for serial fits and by `parLapply()` on a PSOCK cluster for parallel fits. Two alternatives were considered. Recursive single-chain calls of `runOccJSDM()` stitched together by the parent were rejected because the data preparation would repeat per chain, the scalar WAIC cannot be merged exactly, and the serial path would also have to recurse. Partial combination inside the workers was rejected as premature because the parent holds every chain's arrays either way.
- Doug's approvals are provisional until this spec has been through an external review.

## Non-goals

- No redesign of the random streams within a process, and no stream splitting. Independence between chains comes from distinct seeds on two Mersenne-Twister generators.
- No dispersed starting values. The mirror-labelling check in `TODO.md` remains a separate item.
- No change to the sampling order within an iteration, to any sampler, prior or default other than `ncores`.
- No fail-fast across chains beyond what base `parallel` allows.
- No change to the lessons' teaching archives, and no `TODO.md` edit in the pull request.
- No Windows benchmark. The mechanisms are cross-platform, and the report says it was not run there.
- No dependency beyond base `parallel`. `pkgload` is used only when the parent itself was loaded by it.

## Interface

- `MCMCparams$ncores`, read with `get_param()` like `nchain`, default 1. It must be a single whole number of at least 1, otherwise an error. Above `parallel::detectCores()` it warns and proceeds.
- Budget rule: workers = min(nchain, ncores). Each worker gets floor(ncores / workers) RcppParallel threads, at least 1, set in the worker through `RCPP_PARALLEL_NUM_THREADS` and `RcppParallel::setThreadOptions()`, as `dev/simstudy/run_study.R` already does for its own workers. Any remainder of the budget is idle; for example `ncores = 3` with two chains gives two single-threaded workers.
- `ncores = 1` is today's fit: in process, with the session's RcppParallel setting untouched.
- `ncores > 1` with `nchain = 1` also runs in process, with the session's RcppParallel thread count set to `ncores` for the fit and restored afterwards, so the budget means the same thing everywhere except the untouched default.
- Messages. "Running MCMC" is printed as today. For `ncores > 1` the parent then prints one start line, for example "Running 4 chains on 4 workers, 1 RcppParallel thread each. Follow progress with: tail -f /path/occJSDM-progress-1a2b3c.log". The workers' existing "Chain k - Iteration n" and burn-in messages, printed every 100 iterations, go to that file. In both modes the parent prints one line per chain with its elapsed time once the chains have finished, for example "Chain 1 finished in 312.4 s". For `ncores = 1` the chain messages print to the console as today.
- The progress file is `tempfile("occJSDM-progress-", fileext = ".log")`. It is not deleted; the session's temporary directory goes when R exits.
- The returned object's `infos` gains two fields, for provenance and for rerunning one chain; nothing in the package reads them. `MCMCparams` holds the resolved settings: `nchain`, `nburn`, `niter`, `nthin`, `ncores`, `workers`, `threads_per_worker`. `chain_seeds` is a 2 by `nchain` integer matrix with rows `r` and `cpp`.
- The roxygen promise that the calling session's RcppParallel thread setting is preserved stays true in every mode.

## Seeding and reproducibility

- Where seeds come from. After the existing setup and before any chain runs, the parent draws `2 * nchain` distinct integers with `sample.int(.Machine$integer.max, 2L * nchain)` from R's RNG, in whatever state the user left it. Chain k receives an R seed and a C++ seed. The parent then saves `.Random.seed`.
- What a chain does with them. At its start a chain calls `set.seed()` with its R seed and the kinds pinned to Mersenne-Twister, Inversion and Rejection, then `setOccJSDMSeed()` with its C++ seed. Pinning the kinds means a worker, whose fresh session has R's defaults, and the parent, whose session may be in any kind, produce the same R-level draws.
- After the chains. The serial path restores the parent's saved `.Random.seed` by assignment in the global environment, which also restores the RNG kinds. The parallel path never touches the parent's RNG after the draw. So at `ncores = 1` and at any larger value the session ends the fit in the identical state. The parent runs no random code after the chains: factor reparameterisation and diagnostics are deterministic.
- Independence by construction. A chain's randomness is fixed by its own two seeds alone, so it is the same whether it runs first, last, in process or in another process. No chain continues from where another stopped. That is the property the unexplained mirror observation in `TODO.md` ("Checks for Doug") singled out: 2 of 6 chains run one after another in one process entered a mirror labelling, against 0 of 16 run as separate processes. Whether the shared stream caused it is not known and this design does not claim to settle it; it removes the layout. Starting values stay common.
- The documented contract. `set.seed()` before `runOccJSDM()` reproduces a fit for any `ncores` on the same platform and build. Consecutive fits without resetting the seed differ, because the chain seeds come from R's RNG. A seeded fit's numbers differ from those of version 0.1.1 because every chain is now seeded on its own, including chain 1, whose R-level draws no longer continue the session's stream.
- What this changes elsewhere. Everything pinned to old numbers: the shipped `sampleresults` and its 7 October 2026 provenance, handled by the refit below; the quickstart's statement that the same seed reproduces `sampleresults`, true again after the refit, with its software line updated from R 4.5.0 to R 4.6.1; and the opt-in coverage snapshot in `tests/testthat/test-coverage-study.R`, which may need re-recording if a block moves by more than its 0.02 tolerance. No test pins fitted values to a seed; the tests assert that the same seed reproduces, that different seeds differ and that thread counts do not matter, all of which remain true. The lesson archives, the study records under `dev/simstudy` and the MCP pilot's pinned runtime are unaffected. A refit under the new code cannot be expected to match any archive made before it.

## The single-chain function

- Placement. A new file `R/chains.R` holds `run_one_chain()`, the task wrapper, the driver, the combination helpers and the cluster helpers. `R/runOccJSDM.R` loses the loop body and keeps the setup and the output assembly. The names below are indicative; the implementation plan fixes them.
- `run_one_chain(setup, seed_r, seed_cpp, chain)`. `setup` is one plain list built in `runOccJSDM()` from what the loop body reads: the model flags (`model`, `jsdmModel`, `threshold`); the data and index objects (`y`, `y_NA`, `M`, `K`, `P`, `sumM`, `sumP`, `sumK`, `maxP`, `primerId_p`, `idx_w_k`, `idx_p_k`, `idx_z_w`, `list_idx`, and the sizes `n`, `N`, `N3`, `S`, `d`, `g`, `gt`, `ps`, `ncov_theta`, `ncov_psi`); the design matrices (`X_theta`, `X_psi`, `Tr`); the precomputed spatial objects (`list_Xs`, `list_SoRSummaries`); the priors (`list_data`, `list_priors`, `a_theta0`, `b_theta0`, `a_p`, `b_p`, `a_q`, `b_q`, `b_betatheta`, `B_betatheta`, `spatial_sd_prior`); and the MCMC settings (`nburn`, `niter`, `nthin`, `summarisedLatentPresences`). The implementation plan fixes the exact list from a read of the loop. A PSOCK worker enforces it: a missing field is an error there, and the always-on tests run on workers.
- Body. The current starting-values block and iteration loop, moved, with four changes only. The seeds are set at the top. The latent-state means become per-chain sums of the retained `z`, `psi`, `theta` and `w`. The WAIC accumulators become per-chain, with the iteration counter running from 1 within the chain. The dead `theta_output_chain` allocation is dropped, since it would otherwise be `N * S * niter` doubles of `NA` per chain crossing the socket for nothing. Sampling order within an iteration and the message wording are unchanged.
- Output. A list with `draws`, the per-chain arrays as today; `latent_sums`, the four sums, or `NULL` where the model has none; `waic`, the three partial accumulators, `NULL` where the model has none, with their draw count; `elapsed`, in seconds; `chain`; and the two seeds.

## Workers

- Cluster. `parallel::makeCluster(workers, outfile = progress_file)`, PSOCK, with `parallel` added to `Imports` in `DESCRIPTION`. Immediately after creation an `on.exit()` handler registers the shutdown described under failures. No `clusterSetRNGStream()`: the chain seeds are explicit.
- Bootstrap, run once per worker through `clusterCall()`. The worker sets `.libPaths()` to the parent's, because a PSOCK worker inherits environment variables but not library paths changed at run time. It pins its thread count from the budget rule, through the environment variable and `setThreadOptions()`. Then it loads the same occJSDM as the parent. If the parent's namespace was loaded by pkgload, detected from the marker pkgload leaves in the namespace, the worker runs `pkgload::load_all()` on the parent's source path, which `getNamespaceInfo("occJSDM", "path")` gives; the shared object is already compiled by the parent, so this only loads it. Otherwise the worker runs `loadNamespace("occJSDM")`. `pkgload` is not an import: that branch runs only when the parent was loaded by pkgload, which means pkgload is installed on the same library paths.
- Build identity check. Each worker reports its occJSDM namespace path, its version and the md5 of its loaded shared object, taken from `getLoadedDLLs()`. The parent compares all three with its own and stops, naming both paths, on any mismatch. The md5 is the real test of "same build"; the path catches library mix-ups.
- What crosses the socket. The task function is defined at package level in `R/chains.R`, so it serialises as a reference to the namespace, and takes `setup` and the seed matrix as arguments passed through the dots of `parLapply()`. It must not be a closure created inside `runOccJSDM()`: serialising such a closure would ship the whole frame of `runOccJSDM()`, including the preallocated full output arrays. `parLapply()` sends the arguments once per chunk and splits the chains into consecutive chunks, one per worker, so a worker running two chains receives `setup` once. The serial path calls the same task function through `lapply()`.
- Sequence in the parent. Build `setup`; draw the seeds; save the RNG state; if workers > 1, create the cluster, bootstrap, check builds, run `parLapply()`; otherwise run `lapply()` and restore the RNG state; combine; finish as today.

## Combination, in the parent, in chain order

- Draws. Each chain's arrays are copied into the preallocated `[..., niter, nchain]` arrays at index k, as the loop does today. Which worker produced chain k cannot matter.
- Latent means. The per-chain sums are added in chain order and divided by `niter * nchain`. Same quantity as today, a different floating-point path, independent of `ncores`.
- WAIC. For each of the three accumulators the chains are merged pairwise in chain order with the standard formula for combining running means and sums of squared deviations. For parts A and B with counts nA and nB and n = nA + nB: delta is meanB minus meanA; the merged mean is meanA + delta * nB / n; the merged M2 is M2A + M2B + delta squared * nA * nB / n; the likelihood means combine as a count-weighted average. Then `compute_waic()` as today with `niter * nchain` draws. The current check that the counter equals `niter * nchain` becomes a check that every chain reports `niter` draws. The merged value differs from today's single pass in rounding only.
- Afterwards, unchanged: factor reparameterisation on the full arrays, the output lists, `computeDiagnostics()`, `infos`.

## Failures, interrupts and memory

- Error in a chain. The task wraps `run_one_chain()` in `tryCatch()` and returns a structured failure holding the chain number and the message, instead of letting the generic "one node produced an error" of `parallel` surface. Base `parallel` delivers results only once every chunk has returned, so the parent waits for the running chains and then stops with "Chain k failed: message". The worker's own output, including any traceback, is in the progress file.
- Shutdown on every exit, normal, error or interrupt. The `on.exit()` handler calls `stopCluster()` and then `tools::pskill()` on the worker PIDs recorded at bootstrap. The second step matters: a worker that is mid-chain reads the stop message only after its chain finishes, so without it, pressing Escape in the parent would leave the workers computing to the end of their chains. `pskill()` on a worker that has already exited is a harmless no-op.
- Memory. Each worker holds `setup` plus one chain's arrays. The parent holds `setup`, the full arrays and, while assembling, the list of every chain's returned draws, so its peak is briefly about twice the draws, against today's full arrays plus one chain. Each returned element is dropped as soon as it is copied. The benchmark measures the parent and the workers separately, so that the report states the overhead rather than guesses it.

## Documentation changes

- Help page of `runOccJSDM()`: the `MCMCparams` entry gains `ncores` with the budget rule; `@details` is rewritten for the seeding contract, the progress file and the change from 0.1.1; `@return` notes the two new `infos` fields. `man/` is regenerated with roxygen2 8.0.0.
- Help page of `sampleresults` in `R/data.R`: the refit settings, including `ncores` and R 4.6.1.
- `NEWS.md`, development version: the feature in one bullet and the behaviour change in its own bullet, because every seeded fit changes.
- Quickstart `vignettes/occJSDM-00-quickstart.Rmd`: one sentence on `ncores` at the `MCMCparams` bullet; the software line in the "same seed" sentence updated to R 4.6.1; every sentence that describes the shipped fit, such as which blocks it warned about, re-checked against the refitted object and corrected. The rendered `.md` is rebuilt.
- `AGENTS.md`: a dated note marking the "Detail behind Fixed bugs 40-42" reproducibility paragraph, and the Fixed bugs 41 advice to pin one thread, as superseded by Fixed bugs 50. Nothing is rewritten or removed.
- `TODO.md` is not touched in the pull request. After Alex merges, in its own commit on `main`: the parallel item moves to completed work, as the WAIC item did; the mirror item gets a line saying chains are now seeded independently, so the planned repeat runs on the new code; and the "Fully reproducible parallel random draws within a process" item no longer waits on parallel chains.

## The sampleresults refit

- When: after implementation and both test runs, and after the benchmark has confirmed identity across worker counts at production size.
- How: the quickstart call exactly, with `set.seed(20261007)`, two chains of 5,000 burn-in and 5,000 retained iterations, thinning 1, two factors, both covariate sets, the spatial covariates, threshold 1, summarised latent presences, default priors. Fitted with an installation made from a `git archive` export of the branch head into a dedicated library, under R 4.6.1. Fitted twice, at `ncores = 1` and `ncores = 2`, and the two objects must be identical before either is saved; the `ncores = 2` object is what ships.
- Provenance: a new dated folder under `dev/release/` with the refit script, `provenance.rds` holding the source hashes, the seed, the chain seeds, `ncores`, the input and output md5s, a verify script and session info, following the 7 October 2026 layout. The 7 October folder stays as history, and its verify script is expected to fail its source-hash check from now on.
- Size: the new `.rda` is about 52 MB, so the repository grows by that much.

## Tests

A new `tests/testthat/test-parallel-chains.R`, under the tier-1 rules in `helper-fixtures.R`: structural assertions, a 30-second budget for the always-on tier, and at most two workers, because that is what `R CMD check` allows.

- Identity across worker counts, always on: for the two-stage spatial fixture and the binary fixture, `nchain = 2`, the same `set.seed()`: `results_output` and the session's `.Random.seed` identical at `ncores = 1` and `ncores = 2`. Each cluster start costs a few seconds, so the continuous and occupancy models, and an `nchain = 3` case at `ncores = 2` where one worker runs two chains, run under `skip_on_cran()` and in the benchmark.
- Seeds: the same seed reproduces at `ncores = 2`; a different seed differs. Extends the Fixed bugs 28 test in `test-regression-bugs.R`.
- Combination without fitting: the WAIC merge against a single pass of `update_waic_summary()` over concatenated draws, to a relative tolerance, including the final WAIC value. The latent means, by checking under one seed that the mean over all draws of an unsummarised fit's `z` equals the summarised fit's posterior mean, since the option changes storage only.
- Interface: `ncores` of 0, 1.5, a string or a vector errors; above `detectCores()` warns.
- Workers, through the internal driver without a fit: a task that errors produces "Chain k failed: ..." and no worker PID survives; workers report the pinned thread count from the budget rule; the build-identity comparison is unit tested with fabricated inputs, since two builds cannot be made inside a test.
- Progress file: one parallel fit with 100 burn-in and 100 retained iterations on the small fixture, so the every-100 messages fire; the file named in the start message exists and contains the chain lines.
- Both loading branches: `devtools::test()` under pkgload exercises the `load_all()` branch; `R CMD check` from the export exercises `loadNamespace()`. The final evidence includes both runs, with the tier-1 timing recorded.
- The existing suite must pass unchanged, apart from the opt-in coverage snapshot noted above.

## Benchmark

In `dev/simstudy/parallel-chains/`, following the repository's study layout.

- Fit: the quickstart call on `sampledata`, four chains of 5,000 and 5,000, at `ncores = 1` with the session's default threads; `ncores = 1` with one thread pinned, to show what the threaded deterministic parts contribute in process; then 2, 4 and 8. Each configuration twice. Roughly 90 minutes on this laptop (Apple M4, 10 cores, 32 GB), R 4.6.1, from an installation made from a `git archive` export of the branch head. Launched detached with `nohup` and `setsid` from a non-interactive shell, on mains power with the lid open.
- Measured: wall time; per-chain elapsed from the summary; peak resident memory of the parent and of each worker, sampled every second with `ps`; and identity of `results_output` and of `.Random.seed` against the `ncores = 1` fit for every configuration.
- Reported: `REPORT.md` with bullet lists and a CSV of raw timings. The pull request description quotes the speed-up and the memory overhead. Windows is untested and the report says so.

## Process and deliverables

- Worktree `.worktrees/parallel-chains`, branch `claude/parallel-chains` from `main`. This spec is the branch's first commit. All builds and installs come from a `git archive` export into a scratch library; never `R CMD build` a checkout.
- Order: this spec approved after Doug's review; an implementation plan with the writing-plans skill; implementation with tests first; the full suite under pkgload and `R CMD check` from the export, against the recorded baseline of 0 errors, 2 warnings and 2 notes; the benchmark; the refit; the documentation; the pull request against `AlexDiana/occJSDM` `main` with an explanation for Alex in the PR #14 style: what changed, the seeding contract, why every seeded fit's numbers change, the benchmark headline, and how to review the moved loop.
- After merge: the `TODO.md` commit on `main` described above.

## Risks and open points

- The tier-1 time budget: cluster starts in the always-on tests are measured during implementation, and more cases move under `skip_on_cran()` if the 30 seconds are exceeded.
- `pkgload::load_all()` in several workers at once recompiles only if the sources changed after the parent loaded them, which does not happen in a normal session; if it ever did, the workers could race on the compilation. Not guarded.
- The parent's peak memory at assembly, about twice the draws, is the price of returning whole chains through `parLapply()`. If the benchmark shows it matters at realistic sizes, chains could be collected one at a time with the lower-level functions of `parallel` in a later change.
- Windows: `makeCluster()`, `outfile` and `tools::pskill()` are documented to work there; not run here.
- The WAIC merge changes the floating-point path of a value that also changes with the seeds; its correctness is covered by the unit test against the single pass.

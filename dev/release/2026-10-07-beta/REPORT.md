# Beta release refresh, 7 October 2026

The records below describe the initial beta published at `63a4faf`. The [subsequent release refresh](refresh-REPORT.md) includes the pushed MCP and teaching work, its new downloads and fresh package checks; the original fit and verification records are retained.

The beta release refresh uses merged main `a56f5548e68e09d244513d09fcc04ad52f0caa19`, including PR #14. It changes the shipped example and release documentation, not model code, priors or historical scientific evidence. Source tests, the full source-package build, the installed-package check and both installed beta vignettes pass. Independent review found no critical or important issues and judged the change ready to merge. Successful software checks do not establish scientific convergence or remove the limitations below.

## Refitted example

`data/sampleresults.rda` was regenerated from the unchanged `data/sampledata.rda` using a fresh private installation of the merged package. The model matches the quickstart: two-stage detection, 100 sites, 10 species, two latent species factors, both occupancy and collection covariates, species traits, the spatial covariates, threshold 1 and posterior means for latent presences. Two chains each use 5,000 burn-in iterations and 5,000 retained iterations, with thinning 1. The seed is `20261007`, using Mersenne-Twister, Inversion and Rejection.

No prior override was supplied. The model retains the default Normal(0, 1) occupancy baseline, Beta(5, 1) laboratory true-positive rate, Beta(1, 20) laboratory and field false-positive rates, collection-slope prior variance 2 and default inverse-gamma spatial-amplitude prior. The fitted object retains the threshold, trait transformation and applied prior metadata.

The refit requested one RcppParallel thread and set `OMP_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS` and `MKL_NUM_THREADS` to 1. This machine's R Makeconf has an empty `SHLIB_OPENMP_CXXFLAGS`, and the fresh library has no undefined OpenMP runtime symbols. These are macOS checks, not evidence about Linux OpenMP performance. The fit took 164.673 seconds before saving and compression.

The new fit has MD5 `01ed6e3da2347ddd99c2643ffc285180` and is 52,102,192 bytes, compared with the previous 54,063,081 bytes. [provenance.rds](provenance.rds) records the production source hashes, input identity, old/new fit hashes, exact call, RNG, thread settings, times, package versions and fitting warnings. [refit-sampleresults.R](refit-sampleresults.R) reproduces the fit, [verify-release.R](verify-release.R) verifies its identity and installed-package copy, and [session-info.txt](session-info.txt) records the runtime.

## Convergence limits

This is a workflow example, not a fully converged scientific analysis. Of the 130 occupancy/detection rows returned by `returnConvergenceDiagnostics()`, 35 exceed Rhat 1.01 or fall below ESS 400. The largest Rhat is 1.083009 for `OTU_3`'s `theta0`; the smallest ESS is 137.6736 for an occupancy slope. Lab rates mix better in this run: all `p` and `q` Rhats are below 1.01 and their ESS values exceed 2,300.

The internal fitting diagnostics also warn about `Bs_output`, `Gs_output`, `As_output` and `Cs_output`, with low ESS in `Bs` and `As`. Those spatial blocks are not included in the public occupancy/detection diagnostic table. The complete warnings and block summaries are retained in [refit.log](refit.log); all 130 tabular diagnostics are in [diagnostics.csv](diagnostics.csv), with flagged rows in [diagnostics-flagged.csv](diagnostics-flagged.csv).

The quickstart discloses the flags, retains its advice to leave spatial covariates out of real beta fits for now, and explains that flagged estimates require longer chains and reassessment before scientific interpretation. The existing two-chain schedule is kept so the saved result matches the teaching call. No seed was selected to obtain a better-looking result, and no convergence claim is inferred from successful package tests.

## Documentation and verification

The source-test baseline passed with one deliberately opt-in coverage study skipped. The refreshed fit's production hashes, input identity, seed, default priors, dimensions, threshold and finite occupancy summaries pass the release verification script.

The quickstart and simulator guide both render independently to HTML and unwrapped GitHub Markdown from the refreshed installation. The quickstart's diagnostics and two figures were regenerated and visually inspected. The simulator's Markdown and three figures are newly tracked, giving the beta simulator guide the same readable review format as the quickstart; its R Markdown source is unchanged. All five figures were visually inspected. Startup notices report installed dependencies built under R 4.5.2 while this host runs R 4.5.0; these notices are outside the rendered teaching output and do not indicate new occJSDM test warnings.

The final source suite passes 1,050 expectations with zero failures or test warnings and one deliberately opt-in coverage study skipped. The clean source archive builds both vignettes. `R CMD check --no-manual` completes with zero errors, three existing warnings and two existing notes, including successful rebuilding of vignette outputs. Installed-package tests pass 1,014 expectations with zero failures/test warnings and eight standard skips. [source-tests.log](source-tests.log), [installed-tests.log](installed-tests.log), [package-build.log](package-build.log) and [package-check.log](package-check.log) retain the exact results.

The three warnings are the R header's unsupported clang warning option `-Wfixed-enum-extension`, the existing undocumented `predictNewSites(verbose)` argument, and GNU Makefile extensions in `src/Makevars`. The two notes concern `LICENSE` metadata and existing R-code analysis findings, including global-variable/function references and partial argument matching. These are the same diagnostic categories recorded before the refresh. Installed size is 51.9 MB, reported as information rather than a check failure. This is not a CRAN-readiness claim.

[verify-installed.log](verify-installed.log) confirms that the checked installation's `sampleresults` is identical to the saved file, its public diagnostics match the CSV, and the only installed vignettes are `occJSDM` and `simulateOccJSDMData`. The lesson-link regression checks pass and all 10 site/publication tests pass; their logs are retained beside this report. Production files, `sampledata`, historical teaching bundles, lesson publication flags and exclusions match the merged base.

The checked source archive is `/private/tmp/occjsdm-beta-release-20261007/occJSDM_0.1.0.tar.gz`, with SHA-256 `2ee848e7c09c8f656456da4488e35abb5eaec55f6adf13f1ec9fb9c9b3a2ef5f`. It contains the refreshed example and the two built beta vignettes, excluding Lessons 0-7, teaching bundles, development evidence and Git/workflow bookkeeping.

## Independent review

A fresh reviewer inspected commit `fbaa703` against merged main, independently ran the installed-data verification, confirmed the archive hash and compared its data and vignette sources with the worktree. The reviewer found no critical or important issues and judged the convergence disclosures adequate for a workflow demonstration. The reviewer read the logged suites without rerunning them and visually checked the quickstart figures and simulator trait figure.

One minor reproducibility improvement is deferred: the refit script sets RcppParallel threads but only records the four BLAS/OpenMP environment variables, so a future launch must use the recorded environment to reproduce this artifact's conditions. The checked artifact does record all four as 1. This does not undermine its verification.

Review scope was settled as follows: retain and disclose the existing demonstration schedule instead of certifying scientific convergence; retain the documented prior/model limitations rather than re-investigate unchanged sampler and calibration studies; restrict software claims to the checked macOS environment rather than CRAN acceptance or cross-platform numerical equivalence; treat integration, publication and announcement delivery as subsequent actions under the user's release instruction; preserve the primary checkout by isolation, with its same two modified MCP documents still present; and retain trailing spaces in raw logs and generated output to preserve emitted evidence. Authored source and prose pass the whitespace check when these emitted artifacts are excluded.

## Publication scope

The beta includes only the quickstart and simulator guide. Lessons 0-7 and all historical teaching bundles retain their exclusions and recorded evidence. CRAN submission, pkgdown deployment, broader calibration studies and lesson publication remain deferred.

[release-notes.md](release-notes.md) is the proposed GitHub prerelease description for `v0.1.0-beta`; [announcement.md](announcement.md) reconciles the existing announcement with the README and withheld lessons. Neither document's existence means a release was published or an announcement was sent.

Integration is tracked in [PR #30](https://github.com/AlexDiana/occJSDM/pull/30). The [GitHub beta release](https://github.com/AlexDiana/occJSDM/releases/tag/v0.1.0-beta) records publication status and the release revision. Announcement recipients have not been supplied.

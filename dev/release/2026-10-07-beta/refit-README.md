# Launching the beta example refit

Run this from the root of an isolated occJSDM checkout. The script replaces that checkout's `data/sampleresults.rda`, so use a worktree rather than your ordinary working copy. It needs an installed copy of the package from the same source checkout and a fresh evidence directory.

To reproduce the original beta fit, use the scientific source revision and R/package versions recorded in [REPORT.md](REPORT.md), [provenance.rds](provenance.rds) and [session-info.txt](session-info.txt). The recorded R version is 4.5.0. Running the command against newer code is a new refit, whose provenance records that revision; it does not recreate the historical artifact by definition.

Select the intended R installation before starting. The `R` and `Rscript` commands below must use that same installation. The four thread variables must be set in the launching shell, before R starts; `refit-sampleresults.R` itself sets the RcppParallel thread count to one and records these environment variables.

```sh
refit_repo="$(pwd -P)"
mkdir -p "$refit_repo/dev/simstudy/results"
refit_run="$(mktemp -d "$refit_repo/dev/simstudy/results/beta-refit-XXXXXX")"
refit_library="$refit_run/library"
refit_evidence="$refit_run/evidence"
mkdir -p "$refit_library" "$refit_evidence"

R CMD INSTALL --library="$refit_library" "$refit_repo"

env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 MKL_NUM_THREADS=1 \
    Rscript --vanilla \
    "$refit_repo/dev/release/2026-10-07-beta/refit-sampleresults.R" \
    "$refit_repo" "$refit_library" "$refit_evidence" \
    > "$refit_run/refit.log" 2>&1
```

The evidence directory receives the seed, source/input/fit hashes, fitting settings, requested thread counts, timings, warnings, diagnostics and session information. Check `provenance.rds` after the run: all four environment variables should be `"1"`. These requested settings do not establish that OpenMP is enabled, nor that the example has converged; the original beta report describes both limitations.

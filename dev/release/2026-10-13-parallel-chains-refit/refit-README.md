# Launching the parallel-chains example refit

Run this from the root of the isolated worktree for the parallel-chains branch. The script replaces `data/sampleresults.rda` in the checkout it is given, so it must be run in the worktree and never in your ordinary working copy. It needs an installed copy of the package built from the same frozen source, a fresh evidence directory, and a BLAS evidence file collected for that installation by the benchmark's `blas-provenance.R` script.

The refit fits the quickstart example twice under per-chain seeding, once with one core and once with two, requires the two fits to be identical apart from the execution record, and ships the two-core object. Per-chain seeding changes the numbers of every seeded fit, which is why the example had to be refitted.

The `OCCJSDM_REFERENCE_BLAS` variable must be set in the launching shell, before R starts. The script refuses to run without it and refuses an evidence directory that already holds a `provenance.rds`, so make a fresh directory for every run. `OUT` is the benchmark output directory, which holds the frozen `library`, `blas-evidence.rds` and `blas-declaration-basis.txt`, or a directory prepared the same way from a later export, as for the second refit described in `REPORT.md`.

```sh
cd /Users/douglasyu/src/occJSDM/.worktrees/parallel-chains
EV=$(mktemp -d "/Users/douglasyu/src/occJSDM/dev/simstudy/results/parallel-chains-refit-$(date +%Y%m%d)-XXXXXX")
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' env OCCJSDM_REFERENCE_BLAS=1 /usr/local/bin/Rscript --vanilla dev/release/2026-10-13-parallel-chains-refit/refit-sampleresults.R "$(pwd -P)" "$OUT/library" "$EV" "$OUT/blas-evidence.rds" > "$EV/refit.log" 2>&1 < /dev/null &
cp "$OUT/blas-declaration-basis.txt" "$EV/"
```

The script takes four arguments, in this order: `REPO` (the checkout whose `data/sampleresults.rda` is replaced), `LIB` (the frozen library), `EVIDENCE` (the fresh evidence directory) and `BLAS_EVIDENCE_RDS` (the benchmark's `blas-evidence.rds`, copied into the evidence directory as `blas-evidence.rds`).

Wait for `REFIT COMPLETE` in `$EV/refit.log`. Two fits of the quickstart size take about five minutes. Then copy the evidence into this folder and run the verify script against the checkout.

```sh
cp "$EV"/provenance.rds "$EV"/session-info.txt "$EV"/diagnostics*.csv "$EV"/blas-evidence.rds "$EV"/blas-declaration-basis.txt dev/release/2026-10-13-parallel-chains-refit/
Rscript dev/release/2026-10-13-parallel-chains-refit/verify-release.R "$(pwd -P)"
```

The verify script should print `Refit provenance, source hashes, seed, settings and data identity verified.` It also accepts an installed library as a second argument and then checks that the installed `sampleresults` is identical to the shipped file. The frozen library holds the pre-refit data, so that form of the check is meant for a fresh export built after the refit.

The evidence directory receives `provenance.rds`, `session-info.txt`, `diagnostics.csv`, `diagnostics-flagged.csv` and `blas-evidence.rds`. The provenance records the seed, the chain seeds, the source, input and fit hashes, the fitting settings, the execution records of both fits, the BLAS evidence file names, the timings and the warnings. A flagged row in `diagnostics-flagged.csv` has an Rhat above 1.01 or an effective sample size below 400, or a non-finite value of either.

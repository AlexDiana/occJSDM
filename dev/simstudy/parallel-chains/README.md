# Parallel-chains benchmark

This folder holds the scripts that time and measure the quickstart fit at several core budgets (`MCMCparams$ncores`), check that the fit and the random-number stream are identical across worker counts, and record the provenance of the BLAS the run used.

## What a run does

- Freezes the build: `launch.sh` exports `HEAD` with `git archive`, installs the export into a dedicated library under the run directory, and records the commit in `revision.txt`. The benchmark and the later refit both run from that one library, never from the working tree.
- Collects BLAS evidence with `blas-provenance.R` (location, symlink status, linkage, whether the file defines `dgemm_` itself, md5, and the macOS installer receipt). The script makes no verdict.
- Runs `benchmark.R`: four chains at five core budgets (`serial-ambient`, `serial-1thread`, `ncores-2`, `ncores-4`, `ncores-8`), two repeats each, with `sample-rss.sh` sampling resident memory of the parent and of every PSOCK worker once a second.

## The BLAS declaration

The declaration is the operator's decision, made after reading the evidence. For the 13 October 2026 runs the implementation plan prescribed it, and Claude made it once every recorded fact matched. Run the provenance script on its own first, compare the output with what is known about the installation, and only then launch with `OCCJSDM_REFERENCE_BLAS` and `OCCJSDM_REFERENCE_BLAS_BASIS` set in the environment. `launch.sh` refuses to start without both, and `benchmark.R` writes the basis to `blas-declaration-basis.txt`. The identity claims in the report hold under that declaration only.

## How to launch

Use a fresh run directory every time, from a non-interactive shell, on mains power with the lid open, from the worktree root. The run directory goes in the main checkout's results folder, given as an absolute path below, because `dev/simstudy/results/` is git-ignored and does not exist inside a worktree: a relative path there makes `mktemp` fail and leaves `OUT` empty.

```sh
OUT=$(mktemp -d "/Users/douglasyu/src/occJSDM/dev/simstudy/results/parallel-chains-$(date +%Y%m%d)-XXXXXX")
Rscript dev/simstudy/parallel-chains/blas-provenance.R "$OUT/blas-evidence-preview.rds"
OCCJSDM_REFERENCE_BLAS=1 OCCJSDM_REFERENCE_BLAS_BASIS="..." \
nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' \
  sh dev/simstudy/parallel-chains/launch.sh "$OUT" \
  > "$OUT/launch.log" 2>&1 < /dev/null &
```

The full run takes about 35 minutes (36 minutes on 13 October 2026). `launch.sh` refuses to reuse a run directory that already holds a build, results or a reference.

To run only some configurations, name their labels after the run directory, for example `sh dev/simstudy/parallel-chains/launch.sh "$OUT" serial-ambient`; with none named all five run, and an unknown label stops the run before any fit. The configurations run in their fixed order and the first one to run becomes the identity reference that every later fit is compared with, so a subset that omits `serial-ambient` is compared with its own first fit (which may be a worker fit) rather than with the in-process fit; `benchmark.R` does not warn about this, so include `serial-ambient` unless that is intended.

## Where results land

Everything goes under the run directory, which is in the main checkout's git-ignored `dev/simstudy/results/` folder:

- `launcher.pid`, `launch.log`, `install.log`, `benchmark.log`: the launch and its progress. The benchmark log ends with `ALL DONE`.
- `revision.txt`, `src-export/`, `library/`, `scripts/`: the frozen build and a copy of these scripts.
- `blas-evidence.rds`, `blas-evidence.log`, `blas-declaration-basis.txt`: the BLAS evidence and the declaration.
- `results.csv`, `rss/`, `execution/`: timings and identity checks, per-configuration memory samples, and each fit's execution record.
- `reference-fit.rds`, `reference-rng.rds`: the first fit and its random-number state, which every later fit is compared against.
- `session-info.txt`: written when the benchmark finishes.

## After the run

`Rscript dev/simstudy/parallel-chains/summarise.R "$OUT"` prints the bullet lists that `REPORT.md` is written from. The 13 October 2026 results are kept in this folder: `REPORT.md` (what was run, the BLAS evidence and declaration, the results, the rerun and the limits), `results.csv` (the full run), `rerun-results.csv` (the later rerun of the two ambient serial fits, made with the label selector), `blas-evidence.rds` and `blas-declaration-basis.txt`. The run directories themselves stay in the main checkout's ignored `dev/simstudy/results/`.

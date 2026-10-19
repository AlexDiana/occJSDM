#!/bin/sh
# Freeze the build from a git archive export, collect BLAS evidence, run the
# benchmark detached. The BLAS declaration is the operator's, made in the
# launching environment after reading the evidence from a first run of
# blas-provenance.R; this script never sets it. Run from the worktree root,
# from a non-interactive shell:
#   OCCJSDM_REFERENCE_BLAS=1 OCCJSDM_REFERENCE_BLAS_BASIS="..." \
#   nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' \
#     sh dev/simstudy/parallel-chains/launch.sh OUTDIR [LABEL ...] \
#     > OUTDIR/launch.log 2>&1 < /dev/null &
# Labels after OUTDIR are passed to benchmark.R, which then runs only those
# configurations.
set -eu
OUT="$1"
shift
REPO="$(pwd -P)"
RSCRIPT=/usr/local/bin/Rscript
: "${OCCJSDM_REFERENCE_BLAS:?declare the BLAS yourself after reading the evidence}"
: "${OCCJSDM_REFERENCE_BLAS_BASIS:?state the basis of the declaration}"
# A run directory is used once. Refuse anything that already holds a build,
# results or a reference, so an earlier run's provenance is never overwritten
# and a new build is never compared against an old reference.
for stale in "$OUT/revision.txt" "$OUT/results.csv" "$OUT/reference-fit.rds" "$OUT/library/occJSDM" "$OUT/rss"; do
  if [ -e "$stale" ]; then echo "refusing to reuse $OUT: $stale exists; use a fresh run directory" >&2; exit 1; fi
done
mkdir -p "$OUT/src-export" "$OUT/library" "$OUT/scripts" "$OUT/rss" "$OUT/execution"
cp "$REPO"/dev/simstudy/parallel-chains/*.R "$REPO"/dev/simstudy/parallel-chains/*.sh "$OUT/scripts/"
git -C "$REPO" rev-parse HEAD > "$OUT/revision.txt"
git -C "$REPO" archive HEAD | tar -x -C "$OUT/src-export"
R CMD INSTALL --library="$OUT/library" "$OUT/src-export" > "$OUT/install.log" 2>&1
"$RSCRIPT" "$OUT/scripts/blas-provenance.R" "$OUT/blas-evidence.rds" > "$OUT/blas-evidence.log" 2>&1
caffeinate -i -w $$ &
echo $$ > "$OUT/launcher.pid"
"$RSCRIPT" "$OUT/scripts/benchmark.R" "$OUT" "$@" > "$OUT/benchmark.log" 2>&1
echo "ALL DONE $(date '+%F %T')"

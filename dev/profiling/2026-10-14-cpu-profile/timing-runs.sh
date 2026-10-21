#!/bin/sh
# Baseline and prototype timings, run one after another so they do not compete.
# Usage: sh timing-runs.sh <run-dir>   (run from the checkout root; OCCJSDM_PKG may name another checkout)
set -e
OUT=$1; HERE=$(cd "$(dirname "$0")" && pwd); export OCCJSDM_PKG=${OCCJSDM_PKG:-$(pwd)}
mkdir -p "$OUT"; cp "$HERE/proto-logweights.R" "$OUT/"
[ -f "$OUT/binary900.rds" ] || Rscript "$HERE/make-binary900.R" "$OUT/binary900.rds"
NOPROF=1 Rscript "$HERE/profile-fit.R" quickstart 1000 1000 "$OUT/quickstart-noprof" > "$OUT/quickstart-noprof.log" 2>&1
NOPROF=1 Rscript "$HERE/profile-fit.R" binary900 8 8 "$OUT/binary900-noprof" "$OUT/binary900.rds" > "$OUT/binary900-noprof.log" 2>&1
PROTO=1 NOPROF=1 Rscript "$HERE/profile-fit.R" quickstart 1000 1000 "$OUT/quickstart-proto" > "$OUT/quickstart-proto.log" 2>&1
PROTO=1 NOPROF=1 Rscript "$HERE/profile-fit.R" binary900 8 8 "$OUT/binary900-proto" "$OUT/binary900.rds" > "$OUT/binary900-proto.log" 2>&1
echo ALLDONE

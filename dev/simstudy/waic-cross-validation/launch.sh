#!/bin/sh
# Run every fit, then every score, then the summary, detached from any tool
# session (memory: fits spawned under an agent session die with it).
# Launch from a non-interactive shell:
#   nohup perl -MPOSIX -e 'POSIX::setsid() or die; exec @ARGV or die' \
#     sh OUTDIR/scripts/launch.sh OUTDIR RESULTS_DIR \
#     > OUTDIR/logs/launch.log 2>&1 < /dev/null &
# Runs the copies in OUTDIR/scripts, so editing the repository mid-run is safe.
set -eu
OUT="$1"
RESULTS="$2"
SCRIPTS="$OUT/scripts"
RSCRIPT=/usr/local/bin/Rscript
WORKERS="${CV_WORKERS:-8}"
export R_LIBS="$OUT/library"
export RCPP_PARALLEL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
echo $$ > "$OUT/logs/launcher.pid"
# Prevent idle sleep while this launcher runs (lid-close sleep still pauses it).
caffeinate -i -w $$ &

stage() {
  name="$1"; script="$2"
  echo "START $name $(date '+%F %T')"
  xargs -P "$WORKERS" -I{} sh -c \
    '"$0" "$1/scripts/$2" "$1" "$3" > "$1/logs/$4-$3.log" 2>&1 && echo "DONE $4 $3 $(date +%T)" || echo "FAILED $4 $3 $(date +%T)"' \
    "$RSCRIPT" "$OUT" "$script" {} "$name" < "$OUT/jobs.txt"
  missing=0
  while read -r job; do
    if [ ! -f "$OUT/$name/$job.rds" ]; then echo "MISSING $name $job"; missing=1; fi
  done < "$OUT/jobs.txt"
  if [ "$missing" -ne 0 ]; then echo "STOPPED after $name $(date '+%F %T')"; exit 1; fi
  echo "END $name $(date '+%F %T')"
}

stage fits fit-job.R
stage scores score-job.R
"$RSCRIPT" "$SCRIPTS/summarise.R" "$OUT" "$RESULTS" > "$OUT/logs/summarise.log" 2>&1
echo "ALL DONE $(date '+%F %T')"

#!/bin/sh
# One step of the run, named by the first argument; see README.md.
#   draws: replicates 1 to 10 of the four M = 3 arms at 1,000 + 1,000,
#          keeping draws (the both-means run's refit.R).
#   main:  the four M = 3 arms at R = 100 through the study runner.
#   long:  replicates 1 to 5 of the four M = 3 arms at 5,000 + 5,000,
#          keeping draws.
# 5 worker processes and one thread per fit in every step; the steps run one
# at a time, so the laptop never has more than 5 fitting workers. Launched
# detached from a non-interactive shell. Every output goes to a path that
# .gitignore lists, so each step records a clean tree.
set -u
STEP=$1
WT=/Users/douglasyu/src/occJSDM/.claude/worktrees/relaxed-herschel-025463
RUN=$WT/dev/simstudy/ordinary-m3-20261016
REFIT=$WT/dev/simstudy/theta-community-mean-20261016/refit.R
RSCRIPT=/usr/local/bin/Rscript
ALL=base_m3,occupancy_m3,base_cm_m3,occupancy_cm_m3
cd "$RUN"
case "$STEP" in
  draws) OUT=$RUN/draws/short; SETTINGS="replicates 1 to 10 of $ALL, nburn = 1000, niter = 1000, nchain = 2" ;;
  main)  OUT=$RUN/results; SETTINGS="cells $ALL, R = 100 per cell, nburn = 1000, niter = 1000, nchain = 2" ;;
  long)  OUT=$RUN/draws/long; SETTINGS="replicates 1 to 5 of $ALL, nburn = 5000, niter = 5000, nchain = 2" ;;
  *) echo "unknown step: $STEP" >&2; exit 2 ;;
esac
if [ -e "$OUT" ]; then echo "refusing: $OUT exists" >&2; exit 1; fi
{
  echo "step: $STEP"
  echo "launched: $(date '+%F %T')"
  echo "revision: $(git -C "$WT" rev-parse HEAD)"
  echo "branch: $(git -C "$WT" rev-parse --abbrev-ref HEAD)"
  echo "dirty files: $(git -C "$WT" status --porcelain | wc -l | tr -d ' ')"
  echo "R: $("$RSCRIPT" -e 'cat(R.version.string)')"
  echo "BLAS: $("$RSCRIPT" -e 'cat(sessionInfo()$BLAS)')"
  echo "lowpowermode: $(pmset -g | awk '/lowpowermode/ {print $2}')"
  echo "power: $(pmset -g batt | head -1)"
  echo "settings: $SETTINGS, cores = 5, threads = 1"
} > "$STEP-provenance.txt"
echo $$ > "$STEP-launcher.pid"
case "$STEP" in
  draws)
    "$RSCRIPT" "$REFIT" "$WT" "$OUT" 5 "$ALL" 1-10 1000 1000 > draws.log 2>&1 &
    ;;
  main)
    "$RSCRIPT" "$WT/dev/simstudy/run_study.R" --pkg="$WT" --cores=5 --R=100 \
      --scenarios="$ALL" --out="$OUT" --caffeinate > main.log 2>&1 &
    ;;
  long)
    "$RSCRIPT" "$REFIT" "$WT" "$OUT" 5 "$ALL" 1-5 5000 5000 > long.log 2>&1 &
    ;;
esac
RPID=$!
echo $RPID > "$STEP-rscript.pid"
caffeinate -i -w $RPID &
wait $RPID
status=$?
echo "exit status: $status at $(date '+%F %T')" >> "$STEP-provenance.txt"
echo "ALL DONE $(date '+%F %T') status $status" >> "$STEP.log"

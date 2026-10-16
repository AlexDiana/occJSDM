#!/bin/sh
# One step of the run, named by the first argument; see README.md.
#   draws: replicates 1 to 10 of the six both-means arms and 1 to 3 of the
#          eight comparison arms at 1,000 + 1,000, keeping draws (refit.R).
#   main:  the six both-means arms at R = 100 through the study runner.
#   long:  replicates 1 to 5 of the four main both-means arms at
#          5,000 + 5,000, keeping draws (refit.R).
# 5 worker processes and one thread per fit in every step; the steps run one
# at a time, so the laptop never has more than 5 fitting workers. Launched
# detached from a non-interactive shell. Every output goes to a path that
# .gitignore lists, so each step records a clean tree.
set -u
STEP=$1
WT=/Users/douglasyu/src/occJSDM/.claude/worktrees/pensive-chaum-ade47b
RUN=$WT/dev/simstudy/theta-community-mean-20261016
RSCRIPT=/usr/local/bin/Rscript
BOTH=base_cm_both,occupancy_cm_both,base_rare_cm_both,occupancy_rare_cm_both,base_rare_cm_both_wide,occupancy_rare_cm_both_wide
COMPARE=base,occupancy,base_rare,occupancy_rare,base_cm,occupancy_cm,base_rare_cm,occupancy_rare_cm
MAIN4=base_cm_both,occupancy_cm_both,base_rare_cm_both,occupancy_rare_cm_both
cd "$RUN"
case "$STEP" in
  draws) OUT=$RUN/draws/short; SETTINGS="replicates 1 to 10 of $BOTH and 1 to 3 of $COMPARE, nburn = 1000, niter = 1000, nchain = 2" ;;
  main)  OUT=$RUN/results; SETTINGS="cells $BOTH, R = 100 per cell, nburn = 1000, niter = 1000, nchain = 2" ;;
  long)  OUT=$RUN/draws/long; SETTINGS="replicates 1 to 5 of $MAIN4, nburn = 5000, niter = 5000, nchain = 2" ;;
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
    "$RSCRIPT" "$RUN/refit.R" "$WT" "$OUT" 5 "$BOTH" 1-10 1000 1000 > draws.log 2>&1 &&
    "$RSCRIPT" "$RUN/refit.R" "$WT" "$OUT" 5 "$COMPARE" 1-3 1000 1000 >> draws.log 2>&1 &
    ;;
  main)
    "$RSCRIPT" "$WT/dev/simstudy/run_study.R" --pkg="$WT" --cores=5 --R=100 \
      --scenarios="$BOTH" --out="$OUT" --caffeinate > main.log 2>&1 &
    ;;
  long)
    "$RSCRIPT" "$RUN/refit.R" "$WT" "$OUT" 5 "$MAIN4" 1-5 5000 5000 > long.log 2>&1 &
    ;;
esac
RPID=$!
echo $RPID > "$STEP-rscript.pid"
caffeinate -i -w $RPID &
wait $RPID
status=$?
echo "exit status: $status at $(date '+%F %T')" >> "$STEP-provenance.txt"
echo "ALL DONE $(date '+%F %T') status $status" >> "$STEP.log"

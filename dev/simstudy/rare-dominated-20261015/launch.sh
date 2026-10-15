#!/bin/sh
# Rare-dominated communities: base_rare, occupancy_rare and binary_rare at
# R = 100, paired with replicates 1 to 100 of base, occupancy and binary in the
# R = 300 paper study. 5 workers, one thread per fit. Launched detached from a
# non-interactive shell; see README.md.
set -u
RUN=/Users/douglasyu/src/occJSDM/.worktrees/rare-dominated/dev/simstudy/rare-dominated-20261015
WT=/Users/douglasyu/src/occJSDM/.worktrees/rare-dominated
RSCRIPT=/usr/local/bin/Rscript
CELLS=base_rare,binary_rare,occupancy_rare
cd "$RUN"
if [ -e results ]; then echo "refusing: $RUN/results exists; use a fresh run directory" >&2; exit 1; fi
{
  echo "launched: $(date '+%F %T')"
  echo "revision: $(git -C "$WT" rev-parse HEAD)"
  echo "branch: $(git -C "$WT" rev-parse --abbrev-ref HEAD)"
  echo "dirty files: $(git -C "$WT" status --porcelain | wc -l | tr -d ' ')"
  echo "R: $("$RSCRIPT" -e 'cat(R.version.string)')"
  echo "BLAS: $("$RSCRIPT" -e 'cat(sessionInfo()$BLAS)')"
  echo "lowpowermode: $(pmset -g | awk '/lowpowermode/ {print $2}')"
  echo "power: $(pmset -g batt | head -1)"
  echo "cells: $CELLS"
  echo "settings: R = 100 per cell, cores = 5, threads = 1, nburn = 1000, niter = 1000, nchain = 2"
} > provenance.txt
echo $$ > launcher.pid
"$RSCRIPT" "$WT/dev/simstudy/run_study.R" --pkg="$WT" --cores=5 --R=100 \
  --scenarios="$CELLS" --out="$RUN/results" --caffeinate > run.log 2>&1 &
RPID=$!
echo $RPID > rscript.pid
wait $RPID
status=$?
echo "exit status: $status at $(date '+%F %T')" >> provenance.txt
echo "ALL DONE $(date '+%F %T') status $status" >> run.log

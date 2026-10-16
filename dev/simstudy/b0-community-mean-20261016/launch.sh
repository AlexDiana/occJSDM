#!/bin/sh
# Community mean of B0: base_cm, occupancy_cm, binary_cm and their rare
# counterparts at R = 100, all fitted with listPriors$b0_prior = "community".
# Each pairs with a default-prior arm already run on identical datasets; see
# README.md. 5 workers, one thread per fit. Launched detached from a
# non-interactive shell.
set -u
WT=/Users/douglasyu/src/occJSDM/.claude/worktrees/gracious-almeida-d7dbce
RUN=$WT/dev/simstudy/b0-community-mean-20261016
RSCRIPT=/usr/local/bin/Rscript
CELLS=base_cm,occupancy_cm,binary_cm,base_rare_cm,occupancy_rare_cm,binary_rare_cm
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

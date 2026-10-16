#!/bin/sh
# The two levers that moved B in the non-spatial candidates run, together:
# no_spatial's datasets fitted with 4 latent traits and 5000 + 5000 iterations
# per chain, at R = 100. One cell, 5 workers, one thread per fit; the cell's
# own MCMC settings replace the runner's. Launched detached from a
# non-interactive shell; see README.md.
set -u
WT=/Users/douglasyu/src/occJSDM/.claude/worktrees/sharp-wilson-52e053
RUN=$WT/dev/simstudy/b-gt4-long-20261016
RSCRIPT=/usr/local/bin/Rscript
CELLS=no_spatial_gt4_long
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
  echo "load at launch: $(uptime | sed 's/.*load averages*: //')"
  echo "busiest other process: $(ps -Ao pcpu,comm -r | sed -n 2p | sed 's/^ *//')"
  echo "cells: $CELLS"
  echo "settings: R = 100, cores = 5, threads = 1; no_spatial_gt4_long nburn = 5000, niter = 5000, nchain = 2, fit_gt = 4"
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

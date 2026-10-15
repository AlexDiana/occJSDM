#!/bin/sh
# Full paper simulation grid on the combined branch of PRs 31 and 32, R = 300.
# Production ten cells, 5 replicate processes, one thread per fit, no time cap.
# Launched detached from a non-interactive shell; see README.md for the command.
# To stop it: kill $(cat rscript.pid); the PSOCK workers exit with the master.
# To resume after an interruption: add --resume to the Rscript line below and relaunch.
set -u
RUN=/Users/douglasyu/src/occJSDM/.worktrees/paper-study-combined/dev/simstudy/paper-study-20261015
WT=/Users/douglasyu/src/occJSDM/.worktrees/paper-study-combined
RSCRIPT=/usr/local/bin/Rscript
CELLS=base,binary,d_overfit,d_underfit,low_information,occupancy,primers_3,spatial_isolated,species_20,traits_isolated
R_REPS=300
cd "$RUN"
if [ -e results ] && [ -z "${RESUME:-}" ]; then
  echo "refusing: $RUN/results exists; set RESUME=1 to continue it with --resume" >&2
  exit 1
fi
{
  echo "launched: $(date '+%F %T')"
  echo "worktree: $WT"
  echo "revision: $(git -C "$WT" rev-parse HEAD)"
  echo "branch: $(git -C "$WT" rev-parse --abbrev-ref HEAD)"
  echo "dirty files: $(git -C "$WT" status --porcelain | wc -l | tr -d ' ')"
  echo "merged PR heads: 7a47d1b (PR 31, claude/parallel-chains) and 375b387 (PR 32, claude/symmetric-gram) onto main d156492"
  echo "R: $("$RSCRIPT" -e 'cat(R.version.string)')"
  echo "BLAS: $("$RSCRIPT" -e 'cat(sessionInfo()$BLAS)')"
  echo "lowpowermode: $(pmset -g | awk '/lowpowermode/ {print $2}')"
  echo "power: $(pmset -g batt | head -1)"
  echo "cells: $CELLS"
  echo "settings: R = $R_REPS per cell, cores = 5, threads = 1, nburn = 1000, niter = 1000, nchain = 2, no cap${RESUME:+, resumed}"
} >> provenance.txt
echo $$ > launcher.pid
"$RSCRIPT" "$WT/dev/simstudy/run_study.R" --pkg="$WT" --cores=5 --R=$R_REPS \
  --scenarios="$CELLS" --out="$RUN/results" --caffeinate ${RESUME:+--resume} >> run.log 2>&1 &
RPID=$!
echo $RPID > rscript.pid
wait $RPID
status=$?
echo "exit status: $status at $(date '+%F %T')" >> provenance.txt
echo "ALL DONE $(date '+%F %T') status $status" >> run.log

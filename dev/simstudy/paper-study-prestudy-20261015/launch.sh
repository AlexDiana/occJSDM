#!/bin/sh
# Timing pre-study of the paper simulation grid on the combined branch of PRs 31 and 32.
# Runs the production ten cells at R = 7 on 5 cores, one thread per fit, capped at 30 minutes.
# Launched detached from a non-interactive shell; see README.md for the command.
set -u
RUN=/Users/douglasyu/src/occJSDM/.worktrees/paper-study-combined/dev/simstudy/paper-study-prestudy-20261015
WT=/Users/douglasyu/src/occJSDM/.worktrees/paper-study-combined
RSCRIPT=/usr/local/bin/Rscript
CELLS=base,binary,d_overfit,d_underfit,low_information,occupancy,primers_3,spatial_isolated,species_20,traits_isolated
cd "$RUN"
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
  echo "settings: R = 7 per cell, cores = 5, threads = 1, nburn = 1000, niter = 1000, nchain = 2, cap = 1800 s"
} > provenance.txt
echo $$ > launcher.pid
/opt/homebrew/bin/timeout --kill-after=60 1800 \
  "$RSCRIPT" "$WT/dev/simstudy/run_study.R" --pkg="$WT" --cores=5 --R=7 \
  --scenarios="$CELLS" --out="$RUN/results" --caffeinate > run.log 2>&1
status=$?
echo "exit status: $status at $(date '+%F %T') (124 means the 30-minute cap fired)" >> provenance.txt
echo "ALL DONE $(date '+%F %T') status $status" >> run.log

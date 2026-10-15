#!/bin/sh
# Range-inside-grid diagnostic: base and binary copies with l_s = 0.05 and
# 0.075 raw, R = 50, 5 replicate processes, one thread per fit, no time cap.
# Launched detached from a non-interactive shell; see README.md for the command.
# To stop it: kill $(cat rscript.pid); the PSOCK workers exit with the master.
# To resume after an interruption: RESUME=1 sh launch.sh (adds --resume).
set -u
WT=/Users/douglasyu/src/occJSDM/.claude/worktrees/gracious-heyrovsky-2adc4e
RUN=$WT/dev/simstudy/range-inside-grid-20261015
RSCRIPT=/usr/local/bin/Rscript
CELLS=base_ls05,binary_ls05,base_ls075,binary_ls075
R_REPS=50
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
  echo "package code (R/, src/) identical to the production run's f3668e0"
  echo "R: $("$RSCRIPT" -e 'cat(R.version.string)')"
  echo "BLAS: $("$RSCRIPT" -e 'cat(sessionInfo()$BLAS)')"
  echo "lowpowermode: $(pmset -g | awk '/lowpowermode/ {print $2}')"
  echo "power: $(pmset -g batt | head -1)"
  echo "cells: $CELLS"
  echo "settings: R = $R_REPS per cell, cores = 5, threads = 1, nburn = 1000, niter = 1000, nchain = 2, no cap${RESUME:+, resumed}"
} >> provenance.txt
echo $$ > launcher.pid
"$RSCRIPT" "$WT/dev/simstudy/run_study.R" --pkg="$WT" --cores=5 --threads=1 --R=$R_REPS \
  --scenarios="$CELLS" --out="$RUN/results" --caffeinate ${RESUME:+--resume} >> run.log 2>&1 &
RPID=$!
echo $RPID > rscript.pid
wait $RPID
status=$?
echo "exit status: $status at $(date '+%F %T')" >> provenance.txt
echo "ALL DONE $(date '+%F %T') status $status" >> run.log

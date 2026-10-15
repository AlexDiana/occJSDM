#!/bin/sh
# The field-sample ladder on current code: M2, M5, M10, M20 and K30 at R = 100,
# 5 replicate processes, one thread per fit. Refuses to start while any other
# simulation study is running, so it never competes with one for the CPU.
# Launch detached from a non-interactive shell; see README.md.
set -u
RUN=/Users/douglasyu/src/occJSDM/.worktrees/field-sample-ladder/dev/simstudy/field-sample-ladder
WT=/Users/douglasyu/src/occJSDM/.worktrees/field-sample-ladder
RSCRIPT=/usr/local/bin/Rscript
CELLS=M2,M5,M10,M20,K30
cd "$RUN"
busy=$(pgrep -f 'dev/simstudy/run_study.R' || true)
if [ -n "$busy" ]; then
  echo "refusing: another simulation study is running (pids: $(echo $busy | tr '\n' ' '))" >&2
  pgrep -fl 'dev/simstudy/run_study.R' >&2
  exit 1
fi
if [ -e results ] && [ -z "${RESUME:-}" ]; then
  echo "refusing: $RUN/results exists; set RESUME=1 to continue it with --resume" >&2
  exit 1
fi
if [ -n "$(git -C "$WT" status --porcelain)" ]; then
  echo "refusing: the worktree has uncommitted changes, so the run's provenance would be dirty" >&2
  git -C "$WT" status --short >&2
  exit 1
fi
{
  echo "launched: $(date '+%F %T')${RESUME:+ (resumed)}"
  echo "revision: $(git -C "$WT" rev-parse HEAD)"
  echo "branch: $(git -C "$WT" rev-parse --abbrev-ref HEAD)"
  echo "R: $("$RSCRIPT" -e 'cat(R.version.string)')"
  echo "BLAS: $("$RSCRIPT" -e 'cat(sessionInfo()$BLAS)')"
  echo "lowpowermode: $(pmset -g | awk '/lowpowermode/ {print $2}')"
  echo "power: $(pmset -g batt | head -1)"
  echo "cells: $CELLS"
  echo "settings: R = 100 per cell, cores = 5, threads = 1, nburn = 1000, niter = 1000, nchain = 2"
} >> provenance.txt
echo $$ > launcher.pid
"$RSCRIPT" "$WT/dev/simstudy/run_study.R" --pkg="$WT" --cores=5 --R=100 \
  --scenarios="$CELLS" --out="$RUN/results" --caffeinate ${RESUME:+--resume} >> run.log 2>&1 &
RPID=$!
echo $RPID > rscript.pid
wait $RPID
status=$?
echo "exit status: $status at $(date '+%F %T')" >> provenance.txt
echo "ALL DONE $(date '+%F %T') status $status" >> run.log

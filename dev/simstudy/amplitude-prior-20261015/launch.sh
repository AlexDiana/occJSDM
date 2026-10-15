#!/bin/sh
# Free-amplitude diagnostic: base and binary copies fitted with the half-Cauchy
# spatial amplitude prior, at the production range and at 0.075 raw, R = 50,
# 5 replicate processes, one thread per fit. Then amp-diag.R refits
# replicates 1 to 10 of eight cells to record amplitude, range and spatial
# share. Both run in sequence so they do not compete for cores.
# Launched detached from a non-interactive shell; see README.md for the command.
# To stop it: kill $(cat launcher.pid) $(cat rscript.pid).
set -u
WT=/Users/douglasyu/src/occJSDM/.claude/worktrees/gracious-heyrovsky-2adc4e
RUN=$WT/dev/simstudy/amplitude-prior-20261015
RSCRIPT=/usr/local/bin/Rscript
CELLS=base_hc,binary_hc,base_ls075_hc,binary_ls075_hc
R_REPS=50
cd "$RUN"
if [ -e results ]; then
  echo "refusing: $RUN/results exists" >&2
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
  echo "settings: R = $R_REPS per cell, cores = 5, threads = 1, nburn = 1000, niter = 1000, nchain = 2, no cap"
} >> provenance.txt
echo $$ > launcher.pid
"$RSCRIPT" "$WT/dev/simstudy/run_study.R" --pkg="$WT" --cores=5 --threads=1 --R=$R_REPS \
  --scenarios="$CELLS" --out="$RUN/results" --caffeinate >> run.log 2>&1 &
RPID=$!
echo $RPID > rscript.pid
wait $RPID
status=$?
echo "run exit status: $status at $(date '+%F %T')" >> provenance.txt
echo "RUN DONE $(date '+%F %T') status $status" >> run.log
echo "amp-diag start: $(date '+%F %T')" >> provenance.txt
cd "$WT"
caffeinate -i "$RSCRIPT" "$RUN/amp-diag.R" "$RUN/results/amp-diag.rds" > "$RUN/amp-diag.log" 2>&1
status=$?
echo "amp-diag exit status: $status at $(date '+%F %T')" >> "$RUN/provenance.txt"
echo "ALL DONE $(date '+%F %T') status $status" >> "$RUN/run.log"

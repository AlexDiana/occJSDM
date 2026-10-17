#!/bin/sh
# The prior on sigma_b: five arms on no_spatial's datasets at R = 100, each
# changing sigma_b's prior through the new listPriors option, three of them at
# 5000 + 5000 iterations per chain and one with 4 latent traits, so every arm
# pairs with an arm already run (no_spatial, no_spatial_long or
# no_spatial_gt4_long). 5 workers, one thread per fit; a cell's own MCMC
# settings replace the runner's. Launched detached from a non-interactive
# shell; see README.md.
set -u
WT=/Users/douglasyu/src/occJSDM/.worktrees/sigma-b-prior
RUN=$WT/dev/simstudy/sigma-b-prior-20261017
RSCRIPT=/usr/local/bin/Rscript
CELLS=no_spatial_sbhc,no_spatial_sbig2,no_spatial_sbhc_long,no_spatial_sbig2_long,no_spatial_sbhc_gt4_long
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
  echo "other R fits at launch: $(ps -axo command | grep -c '[R]script.*run_study')"
  echo "cells: $CELLS"
  echo "settings: R = 100 per cell, cores = 5, threads = 1, nchain = 2; no_spatial_sbhc and no_spatial_sbig2 nburn = 1000, niter = 1000; no_spatial_sbhc_long, no_spatial_sbig2_long and no_spatial_sbhc_gt4_long nburn = 5000, niter = 5000; no_spatial_sbhc_gt4_long fit_gt = 4; sbhc arms sigma_b_prior = half_cauchy, sigma_b_scale = 1; sbig2 arms a_sigmab = 2, b_sigmab = 1/9"
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

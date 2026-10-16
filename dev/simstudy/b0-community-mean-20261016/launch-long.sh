#!/bin/sh
# Long-chain check: replicates 1 to 5 of base_rare_cm, occupancy_rare_cm,
# base_cm, occupancy_cm, base_rare and occupancy_rare at 5,000 + 5,000,
# keeping draws; see refit_long.R and README.md. Launched detached from a
# non-interactive shell.
set -u
WT=/Users/douglasyu/src/occJSDM/.claude/worktrees/gracious-almeida-d7dbce
RUN=$WT/dev/simstudy/b0-community-mean-20261016
SCRIPT=$1
RSCRIPT=/usr/local/bin/Rscript
# Everything goes under draws/, which git ignores, so the main run, still in
# flight, records a clean tree when it writes.
OUT=$RUN/draws/long
if [ -e "$OUT" ]; then echo "refusing: $OUT exists" >&2; exit 1; fi
mkdir -p "$OUT"
cd "$OUT"
{
  echo "launched: $(date '+%F %T')"
  echo "revision: $(git -C "$WT" rev-parse HEAD)"
  echo "branch: $(git -C "$WT" rev-parse --abbrev-ref HEAD)"
  echo "dirty files: $(git -C "$WT" status --porcelain | wc -l | tr -d ' ')"
  echo "script: $SCRIPT md5 $(md5 -q "$SCRIPT")"
  echo "R: $("$RSCRIPT" -e 'cat(R.version.string)')"
  echo "lowpowermode: $(pmset -g | awk '/lowpowermode/ {print $2}')"
  echo "power: $(pmset -g batt | head -1)"
  echo "settings: replicates 1 to 5, 6 arms, cores = 5, threads = 1, nburn = 5000, niter = 5000, nchain = 2"
} > long-provenance.txt
echo $$ > long-launcher.pid
"$RSCRIPT" "$SCRIPT" "$WT" "$OUT" 5 > long.log 2>&1 &
RPID=$!
echo $RPID > long-rscript.pid
caffeinate -i -w $RPID &
wait $RPID
status=$?
echo "exit status: $status at $(date '+%F %T')" >> long-provenance.txt
echo "ALL DONE $(date '+%F %T') status $status" >> long.log

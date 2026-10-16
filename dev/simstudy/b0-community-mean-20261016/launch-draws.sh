#!/bin/sh
# Mixing check: refits replicates 1 to 10 of the six community arms and their
# six default counterparts, keeping draws; see refit_draws.R and README.md.
# Launched detached from a non-interactive shell.
set -u
WT=/Users/douglasyu/src/occJSDM/.claude/worktrees/gracious-almeida-d7dbce
RUN=$WT/dev/simstudy/b0-community-mean-20261016
RSCRIPT=/usr/local/bin/Rscript
cd "$RUN"
if [ -e draws ]; then echo "refusing: $RUN/draws exists" >&2; exit 1; fi
{
  echo "launched: $(date '+%F %T')"
  echo "revision: $(git -C "$WT" rev-parse HEAD)"
  echo "branch: $(git -C "$WT" rev-parse --abbrev-ref HEAD)"
  echo "dirty files: $(git -C "$WT" status --porcelain | wc -l | tr -d ' ')"
  echo "R: $("$RSCRIPT" -e 'cat(R.version.string)')"
  echo "lowpowermode: $(pmset -g | awk '/lowpowermode/ {print $2}')"
  echo "power: $(pmset -g batt | head -1)"
  echo "settings: replicates 1 to 10, 12 arms, cores = 5, threads = 1, nburn = 1000, niter = 1000, nchain = 2"
} > draws-provenance.txt
echo $$ > draws-launcher.pid
"$RSCRIPT" "$RUN/refit_draws.R" "$WT" "$RUN/draws" 5 > draws.log 2>&1 &
RPID=$!
echo $RPID > draws-rscript.pid
caffeinate -i -w $RPID &
wait $RPID
status=$?
echo "exit status: $status at $(date '+%F %T')" >> draws-provenance.txt
echo "ALL DONE $(date '+%F %T') status $status" >> draws.log

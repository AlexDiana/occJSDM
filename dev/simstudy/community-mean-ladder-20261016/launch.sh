#!/bin/sh
# One step of the run, named by the first argument; see README.md.
#   draws: replicates 1 to 10 of the six community-mean arms at 1,000 + 1,000,
#          keeping draws (the both-means run's refit.R).
#   main:  the nine arms at R = 100 through the study runner.
#   long:  replicates 1 to 5 of the arms in draws/long-arms.txt (written by
#          mu_b0_check.R) at 5,000 + 5,000, keeping draws.
# 5 worker processes and one thread per fit in every step; the steps run one
# at a time, so the laptop never has more than 5 fitting workers. Launched
# detached from a non-interactive shell. Every output goes to a path that
# .gitignore lists, so each step records a clean tree. Exits with the R
# process's status.
set -u
STEP=$1
WT=/Users/douglasyu/src/occJSDM/.claude/worktrees/elastic-cohen-deb516
RUN=$WT/dev/simstudy/community-mean-ladder-20261016
REFIT=$WT/dev/simstudy/theta-community-mean-20261016/refit.R
RSCRIPT=/usr/local/bin/Rscript
CM=M5_cm,M10_cm,M20_cm,base_rare_cm_m5,base_rare_cm_m10,base_rare_cm_m20
ALL=$CM,base_rare_m5,base_rare_m10,base_rare_m20
cd "$RUN"
case "$STEP" in
  draws) OUT=$RUN/draws/short; SETTINGS="replicates 1 to 10 of $CM, nburn = 1000, niter = 1000, nchain = 2" ;;
  main)  OUT=$RUN/results; SETTINGS="cells $ALL, R = 100 per cell, nburn = 1000, niter = 1000, nchain = 2" ;;
  long)
    if [ ! -s "$RUN/draws/long-arms.txt" ]; then echo "refusing: no arms in draws/long-arms.txt" >&2; exit 1; fi
    LONG=$(cat "$RUN/draws/long-arms.txt")
    OUT=$RUN/draws/long; SETTINGS="replicates 1 to 5 of $LONG, nburn = 5000, niter = 5000, nchain = 2" ;;
  *) echo "unknown step: $STEP" >&2; exit 2 ;;
esac
if [ -e "$OUT" ]; then echo "refusing: $OUT exists" >&2; exit 1; fi
# Study runners, refit and long-chain scripts, and PSOCK workers of any of
# them. A minute's grace lets the previous step's workers exit.
PATTERN='dev/simstudy/run_study[.]R|/refit[^/ ]*[.]R |/long[.]R |workRSOCK'
for i in 1 2 3 4 5 6 7 8 9 10 11 12; do
  busy=$(pgrep -f "$PATTERN" || true)
  [ -z "$busy" ] && break
  sleep 5
done
if [ -n "$busy" ]; then
  echo "refusing: another simulation study, refit or R worker is running (pids: $(echo $busy | tr '\n' ' '))" >&2
  pgrep -fl "$PATTERN" >&2
  exit 1
fi
code_dirty=$(git -C "$WT" status --porcelain -- R src tests DESCRIPTION NAMESPACE dev/simstudy/run_study.R "$REFIT")
if [ -n "$code_dirty" ]; then
  echo "refusing: the package code, tests or harness have uncommitted changes" >&2
  echo "$code_dirty" >&2
  exit 1
fi
{
  echo "step: $STEP"
  echo "launched: $(date '+%F %T')"
  echo "revision: $(git -C "$WT" rev-parse HEAD)"
  echo "branch: $(git -C "$WT" rev-parse --abbrev-ref HEAD)"
  echo "dirty files: $(git -C "$WT" status --porcelain | wc -l | tr -d ' ')"
  echo "R: $("$RSCRIPT" -e 'cat(R.version.string)')"
  echo "BLAS: $("$RSCRIPT" -e 'cat(sessionInfo()$BLAS)')"
  echo "lowpowermode: $(pmset -g | awk '/lowpowermode/ {print $2}')"
  echo "power: $(pmset -g batt | head -1)"
  echo "settings: $SETTINGS, cores = 5, threads = 1"
} > "$STEP-provenance.txt"
echo $$ > "$STEP-launcher.pid"
case "$STEP" in
  draws)
    "$RSCRIPT" "$REFIT" "$WT" "$OUT" 5 "$CM" 1-10 1000 1000 > draws.log 2>&1 &
    ;;
  main)
    "$RSCRIPT" "$WT/dev/simstudy/run_study.R" --pkg="$WT" --cores=5 --R=100 \
      --scenarios="$ALL" --out="$OUT" --caffeinate > main.log 2>&1 &
    ;;
  long)
    "$RSCRIPT" "$REFIT" "$WT" "$OUT" 5 "$LONG" 1-5 5000 5000 > long.log 2>&1 &
    ;;
esac
RPID=$!
echo $RPID > "$STEP-rscript.pid"
caffeinate -i -w $RPID &
wait $RPID
status=$?
echo "exit status: $status at $(date '+%F %T')" >> "$STEP-provenance.txt"
echo "ALL DONE $(date '+%F %T') status $status" >> "$STEP.log"
exit $status

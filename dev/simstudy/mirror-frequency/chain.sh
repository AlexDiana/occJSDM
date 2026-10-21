#!/bin/sh
# Run the mirror-frequency check end to end, detached (README.md, "Launch"):
# stage a (25 repeats), its analysis, then stage b if the stage-a analysis
# found any chain outside the near-truth region, its analysis, and a
# CHAIN-DONE file holding the exit status.
#   sh chain.sh REPO STUDY INPUTS MAXPROCS
set -u
REPO="$1"; STUDY="$2"; INPUTS="$3"; MAXPROCS="$4"
RSCRIPT=/usr/local/bin/Rscript
S="$REPO/dev/simstudy/mirror-frequency"
LOG="$STUDY/logs/chain.log"
mkdir -p "$STUDY/logs"
say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }
run_stage() {
  stage="$1"; reps="$2"
  say "stage $stage: launching repeats $reps with at most $MAXPROCS processes"
  "$RSCRIPT" "$S/launch.R" --repo="$REPO" --study="$STUDY" --inputs-root="$INPUTS" \
    --stage="$stage" --repeats="$reps" --max-procs="$MAXPROCS" >> "$LOG" 2>&1 || { say "stage $stage: launcher exit $?"; return 1; }
  say "stage $stage: analysing"
  "$RSCRIPT" "$S/analyse.R" --repo="$REPO" --study="$STUDY" --inputs-root="$INPUTS" \
    --stage="$stage" >> "$LOG" 2>&1 || { say "stage $stage: analysis failed"; return 1; }
  return 0
}
rm -f "$STUDY/CHAIN-DONE"
status=0
run_stage a 1-25 || status=1
if [ "$status" -eq 0 ]; then
  decision=$(cat "$S/results/stage-a/decision.txt" 2>/dev/null || echo "run_stage_b: unreadable")
  say "stage a decision: $decision"
  case "$decision" in
    "run_stage_b: yes") run_stage b 1-25 || status=1 ;;
    "run_stage_b: no") say "stage b not run: every stage-a chain stayed in the near-truth region" ;;
    *) say "stage b not run: decision unreadable"; status=1 ;;
  esac
fi
say "CHAIN DONE status $status"
echo "$status" > "$STUDY/CHAIN-DONE"
exit "$status"

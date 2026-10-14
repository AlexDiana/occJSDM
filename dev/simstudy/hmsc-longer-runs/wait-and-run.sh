#!/bin/sh
# Wait for the mirror-frequency check to finish, then run the Hmsc refits and
# the comparison, detached (README.md, "Launch").
#   sh wait-and-run.sh REPO ROOT STUDY MIRROR_STUDY MAXPROCS
set -u
REPO="$1"; ROOT="$2"; STUDY="$3"; MIRROR="$4"; MAXPROCS="$5"
RSCRIPT=/usr/local/bin/Rscript
S="$REPO/dev/simstudy/hmsc-longer-runs"
LOG="$STUDY/logs/chain.log"
mkdir -p "$STUDY/logs"
say() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }
rm -f "$STUDY/CHAIN-DONE"
say "waiting for $MIRROR/CHAIN-DONE"
while [ ! -f "$MIRROR/CHAIN-DONE" ]; do sleep 60; done
say "mirror check finished with status $(cat "$MIRROR/CHAIN-DONE"); launching the Hmsc refits, at most $MAXPROCS at once"
status=0
"$RSCRIPT" "$S/launch.R" --repo="$REPO" --root="$ROOT" --study="$STUDY" --max-procs="$MAXPROCS" >> "$LOG" 2>&1 || { say "launcher exit $?"; status=1; }
if [ "$status" -eq 0 ]; then
  say "comparing refits with the saved fits"
  "$RSCRIPT" "$S/compare.R" --repo="$REPO" --root="$ROOT" --study="$STUDY" >> "$LOG" 2>&1 || { say "comparison failed"; status=1; }
fi
say "CHAIN DONE status $status"
echo "$status" > "$STUDY/CHAIN-DONE"
exit "$status"

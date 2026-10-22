#!/bin/sh
# Sample resident memory (KB) of a parent R process and of every PSOCK worker
# on the machine once a second, until OUT.csv.stop appears or the parent exits.
# PSOCK workers are not children of the parent (makeCluster() launches them
# through a shell that exits, so they are reparented to PID 1; checked
# 13 October 2026), so they are found by their command line, which contains
# workRSOCK. summarise.R keeps only the PIDs the fit's execution record names.
# Usage: sample-rss.sh PID OUT.csv
PID="$1"; OUT="$2"
echo "time,pid,role,rss_kb" > "$OUT"
while [ ! -f "$OUT.stop" ] && kill -0 "$PID" 2> /dev/null; do
  now=$(date +%s)
  ps -o rss= -p "$PID" | awk -v t="$now" -v p="$PID" '{print t","p",parent,"$1}' >> "$OUT"
  for c in $(pgrep -f workRSOCK 2> /dev/null); do
    ps -o rss= -p "$c" | awk -v t="$now" -v p="$c" '{print t","p",worker,"$1}' >> "$OUT"
  done
  sleep 1
done

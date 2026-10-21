#!/bin/sh
# Usage: probe-gllvm.sh STUDY DATASET_ID N MAXIT LIMIT_GB [VARIANT]
# Runs probe-gllvm.R, samples the process footprint (top's MEM, which includes compressed pages) every few seconds, kills the fit above
# LIMIT_GB, and appends one line to STUDY/logs/gllvm-probe.csv: dataset,n,maxit,limit_gb,peak_gb,seconds,killed,exit,variant.
STUDY=$1; DATASET=$2; N=$3; MAXIT=$4; LIMIT=$5; VARIANT=${6:-unstructured}
D=$(cd "$(dirname "$0")" && pwd)
OUT=$STUDY/logs/gllvm-probe-$DATASET-n$N-$VARIANT.csv
"$D/run-r" -f "$D/probe-gllvm.R" --args "$STUDY" "$DATASET" "$N" "$MAXIT" "$OUT" "$VARIANT" > "$STUDY/logs/gllvm-probe-$DATASET-n$N-$VARIANT.log" 2>&1 &
PID=$!; START=$(date +%s); PEAK=0; KILLED=no
while kill -0 $PID 2>/dev/null; do
  GB=$(top -l 1 -pid $PID -stats mem 2>/dev/null | tail -1 | awk '{v=$1; gsub(/[+-]/,"",v); u=substr(v,length(v),1); x=substr(v,1,length(v)-1)+0; if (u=="G") print x; else if (u=="M") print x/1024; else if (u=="K") print x/1048576; else print 0}')
  PEAK=$(echo "$PEAK ${GB:-0}" | awk '{print ($2>$1)?$2:$1}')
  if [ "$(echo "$PEAK $LIMIT" | awk '{print ($1>=$2)?1:0}')" = 1 ]; then kill -KILL $PID 2>/dev/null; KILLED=yes; break; fi
  sleep 4
done
wait $PID 2>/dev/null; RC=$?
END=$(date +%s)
[ -f "$STUDY/logs/gllvm-probe.csv" ] || echo "dataset,n,maxit,limit_gb,peak_gb,seconds,killed,exit,variant" > "$STUDY/logs/gllvm-probe.csv"
echo "$DATASET,$N,$MAXIT,$LIMIT,$PEAK,$((END-START)),$KILLED,$RC,$VARIANT" >> "$STUDY/logs/gllvm-probe.csv"
tail -1 "$STUDY/logs/gllvm-probe.csv"

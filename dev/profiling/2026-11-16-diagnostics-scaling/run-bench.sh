#!/bin/sh
# Run bench-diagnostics.R at each block size against one installed build and
# append everything to a log, with the conditions the timings were taken under.
#
# Usage:
#   sh run-bench.sh <label> <lib> <logfile> [sizes] [max_minutes] [max_gb]
#
# <lib> holds the occJSDM build to measure. Install an export of the commit,
# never the working tree of the main checkout (see AGENTS.md):
#   git archive <commit> | tar -x -C <src>
#   R CMD INSTALL --no-docs --library=<lib> <src>
#
# [sizes] is a quoted list of block sizes, by default "300 3000 12500".
# [max_minutes] and [max_gb] stop a run that outgrows the laptop: each run is
# polled every 5 seconds and killed once it has run longer than the first or
# has a larger memory footprint than the second. A killed run is logged as
# such, with the time and footprint it had reached. Defaults: 60 minutes,
# 20 GB. The footprint (footprint(1), what Activity Monitor shows) is used and
# not the resident set size, because macOS compresses idle pages: at 6,000
# parameters the old code peaked at 6.7 GB resident but 11.1 GB of footprint.
#
# Each run is a fresh R process, so one size's garbage cannot inflate the
# next size's peak. /usr/bin/time -l reports the process peaks.

set -u
label=$1
lib=$2
log=$3
sizes=${4:-"300 3000 12500"}
max_minutes=${5:-60}
max_gb=${6:-20}

here=$(cd "$(dirname "$0")" && pwd)
max_seconds=$((max_minutes * 60))
max_mb=$((max_gb * 1024))

{
  echo "=== $label: $(date '+%Y-%m-%d %H:%M:%S %Z') ==="
  echo "lowpowermode: $(pmset -g | grep -i lowpowermode | tr -s ' ')"
  echo "power source: $(pmset -g batt | head -1)"
  echo "machine: $(sysctl -n machdep.cpu.brand_string), $(($(sysctl -n hw.memsize) / 1073741824)) GB, macOS $(sw_vers -productVersion)"
  echo "R: $(R --version | head -1)"
  echo "BLAS: $(Rscript --vanilla -e 'cat(extSoftVersion()[["BLAS"]])')"
  echo "lib: $lib"
} >> "$log"

for n in $sizes; do
  echo "--- $label n_params=$n ---" >> "$log"
  out=$(mktemp)
  /usr/bin/time -l Rscript --vanilla "$here/bench-diagnostics.R" "$label" "$n" "$lib" > "$out" 2>&1 &
  pid=$!
  start=$(date +%s)
  peak_mb=0
  killed=""
  while kill -0 "$pid" 2>/dev/null; do
    sleep 5
    # time(1) is the parent of the R process, which is the one to measure.
    rpid=$(pgrep -P "$pid" | head -1)
    mb=0
    if [ -n "$rpid" ]; then
      mb=$(footprint "$rpid" 2>/dev/null | awk '/phys_footprint:/ { v = $2; if ($3 == "GB") v *= 1024; if ($3 == "KB") v /= 1024; print int(v); exit }')
    fi
    mb=${mb:-0}
    [ "$mb" -gt "$peak_mb" ] && peak_mb=$mb
    now=$(date +%s)
    if [ $((now - start)) -gt "$max_seconds" ]; then killed="time limit of $max_minutes min"; fi
    if [ "$mb" -gt "$max_mb" ]; then killed="memory limit of $max_gb GB"; fi
    if [ -n "$killed" ]; then
      pkill -P "$pid" 2>/dev/null
      kill "$pid" 2>/dev/null
      break
    fi
  done
  wait "$pid" 2>/dev/null
  end=$(date +%s)
  if [ -n "$killed" ]; then
    echo "STOPPED label=$label n_params=$n reason=\"$killed\" elapsed_s=$((end - start)) polled_peak_footprint_mb=$peak_mb" >> "$log"
  fi
  grep -v -e 'page reclaims' -e 'page faults' -e 'swaps' -e 'block ' -e 'messages ' \
    -e 'signals' -e 'context switches' -e 'instructions' -e 'cycles' \
    -e 'shared memory' -e 'unshared' "$out" >> "$log"
  rm -f "$out"
done
echo >> "$log"

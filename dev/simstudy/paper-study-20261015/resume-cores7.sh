#!/bin/sh
# Resume the R = 300 paper study with 7 replicate processes instead of 5.
# Kept outside the worktree so the runner's clean-tree check stays clean; its text is copied into provenance.txt.
set -u
RUN=/Users/douglasyu/src/occJSDM/.worktrees/paper-study-combined/dev/simstudy/paper-study-20261015
WT=/Users/douglasyu/src/occJSDM/.worktrees/paper-study-combined
RSCRIPT=/usr/local/bin/Rscript
CELLS=base,binary,d_overfit,d_underfit,low_information,occupancy,primers_3,spatial_isolated,species_20,traits_isolated
cd "$RUN"
{
  echo "resumed: $(date '+%F %T') with --cores=7 --resume, at Doug's request"
  echo "resume revision: $(git -C "$WT" rev-parse HEAD); dirty files: $(git -C "$WT" status --porcelain | wc -l | tr -d ' ')"
  echo "resume power: $(pmset -g batt | head -1); lowpowermode $(pmset -g | awk '/lowpowermode/ {print $2}')"
  echo "resume command: $RSCRIPT $WT/dev/simstudy/run_study.R --pkg=$WT --cores=7 --R=300 --scenarios=$CELLS --out=$RUN/results --caffeinate --resume"
} >> provenance.txt
echo $$ > launcher.pid
"$RSCRIPT" "$WT/dev/simstudy/run_study.R" --pkg="$WT" --cores=7 --R=300 \
  --scenarios="$CELLS" --out="$RUN/results" --caffeinate --resume >> run.log 2>&1 &
RPID=$!
echo $RPID > rscript.pid
wait $RPID
status=$?
echo "exit status: $status at $(date '+%F %T')" >> provenance.txt
echo "ALL DONE $(date '+%F %T') status $status" >> run.log

#!/bin/sh
# Bisection refits for the calibration report: base and binary, R = 20,
# at three main commits between the 10 August run and the October run.
# The second stage used the same script with the commit list replaced by
# 54b5cca f1d2e12 e9fd2d3 (the base and two commits of the spatial branch,
# PR 8) and progress written to progress2.log.
set -u
WT=/Users/douglasyu/src/occJSDM/.claude/worktrees/sharp-elgamal-46c8c5
B=/private/tmp/claude-501/-Users-douglasyu-src-occJSDM--claude-worktrees-sharp-elgamal-46c8c5/1042439f-5d8b-4b4a-8f51-2fea77c4afc7/scratchpad/bisect
cd "$B"
for sha in 5382460 9a784a9 cc1b4e6; do
  mkdir -p "tree-$sha"
  git -C "$WT" archive "$sha" | tar -x -C "tree-$sha"
  echo "start $sha $(date '+%F %T')" >> progress.log
  /usr/local/bin/Rscript "$WT/dev/simstudy/run_study.R" --pkg="$B/tree-$sha" \
    --cores=5 --R=20 --scenarios=base,binary --out="$B/out-$sha" --caffeinate \
    > "run-$sha.log" 2>&1
  echo "done $sha status $? $(date '+%F %T')" >> progress.log
done
echo "ALL DONE $(date '+%F %T')" >> progress.log

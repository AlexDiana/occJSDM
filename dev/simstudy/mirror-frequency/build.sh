#!/bin/sh
# Build the study library from a git archive export of one revision (README.md,
# "Library"):
#   sh build.sh REPO STUDY REVISION
# Refuses a STUDY that already holds a build. Writes STUDY/source (the package
# files only), STUDY/source-revision.txt, STUDY/library, STUDY/install.log and,
# through fingerprint.R, the fingerprint CSV under the study's results directory.
set -eu
REPO="$1"; STUDY="$2"; REVISION="$3"
R=/usr/local/bin/R
RSCRIPT=/usr/local/bin/Rscript
for stale in "$STUDY/source" "$STUDY/library/occJSDM" "$STUDY/source-revision.txt"; do
  if [ -e "$stale" ]; then echo "refusing to reuse $STUDY: $stale exists" >&2; exit 1; fi
done
mkdir -p "$STUDY/source" "$STUDY/library"
git -C "$REPO" rev-parse --verify "$REVISION^{commit}" > "$STUDY/source-revision.txt"
git -C "$REPO" archive "$(cat "$STUDY/source-revision.txt")" R src DESCRIPTION NAMESPACE LICENSE man data | tar -x -C "$STUDY/source"
"$R" CMD INSTALL --library="$STUDY/library" "$STUDY/source" > "$STUDY/install.log" 2>&1
"$RSCRIPT" "$REPO/dev/simstudy/mirror-frequency/fingerprint.R" --repo="$REPO" --study="$STUDY"
echo "built $(cat "$STUDY/source-revision.txt") into $STUDY/library"

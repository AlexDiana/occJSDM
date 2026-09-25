#!/usr/bin/env python3
"""Fit, score, verify/export and render the approved study without a live chat."""
import argparse
import csv
import datetime
import json
from pathlib import Path
import subprocess
import sys

p = argparse.ArgumentParser()
p.add_argument('root', type=Path)
p.add_argument('code', type=Path)
p.add_argument('repo', type=Path)
p.add_argument('--launcher', type=Path, required=True)
p.add_argument('--workers', type=int, default=4)
a = p.parse_args()
root, code, repo = a.root.resolve(), a.code.resolve(), a.repo.resolve()

def run(command):
    print('Running:', ' '.join(map(str, command)), flush=True)
    subprocess.run([str(x) for x in command], check=True, cwd=repo)

try:
    run([sys.executable, code / 'extension/run-study.py', root, code,
         '--launcher', a.launcher, '--workers', a.workers, '--score'])
    run([a.launcher, code / 'extension/export.R', root, repo])
    run([a.launcher, code / 'extension/verify.R', root, repo])
    run([a.launcher, '-e', 'rmarkdown::render("vignettes/occJSDM-lesson-4.Rmd", quiet=TRUE); '
         'rmarkdown::render("vignettes/occJSDM-lesson-4.Rmd", '
         'output_format=rmarkdown::github_document(html_preview=FALSE), quiet=TRUE)'])
    final = dict(state='all_attempted_exported_and_verified', finished=datetime.datetime.now().isoformat())
except Exception as error:
    final = dict(state='needs_attention', error=str(error), finished=datetime.datetime.now().isoformat())
    (root / 'completion.json').write_text(json.dumps(final, indent=2))
    raise
(root / 'completion.json').write_text(json.dumps(final, indent=2))
print(json.dumps(final), flush=True)

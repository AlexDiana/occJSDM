"""Run the complete project pytest collection in three independent processes."""
import json
import hashlib
from pathlib import Path
import subprocess
import sys
import time

root=Path(__file__).resolve().parents[3]
production=[root/'src/tools/quickstart.py',root/'src/tools/runtime.py',root/'src/tools/__init__.py',root/'src/r_scripts/quickstart.R']
def hashes():return {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in production}
before=hashes()
collection=subprocess.run([sys.executable,'-m','pytest','--collect-only','-q'],cwd=root,text=True,capture_output=True)
if collection.returncode:
    print(collection.stdout+collection.stderr);sys.exit(collection.returncode)
nodes=[line for line in collection.stdout.splitlines() if '::test_' in line]
if not nodes:raise SystemExit('No pytest cases collected')
shards=[nodes[i::3] for i in range(3)]
started=time.monotonic()
workers=[]
for i,targets in enumerate(shards):
    log=root/'reports'/f'implementation-green-shard-{i}.log'
    handle=log.open('w')
    command=[sys.executable,'-m','pytest','-q',*targets]
    proc=subprocess.Popen(command,cwd=root,stdout=handle,stderr=subprocess.STDOUT)
    workers.append((proc,handle,command,log))
results=[]
for proc,handle,command,log in workers:
    code=proc.wait();handle.close()
    results.append(dict(command=command,returncode=code,log=str(log.relative_to(root)),tail=log.read_text().splitlines()[-6:]))
after=hashes()
report=dict(production_before=before,production_after=after,production_unchanged=before==after,collected=len(nodes),shards=results,elapsed_seconds=time.monotonic()-started)
(root/'reports/implementation-green-suite.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
sys.exit(0 if before==after and all(x['returncode']==0 for x in results) else 1)

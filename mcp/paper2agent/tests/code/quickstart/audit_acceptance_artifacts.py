"""Verify publication hashes and call the independent native comparator after stdio."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
root=Path(sys.argv[1]).resolve();selector=sys.argv[2];runtime=Path(sys.argv[3]).resolve();report=Path(sys.argv[4]).resolve()
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
producer_hashes=json.loads((root/'reports/verification/quickstart/producer_hashes.json').read_text())
assert all(sha(root/p)==h for p,h in producer_hashes.items()),'Production changed after independent verification'
if Path(selector).is_file():
    previous={Path(p).name for p in json.loads(Path(selector).read_text())}
    selected=[p for p in (root/'artifacts/calls').iterdir() if p.name not in previous]
else:
    cutoff=float(selector)
    selected=[p for p in (root/'artifacts/calls').iterdir() if (p/'request.json').is_file() and (p/'request.json').stat().st_mtime>=cutoff]
checks=[]
for p in selected:
    if not (p/'manifest.json').is_file():continue
    m=json.loads((p/'manifest.json').read_text())
    assert m['source_revision']==m['runtime']['source_revision']=='b7b7001e56cea8e3931b09c917ea0e29e6cef3c6'
    assert Path(m['runtime']['occJSDM_path']).resolve().is_relative_to(runtime)
    assert sha(m['fit_path'])==m['fit_sha256']
    assert sha(m['input_path'])==m['input_sha256']==sha(p/'input.rds')
    registry=json.loads((root/'artifacts/registry'/(m['fit_id']+'.json')).read_text())
    assert sha(registry['manifest_path'])==registry['manifest_sha256']
    checks.append(m['fit_id'])
run=subprocess.run(['/usr/local/bin/Rscript','--vanilla',str(root/'tests/code/quickstart/audit_acceptance_artifacts.R'),str(root),selector,str(runtime),str(report)],capture_output=True,text=True)
report.with_suffix('.log').write_text(run.stdout+run.stderr)
assert run.returncode==0,run.stdout+run.stderr
result=json.loads(report.read_text());result['publication_hash_checks']=checks
report.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
print(json.dumps({'success':result['success'],'fits':len(result['fits']),'tables':len(result['tables']),'publication_hash_checks':len(checks),'report':str(report)}))

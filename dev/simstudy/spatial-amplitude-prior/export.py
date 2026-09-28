#!/usr/bin/env python3
"""Export audited tables/figures without requiring raw fits to render the report."""
import csv
import hashlib
import os
from pathlib import Path
import shutil
import sys

study = Path(sys.argv[1]).resolve()
repo = Path(sys.argv[2]).resolve()
source = study / 'summary-binary-final'
initial = study / 'summary-binary-initial'
scripts = repo / 'dev/simstudy/spatial-amplitude-prior'
out = scripts / 'results'
out.mkdir(exist_ok=True)

def digest_stream(stream):
    h = hashlib.md5()
    for block in iter(lambda: stream.read(1024 * 1024), b''):
        h.update(block)
    return h.hexdigest()

def rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))

def write(path, records):
    assert records
    with path.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(records[0]))
        w.writeheader();w.writerows(records)

fits = rows(source / 'fits.csv')
audit = rows(source / 'independent-audit.csv')
assert len(fits) == len(audit) == 18
by_id = {(a['key'], a['prior']): a for a in audit}
for f in fits:
    a = by_id[(f['key'], f['prior'])]
    assert a['result_md5'] == f['result_md5'] and a['fit_md5'] == f['fit_md5']
    for path_col, hash_col in [('result_file', 'result_md5'), ('fit_file', 'fit_md5')]:
        with open(f[path_col], 'rb') as stream:
            assert digest_stream(stream) == f[hash_col]

files = ['fits.csv','paired-communities.csv','paired-summary.csv','groups.csv',
         'spatial-diagnostics.csv','selection.csv','extension-gate.csv',
         'initialization-selection.csv','initialization-summary.csv',
         'independent-audit.csv','audit-source.csv',
         'paired-field-recovery.png','spatial-amplitude-priors.png',
         'occupancy-bias.png','illustrative-spatial-patterns.png']
for name in files:
    src = source / name
    assert src.exists(), name
    if src.suffix == '.csv':
        records = rows(src)
        for row in records:
            for key, value in row.items():
                if value.startswith(str(study.parent) + '/'):
                    row[key] = os.path.relpath(value, study)
                elif value.startswith(str(repo) + '/'):
                    row[key] = os.path.relpath(value, scripts)
        write(out / name, records)
    else:
        shutil.copy2(src, out / name)
for name in ['fits.csv','paired-communities.csv','paired-summary.csv','spatial-diagnostics.csv']:
    records = rows(initial / name)
    for row in records:
        for key, value in row.items():
            if value.startswith(str(study.parent) + '/'):
                row[key] = os.path.relpath(value, study)
    write(out / ('initial-' + name), records)
shutil.copy2(study / 'source-revision.txt', out / 'source-revision.txt')
shutil.copy2(study / 'default-equivalence.txt', out / 'default-equivalence.txt')
write(out / 'research-source-md5.csv', [dict(file=p.name, md5=hashlib.md5(p.read_bytes()).hexdigest())
      for p in sorted(scripts.iterdir()) if p.is_file() and p.suffix in ['.R','.Rmd','.py','.md']])
long_pairs = sum(f['prior'] == 'half_cauchy' and f['phase'] == 'long' for f in fits)
flags = sum(int(f['spatial_flag_count']) > 0 for f in fits)
(out / 'README.md').write_text(f'''# Compact spatial-amplitude comparison

Nine binary communities, both priors, 100 supports. {long_pairs} pairs use the longer schedule; {flags} selected fits retain spatial diagnostic flags. Read the report for the scientific conclusion and all limitations.

The paired tables give half-Cauchy minus inverse-gamma changes. Probability errors are proportions; multiply by 100 for percentage points. Communities are the replication units, with three fixed generating-range strata. Interval containment is descriptive and does not establish calibration from nine communities.

`fits.csv` identifies every selected result and fit; `independent-audit.csv` checks each one. The `initial-` tables preserve initial-versus-longer sensitivity. `extension-gate.csv` records the prespecified binary gate. Review native Rhat threshold crossings separately before extending. Figures and `initialization-summary.csv` support rendering the HTML report without raw fits. `artifact-md5.csv` protects this compact bundle; `research-source-md5.csv` records the source scripts.

Paths in result/fit manifests are relative to the current raw archive. Paths beginning `../spatial-targeted-20260927/` refer to the immutable original study. The full archives remain local at `{study}` and `{study.parent / 'spatial-targeted-20260927'}`; posterior draws and installed packages are required to rerun the numerical audit. Production fitting revision is recorded in `source-revision.txt`. Defaults were not changed.

Render `../report.Rmd` with its `summary` parameter set to this directory's absolute path. This bundle supports inspection and rendering; it does not contain every posterior draw.
''')
write(out / 'artifact-md5.csv', [dict(file=p.name,md5=hashlib.md5(p.read_bytes()).hexdigest())
      for p in sorted(out.iterdir()) if p.is_file() and p.name != 'artifact-md5.csv'])
print('Exported', len(list(out.iterdir())), 'compact artifacts to', out)

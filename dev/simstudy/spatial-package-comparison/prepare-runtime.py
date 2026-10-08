#!/usr/bin/env python3
"""Create a private R 4.5 home with corrected recursive shell entry points.
No installed R or global library is changed. Dependencies remain read-only links.
"""
import json
import pathlib
import sys
old = pathlib.Path('/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources')
root = pathlib.Path(sys.argv[1]).resolve()
home = root / 'runtime' / 'R-home'
if home.exists():
    raise SystemExit('Runtime already exists; verify it rather than overwrite it')
(home / 'bin').mkdir(parents=True)
for item in old.iterdir():
    if item.name not in ('bin', 'etc'):
        (home / item.name).symlink_to(item, target_is_directory=item.is_dir())
for item in (old / 'bin').iterdir():
    if item.name == 'Rscript':
        continue
    if item.name == 'R':
        text = item.read_text().replace('/Library/Frameworks/R.framework/Resources', str(home))
        (home / 'bin' / 'R').write_text(text)
        (home / 'bin' / 'R').chmod(0o755)
    else:
        (home / 'bin' / item.name).symlink_to(item, target_is_directory=item.is_dir())
# Recursive Rscript used by package Makevars must also avoid the global wrapper.
script = """#!/usr/bin/python3
import os, pathlib, sys
home = pathlib.Path(__file__).resolve().parent.parent
os.environ['R_HOME'] = str(home)
a = sys.argv[1:]
if '-e' not in a:
    position = next((i for i,x in enumerate(a) if not x.startswith('-')), None)
    if position is not None:
        a = a[:position] + ['--file='+a[position], '--args'] + a[position+1:]
os.execv(str(home/'bin'/'exec'/'R'), ['R','--no-echo','--no-restore']+a)
"""
(home/'bin'/'Rscript').write_text(script)
(home/'bin'/'Rscript').chmod(0o755)
(home/'etc').mkdir()
for item in (old/'etc').iterdir():
    dest = home/'etc'/item.name
    if item.name == 'Makeconf':
        dest.write_text(item.read_text().replace('LIBR = -F/Library/Frameworks/R.framework/.. -framework R', 'LIBR = $(LIBR0) $(LIBR1)'))
    else:
        dest.symlink_to(item, target_is_directory=item.is_dir())
(root / 'library').mkdir(exist_ok=True)
(root / 'logs').mkdir(exist_ok=True)
(root / 'runtime' / 'identity.json').write_text(json.dumps(dict(
    original_home=str(old), private_home=str(home),
    reason='System R framework changed to 4.6; old shell wrapper follows global symlink. Private wrapper pins recursive children to 4.5.',
    modification='Private R/Rscript wrappers and Makeconf LIBR pin 4.5; installed binaries and libraries unchanged'), indent=2)+'\n')
print(home)

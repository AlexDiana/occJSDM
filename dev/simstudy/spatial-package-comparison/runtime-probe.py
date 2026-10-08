#!/usr/bin/env python3
import hashlib, json, pathlib, platform, sys
import numpy, torch
paths={pathlib.Path(sys.executable).resolve()}
for module in [numpy,torch]:
    base=pathlib.Path(module.__file__).parent
    paths.update(p for p in base.rglob('*') if p.is_file() and p.suffix in ('.py','.so','.dylib') and '__pycache__' not in p.parts)
print(json.dumps(dict(python=sys.version,executable=sys.executable,platform=platform.platform(),
                     torch=torch.__version__,numpy=numpy.__version__,cuda_available=torch.cuda.is_available(),
                     files={str(p):hashlib.md5(p.read_bytes()).hexdigest() for p in sorted(paths)})))

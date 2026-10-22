"""A resume must reject a modified result, and preserve completed failed fits."""
import hashlib
import importlib.util
import json
from pathlib import Path
import tempfile

spec = importlib.util.spec_from_file_location('runner', Path(__file__).with_name('run-study.py'))
runner = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runner)
with tempfile.TemporaryDirectory() as temporary:
    root = Path(temporary)
    job = dict(job='example', input_md5='original')
    assert not runner.complete(job, root)
    output = root / 'jobs' / 'example'
    output.mkdir(parents=True)
    result = output / 'result.rds'
    result.write_bytes(b'preserved failed fit')
    (output / 'status.json').write_text(json.dumps(dict(ok=False, input_md5='original',
        result_md5=hashlib.md5(b'preserved failed fit').hexdigest())))
    assert runner.complete(job, root)
    result.write_bytes(b'truncated or modified')
    try:
        runner.complete(job, root)
        raise AssertionError('Modified result accepted')
    except RuntimeError as error:
        assert 'modified result' in str(error)
print('Resume preserves recorded failures and rejects modified results.')

"""Project-confined, file-based transport to the fixed native R dispatcher."""
import asyncio
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import signal
import uuid

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
OPERATIONS = frozenset({'validate_data','fit_model','diagnostics','summarise_fit'})


def project_root() -> Path:
    root = Path(os.environ.get('P2A_PROJECT_ROOT', str(PACKAGE_ROOT))).expanduser().resolve(strict=True)
    if not root.is_dir():
        raise ValueError('P2A_PROJECT_ROOT must be an existing directory')
    return root


def confined_path(value: str, *, exists: bool = True) -> Path:
    raw = Path(value).expanduser()
    if '..' in raw.parts:
        raise ValueError('Parent traversal (..) is not allowed; use a path within the project')
    root = project_root()
    path = (raw if raw.is_absolute() else root / raw).resolve(strict=exists)
    if not path.is_relative_to(root):
        raise ValueError('Path escapes the configured project root, including through a symlink')
    if exists and not path.is_file():
        raise ValueError('Input must be an existing regular file')
    return path


def digest(path: Path) -> str:
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def write_json_atomic(path: Path, data: dict) -> None:
    pending = path.with_name(path.name + '.' + uuid.uuid4().hex + '.pending')
    pending.write_text(json.dumps(data, allow_nan=False, indent=2) + '\n')
    os.replace(pending, path)


def artifact(kind: str, path: Path, description: str) -> dict:
    checked = confined_path(str(path))
    return {'kind':kind, 'description':description, 'path':str(checked), 'sha256':digest(checked)}


def fresh_output(operation: str) -> Path:
    parent = confined_path('artifacts/calls', exists=False)
    parent.mkdir(parents=True, exist_ok=True)
    # Resolve again after mkdir to catch an existing symlink escape.
    parent = confined_path(str(parent), exists=False)
    out = parent / (operation + '-' + uuid.uuid4().hex)
    out.mkdir()
    return out


def native_config() -> tuple[str, Path, float]:
    configured = os.environ.get('P2A_RSCRIPT', 'Rscript')
    executable = shutil.which(configured)
    if executable is None:
        raise ValueError('Rscript unavailable; set P2A_RSCRIPT to an executable Rscript path')
    r_project = Path(os.environ.get('P2A_R_PROJECT', str(PACKAGE_ROOT/'r-runtime'))).expanduser().resolve(strict=True)
    if not r_project.is_dir() or not (r_project/'activate.R').is_file() or not (r_project/'renv/activate.R').is_file() or not (r_project/'renv.lock').is_file():
        raise ValueError('P2A_R_PROJECT requires the isolated runtime activate.R, renv/activate.R and renv.lock; user-library fallback is prohibited')
    timeout = float(os.environ.get('P2A_TIMEOUT_SECONDS','600'))
    if not math.isfinite(timeout) or timeout <= 0:
        raise ValueError('P2A_TIMEOUT_SECONDS must be a finite positive deadline')
    return executable, r_project, timeout


async def _reap(process: asyncio.subprocess.Process) -> None:
    if process.returncode is not None:
        return
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        await asyncio.wait_for(process.wait(), 2)
    except asyncio.TimeoutError:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        await process.wait()


async def run_r(operation: str, request: dict, out: Path) -> dict:
    if operation not in OPERATIONS:
        raise ValueError('Unknown native operation')
    executable, r_project, timeout = native_config()
    request_path = out/'request.json'
    write_json_atomic(request_path, {**request, 'operation':operation, 'output_dir':str(out)})
    env = os.environ.copy()
    env.update(P2A_R_PROJECT=str(r_project), P2A_PROJECT_ROOT=str(project_root()))
    argv = [executable,'--vanilla',str(PACKAGE_ROOT/'src/r_scripts/quickstart.R'),str(r_project/'activate.R'),str(request_path),str(out/'result.json')]
    # Direct log files keep unbounded console output out of memory and protocol stdout.
    with (out/'stdout.log').open('wb') as stdout, (out/'stderr.log').open('wb') as stderr:
        process = await asyncio.create_subprocess_exec(*argv, cwd=project_root(), env=env, stdout=stdout, stderr=stderr, start_new_session=True)
        try:
            await asyncio.wait_for(process.wait(), timeout)
        except asyncio.TimeoutError:
            await _reap(process)
            raise RuntimeError(f'Native R operation exceeded {timeout:g} seconds; no completed fit was published. Logs: {out}') from None
        except asyncio.CancelledError:
            await asyncio.shield(_reap(process))
            raise
    result_path = out/'result.json'
    if not result_path.is_file():
        raise RuntimeError(f'Native R startup or operation failed (exit {process.returncode}); inspect local logs in {out}')
    result = json.loads(result_path.read_text(),parse_constant=lambda s: (_ for _ in ()).throw(ValueError('Nonfinite JSON from R')))
    if process.returncode != 0 or result.get('status') == 'error':
        raise RuntimeError('Native R error: ' + result.get('error','unsuccessful execution') + f'. Local logs: {out}')
    return result


def resolve_fit(fit_id: str) -> tuple[Path, dict]:
    # IDs are generated by the server; never accept RDS paths from an MCP caller.
    if len(fit_id) != 36 or not fit_id.startswith('fit_') or any(c not in '0123456789abcdef' for c in fit_id[4:]):
        raise ValueError('Unknown fit_id; use an identifier returned by fit_model')
    index = confined_path('artifacts/registry/'+fit_id+'.json')
    record = json.loads(index.read_text())
    manifest_path = confined_path(record['manifest_path'])
    if digest(manifest_path) != record['manifest_sha256']:
        raise ValueError('Fit manifest integrity check failed')
    manifest = json.loads(manifest_path.read_text())
    if manifest.get('fit_id') != fit_id or manifest.get('status') != 'complete':
        raise ValueError('Fit manifest is not a completed server-owned fit')
    path = confined_path(manifest['fit_path'])
    if digest(path) != manifest['fit_sha256']:
        raise ValueError('Fit artifact integrity check failed')
    return path, manifest

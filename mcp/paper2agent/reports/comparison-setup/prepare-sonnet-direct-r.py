#!/usr/bin/env python3
"""Remove the occJSDM MCP route only from Doug's Condition C Cloud copy.

--apply-copy is operator confirmation, not automatic Cloud identity detection.
Archives are private temporary recovery copies, not durable Cloud backups.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import uuid

PROJECT_ID = '13080637'
REVISION = 'b7b7001e56cea8e3931b09c917ea0e29e6cef3c6'
PINS = {
    'trial/data.rds': '86a183af12fa6c8e7f2d87abd01decc30f6280bb8d4002c8ccc9dee4fef7b584',
    'trial/fit-short.rds': '446556b0ab908541034af3fb5369864e8164531086e25d3256d19f0078fbd5f8',
}
MCP_PATHS = ('mcp/paper2agent/src', 'mcp/paper2agent/occJSDM-env')
SETTINGS = '.posit/assistant/settings.json'


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def inventory(path):
    """Compare regular file bytes and link text without following directory links."""
    result = {}
    for folder, dirs, files in os.walk(path, followlinks=False):
        for name in dirs + files:
            p = Path(folder) / name
            rel = str(p.relative_to(path))
            if p.is_symlink():
                result[rel] = ('link', os.readlink(p))
            elif p.is_file():
                result[rel] = ('file', sha(p))
            elif p.is_dir():
                result[rel] = ('directory',)
            else:
                raise ValueError('Unsupported file type in archive: ' + rel)
    return result


def check_pins(root, pins):
    for rel, digest in pins.items():
        p = root / rel
        if not p.resolve(strict=True).is_relative_to(root) or sha(p) != digest:
            raise ValueError('Reference integrity check failed: ' + rel)


def delete_discarded(path):
    # renv/venv may contain owner-read-only directories. Never chmod link targets.
    for folder, dirs, _ in os.walk(path, followlinks=False):
        p = Path(folder)
        if not p.is_symlink():
            p.chmod(p.stat().st_mode | 0o700)
        for name in dirs:
            p = Path(folder) / name
            if not p.is_symlink():
                p.chmod(p.stat().st_mode | 0o700)
    shutil.rmtree(path)


def prepare_direct(root, backup, pins, verify=None):
    root, backup = root.resolve(), backup.resolve()
    if backup.is_relative_to(root) or root.is_relative_to(backup):
        raise ValueError('Archive must be outside the project')
    if any(backup.iterdir()):
        raise ValueError('Archive must be new and empty')
    backup.chmod(0o700)
    check_pins(root, pins)
    settings = root / SETTINGS
    for rel in (*MCP_PATHS, SETTINGS):
        p = root / rel
        if p.is_symlink() or not p.resolve(strict=True).is_relative_to(root):
            raise ValueError('Setup paths must be real paths inside the project')
    original = settings.read_bytes()
    config = json.loads(original)
    servers = config.get('mcpServers', {})
    if not isinstance(servers, dict) or 'occJSDM' not in servers:
        raise ValueError('Expected occJSDM registration missing; no changes made')
    # Archive and verify every removed byte before touching the project.
    saved_settings = backup / SETTINGS
    saved_settings.parent.mkdir(parents=True)
    saved_settings.write_bytes(original)
    saved_settings.chmod(0o600)
    for rel in MCP_PATHS:
        target = backup / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(root / rel, target, symlinks=True)
        if inventory(root / rel) != inventory(target):
            raise ValueError('MCP archive integrity check failed')
    before = inventory(root)
    token = uuid.uuid4().hex
    holding = root / ('.direct-r-discard-' + token)
    stage = settings.parent / ('.direct-r-settings-' + token)
    servers.pop('occJSDM')
    stage.write_text(json.dumps(config, indent=2) + '\n')
    stage.chmod(0o600)
    moved = []
    modes = {}
    replaced = False
    holding.mkdir(mode=0o700)
    try:
        for rel in MCP_PATHS:
            source = root / rel
            target = holding / source.name
            modes[source] = source.stat().st_mode
            source.chmod(modes[source] | 0o700)
            os.rename(source, target)
            moved.append((source, target))
        os.replace(stage, settings)
        replaced = True
        check_pins(root, pins)
        if 'occJSDM' in json.loads(settings.read_text()).get('mcpServers', {}):
            raise ValueError('MCP registration remains')
        after = inventory(root)
        for rel, value in before.items():
            if rel == SETTINGS or any(rel == p or rel.startswith(p + '/') for p in MCP_PATHS):
                continue
            if after.get(rel) != value:
                raise ValueError('Retained project material changed: ' + rel)
        if verify is not None:
            verify()
    except BaseException:
        if replaced:
            stage.write_bytes(original)
            stage.chmod(0o600)
            os.replace(stage, settings)
        for source, target in reversed(moved):
            os.rename(target, source)
        for source, mode in modes.items():
            source.chmod(mode)
        stage.unlink(missing_ok=True)
        holding.rmdir()
        raise
    try:
        delete_discarded(holding)
    except OSError as exc:
        raise RuntimeError('Do not start trials: discard cleanup incomplete. Recovery archive: '
                           + str(backup) + '; residual directory: ' + str(holding)) from exc
    return {'MCP_registration_removed': True, 'MCP_implementation_archived': True,
            'preserved_hashes': pins, 'private_temporary_archive': str(backup)}


def native_check(root, executable, expected_r, revision):
    script = r'''
root <- Sys.getenv("DIRECT_R_ROOT")
runtime <- file.path(root, "mcp/paper2agent/cloud-r-runtime")
Sys.setenv(P2A_R_PROJECT=runtime)
source(file.path(runtime, "activate.R"))
library(occJSDM)
origin <- normalizePath(find.package("occJSDM"))
stopifnot(startsWith(origin, paste0(normalizePath(runtime), "/library/")),
          identical(R.version.string, Sys.getenv("DIRECT_R_VERSION")),
          identical(as.character(packageVersion("occJSDM")), "0.1.0"))
rev <- as.character(read.dcf(file.path(origin,"DESCRIPTION"),fields="RemoteSha")[,"RemoteSha"])
stopifnot(identical(unique(rev[!is.na(rev)]), Sys.getenv("DIRECT_R_REVISION")))
fit <- readRDS(file.path(root,"trial/fit-short.rds"))
d <- occJSDM::returnConvergenceDiagnostics(fit)
cat(jsonlite::toJSON(list(R=R.version.string, package_path=origin,
    source_revision=unique(rev[!is.na(rev)]), rows=nrow(d),
    rhat_above_1_01=sum(d$rhat>1.01,na.rm=TRUE),
    ess_below_400=sum(d$ess<400,na.rm=TRUE),
    rhat_unavailable=sum(!is.finite(d$rhat)),
    ess_unavailable=sum(!is.finite(d$ess))),auto_unbox=TRUE),"\n")
'''
    env = dict(os.environ, DIRECT_R_ROOT=str(root), DIRECT_R_VERSION=expected_r,
               DIRECT_R_REVISION=revision)
    completed = subprocess.run([executable, '--vanilla', '-'], input=script,
                               text=True, capture_output=True, env=env, cwd=root,
                               timeout=180)
    if completed.returncode:
        raise RuntimeError('Native read-only verification failed; no trials should start.\n'
                           + completed.stdout + completed.stderr)
    line = next((s for s in reversed(completed.stdout.splitlines()) if s.startswith('{')), None)
    if line is None:
        raise RuntimeError('Native verification returned no report')
    return json.loads(line)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--apply-copy', required=True, choices=[PROJECT_ID])
    args = parser.parse_args()
    root = Path('/cloud/project')
    backup = Path(tempfile.mkdtemp(prefix='occjsdm-direct-r-', dir='/tmp'))
    native = {}
    def verify():
        native.update(native_check(root, '/opt/R/4.6.1/lib/R/bin/Rscript',
                                   'R version 4.6.1 (2026-06-24)', REVISION))
        expected = dict(rows=32, rhat_above_1_01=29, ess_below_400=32,
                        rhat_unavailable=0, ess_unavailable=0)
        if any(native[k] != v for k, v in expected.items()):
            raise ValueError('Native reference counts differ from the frozen Cloud fit')
    try:
        verify()
        report = prepare_direct(root, backup, PINS, verify)
        # Setup/grading material must not remain in the student project.
        script = Path(__file__).resolve()
        if script.is_relative_to(root):
            saved = backup / 'prepare-direct-r.py'
            shutil.copy2(script, saved)
            if sha(script) != sha(saved):
                raise ValueError('Setup-script archive integrity failed')
            script.unlink()
        report.update(condition='C -- Sonnet direct R',
                      project_url='https://posit.cloud/spaces/828411/content/' + args.apply_copy,
                      native_read_only_check=native,
                      status='prepared_restart_R_and_Assistant_required',
                      scored_trials='not_started; teaching review still pending',
                      target_confirmation='operator supplied; Cloud identity not autodetected')
        (backup / 'report.json').write_text(json.dumps(report, indent=2) + '\n')
        print(json.dumps(report, indent=2))
    except Exception as exc:
        raise SystemExit('Setup failed. Do not start trials. Private archive: '
                         + str(backup) + '\n' + str(exc))


if __name__ == '__main__':
    main()

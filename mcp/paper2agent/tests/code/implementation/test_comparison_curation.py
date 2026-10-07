"""Exercise isolation, preserved evidence and rollback for copied Cloud projects."""
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile

import pytest

ROOT = Path(__file__).resolve().parents[3]
HELPER = ROOT / 'reports/comparison-setup/prepare-comparison-base.py'


def module():
    assert HELPER.is_file(), 'Comparison curation helper has not been implemented'
    spec = importlib.util.spec_from_file_location('curation', HELPER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def copied_project(tmp_path):
    root = tmp_path / 'project'
    root.mkdir()
    runtime = root / 'mcp/paper2agent/cloud-r-runtime'
    for name in ('activate.R', 'renv/activate.R', 'renv.lock',
                 'bootstrap-library/renv/DESCRIPTION', 'library/package/lib.so',
                 'sandbox/package/DESCRIPTION', 'state/runtime-state'):
        path = runtime / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text('preserved runtime: ' + name)
    archive = runtime / 'renv/cellar/occJSDM/occJSDM_0.1.0.tar.gz'
    archive.parent.mkdir(parents=True)
    entries = {'occJSDM/DESCRIPTION': b'Package: occJSDM\n',
               'occJSDM/R/diagnostics.R': b'native package source',
               'occJSDM/man/runOccJSDM.Rd': b'package help',
               'occJSDM/vignettes/occJSDM.Rmd': b'lesson',
               'occJSDM/tests/answers.txt': b'UPSTREAM TEST ANSWER',
               'occJSDM/dev/notes.md': b'DEVELOPMENT NOTES'}
    with tarfile.open(archive, 'w:gz') as tar:
        for name, data in entries.items():
            item = tarfile.TarInfo(name)
            item.size = len(data)
            tar.addfile(item, io.BytesIO(data))
    fit_id = 'fit_' + 'a' * 32
    folder = root / 'artifacts/calls/fit_model-reference'
    folder.mkdir(parents=True)
    fit = folder / 'fit.rds'
    data = folder / 'input.rds'
    fit.write_bytes(b'EXISTING FIT')
    data.write_bytes(b'EXISTING INPUT')
    manifest = {'status':'complete', 'fit_id':fit_id, 'fit_path':str(fit),
                'fit_sha256':digest(fit), 'input_sha256':digest(data),
                'source_revision':'pinned', 'seed':1702,
                'runtime':{'R':'example R', 'Rscript':'Rscript'}}
    manifest_path = folder / 'manifest.json'
    manifest_path.write_text(json.dumps(manifest))
    registry = root / 'artifacts/registry'
    registry.mkdir()
    (registry / (fit_id + '.json')).write_text(json.dumps({
        'manifest_path':str(manifest_path), 'manifest_sha256':digest(manifest_path)}))
    (folder / 'stdout.log').write_text('OLD NUMERICAL ANSWERS')
    for name in ('assistant-instructions.md', 'glossary.md', 'model-trial.md'):
        p = root / 'mcp/paper2agent/class' / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text('teaching' if name != 'model-trial.md' else 'GRADER')
    for name in ('src/occJSDM_mcp.py', 'occJSDM-env/bin/python'):
        p = root / 'mcp/paper2agent' / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text('MCP runtime')
    for name in ('AGENTS.md', 'dev/grader.txt', '.git/history', '.Rhistory',
                 '.posit/assistant/old-conversation.json'):
        p = root / name
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text('GRADER / OLD ANSWER')
    settings = root / '.posit/assistant/settings.json'
    settings.write_text(json.dumps({'provider':'preserve private connection',
                                   'mcpServers':{'occJSDM':{'command':['python','server']}}}))
    return root, dict(source_sha256=digest(archive), source_revision='pinned',
                      fit_id=fit_id, fit_sha256=digest(fit), input_sha256=digest(data))


def test_stage_has_only_student_material_and_exact_preserved_artifacts(tmp_path):
    mod = module()
    root, pins = copied_project(tmp_path)
    stage = root / '.stage'
    mod.stage_project(root, stage, pins, 'base')
    assert (stage / 'R/diagnostics.R').read_bytes() == b'native package source'
    assert digest(stage / 'trial/fit-short.rds') == pins['fit_sha256']
    assert digest(stage / 'trial/data.rds') == pins['input_sha256']
    assert (stage / 'mcp/paper2agent/cloud-r-runtime/library/package/lib.so').is_file()
    assert (stage / 'mcp/paper2agent/occJSDM-env/bin/python').is_file()
    assert not (stage / 'mcp/paper2agent/cloud-r-runtime/renv/cellar').exists()
    assert not (stage / 'artifacts/calls/fit_model-reference/stdout.log').exists()
    for name in ('tests', 'dev', '.git', 'AGENTS.md', '.Rhistory',
                 'mcp/paper2agent/class/model-trial.md',
                 '.posit/assistant/old-conversation.json'):
        assert not (stage / name).exists(), name
    assert (root / 'dev/grader.txt').is_file(), 'Staging must not alter the copied project'
    assert 'private connection' in (stage / '.posit/assistant/settings.json').read_text()


def test_direct_route_excludes_mcp_implementation_and_disables_only_its_registration(tmp_path):
    mod = module()
    root, pins = copied_project(tmp_path)
    stage = root / '.stage'
    mod.stage_project(root, stage, pins, 'direct-R')
    assert not (stage / 'mcp/paper2agent/src').exists()
    assert not (stage / 'mcp/paper2agent/occJSDM-env').exists()
    settings = json.loads((stage / '.posit/assistant/settings.json').read_text())
    assert 'occJSDM' not in settings['mcpServers']
    assert settings['provider'] == 'preserve private connection'


def test_corrupt_saved_fit_stops_before_staging_or_moving_anything(tmp_path):
    mod = module()
    root, pins = copied_project(tmp_path)
    (root / 'artifacts/calls/fit_model-reference/fit.rds').write_bytes(b'CHANGED')
    with pytest.raises(ValueError, match='fit.*integrity'):
        mod.stage_project(root, root / '.stage', pins, 'base')
    assert not (root / '.stage').exists()
    assert (root / 'dev/grader.txt').is_file()


def test_publish_archives_excluded_material_outside_project_and_preserves_reference(tmp_path):
    mod = module()
    root, pins = copied_project(tmp_path)
    stage = root / '.stage'
    mod.stage_project(root, stage, pins, 'base')
    backup = tmp_path / 'private-backup'
    backup.mkdir(mode=0o700)
    mod.publish_stage(root, stage, backup)
    assert not (root / 'dev').exists()
    assert (backup / 'original/dev/grader.txt').is_file()
    assert (backup / 'original/.posit/assistant/old-conversation.json').is_file()
    assert digest(root / 'trial/fit-short.rds') == pins['fit_sha256']
    assert (root / '.posit/assistant/settings.json').is_file()


def test_publish_failure_rolls_back_original_project(tmp_path, monkeypatch):
    mod = module()
    root, pins = copied_project(tmp_path)
    stage = root / '.stage'
    mod.stage_project(root, stage, pins, 'base')
    backup = tmp_path / 'private-backup'
    backup.mkdir(mode=0o700)
    original_move = mod.os.rename
    def fail_new_r(source, target):
        if Path(source) == stage / 'R':
            raise OSError('simulated interrupted publication')
        return original_move(source, target)
    monkeypatch.setattr(mod.os, 'rename', fail_new_r)
    with pytest.raises(OSError, match='interrupted'):
        mod.publish_stage(root, stage, backup)
    assert (root / 'dev/grader.txt').read_text() == 'GRADER / OLD ANSWER'
    assert (root / '.git/history').is_file()
    assert not (root / 'trial').exists()


def test_cross_filesystem_archive_never_deletes_source_during_copy(tmp_path, monkeypatch):
    mod = module()
    root, pins = copied_project(tmp_path)
    stage = root / '.stage'
    mod.stage_project(root, stage, pins, 'base')
    backup = tmp_path / 'private-backup'
    backup.mkdir(mode=0o700)
    def unsafe_move(*args, **kwargs):
        raise AssertionError('Cross-filesystem move must not be used')
    monkeypatch.setattr(mod.shutil, 'move', unsafe_move)
    mod.publish_stage(root, stage, backup)
    assert (backup / 'original/dev/grader.txt').read_text() == 'GRADER / OLD ANSWER'
    assert digest(root / 'trial/fit-short.rds') == pins['fit_sha256']


def test_incomplete_backup_aborts_before_original_changes(tmp_path, monkeypatch):
    mod = module()
    root, pins = copied_project(tmp_path)
    stage = root / '.stage'
    mod.stage_project(root, stage, pins, 'base')
    backup = tmp_path / 'private-backup'
    backup.mkdir(mode=0o700)
    copytree = mod.shutil.copytree
    def fail_backup(source, target, *args, **kwargs):
        if Path(source) == root / 'dev':
            raise PermissionError('simulated partial archive copy')
        return copytree(source, target, *args, **kwargs)
    monkeypatch.setattr(mod.shutil, 'copytree', fail_backup)
    with pytest.raises(PermissionError, match='partial archive'):
        mod.publish_stage(root, stage, backup)
    assert (root / 'dev/grader.txt').read_text() == 'GRADER / OLD ANSWER'
    assert (root / '.git/history').is_file()
    assert not (root / 'trial').exists()


def test_runtime_links_must_resolve_in_the_retained_project(tmp_path):
    mod = module()
    root, pins = copied_project(tmp_path)
    lib = root / 'mcp/paper2agent/cloud-r-runtime/library'
    (lib / 'retained-link').symlink_to(lib / 'package', target_is_directory=True)
    stage = root / '.stage'
    mod.stage_project(root, stage, pins, 'base')
    mod.check_retained_links(root, stage)
    (stage / 'mcp/paper2agent/cloud-r-runtime/library/leak').symlink_to(root / 'dev')
    with pytest.raises(ValueError, match='excluded'):
        mod.check_retained_links(root, stage)


def test_partial_cleanup_failure_retains_verified_complete_archive(tmp_path, monkeypatch):
    mod = module()
    root, pins = copied_project(tmp_path)
    stage = root / '.stage'
    mod.stage_project(root, stage, pins, 'base')
    backup = tmp_path / 'private-backup'
    backup.mkdir(mode=0o700)
    def fail_cleanup(path):
        (path / 'dev/grader.txt').unlink()
        raise PermissionError('simulated partial cleanup')
    monkeypatch.setattr(mod.shutil, 'rmtree', fail_cleanup)
    with pytest.raises(RuntimeError, match='Do not start trials'):
        mod.publish_stage(root, stage, backup)
    assert (backup / 'original/dev/grader.txt').read_text() == 'GRADER / OLD ANSWER'
    assert digest(root / 'trial/fit-short.rds') == pins['fit_sha256']


def test_read_only_renv_sandbox_can_be_archived_without_changing_retained_modes(tmp_path):
    mod = module()
    root, pins = copied_project(tmp_path)
    sandbox = root / 'mcp/paper2agent/cloud-r-runtime/sandbox/package'
    sandbox.chmod(0o555)
    stage = root / '.stage'
    mod.stage_project(root, stage, pins, 'base')
    backup = tmp_path / 'private-backup'
    backup.mkdir(mode=0o700)
    mod.publish_stage(root, stage, backup)
    rel = 'mcp/paper2agent/cloud-r-runtime/sandbox/package'
    assert (root / rel).stat().st_mode & 0o777 == 0o555
    assert (backup / 'original' / rel / 'DESCRIPTION').is_file()
    assert not list(root.glob('.comparison-original-*'))

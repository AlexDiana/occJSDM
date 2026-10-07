"""Preserve teaching inputs and runtime while isolating a direct-R copy."""
import hashlib
import importlib.util
import json
from pathlib import Path
import pytest

HELPER = Path(__file__).resolve().parents[3] / 'reports/comparison-setup/prepare-direct-r.py'


def load():
    assert HELPER.exists(), 'Direct-R setup is not implemented'
    spec = importlib.util.spec_from_file_location('direct', HELPER)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def project(tmp_path):
    root = tmp_path / 'project'
    root.mkdir()
    files = {'trial/data.rds':b'input','trial/fit-short.rds':b'fit',
             'trial/no-traits.rds':b'no traits','R/diagnostics.R':b'science',
             'teaching/glossary.md':b'teaching',
             'mcp/paper2agent/cloud-r-runtime/library/occJSDM/library.so':b'runtime',
             'mcp/paper2agent/src/server.py':b'MCP code',
             'mcp/paper2agent/occJSDM-env/bin/python':b'MCP python'}
    for rel, data in files.items():
        p = root / rel
        p.parent.mkdir(parents=True,exist_ok=True)
        p.write_bytes(data)
    p = root / '.posit/assistant/settings.json'
    p.parent.mkdir(parents=True)
    p.write_text(json.dumps({'provider':'private connection','mcpServers':{
        'occJSDM':{'command':['python','server']},'other':{'command':['other']}}}))
    pins = {rel:hashlib.sha256(data).hexdigest() for rel,data in files.items()
            if rel.startswith('trial/') and not rel.endswith('no-traits.rds')}
    backup = tmp_path / 'private'
    backup.mkdir(mode=0o700)
    return root,backup,pins


def test_direct_copy_has_no_occjsdm_mcp_and_preserves_all_other_material(tmp_path):
    mod=load();root,backup,pins=project(tmp_path)
    before={str(p.relative_to(root)):p.read_bytes() for p in root.rglob('*') if p.is_file()}
    mod.prepare_direct(root,backup,pins)
    assert not (root/'mcp/paper2agent/src').exists()
    assert not (root/'mcp/paper2agent/occJSDM-env').exists()
    config=json.loads((root/'.posit/assistant/settings.json').read_text())
    assert config=={'provider':'private connection','mcpServers':{'other':{'command':['other']}}}
    for rel,data in before.items():
        if rel.startswith(('mcp/paper2agent/src/','mcp/paper2agent/occJSDM-env/')):
            assert (backup/rel).read_bytes()==data
        elif rel!='.posit/assistant/settings.json':
            assert (root/rel).read_bytes()==data
    assert (backup/'.posit/assistant/settings.json').read_bytes()==before['.posit/assistant/settings.json']
    assert (root/'.posit/assistant/settings.json').stat().st_mode & 0o777==0o600


def test_reference_integrity_failure_leaves_copy_unchanged(tmp_path):
    mod=load();root,backup,pins=project(tmp_path)
    (root/'trial/fit-short.rds').write_bytes(b'changed')
    settings=(root/'.posit/assistant/settings.json').read_bytes()
    with pytest.raises(ValueError,match='integrity'):
        mod.prepare_direct(root,backup,pins)
    assert (root/'mcp/paper2agent/src/server.py').is_file()
    assert (root/'.posit/assistant/settings.json').read_bytes()==settings


def test_failed_publication_restores_settings_and_mcp_directories(tmp_path,monkeypatch):
    mod=load();root,backup,pins=project(tmp_path)
    before=(root/'.posit/assistant/settings.json').read_bytes()
    rename=mod.os.rename
    def fail_venv(source,target):
        if Path(source)==root/'mcp/paper2agent/occJSDM-env':
            raise OSError('simulated interrupted setup')
        return rename(source,target)
    monkeypatch.setattr(mod.os,'rename',fail_venv)
    with pytest.raises(OSError,match='interrupted'):
        mod.prepare_direct(root,backup,pins)
    assert (root/'.posit/assistant/settings.json').read_bytes()==before
    assert (root/'mcp/paper2agent/src/server.py').is_file()
    assert (root/'mcp/paper2agent/occJSDM-env/bin/python').is_file()


def test_archive_copy_failure_precedes_any_project_mutation(tmp_path,monkeypatch):
    mod=load();root,backup,pins=project(tmp_path)
    settings=(root/'.posit/assistant/settings.json').read_bytes()
    def fail(*args,**kwargs): raise OSError('archive failed')
    monkeypatch.setattr(mod.shutil,'copytree',fail)
    with pytest.raises(OSError,match='archive failed'):
        mod.prepare_direct(root,backup,pins)
    assert (root/'.posit/assistant/settings.json').read_bytes()==settings
    assert (root/'mcp/paper2agent/src/server.py').is_file()


def test_readonly_discarded_dirs_do_not_change_external_link_targets(tmp_path):
    mod=load();root,backup,pins=project(tmp_path)
    external=tmp_path/'external';external.mkdir();external.chmod(0o555)
    server=root/'mcp/paper2agent/src';(server/'external').symlink_to(external)
    server.chmod(0o555)
    mod.prepare_direct(root,backup,pins)
    assert external.stat().st_mode & 0o777==0o555
    assert not list(root.glob('.direct-r-discard-*'))


def test_postpublication_verification_failure_restores_original_route(tmp_path):
    mod=load();root,backup,pins=project(tmp_path)
    settings=(root/'.posit/assistant/settings.json').read_bytes()
    def fail_native():
        assert not (root/'mcp/paper2agent/src').exists()
        raise ValueError('native runtime mismatch')
    with pytest.raises(ValueError,match='runtime mismatch'):
        mod.prepare_direct(root,backup,pins,verify=fail_native)
    assert (root/'.posit/assistant/settings.json').read_bytes()==settings
    assert (root/'mcp/paper2agent/src/server.py').is_file()
    assert (root/'mcp/paper2agent/occJSDM-env/bin/python').is_file()
    assert not list(root.glob('.direct-r-discard-*'))

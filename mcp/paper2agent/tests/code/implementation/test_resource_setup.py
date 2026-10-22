import importlib.util
import json
from pathlib import Path

SETUP = Path(__file__).resolve().parents[3] / "reports/resource-benchmark/setup-resource-benchmark.py"


def load_setup():
    spec = importlib.util.spec_from_file_location("resource_setup", SETUP)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_install_and_remove_preserve_other_settings_and_prior_override(tmp_path):
    setup = load_setup()
    settings = tmp_path / ".posit/assistant/settings.json"
    settings.parent.mkdir(parents=True)
    config = {"theme": "test", "mcpServers": {
        "other": {"command": ["other"]},
        "occJSDM": {"command": ["python", "server"], "timeout": 30000,
                    "environment": {"P2A_RSCRIPT": "/old/Rscript", "keep": "value"}}
    }}
    settings.write_text(json.dumps(config))
    original = settings.read_bytes()
    setup.install(tmp_path)
    updated = json.loads(settings.read_text())
    launcher = tmp_path / "mcp/paper2agent/measure-rscript.py"
    assert launcher.read_text() == setup.LAUNCHER_SOURCE
    assert launcher.stat().st_mode & 0o111
    assert updated["mcpServers"]["occJSDM"]["environment"]["P2A_RSCRIPT"] == str(launcher)
    assert updated["mcpServers"]["other"] == config["mcpServers"]["other"]
    assert updated["theme"] == "test"
    assert updated["mcpServers"]["occJSDM"]["timeout"] == 30000
    assert list(settings.parent.glob("settings.json.benchmark-*.bak"))[0].read_bytes() == original
    # Removal must merge into current settings, preserving edits made meanwhile.
    updated["theme"] = "changed meanwhile"
    settings.write_text(json.dumps(updated))
    setup.remove(tmp_path)
    restored = json.loads(settings.read_text())
    assert restored["theme"] == "changed meanwhile"
    assert restored["mcpServers"]["occJSDM"]["environment"] == config["mcpServers"]["occJSDM"]["environment"]


def test_missing_server_does_not_mutate_settings(tmp_path):
    setup = load_setup()
    settings = tmp_path / ".posit/assistant/settings.json"
    settings.parent.mkdir(parents=True)
    settings.write_text('{"theme":"keep"}')
    original = settings.read_bytes()
    try:
        setup.install(tmp_path)
    except ValueError:
        pass
    else:
        raise AssertionError("Expected missing-server failure")
    assert settings.read_bytes() == original
    assert not (tmp_path / "mcp").exists()

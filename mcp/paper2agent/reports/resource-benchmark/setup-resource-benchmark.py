#!/usr/bin/env python3
"""Run in Cloud Terminal; --remove restores only the previous R launcher setting."""
import json
import os
from pathlib import Path
import stat
import sys
import uuid

LAUNCHER_SOURCE = '#!/usr/bin/env python3\n"""Temporary P2A_RSCRIPT observer; native code and arguments pass through."""\nimport json\nimport os\nfrom pathlib import Path\nimport resource\nimport shutil\nimport signal\nimport subprocess\nimport sys\nimport time\n\n\ndef main():\n    args = sys.argv[1:]\n    if len(args) != 5 or args[0] != "--vanilla":\n        raise ValueError("Expected the MCP dispatcher\'s five Rscript arguments")\n    if sys.platform not in ("linux", "darwin"):\n        raise RuntimeError("Peak RSS units are verified only for Linux and macOS")\n    native = shutil.which("Rscript")\n    if native is None or Path(native).resolve() == Path(__file__).resolve():\n        raise RuntimeError("A separate native Rscript executable is required on PATH")\n    report = Path(args[-1]).parent / "resource-usage.json"\n    if report.exists():\n        raise FileExistsError(f"Preserving existing resource report: {report}")\n    # The MCP sends termination to the whole process group, including native R.\n    # Let this observer reap R and record its exit; it stays in that same group.\n    termination_signal = None\n\n    def remember_termination(signum, frame):\n        nonlocal termination_signal\n        termination_signal = signum\n\n    signal.signal(signal.SIGTERM, remember_termination)\n    signal.signal(signal.SIGINT, remember_termination)\n    started = time.monotonic()\n    process = subprocess.run([native, *args], check=False)\n    elapsed = time.monotonic() - started\n    usage = resource.getrusage(resource.RUSAGE_CHILDREN)\n    payload = {\n        "native_executable": native,\n        "native_returncode": process.returncode,\n        "termination_signal": termination_signal,\n        "elapsed_seconds": elapsed,\n        "user_cpu_seconds": usage.ru_utime,\n        "system_cpu_seconds": usage.ru_stime,\n        "peak_individual_child_rss_bytes": usage.ru_maxrss * (1024 if sys.platform == "linux" else 1),\n        "platform": sys.platform,\n        "cpu_scope": "Native R dispatcher, including startup and waited descendants; excludes MCP and model latency",\n        "memory_scope": "Largest individual waited child peak RSS; not simultaneous process-tree or Cloud-project memory",\n    }\n    with report.open("x") as stream:\n        json.dump(payload, stream, indent=2, allow_nan=False)\n        stream.write("\\n")\n    # run_r creates this private group. Finish cleanup even if a descendant\n    # ignored TERM after native R exited; never kill a shared caller group.\n    if termination_signal is not None and os.getpgrp() == os.getpid():\n        os.killpg(os.getpgrp(), signal.SIGKILL)\n    return process.returncode if process.returncode >= 0 else 128 - process.returncode\n\n\nif __name__ == "__main__":\n    sys.exit(main())\n'


def configuration(root):
    settings = root / ".posit/assistant/settings.json"
    config = json.loads(settings.read_text())
    try:
        server = config["mcpServers"]["occJSDM"]
    except (KeyError, TypeError):
        raise ValueError("The existing occJSDM MCP configuration was not found") from None
    environment = server.setdefault("environment", {})
    if not isinstance(environment, dict):
        raise ValueError("Expected an environment object in the existing configuration")
    return settings, config, environment


def save_settings(settings, config):
    mode = stat.S_IMODE(settings.stat().st_mode)
    backup = settings.with_name(settings.name + ".benchmark-" + uuid.uuid4().hex + ".bak")
    with backup.open("xb") as stream:
        stream.write(settings.read_bytes())
    backup.chmod(mode)
    pending = settings.with_name(settings.name + "." + uuid.uuid4().hex + ".pending")
    with pending.open("x") as stream:
        json.dump(config, stream, indent=2)
        stream.write("\n")
    pending.chmod(mode)
    os.replace(pending, settings)
    print("Settings backup:", backup)


def install(root):
    settings, config, environment = configuration(root)
    folder = root / "mcp/paper2agent"
    launcher = folder / "measure-rscript.py"
    state = folder / "resource-benchmark-setup.json"
    if state.exists():
        raise FileExistsError("Benchmark setup already exists; use --remove before a new setup")
    if launcher.exists() and launcher.read_text() != LAUNCHER_SOURCE:
        raise FileExistsError(f"Preserving a different existing file: {launcher}")
    previous = {"had_override": "P2A_RSCRIPT" in environment,
                "previous_override": environment.get("P2A_RSCRIPT"),
                "installed_override": str(launcher)}
    folder.mkdir(parents=True, exist_ok=True)
    if not launcher.exists():
        with launcher.open("x") as stream:
            stream.write(LAUNCHER_SOURCE)
    launcher.chmod(0o755)
    with state.open("x") as stream:
        json.dump(previous, stream, indent=2)
        stream.write("\n")
    environment["P2A_RSCRIPT"] = str(launcher)
    save_settings(settings, config)
    print("Resource measurement ready. Restart occJSDM through /mcp in Posit Assistant.")


def remove(root):
    settings, config, environment = configuration(root)
    state = root / "mcp/paper2agent/resource-benchmark-setup.json"
    previous = json.loads(state.read_text())
    if environment.get("P2A_RSCRIPT") != previous["installed_override"]:
        raise ValueError("The launcher setting changed meanwhile; preserving it")
    if previous["had_override"]:
        environment["P2A_RSCRIPT"] = previous["previous_override"]
    else:
        environment.pop("P2A_RSCRIPT", None)
    save_settings(settings, config)
    state.rename(state.with_name("resource-benchmark-removed-" + uuid.uuid4().hex + ".json"))
    print("Original launcher setting restored. Restart occJSDM through /mcp.")


if __name__ == "__main__":
    if sys.argv[1:] not in ([], ["--remove"]):
        raise SystemExit("Usage: python3 setup-resource-benchmark.py [--remove]")
    (remove if sys.argv[1:] else install)(Path("/cloud/project"))

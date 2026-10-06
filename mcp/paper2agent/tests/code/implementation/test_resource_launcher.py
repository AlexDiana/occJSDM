"""Exercise the optional observer as a real subprocess, without native refits."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

LAUNCHER = Path(__file__).resolve().parents[3] / "reports/resource-benchmark/measure_rscript.py"


def launch(tmp_path, exit_code=0):
    native = tmp_path / "Rscript"
    native.write_text(
        f"#!{sys.executable}\n"
        "import json, sys\n"
        "from pathlib import Path\n"
        "print('native stdout')\n"
        "print('native stderr', file=sys.stderr)\n"
        "Path(sys.argv[-1]).write_text(json.dumps(sys.argv[1:]))\n"
        f"sys.exit({exit_code})\n"
    )
    native.chmod(0o755)
    args = ["--vanilla", "dispatcher with spaces.R", "activate.R",
            str(tmp_path / "request.json"), str(tmp_path / "result.json")]
    env = dict(os.environ, PATH=str(tmp_path) + os.pathsep + os.environ["PATH"])
    result = subprocess.run([sys.executable, str(LAUNCHER), *args], env=env,
                            capture_output=True, text=True, timeout=10)
    return result, args


def test_success_forwards_arguments_and_logs_and_records_resources(tmp_path):
    result, args = launch(tmp_path)
    assert result.returncode == 0, result.stderr
    assert result.stdout == "native stdout\n"
    assert result.stderr == "native stderr\n"
    assert json.loads((tmp_path / "result.json").read_text()) == args
    report = json.loads((tmp_path / "resource-usage.json").read_text())
    assert report["native_returncode"] == 0
    assert report["elapsed_seconds"] > 0
    assert report["user_cpu_seconds"] >= 0
    assert report["system_cpu_seconds"] >= 0
    assert report["peak_individual_child_rss_bytes"] > 0
    assert report["native_executable"] == str(tmp_path / "Rscript")
    assert "simultaneous" in report["memory_scope"]


def test_native_failure_remains_failure_with_resource_report(tmp_path):
    result, _ = launch(tmp_path, exit_code=7)
    assert result.returncode == 7
    report = json.loads((tmp_path / "resource-usage.json").read_text())
    assert report["native_returncode"] == 7


def test_existing_report_is_preserved_and_native_is_not_run(tmp_path):
    report = tmp_path / "resource-usage.json"
    report.write_text("preserve me")
    result, _ = launch(tmp_path)
    assert result.returncode != 0
    assert report.read_text() == "preserve me"
    assert not (tmp_path / "result.json").exists()


def test_termination_reaps_native_and_kills_remaining_private_group(tmp_path):
    native = tmp_path / "Rscript"
    native.write_text(
        f"#!{sys.executable}\n"
        "import subprocess, sys, time\n"
        "from pathlib import Path\n"
        "child = subprocess.Popen([sys.executable, '-c', "
        "'import signal,time; signal.signal(signal.SIGTERM, signal.SIG_IGN); time.sleep(30)'])\n"
        "time.sleep(0.2)\n"
        "Path(sys.argv[-1]).with_name('ready').write_text(str(child.pid))\n"
        "time.sleep(30)\n"
    )
    native.chmod(0o755)
    env = dict(os.environ, PATH=str(tmp_path) + os.pathsep + os.environ["PATH"])
    process = subprocess.Popen(
        [sys.executable, str(LAUNCHER), "--vanilla", "dispatch.R", "activate.R",
         str(tmp_path / "request.json"), str(tmp_path / "result.json")],
        env=env, start_new_session=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    try:
        deadline = time.monotonic() + 5
        while not (tmp_path / "ready").exists():
            assert time.monotonic() < deadline
            time.sleep(0.02)
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=5)
        assert process.returncode == -signal.SIGKILL
        report = json.loads((tmp_path / "resource-usage.json").read_text())
        assert report["native_returncode"] == -signal.SIGTERM
    finally:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        process.wait(timeout=5)

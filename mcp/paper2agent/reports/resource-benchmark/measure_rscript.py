#!/usr/bin/env python3
"""Temporary P2A_RSCRIPT observer; native code and arguments pass through."""
import json
import os
from pathlib import Path
import resource
import shutil
import signal
import subprocess
import sys
import time


def main():
    args = sys.argv[1:]
    if len(args) != 5 or args[0] != "--vanilla":
        raise ValueError("Expected the MCP dispatcher's five Rscript arguments")
    if sys.platform not in ("linux", "darwin"):
        raise RuntimeError("Peak RSS units are verified only for Linux and macOS")
    native = shutil.which("Rscript")
    if native is None or Path(native).resolve() == Path(__file__).resolve():
        raise RuntimeError("A separate native Rscript executable is required on PATH")
    report = Path(args[-1]).parent / "resource-usage.json"
    if report.exists():
        raise FileExistsError(f"Preserving existing resource report: {report}")
    # The MCP sends termination to the whole process group, including native R.
    # Let this observer reap R and record its exit; it stays in that same group.
    termination_signal = None

    def remember_termination(signum, frame):
        nonlocal termination_signal
        termination_signal = signum

    signal.signal(signal.SIGTERM, remember_termination)
    signal.signal(signal.SIGINT, remember_termination)
    started = time.monotonic()
    process = subprocess.run([native, *args], check=False)
    elapsed = time.monotonic() - started
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    payload = {
        "native_executable": native,
        "native_returncode": process.returncode,
        "termination_signal": termination_signal,
        "elapsed_seconds": elapsed,
        "user_cpu_seconds": usage.ru_utime,
        "system_cpu_seconds": usage.ru_stime,
        "peak_individual_child_rss_bytes": usage.ru_maxrss * (1024 if sys.platform == "linux" else 1),
        "platform": sys.platform,
        "cpu_scope": "Native R dispatcher, including startup and waited descendants; excludes MCP and model latency",
        "memory_scope": "Largest individual waited child peak RSS; not simultaneous process-tree or Cloud-project memory",
    }
    with report.open("x") as stream:
        json.dump(payload, stream, indent=2, allow_nan=False)
        stream.write("\n")
    # run_r creates this private group. Finish cleanup even if a descendant
    # ignored TERM after native R exited; never kill a shared caller group.
    if termination_signal is not None and os.getpgrp() == os.getpid():
        os.killpg(os.getpgrp(), signal.SIGKILL)
    return process.returncode if process.returncode >= 0 else 128 - process.returncode


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""Targeted read-only checks of existing archive manifests, not numerical audits."""
import argparse
import csv
import hashlib
import json
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--bulk", type=Path, required=True)
    args = parser.parse_args()
    raw, repo, output = args.archive.resolve(), args.repo.resolve(), args.output.resolve()
    bulk = args.bulk.resolve()
    if output == raw or raw in output.parents or bulk == raw or raw in bulk.parents:
        parser.error("Outputs must be outside the raw archive")
    output.mkdir(parents=True, exist_ok=True)
    bulk.mkdir(parents=True, exist_ok=True)
    ledger, cache, sizes = [], {}, {}

    def digest(path, algorithm="md5"):
        path = path.resolve()
        key = (str(path), algorithm)
        if key not in cache:
            hasher = hashlib.new(algorithm)
            with path.open("rb") as handle:
                for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
                    hasher.update(block)
            cache[key] = hasher.hexdigest()
            sizes[str(path)] = path.stat().st_size
        return cache[key]

    def check(group, path, expected="", algorithm="md5"):
        path = Path(path)
        found = path.is_file()
        actual = digest(path, algorithm) if found and expected else ""
        status = "missing" if not found else "match" if expected and actual == expected else "mismatch" if expected else "exists_only"
        ledger.append(dict(group=group, path=str(path), algorithm=algorithm if expected else "existence",
                           expected=expected, actual=actual, status=status,
                           bytes=path.stat().st_size if found else ""))

    def manifest(group, path, base):
        for row in rows(path):
            check(group, base / row["file"], row["md5"])

    extension = raw / "lesson-4-extension-20260923"
    for name in ["source-manifest.json", "source-manifest-normal-integral-v1.json"]:
        data = json.loads((extension / name).read_text())
        for file, expected in data["files"].items():
            check("extension-frozen-source", extension / file, expected, "sha256")
    jobs = rows(extension / "manifest.csv")
    for row in jobs:
        check("extension-inputs", extension / "inputs" / (row["input"] + ".rds"), row["input_md5"])
    states = []
    samples = {}
    for row in jobs:
        base = extension / "jobs" / row["job"]
        check("extension-status-files", base / "status.json")
        status = json.loads((base / "status.json").read_text())
        states.append((status.get("ok"), status.get("diagnostic_pass")))
        check("extension-result-presence", base / "result.rds")
        if status.get("ok"):
            check("extension-score-presence", base / "score.rds")
        if status.get("ok") and row["package"] not in samples:
            samples[row["package"]] = row["job"]
            check("extension-result-sample", base / "result.rds", status["result_md5"])
    for file, expected in json.loads((raw / "lesson-4-calibration-20260925/baseline-sha256.json").read_text()).items():
        check("calibration-baseline", raw / "lesson-4-calibration-20260925/baseline" / file, expected, "sha256")
    selected = {r["key"]: r for r in rows(raw / "pr11-current-20260927/selected-manifest.csv")}
    for row in rows(raw / "pr11-current-20260927/verification-selected.csv"):
        entry = selected[row["key"]]
        check("pr11-selected-results", entry["result_file"], row["result_md5"])
        check("pr11-inputs-external", entry["input_file"], row["input_md5"])
        fit = Path(entry["result_file"].replace("-result.rds", "-fit.rds"))
        check("pr11-selected-fit-presence", fit)
        if row["key"] in ["jsdm-n0100-01", "design-qfar_K6-sites300-05"]:
            check("pr11-fit-sample", fit, row["fit_md5"])
    for row in rows(raw / "spatial-targeted-20260927/input-manifest.csv"):
        check("spatial-inputs", raw / "spatial-targeted-20260927/inputs" / (row["community"] + ".rds"), row["input_md5"])
    for row in rows(raw / "spatial-targeted-20260927/verification-selected-independent.csv"):
        check("spatial-selected-results", row["file"], row["result_md5"])
        fit = Path(row["file"].replace("-result.rds", "-fit.rds"))
        check("spatial-selected-fit-presence", fit)
        if row["key"] == "range4-rep01-binary-k020":
            check("spatial-fit-sample", fit, row["fit_md5"])
    for row in rows(repo / "dev/simstudy/spatial-targeted-recheck/results/fit-source-hashes.csv"):
        file = Path(row["file"])
        path = raw / "spatial-targeted-20260927" / file
        if str(file).startswith("dev/"):
            path = raw / "spatial-targeted-20260927/research-frozen" / file.name
        check("spatial-frozen-source-library", path, row["md5"])
    manifest("intercept-frozen-source-library", raw / "intercept-prior-20260929/library-fingerprint.csv", raw / "intercept-prior-20260929")
    manifest("convergence-frozen-source-library", repo / "dev/simstudy/convergence-flag-diagnosis/results/library-fingerprint.csv", raw / "convergence-diagnosis-20261001")
    evidence = repo / "dev/simstudy/occupancy-intercept-prior/results"
    for row in rows(evidence / "evidence-manifest.csv"):
        check("intercept-compact-evidence", evidence / row["file"], row["md5"])
        check("intercept-evidence-original", raw / "intercept-prior-20260929" / row["source"], row["md5"])
    for row in rows(raw / "intercept-prior-20260929/controls/control-provenance.csv"):
        check("intercept-control-fit-presence", raw / row["control_fit"])
    for sub in ["spatial-amplitude-prior/results", "spatial-amplitude-prior/diagnosis/results"]:
        base = repo / "dev/simstudy" / sub
        manifest("amplitude-compact-evidence", base / "artifact-md5.csv", base)
    for row in rows(repo / "dev/simstudy/spatial-amplitude-prior/results/fits.csv"):
        check("amplitude-selected-results", raw / "spatial-amplitude-20260928" / row["result_file"], row["result_md5"])
        check("amplitude-selected-fit-presence", raw / "spatial-amplitude-20260928" / row["fit_file"])
    for row in rows(repo / "dev/simstudy/spatial-amplitude-prior/diagnosis/results/all-selected.csv"):
        check("amplitude-diagnosis-presence", raw / "spatial-amplitude-20260928/diagnosis-v1" / row["file"])
    for category in ["extended", "variants/a", "variants/b", "variants/c"]:
        for row in rows(repo / "dev/simstudy/convergence-flag-diagnosis/results" / category / "provenance.csv"):
            check("convergence-diagnostic-fit-presence", raw / row["fit_file"])
            if category == "extended" and row["chain"] == "1":
                check("convergence-fit-sample", raw / row["fit_file"], row["fit_md5"])
    for path in sorted((raw / "intercept-prior-inputs").rglob("*.rds")):
        matches = [r for r in selected.values() if Path(r["input_file"]).name == path.name and path.parent.parent.name in Path(r["input_file"]).parts]
        if matches:
            check("intercept-copied-inputs", path, matches[0]["input_md5"])
        else:
            check("intercept-copied-inputs-unmapped", path)
    pairs = []

    def pair(label, first, second):
        first, second = Path(first), Path(second)
        a, b = digest(first, "sha256"), digest(second, "sha256")
        pairs.append(dict(label=label, first=str(first), second=str(second), first_bytes=first.stat().st_size,
                          second_bytes=second.stat().st_size, first_sha256=a, second_sha256=b,
                          status="identical" if a == b else "distinct"))
    first = raw / "spatial-targeted-20260927/initial-review"
    second = raw / "spatial-targeted-20260927/summary-initial"
    for path in sorted(first.iterdir()):
        if path.is_file() and (second / path.name).is_file():
            pair("spatial-initial-review-copy", path, second / path.name)
    pair("native-calibration-shipped-copy", raw / "lesson-4-native-intervals-20260925/full/calibration.rds", repo / "vignettes/teaching-data/lesson-4-calibration.rds")
    for key in ["design-qfar_K6-sites300-05"]:
        pair("pr11-initial-long-distinct", raw / "pr11-current-20260927/initial" / (key + "-fit.rds"), raw / "pr11-current-20260927/long" / (key + "-fit.rds"))
    groups = defaultdict(list)
    for row in ledger:
        groups[row["group"]].append(row)
    summary = []
    for group, items in sorted(groups.items()):
        summary.append(dict(group=group, checks=len(items), unique_paths=len({r["path"] for r in items}), **{status: sum(r["status"] == status for r in items) for status in ["match", "mismatch", "missing", "exists_only"]}))
    for filename, data in [("manifest-checks.csv", ledger), ("check-summary.csv", summary), ("duplicate-checks.csv", pairs)]:
        destination = bulk if filename == "manifest-checks.csv" else output
        with (destination / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(data[0]), lineterminator="\n")
            writer.writeheader()
            writer.writerows(data)
    findings = dict(completed_utc=datetime.now(timezone.utc).isoformat(),
                    checksum_unique_files=len(sizes), checksum_unique_logical_bytes=sum(sizes.values()),
                    raw_checksum_files=sum(raw == Path(p) or raw in Path(p).parents for p in sizes),
                    raw_checksum_logical_bytes=sum(n for p,n in sizes.items() if raw == Path(p) or raw in Path(p).parents),
                    states={str(key): value for key,value in Counter(states).items()}, extension_sample_jobs=samples,
                    issues=[r for r in ledger if r["status"] in ["missing", "mismatch"]],
                    method="Existing MD5/SHA256 manifests, full reads of selected payloads, existence-only checks for other named fits; no fit loading or numerical reanalysis")
    (output / "check-totals.json").write_text(json.dumps(findings, indent=2) + "\n")
    print(json.dumps({**findings, "issues": len(findings["issues"])}, indent=2))
    print("Manifest checks:", len(ledger), "matches:", sum(r["status"] == "match" for r in ledger),
          "existence only:", sum(r["status"] == "exists_only" for r in ledger))
    if findings["issues"]:
        sys.exit(1)


if __name__ == "__main__":
    main()

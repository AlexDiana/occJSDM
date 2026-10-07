#!/usr/bin/env python3
"""Read-only archive metadata inventory; never follows symlinks or loads fits."""
import argparse
import csv
import json
import os
import stat
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path


def write_csv(path, rows, fields):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--compact", type=Path, required=True)
    parser.add_argument("--bulk", type=Path, required=True)
    args = parser.parse_args()
    root = args.archive.resolve()
    compact, bulk = args.compact.resolve(), args.bulk.resolve()
    if not root.is_dir():
        parser.error("Archive must be an existing directory")
    if root == bulk or root in bulk.parents or root == compact or root in compact.parents:
        parser.error("Output must be outside the archive being inventoried")
    compact.mkdir(parents=True, exist_ok=True)
    bulk.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat()
    records, errors = [], []
    fields = ["path", "study", "category", "kind", "extension", "bytes",
              "allocated_bytes", "mtime_ns", "device", "inode", "symlink_target"]

    def record_walk_error(exc):
        path = Path(exc.filename) if exc.filename else root
        errors.append(dict(path=os.path.relpath(path, root), error=str(exc)))

    for base, dirs, files in os.walk(root, followlinks=False, onerror=record_walk_error):
        dirs.sort()
        names = sorted(files + [name for name in dirs if (Path(base) / name).is_symlink()])
        for name in names:
            path = Path(base) / name
            rel = path.relative_to(root)
            try:
                info = path.lstat()
                kind = "file" if stat.S_ISREG(info.st_mode) else "symlink" if stat.S_ISLNK(info.st_mode) else "other"
                study = rel.parts[0] if len(rel.parts) > 1 else "_loose"
                category = rel.parts[1] if len(rel.parts) > 2 else "_top_level"
                records.append(dict(path=str(rel), study=study, category=category,
                                    kind=kind, extension=path.suffix.lower() or "[none]",
                                    bytes=info.st_size, allocated_bytes=info.st_blocks * 512,
                                    mtime_ns=info.st_mtime_ns, device=info.st_dev, inode=info.st_ino,
                                    symlink_target=os.readlink(path) if kind == "symlink" else ""))
            except OSError as exc:
                errors.append(dict(path=str(rel), error=str(exc)))
    records.sort(key=lambda row: row["path"])
    write_csv(bulk / "files.csv", records, fields)
    regular = [row for row in records if row["kind"] == "file"]
    try:
        immediate_directories = sum(p.is_dir() and not p.is_symlink() for p in root.iterdir())
    except OSError as exc:
        record_walk_error(exc)
        immediate_directories = None
    totals = dict(archive=str(root), started_utc=started,
                  completed_utc=datetime.now(timezone.utc).isoformat(),
                  immediate_directories=immediate_directories,
                  regular_files=len(regular), symlinks=sum(row["kind"] == "symlink" for row in records),
                  logical_bytes=sum(row["bytes"] for row in regular),
                  allocated_bytes=sum(row["allocated_bytes"] for row in regular), errors=errors,
                  method="os.walk without following symlinks; lstat sizes, allocation, mtime, device and inode; no fit deserialization or payload hash")
    (compact / "totals.json").write_text(json.dumps(totals, indent=2) + "\n")
    for name, keys in [("studies", ["study"]), ("categories", ["study", "category"]), ("file-types", ["extension"])]:
        groups = defaultdict(list)
        for row in regular:
            groups[tuple(row[key] for key in keys)].append(row)
        rows = []
        for key, members in sorted(groups.items()):
            row = dict(zip(keys, key))
            row.update(files=len(members), logical_bytes=sum(r["bytes"] for r in members),
                       allocated_bytes=sum(r["allocated_bytes"] for r in members))
            rows.append(row)
        write_csv(compact / (name + ".csv"), rows, keys + ["files", "logical_bytes", "allocated_bytes"])
    write_csv(compact / "largest-files.csv", sorted(regular, key=lambda row: -row["bytes"])[:30], fields)
    print(json.dumps(totals, indent=2))
    for row in sorted(regular, key=lambda row: -row["bytes"])[:5]:
        print(row["bytes"], row["path"])
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())

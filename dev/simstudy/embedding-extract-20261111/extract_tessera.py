"""Extract Tessera embeddings at the points written by extract_alphaearth.py.

DRAFT: this script compiles but has never been run. It is a starting point for
the extract, not a tested tool. See "Real embedding extract" in
dev/superpowers/specs/2026-11-11-grouped-covariate-prior-design.md.

Uses the geotessera package, which downloads the 0.1-degree tiles it needs from
Source Cooperative into --cache-dir (roughly 80 MB per tile at this latitude).

Run with (one line), after extract_alphaearth.py:
  uv run --no-project --with geotessera --with pandas --with numpy \
    python -I extract_tessera.py --cache-dir <empty dir for downloads> --out-dir .

Outputs, in --out-dir:
  tessera_<year>.csv.gz           id, T000..T127
  tessera_<year>_provenance.json

Never mix Tessera dataset versions or variants in one analysis: each is a
separately learned feature space.
"""

import argparse
import json
import sys
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pandas as pd
from geotessera import GeoTessera


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--year", type=int, default=2024)
    p.add_argument("--cache-dir", required=True)
    p.add_argument("--out-dir", required=True)
    args = p.parse_args()
    out = Path(args.out_dir)

    points = pd.read_csv(out / "points.csv")
    client = GeoTessera(cache_dir=args.cache_dir)
    lonlat = list(zip(points["lon"].tolist(), points["lat"].tolist()))
    values = np.asarray(client.sample_embeddings_at_points(lonlat, year=args.year), dtype="float32")
    if values.shape != (len(points), 128):
        raise SystemExit(f"unexpected shape {values.shape}; expected ({len(points)}, 128)")

    frame = pd.DataFrame(values, columns=[f"T{k:03d}" for k in range(128)])
    frame.insert(0, "id", points["id"].to_numpy())
    frame.to_csv(out / f"tessera_{args.year}.csv.gz", index=False, float_format="%.6f")

    provenance = {
        "product": "Tessera embeddings", "year": args.year,
        "geotessera": version("geotessera"),
        "dataset_version": str(getattr(client, "dataset_version", "package default")),
        "dataset_variant": str(getattr(client, "dataset_variant", "package default")),
        "n_points": int(len(points)), "n_missing": int(np.isnan(values).any(axis=1).sum()),
        "python": sys.version.split()[0],
    }
    (out / f"tessera_{args.year}_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps(provenance))


if __name__ == "__main__":
    main()

"""Extract AlphaEarth Foundations embeddings at random points in a bounding box.

DRAFT: this script compiles but has never completed a run. It is a starting
point for the extract, not a tested tool. The default bounding box and year are
guesses. See "Real embedding extract" in
dev/superpowers/specs/2026-11-11-grouped-covariate-prior-design.md.

Reads the public cloud-optimised GeoTIFF copy of the annual Satellite Embedding
dataset on Source Cooperative with windowed reads, so only the bounding box is
transferred. Also reads the ESA WorldCover 2021 class at each point, so forest
points can be selected later.

Run with (one line):
  uv run --no-project --with rasterio --with pandas --with pyarrow --with numpy \
    python -I extract_alphaearth.py --raw-dir <empty dir for downloads> --out-dir .

Outputs, in --out-dir:
  points.csv                      id, lon, lat, worldcover, forest
  alphaearth_<year>.csv.gz        id, A00..A63 at the point's own 10 m pixel
  alphaearth_<year>_mean3.csv.gz  id, A00..A63 averaged over the 3 x 3 pixels
  alphaearth_<year>_provenance.json

Attribution required by the licence (CC-BY 4.0): "The AlphaEarth Foundations
Satellite Embedding dataset is produced by Google and Google DeepMind."
"""

import argparse
import json
import math
import os
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.warp import transform as warp_transform
from rasterio.windows import Window

BASE = "https://data.source.coop/tge-labs/aef/v1/annual"
WORLDCOVER = ("https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map/"
              "ESA_WorldCover_10m_2021_v200_{tile}_Map.tif")
NODATA = -128


def parse_args():
    p = argparse.ArgumentParser()
    # West and north of Lake Hurdal. Claude's guess at the study landscape,
    # 11 November 2026; replace with the real extent of the plots.
    p.add_argument("--bbox", nargs=4, type=float, default=[10.80, 60.30, 11.10, 60.55],
                   metavar=("WEST", "SOUTH", "EAST", "NORTH"))
    p.add_argument("--year", type=int, default=2024)
    p.add_argument("--n-candidates", type=int, default=6000)
    p.add_argument("--seed", type=int, default=20261111)
    p.add_argument("--raw-dir", required=True, help="empty directory for downloaded index files")
    p.add_argument("--out-dir", required=True)
    return p.parse_args()


def fetch(url, dest):
    if not dest.exists():
        print(f"downloading {url}", flush=True)
        urllib.request.urlretrieve(url, dest)
    return dest


def find_path_column(index):
    """The index's file column is not documented by name; find it by content."""
    for col in index.columns:
        values = index[col].dropna()
        if len(values) and values.astype(str).str.endswith((".tiff", ".tif")).all():
            return col
    raise SystemExit(f"no file-path column found in the index; columns: {list(index.columns)}")


def to_https(path, year, zone):
    """Index paths may be gs:// or s3:// URLs or bare names; keep only the file name."""
    return f"{BASE}/{year}/{zone}/{Path(str(path)).name}"


def open_raster(tiff_url):
    """The GeoTIFFs are stored bottom-up; the provider says to read through the VRT."""
    vrt_url = tiff_url.rsplit(".", 1)[0] + ".vrt"
    try:
        return rasterio.open("/vsicurl/" + vrt_url), vrt_url
    except Exception as err:  # noqa: BLE001
        print(f"  VRT failed ({err}); reading the GeoTIFF directly", flush=True)
        return rasterio.open("/vsicurl/" + tiff_url), tiff_url


def dequantise(values):
    v = values.astype("float32")
    return ((v / 127.5) ** 2) * np.sign(v)


def main():
    args = parse_args()
    west, south, east, north = args.bbox
    raw = Path(args.raw_dir)
    out = Path(args.out_dir)
    raw.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)

    os.environ.setdefault("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
    os.environ.setdefault("CPL_VSIL_CURL_ALLOWED_EXTENSIONS", ".tif,.tiff,.vrt")
    os.environ.setdefault("GDAL_HTTP_MULTIRANGE", "YES")
    os.environ.setdefault("VSI_CACHE", "TRUE")

    rng = np.random.default_rng(args.seed)
    lon = rng.uniform(west, east, args.n_candidates)
    lat = rng.uniform(south, north, args.n_candidates)
    n = args.n_candidates

    index = pd.read_parquet(fetch(f"{BASE}/aef_index.parquet", raw / "aef_index.parquet"))
    path_col = find_path_column(index)
    hits = index[(index["year"].astype(int) == args.year)
                 & (index["wgs84_east"] > west) & (index["wgs84_west"] < east)
                 & (index["wgs84_north"] > south) & (index["wgs84_south"] < north)]
    if hits.empty:
        raise SystemExit("no AlphaEarth file intersects the bounding box for that year")
    print(f"{len(hits)} file(s) intersect the box", flush=True)

    pixel = np.full((n, 64), np.nan, dtype="float32")
    mean3 = np.full((n, 64), np.nan, dtype="float32")
    source = np.full(n, "", dtype=object)

    for _, row in hits.iterrows():
        url = to_https(row[path_col], args.year, row["utm_zone"])
        ds, opened = open_raster(url)
        print(f"reading {opened}", flush=True)
        xs, ys = warp_transform("EPSG:4326", ds.crs, lon.tolist(), lat.tolist())
        inv = ~ds.transform
        cols, rows = inv * (np.asarray(xs), np.asarray(ys))
        cols = np.floor(cols).astype(int)
        rows = np.floor(rows).astype(int)
        todo = (np.isnan(pixel[:, 0]) & (cols >= 1) & (cols < ds.width - 1)
                & (rows >= 1) & (rows < ds.height - 1))
        if not todo.any():
            ds.close()
            continue
        c0, c1 = cols[todo].min() - 1, cols[todo].max() + 2
        r0, r1 = rows[todo].min() - 1, rows[todo].max() + 2
        block = ds.read(window=Window(c0, r0, c1 - c0, r1 - r0))  # (64, h, w) int8
        ds.close()
        for i in np.flatnonzero(todo):
            r, c = rows[i] - r0, cols[i] - c0
            centre = block[:, r, c]
            if (centre == NODATA).any():
                continue
            pixel[i] = dequantise(centre)
            patch = block[:, r - 1:r + 2, c - 1:c + 2].reshape(64, -1)
            valid = (patch != NODATA).all(axis=0)
            mean3[i] = dequantise(patch[:, valid]).mean(axis=1)
            source[i] = Path(url).name

    # WorldCover tiles are 3 by 3 degrees, named by their south-west corner.
    cover = np.zeros(n, dtype="int16")
    tiles = {(3 * math.floor(a / 3), 3 * math.floor(o / 3)) for a, o in zip(lat, lon)}
    for tlat, tlon in sorted(tiles):
        name = f"{'N' if tlat >= 0 else 'S'}{abs(tlat):02d}{'E' if tlon >= 0 else 'W'}{abs(tlon):03d}"
        with rasterio.open("/vsicurl/" + WORLDCOVER.format(tile=name)) as wc:
            inside = (lat >= tlat) & (lat < tlat + 3) & (lon >= tlon) & (lon < tlon + 3)
            cols, rows = ~wc.transform * (lon[inside], lat[inside])
            cols = np.floor(cols).astype(int)
            rows = np.floor(rows).astype(int)
            c0, r0 = cols.min(), rows.min()
            block = wc.read(1, window=Window(c0, r0, cols.max() - c0 + 1, rows.max() - r0 + 1))
            cover[inside] = block[rows - r0, cols - c0]

    keep = ~np.isnan(pixel[:, 0])
    ids = np.arange(1, n + 1)[keep]
    names = [f"A{k:02d}" for k in range(64)]
    pd.DataFrame({"id": ids, "lon": lon[keep], "lat": lat[keep],
                  "worldcover": cover[keep], "forest": cover[keep] == 10}
                 ).to_csv(out / "points.csv", index=False)
    for label, values in (("", pixel), ("_mean3", mean3)):
        frame = pd.DataFrame(values[keep], columns=names)
        frame.insert(0, "id", ids)
        frame.to_csv(out / f"alphaearth_{args.year}{label}.csv.gz", index=False, float_format="%.6f")

    provenance = {
        "product": "AlphaEarth Foundations Satellite Embedding V1, annual",
        "source": BASE, "year": args.year, "bbox_wsen": args.bbox, "seed": args.seed,
        "n_candidates": n, "n_kept": int(keep.sum()), "n_forest": int((cover[keep] == 10).sum()),
        "files": sorted(set(source[keep])),
        "worldcover": "ESA WorldCover 2021 v200; class 10 is tree cover",
        "rasterio": rasterio.__version__, "python": sys.version.split()[0],
        "attribution": "The AlphaEarth Foundations Satellite Embedding dataset is produced "
                       "by Google and Google DeepMind.",
    }
    (out / f"alphaearth_{args.year}_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    print(json.dumps({k: provenance[k] for k in ("n_candidates", "n_kept", "n_forest")}))


if __name__ == "__main__":
    main()

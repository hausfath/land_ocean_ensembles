"""JRA-3Q monthly 2 m temperature over land, from NSF NCAR GDEX (d640002 monthly
statistics, 1.25°; d640000 fractional land mask).

JMA's own global mean uses the 2 m temperature analysis over land and the
background (first guess) over the ocean (Kosaka et al. 2024, JMSJ 102,
doi:10.2151/jmsj.2024-004, Sect. 9.1a). The land series here is therefore the
land-fraction-weighted mean of the ANALYSIS; the BACKGROUND land mean is kept
for the analysis-increment diagnostic and sensitivity test.

Per month (1947-09 onward) it writes
  CSV:  land_an_K, land_fc_K  = sum(T·cos(lat)·f) / sum(cos(lat)·f), all land
        incl. Antarctica and Greenland (f = fractional land mask)
  NetCDF: the same quantities as 5°x5° box means on the CRUTEM5/ERA5-5° grid
        (lat -87.5..87.5, lon -177.5..177.5) plus each box's cos·f weight, for
        coverage-matched sensitivity tests.
Each month's two ~130 KB files are downloaded to a temp dir, reduced and deleted.
Reuses the GDEX reader in ../ipcc_gmst_ensemble/src/gmst_merge/jra3q_gdex.py.

Outputs: Land Data/JRA-3Q/jra3q_land_monthly.csv, jra3q_land_5deg_monthly.nc
Usage:   python prepare_jra3q_land.py [--through 2025-12]
"""
from __future__ import annotations

import argparse
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import xarray as xr

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent / "ipcc_gmst_ensemble" / "src"))
from gmst_merge import jra3q_gdex as J  # noqa: E402

OUT_DIR = ROOT / "Land Data" / "JRA-3Q"
OUT_CSV = OUT_DIR / "jra3q_land_monthly.csv"
OUT_NC = OUT_DIR / "jra3q_land_5deg_monthly.nc"


def box_index(lat, lon):
    lon180 = np.where(lon >= 180, lon - 360, lon)
    bi = np.clip(np.floor((lat + 90) / 5).astype(int), 0, 35)
    bj = np.clip(np.floor((lon180 + 180) / 5).astype(int), 0, 71)
    return bi, bj


def make_reducer(f, lat, lon):
    w = np.cos(np.deg2rad(lat))[:, None] * f                       # cos·f weights
    bi, bj = box_index(lat, lon)
    B = (bi[:, None] * 72 + bj[None, :]).ravel()
    wb = np.bincount(B, weights=w.ravel(), minlength=36 * 72)

    def reduce(T):
        if not np.isfinite(T).all():
            raise ValueError("non-finite JRA-3Q 2 m temperature")
        glob = float((T * w).sum() / w.sum())
        box = np.bincount(B, weights=(T * w).ravel(), minlength=36 * 72)
        with np.errstate(invalid="ignore", divide="ignore"):
            box = np.where(wb > 0, box / wb, np.nan)
        return glob, box.reshape(36, 72)
    return reduce, wb.reshape(36, 72)


def fetch(y, m, reduce, session):
    with tempfile.TemporaryDirectory() as td:
        out = {}
        for k in ("an", "fc"):
            p = Path(td) / f"{k}.nc"
            if not J.download(J.month_url(k, y, m), p, session=session):
                return None
            T, _ = J._read(p, J.STREAMS[k][1])
            out[k] = reduce(T)
    return y, m, out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--through", default="2025-12")
    a = ap.parse_args()
    ty, tm = map(int, a.through.split("-"))
    s = requests.Session()
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "land.nc"
        if not J.download(J.LAND_URL, p, session=s):
            raise RuntimeError("land mask not found")
        with xr.open_dataset(p) as ds:
            f = ds[J.LAND_VAR].values.squeeze().astype(np.float64)
            lat = ds["lat"].values.astype(np.float64)
            lon = ds["lon"].values.astype(np.float64)
    if not (f.min() >= 0 and f.max() <= 1):
        raise ValueError("land fraction outside [0, 1]")
    reduce, wbox = make_reducer(f, lat, lon)
    months, (y, m) = [], J.FIRST_MONTH
    while (y, m) <= (ty, tm):
        months.append((y, m))
        y, m = J.next_month(y, m)
    print(f"{len(months)} months {months[0]}..{months[-1]}; grid {f.shape}; "
          f"land fraction {(np.cos(np.deg2rad(lat))[:, None] * f).sum() / (np.cos(np.deg2rad(lat))[:, None] * np.ones_like(f)).sum():.3f}", flush=True)
    with ThreadPoolExecutor(8) as ex:
        res = list(ex.map(lambda ym: fetch(ym[0], ym[1], reduce, s), months))
    missing = [ym for ym, r in zip(months, res) if r is None]
    if missing:
        raise RuntimeError(f"months not published: {missing[:5]}")
    rows = [{"year": y, "month": m, "land_an_K": o["an"][0], "land_fc_K": o["fc"][0]} for y, m, o in res]
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT_CSV, index=False, float_format="%.5f")
    t = pd.to_datetime([f"{y}-{m:02d}-01" for y, m, _ in res])
    xr.Dataset(
        {"an_land": (("time", "latitude", "longitude"), np.stack([o["an"][1] for *_, o in res]).astype(np.float32)),
         "fc_land": (("time", "latitude", "longitude"), np.stack([o["fc"][1] for *_, o in res]).astype(np.float32)),
         "landweight": (("latitude", "longitude"), wbox.astype(np.float32))},
        coords={"time": t, "latitude": np.arange(-87.5, 90, 5.0), "longitude": np.arange(-177.5, 180, 5.0)},
        attrs={"source": "NSF NCAR GDEX d640002 (doi:10.5065/PH0D-MH18) + d640000 land mask",
               "method": "box mean of 1.25deg points weighted by cos(lat) x land fraction; absolute K",
               "script": "land_ocean_temps/prepare_jra3q_land.py"},
    ).to_netcdf(OUT_NC)
    print(f"wrote {OUT_CSV.name} and {OUT_NC.name} ({OUT_NC.stat().st_size / 1e6:.1f} MB)", flush=True)


if __name__ == "__main__":
    main()

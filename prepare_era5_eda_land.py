"""ERA5 EDA (10-member ensemble of data assimilations) monthly 2 m temperature,
reduced to global-mean LAND and global-mean ALL-surface series per member.

Source: Copernicus CDS `reanalysis-era5-single-levels-monthly-means`,
product_type `monthly_averaged_ensemble_members`, 2m_temperature plus the EDA
land-sea mask (lsm, fractional), native 0.5° grid.

Weighting matches prepare_era5_5deg.py: land mean = sum(T·cos(lat)·lsm) /
sum(cos(lat)·lsm) over all grid points (fractional mask, Antarctica and
Greenland included); global mean = cos(lat)-weighted mean of all points.

The CDS queues one request per dataset per user, so the whole period and the
land-sea mask go in a single request (cost 4,416 of the 20,000 limit for
1980-2025; a few GB). The raw file is cached under temp_files/era5_eda_cache/,
reduced in 12-month chunks, and deleted. The global means are a validation
against the GSAT ensemble's vendored ERA5 ensemble summary
(../ipcc_gmst_ensemble/data/summaries/ERA5 ensemble/ensemble_time_series.csv).

Output: Land Data/ERA5/era5_eda_land_global_monthly.csv
  columns: year, month, member (0-9), land_K, global_K
Credentials: ~/.cdsapirc.  Usage: python prepare_era5_eda_land.py [first_year last_year]
"""
from __future__ import annotations

import shutil
import sys
import zipfile
from pathlib import Path

import cdsapi
import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).parent
OUT = ROOT / "Land Data" / "ERA5" / "era5_eda_land_global_monthly.csv"
DATASET = "reanalysis-era5-single-levels-monthly-means"
PRODUCT = "monthly_averaged_ensemble_members"
CACHE = ROOT / "temp_files" / "era5_eda_cache"   # raw file, deleted once reduced


def retrieve(years, target: Path) -> list[Path]:
    cdsapi.Client(quiet=True).retrieve(DATASET, {
        "product_type": [PRODUCT], "variable": ["2m_temperature", "land_sea_mask"],
        "year": [str(y) for y in years], "month": [f"{m:02d}" for m in range(1, 13)],
        "time": ["00:00"], "data_format": "netcdf", "download_format": "unarchived",
    }, str(target))
    return unpack(target)


def unpack(target: Path) -> list[Path]:
    if zipfile.is_zipfile(target):
        dest = target.parent / (target.stem + "_x")
        with zipfile.ZipFile(target) as z:
            z.extractall(dest)
        target.unlink()
        return sorted(dest.glob("*.nc"))
    return [target]


def open_var(paths, name):
    for p in paths:
        ds = xr.open_dataset(p)
        if name in ds:
            return ds
        ds.close()
    raise KeyError(f"{name} not in {paths}")


def land_mask(paths):
    ds = open_var(paths, "lsm")
    lsm = ds["lsm"]
    for d in list(lsm.dims):
        if d not in ("latitude", "longitude"):
            lsm = lsm.isel({d: 0})
    out = lsm.transpose("latitude", "longitude").values.astype(np.float64), ds["latitude"].values
    ds.close()
    return out


def reduce_file(paths, lsm, lat_lsm, chunk=12) -> pd.DataFrame:
    ds = open_var(paths, "t2m")
    if not np.allclose(ds["latitude"].values, lat_lsm):
        raise ValueError("EDA grid differs from lsm grid")
    tname = "valid_time" if "valid_time" in ds.dims else "time"
    t2m = ds["t2m"].transpose(tname, "number", "latitude", "longitude")
    cw = np.cos(np.deg2rad(ds["latitude"].values))[:, None] * np.ones((1, t2m.shape[-1]))
    wl = cw * lsm
    times = pd.DatetimeIndex(ds[tname].values)
    members = ds["number"].values
    rows = []
    for i0 in range(0, len(times), chunk):
        X = t2m.isel({tname: slice(i0, i0 + chunk)}).values.astype(np.float64)
        ok = np.isfinite(X)
        land = np.where(ok, X * wl, 0).sum((2, 3)) / np.where(ok, wl, 0).sum((2, 3))
        glob = np.where(ok, X * cw, 0).sum((2, 3)) / np.where(ok, cw, 0).sum((2, 3))
        for i, tt in enumerate(times[i0:i0 + chunk]):
            for j, m in enumerate(members):
                rows.append((tt.year, tt.month, int(m), land[i, j], glob[i, j]))
        print(f"  reduced through {times[min(i0 + chunk, len(times)) - 1]:%Y-%m}", flush=True)
    ds.close()
    return pd.DataFrame(rows, columns=["year", "month", "member", "land_K", "global_K"])


def main():
    y0, y1 = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) == 3 else (1980, 2025)
    CACHE.mkdir(parents=True, exist_ok=True)
    name = f"eda_{y0}-{y1}"
    target, xdir = CACHE / f"{name}.nc", CACHE / f"{name}_x"
    if xdir.exists():
        paths = sorted(xdir.glob("*.nc"))
    elif target.exists():
        paths = unpack(target)
    else:
        part = CACHE / f"{name}.part.nc"
        print(f"submitting {name} (one request: t2m + lsm, {y1 - y0 + 1} years)", flush=True)
        paths = retrieve(range(y0, y1 + 1), part)
        if paths == [part]:
            part.rename(target)
            paths = [target]
    print(f"downloaded {name} ({sum(q.stat().st_size for q in paths) / 1e6:.0f} MB)", flush=True)
    lsm, lat = land_mask(paths)
    df = reduce_file(paths, lsm, lat)
    n = df.groupby("year").size()
    if not (n == 120).all():
        raise ValueError(f"incomplete years: {n[n != 120].to_dict()}")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    df.sort_values(["year", "month", "member"]).to_csv(OUT, index=False, float_format="%.5f")
    for q in paths:
        q.unlink()
    shutil.rmtree(xdir, ignore_errors=True)
    print(f"wrote {OUT.name}: {len(df)} rows; raw deleted", flush=True)


if __name__ == "__main__":
    main()

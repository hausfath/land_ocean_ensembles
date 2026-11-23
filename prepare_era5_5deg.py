"""Reduce ERA5 monthly 2 m temperature (0.25°) to CRUTEM5's 5°x5° grid.

Used only to estimate the year-to-year persistence of CRUTEM5's coverage error
(estimate_coverage_persistence.py). For every 5° box and month it writes

  t2m_land : cos(lat) x land-fraction weighted mean over the box (land points only)
  t2m_all  : cos(lat) weighted mean over the whole box
  landfrac : cos(lat) weighted land fraction of the box (time-invariant)

Inputs (not vendored; ~1.4 GB):
  ERA5 monthly averaged reanalysis, 2m_temperature, 0.25°, 1940-01 onward
  ERA5 land-sea mask (lsm), 0.25°
Both from the Copernicus Climate Data Store (reanalysis-era5-single-levels-monthly-means).

Output: Land Data/ERA5/era5_t2m_5deg_monthly.nc (~25 MB), absolute K.

Usage: python prepare_era5_5deg.py T2M_FILE LSM_FILE
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import xarray as xr

ROOT = Path(__file__).parent
OUT = ROOT / "Land Data" / "ERA5" / "era5_t2m_5deg_monthly.nc"


def to_boxes(field: np.ndarray, w: np.ndarray) -> np.ndarray:
    """field, w: (..., 720, 1440) on rows 90..-89.75, cols -180..179.75.
    Returns weighted box means (..., 36, 72), north-first rows."""
    sh = field.shape[:-2]
    f = field.reshape(*sh, 36, 20, 72, 20)
    ww = np.broadcast_to(w, field.shape).reshape(*sh, 36, 20, 72, 20)
    num = np.nansum(f * ww, axis=(-3, -1))
    den = np.sum(np.where(np.isfinite(f), ww, 0.0), axis=(-3, -1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, num / den, np.nan)


def main(t2m_file: str, lsm_file: str) -> None:
    lsm = xr.open_dataset(lsm_file)["lsm"].squeeze(drop=True)
    ds = xr.open_dataset(t2m_file)
    assert np.allclose(ds.latitude.values, lsm.latitude.values)
    lat = ds.latitude.values[:720]                       # drop the -90 row
    roll = 720                                            # lon 180.0 -> first column (-180)
    lon = np.roll(ds.longitude.values, -roll)
    lon = np.where(lon >= 180, lon - 360, lon)
    assert lon[0] == -180.0 and lat[0] == 90.0
    coslat = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, 1440))
    land = np.roll(lsm.values[:720], -roll, axis=1)
    w_land = coslat * land
    w_all = coslat
    landfrac = to_boxes(land, coslat)

    nt = ds.sizes["time"]
    out_land = np.full((nt, 36, 72), np.nan, dtype=np.float32)
    out_all = np.full((nt, 36, 72), np.nan, dtype=np.float32)
    for i0 in range(0, nt, 12):
        i1 = min(nt, i0 + 12)
        t = ds["t2m"].isel(time=slice(i0, i1)).values[:, :720]
        t = np.roll(t, -roll, axis=2).astype(np.float64)
        out_land[i0:i1] = to_boxes(t, w_land)
        out_all[i0:i1] = to_boxes(t, w_all)
        if i0 % 120 == 0:
            print(" ", str(ds.time.values[i0])[:7], flush=True)

    # north-first -> south-first rows to match CRUTEM5 (-87.5 .. 87.5)
    blat = np.arange(-87.5, 90, 5.0)
    blon = np.arange(-177.5, 180, 5.0)
    res = xr.Dataset(
        {
            "t2m_land": (("time", "latitude", "longitude"), out_land[:, ::-1]),
            "t2m_all": (("time", "latitude", "longitude"), out_all[:, ::-1]),
            "landfrac": (("latitude", "longitude"), landfrac[::-1].astype(np.float32)),
        },
        coords={"time": ds.time.values, "latitude": blat, "longitude": blon},
        attrs={
            "source": "ERA5 monthly means 2m_temperature 0.25deg + lsm, Copernicus CDS",
            "method": "cos(lat) (x lsm for t2m_land) weighted means over 20x20 0.25deg points per 5deg box; -90 row dropped",
            "script": "land_ocean_temps/prepare_era5_5deg.py",
        },
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    res.to_netcdf(OUT)
    print("wrote", OUT, f"{OUT.stat().st_size/1e6:.1f} MB")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

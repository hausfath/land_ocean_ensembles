"""Reduce HadISST1 monthly SST (1°) to HadSST4's 5°x5° grid.

Used only to estimate the year-to-year persistence of HadSST4's coverage error
(estimate_sst_coverage_persistence.py). HadSST4's coverage uncertainty is itself
estimated by subsampling HadISST.1.1 (HadSST.4.2.0.0 Product User Guide).

Downloads HadISST_sst.nc.gz (243 MB) from the Met Office to a temporary
directory, writes the 5° box means, and deletes the raw file.
Sea-ice cells (-1000 in HadISST) are set to -1.8 °C; land is excluded.
Each 5° box mean is the cos(lat)-weighted mean of its ocean 1° cells.

Output: Ocean Data/HadISST/HadISST1_sst_5deg_monthly.nc (~20 MB)
"""
from __future__ import annotations

import gzip
import shutil
import tempfile
from pathlib import Path

import numpy as np
import requests
import xarray as xr

ROOT = Path(__file__).parent
OUT = ROOT / "Ocean Data" / "HadISST" / "HadISST1_sst_5deg_monthly.nc"
URL = "https://www.metoffice.gov.uk/hadobs/hadisst/data/HadISST_sst.nc.gz"


def main() -> None:
    with tempfile.TemporaryDirectory() as td:
        gz = Path(td) / "HadISST_sst.nc.gz"
        nc = Path(td) / "HadISST_sst.nc"
        with requests.get(URL, stream=True, timeout=600, headers={"User-Agent": "Mozilla/5.0"}) as r:
            r.raise_for_status()
            with open(gz, "wb") as f:
                for chunk in r.iter_content(1 << 22):
                    f.write(chunk)
        with gzip.open(gz, "rb") as fi, open(nc, "wb") as fo:
            shutil.copyfileobj(fi, fo, 1 << 22)
        gz.unlink()
        ds = xr.open_dataset(nc)
        lat = ds.latitude.values            # 89.5 .. -89.5
        lon = ds.longitude.values           # -179.5 .. 179.5
        assert lat[0] > 0 and lon[0] < 0 and len(lat) == 180 and len(lon) == 360
        w = np.cos(np.deg2rad(lat))[:, None] * np.ones((1, 360))
        nt = ds.sizes["time"]
        out = np.full((nt, 36, 72), np.nan, dtype=np.float32)
        for i0 in range(0, nt, 120):
            x = ds["sst"].isel(time=slice(i0, i0 + 120)).values.astype(np.float64)
            x = np.where(x < -100, -1.8, x)                  # sea ice
            fin = np.isfinite(x)
            num = np.where(fin, x, 0.0) * w
            den = fin * w
            num = num.reshape(-1, 36, 5, 72, 5).sum(axis=(2, 4))
            den = den.reshape(-1, 36, 5, 72, 5).sum(axis=(2, 4))
            with np.errstate(invalid="ignore"):
                out[i0:i0 + 120] = np.where(den > 0, num / den, np.nan)
        times = ds.time.values
        ds.close()
    res = xr.Dataset(
        {"sst": (("time", "latitude", "longitude"), out[:, ::-1])},   # south-first like HadSST4
        coords={"time": times, "latitude": np.arange(-87.5, 90, 5.0), "longitude": np.arange(-177.5, 180, 5.0)},
        attrs={"source": URL, "script": "land_ocean_temps/prepare_hadisst_5deg.py",
               "method": "cos(lat)-weighted mean of ocean 1deg cells per 5deg box; sea ice set to -1.8C"},
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    res.to_netcdf(OUT)
    print("wrote", OUT, f"{OUT.stat().st_size / 1e6:.1f} MB")


if __name__ == "__main__":
    main()

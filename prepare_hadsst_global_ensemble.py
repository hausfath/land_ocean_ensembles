"""
Build a 200-member global-mean monthly SST anomaly ensemble from the
HadSST 4.2.0.0 gridded ensemble files.

Inputs:
    Ocean Data/HadSST4/HadSST.4.2.0.0_ensemble_members_{1-50,51-100,101-150,151-200}.zip
    (each zip contains 50 NetCDFs named HadSST.4.2.0.0_ensemble_member_N.nc)

For each member we compute the area-weighted (cos(lat)) global ocean mean
of the `tos` field per month (NaN cells, i.e. land/ice-covered, are
skipped from both numerator and denominator).

Output:
    Ocean Data/HadSST4/HadSST.4.2.0.0_global_ensemble_monthly.csv
        (n_months × 202: year, month, m001..m200)
"""
from __future__ import annotations

import shutil
import subprocess
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).parent
HADSST = ROOT / "Ocean Data" / "HadSST4"
TMP = ROOT / "temp_files" / "hadsst_members"
TMP.mkdir(parents=True, exist_ok=True)
OUT = HADSST / "HadSST.4.2.0.0_global_ensemble_monthly.csv"

ZIPS = [
    (HADSST / "HadSST.4.2.0.0_ensemble_members_1-50.zip", range(1, 51)),
    (HADSST / "HadSST.4.2.0.0_ensemble_members_51-100.zip", range(51, 101)),
    (HADSST / "HadSST.4.2.0.0_ensemble_members_101-150.zip", range(101, 151)),
    (HADSST / "HadSST.4.2.0.0_ensemble_members_151-200.zip", range(151, 201)),
]


def area_weighted_ocean_mean(da: xr.DataArray) -> np.ndarray:
    lat_name = "latitude" if "latitude" in da.dims else "lat"
    lon_name = "longitude" if "longitude" in da.dims else "lon"
    weights = np.cos(np.deg2rad(da[lat_name]))
    weights_b = xr.broadcast(weights, da)[0]
    mask = np.isfinite(da)
    w = weights_b.where(mask, 0.0)
    num = (da.where(mask, 0.0) * w).sum(dim=(lat_name, lon_name))
    den = w.sum(dim=(lat_name, lon_name))
    return (num / den).values


def main() -> None:
    rows = None
    member_cols: dict[str, np.ndarray] = {}

    t0 = time.time()
    total_done = 0
    for zip_path, member_range in ZIPS:
        print(f"\nUnzipping {zip_path.name} ...")
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(TMP)

        for m in member_range:
            nc = TMP / f"HadSST.4.2.0.0_ensemble_member_{m}.nc"
            with xr.open_dataset(nc) as ds:
                series = area_weighted_ocean_mean(ds["tos"])
                time_index = ds["time"].values

            df = pd.DataFrame({
                "year": pd.DatetimeIndex(time_index).year,
                "month": pd.DatetimeIndex(time_index).month,
            })
            if rows is None:
                rows = df.copy()
            member_cols[f"m{m:03d}"] = series

            nc.unlink()
            total_done += 1
            elapsed = time.time() - t0
            print(
                f"[{total_done:3d}/200] m{m:03d} mean={series.mean():+.3f}  "
                f"elapsed={elapsed:.1f}s"
            )

    # safety: remove any leftover files in TMP, but keep the directory
    for p in TMP.iterdir():
        if p.is_file():
            p.unlink()

    assert rows is not None
    out = pd.concat([rows, pd.DataFrame(member_cols)], axis=1)
    out.to_csv(OUT, index=False, float_format="%.6f")
    print(f"\nSaved {OUT}  ({len(out)} rows × {out.shape[1]} cols)")


if __name__ == "__main__":
    main()

"""
COBE-SST2 ships as absolute monthly SST on a 1°×1° grid; we need
anomalies vs. the 1991-2020 climatology and a global ocean area-weighted
monthly mean.

This script:
  1. Subtracts the 1991-2020 monthly climatology from each month's grid.
  2. Computes the cos(lat)-weighted global ocean mean (ocean cells are those
     where SST is finite — land/sea-ice cells are NaN).
  3. Saves a CSV of monthly anomaly time series 1850-2026.

Output:
    Ocean Data/COBE-SST2/COBE-SST2_global_monthly.csv
"""
from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).parent
COBE = ROOT / "Ocean Data" / "COBE-SST2"
OUT = COBE / "COBE-SST2_global_monthly.csv"


def main() -> None:
    print("Loading COBE-SST2 monthly mean (549 MB) ...")
    sst = xr.open_dataset(COBE / "sst.mon.mean.nc")["sst"]
    print(f"  time: {sst.time.values[[0,-1]]}  shape: {sst.shape}")

    print("Loading 1991-2020 climatology ...")
    ltm = xr.open_dataset(COBE / "sst.mon.ltm.1991-2020.nc")["sst"]
    # ltm has 12 months indexed by cftime year-0001
    # match by calendar month
    ltm_arr = ltm.values   # (12, lat, lon)

    print("Subtracting climatology and computing ocean global mean ...")
    months = sst["time"].dt.month.values
    # Compute weights once
    weights = np.cos(np.deg2rad(sst["lat"].values))[:, None]  # (lat, 1)

    monthly_anom_mean = np.empty(sst.shape[0], dtype=np.float64)
    for i in range(sst.shape[0]):
        grid = sst.isel(time=i).values - ltm_arr[months[i] - 1]
        mask = np.isfinite(grid)
        w = np.where(mask, weights, 0.0)
        num = np.where(mask, grid * w, 0.0).sum()
        den = w.sum()
        monthly_anom_mean[i] = num / den if den > 0 else np.nan
        if i % 240 == 0:
            print(f"  {sst.time.values[i]}: anom={monthly_anom_mean[i]:+.3f}")

    df = pd.DataFrame({
        "year": pd.DatetimeIndex(sst.time.values).year,
        "month": pd.DatetimeIndex(sst.time.values).month,
        "anomaly": monthly_anom_mean,
    })
    df.to_csv(OUT, index=False, float_format="%.6f")
    print(f"\nSaved {OUT} ({len(df)} rows)")


if __name__ == "__main__":
    main()

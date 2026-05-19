"""C-LSAT 2.1 (Sun / Li, Sun Yat-sen University) ships as a 5°×5° monthly
land surface air temperature anomaly NetCDF with non-CF time encoded as
integer YYYYMM. This script:

  1. Decodes the YYYYMM time axis.
  2. Computes the cos(lat)-weighted global land-area mean per month
     (land cells are those where the anomaly field is finite).
  3. Saves a CSV of monthly anomaly time series 1850-01 through 2025-12.

Native climatology of the source NetCDF is 1961-1990; we leave anomalies
on that baseline here and re-baseline to 1981-2010 inside the build.

Source: http://www.gwpu.net/en/h-col-103.html
File:   Land Data/C-LSAT/China-LSAT2.1_tavg.nc

Output:
    Land Data/C-LSAT/C-LSAT2.1_global_monthly.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).parent
CLS = ROOT / "Land Data" / "C-LSAT"
SRC = CLS / "China-LSAT2.1_tavg.nc"
OUT = CLS / "C-LSAT2.1_global_monthly.csv"


def main() -> None:
    print(f"Loading {SRC.name} ...")
    ds = xr.open_dataset(SRC)
    da = ds["tavg_anomaly"]
    tcode = da["time"].values.astype(np.int64)
    years = (tcode // 100).astype(int)
    months = (tcode % 100).astype(int)
    print(f"  shape {tuple(da.sizes.values())}  "
          f"time {int(tcode[0])} → {int(tcode[-1])}  "
          f"(n={len(tcode)} months)")

    lat = da["lat"].values
    w = np.cos(np.deg2rad(lat))[None, :, None]   # (1, lat, 1) broadcastable
    mask = np.isfinite(da.values)
    num = np.nansum(np.where(mask, da.values, 0.0) * w, axis=(1, 2))
    den = np.nansum(mask * w, axis=(1, 2))
    gm = np.where(den > 0, num / den, np.nan)
    print(f"  finite monthly means: {np.isfinite(gm).sum()} / {len(gm)}")

    df = pd.DataFrame({"year": years, "month": months, "anomaly": gm})
    df.to_csv(OUT, index=False, float_format="%.6f")
    print(f"\nSaved {OUT} ({len(df)} rows)")
    print("\nAnnual sanity (ref 1961-1990, °C):")
    ann = df.groupby("year")["anomaly"].mean()
    for y in [1850, 1900, 1950, 2000, 2024, 2025]:
        if y in ann.index:
            print(f"  {y}: {ann.loc[y]:+.3f}")


if __name__ == "__main__":
    main()

"""
COBE-SST3 best-estimate global-mean monthly anomaly.

The JMA archive ships per-year monthly 1° gridded NetCDFs of absolute
SST (`sst`, °C) along with `err`, `sic`, `da`. We:

  1. Download each year file 1850-2024 to a streaming buffer (~2 MB
     each, ~367 MB total).
  2. Concatenate along time and compute a 1991-2020 monthly climatology
     per grid cell.
  3. Subtract climatology and area-weight (`cos(lat)`) over ocean cells
     (where `sst` is finite — SST3 distributes finite SST under sea
     ice with the analysis-interpolated value, so ice cells *are*
     included; this matches the COBE-SST2 / DCENT-I convention of
     including ice cells in the SST mean).
  4. Save monthly anomaly time series 1850-2024.
  5. Cross-check against JMA's pre-computed `gm/gm_cobe-sst3` by
     1991-2020-anchoring both and reporting RMSE.

Output:
    Ocean Data/COBE-SST3/COBE-SST3_global_monthly.csv

The downstream `build_ocean_ensemble.py` re-baselines this to 1981-2010
member-wise; the 1991-2020 internal climatology cancels out then.
"""
from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).parent
DATA = ROOT / "Ocean Data" / "COBE-SST3"
TMP = ROOT / "temp_files" / "cobe_sst3_central"
DATA.mkdir(parents=True, exist_ok=True)
TMP.mkdir(parents=True, exist_ok=True)

OUT_CSV = DATA / "COBE-SST3_global_monthly.csv"
GM_REF = DATA / "gm_cobe-sst3"  # JMA pre-computed (downloaded for cross-check)

URL_FMT = (
    "https://climate.mri-jma.go.jp/pub/archives/Ishii-et-al_COBE-SST3/"
    "cobe-sst3/1x1/monthly/cobe-sst3.monthly.{year}.nc"
)

YEAR_START = 1850
YEAR_END = 2024
CLIM_START, CLIM_END = 1991, 2020


def download(url: str, dest: Path, attempts: int = 4) -> None:
    for k in range(attempts):
        try:
            # --connect-timeout / --max-time guard against stale TCP routes
            # after a system sleep (otherwise curl hangs indefinitely on a
            # half-open connection). 300 s is well above typical transfer
            # times for any file we pull.
            subprocess.run(
                ["curl", "-L", "--fail", "-sS",
                 "--connect-timeout", "30", "--max-time", "300",
                 "-o", str(dest), url],
                check=True,
            )
            if dest.stat().st_size < 500_000:
                raise RuntimeError(f"unexpectedly small file ({dest.stat().st_size} B)")
            return
        except (subprocess.CalledProcessError, RuntimeError) as exc:
            wait = 2 ** k
            print(f"  retry {k+1}/{attempts} after {wait}s: {exc}")
            if dest.exists():
                dest.unlink()
            time.sleep(wait)
    raise RuntimeError(f"failed to download {url} after {attempts} attempts")


def load_year(year: int, keep_local: bool) -> xr.DataArray:
    path = TMP / f"cobe-sst3.monthly.{year}.nc"
    if not path.exists():
        download(URL_FMT.format(year=year), path)
    ds = xr.open_dataset(path)
    sst = ds["sst"].load()
    ds.close()
    if not keep_local:
        path.unlink()
    return sst


def compute_global_mean(sst: xr.DataArray, climatology: np.ndarray) -> np.ndarray:
    """Return (12,) monthly global ocean-mean anomalies for a single year.

    `sst` has shape (12, lat, lon); `climatology` shape (12, lat, lon).
    Ocean cells are wherever `sst` is finite (SST3 includes ice cells
    with the analysis value; land is NaN).
    """
    lats = sst["latitude"].values
    weights = np.cos(np.deg2rad(lats))[:, None]  # (lat, 1)
    out = np.empty(12, dtype=np.float64)
    for m in range(12):
        grid = sst.isel(time=m).values - climatology[m]
        mask = np.isfinite(grid)
        w = np.where(mask, weights, 0.0)
        num = np.where(mask, grid * w, 0.0).sum()
        den = w.sum()
        out[m] = num / den if den > 0 else np.nan
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep-local", action="store_true",
                    help="retain yearly NC files in temp_files/ after processing")
    args = ap.parse_args()

    print(f"Loading climatology window {CLIM_START}-{CLIM_END} ...")
    # Pass 1: build the climatology from CLIM_START..CLIM_END
    clim_sum = None
    clim_cnt = None
    for y in range(CLIM_START, CLIM_END + 1):
        sst = load_year(y, keep_local=True)  # always keep clim years; we re-open below if streaming
        arr = sst.values  # (12, lat, lon)
        if clim_sum is None:
            clim_sum = np.where(np.isfinite(arr), arr, 0.0)
            clim_cnt = np.isfinite(arr).astype(np.int32)
        else:
            clim_sum += np.where(np.isfinite(arr), arr, 0.0)
            clim_cnt += np.isfinite(arr).astype(np.int32)
        print(f"  clim year {y} ingested")
    climatology = np.where(clim_cnt > 0, clim_sum / np.maximum(clim_cnt, 1), np.nan)
    print(f"  climatology shape: {climatology.shape}  finite cells (Jan): "
          f"{np.isfinite(climatology[0]).sum()}")

    # Pass 2: compute anomalies year-by-year
    print(f"\nComputing monthly anomalies {YEAR_START}-{YEAR_END} ...")
    years, months, anoms = [], [], []
    for y in range(YEAR_START, YEAR_END + 1):
        # within climatology window we already have the file in tmp
        in_clim_window = CLIM_START <= y <= CLIM_END
        sst = load_year(y, keep_local=(in_clim_window or args.keep_local))
        monthly_anom = compute_global_mean(sst, climatology)
        years.extend([y] * 12)
        months.extend(list(range(1, 13)))
        anoms.extend(monthly_anom.tolist())
        if y % 20 == 0 or y in (YEAR_START, YEAR_END):
            print(f"  {y}: jan={monthly_anom[0]:+.3f}  jul={monthly_anom[6]:+.3f}  "
                  f"annual={np.mean(monthly_anom):+.3f}")

    df = pd.DataFrame({"year": years, "month": months, "anomaly": anoms})
    df.to_csv(OUT_CSV, index=False, float_format="%.6f")
    print(f"\nSaved {OUT_CSV} ({len(df)} rows)")

    # cross-check vs JMA gm_cobe-sst3 if it's been downloaded
    if GM_REF.exists():
        print("\nCross-checking vs JMA gm_cobe-sst3 ...")
        ref = pd.read_csv(GM_REF, sep=r"\s+", header=None,
                          names=["year", "month", "value", "uncertainty"])
        ref = ref[(ref["year"] >= YEAR_START) & (ref["year"] <= YEAR_END)]
        merged = df.merge(ref, on=["year", "month"], how="inner")
        # both should be re-zeroed on 1991-2020 mean before comparison
        ours = merged["anomaly"].values
        theirs = merged["value"].values
        clim_mask = (merged["year"] >= CLIM_START) & (merged["year"] <= CLIM_END)
        ours -= ours[clim_mask].mean()
        theirs -= theirs[clim_mask].mean()
        diff = ours - theirs
        print(f"  {len(merged)} months overlapping")
        print(f"  RMSE: {np.sqrt(np.mean(diff ** 2)):.4f} °C")
        print(f"  max abs diff: {np.max(np.abs(diff)):.4f} °C")
        print(f"  mean bias (ours - theirs): {diff.mean():+.4f} °C")


if __name__ == "__main__":
    main()

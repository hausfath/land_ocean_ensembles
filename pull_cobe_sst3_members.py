"""
COBE-SST3 perturbation ensemble (Ishii et al. 2025; doi:10.2151/jmsj.2024-010).

The MRI/JMA archive distributes the SST3 perturbation ensemble as
per-member-per-year monthly NetCDF files on a 1°×1° grid at
`https://climate.mri-jma.go.jp/pub/archives/Ishii-et-al_COBE-SST3/
 cobe-sst3/1x1/perturb/m{NNN}/cobe-sst3.monthly.m{NNN}.{YYYY}.nc`.

As verified 2026-05-22, every sampled member (m001, m099, m150, m300)
carries 155 yearly files spanning **1870–2024** (the `00README`'s claim
of "1890–2020, post-2020 by request" appears to be outdated relative
to the 2025-10-01 archive contents). We adopt 1870–2024 as the native
perturbation window.

Algorithm — streaming per-member, peak transient disk ~155 MB:

  For each member m:
    1. Download all 155 yearly files into TMP/m{NNN}/.
    2. Compute the per-cell 1991–2020 monthly climatology from this
       member's own gridded series.
    3. Subtract climatology and compute the cos(lat)-weighted global
       ocean mean per month (mask = `sst` is finite, matching the
       SST3 sea-ice-cells-included convention).
    4. Append the 1860-month time series to the in-memory ensemble.
    5. Delete the member's TMP directory.

Output:
    Ocean Data/COBE-SST3/COBE-SST3_monthly_ensemble.csv
        rows: 1860 (Jan 1870 .. Dec 2024)
        cols: year, month, m001..m300

Total transfer: 300 members × ~155 MB = ~46.5 GB. Wall time depends on
the link to MRI; expect 2–4 h on a typical home connection. Re-running
with `--first M --last N` resumes a partial population.
"""
from __future__ import annotations

import argparse
import shutil
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).parent
OUT = ROOT / "Ocean Data" / "COBE-SST3" / "COBE-SST3_monthly_ensemble.csv"
TMP = ROOT / "temp_files" / "cobe_sst3_members"
TMP.mkdir(parents=True, exist_ok=True)
OUT.parent.mkdir(parents=True, exist_ok=True)

URL_FMT = (
    "https://climate.mri-jma.go.jp/pub/archives/Ishii-et-al_COBE-SST3/"
    "cobe-sst3/1x1/perturb/m{member:03d}/cobe-sst3.monthly.m{member:03d}.{year}.nc"
)

YEAR_START = 1870
YEAR_END = 2024
N_YEARS = YEAR_END - YEAR_START + 1
N_MEMBERS_TOTAL = 300
CLIM_START, CLIM_END = 1991, 2020


def download(url: str, dest: Path, attempts: int = 4) -> None:
    for k in range(attempts):
        try:
            subprocess.run(
                ["curl", "-L", "--fail", "-sS", "-o", str(dest), url],
                check=True,
            )
            if dest.stat().st_size < 500_000:
                raise RuntimeError(f"unexpectedly small file ({dest.stat().st_size} B)")
            return
        except (subprocess.CalledProcessError, RuntimeError) as exc:
            wait = 2 ** k
            print(f"    retry {k+1}/{attempts} after {wait}s: {exc}")
            if dest.exists():
                dest.unlink()
            time.sleep(wait)
    raise RuntimeError(f"failed to download {url} after {attempts} attempts")


def fetch_member_yearfiles(member: int, member_dir: Path) -> None:
    member_dir.mkdir(parents=True, exist_ok=True)
    for year in range(YEAR_START, YEAR_END + 1):
        path = member_dir / f"cobe-sst3.monthly.m{member:03d}.{year}.nc"
        if path.exists() and path.stat().st_size > 500_000:
            continue
        download(URL_FMT.format(member=member, year=year), path)


def member_global_monthly(member_dir: Path, member: int) -> np.ndarray:
    """Returns (N_YEARS*12,) monthly global ocean anomaly time series
    1870-01 .. 2024-12. Climatology is this member's own 1991–2020 mean."""
    # Pass 1: per-cell 1991-2020 climatology from this member's files
    clim_sum = None
    clim_cnt = None
    weights = None
    for y in range(CLIM_START, CLIM_END + 1):
        path = member_dir / f"cobe-sst3.monthly.m{member:03d}.{y}.nc"
        with xr.open_dataset(path) as ds:
            arr = ds["sst"].values  # (12, lat, lon)
            if clim_sum is None:
                clim_sum = np.zeros_like(arr, dtype=np.float64)
                clim_cnt = np.zeros_like(arr, dtype=np.int32)
                lats = ds["latitude"].values
                weights = np.cos(np.deg2rad(lats))[:, None].astype(np.float64)
        finite = np.isfinite(arr)
        clim_sum += np.where(finite, arr, 0.0)
        clim_cnt += finite.astype(np.int32)
    clim = np.where(clim_cnt > 0, clim_sum / np.maximum(clim_cnt, 1), np.nan)

    # Pass 2: monthly global mean anomaly
    out = np.empty(N_YEARS * 12, dtype=np.float64)
    for iy, y in enumerate(range(YEAR_START, YEAR_END + 1)):
        path = member_dir / f"cobe-sst3.monthly.m{member:03d}.{y}.nc"
        with xr.open_dataset(path) as ds:
            arr = ds["sst"].values  # (12, lat, lon)
        for m in range(12):
            grid = arr[m] - clim[m]
            mask = np.isfinite(grid)
            w = np.where(mask, weights, 0.0)
            denom = w.sum()
            num = np.where(mask, grid * w, 0.0).sum()
            out[iy * 12 + m] = num / denom if denom > 0 else np.nan
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--first", type=int, default=1)
    ap.add_argument("--last", type=int, default=N_MEMBERS_TOTAL)
    ap.add_argument("--out-suffix", type=str, default=None,
                    help="suffix appended to the output CSV (for partial runs)")
    ap.add_argument("--keep-member-files", action="store_true",
                    help="retain per-member yearly files after processing")
    args = ap.parse_args()
    assert 1 <= args.first <= args.last <= N_MEMBERS_TOTAL

    years = np.repeat(np.arange(YEAR_START, YEAR_END + 1), 12)
    months = np.tile(np.arange(1, 13), N_YEARS)

    member_cols: dict[str, np.ndarray] = {}
    t0 = time.time()
    n_members = args.last - args.first + 1
    for k, member in enumerate(range(args.first, args.last + 1), start=1):
        member_dir = TMP / f"m{member:03d}"
        print(f"\n[{k:3d}/{n_members}] member m{member:03d}: downloading "
              f"{N_YEARS} year files ...")
        fetch_member_yearfiles(member, member_dir)

        series = member_global_monthly(member_dir, member)
        member_cols[f"m{member:03d}"] = series

        if not args.keep_member_files:
            shutil.rmtree(member_dir)

        elapsed = time.time() - t0
        rate = k / elapsed if elapsed > 0 else 0.0
        eta = (n_members - k) / rate if rate > 0 else float("nan")
        ann_2024 = np.nanmean(series[years == 2024])
        ann_1900 = np.nanmean(series[years == 1900])
        print(
            f"    m{member:03d}: 1900={ann_1900:+.3f}  2024={ann_2024:+.3f} °C  "
            f"elapsed={elapsed/60:.1f}m  eta={eta/60:.1f}m"
        )

    suf = args.out_suffix or ""
    out_path = OUT if not suf else OUT.with_stem(OUT.stem + suf)
    base = pd.DataFrame({"year": years, "month": months})
    df = pd.concat([base, pd.DataFrame(member_cols)], axis=1)
    df.to_csv(out_path, index=False, float_format="%.6f")
    print(f"\nSaved {out_path} ({len(df)} rows × {df.shape[1]} cols)")


if __name__ == "__main__":
    main()

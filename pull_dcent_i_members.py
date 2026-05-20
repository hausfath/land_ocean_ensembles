"""DCENT-I v1.1.0.0 (Chan et al. 2026, Geoscience Data Journal) —
Dynamically Consistent ENsemble of Temperature, Infilled.

Iteratively pulls each of 200 member NetCDFs from the Harvard Dataverse
deposit at doi:10.7910/DVN/ROG38Q. Each member file contains separate
`sst` (sea_surface_temperature_anomaly), `lsat` (land near-surface air
temperature anomaly), AND `ts` (merged surface temperature anomaly) on a
5°×5° monthly grid from 1850-01 through 2025-12.

For our build we only need `sst` and `lsat`. Each is finite over its
domain (sst over open ocean + sea-ice cells; lsat over land + sea-ice
cells), with sea-ice cells appearing in BOTH fields via the DCENT-I
kriging-with-air-temp-blending design. We compute area-weighted
(cos lat) global means over each field's finite cells per member per
month, then delete the raw NetCDF.

Outputs:
    Land Data/DCLSAT/DCLSAT_I_monthly_ensemble.csv          (200 members × 2112 months)
    Ocean Data/DCENT SST/DCENT_I_SST_monthly_ensemble.csv   (200 members × 2112 months)

Peak transient disk: ~42 MB. Total transfer: ~8.4 GB. Wall time ~17 min.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).parent
LAND_OUT = ROOT / "Land Data" / "DCLSAT" / "DCLSAT_I_monthly_ensemble.csv"
SST_OUT = ROOT / "Ocean Data" / "DCENT SST" / "DCENT_I_SST_monthly_ensemble.csv"
TMP_DIR = ROOT / "temp_files" / "dcent_i_members"
TMP_DIR.mkdir(parents=True, exist_ok=True)

MANIFEST = ROOT / "dcent_i_member_fileids.json"
DATAVERSE_URL = "https://dataverse.harvard.edu/api/access/datafile/{fid}"


def download(fid: int, dest: Path, attempts: int = 4) -> None:
    url = DATAVERSE_URL.format(fid=fid)
    for k in range(attempts):
        try:
            subprocess.run(
                ["curl", "-L", "--fail", "-sS", "-o", str(dest), url],
                check=True,
            )
            if dest.stat().st_size < 1_000_000:
                raise RuntimeError(f"unexpectedly small file ({dest.stat().st_size} B)")
            return
        except (subprocess.CalledProcessError, RuntimeError) as exc:
            wait = 2 ** k
            print(f"  retry {k+1}/{attempts} after {wait}s: {exc}")
            if dest.exists():
                dest.unlink()
            time.sleep(wait)
    raise RuntimeError(f"failed to download fileId {fid} after {attempts} attempts")


def area_weighted_global_mean(field: np.ndarray, cos_lat: np.ndarray) -> np.ndarray:
    """field: (n_months, lat, lon) with NaN over non-domain cells.
       cos_lat: (lat, lon). Returns (n_months,) global mean."""
    out = np.empty(field.shape[0], dtype=np.float64)
    for k in range(field.shape[0]):
        m = np.isfinite(field[k])
        contrib = np.where(m, field[k] * cos_lat, 0.0).sum()
        denom = (m * cos_lat).sum()
        out[k] = contrib / denom if denom > 0 else np.nan
    return out


def member_global_means(
    nc_path: Path, cos_lat: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Read one member's `sst` and `lsat` fields, compute their separate
    area-weighted global monthly time series. Returns (lsat_global, sst_global)
    each shape (n_months,)."""
    with xr.open_dataset(nc_path) as ds:
        sst = ds["sst"].values
        lsat = ds["lsat"].values
        if sst.ndim == 4:
            sst = sst.squeeze()
            lsat = lsat.squeeze()
    return area_weighted_global_mean(lsat, cos_lat), area_weighted_global_mean(sst, cos_lat)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--first", type=int, default=1)
    ap.add_argument("--last", type=int, default=200)
    ap.add_argument("--out-suffix", type=str, default=None,
                    help="if given, suffix appended to output CSV names (for partial runs)")
    args = ap.parse_args()

    with open(MANIFEST) as f:
        mapping = {int(k): int(v) for k, v in json.load(f).items()}
    assert sorted(mapping) == list(range(1, 201))

    # Time + lat axis from one member (downloaded as part of the run)
    first_fid = mapping[args.first]
    first_tmp = TMP_DIR / f"m{args.first:03d}.nc"
    if not first_tmp.exists():
        print(f"Downloading member {args.first:03d} to read axes ...")
        download(first_fid, first_tmp)
    with xr.open_dataset(first_tmp) as ds0:
        time_axis = pd.DatetimeIndex(ds0["time"].values)
        lat = ds0["lat"].values
        lon = ds0["lon"].values
    years = time_axis.year
    months = time_axis.month
    n_months = len(time_axis)
    cos_lat_2d = np.broadcast_to(
        np.cos(np.deg2rad(lat))[:, None], (len(lat), len(lon))
    ).astype(np.float64).copy()
    print(f"time: {time_axis[0]} to {time_axis[-1]}  (n={n_months} months)")

    lsat_cols: dict[str, np.ndarray] = {}
    sst_cols: dict[str, np.ndarray] = {}

    t0 = time.time()
    n_members = args.last - args.first + 1
    for k, member in enumerate(range(args.first, args.last + 1), start=1):
        tmp_nc = TMP_DIR / f"m{member:03d}.nc"
        if not tmp_nc.exists():
            download(mapping[member], tmp_nc)

        lsat_m, sst_m = member_global_means(tmp_nc, cos_lat_2d)
        lsat_cols[f"m{member:03d}"] = lsat_m
        sst_cols[f"m{member:03d}"] = sst_m

        tmp_nc.unlink()

        elapsed = time.time() - t0
        rate = k / elapsed if elapsed > 0 else 0.0
        eta = (n_members - k) / rate if rate > 0 else float("nan")
        ann2024_lsat = np.nanmean(lsat_m[years == 2024])
        ann2024_sst = np.nanmean(sst_m[years == 2024])
        ann2025_lsat = np.nanmean(lsat_m[years == 2025])
        ann2025_sst = np.nanmean(sst_m[years == 2025])
        print(
            f"[{k:3d}/{n_members}] m{member:03d}  "
            f"2024 L={ann2024_lsat:+.3f} S={ann2024_sst:+.3f}  "
            f"2025 L={ann2025_lsat:+.3f} S={ann2025_sst:+.3f}  "
            f"elapsed={elapsed/60:.1f}m  eta={eta/60:.1f}m"
        )

    suf = args.out_suffix or ""
    lsat_out = LAND_OUT if not suf else LAND_OUT.with_stem(LAND_OUT.stem + suf)
    sst_out = SST_OUT if not suf else SST_OUT.with_stem(SST_OUT.stem + suf)

    base = pd.DataFrame({"year": years, "month": months})
    lsat_df = pd.concat([base, pd.DataFrame(lsat_cols)], axis=1)
    sst_df = pd.concat([base, pd.DataFrame(sst_cols)], axis=1)
    lsat_df.to_csv(lsat_out, index=False, float_format="%.6f")
    sst_df.to_csv(sst_out, index=False, float_format="%.6f")
    print(
        f"\nSaved:\n  {lsat_out}  ({len(lsat_df)} rows × {lsat_df.shape[1]} cols)\n  "
        f"{sst_out}  ({len(sst_df)} rows × {sst_df.shape[1]} cols)"
    )


if __name__ == "__main__":
    main()

"""
Pull all 200 DCENT v3.0 member NetCDFs from Harvard Dataverse, compute the
ocean-area-weighted global-mean SST anomaly per member per month, and save
a single CSV holding the 200-member monthly ensemble. Raw NetCDFs are
deleted after each successful global-mean computation to keep disk use
bounded.

Output:
    Ocean Data/DCENT SST/DCENT_SST_monthly_ensemble.csv
        (time × 202: year, month, m001..m200)

Mirrors `pull_dcent_lsat_members.py` but operates on the `sst` field. If
the LSAT pull was previously run on this machine, the raw NetCDFs were
already deleted; we re-download each one.
"""
from __future__ import annotations

import json
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "Ocean Data" / "DCENT SST"
TMP_DIR = ROOT / "temp_files" / "dcent_members_sst"
TMP_DIR.mkdir(parents=True, exist_ok=True)
OUT_CSV = OUT_DIR / "DCENT_SST_monthly_ensemble.csv"

MANIFEST = ROOT / "dcent_member_fileids.json"
DATAVERSE_URL = "https://dataverse.harvard.edu/api/access/datafile/{fid}"


def area_weighted_ocean_mean(da: xr.DataArray) -> np.ndarray:
    """Area-weighted (cos(lat)) global mean of `da`, restricted to cells where
    `da` is finite (i.e. ocean). Returns 1-D array indexed by time."""
    weights = np.cos(np.deg2rad(da["lat"]))
    weights_b = xr.broadcast(weights, da)[0]
    mask = np.isfinite(da)
    w = weights_b.where(mask, 0.0)
    num = (da.where(mask, 0.0) * w).sum(dim=("lat", "lon"))
    den = w.sum(dim=("lat", "lon"))
    return (num / den).values


def download(fid: int, dest: Path, attempts: int = 4) -> None:
    url = DATAVERSE_URL.format(fid=fid)
    for k in range(attempts):
        try:
            subprocess.run(
                ["curl", "-L", "--fail", "-sS", "-o", str(dest), url],
                check=True,
            )
            return
        except subprocess.CalledProcessError as exc:
            wait = 2 ** k
            print(f"  download retry {k+1}/{attempts} after {wait}s: {exc}")
            time.sleep(wait)
    raise RuntimeError(f"failed to download fileId {fid} after {attempts} attempts")


def main() -> None:
    with open(MANIFEST) as f:
        mapping = {int(k): int(v) for k, v in json.load(f).items()}
    assert sorted(mapping) == list(range(1, 201))

    rows: pd.DataFrame | None = None
    member_cols: dict[str, np.ndarray] = {}

    t0 = time.time()
    for member in range(1, 201):
        fid = mapping[member]
        tmp_nc = TMP_DIR / f"dcent_m{member:03d}.nc"

        if not tmp_nc.exists():
            download(fid, tmp_nc)

        with xr.open_dataset(tmp_nc) as ds:
            series = area_weighted_ocean_mean(ds["sst"])
            time_index = ds["time"].values

        df = pd.DataFrame({
            "year": pd.DatetimeIndex(time_index).year,
            "month": pd.DatetimeIndex(time_index).month,
        })
        if rows is None:
            rows = df.copy()
        member_cols[f"m{member:03d}"] = series

        tmp_nc.unlink()

        elapsed = time.time() - t0
        rate = member / elapsed if elapsed > 0 else 0
        eta = (200 - member) / rate if rate > 0 else float("nan")
        print(
            f"[{member:3d}/200] mean={series.mean():+.3f}  "
            f"elapsed={elapsed/60:.1f}m  eta={eta/60:.1f}m"
        )

    assert rows is not None
    out = pd.concat([rows, pd.DataFrame(member_cols)], axis=1)
    out.to_csv(OUT_CSV, index=False, float_format="%.6f")
    print(f"\nSaved {OUT_CSV} ({len(out)} rows × {out.shape[1]} cols)")


if __name__ == "__main__":
    main()

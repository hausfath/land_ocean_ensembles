"""
Pull the ERSSTv6 native ensemble (Huang et al., NCEI pre-release).
Iteratively: download one member's gridded Fortran-sequential big-endian file,
compute the ocean-area-weighted global-mean monthly SSTA, append to an
in-memory ensemble, then delete the raw .dat file. Final write is a single
CSV with monthly global-mean anomalies for all members.

Source:
    https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/tmp/ersstv6.ensemble/
    sst2d.ano.1850.2024.1991-2020.ens.NNNN.dat  (NNNN = 0001..1000)

Format (per Boyin Huang's README and the GrADS .ctl):
    Fortran sequential unformatted, big_endian
    real ssta(180,89)   -- 2-deg lon (0..358E), 2-deg lat (-88..+88)
    do ny=1850,2024 ; do nm=1,12 ; read(11) ssta ; enddo ; enddo
    Undef sentinel: -999.9 °C

Output:
    Ocean Data/ERSSTv6/ERSSTv6_monthly_ensemble.csv
        rows: 2100 (Jan 1850 .. Dec 2024)
        cols: year, month, m0001..mNNNN
"""
from __future__ import annotations

import argparse
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import FortranFile

ROOT = Path(__file__).parent
OUT_DIR = ROOT / "Ocean Data" / "ERSSTv6"
TMP_DIR = ROOT / "temp_files" / "ersstv6_members"
TMP_DIR.mkdir(parents=True, exist_ok=True)

NLON, NLAT = 180, 89
NYEARS = 175                # 1850..2024 inclusive
NMON = NYEARS * 12
UNDEF = -999.9

URL_FMT = (
    "https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/tmp/ersstv6.ensemble/"
    "sst2d.ano.1850.2024.1991-2020.ens.{member:04d}.dat"
)

# Latitude grid: -88..+88 in 2° steps (89 cells). Cosine-area weights.
LAT = -88.0 + 2.0 * np.arange(NLAT)
LAT_WEIGHTS = np.cos(np.deg2rad(LAT))[:, None]  # (89, 1) broadcastable


def download(url: str, dest: Path, attempts: int = 4) -> None:
    for k in range(attempts):
        try:
            # --connect-timeout / --max-time guard against stale TCP routes
            # after a system sleep (otherwise curl hangs indefinitely on a
            # half-open connection). 300 s is well above typical transfer
            # times for the ~135 MB members (needs ≥0.45 MB/s).
            subprocess.run(
                ["curl", "-L", "--fail", "-sS",
                 "--connect-timeout", "30", "--max-time", "300",
                 "-o", str(dest), url],
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
    raise RuntimeError(f"failed to download {url} after {attempts} attempts")


def member_global_monthly(path: Path) -> np.ndarray:
    """Read one member's .dat file and return (NMON,) global ocean SSTA."""
    ff = FortranFile(path, "r", header_dtype=np.dtype(">u4"))
    out = np.empty(NMON, dtype=np.float64)
    for k in range(NMON):
        grid = ff.read_reals(dtype=">f4").reshape(NLAT, NLON)
        # mask the GrADS undef sentinel (-999.9) and any other obvious fills
        m = (grid > -50.0) & np.isfinite(grid)
        w = LAT_WEIGHTS * m  # (89,180)
        denom = w.sum()
        if denom > 0:
            out[k] = float((grid * w).sum()) / denom
        else:
            out[k] = np.nan
    ff.close()
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--first", type=int, default=1, help="first member (1-based)")
    ap.add_argument("--last", type=int, default=1000, help="last member (inclusive)")
    ap.add_argument(
        "--out", type=str, default=None,
        help="output CSV (default: ERSSTv6_monthly_ensemble.csv, or _firstN.csv for partial runs)",
    )
    args = ap.parse_args()

    first, last = args.first, args.last
    n_members = last - first + 1
    assert 1 <= first <= last <= 1000

    if args.out is None:
        if first == 1 and last == 1000:
            out_csv = OUT_DIR / "ERSSTv6_monthly_ensemble.csv"
        else:
            out_csv = OUT_DIR / f"ERSSTv6_monthly_ensemble_m{first:04d}_m{last:04d}.csv"
    else:
        out_csv = Path(args.out)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    # year/month columns
    years = np.repeat(np.arange(1850, 2025), 12)
    months = np.tile(np.arange(1, 13), NYEARS)

    member_cols: dict[str, np.ndarray] = {}
    t0 = time.time()
    for k, member in enumerate(range(first, last + 1), start=1):
        url = URL_FMT.format(member=member)
        tmp = TMP_DIR / f"ens.{member:04d}.dat"
        if not tmp.exists():
            download(url, tmp)
        series = member_global_monthly(tmp)
        member_cols[f"m{member:04d}"] = series
        tmp.unlink()

        elapsed = time.time() - t0
        rate = k / elapsed if elapsed > 0 else 0.0
        eta = (n_members - k) / rate if rate > 0 else float("nan")
        ann_2024 = np.nanmean(series[years == 2024])
        print(
            f"[{k:4d}/{n_members}] member {member:04d}  2024={ann_2024:+.3f} °C  "
            f"elapsed={elapsed/60:.1f}m  eta={eta/60:.1f}m"
        )

    out = pd.DataFrame({"year": years, "month": months})
    for name, series in member_cols.items():
        out[name] = series
    out.to_csv(out_csv, index=False, float_format="%.6f")
    print(f"\nSaved {out_csv}  ({len(out)} rows × {out.shape[1]} cols)")


if __name__ == "__main__":
    main()

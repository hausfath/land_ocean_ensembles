"""
Combined land + ocean ensemble figure, Berkeley-Earth-style.
Pattern: `reference_image.png` in this folder.

Inputs:
    land_ensemble.csv
    ocean_ensemble.csv

Output:
    land_ocean_ensemble.png

Baseline: percentiles are computed on the **1981-2010 baseline** (which
is how `land_ensemble.csv` / `ocean_ensemble.csv` are stored: each
member has mean ~0 over 1981-2010). For display we then add a
constant offset so the **median trace** has zero mean over 1850-1900,
i.e. anomalies are plotted relative to a 1850-1900 reference frame
*without* collapsing the uncertainty band in the well-sampled modern
era. This is the "modern baseline plus offset" approach used in
Thorne et al. (2026) and a number of Berkeley Earth GSAT figures.
Bands: 5th-95th percentile.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).parent

DISPLAY_BASELINE = (1850, 1900)   # plot anomalies relative to this period
SPREAD_BASELINE  = (1981, 2010)   # but compute percentiles on the modern baseline


def rebaseline_and_percentiles(csv_path: Path) -> pd.DataFrame:
    """Percentiles on 1981-2010 baseline (members already centred there in CSV),
    plus a constant offset so the median has zero mean over 1850-1900.
    The offset is reported separately for transparency."""
    df = pd.read_csv(csv_path)
    years = df["year"].values
    mems = df.iloc[:, 1:].values    # (n_years, n_members), already on 1981-2010 baseline

    p05 = np.nanpercentile(mems,  5, axis=1)
    p50 = np.nanpercentile(mems, 50, axis=1)
    p95 = np.nanpercentile(mems, 95, axis=1)

    # constant offset: shift so the median's 1850-1900 mean is zero
    pre_mask = (years >= DISPLAY_BASELINE[0]) & (years <= DISPLAY_BASELINE[1])
    offset = -float(np.nanmean(p50[pre_mask]))

    return pd.DataFrame({
        "year": years,
        "p05":  p05 + offset,
        "p50":  p50 + offset,
        "p95":  p95 + offset,
    }), offset


def main() -> None:
    land, land_offset = rebaseline_and_percentiles(ROOT / "land_ensemble.csv")
    ocean, ocean_offset = rebaseline_and_percentiles(ROOT / "ocean_ensemble.csv")
    print(f"Land  1981-2010 → 1850-1900 offset: {land_offset:+.3f} °C")
    print(f"Ocean 1981-2010 → 1850-1900 offset: {ocean_offset:+.3f} °C")

    LAND_COLOR = "#c0392b"
    OCEAN_COLOR = "#2455a3"

    fig, ax = plt.subplots(figsize=(9.5, 5.2))

    # Shaded 5-95% bands
    ax.fill_between(ocean["year"], ocean["p05"], ocean["p95"],
                    color=OCEAN_COLOR, alpha=0.22, linewidth=0)
    ax.fill_between(land["year"], land["p05"], land["p95"],
                    color=LAND_COLOR, alpha=0.22, linewidth=0)

    # Dotted-marker + thin connecting line for the medians
    ax.plot(ocean["year"], ocean["p50"], color=OCEAN_COLOR, lw=0.7,
            marker="o", markersize=2.2, markerfacecolor=OCEAN_COLOR,
            markeredgecolor=OCEAN_COLOR, zorder=3)
    ax.plot(land["year"], land["p50"], color=LAND_COLOR, lw=0.7,
            marker="o", markersize=2.2, markerfacecolor=LAND_COLOR,
            markeredgecolor=LAND_COLOR, zorder=4)

    ax.axhline(0, color="0.4", lw=0.5, ls="-", alpha=0.6)

    ax.set_xlim(1850, 2025)
    # set y-limits with a small margin around data
    ymin = min(land["p05"].min(), ocean["p05"].min()) - 0.15
    ymax = max(land["p95"].max(), ocean["p95"].max()) + 0.20
    ax.set_ylim(ymin, ymax)

    ax.set_title("Land and Ocean Temperatures 1850-2025",
                 fontsize=15, loc="left", pad=8)
    ax.set_ylabel("Temperature Anomaly (°C)", fontsize=11)
    ax.set_xlabel("")
    ax.tick_params(axis="both", labelsize=10)
    ax.grid(False)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color("0.4")
    ax.spines["bottom"].set_color("0.4")

    # Inline labels near the right end of each trace
    last_year = int(land["year"].iloc[-1])
    land_y = float(land["p50"].iloc[-1])
    ocean_y = float(ocean["p50"].iloc[-1])
    # place labels offset from the latest point, on the curve itself a few years back
    label_year = 2000
    iL = int(land[land["year"] == label_year].index[0])
    iO = int(ocean[ocean["year"] == label_year].index[0])
    ax.annotate("Land Average",
                xy=(label_year, land["p50"].iloc[iL]),
                xytext=(label_year - 8, land["p50"].iloc[iL] + 0.35),
                color=LAND_COLOR, fontsize=12,
                ha="center", fontweight="normal")
    ax.annotate("Ocean Average",
                xy=(label_year, ocean["p50"].iloc[iO]),
                xytext=(label_year - 8, ocean["p50"].iloc[iO] - 0.35),
                color=OCEAN_COLOR, fontsize=12,
                ha="center", fontweight="normal")

    # Footer annotation
    footer = (
        "Temperature anomalies relative to 1850-1900 average; 5th-95th percentile range computed on 1981-2010 baseline\n"
        "Land ensemble from 6 LSAT datasets; ocean ensemble from 4 SST datasets (family-tree weighted, 10,000 members each)"
    )
    ax.text(0.99, 0.04, footer, transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8, color="0.30",
            linespacing=1.3)

    fig.tight_layout()
    out = ROOT / "land_ocean_ensemble.png"
    fig.savefig(out, dpi=200)
    print(f"Saved {out}")


if __name__ == "__main__":
    main()

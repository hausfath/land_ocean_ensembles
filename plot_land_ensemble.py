"""
Plot the 10,000-member land LSAT ensemble: central estimate + 90% band,
with the five per-dataset best estimates overlaid for cross-check.

Inputs:
    land_ensemble_summary.csv
    land_ensemble_perdataset.csv

Output:
    land_ensemble.png
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).parent


def main() -> None:
    summary = pd.read_csv(ROOT / "land_ensemble_summary.csv")
    perds = pd.read_csv(ROOT / "land_ensemble_perdataset.csv")

    fig, ax = plt.subplots(figsize=(10, 5.5))

    ax.fill_between(
        summary["year"], summary["p025"], summary["p975"],
        color="0.7", alpha=0.6, label="Ensemble 2.5–97.5 percentile",
    )
    ax.plot(summary["year"], summary["mean"], "k-", lw=1.6, label="Ensemble mean")

    colors = {
        "berkeley_earth": "#d73027",
        "crutem5":        "#4575b4",
        "glosatlat":      "#74add1",
        "dclsat":         "#5aae61",
        "noaa_land":      "#fdae61",
    }
    labels = {
        "berkeley_earth": "Berkeley Earth Highres",
        "crutem5":        "CRUTEM5",
        "glosatlat":      "GloSATLAT (1850-2021)",
        "dclsat":         "DCLSAT (DCENT v3.0 land)",
        "noaa_land":      "NOAA Land v6.1",
    }
    for name, color in colors.items():
        ax.plot(
            perds["year"], perds[f"{name}_mean"],
            color=color, lw=0.9, alpha=0.85, label=labels[name],
        )

    ax.axhline(0, color="k", lw=0.4, ls="--", alpha=0.4)
    ax.set_xlim(1850, 2025)
    ax.set_xlabel("Year")
    ax.set_ylabel("Land surface air temperature anomaly (°C, 1981–2010 baseline)")
    ax.set_title("Global LSAT ensemble (10,000 members, family-tree weighted)")
    ax.legend(loc="upper left", fontsize=8, ncol=2, frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(ROOT / "land_ensemble.png", dpi=180)
    print(f"Saved {ROOT / 'land_ensemble.png'}")


if __name__ == "__main__":
    main()

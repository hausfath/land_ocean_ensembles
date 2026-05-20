"""Plot the 10,000-member ocean SST ensemble + per-dataset best estimates."""
from __future__ import annotations
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).parent


def main() -> None:
    summary = pd.read_csv(ROOT / "ocean_ensemble_summary.csv")
    perds = pd.read_csv(ROOT / "ocean_ensemble_perdataset.csv")

    fig, ax = plt.subplots(figsize=(10, 5.5))
    ax.fill_between(summary["year"], summary["p025"], summary["p975"],
                    color="0.7", alpha=0.6, label="Ensemble 2.5–97.5 percentile")
    ax.plot(summary["year"], summary["mean"], "k-", lw=1.6, label="Ensemble mean")

    colors = {
        "hadsst4":   "#1f78b4",
        "ersstv6":   "#fdae61",
        "cobe_sst2": "#74add1",
        "dcent_sst": "#5aae61",
    }
    labels = {
        "hadsst4":   "HadSST 4.2.0.0",
        "ersstv6":   "ERSSTv6",
        "cobe_sst2": "COBE-SST2",
        "dcent_sst": "DCENT-I SST (v1.1.0.0)",
    }
    for name, color in colors.items():
        ax.plot(perds["year"], perds[f"{name}_mean"], color=color, lw=0.9, alpha=0.85, label=labels[name])

    ax.axhline(0, color="k", lw=0.4, ls="--", alpha=0.4)
    ax.set_xlim(1850, 2025)
    ax.set_xlabel("Year")
    ax.set_ylabel("Sea-surface temperature anomaly (°C, 1981–2010 baseline)")
    ax.set_title("Global SST ensemble (10,000 members, family-tree weighted)")
    ax.legend(loc="upper left", fontsize=8, ncol=2, frameon=False)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(ROOT / "ocean_ensemble.png", dpi=180)
    print(f"Saved {ROOT / 'ocean_ensemble.png'}")


if __name__ == "__main__":
    main()

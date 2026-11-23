#!/usr/bin/env python3
"""
Land/ocean warming ratio from the AR7 land and ocean GMST ensembles.

Computes, per ensemble member, the warming between 1850-1900 and a recent
period, then forms the land-to-ocean ratio by Monte Carlo pairing of the two
(independent) 10,000-member ensembles. Reports median + 2.5-97.5% range.
Note: independent pairing is conservative (wide) for the ratio/difference,
since land and ocean errors share structural components (positive correlation
would narrow both intervals). Quoted central ratio is the ratio of medians.

Validation: the ocean warming reproduces the assessed SST values
(0.90 [0.70 to 1.09] degC, 5-95%, for 1850-1900 -> 2016-2025).

Inputs : land_ensemble.csv, ocean_ensemble.csv (year + 10,000 member cols,
         annual anomalies 1850-2025, each member zero-mean over its baseline).
Output : prints summary; writes land_ocean_ratio_summary.csv
"""
import pandas as pd, numpy as np

SEED = 42
N_MC = 2_000_000
BASE = (1850, 1900)
PERIODS = {"2016-2025": (2016, 2025), "2011-2020": (2011, 2020)}


def warming(df, recent, base=BASE):
    """Per-member warming = mean(recent) - mean(base), shape (n_members,)."""
    yr = df["year"].values
    M = df.drop(columns="year").values
    return (M[(yr >= recent[0]) & (yr <= recent[1])].mean(0)
            - M[(yr >= base[0]) & (yr <= base[1])].mean(0))


def summ(x):
    return np.percentile(x, 50), np.percentile(x, 2.5), np.percentile(x, 97.5)


def main():
    rng = np.random.default_rng(SEED)
    land = pd.read_csv("land_ensemble.csv")
    ocean = pd.read_csv("ocean_ensemble.csv")
    rows = []
    for label, period in PERIODS.items():
        L = warming(land, period)
        O = warming(ocean, period)
        # Independent random pairing with replacement -> ratio distribution.
        li = rng.integers(0, len(L), N_MC)
        oi = rng.integers(0, len(O), N_MC)
        ratio = L[li] / O[oi]
        diff = L[li] - O[oi]
        lm, lo, lh = summ(L); om, oo, oh = summ(O)
        rm, ro, rh = summ(ratio); dm, do, dh = summ(diff)
        print(f"=== 1850-1900 -> {label} ===")
        print(f"  land   {lm:.2f} [{lo:.2f} to {lh:.2f}] degC")
        print(f"  ocean  {om:.2f} [{oo:.2f} to {oh:.2f}] degC")
        print(f"  ratio  {rm:.2f} [{ro:.2f} to {rh:.2f}]  (ratio of medians {lm/om:.2f})")
        print(f"  diff   {dm:.2f} [{do:.2f} to {dh:.2f}] degC")
        rows.append(dict(period=label, land_med=lm, land_lo=lo, land_hi=lh,
                         ocean_med=om, ocean_lo=oo, ocean_hi=oh,
                         ratio_med=rm, ratio_lo=ro, ratio_hi=rh,
                         ratio_of_medians=lm/om,
                         diff_med=dm, diff_lo=do, diff_hi=dh))
    pd.DataFrame(rows).to_csv("land_ocean_ratio_summary.csv", index=False)
    print("\nwrote land_ocean_ratio_summary.csv")


if __name__ == "__main__":
    main()

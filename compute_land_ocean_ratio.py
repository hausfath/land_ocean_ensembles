#!/usr/bin/env python3
"""
Land/ocean warming ratio from the AR7 land and ocean GMST ensembles.

Computes, per ensemble member, the warming between 1850-1900 and a recent
period, then forms the land-to-ocean ratio by Monte Carlo pairing of the two
(independent) 10,000-member ensembles. Reports median + 5-95% range.
Also reports land/ocean ratios of per-member linear trends (OLS, degC/decade).
Independent pairing ignores land-ocean error dependence, whose sign is
uncertain (see land_ocean_ratio_methodology.md). Quoted central ratio is the
ratio of medians.

Validation: the ocean warming reproduces the assessed SST values
(0.90 [0.70 to 1.09] degC, 5-95%, for 1850-1900 -> 2016-2025).

Inputs : land_ensemble.csv, ocean_ensemble.csv (year + 10,000 member cols,
         annual anomalies 1850-2025, each member zero-mean over its baseline).
Output : prints summary; writes land_ocean_ratio_summary.csv (headline,
         independent pairing) and land_ocean_ratio_sensitivity.csv (rank-
         correlated pairing at RHOS; producer-matched land/ocean pairs from
         the per-dataset ensemble means).
"""
import pandas as pd, numpy as np
from scipy.stats import norm

SEED = 42
N_MC = 2_000_000
BASE = (1850, 1900)
PERIODS = {"2016-2025": (2016, 2025), "2011-2020": (2011, 2020)}
RHOS = (0.3, 0.5, 0.7)
SENS_PERIODS = {"warming 2016-2025": ("warming", (2016, 2025)),
                "trend 1950-2025": ("trend", (1950, 2025)),
                "trend 1981-2025": ("trend", (1981, 2025))}
# Producer's own land and ocean series (columns of the *_perdataset.csv files).
PRODUCERS = {"HadCRUT5": ("crutem5", "hadsst4"),
             "Berkeley Earth": ("berkeley_earth", "hadsst4"),
             "DCENT": ("dclsat", "dcent_sst"),
             "NOAAGlobalTemp": ("noaa_land", "ersstv6"),
             "CMST3": ("c_lsat", "ersstv6")}
TREND_PERIODS = {f"{a}-{b}": (a, b) for a, b in
                 [(1850, 2025), (1901, 2025), (1950, 2025), (1960, 2025),
                  (1979, 2025), (1981, 2025), (1991, 2025)]}


def warming(df, recent, base=BASE):
    """Per-member warming = mean(recent) - mean(base), shape (n_members,)."""
    yr = df["year"].values
    M = df.drop(columns="year").values
    return (M[(yr >= recent[0]) & (yr <= recent[1])].mean(0)
            - M[(yr >= base[0]) & (yr <= base[1])].mean(0))


def trend(df, period):
    """Per-member OLS trend over `period` (inclusive), degC/decade."""
    yr = df["year"].values
    M = df.drop(columns="year").values
    s = (yr >= period[0]) & (yr <= period[1])
    x = yr[s] - yr[s].mean()
    return (x[:, None] * (M[s] - M[s].mean(0))).sum(0) / (x ** 2).sum() * 10


def summ(x):
    return np.percentile(x, 50), np.percentile(x, 5), np.percentile(x, 95)


def copula_draws(L, O, rho, n, rng):
    """Pair land and ocean members with Gaussian-copula rank correlation rho."""
    z1 = rng.standard_normal(n)
    z2 = rho * z1 + np.sqrt(1 - rho ** 2) * rng.standard_normal(n)
    li = np.minimum((norm.cdf(z1) * len(L)).astype(int), len(L) - 1)
    oi = np.minimum((norm.cdf(z2) * len(O)).astype(int), len(O) - 1)
    return np.sort(L)[li], np.sort(O)[oi]


def series_stat(s, metric, period):
    s = s.dropna()
    if metric == "warming":
        return s.loc[period[0]:period[1]].mean() - s.loc[BASE[0]:BASE[1]].mean()
    s = s.loc[period[0]:period[1]]
    x = s.index.values - s.index.values.mean()
    return (x * (s.values - s.values.mean())).sum() / (x ** 2).sum() * 10


def sensitivity(land, ocean, rng):
    rows = []
    for name, (metric, period) in SENS_PERIODS.items():
        fn = warming if metric == "warming" else trend
        L, O = fn(land, period), fn(ocean, period)
        for rho in RHOS:
            l, o = copula_draws(L, O, rho, N_MC, rng)
            (rm, ro, rh), (dm, do, dh) = summ(l / o), summ(l - o)
            rows.append(dict(test=f"rank correlation {rho}", quantity=name,
                             ratio_med=rm, ratio_lo=ro, ratio_hi=rh,
                             diff_med=dm, diff_lo=do, diff_hi=dh))
    lp = pd.read_csv("land_ensemble_perdataset.csv").set_index("year")
    op = pd.read_csv("ocean_ensemble_perdataset.csv").set_index("year")
    for prod, (lc, oc) in PRODUCERS.items():
        for name, (metric, period) in SENS_PERIODS.items():
            lv = series_stat(lp[f"{lc}_mean"], metric, period)
            ov = series_stat(op[f"{oc}_mean"], metric, period)
            rows.append(dict(test=f"producer pair: {prod} ({lc}/{oc})", quantity=name,
                             land=lv, ocean=ov, ratio_med=lv / ov, diff_med=lv - ov))
    out = pd.DataFrame(rows)
    out.to_csv("land_ocean_ratio_sensitivity.csv", index=False, float_format="%.4f")
    print("\n=== sensitivity (5-95%) ===")
    for r in out.itertuples():
        rng_txt = "" if pd.isna(r.ratio_lo) else f" [{r.ratio_lo:.2f}-{r.ratio_hi:.2f}]"
        print(f"  {r.test:45s} {r.quantity:18s} ratio {r.ratio_med:.2f}{rng_txt}  diff {r.diff_med:.3f}")
    print("wrote land_ocean_ratio_sensitivity.csv")


def main():
    rng = np.random.default_rng(SEED)
    land = pd.read_csv("land_ensemble.csv")
    ocean = pd.read_csv("ocean_ensemble.csv")
    rows = []
    jobs = [("warming", lbl, warming(land, p), warming(ocean, p)) for lbl, p in PERIODS.items()]
    jobs += [("trend", lbl, trend(land, p), trend(ocean, p)) for lbl, p in TREND_PERIODS.items()]
    for metric, label, L, O in jobs:
        # Independent random pairing with replacement -> ratio distribution.
        li = rng.integers(0, len(L), N_MC)
        oi = rng.integers(0, len(O), N_MC)
        ratio = L[li] / O[oi]
        diff = L[li] - O[oi]
        lm, lo, lh = summ(L); om, oo, oh = summ(O)
        rm, ro, rh = summ(ratio); dm, do, dh = summ(diff)
        unit = "degC" if metric == "warming" else "degC/decade"
        head = f"1850-1900 -> {label}" if metric == "warming" else f"trend {label}"
        print(f"=== {head} ===")
        f = (lambda v: f"{v:.2f}") if metric == "warming" else (lambda v: f"{v:.3f}")
        print(f"  land   {f(lm)} [{f(lo)} to {f(lh)}] {unit}")
        print(f"  ocean  {f(om)} [{f(oo)} to {f(oh)}] {unit}")
        print(f"  ratio  {rm:.2f} [{ro:.2f} to {rh:.2f}]  (ratio of medians {lm/om:.2f})")
        print(f"  diff   {f(dm)} [{f(do)} to {f(dh)}] {unit}")
        rows.append(dict(metric=metric, period=label, land_med=lm, land_lo=lo, land_hi=lh,
                         ocean_med=om, ocean_lo=oo, ocean_hi=oh,
                         ratio_med=rm, ratio_lo=ro, ratio_hi=rh,
                         ratio_of_medians=lm/om,
                         diff_med=dm, diff_lo=do, diff_hi=dh))
    pd.DataFrame(rows).to_csv("land_ocean_ratio_summary.csv", index=False)
    print("\nwrote land_ocean_ratio_summary.csv")
    sensitivity(land, ocean, rng)


if __name__ == "__main__":
    main()

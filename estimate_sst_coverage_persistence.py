"""Coverage error of HadSST4's global-mean SST, estimated as the Met Office does
(HadSST.4.2.0.0 Product User Guide: subsample the complete HadISST.1.1 field
with HadSST4's coverage over a fixed historical period; coverage uncertainty is
the RMS difference), extended to lags to give the year-to-year persistence.

For mask year t (HadSST4 monthly coverage, from HadSST.4.2.0.0_median.nc) and
HadISST year j: e(t, j) = annual mean over months of
[mean over observed boxes - mean over all ocean boxes] of HadISST anomalies.

Step A reproduces the published monthly and annual coverage σ (RMS over j) and
their ratio, for a few global-mean definitions and HadISST windows.
Step B: uncentred lag correlation
    r(L) = Σ_j e(t,j) e(t+L,j+L) / sqrt(Σ_j e(t,j)^2 · Σ_j e(t+L,j+L)^2),
averaged over mask years, fitted with corr(L) = (1-p) φ^L + p.

Outputs: sst_coverage_persistence.csv, sst_coverage_persistence_fit.json
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.optimize import least_squares

ROOT = Path(__file__).parent
DATA = ROOT / "Ocean Data"
FIT_LAGS = (1, 10)
MAX_LAG = 20
WINDOWS = {"1961-1990": (1961, 1990), "1979-2018": (1979, 2018), "1950-2020": (1950, 2020)}


def load_field():
    ds = xr.open_dataset(DATA / "HadISST" / "HadISST1_sst_5deg_monthly.nc")
    t = pd.DatetimeIndex(ds.time.values)
    yrs = np.unique(t.year)
    yrs = np.array([y for y in yrs if (t.year == y).sum() == 12 and y <= 2025])
    keep = np.isin(t.year, yrs)
    X = ds.sst.values[keep].reshape(len(yrs), 12, -1).astype(np.float64)
    ocean = np.isfinite(X).all(axis=(0, 1))
    clim = np.nanmean(X[(yrs >= 1961) & (yrs <= 1990)], axis=0, keepdims=True)
    X = np.where(ocean[None, None, :], X - clim, 0.0)
    return yrs, X, ocean, ds.latitude.values


def load_masks(ocean):
    ds = xr.open_dataset(DATA / "HadSST4" / "HadSST.4.2.0.0_median.nc")
    t = pd.DatetimeIndex(ds.time.values)
    ok = np.isfinite(ds.tos.values)
    myrs = np.arange(1850, 2026)
    M = np.zeros((len(myrs), 12, 36 * 72), bool)
    for i, y in enumerate(myrs):
        sel = t.year == y
        M[i, t.month[sel] - 1] = ok[sel].reshape(sel.sum(), -1)
    return myrs, M & ocean[None, None, :]


def errors(X, M, ocean, lat, variant):
    """e[t, j, m]: monthly coverage error for mask (t, m) on field year j."""
    coslat = np.repeat(np.cos(np.deg2rad(lat)), 72)
    nh = np.repeat(lat > 0, 72)
    parts = [(np.ones_like(nh), 1.0)] if variant == "global" else [(nh, 0.5), (~nh, 0.5)]
    e = 0.0
    for region, wgt in parts:
        wt = coslat * ocean * region
        truth = np.einsum("jmb,b->jm", X, wt) / wt.sum()
        W = M * (coslat * region)[None, None, :]
        den = W.sum(axis=2)
        masked = np.einsum("jmb,tmb->tjm", X, W) / np.where(den > 0, den, np.nan)[:, None, :]
        e = e + wgt * (masked - truth[None])
    return e


def main():
    yrs, X, ocean, lat = load_field()
    myrs, M = load_masks(ocean)
    pub_m = pd.read_csv(DATA / "HadSST4" / "HadSST.4.2.0.0_monthly_GLOBE.csv")
    pub_m = pub_m[pub_m.year <= 2025].groupby("year").coverage_uncertainty.apply(lambda s: np.sqrt(np.mean(s ** 2)))
    pub_a = pd.read_csv(DATA / "HadSST4" / "HadSST.4.2.0.0_annual_GLOBE.csv").set_index("year").coverage_uncertainty
    chk = [1860, 1880, 1900, 1920, 1950, 1975, 1995, 2010, 2025]
    print("Step A: coverage RMS vs published (monthly RMS-of-year / annual / ratio)")
    best, best_score, E_best = None, np.inf, None
    for variant in ("global", "hemispheric"):
        e = errors(X, M, ocean, lat, variant)
        for wname, (a, b) in WINDOWS.items():
            sel = (yrs >= a) & (yrs <= b)
            em = e[:, sel, :]
            mon = np.sqrt(np.nanmean(em ** 2, axis=(1, 2)))
            ann = np.sqrt(np.nanmean(np.nanmean(em, axis=2) ** 2, axis=1))
            r_m = mon / pub_m.reindex(myrs).values
            r_a = ann / pub_a.reindex(myrs).values
            ok = myrs >= 1850
            score = np.nanmedian(np.abs(np.log(r_a[ok]))) + np.nanmedian(np.abs(np.log(r_m[ok])))
            print(f"\n  {variant}, HadISST {wname}: ratio to published monthly median {np.nanmedian(r_m):.2f}, "
                  f"annual {np.nanmedian(r_a):.2f}; annual/monthly ratio median {np.nanmedian(ann / mon):.2f} (published 0.81)")
            for y in chk:
                i = y - 1850
                print(f"    {y}: monthly {mon[i]:.3f} vs {pub_m.get(y, np.nan):.3f}   annual {ann[i]:.3f} vs {pub_a.get(y, np.nan):.3f}")
            if score < best_score:
                best, best_score, E_best = (variant, wname), score, np.nanmean(e, axis=2)
    print(f"\nbest reproduction: {best}")
    # lag correlations over the same HadISST window that reproduces the published sigma
    lw = WINDOWS[best[1]]
    jsel = (yrs >= lw[0]) & (yrs <= lw[1])
    yrs, E_best = yrs[jsel], E_best[:, jsel]

    print("\nStep B: uncentred lag correlation of annual coverage error over HadISST years")
    rows = []
    epochs = {"1850-1899": (1850, 1899), "1900-1949": (1900, 1949), "1950-2005": (1950, 2005), "all": (1850, 2005)}
    for detrend in (False, True):
        E = E_best.copy()
        if detrend:
            jj = yrs - yrs.mean()
            for i in range(E.shape[0]):
                if np.isfinite(E[i]).all():
                    E[i] = E[i] - np.polyval(np.polyfit(jj, E[i], 1), jj)
        for ep, (a, b) in epochs.items():
            ts = [i for i, y in enumerate(myrs) if a <= y <= b]
            for L in range(0, MAX_LAG + 1):
                cs = []
                for i in ts:
                    if i + L >= len(myrs):
                        continue
                    x, y = E[i, : len(yrs) - L], E[i + L, L:]
                    if np.isfinite(x).all() and np.isfinite(y).all():
                        cs.append(np.sum(x * y) / np.sqrt(np.sum(x * x) * np.sum(y * y)))
                rows.append({"detrended": detrend, "epoch": ep, "lag": L, "corr": np.mean(cs), "n_masks": len(cs)})
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "sst_coverage_persistence.csv", index=False, float_format="%.4f")
    print(df.pivot_table(index="lag", columns=["detrended", "epoch"], values="corr").loc[[1, 2, 3, 5, 10, 15, 20]].round(3).to_string())
    fits = {}
    for detrend in (False, True):
        c = df[(df.detrended == detrend) & (df.epoch == "all") & (df.lag >= FIT_LAGS[0]) & (df.lag <= FIT_LAGS[1])]
        L, r = c["lag"].values, c["corr"].values
        sol = least_squares(lambda q: (1 - q[0]) * q[1] ** L + q[0] - r, x0=[0.3, 0.5], bounds=([0, 0], [1, 0.999]))
        fits["detrended" if detrend else "raw"] = {"p": float(sol.x[0]), "phi": float(sol.x[1])}
    fits.update({"definition": best[0], "hadisst_window_sigma_check": best[1],
                 "hadisst_years_lag": [int(yrs[0]), int(yrs[-1])], "note": "HadISST before ~1950 is a smooth reconstruction; lags use the sigma-check window", "moments": "uncentred (RMS), as HadSST4"})
    print("\nfit:", json.dumps(fits, indent=1))
    (ROOT / "sst_coverage_persistence_fit.json").write_text(json.dumps(fits, indent=1))


if __name__ == "__main__":
    main()

"""Coverage error of CRUTEM5's global land mean, estimated the Met Office way
(Morice et al. 2012 §5.3; Osborn et al. 2021): apply CRUTEM5's monthly coverage
mask for year t to every year j of a spatially complete reanalysis (ERA5,
reduced to 5° by prepare_era5_5deg.py) and take masked minus complete means,
hemispheres combined as (2 NH + SH)/3. The SD over j is the coverage sigma for
year t. Extended to lags: corr over j of e(t, j) and e(t+L, j+L) gives the
year-to-year persistence of the coverage error.

Step A reproduces the published annual coverage sigma (land-definition check).
Step B estimates the lag correlation and fits corr(L) = (1-p) phi^L + p.

Output: coverage_persistence.csv (lag correlations by epoch, raw and detrended)
        and coverage_persistence_fit.json (p, phi for the synthesis).
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from scipy.optimize import least_squares

ROOT = Path(__file__).parent
DATA = ROOT / "Land Data"
CRU_VERSION = "5.2.0.0"
REPRO_YEARS = (1979, 2018)     # ERA5 years used by HadCRUT5/CRUTEM5 for coverage sigma
LAG_YEARS = (1940, 2025)       # ERA5 years for the lag correlations
MAX_LAG = 20
FIT_LAGS = (1, 10)            # lags 15-16 carry an ERA5 sampling bump (present with fixed masks)


def load_era5():
    ds = xr.open_dataset(DATA / "ERA5" / "era5_t2m_5deg_monthly.nc")
    lf = ds.landfrac.values
    box = np.where(np.isfinite(ds.t2m_land.values), ds.t2m_land.values, ds.t2m_all.values)
    t = pd.DatetimeIndex(ds.time.values)
    keep = (t.year >= LAG_YEARS[0]) & (t.year <= LAG_YEARS[1])
    box, t = box[keep], t[keep]
    years = np.unique(t.year)
    X = box.reshape(len(years), 12, 36 * 72)
    X = X - X.mean(axis=0, keepdims=True)              # anomalies (baseline irrelevant)
    return years, X, lf.reshape(-1), ds.latitude.values


def load_masks():
    ds = xr.open_dataset(DATA / "CRUTEM5" / f"CRUTEM.{CRU_VERSION}.anomalies.nc")
    t = pd.DatetimeIndex(ds.time.values)
    ok = np.isfinite(ds.tas.values)
    yrs = np.arange(1850, 2026)
    M = np.zeros((len(yrs), 12, 36 * 72), bool)
    for i, y in enumerate(yrs):
        sel = t.year == y
        M[i, t.month[sel] - 1] = ok[sel].reshape(sel.sum(), -1)
    return yrs, M


def residuals(X, M, lf, lat, variant):
    """e[t, j] annual coverage error for mask-year t applied to ERA5 year j."""
    coslat = np.repeat(np.cos(np.deg2rad(lat)), 72)
    nh = np.repeat(lat > 0, 72)
    land = lf > 0
    if variant == "cos/cos":            # masked cos-weighted, truth cos over land boxes
        wm_base, wt = coslat, coslat * land
    elif variant == "cos/cos*lf":       # masked cos-weighted, truth land-area weighted
        wm_base, wt = coslat, coslat * lf
    elif variant == "cos*lf/cos*lf":
        wm_base, wt = coslat * np.maximum(lf, 0.05), coslat * lf
    out = []
    for hemi in (nh, ~nh):
        wt_h = wt * hemi
        truth = np.einsum("jmb,b->jm", X, wt_h) / wt_h.sum()             # (j, m)
        W = M * (wm_base * hemi)[None, None, :]                           # (t, m, b)
        den = W.sum(axis=2)                                               # (t, m)
        masked = np.einsum("jmb,tmb->tjm", X, W) / np.where(den > 0, den, np.nan)[:, None, :]
        out.append(masked - truth[None])                                  # (t, j, m)
    e_m = (2 * out[0] + out[1]) / 3
    return e_m.mean(axis=2), out[0].mean(axis=2), out[1].mean(axis=2), e_m


def main():
    years, X, lf, lat = load_era5()
    myrs, M = load_masks()
    pub = pd.read_csv(DATA / "CRUTEM5" / f"CRUTEM.{CRU_VERSION}.component_series.global.annual.csv").set_index("Time").iloc[:, 3]
    pubn = pd.read_csv(DATA / "CRUTEM5" / f"CRUTEM.{CRU_VERSION}.component_series.northern_hemisphere.annual.csv").set_index("Time").iloc[:, 3]
    pubs = pd.read_csv(DATA / "CRUTEM5" / f"CRUTEM.{CRU_VERSION}.component_series.southern_hemisphere.annual.csv").set_index("Time").iloc[:, 3]
    rsel = (years >= REPRO_YEARS[0]) & (years <= REPRO_YEARS[1])

    print("Step A: annual coverage sigma, ERA5 %d-%d, vs published (global / NH / SH)" % REPRO_YEARS)
    chk = [1860, 1880, 1900, 1930, 1950, 1970, 1990, 2010, 2020]
    results = {}
    for variant in ("cos/cos", "cos/cos*lf", "cos*lf/cos*lf"):
        e, en, es, _ = residuals(X, M, lf, lat, variant)
        sd = np.nanstd(e[:, rsel], axis=1, ddof=1)
        sdn = np.nanstd(en[:, rsel], axis=1, ddof=1)
        sds = np.nanstd(es[:, rsel], axis=1, ddof=1)
        ii = [list(myrs).index(y) for y in chk]
        ratio = sd[[list(myrs).index(y) for y in range(1857, 2026)]] / pub.loc[1857:2025].values
        print(f"\n  {variant}: ratio ERA5/published global median {np.nanmedian(ratio):.2f} "
              f"(5-95%: {np.nanpercentile(ratio, 5):.2f}-{np.nanpercentile(ratio, 95):.2f})")
        for y, i in zip(chk, ii):
            print(f"    {y}: {sd[i]:.3f} / {sdn[i]:.3f} / {sds[i]:.3f}   published {pub[y]:.3f} / {pubn[y]:.3f} / {pubs[y]:.3f}")
        results[variant] = (e, np.nanmedian(np.abs(np.log(ratio))))
    best = min(results, key=lambda k: results[k][1])
    print(f"\nclosest land definition: {best}")
    e = results[best][0]

    print("\nStep B: lag correlation of annual coverage error, ERA5 %d-%d" % LAG_YEARS)
    rows = []
    epochs = {"1857-1899": (1857, 1899), "1900-1949": (1900, 1949), "1950-2010": (1950, 2010), "all": (1857, 2010)}
    for detrend in (False, True):
        E = e.copy()
        if detrend:
            jj = years - years.mean()
            for i in range(E.shape[0]):
                if np.isfinite(E[i]).all():
                    b = np.polyfit(jj, E[i], 1)
                    E[i] = E[i] - np.polyval(b, jj)
        for ep, (a, b) in epochs.items():
            ts = [i for i, y in enumerate(myrs) if a <= y <= b]
            for L in range(0, MAX_LAG + 1):
                cs = []
                for i in ts:
                    if i + L >= len(myrs):
                        continue
                    x = E[i, : len(years) - L]
                    y = E[i + L, L:]
                    if np.isfinite(x).all() and np.isfinite(y).all():
                        cs.append(np.corrcoef(x, y)[0, 1])
                rows.append({"detrended": detrend, "epoch": ep, "lag": L, "corr": np.mean(cs), "n_masks": len(cs)})
    df = pd.DataFrame(rows)
    df.to_csv(ROOT / "coverage_persistence.csv", index=False, float_format="%.4f")
    piv = df.pivot_table(index="lag", columns=["detrended", "epoch"], values="corr")
    print(piv.loc[[1, 2, 3, 5, 10, 15, 20]].round(3).to_string())

    fits = {}
    for detrend in (False, True):
        c = df[(df.detrended == detrend) & (df.epoch == "all") & (df.lag >= FIT_LAGS[0]) & (df.lag <= FIT_LAGS[1])]
        L, r = c["lag"].values, c["corr"].values

        def f(par):
            p, phi = par
            return (1 - p) * phi ** L + p - r
        sol = least_squares(f, x0=[0.1, 0.5], bounds=([0, 0], [1, 0.999]))
        fits["detrended" if detrend else "raw"] = {"p": float(sol.x[0]), "phi": float(sol.x[1])}
    fits["land_definition"] = best
    fits["era5_years_lag"] = list(LAG_YEARS)
    fits["era5_years_sigma_check"] = list(REPRO_YEARS)
    print("\nfit corr(L) = (1-p) phi^L + p:", json.dumps(fits, indent=1))
    (ROOT / "coverage_persistence_fit.json").write_text(json.dumps(fits, indent=1))


if __name__ == "__main__":
    main()

"""Reanalysis land leaves (ERA5, JRA-3Q) and their 1995/96 splice into the LSAT ensemble.

Land series (annual means of monthly land means; absolute K, rebaselined later):
  ERA5    HRES monthly 2 m T, cos(lat) x fractional land-sea mask, from
          Land Data/ERA5/era5_t2m_5deg_monthly.nc (prepare_era5_5deg.py).
  JRA-3Q  2 m T ANALYSIS (JMA's own choice over land), cos(lat) x fractional
          land mask, from Land Data/JRA-3Q/ (prepare_jra3q_land.py).
          field="fc" uses the background instead (sensitivity test).
Both cover all land, Antarctica and Greenland included (domain="all"), or
exclude Antarctica (domain="ex_antarctica", boxes south of 60°S dropped), or are
subsampled to CRUTEM5's monthly coverage (mask="crutem5"; boxes are first converted
to anomalies vs their own 1981-2010 monthly climatology, since coverage varies).

Uncertainty: ERA5 EDA (10 members incl. the control) land means from
Land Data/ERA5/era5_eda_land_global_monthly.csv (prepare_era5_eda_land.py).
Perturbation_k = EDA_k - mean over members (annual). ERA5 member k = HRES +
perturbation_k; JRA-3Q member k = JRA-3Q + perturbation_k (ERA5-EDA donor, unscaled,
as in the GSAT ensemble). Optional persistent structural term for sensitivity tests.

Splice (member replacement): the existing 10,000 members are built unchanged; a
stratified 1/6 of them (exactly 1/6 of each primary leaf's members, largest-remainder
rounding; systematic sample ordered by each member's 1850-1900 offset) have their
values from SPLICE_YEAR onward replaced by a reanalysis member
(ERA5 and JRA-3Q split as evenly as possible, EDA member index balanced). Both
segments are zero-mean over the alignment window (1981-2010) before the splice; the
spliced member is then re-zeroed over 1981-2010, so every member keeps the file's
1981-2010 convention and pre-splice values relative to 1850-1900 are unchanged.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

ROOT = Path(__file__).resolve().parent
ERA5_5DEG = ROOT / "Land Data" / "ERA5" / "era5_t2m_5deg_monthly.nc"
EDA_CSV = ROOT / "Land Data" / "ERA5" / "era5_eda_land_global_monthly.csv"
JRA_CSV = ROOT / "Land Data" / "JRA-3Q" / "jra3q_land_monthly.csv"
JRA_5DEG = ROOT / "Land Data" / "JRA-3Q" / "jra3q_land_5deg_monthly.nc"
CRUTEM_NC = ROOT / "Land Data" / "CRUTEM5" / "CRUTEM.5.2.0.0.anomalies.nc"
N_EDA = 10
REAN_SEED = 20261125
FIRST_FULL_YEAR = {"era5": 1940, "jra3q": 1948}


def _annual(monthly: pd.Series) -> pd.Series:
    g = monthly.groupby(monthly.index.year)
    a = g.mean()
    return a[g.count() == 12]


def _crutem_mask(t: pd.DatetimeIndex) -> np.ndarray:
    ds = xr.open_dataset(CRUTEM_NC)
    tc = pd.DatetimeIndex(ds.time.values)
    ok = np.isfinite(ds.tas.values)
    lat = ds.latitude.values
    ds.close()
    if lat[0] > lat[-1]:
        ok = ok[:, ::-1]
    idx = {(d.year, d.month): i for i, d in enumerate(tc)}
    return np.stack([ok[idx[(d.year, d.month)]] if (d.year, d.month) in idx else np.zeros(ok.shape[1:], bool)
                     for d in t])


def _box_anomalies(T, t: pd.DatetimeIndex, base=(1981, 2010)):
    """Per-box, per-calendar-month anomalies vs `base`. Needed whenever the set of
    boxes changes over time (coverage masks); absolute-K means would jump with it."""
    A = np.full_like(T, np.nan, dtype=np.float64)
    inb = (t.year >= base[0]) & (t.year <= base[1])
    for m in range(1, 13):
        sel = t.month == m
        A[sel] = T[sel] - np.nanmean(T[sel & inb], axis=0)
    return A


def _box_mean(T, W, lat, domain, mask):
    W = np.broadcast_to(W, T.shape).copy()
    if domain == "ex_antarctica":
        W[:, lat < -60, :] = 0
    if mask is not None:
        W = W * mask
    ok = np.isfinite(T) & (W > 0)
    return np.where(ok, T * W, 0).sum((1, 2)) / np.where(ok, W, 0).sum((1, 2))


def era5_land_annual(domain="all", mask=None) -> pd.Series:
    ds = xr.open_dataset(ERA5_5DEG)
    t = pd.DatetimeIndex(ds.time.values)
    lat = ds.latitude.values
    W = (np.cos(np.deg2rad(lat))[:, None] * ds.landfrac.values)[None]
    m = _crutem_mask(t) if mask == "crutem5" else None
    T = ds.t2m_land.values if m is None else _box_anomalies(ds.t2m_land.values, t)
    s = pd.Series(_box_mean(T, W, lat, domain, m), index=t)
    ds.close()
    return _annual(s)


def jra3q_land_annual(field="an", domain="all", mask=None) -> pd.Series:
    if domain == "all" and mask is None:
        d = pd.read_csv(JRA_CSV)
        s = pd.Series(d[f"land_{field}_K"].values,
                      index=pd.to_datetime(dict(year=d.year, month=d.month, day=1)))
        return _annual(s)
    ds = xr.open_dataset(JRA_5DEG)
    t = pd.DatetimeIndex(ds.time.values)
    lat = ds.latitude.values
    m = _crutem_mask(t) if mask == "crutem5" else None
    T = ds[f"{field}_land"].values if m is None else _box_anomalies(ds[f"{field}_land"].values, t)
    s = pd.Series(_box_mean(T, ds.landweight.values[None], lat, domain, m), index=t)
    ds.close()
    return _annual(s)


def eda_annual(col="land_K") -> pd.DataFrame:
    """(year x member) annual means of the EDA land (or global) series."""
    d = pd.read_csv(EDA_CSV)
    p = d.pivot_table(index=["year", "month"], columns="member", values=col)
    n = p.groupby(level=0).count().min(axis=1)
    a = p.groupby(level=0).mean()
    return a[n == 12]


def eda_perturbations(years: np.ndarray) -> np.ndarray:
    a = eda_annual("land_K")
    if a.shape[1] != N_EDA:
        raise ValueError(f"expected {N_EDA} EDA members, got {a.shape[1]}")
    pert = a.sub(a.mean(axis=1), axis=0).reindex(years)
    return pert.values                                   # (n_years, 10); NaN before 1940


def structural_noise(n_years, n, sigma, tau, rng) -> np.ndarray:
    """Persistent AR(1) noise with stationary SD sigma and e-folding time tau (yr)."""
    phi = np.exp(-1.0 / tau)
    z = np.empty((n_years, n))
    z[0] = rng.standard_normal(n)
    for i in range(1, n_years):
        z[i] = phi * z[i - 1] + np.sqrt(1 - phi ** 2) * rng.standard_normal(n)
    return sigma * z


def reanalysis_leaves(years, baseline=(1981, 2010), jra_field="an", domain="all", mask=None,
                      struct_sigma=0.0, struct_tau=30.0, n_struct=100, which=("era5", "jra3q")):
    """Dict name -> (n_years, n_members) anomalies, zero-mean over `baseline`.
    Without a structural term each leaf has the 10 EDA-perturbation members."""
    pert = eda_perturbations(years)
    best = {"era5": era5_land_annual(domain, mask), "jra3q": jra3q_land_annual(jra_field, domain, mask)}
    rng = np.random.default_rng(REAN_SEED + 1)
    out = {}
    for name in which:
        b = best[name].reindex(years).values
        b[years < FIRST_FULL_YEAR[name]] = np.nan
        ens = b[:, None] + pert
        if struct_sigma > 0:
            k = np.arange(n_struct) % N_EDA
            ens = ens[:, k] + structural_noise(len(years), n_struct, struct_sigma, struct_tau, rng)
        base = (years >= baseline[0]) & (years <= baseline[1])
        out[name] = ens - np.nanmean(ens[base], axis=0, keepdims=True)
    return out


def _largest_remainder(counts: np.ndarray, frac: float) -> np.ndarray:
    target = int(round(counts.sum() * frac))
    raw = counts * frac
    n = np.floor(raw).astype(int)
    order = np.argsort(-(raw - n))
    n[order[: target - n.sum()]] += 1
    return n


def splice(final, years, primary_leaf, rean, weight, split_year=1996, baseline=(1981, 2010),
           align=(1981, 2010), seed=REAN_SEED):
    """Replace years >= split_year of a stratified `weight` fraction of members.
    Returns (spliced copy of final, manifest DataFrame of replaced members)."""
    out = final.copy()
    names = list(rean)
    if weight <= 0 or not names:
        return out, pd.DataFrame(columns=["member", "head_leaf", "head_member"])
    rng = np.random.default_rng(seed)
    leaves_u, counts = np.unique(primary_leaf, return_counts=True)
    n_take = _largest_remainder(counts, weight)
    # Systematic sample within each leaf, ordered by the member's early-baseline offset
    # (mean over 1850-1900 on the 1981-2010 baseline), so the replaced members span the
    # baseline distribution evenly and the headline quantiles are stable across seeds.
    key = np.nanmean(out[(years >= 1850) & (years <= 1900)], axis=0)
    chosen = []
    for l, k in zip(leaves_u, n_take):
        idx = np.flatnonzero(primary_leaf == l)
        idx = idx[np.argsort(key[idx], kind="stable")]
        pos = np.floor((np.arange(k) + rng.random()) * len(idx) / k).astype(int)
        chosen.append(idx[pos])
    chosen = np.concatenate(chosen)
    chosen = rng.permutation(chosen)
    head = np.array([names[i % len(names)] for i in range(len(chosen))])
    head = rng.permutation(head)
    head_k = np.empty(len(chosen), int)
    for nm in names:
        sel = np.flatnonzero(head == nm)
        nmem = rean[nm].shape[1]
        head_k[sel] = rng.permutation(np.arange(len(sel)) % nmem)
    a = (years >= align[0]) & (years <= align[1])
    b = (years >= baseline[0]) & (years <= baseline[1])
    late = years >= split_year
    for m, nm, k in zip(chosen, head, head_k):
        tail = out[:, m] - np.nanmean(out[a, m])
        r = rean[nm][:, k] - np.nanmean(rean[nm][a, k])
        x = np.where(late, r, tail)
        if not np.isfinite(x[late]).all():
            raise ValueError(f"{nm} member {k} has NaN after {split_year}")
        out[:, m] = x - np.nanmean(x[b])
    man = pd.DataFrame({"member": chosen, "head_leaf": head, "head_member": head_k}).sort_values("member")
    return out, man

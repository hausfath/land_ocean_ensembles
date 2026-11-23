"""Annual pseudo-ensembles for CRU-lineage LSAT products (CRUTEM5, GloSATLAT)
from their published ANNUAL global uncertainty components.

Error model per member k and year t (all terms independent):

  uncorrelated  u(t) * w_k(t)                 w white in time
  correlated    c(t) * x_k(t)                 x AR(1), e-folding TAU_CORR years
                                              (station homogenisation breaks,
                                              Morice et al. 2012: 40-yr repeat rate)
  coverage      v(t) * y_k(t)                 y standardised, corr(L) = (1-p) phi^L + p
                                              (p, phi from estimate_coverage_persistence.py)
  urbanisation  -max(0, Z_k) * R * d(t)       R = 0.0055 C/decade (1 sigma),
                                              d(t) = (max(t, 1900) - 1975.5)/10 decades;
                                              Brohan et al. 2006 / Morice et al. 2012:
                                              one-sided, negative draws set to zero
  exposure      split-normal(s_lo(t), s_hi(t)) * one draw per member
                                              s_lo, s_hi fitted so that urbanisation +
                                              exposure reproduces the published bias
                                              2.5/97.5% limits each year

Standardised processes have unit variance, so each component's annual SD on the
product's native baseline equals the published annual sigma.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import norm

URBAN_RATE = 0.0055        # C/decade, 1 sigma (Brohan et al. 2006)
URBAN_PIVOT = 1975.5       # centre of the 1961-1990 normal period
URBAN_START = 1900         # ramp flat before 1900
TAU_CORR = 40.0            # years


def urban_ramp(years: np.ndarray) -> np.ndarray:
    return (np.maximum(years, URBAN_START) - URBAN_PIVOT) / 10.0


def fit_exposure(years, lo, hi, n=400_000, seed=0, iters=60):
    """Per-year split-normal exposure scales (s_lo, s_hi) such that the 2.5/97.5%
    quantiles of urban + exposure offsets equal -lo and +hi (published bias limits
    relative to the best estimate). Returns (s_lo, s_hi, achieved_lo, achieved_hi)."""
    rng = np.random.default_rng(seed)
    z1 = np.maximum(0.0, rng.standard_normal(n))
    z2 = rng.standard_normal(n)
    d = urban_ramp(years)
    s_lo = np.zeros(len(years)); s_hi = np.zeros(len(years))
    a_lo = np.zeros(len(years)); a_hi = np.zeros(len(years))
    for i in range(len(years)):
        if not (np.isfinite(lo[i]) and np.isfinite(hi[i])):
            s_lo[i] = s_hi[i] = a_lo[i] = a_hi[i] = np.nan
            continue
        urb = -z1 * URBAN_RATE * d[i]
        ql, qh = np.quantile(urb, [0.025, 0.975])
        # exposure needed only where the urban term alone does not reach the limits
        sl = max(0.0, lo[i] ** 2 - ql ** 2) ** 0.5 / 1.96 if -ql < lo[i] - 1e-4 else 0.0
        sh = max(0.0, hi[i] ** 2 - qh ** 2) ** 0.5 / 1.96 if qh < hi[i] - 1e-4 else 0.0
        for _ in range(iters if (sl > 0 or sh > 0) else 0):
            e = np.where(z2 < 0, z2 * sl, z2 * sh)
            q = np.quantile(urb + e, [0.025, 0.975])
            if sl > 0:
                sl *= (lo[i] / max(-q[0], 1e-9)) ** 0.7
            if sh > 0:
                sh *= (hi[i] / max(q[1], 1e-9)) ** 0.7
        e = np.where(z2 < 0, z2 * sl, z2 * sh)
        q = np.quantile(urb + e, [0.025, 0.975])
        s_lo[i], s_hi[i], a_lo[i], a_hi[i] = sl, sh, -q[0], q[1]
    return s_lo, s_hi, a_lo, a_hi


def stratified_normal(n, rng):
    """n standard-normal draws at evenly spaced quantiles, randomly permuted.
    Used for the one-draw-per-member terms (urbanisation, exposure), so a
    200-member leaf has the intended mean offset instead of a sampled one."""
    return rng.permutation(norm.ppf((np.arange(n) + 0.5) / n))


def standard_ar1(n_years, n_members, phi, rng):
    out = np.empty((n_years, n_members))
    out[0] = rng.standard_normal(n_members)
    s = np.sqrt(1 - phi ** 2)
    for t in range(1, n_years):
        out[t] = phi * out[t - 1] + s * rng.standard_normal(n_members)
    return out


def standard_cov_process(n_years, n_members, p, phi, rng):
    """unit-variance process with corr(L) = (1-p) phi^L + p"""
    persist = rng.standard_normal(n_members)[None, :]
    if phi >= 1.0:
        ar = np.repeat(rng.standard_normal(n_members)[None, :], n_years, axis=0)
    else:
        ar = standard_ar1(n_years, n_members, phi, rng)
    return np.sqrt(p) * persist + np.sqrt(1 - p) * ar


def synthesize(comp: pd.DataFrame, n_members: int, rng: np.random.Generator,
               cov_p: float, cov_phi: float, tau_corr: float = TAU_CORR,
               exposure=None, recentre: bool = False, return_parts: bool = False):
    """comp: index = year (contiguous), columns a, u, c, v, lo, hi (annual, published;
    lo/hi = distance of bias 2.5/97.5% limits from the best estimate, >= 0).
    Returns (n_years, n_members) members on the product's native baseline."""
    years = comp.index.values
    ny = len(years)
    fin = np.isfinite(comp.a.values)
    u = np.nan_to_num(comp.u.values); c = np.nan_to_num(comp.c.values); v = np.nan_to_num(comp.v.values)
    if exposure is None:
        exposure = fit_exposure(years, comp.lo.values, comp.hi.values)
    s_lo, s_hi = np.nan_to_num(exposure[0]), np.nan_to_num(exposure[1])

    parts = {}
    parts["uncorrelated"] = u[:, None] * rng.standard_normal((ny, n_members))
    parts["correlated"] = c[:, None] * standard_ar1(ny, n_members, np.exp(-1.0 / tau_corr), rng)
    parts["coverage"] = v[:, None] * standard_cov_process(ny, n_members, cov_p, cov_phi, rng)
    z1 = np.maximum(0.0, stratified_normal(n_members, rng))
    parts["urban"] = -(URBAN_RATE * urban_ramp(years))[:, None] * z1[None, :]
    z2 = stratified_normal(n_members, rng)[None, :]
    parts["exposure"] = np.where(z2 < 0, z2 * s_lo[:, None], z2 * s_hi[:, None])
    err = sum(parts.values())
    if recentre:
        err = err - err.mean(axis=1, keepdims=True)
    members = comp.a.values[:, None] + err
    members[~fin] = np.nan
    if return_parts:
        return members, parts
    return members


def with_central(comp: pd.DataFrame, central: pd.Series) -> pd.DataFrame:
    """Another product's best estimate with CRUTEM5's error model: sigma columns
    taken from `comp` (back/forward-filled where CRUTEM5 has no value)."""
    out = comp.reindex(central.index).copy()
    out[["u", "c", "v", "lo", "hi"]] = out[["u", "c", "v", "lo", "hi"]].bfill().ffill()
    out["a"] = central.values
    return out


def load_crutem5_annual(data_dir, version="5.2.0.0") -> pd.DataFrame:
    df = pd.read_csv(data_dir / "CRUTEM5" / f"CRUTEM.{version}.component_series.global.annual.csv")
    df.columns = ["year", "a", "u", "c", "v", "blo", "bhi"]
    df = df.set_index("year")
    df["lo"] = (df.a - df.blo).clip(lower=0)
    df["hi"] = (df.bhi - df.a).clip(lower=0)
    return df[["a", "u", "c", "v", "lo", "hi"]]


def load_glosatlat_annual(data_dir) -> pd.DataFrame:
    import xarray as xr
    ds = xr.open_dataset(data_dir / "GloSATLAT" / "GloSATLAT-1-0-0-0_component-series_global_annual.nc")
    df = pd.DataFrame({
        "a": ds.tas.values, "u": ds.tas_unc.values, "c": ds.tas_corr.values,
        "v": ds.coverage_unc.values,
        "lo": ds.tas.values - ds.tas_lower_bias.values,
        "hi": ds.tas_upper_bias.values - ds.tas.values,
    }, index=pd.DatetimeIndex(ds.time.values).year)
    df.index.name = "year"
    df[["lo", "hi"]] = df[["lo", "hi"]].clip(lower=0)
    return df

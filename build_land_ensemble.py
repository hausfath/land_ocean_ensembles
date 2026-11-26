"""
Build a 10,000-member ensemble of annual global-mean land surface air
temperature (LSAT) anomaly, 1850-2025, by combining six LSAT products
through a structural-method family tree with donor-imputed uncertainty
where native ensembles are not available.

See `land_ensemble_methodology.md` for full description.

Inputs (in `Land Data/`):
    Berkeley Earth Highres/Land_TAVG_ensemble.txt          (10-member monthly)
    CRUTEM5/CRUTEM.5.2.0.0.component_series.global.annual.csv   (annual σ components)
    GloSATLAT/GloSATLAT-1-0-0-0_component-series_global_annual.nc
    (monthly component files: --method legacy only)
    ../coverage_persistence_fit.json                       (estimate_coverage_persistence.py)
CRU-lineage synthesis: lsat_synthesis.py. NOAA Land and C-LSAT 2.1 get CRUTEM5's
error model around their own best estimates.
    DCLSAT/DCLSAT_I_monthly_ensemble.csv                   (200-member monthly,
                                                            from pull_dcent_lsat_members.py)
    NOAA Land/aravg.mon.land.90S.90N.v6.1.0.202604.asc

Outputs (in the project root, alongside this script):
    land_ensemble.csv                  176 rows × 10,001 cols
                                         (year, m00001 .. m10000)
    land_ensemble_summary.csv          per-year mean, 2.5/50/97.5 pct, SD
    land_ensemble_perdataset.csv       per-dataset 200-or-10-member ensembles
                                         used as inputs (annual, baselined)
    temp_files/                        intermediate per-dataset CSVs
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

import reanalysis_land as R

ROOT = Path(__file__).parent
DATA = ROOT / "Land Data"
TMP = ROOT / "temp_files"
TMP.mkdir(exist_ok=True)

# ---------- configuration ----------
START_YEAR = 1850
END_YEAR = 2025                         # full year completed
YEARS = np.arange(START_YEAR, END_YEAR + 1)
N_YEARS = len(YEARS)
BASELINE = (1981, 2010)
N_FINAL = 10_000
N_SYNTH = 200                           # per-dataset synthesized members
RHO_MONTHLY = 0.3                       # AR(1) coefficient for uncorrelated term
RNG = np.random.default_rng(20260512)   # mixture sampling (and the legacy synthesis)
SYNTH_SEED = 20261123                   # own stream for the CRU-lineage synthesis
CRU_VERSION = "5.2.0.0"
SMOOTH_YEARS = 21                       # running window for the donor scale factor
# coverage-error persistence variants: corr(L) = (1-p) phi^L + p
COVERAGE_VARIANTS = {
    "era5":       None,                 # fitted (raw) by estimate_coverage_persistence.py
    "detrended":  None,                 # fitted to detrended errors
    "white":      (0.0, 0.0),
    "ar1_4yr":    (0.0, float(np.exp(-1 / 4))),
    "persistent": (1.0, 0.0),
}


# ---------- helpers ----------
def reduce_monthly_to_annual(months_array: np.ndarray, year_idx: np.ndarray) -> np.ndarray:
    """Average monthly values within each calendar year.
    months_array: shape (n_months, n_members)
    year_idx:     shape (n_months,) integer year per month
    Returns: shape (N_YEARS, n_members) on the common YEARS axis. Years
             with <12 valid months are NaN.
    """
    n_members = months_array.shape[1] if months_array.ndim == 2 else 1
    arr2d = months_array.reshape(-1, n_members) if months_array.ndim == 1 else months_array
    out = np.full((N_YEARS, n_members), np.nan)
    for i, y in enumerate(YEARS):
        mask = year_idx == y
        if mask.sum() == 12:
            blk = arr2d[mask]
            if np.isfinite(blk).all():
                out[i] = blk.mean(axis=0)
    return out


def rebaseline_to_1981_2010(ensemble: np.ndarray) -> np.ndarray:
    """Subtract each member's 1981-2010 mean."""
    base_mask = (YEARS >= BASELINE[0]) & (YEARS <= BASELINE[1])
    means = np.nanmean(ensemble[base_mask], axis=0, keepdims=True)
    return ensemble - means


def ar1_noise(sigma_t: np.ndarray, n_members: int, rho: float, rng: np.random.Generator) -> np.ndarray:
    """AR(1) noise with time-varying marginal SD sigma_t.
    sigma_t: shape (n_months,)
    returns: shape (n_months, n_members)
    """
    n = sigma_t.shape[0]
    out = np.zeros((n, n_members))
    # initialize at marginal distribution
    out[0] = rng.standard_normal(n_members) * sigma_t[0]
    rho2 = np.sqrt(max(0.0, 1.0 - rho ** 2))
    for t in range(1, n):
        # scale previous to current sigma, plus innovation
        scale = sigma_t[t] / sigma_t[t - 1] if sigma_t[t - 1] > 0 else 0.0
        innov = rng.standard_normal(n_members) * sigma_t[t]
        out[t] = rho * scale * out[t - 1] + rho2 * innov
    return out


# ---------- per-dataset loaders ----------
def load_berkeley_earth() -> np.ndarray:
    """Returns (N_YEARS, 10) annual ensemble. Native baseline 1951-1980."""
    path = DATA / "Berkeley Earth Highres" / "Land_TAVG_ensemble.txt"
    rows: list[list[float]] = []
    with open(path) as f:
        for line in f:
            if line.startswith("%") or not line.strip():
                continue
            parts = line.split()
            if len(parts) >= 12:
                rows.append([float(p) for p in parts[:12]])
    arr = np.array(rows)
    years = arr[:, 0].astype(int)
    months = arr[:, 1].astype(int)
    members = arr[:, 2:12]  # 10 members
    # build a (n_months_total, 10) consistent matrix aligned by year/month
    annual = reduce_monthly_to_annual(members, years)
    return annual


def load_crutem5_components(version: str = "5.1.0.0") -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load CRUTEM5 monthly anomaly + four uncertainty components on the 1850-2025
    monthly axis. Returns (years_idx, anomaly, u_unc, u_corr, u_cov, bias_sigma)."""
    df = pd.read_csv(DATA / "CRUTEM5" / f"CRUTEM.{version}.component_series.global.monthly.csv")
    df["Time"] = pd.to_datetime(df["Time"])
    df["year"] = df["Time"].dt.year
    df["month"] = df["Time"].dt.month
    df = df[(df["year"] >= START_YEAR) & (df["year"] <= END_YEAR)]
    bias_sigma = (df["Bias component upper confidence limit (97.5%)"]
                  - df["Bias component lower confidence limit (2.5%)"]) / (2 * 1.96)
    return (
        df["year"].values,
        df["Anomaly (deg C)"].values,
        df["Uncorrelated uncertainty (1 sigma)"].values,
        df["Correlated uncertainty (1 sigma)"].values,
        df["Coverage uncertainty (1 sigma)"].values,
        bias_sigma.values,
    )


def load_glosat_components() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load GloSATLAT monthly anomaly + four uncertainty components on the
    1850-2025 monthly axis. Years beyond GloSAT's 2021-12 end are filled NaN.
    Returns (year_idx, anomaly, u_unc, u_corr, u_cov, bias_sigma)."""
    ds = xr.open_dataset(
        DATA / "GloSATLAT" / "GloSATLAT-1-0-0-0_component-series_global_monthly.nc"
    )
    t = pd.DatetimeIndex(ds["time"].values)
    yrs = t.year.values
    mask = (yrs >= START_YEAR) & (yrs <= END_YEAR)
    bias_sigma = (ds["tas_upper_bias"].values - ds["tas_lower_bias"].values) / (2 * 1.96)
    return (
        yrs[mask],
        ds["tas"].values[mask],
        ds["tas_unc"].values[mask],
        ds["tas_corr"].values[mask],
        ds["coverage_unc"].values[mask],
        bias_sigma[mask],
    )


def synthesize_component_ensemble(
    year_idx: np.ndarray,
    anomaly: np.ndarray,
    u_unc: np.ndarray,
    u_corr: np.ndarray,
    u_cov: np.ndarray,
    bias_sigma: np.ndarray,
    n_members: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Build (N_YEARS, n_members) annual ensemble from monthly central + 4 σ.
    Missing months (NaN in anomaly) are skipped on annual aggregation.
    """
    n_months = anomaly.shape[0]
    # uncorrelated AR(1) noise
    u_unc_safe = np.where(np.isfinite(u_unc), u_unc, 0.0)
    e_unc = ar1_noise(u_unc_safe, n_members, RHO_MONTHLY, rng)
    # correlated, coverage, bias = scalar × time-varying σ
    alpha = rng.standard_normal(n_members)
    beta = rng.standard_normal(n_members)
    gamma = rng.standard_normal(n_members)
    e_corr = np.where(np.isfinite(u_corr), u_corr, 0.0)[:, None] * alpha[None, :]
    e_cov = np.where(np.isfinite(u_cov), u_cov, 0.0)[:, None] * beta[None, :]
    e_bias = np.where(np.isfinite(bias_sigma), bias_sigma, 0.0)[:, None] * gamma[None, :]
    central = np.where(np.isfinite(anomaly), anomaly, 0.0)[:, None]
    monthly_members = central + e_unc + e_corr + e_cov + e_bias
    # mask months that were NaN in central back to NaN
    bad = ~np.isfinite(anomaly)
    monthly_members[bad, :] = np.nan
    annual = reduce_monthly_to_annual(monthly_members, year_idx)
    return annual


def load_dclsat_ensemble() -> np.ndarray:
    """Returns (N_YEARS, 200) annual ensemble. Native baseline 1982-2014.
    Backed by DCENT-I v1.1.0.0 (Chan et al. 2026, GDJ) — the spatially
    complete kriging-infilled extension of DCENT. Each member's native
    `lsat` field has been area-weighted to a global mean by
    `pull_dcent_i_members.py` (see methodology §2.4).
    """
    path = DATA / "DCLSAT" / "DCLSAT_I_monthly_ensemble.csv"
    df = pd.read_csv(path)
    member_cols = [c for c in df.columns if c.startswith("m") and c != "month"]
    assert len(member_cols) == 200, f"expected 200 DCLSAT members, got {len(member_cols)}"
    mask = (df["year"] >= START_YEAR) & (df["year"] <= END_YEAR)
    df = df[mask]
    monthly = df[member_cols].values  # (months, 200)
    return reduce_monthly_to_annual(monthly, df["year"].values)


def donor_rescale(target_sigma: np.ndarray, donor_sigma: np.ndarray, smooth: int = 0) -> np.ndarray:
    """Per-year factor that scales the donor's member deviations to the target σ.
    smooth > 1: centred running mean (in log space) over `smooth` years, so the
    donor's year-to-year shape is not distorted by year-specific factors."""
    r = np.where(donor_sigma > 0, target_sigma / donor_sigma, np.nan)
    if smooth and smooth > 1:
        lr = pd.Series(np.log(r)).rolling(smooth, center=True, min_periods=smooth // 2 + 1).mean()
        r = np.exp(lr.values)
    return np.nan_to_num(r, nan=0.0)


def load_noaa_land(donor_ensemble: np.ndarray, target_sigma: np.ndarray, rng: np.random.Generator,
                   smooth: int = 0) -> np.ndarray:
    """NOAA Land best estimate + donor-imputed uncertainty.
    donor_ensemble: (N_YEARS, n_donor) annual, already baselined.
    target_sigma: (N_YEARS,) target annual SD time series (from GloSATLAT, extended).
    Returns (N_YEARS, n_donor) annual ensemble.
    """
    path = DATA / "NOAA Land" / "aravg.mon.land.90S.90N.v6.1.0.202604.asc"
    rows = []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            yr = int(parts[0])
            mo = int(parts[1])
            anom = float(parts[2])
            rows.append((yr, mo, anom))
    df = pd.DataFrame(rows, columns=["year", "month", "anom"])
    df = df[(df["year"] >= START_YEAR) & (df["year"] <= END_YEAR)]
    annual_central = reduce_monthly_to_annual(df["anom"].values[:, None], df["year"].values)[:, 0]
    # demean donor at each year and rescale to target_sigma / donor's own SD
    donor_demeaned = donor_ensemble - np.nanmean(donor_ensemble, axis=1, keepdims=True)
    donor_sigma = np.nanstd(donor_ensemble, axis=1, ddof=1)
    rescale = donor_rescale(target_sigma, donor_sigma, smooth)[:, None]
    donor_scaled = donor_demeaned * rescale
    return annual_central[:, None] + donor_scaled


def load_c_lsat(donor_ensemble: np.ndarray, target_sigma: np.ndarray, rng: np.random.Generator,
                smooth: int = 0) -> np.ndarray:
    """C-LSAT 2.1 (Sun et al. / Li, Sun Yat-sen University; CMA-homogenization
    lineage) best estimate + DCLSAT-donor-imputed uncertainty.

    Mirrors `load_noaa_land`. The caller MUST pass a DISJOINT slice of the
    DCLSAT donor relative to NOAA Land (NOAA Land draws from members 1-100,
    C-LSAT 2.1 from 101-200), so the two donor-imputed LSAT leaves do not
    share per-member trajectory shape — addresses the stats-review concern
    that two of six leaves would otherwise carry shape-correlated envelopes.
    """
    path = DATA / "C-LSAT" / "C-LSAT2.1_global_monthly.csv"
    df = pd.read_csv(path)
    df = df[(df["year"] >= START_YEAR) & (df["year"] <= END_YEAR)]
    annual_central = reduce_monthly_to_annual(df["anomaly"].values[:, None], df["year"].values)[:, 0]
    donor_demeaned = donor_ensemble - np.nanmean(donor_ensemble, axis=1, keepdims=True)
    donor_sigma = np.nanstd(donor_ensemble, axis=1, ddof=1)
    rescale = donor_rescale(target_sigma, donor_sigma, smooth)[:, None]
    donor_scaled = donor_demeaned * rescale
    return annual_central[:, None] + donor_scaled


def target_sigma_from_glosat(glosat_ensemble: np.ndarray, crutem_ensemble: np.ndarray) -> np.ndarray:
    """Build a 1850-2025 annual target σ time series, using GloSATLAT's spread
    where available (1850-2021) and CRUTEM5's spread for 2022-2025. For the
    pre-coverage years where both are NaN, fall back to the earliest finite
    σ (so the donor envelope persists rather than collapsing to zero)."""
    sigma_g = np.nanstd(glosat_ensemble, axis=1, ddof=1)
    sigma_c = np.nanstd(crutem_ensemble, axis=1, ddof=1)
    out = np.where(np.isfinite(sigma_g) & (sigma_g > 0), sigma_g, sigma_c)
    # back/forward fill any remaining NaN with nearest finite value
    s = pd.Series(out)
    out = s.bfill().ffill().values
    return out


# ---------- family tree ----------
# Five method families, each at 1/5 weight (Thorne et al. 2026 "equal-weight
# method families" principle):
#   - Berkeley Earth   (HOMOG/SCALPEL)
#   - NOAA Land        (PHA)
#   - CRU-lineage      (CRUTEM5 1/10 + GloSATLAT 1/10 = 1/5)
#   - DCLSAT           (Dynamical-Constraint)
#   - C-LSAT 2.1       (CMA-homogenization)
TREE_FULL = {
    "berkeley_earth": 0.20,
    "noaa_land":      0.20,
    "dclsat":         0.20,
    "crutem5":        0.10,
    "glosatlat":      0.10,
    "c_lsat":         0.20,
}
# For 2022-2025 (GloSATLAT absent), redistribute its 0.10 to CRUTEM5 (the
# other CRU-lineage leaf) so CRU-lineage family total stays 0.20.
TREE_NO_GLOSAT = {
    "berkeley_earth": 0.20,
    "noaa_land":      0.20,
    "dclsat":         0.20,
    "crutem5":        0.20,
    "glosatlat":      0.0,
    "c_lsat":         0.20,
}


def leaf_probs_for_year(year: int) -> dict[str, float]:
    return TREE_NO_GLOSAT if year > 2021 else TREE_FULL


# ---------- main ----------
def coverage_params(variant: str) -> tuple[float, float]:
    if COVERAGE_VARIANTS[variant] is not None:
        return COVERAGE_VARIANTS[variant]
    fit = json.loads((ROOT / "coverage_persistence_fit.json").read_text())
    f = fit["raw" if variant == "era5" else "detrended"]
    return f["p"], f["phi"]


def cru_lineage_annual(cov_variant: str, recentre: bool, n: int = N_SYNTH):
    """CRUTEM5 and GloSATLAT members from their published ANNUAL components
    (lsat_synthesis.py). Both use the same seed on the same 1850-2025 axis, so
    member k of each leaf carries the same standardised draws: a mixture member
    that switches from GloSATLAT to CRUTEM5 after 2021 keeps its error draws."""
    import lsat_synthesis as S
    p, phi = coverage_params(cov_variant)
    out = {}
    for name, comp in (("crutem5", S.load_crutem5_annual(DATA, CRU_VERSION)),
                       ("glosatlat", S.load_glosatlat_annual(DATA))):
        comp = comp.reindex(YEARS)
        out[name] = S.synthesize(comp, n, np.random.default_rng(SYNTH_SEED),
                                 cov_p=p, cov_phi=phi, recentre=recentre)
    return out["crutem5"], out["glosatlat"], (p, phi)


def donor_central(name: str) -> np.ndarray:
    """Annual best estimate (native baseline) of NOAA Land or C-LSAT 2.1."""
    if name == "noaa_land":
        path = DATA / "NOAA Land" / "aravg.mon.land.90S.90N.v6.1.0.202604.asc"
        rows = [tuple(map(float, ln.split()[:3])) for ln in open(path) if ln.split()]
        df = pd.DataFrame(rows, columns=["year", "month", "anom"])
        col = "anom"
    else:
        df = pd.read_csv(DATA / "C-LSAT" / "C-LSAT2.1_global_monthly.csv")
        col = "anomaly"
    df = df[(df["year"] >= START_YEAR) & (df["year"] <= END_YEAR)]
    return reduce_monthly_to_annual(df[col].values[:, None], df["year"].values.astype(int))[:, 0]


def cru_model_leaves(cov_variant: str):
    """NOAA Land and C-LSAT 2.1 (no published uncertainty components): their own best
    estimate plus CRUTEM5's annual error model (lsat_synthesis.py), each with an
    independent seed and recentred so the leaf mean equals the product's best estimate
    (CRUTEM5's one-sided urbanisation term adds spread, not an offset, to products
    with their own urban/homogenisation treatment)."""
    import lsat_synthesis as S
    p, phi = coverage_params(cov_variant)
    cru = S.load_crutem5_annual(DATA, CRU_VERSION)
    out = {}
    for k, name in enumerate(("noaa_land", "c_lsat"), start=1):
        comp = S.with_central(cru, pd.Series(donor_central(name), index=YEARS))
        out[name] = S.synthesize(comp, N_SYNTH, np.random.default_rng(SYNTH_SEED + k),
                                 cov_p=p, cov_phi=phi, recentre=True)
    return out


def parse_args(argv=None):
    ap = argparse.ArgumentParser(description="Build the LSAT ensemble.")
    ap.add_argument("--method", choices=["annual", "legacy"], default="annual",
                    help="CRU-lineage synthesis: annual components (default) or the "
                         "pre-2026-11 monthly whole-record-scalar recipe")
    ap.add_argument("--cru-version", default=CRU_VERSION, help="CRUTEM5 version for --method legacy")
    ap.add_argument("--coverage", choices=list(COVERAGE_VARIANTS), default="era5")
    ap.add_argument("--recentre", action="store_true",
                    help="subtract the member-mean error so the leaf mean equals the best estimate")
    ap.add_argument("--donor", choices=["cru_model", "dclsat"], default="cru_model",
                    help="NOAA Land / C-LSAT uncertainty: CRUTEM5 error model (default) or "
                         "DCLSAT member shapes rescaled to the CRU-lineage sigma")
    ap.add_argument("--smooth", type=int, default=SMOOTH_YEARS,
                    help="running window (yr) for the donor scale factor; 0 = per-year (legacy)")
    ap.add_argument("--reanalysis", choices=["era5_jra3q", "era5", "jra3q", "none"], default="era5_jra3q",
                    help="reanalysis family spliced in from --splice-year (none = pre-2026-11-25 build)")
    ap.add_argument("--rean-weight", type=float, default=1 / 6,
                    help="reanalysis family weight from the splice year on (default 1/6)")
    ap.add_argument("--splice-year", type=int, default=1996, help="first year taken from the reanalysis")
    ap.add_argument("--align", default="1981-2010", help="window on which tail and head are aligned")
    ap.add_argument("--jra-field", choices=["an", "fc"], default="an",
                    help="JRA-3Q land field: analysis (JMA method, default) or background")
    ap.add_argument("--rean-domain", choices=["all", "ex_antarctica"], default="all")
    ap.add_argument("--rean-mask", choices=["none", "crutem5"], default="none",
                    help="subsample the reanalyses to CRUTEM5 monthly coverage (sensitivity)")
    ap.add_argument("--rean-struct", type=float, default=0.0,
                    help="SD (degC) of a persistent AR(1) structural term added to the reanalysis leaves")
    ap.add_argument("--rean-struct-tau", type=float, default=30.0, help="e-folding time (yr) of that term")
    ap.add_argument("--rean-seed", type=int, default=None, help="replacement RNG seed (default reanalysis_land.REAN_SEED)")
    ap.add_argument("--out", type=Path, default=ROOT, help="output directory")
    return ap.parse_args(argv)


def build_leaves(args):
    """Per-dataset ensembles, each rebaselined to 1981-2010.
    Returns (leaves dict, donor target sigma, coverage (p, phi) or None, donor smoothing)."""
    legacy = args.method == "legacy"

    print("Loading Berkeley Earth (native 10 members)...")
    be = load_berkeley_earth()
    print(f"  shape: {be.shape}  finite years: {np.isfinite(be).all(axis=1).sum()}")

    if legacy:
        print(f"Legacy synthesis: CRUTEM {args.cru_version} + GloSATLAT, monthly components...")
        cru_in = load_crutem5_components(args.cru_version)
        cru = synthesize_component_ensemble(*cru_in, N_SYNTH, RNG)
        glo_in = load_glosat_components()
        glo = synthesize_component_ensemble(*glo_in, N_SYNTH, RNG)
        cov_pp = None
    else:
        # The legacy synthesis drew from RNG; consume the same number of draws so the
        # mixture's leaf/member assignment below is identical to the legacy build.
        n_months_cru = len(load_crutem5_components(args.cru_version)[0])
        n_months_glo = len(load_glosat_components()[0])
        RNG.standard_normal((n_months_cru + n_months_glo + 6) * N_SYNTH)
        print(f"Annual synthesis: CRUTEM {CRU_VERSION} + GloSATLAT, coverage={args.coverage}, "
              f"recentre={args.recentre}...")
        cru, glo, cov_pp = cru_lineage_annual(args.coverage, args.recentre)
    print(f"  CRUTEM5 {cru.shape} finite years {np.isfinite(cru).all(axis=1).sum()}; "
          f"GloSATLAT finite years {np.isfinite(glo).all(axis=1).sum()}")

    print("Loading DCLSAT (200 native members)...")
    dcl = load_dclsat_ensemble()
    print(f"  shape: {dcl.shape}  finite years: {np.isfinite(dcl).all(axis=1).sum()}")

    # baseline DCLSAT first so the donor demeaning is sensible; also baseline
    # CRUTEM/GloSAT so target_sigma_from_glosat sees them on the common axis.
    dcl_baselined_for_donor = rebaseline_to_1981_2010(dcl)
    cru_baselined_for_sigma = rebaseline_to_1981_2010(cru)
    glo_baselined_for_sigma = rebaseline_to_1981_2010(glo)
    tgt_sigma = target_sigma_from_glosat(glo_baselined_for_sigma, cru_baselined_for_sigma)

    # Two donor-imputed LSAT leaves (NOAA Land and C-LSAT 2.1) draw from
    # DISJOINT slices of DCLSAT's 200-member donor so they don't share per-
    # member trajectory shape. NOAA gets members 1-100; C-LSAT gets 101-200.
    print("Building NOAA Land and C-LSAT 2.1 leaves...")
    if legacy or args.donor == "dclsat":
        smooth = 0 if legacy else args.smooth
        noaa = load_noaa_land(dcl_baselined_for_donor[:, :100], tgt_sigma, RNG, smooth)
        c_lsat = load_c_lsat(dcl_baselined_for_donor[:, 100:], tgt_sigma, RNG, smooth)
    else:
        print("NOAA Land / C-LSAT 2.1: own best estimate + CRUTEM5 annual error model...")
        smooth = None
        d = cru_model_leaves(args.coverage)
        noaa, c_lsat = d["noaa_land"], d["c_lsat"]
    print(f"  NOAA Land finite years {np.isfinite(noaa).all(axis=1).sum()}; "
          f"C-LSAT finite years {np.isfinite(c_lsat).all(axis=1).sum()}")

    print("Re-baselining all per-dataset ensembles to 1981-2010...")
    be_b = rebaseline_to_1981_2010(be)
    cru_b = rebaseline_to_1981_2010(cru)
    glo_b = rebaseline_to_1981_2010(glo)
    dcl_b = rebaseline_to_1981_2010(dcl)
    noaa_b = rebaseline_to_1981_2010(noaa)
    c_lsat_b = rebaseline_to_1981_2010(c_lsat)

    leaves: dict[str, np.ndarray] = {
        "berkeley_earth": be_b,
        "crutem5": cru_b,
        "glosatlat": glo_b,
        "dclsat": dcl_b,
        "noaa_land": noaa_b,
        "c_lsat": c_lsat_b,
    }
    return leaves, tgt_sigma, cov_pp, smooth


def main(argv=None) -> None:
    args = parse_args(argv)
    out_dir = args.out
    out_dir.mkdir(parents=True, exist_ok=True)
    legacy = args.method == "legacy"
    leaves, tgt_sigma, cov_pp, smooth = build_leaves(args)

    # Sanity print: ensemble-mean per dataset for the latest year
    print("\nSanity: ensemble mean (best estimate) for 2024 (re-baselined to 1981-2010):")
    iy = int(np.where(YEARS == 2024)[0][0])
    for name, ens in leaves.items():
        print(f"  {name:16s}: {np.nanmean(ens[iy]):+.3f} °C  spread(1σ)={np.nanstd(ens[iy], ddof=1):.3f}")

    # Save per-dataset best estimates and SDs to a single CSV for diagnostics
    diag = pd.DataFrame({"year": YEARS})
    for name, ens in leaves.items():
        diag[f"{name}_mean"] = np.nanmean(ens, axis=1)
        diag[f"{name}_sd"]   = np.nanstd(ens, axis=1, ddof=1)
    diag["donor_target_sd"] = tgt_sigma
    diag.to_csv(out_dir / "land_ensemble_perdataset.csv", index=False, float_format="%.5f")
    print(f"\nWrote land_ensemble_perdataset.csv")

    # ---------- family-tree sampling to 10000 ----------
    print(f"\nSampling {N_FINAL} ensemble members via year-dependent leaf weights...")
    # Pre-build a per-year leaf-choice array of shape (N_FINAL,). We assign each
    # final member a single leaf for the WHOLE record (so a member is a coherent
    # dataset trajectory), but we still need to swap leaves for years where the
    # assigned leaf is GloSATLAT and the year is >2021 (GloSATLAT is absent).
    leaf_names = list(leaves.keys())
    leaf_to_idx = {n: i for i, n in enumerate(leaf_names)}
    probs_full = np.array([TREE_FULL[n] for n in leaf_names])

    final = np.full((N_YEARS, N_FINAL), np.nan)
    member_leaf_log = np.zeros((N_YEARS, N_FINAL), dtype=np.int8)

    # Choose primary leaf per final member, using the full-tree weights
    primary_leaf = RNG.choice(len(leaf_names), size=N_FINAL, p=probs_full)
    # Choose member index within each chosen leaf
    leaf_member_counts = {n: leaves[n].shape[1] for n in leaf_names}
    primary_member_idx = np.array(
        [RNG.integers(0, leaf_member_counts[leaf_names[primary_leaf[m]]]) for m in range(N_FINAL)]
    )

    # For GloSATLAT-assigned members on years >2021, draw a fallback
    # leaf from {crutem5} (the redistribute target). Cache fallback samples.
    glosat_idx = leaf_to_idx["glosatlat"]
    crutem_idx = leaf_to_idx["crutem5"]
    crutem_count = leaf_member_counts["crutem5"]
    fallback_member_idx = RNG.integers(0, crutem_count, size=N_FINAL)

    # Precompute, per year, the set of leaves whose ensemble has finite values
    # that year. If a member's primary leaf is NaN for a given year, we
    # redistribute among the finite leaves (renormalized weights), and draw a
    # fresh member index from a chosen finite leaf.
    finite_by_year = {}
    for iy, yr in enumerate(YEARS):
        finite_by_year[iy] = [
            i for i, n in enumerate(leaf_names) if np.isfinite(leaves[n][iy]).any()
        ]

    for iy, yr in enumerate(YEARS):
        finite_leaves = finite_by_year[iy]
        # renormalized weights over finite leaves only
        renorm_p = np.array([TREE_FULL[leaf_names[i]] for i in finite_leaves])
        if yr > 2021:
            # GloSATLAT effectively absent: rebuild via TREE_NO_GLOSAT
            renorm_p = np.array([TREE_NO_GLOSAT[leaf_names[i]] for i in finite_leaves])
        renorm_p = renorm_p / renorm_p.sum() if renorm_p.sum() > 0 else renorm_p

        for m in range(N_FINAL):
            pl = primary_leaf[m]
            mi = primary_member_idx[m]
            # year-dependent vetoes: GloSAT after 2021
            if (pl == glosat_idx and yr > 2021):
                pl = crutem_idx
                # annual method: same member index (shared draws, see cru_lineage_annual)
                mi = fallback_member_idx[m] if legacy else mi
            # NaN veto: if the chosen (leaf, member) is NaN, resample leaf
            if not np.isfinite(leaves[leaf_names[pl]][iy, mi]):
                pl = RNG.choice(finite_leaves, p=renorm_p)
                mi = RNG.integers(0, leaf_member_counts[leaf_names[pl]])
                # in case that drew a NaN (shouldn't if leaf has any finite), keep trying
                tries = 0
                while not np.isfinite(leaves[leaf_names[pl]][iy, mi]) and tries < 20:
                    mi = RNG.integers(0, leaf_member_counts[leaf_names[pl]])
                    tries += 1
            final[iy, m] = leaves[leaf_names[pl]][iy, mi]
            member_leaf_log[iy, m] = pl

    # ---------- reanalysis family: splice into a stratified share of members ----------
    manifest = pd.DataFrame({"member": np.arange(N_FINAL),
                             "primary_leaf": [leaf_names[i] for i in primary_leaf],
                             "primary_member": primary_member_idx})
    rean_cfg = None
    if args.reanalysis != "none" and args.rean_weight > 0:
        which = ("era5", "jra3q") if args.reanalysis == "era5_jra3q" else (args.reanalysis,)
        al = tuple(int(x) for x in args.align.split("-"))
        rean = R.reanalysis_leaves(YEARS, BASELINE, args.jra_field, args.rean_domain,
                                   None if args.rean_mask == "none" else args.rean_mask,
                                   args.rean_struct, args.rean_struct_tau, which=which)
        for nm, ens in rean.items():
            diag[f"{nm}_mean"] = np.nanmean(ens, axis=1)
            diag[f"{nm}_sd"] = np.nanstd(ens, axis=1, ddof=1)
        diag.to_csv(out_dir / "land_ensemble_perdataset.csv", index=False, float_format="%.5f")
        final, man_r = R.splice(final, YEARS, primary_leaf, rean, args.rean_weight, args.splice_year,
                                BASELINE, al, seed=args.rean_seed or R.REAN_SEED)
        manifest = manifest.merge(man_r, on="member", how="left")
        late = YEARS >= args.splice_year
        for nm in rean:
            member_leaf_log[np.ix_(late, man_r.member[man_r.head_leaf == nm].values)] = -1 - list(rean).index(nm)
        rean_cfg = {"leaves": list(which), "weight": args.rean_weight, "splice_year": args.splice_year,
                    "align": list(al), "jra_field": args.jra_field, "domain": args.rean_domain,
                    "mask": args.rean_mask, "struct_sigma": args.rean_struct,
                    "struct_tau": args.rean_struct_tau, "seed": args.rean_seed or R.REAN_SEED,
                    "n_replaced": int(len(man_r))}
        print(f"\nReanalysis splice: {len(man_r)} members from {args.splice_year} "
              f"({man_r.head_leaf.value_counts().to_dict()})")
        for nm, ens in rean.items():
            print(f"  {nm}: 2016-2025 mean {np.nanmean(ens[YEARS >= 2016]):+.3f} (rel. 1981-2010), "
                  f"EDA SD 2024 {np.nanstd(ens[YEARS == 2024], ddof=1):.4f}")
    manifest.to_csv(out_dir / "land_ensemble_manifest.csv", index=False)

    # ---------- save outputs ----------
    print("Saving land_ensemble.csv ...")
    cols = ["year"] + [f"m{i:05d}" for i in range(1, N_FINAL + 1)]
    out = pd.DataFrame(np.column_stack([YEARS, final]), columns=cols)
    out["year"] = out["year"].astype(int)
    out.to_csv(out_dir / "land_ensemble.csv", index=False, float_format="%.4f")

    summary = pd.DataFrame({
        "year": YEARS,
        "mean": np.nanmean(final, axis=1),
        "p025": np.nanpercentile(final, 2.5, axis=1),
        "p050": np.nanpercentile(final, 50, axis=1),
        "p500": np.nanpercentile(final, 50, axis=1),
        "p975": np.nanpercentile(final, 97.5, axis=1),
        "sd":   np.nanstd(final, axis=1, ddof=1),
    })
    summary.to_csv(out_dir / "land_ensemble_summary.csv", index=False, float_format="%.4f")
    print(f"Saved land_ensemble.csv ({out.shape}) and land_ensemble_summary.csv")

    # leaf-contribution summary
    leaf_counts = pd.DataFrame({"leaf": leaf_names})
    leaf_counts["share_2024"] = [
        (member_leaf_log[int(np.where(YEARS == 2024)[0][0])] == leaf_to_idx[n]).mean()
        for n in leaf_names
    ]
    leaf_counts["share_1900"] = [
        (member_leaf_log[int(np.where(YEARS == 1900)[0][0])] == leaf_to_idx[n]).mean()
        for n in leaf_names
    ]
    if rean_cfg:
        iy24 = int(np.where(YEARS == 2024)[0][0])
        extra = pd.DataFrame({"leaf": rean_cfg["leaves"],
                              "share_2024": [(member_leaf_log[iy24] == -1 - i).mean() for i in range(len(rean_cfg["leaves"]))],
                              "share_1900": 0.0})
        leaf_counts = pd.concat([leaf_counts, extra], ignore_index=True)
    print("\nLeaf shares of the final 10,000 in 1900 and 2024 (sanity vs. tree weights):")
    print(leaf_counts.to_string(index=False))
    (out_dir / "build_config.json").write_text(json.dumps({
        "method": args.method, "cru_version": CRU_VERSION if not legacy else args.cru_version,
        "coverage": None if legacy else args.coverage, "coverage_p_phi": cov_pp,
        "recentre": args.recentre, "donor": None if legacy else args.donor,
        "donor_smooth_years": smooth, "reanalysis": rean_cfg,
    }, indent=1))


if __name__ == "__main__":
    main()

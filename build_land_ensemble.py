"""
Build a 10,000-member ensemble of annual global-mean land surface air
temperature (LSAT) anomaly, 1850-2025, by combining five LSAT products
through a structural-method family tree with donor-imputed uncertainty
where native ensembles are not available.

See `land_ensemble_methodology.md` for full description.

Inputs (in `Land Data/`):
    Berkeley Earth Highres/Land_TAVG_ensemble.txt          (10-member monthly)
    CRUTEM5/CRUTEM.5.1.0.0.component_series.global.monthly.csv
    GloSATLAT/GloSATLAT-1-0-0-0_component-series_global_monthly.nc
    DCLSAT/DCLSAT_monthly_ensemble.csv                     (200-member monthly,
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

from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

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
RNG = np.random.default_rng(20260512)


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


def load_crutem5_components() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Load CRUTEM5 monthly anomaly + four uncertainty components on the 1850-2025
    monthly axis. Returns (years_idx, anomaly, u_unc, u_corr, u_cov, bias_sigma)."""
    df = pd.read_csv(DATA / "CRUTEM5" / "CRUTEM.5.1.0.0.component_series.global.monthly.csv")
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
    """Returns (N_YEARS, 200) annual ensemble. Native baseline 1982-2014."""
    path = DATA / "DCLSAT" / "DCLSAT_monthly_ensemble.csv"
    df = pd.read_csv(path)
    member_cols = [c for c in df.columns if c.startswith("m") and c != "month"]
    assert len(member_cols) == 200, f"expected 200 DCLSAT members, got {len(member_cols)}"
    mask = (df["year"] >= START_YEAR) & (df["year"] <= END_YEAR)
    df = df[mask]
    monthly = df[member_cols].values  # (months, 200)
    return reduce_monthly_to_annual(monthly, df["year"].values)


def load_noaa_land(donor_ensemble: np.ndarray, target_sigma: np.ndarray, rng: np.random.Generator) -> np.ndarray:
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
    rescale = np.where(donor_sigma > 0, target_sigma / donor_sigma, 0.0)[:, None]
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
# Leaves and their P-of-leaf for years where all 5 are available.
TREE_FULL = {
    "berkeley_earth": 0.25,
    "noaa_land":      0.25,
    "dclsat":         0.25,
    "crutem5":        0.125,
    "glosatlat":      0.125,
}
# For 2022-2025 (GloSATLAT absent), redistribute its 0.125 to CRUTEM5 (the
# other CRU-lineage leaf) so CRU-lineage family total stays 0.25.
TREE_NO_GLOSAT = {
    "berkeley_earth": 0.25,
    "noaa_land":      0.25,
    "dclsat":         0.25,
    "crutem5":        0.25,
    "glosatlat":      0.0,
}


def leaf_probs_for_year(year: int) -> dict[str, float]:
    return TREE_NO_GLOSAT if year > 2021 else TREE_FULL


# ---------- main ----------
def main() -> None:
    print("Loading Berkeley Earth (native 10 members)...")
    be = load_berkeley_earth()
    print(f"  shape: {be.shape}  finite years: {np.isfinite(be).all(axis=1).sum()}")

    print("Loading + synthesizing CRUTEM5 (200 synthetic members)...")
    cru_in = load_crutem5_components()
    cru = synthesize_component_ensemble(*cru_in, N_SYNTH, RNG)
    print(f"  shape: {cru.shape}  finite years: {np.isfinite(cru).all(axis=1).sum()}")

    print("Loading + synthesizing GloSATLAT (200 synthetic members)...")
    glo_in = load_glosat_components()
    glo = synthesize_component_ensemble(*glo_in, N_SYNTH, RNG)
    print(f"  shape: {glo.shape}  finite years: {np.isfinite(glo).all(axis=1).sum()}")

    print("Loading DCLSAT (200 native members)...")
    dcl = load_dclsat_ensemble()
    print(f"  shape: {dcl.shape}  finite years: {np.isfinite(dcl).all(axis=1).sum()}")

    print("Building NOAA Land donor ensemble (200 members; DCLSAT donor rescaled to GloSAT σ)...")
    # baseline DCLSAT first so the donor demeaning is sensible
    dcl_baselined_for_donor = rebaseline_to_1981_2010(dcl)
    cru_baselined_for_sigma = rebaseline_to_1981_2010(cru)
    glo_baselined_for_sigma = rebaseline_to_1981_2010(glo)
    tgt_sigma = target_sigma_from_glosat(glo_baselined_for_sigma, cru_baselined_for_sigma)
    noaa = load_noaa_land(dcl_baselined_for_donor, tgt_sigma, RNG)
    print(f"  shape: {noaa.shape}  finite years: {np.isfinite(noaa).all(axis=1).sum()}")

    print("Re-baselining all per-dataset ensembles to 1981-2010...")
    be_b = rebaseline_to_1981_2010(be)
    cru_b = rebaseline_to_1981_2010(cru)
    glo_b = rebaseline_to_1981_2010(glo)
    dcl_b = rebaseline_to_1981_2010(dcl)
    noaa_b = rebaseline_to_1981_2010(noaa)

    leaves: dict[str, np.ndarray] = {
        "berkeley_earth": be_b,
        "crutem5": cru_b,
        "glosatlat": glo_b,
        "dclsat": dcl_b,
        "noaa_land": noaa_b,
    }

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
    diag.to_csv(ROOT / "land_ensemble_perdataset.csv", index=False, float_format="%.5f")
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
                mi = fallback_member_idx[m]
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

    # ---------- save outputs ----------
    print("Saving land_ensemble.csv ...")
    cols = ["year"] + [f"m{i:05d}" for i in range(1, N_FINAL + 1)]
    out = pd.DataFrame(np.column_stack([YEARS, final]), columns=cols)
    out["year"] = out["year"].astype(int)
    out.to_csv(ROOT / "land_ensemble.csv", index=False, float_format="%.4f")

    summary = pd.DataFrame({
        "year": YEARS,
        "mean": np.nanmean(final, axis=1),
        "p025": np.nanpercentile(final, 2.5, axis=1),
        "p050": np.nanpercentile(final, 50, axis=1),
        "p500": np.nanpercentile(final, 50, axis=1),
        "p975": np.nanpercentile(final, 97.5, axis=1),
        "sd":   np.nanstd(final, axis=1, ddof=1),
    })
    summary.to_csv(ROOT / "land_ensemble_summary.csv", index=False, float_format="%.4f")
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
    print("\nLeaf shares of the final 10,000 in 1900 and 2024 (sanity vs. tree weights):")
    print(leaf_counts.to_string(index=False))


if __name__ == "__main__":
    main()

"""
Build a 10,000-member ensemble of annual global-mean sea-surface
temperature (SST) anomaly, 1850-2025, by combining four SST product
families through a structural-method family tree. As of v3 (May 2026)
the COBE family branch hosts two sibling variants of COBE-SST3 (Ishii
et al. 2025): a native 300-member perturbation ensemble covering
1870-2024, and a HadSST4-donor-imputed variant wrapping SST3's
best-estimate across its full 1850-2024 best-estimate window. Each
sibling carries P=1/8 of the tree inside 1870-2024; outside that
window the sibling-redirect routes draws to the donor variant
(1850-1869); beyond 2024 both COBE variants are NaN, the COBE family
contributes 0 to that year, and the other three families renormalise.

See `ocean_ensemble_methodology.md` for full description.

Inputs (in `Ocean Data/`):
    HadSST4/HadSST.4.2.0.0_global_ensemble_monthly.csv      (200 native bias members,
                                                              from prepare_hadsst_global_ensemble.py)
    HadSST4/HadSST.4.2.0.0_monthly_GLOBE.csv                (decomposed σ for noise terms)
    ERSSTv6/ERSSTv6_monthly_ensemble.csv                    (1000 native members 1850-2024,
                                                              from pull_ersstv6_members.py)
    ERSSTv6/aravg.mon.ocean.90S.90N.v6.1.0.202604.asc       (best-estimate central; used as
                                                              the post-native-end anchor)
    COBE-SST3/COBE-SST3_global_monthly.csv                  (deterministic best estimate
                                                              1850-2024, from
                                                              prepare_cobe_sst3_global_mean.py)
    COBE-SST3/COBE-SST3_monthly_ensemble.csv                (300 native perturbation members
                                                              1870-2024, from
                                                              pull_cobe_sst3_members.py)
    DCENT SST/DCENT_I_SST_monthly_ensemble.csv              (200 native members,
                                                              DCENT-I sea-fraction-weighted SST)

Outputs (in this folder):
    ocean_ensemble.csv                  176 × 10,001
    ocean_ensemble_summary.csv
    ocean_ensemble_perdataset.csv
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).parent
DATA = ROOT / "Ocean Data"
TMP = ROOT / "temp_files"
TMP.mkdir(exist_ok=True)

# ---------- configuration ----------
START_YEAR = 1850
END_YEAR = 2025
YEARS = np.arange(START_YEAR, END_YEAR + 1)
N_YEARS = len(YEARS)
BASELINE = (1981, 2010)
N_FINAL = 10_000
N_SYNTH = 200
RHO_MONTHLY = 0.3
RNG = np.random.default_rng(20260512)

# ERSSTv6's native 1000-member ensemble (Huang et al., NCEI pre-release)
# covers 1850-01 through 2024-12. For years beyond that we apply Option B
# in `extend_with_frozen_offset`: each member's deviation from the last
# native year's median is propagated forward (scaled by HadSST4's σ delta
# to capture any real coverage-driven σ change). This is only justified
# for a single year of extrapolation; the gap is hard-failed at >= 2 to
# force a revisit of the fallback design if the native ensemble falls
# multiple years behind.
ERSSTV6_NATIVE_END = 2024
FROZEN_OFFSET_HARD_FAIL_GAP = 2


# ---------- helpers (mirror build_land_ensemble.py) ----------
def reduce_monthly_to_annual(months_array: np.ndarray, year_idx: np.ndarray) -> np.ndarray:
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
    base_mask = (YEARS >= BASELINE[0]) & (YEARS <= BASELINE[1])
    means = np.nanmean(ensemble[base_mask], axis=0, keepdims=True)
    return ensemble - means


def ar1_noise(sigma_t: np.ndarray, n_members: int, rho: float, rng: np.random.Generator) -> np.ndarray:
    n = sigma_t.shape[0]
    out = np.zeros((n, n_members))
    out[0] = rng.standard_normal(n_members) * sigma_t[0]
    rho2 = np.sqrt(max(0.0, 1.0 - rho ** 2))
    for t in range(1, n):
        scale = sigma_t[t] / sigma_t[t - 1] if sigma_t[t - 1] > 0 else 0.0
        innov = rng.standard_normal(n_members) * sigma_t[t]
        out[t] = rho * scale * out[t - 1] + rho2 * innov
    return out


# ---------- per-dataset loaders ----------
def load_hadsst4_ensemble() -> np.ndarray:
    """HadSST4: 200 native bias members + AR(1) uncorrelated + scalar correlated +
    scalar coverage noise added from CSV σ. Returns (N_YEARS, 200) annual."""
    # native bias ensemble
    nat = pd.read_csv(DATA / "HadSST4" / "HadSST.4.2.0.0_global_ensemble_monthly.csv")
    mem_cols = [c for c in nat.columns if c.startswith("m") and c != "month"]
    assert len(mem_cols) == 200, f"expected 200 HadSST4 members, got {len(mem_cols)}"
    mask = (nat["year"] >= START_YEAR) & (nat["year"] <= END_YEAR)
    nat = nat[mask].reset_index(drop=True)
    bias_monthly = nat[mem_cols].values    # (months, 200)

    # σ components from CSV summary
    csv = pd.read_csv(DATA / "HadSST4" / "HadSST.4.2.0.0_monthly_GLOBE.csv")
    csv["yearmo"] = csv["year"] * 100 + csv["month"]
    nat["yearmo"] = nat["year"] * 100 + nat["month"]
    sig = csv.set_index("yearmo").loc[nat["yearmo"]].reset_index(drop=True)
    u_unc = sig["uncorrelated_uncertainty"].values
    u_corr = sig["correlated_uncertainty"].values
    u_cov = sig["coverage_uncertainty"].values

    # noise additions
    e_unc = ar1_noise(np.where(np.isfinite(u_unc), u_unc, 0.0), 200, RHO_MONTHLY, RNG)
    alpha = RNG.standard_normal(200)
    beta = RNG.standard_normal(200)
    e_corr = np.where(np.isfinite(u_corr), u_corr, 0.0)[:, None] * alpha[None, :]
    e_cov = np.where(np.isfinite(u_cov), u_cov, 0.0)[:, None] * beta[None, :]
    monthly_members = bias_monthly + e_unc + e_corr + e_cov
    bad = ~np.isfinite(bias_monthly[:, 0])
    monthly_members[bad, :] = np.nan
    return reduce_monthly_to_annual(monthly_members, nat["year"].values)


def load_dcent_sst_ensemble() -> np.ndarray:
    """Returns (N_YEARS, 200) annual ensemble. Native baseline 1982-2014.
    Backed by DCENT-I v1.1.0.0 (Chan et al. 2026, GDJ) — the spatially
    complete kriging-infilled extension of DCENT. Each member's native
    `sst` field has been area-weighted to a global mean by
    `pull_dcent_i_members.py` (see methodology §2.2 + sea-ice note in Step 1).
    """
    path = DATA / "DCENT SST" / "DCENT_I_SST_monthly_ensemble.csv"
    df = pd.read_csv(path)
    mem_cols = [c for c in df.columns if c.startswith("m") and c != "month"]
    assert len(mem_cols) == 200, f"expected 200 DCENT SST members, got {len(mem_cols)}"
    mask = (df["year"] >= START_YEAR) & (df["year"] <= END_YEAR)
    df = df[mask]
    monthly = df[mem_cols].values
    return reduce_monthly_to_annual(monthly, df["year"].values)


def load_ersstv6_native_ensemble() -> np.ndarray:
    """ERSSTv6: 1000 native members (Huang et al., NCEI pre-release), 1850-2024.
    Returns (N_YEARS, 1000) annual; years beyond ERSSTV6_NATIVE_END are NaN."""
    path = DATA / "ERSSTv6" / "ERSSTv6_monthly_ensemble.csv"
    df = pd.read_csv(path)
    mem_cols = [c for c in df.columns if c.startswith("m") and c != "month"]
    assert len(mem_cols) == 1000, f"expected 1000 ERSSTv6 members, got {len(mem_cols)}"
    mask = (df["year"] >= START_YEAR) & (df["year"] <= END_YEAR)
    df = df[mask]
    monthly = df[mem_cols].values
    return reduce_monthly_to_annual(monthly, df["year"].values)


def load_ersst_central_aravg() -> np.ndarray:
    """ERSSTv6 best-estimate annual anomaly from NOAAGlobalTemp v6.1 aravg.
    Used as the central anchor for the Option B post-native-end fallback."""
    rows = []
    with open(DATA / "ERSSTv6" / "aravg.mon.ocean.90S.90N.v6.1.0.202604.asc") as f:
        for line in f:
            parts = line.split()
            if not parts:
                continue
            yr = int(parts[0]); mo = int(parts[1]); anom = float(parts[2])
            rows.append((yr, mo, anom))
    df = pd.DataFrame(rows, columns=["year", "month", "anom"])
    df = df[(df["year"] >= START_YEAR) & (df["year"] <= END_YEAR)]
    return reduce_monthly_to_annual(df["anom"].values[:, None], df["year"].values)[:, 0]


def extend_with_frozen_offset(
    native: np.ndarray,
    central_baselined: np.ndarray,
    sigma_donor: np.ndarray,
    label: str,
) -> np.ndarray:
    """Option B fallback for any year where `native` is all-NaN but the build
    end is reached. For year y past the last native year y*:

        σ_scale = σ_donor[y] / σ_donor[y*]
        ensemble_median[y*] = nanmedian(native[y*, :])
        member_m[y]  = central_baselined[y]
                        + (native[y*, m] - ensemble_median[y*]) * σ_scale

    Preserves per-member trajectory continuity at the boundary year while
    letting HadSST4's late-year σ delta scale the spread for any real
    coverage-driven uncertainty growth. Hard-fails if the gap to END_YEAR
    is >= FROZEN_OFFSET_HARD_FAIL_GAP — frozen-offset is only defensible
    for a single year.
    """
    n_yrs, n_mem = native.shape
    finite_year_mask = np.isfinite(native).all(axis=1)
    last_native_iy = int(np.where(finite_year_mask)[0].max())
    last_native_year = int(YEARS[last_native_iy])
    gap = (n_yrs - 1) - last_native_iy
    if gap >= FROZEN_OFFSET_HARD_FAIL_GAP:
        raise RuntimeError(
            f"{label}: native ensemble ends at {last_native_year}; build end "
            f"is {END_YEAR} (gap = {gap} years). Frozen-offset (Option B) is "
            f"justified for one year only — REVISIT THE FALLBACK DESIGN before "
            f"shipping. Either obtain an updated native ensemble, or replace "
            f"`extend_with_frozen_offset` with a model that handles multi-year "
            f"extrapolation (e.g. AR(1) decay of the per-member offsets, or a "
            f"donor-scatter approach reintroduced for years > last_native + 1)."
        )
    if gap == 0:
        return native
    print(
        f"  {label}: native ensemble ends {last_native_year}; "
        f"frozen-offsetting {gap} year(s) to {END_YEAR} using Option B."
    )
    out = native.copy()
    median_last = np.nanmedian(native[last_native_iy])
    sigma_last = sigma_donor[last_native_iy]
    if sigma_last <= 0 or not np.isfinite(sigma_last):
        raise RuntimeError(
            f"{label}: σ_donor at last native year {last_native_year} is "
            f"non-positive ({sigma_last}); cannot scale frozen offset."
        )
    for iy in range(last_native_iy + 1, n_yrs):
        sigma_scale = sigma_donor[iy] / sigma_last if np.isfinite(sigma_donor[iy]) else 1.0
        offsets = native[last_native_iy] - median_last
        out[iy] = central_baselined[iy] + offsets * sigma_scale
        print(
            f"    {int(YEARS[iy])}: central={central_baselined[iy]:+.3f}  "
            f"σ_scale={sigma_scale:.3f}  resulting σ={np.nanstd(out[iy], ddof=1):.4f} °C"
        )
    return out


def load_cobe_sst3_central_monthly() -> tuple[np.ndarray, np.ndarray]:
    """COBE-SST3 deterministic best-estimate monthly anomaly (for diagnostics)."""
    df = pd.read_csv(DATA / "COBE-SST3" / "COBE-SST3_global_monthly.csv")
    df = df[(df["year"] >= START_YEAR) & (df["year"] <= END_YEAR)]
    return df["year"].values, df["anomaly"].values


def load_cobe_sst3_ensemble() -> np.ndarray:
    """COBE-SST3: 300 native perturbation members, monthly, 1870-2024.
    Returns (N_YEARS, 300) annual; years outside the native window are NaN."""
    path = DATA / "COBE-SST3" / "COBE-SST3_monthly_ensemble.csv"
    df = pd.read_csv(path)
    mem_cols = [c for c in df.columns if c.startswith("m") and c != "month"]
    assert len(mem_cols) == 300, f"expected 300 COBE-SST3 members, got {len(mem_cols)}"
    # The file natively spans 1870-2024; reindex onto the build year-month grid
    # so out-of-window months are NaN.
    full_idx = pd.MultiIndex.from_product(
        [np.arange(START_YEAR, END_YEAR + 1), np.arange(1, 13)],
        names=["year", "month"],
    )
    df = df.set_index(["year", "month"]).reindex(full_idx)
    monthly = df[mem_cols].values  # (N_YEARS*12, 300)
    year_idx = full_idx.get_level_values("year").values
    return reduce_monthly_to_annual(monthly, year_idx)


def sparse_era_inflation(years: np.ndarray) -> np.ndarray:
    """Inflation factor: 1.3× pre-1920, linear decay to 1.0× by 1960, 1.0× after."""
    f = np.where(years < 1920, 1.3,
        np.where(years < 1960, 1.3 - 0.3 * (years - 1920) / 40, 1.0))
    return f


def build_donor_ensemble(
    central_annual: np.ndarray,
    donor_ensemble_baselined: np.ndarray,
    target_sigma_base: np.ndarray,
    inflate_factor: np.ndarray,
) -> np.ndarray:
    """Wrap rescaled donor uncertainty around `central_annual`.
       central_annual:    (N_YEARS,)
       donor_ensemble:    (N_YEARS, n_donor)  -- already baselined
       target_sigma_base: (N_YEARS,)
       inflate_factor:    (N_YEARS,)
       Returns (N_YEARS, n_donor) annual ensemble.
    """
    donor_demeaned = donor_ensemble_baselined - np.nanmean(donor_ensemble_baselined, axis=1, keepdims=True)
    donor_sigma = np.nanstd(donor_ensemble_baselined, axis=1, ddof=1)
    target = target_sigma_base * inflate_factor
    rescale = np.where(donor_sigma > 0, target / donor_sigma, 0.0)[:, None]
    donor_scaled = donor_demeaned * rescale
    return central_annual[:, None] + donor_scaled


# ---------- family tree ----------
# v3 (May 2026) family tree: 5 leaves under a 4-branch root. Three branches
# carry one native-ensemble leaf each at P=1/4; the COBE branch hosts two
# siblings (SST3 native + SST2 donor-imputed) at P=1/8 each. In years
# outside SST3's 1870-2024 perturbation window the NaN-aware sampler
# renormalises and COBE-SST2 picks up the full 1/4 COBE-branch share.
#
# Effective independent uncertainty templates:
#   - 1850-1869: 3.5 (HadSST4, ERSSTv6, DCENT; COBE-SST2 shares HadSST4)
#   - 1870-2024: 4.0 (above + SST3's native perturbation template)
#   - 2025+:     3.5 (back to v2 effective count; SST3 leaf NaN)
TREE = {
    "hadsst4":         0.25,
    "dcent_sst":       0.25,
    "ersstv6":         0.25,
    "cobe_sst3":       0.125,   # native 300-member perturbation ensemble (1870-2024)
    "cobe_sst3_donor": 0.125,   # HadSST4-donor wrap around SST3's best estimate (1850-2024)
}

# Sibling-redirect rules for the NaN-aware sampler. When a leaf is NaN
# in a given year, the *generic* renormalisation would proportionally
# spread its weight across all finite leaves — which for the COBE branch
# would shrink the COBE family share in years where one sibling is NaN.
# Instead, a primary draw landing on a leaf listed here is redirected
# to its named sibling first; only if that sibling is *also* NaN
# (e.g. 2025, where both COBE variants lack data) do we fall back to
# the generic finite-leaves renormalisation. This preserves the COBE
# family at its full 1/4 of the tree wherever at least one sibling has
# data, and lets the COBE family naturally drop to 0 in years neither
# does.
LEAF_NAN_SIBLING = {
    "cobe_sst3":       "cobe_sst3_donor",  # 1850-1869: native is NaN; donor carries the branch
    "cobe_sst3_donor": "cobe_sst3",        # symmetric (rare/never used in current windows)
}


# ---------- main ----------
def main() -> None:
    print("Loading HadSST4 (200 native bias members + AR(1)/scalar noise)...")
    had = load_hadsst4_ensemble()
    print(f"  shape: {had.shape}  finite years: {np.isfinite(had).all(axis=1).sum()}")

    print("Loading DCENT SST (200 native members)...")
    dce = load_dcent_sst_ensemble()
    print(f"  shape: {dce.shape}  finite years: {np.isfinite(dce).all(axis=1).sum()}")

    print("Loading ERSSTv6 native ensemble (1000 members, 1850-2024) ...")
    ersst_native = load_ersstv6_native_ensemble()
    n_native_yrs = int(np.isfinite(ersst_native).all(axis=1).sum())
    print(f"  shape: {ersst_native.shape}  finite years: {n_native_yrs}")
    ersst_central_annual = load_ersst_central_aravg()
    print(f"  aravg central anchor finite years: {np.isfinite(ersst_central_annual).sum()}")

    print("Loading COBE-SST3 best-estimate monthly anomaly (deterministic) ...")
    yr_c3, anom_c3 = load_cobe_sst3_central_monthly()
    cobe_sst3_central_annual = reduce_monthly_to_annual(anom_c3[:, None], yr_c3)[:, 0]
    print(f"  finite years: {np.isfinite(cobe_sst3_central_annual).sum()}")

    print("Loading COBE-SST3 native perturbation ensemble (300 members, 1870-2024) ...")
    cobe_sst3 = load_cobe_sst3_ensemble()
    n3_native = int(np.isfinite(cobe_sst3).all(axis=1).sum())
    print(f"  shape: {cobe_sst3.shape}  finite years: {n3_native}")

    # Re-baseline native ensembles to 1981-2010 (member-wise)
    print("Re-baselining HadSST4, DCENT SST, ERSSTv6, COBE-SST3 native to 1981-2010 ...")
    had_b      = rebaseline_to_1981_2010(had)
    dce_b      = rebaseline_to_1981_2010(dce)
    ersst_b    = rebaseline_to_1981_2010(ersst_native)
    cobe_sst3_ens_b = rebaseline_to_1981_2010(cobe_sst3)
    # Re-baseline the SST3 best-estimate (the anchor for the donor variant and
    # the diagnostic column) by subtracting its own 1981-2010 mean.
    def _rebase_central(a: np.ndarray) -> np.ndarray:
        mask = (YEARS >= BASELINE[0]) & (YEARS <= BASELINE[1])
        return a - np.nanmean(a[mask])
    ersst_central_b   = _rebase_central(ersst_central_annual)
    cobe_sst3_cent_b  = _rebase_central(cobe_sst3_central_annual)

    # Apply Option B frozen-offset fallback for any year past ERSSTv6's native end
    print("Extending ERSSTv6 past native end with Option B frozen-offset ...")
    had_sigma = np.nanstd(had_b, axis=1, ddof=1)
    ersst_ens = extend_with_frozen_offset(
        native=ersst_b,
        central_baselined=ersst_central_b,
        sigma_donor=had_sigma,
        label="ERSSTv6",
    )

    # Build the SST3-donor sibling: HadSST4 spread wrapped around the SST3
    # best-estimate central, sparse-era inflated. The result is NaN in any
    # year where the SST3 central is NaN (2025+), so this leaf naturally
    # falls out of the sampler past the SST3 best-estimate end.
    print("Building COBE-SST3-donor sibling (HadSST4 donor, sparse-era inflated) ...")
    infl = sparse_era_inflation(YEARS)
    cobe_sst3_donor_ens = build_donor_ensemble(cobe_sst3_cent_b, had_b, had_sigma, infl)

    leaves = {
        "hadsst4":         had_b,
        "dcent_sst":       dce_b,
        "ersstv6":         ersst_ens,
        "cobe_sst3":       cobe_sst3_ens_b,
        "cobe_sst3_donor": cobe_sst3_donor_ens,
    }

    # Sanity print at several diagnostic years (covering pre-SST3-window,
    # mid-window, and post-window for the COBE branch)
    print("\nSanity: ensemble mean (best estimate) by year (re-baselined to 1981-2010):")
    print(f"  {'leaf':<10s}  {'1880':>14s}  {'1900':>14s}  {'2020':>14s}  {'2024':>14s}  {'2025':>14s}")
    for name, ens in leaves.items():
        row = []
        for y in (1880, 1900, 2020, 2024, 2025):
            iy = int(np.where(YEARS == y)[0][0])
            if np.isfinite(ens[iy]).any():
                mean = np.nanmean(ens[iy])
                sd = np.nanstd(ens[iy], ddof=1)
                row.append(f"{mean:+.3f}±{sd:.3f}")
            else:
                row.append("    NaN     ")
        print(f"  {name:<10s}  " + "  ".join(f"{x:>14s}" for x in row))

    # diagnostics CSV (per-dataset best estimate + 1σ). For COBE branch,
    # also emit the SST3 deterministic best estimate as a separate column
    # so we can visualise how the SST3 best-estimate compares to SST2 and
    # to the SST3-ensemble mean across the SST3 record.
    diag = pd.DataFrame({"year": YEARS})
    for name, ens in leaves.items():
        diag[f"{name}_mean"] = np.nanmean(ens, axis=1)
        diag[f"{name}_sd"]   = np.nanstd(ens, axis=1, ddof=1)
    # The SST3 best-estimate central is logged separately as a diagnostic;
    # the SST3 native and SST3 donor leaves both inherit it (the former
    # via member trajectories that center on it; the latter as the explicit
    # anchor of the HadSST4 donor wrap).
    diag["cobe_sst3_central"] = cobe_sst3_cent_b
    diag.to_csv(ROOT / "ocean_ensemble_perdataset.csv", index=False, float_format="%.5f")
    print(f"\nWrote ocean_ensemble_perdataset.csv")

    # ---------- family-tree sampling to 10000 ----------
    print(f"\nSampling {N_FINAL} ensemble members ...")
    leaf_names = list(leaves.keys())
    probs = np.array([TREE[n] for n in leaf_names])
    leaf_to_idx = {n: i for i, n in enumerate(leaf_names)}
    primary_leaf = RNG.choice(len(leaf_names), size=N_FINAL, p=probs)
    leaf_member_counts = {n: leaves[n].shape[1] for n in leaf_names}
    primary_member_idx = np.array([
        RNG.integers(0, leaf_member_counts[leaf_names[primary_leaf[m]]]) for m in range(N_FINAL)
    ])

    # NaN-aware per-year sampler (same pattern as build_land_ensemble.py)
    # with explicit sibling redirect for LEAF_NAN_SIBLING entries (see TREE comment).
    finite_by_year = {
        iy: [i for i, n in enumerate(leaf_names) if np.isfinite(leaves[n][iy]).any()]
        for iy, _ in enumerate(YEARS)
    }
    final = np.full((N_YEARS, N_FINAL), np.nan)
    member_leaf_log = np.zeros((N_YEARS, N_FINAL), dtype=np.int8)
    for iy, yr in enumerate(YEARS):
        finite_leaves = finite_by_year[iy]
        renorm = np.array([TREE[leaf_names[i]] for i in finite_leaves])
        renorm = renorm / renorm.sum() if renorm.sum() > 0 else renorm
        for m in range(N_FINAL):
            pl = primary_leaf[m]; mi = primary_member_idx[m]
            if not np.isfinite(leaves[leaf_names[pl]][iy, mi]):
                # Step 1: sibling-redirect if rule exists for the primary leaf
                sibling_name = LEAF_NAN_SIBLING.get(leaf_names[pl])
                redirected = False
                if sibling_name is not None:
                    sib_idx = leaf_to_idx[sibling_name]
                    if np.isfinite(leaves[sibling_name][iy]).any():
                        pl = sib_idx
                        mi = RNG.integers(0, leaf_member_counts[sibling_name])
                        redirected = True
                # Step 2: fall back to generic finite-leaves renormalisation
                if not redirected:
                    pl = RNG.choice(finite_leaves, p=renorm)
                    mi = RNG.integers(0, leaf_member_counts[leaf_names[pl]])
                # Step 3: redraw on rare per-member NaN within the chosen leaf
                tries = 0
                while not np.isfinite(leaves[leaf_names[pl]][iy, mi]) and tries < 20:
                    mi = RNG.integers(0, leaf_member_counts[leaf_names[pl]])
                    tries += 1
            final[iy, m] = leaves[leaf_names[pl]][iy, mi]
            member_leaf_log[iy, m] = pl

    print("Saving ocean_ensemble.csv ...")
    cols = ["year"] + [f"m{i:05d}" for i in range(1, N_FINAL + 1)]
    out = pd.DataFrame(np.column_stack([YEARS, final]), columns=cols)
    out["year"] = out["year"].astype(int)
    out.to_csv(ROOT / "ocean_ensemble.csv", index=False, float_format="%.4f")

    summary = pd.DataFrame({
        "year": YEARS,
        "mean": np.nanmean(final, axis=1),
        "p025": np.nanpercentile(final, 2.5, axis=1),
        "p500": np.nanpercentile(final, 50, axis=1),
        "p975": np.nanpercentile(final, 97.5, axis=1),
        "sd":   np.nanstd(final, axis=1, ddof=1),
    })
    summary.to_csv(ROOT / "ocean_ensemble_summary.csv", index=False, float_format="%.4f")
    print(f"Saved ocean_ensemble.csv ({out.shape}) and ocean_ensemble_summary.csv")

    leaf_counts = pd.DataFrame({"leaf": leaf_names})
    for y in (1860, 1900, 2020, 2024, 2025):
        iy = int(np.where(YEARS == y)[0][0])
        leaf_counts[f"share_{y}"] = [
            (member_leaf_log[iy] == leaf_to_idx[n]).mean() for n in leaf_names
        ]
    print("\nLeaf shares of the final 10,000 (COBE branch sums to ~25% in 1850–2024;\n"
          "drops to ~0% in 2025 where neither SST3 variant has data):")
    print(leaf_counts.to_string(index=False))

    # ---------- family-tree invariants (hard-fail if violated) ----------
    # COBE family = cobe_sst3 + cobe_sst3_donor. ±0.020 tolerance is generous
    # against 10k-Bernoulli 3σ at p=0.25 (~0.013).
    print("\nFamily-tree invariant check ...")
    SHARE_TOL = 0.020
    invariant_failures: list[str] = []

    def leaf_share(year: int, leaf: str) -> float:
        iy = int(np.where(YEARS == year)[0][0])
        return float((member_leaf_log[iy] == leaf_to_idx[leaf]).mean())

    for y in (1860, 1900, 2020, 2024, 2025):
        cobe = leaf_share(y, "cobe_sst3") + leaf_share(y, "cobe_sst3_donor")
        hsst = leaf_share(y, "hadsst4")
        erst = leaf_share(y, "ersstv6")
        dcnt = leaf_share(y, "dcent_sst")

        # The COBE family is expected to carry 1/4 in 1850–2024 (where at
        # least one COBE sibling is finite) and 0 in 2025 (where both are
        # NaN and the other three families renormalise to ~1/3 each).
        if y <= 2024:
            for fam, val in (("COBE", cobe), ("HadSST", hsst), ("ERSST", erst), ("DCENT", dcnt)):
                if abs(val - 0.25) > SHARE_TOL:
                    invariant_failures.append(
                        f"{y}: family {fam} share = {val:.4f} (expected 0.25 ± {SHARE_TOL})"
                    )
        else:
            # 2025: COBE family contributes 0; other three renormalise to ~1/3 each
            if cobe != 0.0:
                invariant_failures.append(
                    f"{y}: COBE family share = {cobe:.4f} (expected 0; both variants NaN)"
                )
            for fam, val in (("HadSST", hsst), ("ERSST", erst), ("DCENT", dcnt)):
                if abs(val - 1.0 / 3.0) > SHARE_TOL:
                    invariant_failures.append(
                        f"{y}: family {fam} share = {val:.4f} (expected 0.333 ± {SHARE_TOL}; "
                        f"COBE family absent, renorm over remaining three)"
                    )

        # Within-COBE branch: 1/8 + 1/8 inside 1870–2024, and donor carries
        # the full 1/4 in 1850–1869 (sibling-redirect from the NaN native).
        if 1870 <= y <= 2024:
            n3 = leaf_share(y, "cobe_sst3")
            nd = leaf_share(y, "cobe_sst3_donor")
            if abs(n3 - 0.125) > SHARE_TOL:
                invariant_failures.append(
                    f"{y}: cobe_sst3 native share = {n3:.4f} inside window (expected 0.125)")
            if abs(nd - 0.125) > SHARE_TOL:
                invariant_failures.append(
                    f"{y}: cobe_sst3_donor share = {nd:.4f} inside window (expected 0.125)")
        elif y <= 1869:
            n3 = leaf_share(y, "cobe_sst3")
            nd = leaf_share(y, "cobe_sst3_donor")
            if n3 != 0.0:
                invariant_failures.append(
                    f"{y}: cobe_sst3 native share = {n3:.4f} before 1870 (expected 0)")
            if abs(nd - 0.25) > SHARE_TOL:
                invariant_failures.append(
                    f"{y}: cobe_sst3_donor share = {nd:.4f} before 1870 "
                    f"(expected 0.25 via sibling-redirect from NaN native)")

    if invariant_failures:
        for f in invariant_failures:
            print(f"  FAIL: {f}")
        raise RuntimeError(
            "Family-tree invariant violated — sibling-redirect or TREE weights "
            "are mis-specified. See failures above."
        )
    print("  All family-share invariants satisfied at 1860/1900/2020/2024/2025.")

    # ---------- boundary σ-step diagnostics (soft) ----------
    # Logs the COBE-branch ensemble σ on either side of the COBE-window
    # boundaries (1869→1870 and 2024→2025) so wildly mis-scaled spreads
    # would be visible at run time.
    cobe_branch_member_mask = (
        (member_leaf_log == leaf_to_idx["cobe_sst3"]) |
        (member_leaf_log == leaf_to_idx["cobe_sst3_donor"])
    )
    for y_in, y_out in ((1870, 1869), (2024, 2025)):
        i_in = int(np.where(YEARS == y_in)[0][0])
        i_out = int(np.where(YEARS == y_out)[0][0])
        cobe_in = final[i_in][cobe_branch_member_mask[i_in]]
        cobe_out = final[i_out][cobe_branch_member_mask[i_out]]
        if cobe_in.size > 0 and cobe_out.size > 0:
            s_in = np.nanstd(cobe_in, ddof=1)
            s_out = np.nanstd(cobe_out, ddof=1)
            ratio = s_in / s_out if s_out > 0 else float("nan")
            flag = " <-- WARNING: >50% jump" if (s_out > 0 and abs(ratio - 1.0) > 0.5) else ""
            print(f"  COBE-branch σ {y_out}→{y_in}: {s_out:.4f} → {s_in:.4f} "
                  f"(ratio {ratio:.3f}){flag}")
        else:
            # 2024→2025: COBE branch has no members in 2025 (both NaN)
            print(f"  COBE-branch σ {y_out}→{y_in}: "
                  f"{'(absent)' if cobe_out.size == 0 else f'{np.nanstd(cobe_out, ddof=1):.4f}'} "
                  f"→ "
                  f"{'(absent)' if cobe_in.size == 0 else f'{np.nanstd(cobe_in, ddof=1):.4f}'}")

    # ---------- convention log ----------
    print("\nCOBE-branch sea-ice convention composition:")
    print("  1850-1869:  100% SST3-donor (SST3 central interpolated near-freezing,"
          " HadSST4 spread template)")
    print("  1870-2024:  50% SST3-native / 50% SST3-donor"
          " (both use SST3's interpolated-near-freezing convention)")
    print("  2025+:      COBE family contributes 0")


if __name__ == "__main__":
    main()

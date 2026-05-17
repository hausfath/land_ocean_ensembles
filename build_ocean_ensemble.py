"""
Build a 10,000-member ensemble of annual global-mean sea-surface
temperature (SST) anomaly, 1850-2025, by combining four SST products
through a structural-method family tree with donor-imputed uncertainty
for ERSSTv6 and COBE-SST2.

See `ocean_ensemble_methodology.md` for full description.

Inputs (in `Ocean Data/`):
    HadSST4/HadSST.4.2.0.0_global_ensemble_monthly.csv      (200 native bias members,
                                                              from prepare_hadsst_global_ensemble.py)
    HadSST4/HadSST.4.2.0.0_monthly_GLOBE.csv                (decomposed σ for noise terms)
    ERSSTv6/ERSSTv6_monthly_ensemble.csv                    (1000 native members 1850-2024,
                                                              from pull_ersstv6_members.py)
    ERSSTv6/aravg.mon.ocean.90S.90N.v6.1.0.202604.asc       (best-estimate central; used as
                                                              the post-native-end anchor)
    COBE-SST2/COBE-SST2_global_monthly.csv                  (deterministic, derived)
    DCENT SST/DCENT_SST_monthly_ensemble.csv                (200 native members)

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
    path = DATA / "DCENT SST" / "DCENT_SST_monthly_ensemble.csv"
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


def load_cobe_monthly() -> tuple[np.ndarray, np.ndarray]:
    """COBE-SST2 monthly anomaly from our prepared global-mean CSV."""
    df = pd.read_csv(DATA / "COBE-SST2" / "COBE-SST2_global_monthly.csv")
    df = df[(df["year"] >= START_YEAR) & (df["year"] <= END_YEAR)]
    return df["year"].values, df["anomaly"].values


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
# Equal-weight 4-leaf tree. ERSSTv6 carries its own native 1000-member
# ensemble (Huang et al., NCEI pre-release) as of v2; only COBE-SST2
# remains donor-imputed from HadSST4. Three of the four leaves now
# contribute independent uncertainty templates (HadSST4 bias-perturbation
# family, ERSSTv6 parameter family, DCENT joint-constraint family);
# effective independent leaves ≈ 3.5.
TREE = {
    "hadsst4":   0.25,
    "dcent_sst": 0.25,
    "ersstv6":   0.25,
    "cobe_sst2": 0.25,
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

    print("Loading COBE-SST2 monthly anomaly (deterministic) ...")
    yr_c, anom_c = load_cobe_monthly()
    cobe_annual = reduce_monthly_to_annual(anom_c[:, None], yr_c)[:, 0]
    print(f"  finite years: {np.isfinite(cobe_annual).sum()}")

    # Re-baseline native ensembles to 1981-2010 (member-wise)
    print("Re-baselining HadSST4, DCENT SST, ERSSTv6 native to 1981-2010 ...")
    had_b   = rebaseline_to_1981_2010(had)
    dce_b   = rebaseline_to_1981_2010(dce)
    ersst_b = rebaseline_to_1981_2010(ersst_native)
    # Re-baseline deterministic series (their best estimates) by subtracting their
    # 1981-2010 means
    def _rebase_central(a: np.ndarray) -> np.ndarray:
        mask = (YEARS >= BASELINE[0]) & (YEARS <= BASELINE[1])
        return a - np.nanmean(a[mask])
    ersst_central_b = _rebase_central(ersst_central_annual)
    cobe_b          = _rebase_central(cobe_annual)

    # Apply Option B frozen-offset fallback for any year past ERSSTv6's native end
    print("Extending ERSSTv6 past native end with Option B frozen-offset ...")
    had_sigma = np.nanstd(had_b, axis=1, ddof=1)
    ersst_ens = extend_with_frozen_offset(
        native=ersst_b,
        central_baselined=ersst_central_b,
        sigma_donor=had_sigma,
        label="ERSSTv6",
    )

    # Build donor uncertainty for COBE only (ERSSTv6 is now native)
    print("Building COBE-SST2 donor ensemble (HadSST4 donor, sparse-era inflated) ...")
    infl = sparse_era_inflation(YEARS)
    cobe_ens = build_donor_ensemble(cobe_b, had_b, had_sigma, infl)

    leaves = {
        "hadsst4":   had_b,
        "dcent_sst": dce_b,
        "ersstv6":   ersst_ens,
        "cobe_sst2": cobe_ens,
    }

    # Sanity print
    print("\nSanity: ensemble mean (best estimate) for 2024 (re-baselined to 1981-2010):")
    iy = int(np.where(YEARS == 2024)[0][0])
    for name, ens in leaves.items():
        print(f"  {name:10s}: {np.nanmean(ens[iy]):+.3f} °C  spread(1σ)={np.nanstd(ens[iy], ddof=1):.3f}")

    # diagnostics CSV
    diag = pd.DataFrame({"year": YEARS})
    for name, ens in leaves.items():
        diag[f"{name}_mean"] = np.nanmean(ens, axis=1)
        diag[f"{name}_sd"]   = np.nanstd(ens, axis=1, ddof=1)
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
                pl = RNG.choice(finite_leaves, p=renorm)
                mi = RNG.integers(0, leaf_member_counts[leaf_names[pl]])
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
    leaf_counts["share_2024"] = [
        (member_leaf_log[int(np.where(YEARS == 2024)[0][0])] == leaf_to_idx[n]).mean()
        for n in leaf_names
    ]
    leaf_counts["share_1900"] = [
        (member_leaf_log[int(np.where(YEARS == 1900)[0][0])] == leaf_to_idx[n]).mean()
        for n in leaf_names
    ]
    print("\nLeaf shares of the final 10,000 in 1900 and 2024:")
    print(leaf_counts.to_string(index=False))


if __name__ == "__main__":
    main()

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
    ERSSTv6/aravg.mon.ocean.90S.90N.v6.1.0.202604.asc       (deterministic only)
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


def load_ersst_monthly() -> tuple[np.ndarray, np.ndarray]:
    """ERSSTv6 monthly anomaly from NOAAGlobalTemp v6.1 aravg.mon.ocean ASCII."""
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
    return df["year"].values, df["anom"].values


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
# Equal-weight 4-leaf tree. See methodology MD Step 4 for rationale —
# in v1, ERSSTv6 and COBE-SST2 are donor-imputed from HadSST4 (no public
# ensembles available), but we use equal 1/4 weights in anticipation of
# acquiring a native ERSSTv6 (or NOAA) ensemble in a future revision.
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

    print("Loading ERSSTv6 monthly anomaly (deterministic) ...")
    yr_e, anom_e = load_ersst_monthly()
    ersst_annual = reduce_monthly_to_annual(anom_e[:, None], yr_e)[:, 0]
    print(f"  finite years: {np.isfinite(ersst_annual).sum()}")

    print("Loading COBE-SST2 monthly anomaly (deterministic) ...")
    yr_c, anom_c = load_cobe_monthly()
    cobe_annual = reduce_monthly_to_annual(anom_c[:, None], yr_c)[:, 0]
    print(f"  finite years: {np.isfinite(cobe_annual).sum()}")

    # Re-baseline native ensembles to 1981-2010 (member-wise)
    print("Re-baselining HadSST4 and DCENT SST to 1981-2010 ...")
    had_b = rebaseline_to_1981_2010(had)
    dce_b = rebaseline_to_1981_2010(dce)
    # Re-baseline deterministic series (their best estimates) by subtracting their
    # 1981-2010 means
    def _rebase_central(a: np.ndarray) -> np.ndarray:
        mask = (YEARS >= BASELINE[0]) & (YEARS <= BASELINE[1])
        return a - np.nanmean(a[mask])
    ersst_b = _rebase_central(ersst_annual)
    cobe_b  = _rebase_central(cobe_annual)

    # Build donor uncertainty for ERSST and COBE using HadSST4 as donor
    print("Building ERSSTv6 / COBE-SST2 donor ensembles (HadSST4 donor, sparse-era inflated)...")
    had_sigma = np.nanstd(had_b, axis=1, ddof=1)
    infl = sparse_era_inflation(YEARS)
    ersst_ens = build_donor_ensemble(ersst_b, had_b, had_sigma, infl)
    cobe_ens  = build_donor_ensemble(cobe_b,  had_b, had_sigma, infl)

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

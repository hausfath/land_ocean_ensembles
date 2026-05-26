# Land-only LSAT ensemble — methodology

This note describes how we build a 10,000-member ensemble of annual global
mean **land surface air temperature (LSAT)** anomalies, approximately
replicating the family-tree-with-imputation approach used for GMST by
Thorne et al. (2026, ESSD discussion paper 2025-825), but applied to
land-only datasets and with a much smaller catalogue.

The goal is to produce an LSAT analogue of the Thorne SST-grouped
ensemble — a component-domain ensemble of global LSAT with
structural-method uncertainty captured by the family tree. This
document concerns *only* the LSAT side; the sister ocean methodology
is described in `ocean_ensemble_methodology.md`.

## Input datasets (LSAT, qualifying products in May 2026)

| # | Dataset | Native ensemble? | Uncertainty form | Time span | Native baseline | File(s) |
|---|---------|------------------|------------------|-----------|-----------------|---------|
| 1 | Berkeley Earth High-Resolution Land | 10 members | ensemble | 1750-01 – 2026-04 | 1951-1980 | `Land Data/Berkeley Earth Highres/Land_TAVG_ensemble.txt` |
| 2 | CRUTEM5 v5.1.0.0 | No | uncorrelated σ + correlated σ + coverage σ + bias 2.5/97.5 CI | 1850-01 – 2026-XX | 1961-1990 | `Land Data/CRUTEM5/CRUTEM.5.1.0.0.component_series.global.monthly.csv` |
| 3 | GloSATLAT v1.0.0.0 | No | best estimate + 2.5/97.5 CI (and decomposed σ components) | 1781-01 – 2021-12 | 1961-1990 | `Land Data/GloSATLAT/GloSATLAT-1-0-0-0_{summary,component}-series_global_monthly.nc` |
| 4 | DCLSAT (DCENT v3.0 `lsat` field) | 200 members (gridded; we will reduce to global means) | ensemble of global means we compute ourselves | 1850-01 – 2025-12 | 1982-2014 | `Land Data/DCLSAT/DCENT_ensemble_mean_1850_2025.nc` (mean only currently; 200 members to be pulled) |
| 5 | NOAA Land v6.1.0 | No | deterministic best estimate only — v6.1 aravg text file has variance columns set to `-999`; merged-product gridded NetCDF has no separate land uncertainty field | 1850-01 – 2026-04 | 1971-2000 (aravg text) / 1991-2020 (gridded) | `Land Data/NOAA Land/aravg.mon.land.90S.90N.v6.1.0.202604.asc` |
| 6 | C-LSAT 2.1 (Sun Yat-sen University / CMA-homogenization lineage) | No | deterministic best estimate only | 1850-01 – 2025-12 | 1961-1990 | `Land Data/C-LSAT/China-LSAT2.1_tavg.nc` (raw) → `Land Data/C-LSAT/C-LSAT2.1_global_monthly.csv` (derived via `prepare_c_lsat_global_mean.py`) |

### Time-coverage note

GloSATLAT ends in 2021-12. The Thorne split-and-splice scheme has a
**tail** (pre-1996) and **head** (1996–present) sub-ensemble. GloSATLAT
qualifies for the tail. For the head it does not. Berkeley Earth, CRUTEM5,
DCLSAT, NOAA Land all cover both tail and head.

## Step 1 — Common time/baseline grid

1. Reduce every dataset to **annual means** (Jan-Dec) of the global-mean
   monthly anomaly. For monthly-uncertainty datasets without ensemble we
   propagate uncertainty to the annual mean explicitly (see Step 2).
2. Truncate to **1850–2025**. 2025 is the latest complete calendar year
   as of the working date 2026-05-12 for Berkeley Earth, CRUTEM5,
   DCLSAT and NOAA Land. GloSATLAT ends in 2021-12 and so does not
   contribute to 2022–2025; in those years its 0.125 weight is
   redistributed to CRUTEM5 (the only other CRU-lineage leaf), as
   defined in Step 4. The information loss is minimal because
   CRUTEM5 and GloSATLAT have r=0.999 annual correlation and
   0.026 °C RMS difference over their overlap 1858–2021, i.e. they
   are effectively the same product on the 1850-onward axis.
3. Re-baseline each dataset (and each ensemble member) to a **zero mean
   over 1981–2010**. This matches the Thorne choice and is the modern
   reference period used in AR7 WG1. The 1981–2010 mean is computed on
   the *member* (not just the central estimate) so cross-member spread is
   preserved on the baseline.

## Step 2 — Build a per-dataset ensemble

We target a per-dataset ensemble of native size (no synthetic
upsampling); the family tree (Step 4) handles the expansion to the
final 10,000 draws. Per-dataset details:

### 2.1 Berkeley Earth — keep 10 native members (no bootstrap)

The Berkeley Earth High-Resolution release ships 10 members. We keep
all 10 at native count and let the family-tree sampler (Step 5) draw
uniformly among them. Result: roughly 2,500 of the final 10,000
ensemble draws will be one of 10 distinct BE trajectories.

We **do not** bootstrap or parametric-augment to 200. Bootstrapping
would pretend to higher effective N without adding information, and
parametric augmentation would treat BE's structural perturbations as
Gaussian noise, which they are not. The 10-trajectory granularity will
appear as visible step-clustering in tail percentiles — this is
informative (it advertises BE's under-sampling) rather than a bug.

### 2.2 CRUTEM5 — synthesize from variance components

CRUTEM5 publishes monthly:
- Best-estimate anomaly *a(t)*
- Uncorrelated 1σ uncertainty *u_u(t)* (independent month-to-month)
- Correlated 1σ uncertainty *u_c(t)* (sampling/method correlated across time)
- Coverage 1σ uncertainty *u_v(t)* (correlated across time)
- Bias component 2.5/97.5 CI (derived from a set of bias realisations)

We synthesize a 200-member monthly ensemble as
*a(t) + e_u(t) + α · u_c(t) + β · u_v(t) + γ · b(t)*
where:
- *e_u(t)* — monthly noise term with AR(1) coefficient ρ = 0.3 month-to-month
  (rather than fully independent draws). Each per-member trajectory
  starts from N(0, u_u(0)) and evolves as
  *e_u(t) = ρ · e_u(t−1) · (u_u(t)/u_u(t−1)) + √(1−ρ²) · N(0, u_u(t))*.
  The ρ = 0.3 value is a reasonable compromise between
  "fully independent" (which under-represents sub-annual persistence and
  is wiped out by annual averaging) and "fully correlated" (which is
  already handled by the *u_c* term).
- *α, β* — per-member draws from N(0,1), constant across time, applied
  with the time-varying σ. This is the standard "correlated nuisance
  scalar" approximation Morice et al. (2021) use. It **overstates**
  low-frequency uncertainty (it treats the correlated/coverage σ as
  having an infinite decorrelation timescale, when the real timescale
  is years to decades). We retain it as the simplest defensible
  approximation; full covariance matrices would be cleaner if available.
- *γ · b(t)* — the published bias CI translated into a 1σ envelope
  ( (upper − lower)/(2·1.96) ); *γ ~ N(0,1)* per member, constant in
  time. This is a **fallback**: the bias CI is generated from
  underlying multi-realisation bias trajectories, not from a single σ.
  Where the bias realisations are exposed (HadCRUT5 publishes them;
  CRUTEM5 v5.1 does **not**) we should sample realisations directly.
  Since CRUTEM5 only ships the CI, the Gaussian scalar is what we use,
  flagged here as the weakest link in this dataset's synthesis.

Annual means are computed from each 12-month sequence per member.
This is a deliberate simplification of HadCRUT5's full ensemble
construction — we are not trying to replicate it bit-for-bit, only to
generate plausible LSAT realisations with the correct uncertainty
magnitude and approximate temporal correlation.

After synthesis we verify empirically that the resulting annual σ
matches CRUTEM5's published annual CI to within ~10%; if it
under-shoots we raise ρ from 0.3 toward 0.6 to extend the persistence
of the uncorrelated component.

### 2.3 GloSATLAT — Gaussian draws scaled to native CI, donor temporal structure

GloSATLAT's `summary-series` file provides `tas`, `tas_lower`, `tas_upper`
at monthly resolution. We treat `(upper − lower)/(2·1.96)` as the
monthly 1σ. The `component-series` decomposes that σ into uncorrelated,
correlated, coverage, and bias parts — we use these like CRUTEM5
(Section 2.2) with the same correlated/independent/bias-constant
structure. This guarantees that the resulting 200 ensemble members
recover GloSATLAT's published 95% CI.

GloSATLAT ends in 2021. We do **not** extend it; instead, the family
tree assigns zero weight to GloSATLAT for 2022–present (see Step 4).

### 2.4 DCLSAT — pull 200 native DCENT-I members, extract LSAT global means

We use **DCENT-I v1.1.0.0** (Chan et al. 2026, Geoscience Data
Journal; doi:10.7910/DVN/ROG38Q, dataset version updated 2026-03-30
through end of 2025), the spatially complete kriging-infilled
extension of DCENT. The original unfilled DCENT v3.0 derived data
are archived under `Land Data/DCLSAT/v3_archive/` for reproducibility.

The DCENT-I v1.1.0.0 Harvard Dataverse release exposes 200 gridded
NetCDFs (~42 MB each, ~8.4 GB total). Each member NetCDF contains
**three separate fields** — `sst` (sea_surface_temperature_anomaly),
`lsat` (land near-surface air temperature anomaly), and `ts` (merged
Surface Temperature Anomaly) — on a 5°×5° monthly grid from 1850-01
through 2025-12. The decomposition is performed by the DCENT-I authors
themselves via their ordinary-kriging infill (with anisotropic
heterogeneous kernels) and air-temperature blending over sea-ice cells.

We download each member NetCDF, compute the area-weighted (cos lat)
global mean of the `lsat` field directly (skipping NaN cells), then
delete the raw NetCDF. The same script writes the parallel SST CSV
from each member's `sst` field (ocean methodology §2.2). Peak
transient disk: ~42 MB.

**Sea-ice cells in the DCENT-I decomposition.** Both `sst` and `lsat`
are finite over sea-ice cells (overlap ~20% of the globe), reflecting
DCENT-I's design choice to assign kriged values to BOTH fields at ice
locations. The `sst` value at an ice cell is the kriged surface
temperature including the air-temperature blending; the `lsat` value
similarly reflects the kriged 2-m air temperature. We accept these
as the published DCENT-I quantities rather than re-mask.

**Why we switched from DCENT v3.0 to DCENT-I:**
- DCENT-I is the published, peer-reviewed infilled version.
- Resolves the "no extrapolation beyond observations" coverage gap
  in DCENT v3.0 — modern-era polar regions are now spatially complete.
- 200-member ensemble explicitly samples kriging uncertainty alongside
  the bias-adjustment uncertainty already in DCENT.
- Empirical impact on the headline: 2024 LSAT median shifts by
  −0.03 °C (slight cooling — kriging dampens member-to-member spread
  in well-observed cells), 2024 SST shifts by +0.02 °C (modest polar
  warming uplift); both changes within Monte-Carlo uncertainty.

**Important caveat #1 (carried over from DCENT v3.0):** each DCENT-I
member's LSAT-over-land field is paired with a specific SST field
under the underlying DCENT dynamical-consistency constraint. The
200-member spread of LSAT marginals *under-represents* what a
standalone-LSAT product would show by ~10–30%.

**Important caveat #2 (new for DCENT-I):** the kriging-infilled
modern-era spread is ~40% tighter than DCENT v3.0's was in the same
era (e.g., 2024 SST σ across 200 members: 0.0145 → 0.0088 °C). The
pre-1900 spread is *wider* (1850 LSAT σ: 0.065 → 0.11 °C) because
kriging correctly carries large uncertainty into sparse-data
extrapolation. Net: DCENT-I correctly redistributes σ contribution
across the record — less spread where data is dense, more where it's
sparse. The leaf still carries appropriate sparse-era uncertainty.

(Previous v1.0.0.0 of DCENT-I covered only 1850–2024 and required
weight-based decomposition of a single merged `ts` field via a static
land/sea-fraction diagnostic. The v1.1.0.0 release published 2026-03-30
exposes `sst`, `lsat`, and `ts` as separate variables per member and
extends the time axis through 2025-12, eliminating both prior
complications. The build now uses v1.1.0.0 directly.)

### 2.5 NOAA Land — deterministic best estimate + DCLSAT-donor uncertainty (members 1–100)

NOAA's v6.1 aravg text drops parametric uncertainty, and the gridded
NetCDF has only the merged-product anomaly with no land-only
uncertainty field. We therefore:
- Use the aravg text `anomaly` column as the deterministic central estimate.
- Add a donor uncertainty ensemble to it. The donor is the **first 100
  members of DCLSAT's 200-member ensemble** (Section 2.4), demeaned (so
  the donor's ensemble-mean is removed at each time step), and rescaled
  to match a target 1σ envelope informed by GloSATLAT's σ time series.
  Rationale: NOAA Land uses GHCN-M v4 stations and pairwise
  homogenization (PHA), which is far closer to DCLSAT's station base
  and method than to CRUTEM5's reduced-space optimal-interpolation
  gridding. Using DCLSAT as the donor preserves NOAA-family temporal
  structure of coverage uncertainty (especially pre-1900) better than
  CRUTEM5 would.
- **Why members 1–100 and not all 200?** Two of the six LSAT leaves
  (NOAA Land and C-LSAT 2.1, Section 2.6) are donor-imputed from
  DCLSAT. To prevent the two leaves' per-member trajectories from
  sharing identical DCLSAT-shape correlations, we partition the 200
  DCLSAT members into disjoint halves — NOAA Land draws from members
  1–100, C-LSAT 2.1 from 101–200. Per-member shape correlation between
  the two donor-imputed leaves is zero by construction, while each
  leaf's ensemble still represents the full DCLSAT-shape family at
  half the member count.

NOAA Land contributes **no independent uncertainty information** to
the final ensemble — only its best-estimate trajectory adds signal.
We keep NOAA in the family tree at full method-family weight (1/5);
its donor-imputed uncertainty is flagged in the limitations section.

### 2.6 C-LSAT 2.1 — deterministic best estimate + DCLSAT-donor uncertainty (members 101–200)

C-LSAT 2.1 (Sun, Li, et al., Sun Yat-sen University; updates of the
Xu et al. 2018 "C-LSAT" pipeline; hosted at `gwpu.net` and figshare
[doi:10.6084/m9.figshare.28255394](https://doi.org/10.6084/m9.figshare.28255394))
contributes a distinct method-family lineage to the LSAT side:
station-level homogenization on a merged station network that
substantially expands on GHCN-monthly's Asian and African coverage
(~25,000 stations across the updated network including non-GHCN
Chinese / Russian / Latin-American collections). The structural
distinction from NOAA PHA is real: C-LSAT uses MASH-style penalized
multiple-reference tests rather than NOAA's pairwise SNHT cascade,
and the input network differs by ~20% pre-1950.

**Pre-integration diagnostic.** Decadal r(C-LSAT, CRUTEM5) over
1860–1900 is 0.31–0.86 (mean 0.72), comparable to the pairwise r
between any two other LSAT products over the same era (range 0.63–0.76).
C-LSAT's pre-1900 information is therefore *not* a duplicate of CRUTEM5,
and the reviewer's "shared-input hollow-diversity" threshold (r > 0.95
pre-1900) is comfortably avoided.

**Construction (mirrors NOAA Land):**
- Read the deterministic monthly anomaly from
  `C-LSAT2.1_global_monthly.csv` (derived by
  `prepare_c_lsat_global_mean.py` as the cos(lat)-weighted global
  land-area mean of the 5°×5° gridded `tavg_anomaly` field; native
  baseline 1961-1990 is removed during the 1981-2010 re-baselining of
  Step 3).
- Add the **DCLSAT members 101–200** donor ensemble (Section 2.5
  rationale) rescaled to the same GloSATLAT-derived target σ envelope
  used for NOAA Land.

C-LSAT 2.1's contribution to the **1850–1900 reference period** shifts
the ensemble-mean baseline by only +0.007 °C in practice (measured
empirically by re-running the build with and without C-LSAT). The
realised shift is well below the 0.05 °C concern threshold flagged
at stats review.

## Step 3 — Re-baseline all six per-dataset ensembles to 1981–2010

For every member of every dataset:
- subtract that member's own 1981–2010 mean.

After this step all six ensembles share zero mean over 1981–2010 and
are directly comparable.

## Step 4 — Family tree weighting

We use a 5-method-family tree with equal P=1/5 at each top-level
branch (Thorne 2026 §3.2.5 "equal-weight method families" principle).
The five method families are:

- **HOMOGENIZATION / SCALPEL** — Berkeley Earth's iterative scalpel
  homogenization.
- **PAIRWISE PHA** — NOAA Land's pairwise SNHT cascade on GHCN-M v4.
- **CRU-LINEAGE** — CRUTEM5 + GloSATLAT (reduced-space optimal
  interpolation; common bias-correction lineage). Split internally
  to two leaves at P=1/10 each.
- **DYNAMICAL-CONSTRAINT** — DCLSAT (DCENT v3.0 land field; joint
  SST/LSAT energy-balance correction).
- **CMA-HOMOGENIZATION** — C-LSAT 2.1 (Sun et al., MASH-style penalised
  multiple-reference homogenization on a merged station network that
  includes a distinct ~20% pre-1950 non-GHCN input).

```
LSAT root (P=1)
├── HOMOGENIZATION / SCALPEL family    (P=1/5)
│   └── Berkeley Earth                              (P=1/5)
├── PAIRWISE PHA family                 (P=1/5)
│   └── NOAA Land                                   (P=1/5)
├── CRU-LINEAGE family                  (P=1/5)
│   ├── CRUTEM5                                     (P=1/10)
│   └── GloSATLAT                                   (P=1/10)
├── DYNAMICAL-CONSTRAINT family         (P=1/5)
│   └── DCLSAT                                      (P=1/5)
└── CMA-HOMOGENIZATION family           (P=1/5)
    └── C-LSAT 2.1                                  (P=1/5)
```

Resulting leaf probabilities (sum to 1):
- Berkeley Earth: **0.20**
- NOAA Land: **0.20**
- DCLSAT: **0.20**
- C-LSAT 2.1: **0.20**
- CRUTEM5: **0.10**
- GloSATLAT: **0.10**

**Time-varying weights for 2022–2025:** GloSATLAT does not cover 2022
onward, so for those years we redistribute its 0.10 probability mass
to CRUTEM5 (the only other leaf in the CRU-lineage family) — i.e.
CRUTEM5 becomes 0.20 for years 2022–2025 and GloSATLAT 0.0. All other
weights unchanged. This preserves the family-level weights (CRU
lineage remains P=0.20 throughout).

**Planned sensitivity tests (out of scope for v2):**
- Collapse PAIRWISE-PHA and CMA-HOMOGENIZATION into a single
  "station-network homogenization" super-family (joint P=1/4, NOAA &
  C-LSAT each P=1/8). The reviewer's most important recommended test —
  surfaces whether the C-LSAT 1/5 weight is buying genuine structural
  diversity or partially duplicating the NOAA PHA family.
- Half-weight the two donor-imputed leaves (NOAA Land and C-LSAT 2.1
  each at 0.10 instead of 0.20) to reflect their no-independent-σ
  status; redistribute the freed weight equally to the other three
  families.
- Place DCLSAT in PHA-family with NOAA (a strict reading of
  station-base overlap).
- Run with a 1995/96 splice to quantify the spread inflation we forgo
  by not splicing.

## Step 5 — Generate 10,000-member ensemble

For each of 10,000 draws:
1. Pick a leaf dataset using year-dependent weights from Step 4.
2. From that leaf's already-rebaselined ensemble (Berkeley Earth 10,
   CRUTEM5 / GloSATLAT / DCLSAT 200 each, NOAA Land 100, C-LSAT 2.1
   100), sample one member uniformly at random.
3. Take that member's full 1850–2025 annual time series.
4. Record it as one ensemble member.

The result is a 176 × 10,001 CSV (year + 10,000 anomaly columns), with
each column being a complete-period draw. Year-by-year statistics
(mean, 2.5/97.5 percentiles, SD) are then computed for plotting.

Note: unlike Thorne's GMST pipeline we do **not** splice tail and head
at 1995/96 — for LSAT alone the head/tail distinction matters less
since all leaves cover both, and splicing artificially decorrelates
pre-1996 and post-1995 uncertainty (Thorne flags this as a limitation
of his own approach in his Supplement §3).

## Known limitations (resolved from stats review)

1. **Berkeley Earth's 10 members propagate into 25% of the final
   ensemble.** Tail percentiles (especially 2.5/97.5) will show
   step-clustering at fewer than 10 distinct values in regimes where
   BE dominates. This is informative, not a bug — it advertises BE's
   under-sampling.
2. **CRUTEM5/GloSATLAT bias term is a Gaussian-scalar fallback.**
   CRUTEM5 v5.1 publishes only the bias CI, not realisations, so we
   cannot sample realisations directly. This likely overstates
   low-frequency bias uncertainty (treats it as perfectly correlated
   in time) while understating sub-annual structure.
3. **DCLSAT marginal-LSAT spread is conditioned on SST/LSAT joint
   consistency** and is plausibly ~10–30% narrower than a standalone
   DCLSAT product's uncertainty would be.
4. **Two of six leaves contribute no independent uncertainty.** NOAA
   Land and C-LSAT 2.1 are donor-imputed from DCLSAT (members 1–100 and
   101–200 respectively, so the two leaves do not share per-member
   shape correlations). We keep their best-estimate signals at full
   1/5 family weight; this means ~40% of the final ensemble has
   imputed-DCLSAT-shape uncertainty wrapped around two different best
   estimates. Combined with DCLSAT's own 1/5 leaf, ~60% of the LSAT
   ensemble's uncertainty traces in some form back to DCLSAT — flagged
   here and addressed in the planned half-weight sensitivity test.
5. **Effective N in the tail is far below 10,000.** With Berkeley
   Earth's 10 native trajectories sitting under P=0.25 of the tree,
   the effective independent draws in the tail (1850–1900) is
   plausibly closer to 200–500 once all leaves are accounted for.
   Percentile-based p-values should not be used for trend
   significance — the spread is structural, not stochastic.
6. **Bias toward over-coverage low-frequency / under-coverage
   high-frequency** in CRUTEM5/GloSATLAT synthesis (Section 2.2).

## Safe vs unsafe uses

**Safe:**
- Plotting central estimate + 90% band.
- Reporting the ensemble mean as the LSAT best estimate.
- Comparing with the sister ocean SST ensemble for context.

**Unsafe:**
- Reporting tail-period 2.5/97.5 percentiles as if calibrated
  frequencies.
- Trend-significance p-values from member-percentile counts.
- Per-leaf attribution without documenting the tree weights and
  running the sensitivity tests listed in Step 4.

## File outputs

- `land_ensemble.csv` — 176 × 10,001 (year + 10,000 members), 1850–2025,
  annual anomalies on a 1981–2010 baseline (each member centred at zero
  over 1981–2010 by construction).
- `land_ensemble_summary.csv` — annual mean, 2.5 / 50 / 97.5 percentiles
  and 1σ from the 10,000-member ensemble.
- `land_ensemble_perdataset.csv` — per-dataset annual best estimate and
  1σ across that dataset's own (200-or-10-member) ensemble.
- `land_ensemble.png` — central estimate + 95% band, with the five
  individual best estimates overlaid for sanity-checking.
- `family_tree_land.png` — visual summary of the family tree.

## References

- Thorne, P. W. *et al.* (2026), *A framework for an operational
  assessment of observed global warming-to-date*, ESSD preprint
  https://doi.org/10.5194/essd-2025-825 — Sections 3.2.5, and esp.
  Figures 24–25 for the family-tree GMST construction.
- Morice, C. P. *et al.* (2021), *An updated assessment of near-surface
  temperature change from 1850: the HadCRUT5 data set*, JGR-Atmos
  126(3), e2019JD032361 — for the correlated/uncorrelated/bias σ
  decomposition we re-apply to CRUTEM5 and GloSATLAT.
- Chan, D., Gebbie, G., Huybers, P. & Kent, E. C. (2024), *DCENT:
  Dynamically Consistent ENsemble of Temperature*, Scientific Data 11,
  article 612 — for the dynamical-consistency LSAT/SST joint correction
  underlying DCLSAT.
- Taylor, M. *et al.* (2025), *GloSAT land air temperature*, GDJ —
  for the LATsdb / LEK normals / exposure-bias adjustment used by
  GloSATLAT.
- Rohde, R. A. *et al.* (2013, 2020 updates), Berkeley Earth methodology
  papers — for the scalpel-based homogenization and kriging.
- Yin, X. *et al.* (2024, 2025), NOAAGlobalTemp v6/v6.1 documentation —
  for the LSAT analysis used in NOAA Land.

# Ocean-only SST ensemble — methodology

This note describes how we build a 10,000-member ensemble of annual
global-mean **sea-surface temperature (SST)** anomalies, approximately
replicating the family-tree-with-imputation approach used for GMST by
Thorne et al. (2026, ESSD discussion paper 2025-825) but applied to
ocean-only datasets and with a smaller catalogue. It mirrors the
sister methodology document `land_ensemble_methodology.md` for LSAT
and should be read alongside it.

The goal is to produce an SST analogue of the SST-grouped Thorne
ensemble so that, downstream, an SST ensemble × LSAT ensemble ×
land-fraction weighting can reconstruct a GMST ensemble whose
structural uncertainty is decomposed by component. This memo concerns
*only* the SST side.

## Input datasets (SST, qualifying products in May 2026)

| # | Dataset | Native ensemble? | Uncertainty form | Time span | Native baseline | File(s) |
|---|---------|------------------|------------------|-----------|-----------------|---------|
| 1 | HadSST 4.2.0.0 | 200 native members (gridded; reduced to 200 global means in `prepare_hadsst_global_ensemble.py`) | ensemble of global means we compute ourselves; central + decomposed σ also available in CSV | 1850-01 – 2026-03 | 1961-1990 | `Ocean Data/HadSST4/HadSST.4.2.0.0_monthly_GLOBE.csv` (central) + 4 ensemble zips (members 1-50, 51-100, 101-150, 151-200) totalling ~1.7 GB |
| 2 | ERSSTv6 (NOAA pre-release ensemble + aravg central) | 1000 native members (gridded; reduced to 1000 global means in `pull_ersstv6_members.py`) | ensemble of global means we compute ourselves; aravg deterministic retained as central anchor for years past native end | 1850-01 – 2024-12 (native); 1854-01 – 2026-04 (aravg central) | 1991-2020 | `Ocean Data/ERSSTv6/ERSSTv6_monthly_ensemble.csv` (1000 members) + `Ocean Data/ERSSTv6/aravg.mon.ocean.90S.90N.v6.1.0.202604.asc` (central) |
| 3 | COBE-SST 2 | No | deterministic (analysis-error field exists in gridded but is not in the global mean) | 1850-01 – 2026-04 | 1991-2020 | `Ocean Data/COBE-SST2/COBE-SST2_global_monthly.csv` (derived by `prepare_cobe_global_mean.py` from `sst.mon.mean.nc` − `sst.mon.ltm.1991-2020.nc`) |
| 4 | DCENT SST (DCENT v3.0 `sst` field) | 200 members (gridded; we reduce to 200 global means in `pull_dcent_sst_members.py`) | ensemble of global means we compute ourselves | 1850-01 – 2025-12 | 1982-2014 | `Ocean Data/DCENT SST/DCENT_SST_monthly_ensemble.csv` |

### Note on COBE versions

JMA's operational product is COBE-SST2 (released Dec 2021, updated
monthly). COBE-SST3 (Ishii et al. 2025) exists as a research product
but is not yet operationally distributed for monthly climate
monitoring; we use COBE-SST2.

### Note on HadSST4 ensemble provenance

The Met Office distributes the HadSST4 200-member ensemble as gridded
5°×5° NetCDF in four 50-member zips totalling ~1.7 GB. We download
all four, area-mean each member's monthly grids to produce 200
global-mean monthly series, and save a single CSV
(`HadSST.4.2.0.0_global_ensemble_monthly.csv`). This preserves the
native inter-component covariance structure of the 200 trajectories,
which would be lost in a synthesize-from-σ approach. We then **add**
measurement+sampling uncertainty (the uncorrelated and correlated σ
columns from the CSV summary, which are not represented in the bias
ensemble) by AR(1) noise per member — see Step 2.1 for the recipe.

## Step 1 — Common time/baseline grid

1. Reduce every dataset to **annual means** (Jan-Dec) of the
   global-mean monthly anomaly. For monthly-uncertainty datasets
   without native ensemble we propagate uncertainty to the annual mean
   explicitly (see Step 2).
2. Truncate to **1850–2025** (the latest complete calendar year as of
   working date 2026-05-12).
3. Re-baseline each dataset (and each ensemble member) to a **zero
   mean over 1981–2010** — matching the LSAT methodology and the
   Thorne choice.

### Note on sea-ice region handling

We do **not** impose a consistent treatment of sea-ice-covered cells
across the four SST products. Each product's native convention
propagates through the area-weighted global-mean step in the prep
scripts (`prepare_hadsst_global_ensemble.py`, `pull_dcent_sst_members.py`,
`pull_ersstv6_members.py`, `prepare_cobe_global_mean.py`), all of which
use a `cos(lat)` weighting over `isfinite` cells (plus an explicit
`> -50 °C` undef filter for ERSSTv6's GrADS `-999.9` sentinel).

| Product | Sea-ice handling in source | Effective handling in our build |
|---|---|---|
| **HadSST 4.2.0.0** | NaN over ice (in-situ only — no obs there) | **ice cells excluded** from global mean |
| **ERSSTv6** | `-999.9` undef sentinel over fully-ice cells; EOF reconstruction extends into partial-ice cells | **fully-ice excluded; partial-ice included as reconstructed SST** |
| **COBE-SST 2** | Fills with **−1.8 °C** (freezing point of seawater) wherever ice is present (verified empirically: 53 Arctic cells exactly at −1.8 °C in Jan 2024) | **ice cells INCLUDED, with the −1.8 °C placeholder** |
| **DCENT-I SST** | Full-globe kriging infill; sea-ice cells filled with 2-m air-temperature anomalies blended in (DCENT-I infilling design) | **ice cells INCLUDED via weight=0 → contribute their blended air-temp value to the sea-fraction-weighted SST mean** |

The post-v2 configuration (HadSST4 + ERSSTv6 exclude, COBE-SST 2 +
DCENT-I include) gives 2-of-4 symmetry on this axis — improved over
the v1 configuration (1-of-4 include). This creates two systematic
asymmetries that are not nuisance — they are part of the structural
diversity our family tree is meant to represent — but worth surfacing:

1. **Different denominators.** COBE's and DCENT-I's global-ocean
   denominators are constant in time (they include ice area at
   −1.8 °C placeholder for COBE, blended air-temp for DCENT-I); the
   other two products' denominators *grow* as polar ice retreats
   (previously-NaN cells become observable open ocean). As Arctic ice
   has receded since ~1980, the modern-era global means of HadSST4 and
   ERSSTv6 pick up newly-emerged cold-water cells that COBE/DCENT-I
   were already counting at their respective placeholders — exerting
   a small modern-era cooling bias on those two relative to COBE/DCENT-I.
2. **Edge-of-ice-retreat sign reversal.** Cells transitioning from
   ice-covered to ice-free move from {placeholder | NaN} to their
   actual seasonal-mean SST. For COBE/DCENT-I the transition appears
   as a real anomaly relative to the placeholder; for HadSST4/ERSSTv6
   the transition appears as a new cell entering the denominator.
   These produce different signals at the same geographic location.

**Magnitude:** in practice this contributes ~0.01–0.05 °C of
cross-product disagreement in the modern era, comparable in size to
the bias-correction uncertainty within HadSST4's own native bias
ensemble. We do not attempt to harmonise this asymmetry; the
family-tree treats it as structural diversity. If a future revision
wanted to *remove* it, the cleanest route would be to re-aggregate
each product on a common ice-aware mask (e.g., HadISST ice fraction)
and apply a consistent rule for ice-cell treatment — that would
collapse this dimension of cross-product spread but at the cost of
suppressing genuine methodological diversity.

## Step 2 — Build a per-dataset ensemble

We use native counts where available (HadSST4 200 native bias members
+ noise terms, DCENT 200, ERSSTv6 1000 native through 2024). The one
remaining deterministic-only product (COBE-SST2) is handled via a
HadSST4 donor (Section 2.4).

### 2.1 HadSST4 — area-mean 200 native gridded members + add measurement/sampling noise

The 200 HadSST4 ensemble members published as gridded 5°×5° monthly
NetCDFs sample the **bias-adjustment** dimension of uncertainty (bucket
vs. engine-room vs. buoy bias model parameters). They do NOT include
measurement/sampling uncertainty (`uncorrelated` and `correlated` σ
in the CSV summary), nor coverage uncertainty.

We construct the per-dataset ensemble in two parts:
1. **Native bias ensemble.** Download the 4 zips of 50 members each,
   unzip, and for each of the 200 NetCDFs compute the area-weighted
   (cos lat) global ocean mean per month (skipping NaN cells). This
   yields 200 monthly global-mean anomaly time series 1850-01 to
   2026-03 directly from the published HadSST4 ensemble — preserving
   the native temporal covariance of the bias dimension.
2. **Add measurement + sampling + coverage noise.** To each member's
   monthly series, add Gaussian noise terms reproducing the
   uncorrelated, correlated, and coverage σ from the CSV:
   - uncorrelated: AR(1) noise (ρ = 0.3 month-to-month) with marginal
     σ = `u_unc(t)`,
   - correlated: per-member scalar α ~ N(0,1) × `u_corr(t)`,
   - coverage:  per-member scalar β ~ N(0,1) × `u_cov(t)`.
   These additions match the Morice-style decomposition used for
   CRUTEM5/GloSATLAT/HadSST4 in the literature; together with the
   native bias ensemble they reconstruct HadSST4's full uncertainty
   envelope as documented by Kennedy et al. (2019).

Annual means are computed from each 12-month member sequence after
both ingredients are summed.

After this step the resulting envelope should match HadSST4's
published annual total σ to within a few percent (the only
approximations are the ρ=0.3 AR(1) choice and the perfect-time-
correlation assumption for `correlated` and `coverage`).

### 2.2 DCENT SST — pull 200 native DCENT-I members, extract SST global means

We use **DCENT-I v1.1.0.0** (Chan et al. 2026, Geoscience Data Journal;
doi:10.7910/DVN/ROG38Q, dataset version updated 2026-03-30 through end
of 2025) — the spatially complete kriging-infilled extension of DCENT.
The original unfilled DCENT v3.0 derived data are archived under
`Ocean Data/DCENT SST/v3_archive/` for reproducibility and reversion.

The DCENT-I v1.1.0.0 Harvard Dataverse release exposes 200 gridded
NetCDFs (~42 MB each, ~8.4 GB total). Each member NetCDF contains
**three separate fields** — `sst` (sea_surface_temperature_anomaly),
`lsat` (land near-surface air temperature anomaly), and `ts` (merged
Surface Temperature Anomaly) — on a 5°×5° monthly grid from 1850-01
through 2025-12.

We download each member NetCDF, compute the area-weighted (cos lat)
global mean of the `sst` field directly (skipping NaN cells), then
delete the raw NetCDF. The same script writes the parallel LSAT CSV
from each member's `lsat` field (land methodology §2.4). Peak
transient disk: ~42 MB.

**Sea-ice cells in DCENT-I's `sst` field.** DCENT-I's `sst` is finite
over both open ocean AND sea-ice cells; the value at an ice cell is
the kriged surface temperature including air-temperature blending.
This is methodologically the same compromise COBE-SST 2 makes via its
−1.8 °C freezing-point placeholder, just smarter — and it improves
the sea-ice asymmetry across SST products (see Step 1 sea-ice note:
DCENT-I now joins COBE-SST 2 as a product that includes ice cells in
the SST mean, while HadSST4 and ERSSTv6 remain ice-excluding).

**Why we switched from DCENT v3.0 to DCENT-I:**
- DCENT-I is the published, peer-reviewed infilled version.
- Resolves the "no extrapolation beyond observations" coverage gap
  in DCENT v3.0 — modern-era polar regions are now spatially complete.
- 200-member ensemble explicitly samples kriging uncertainty alongside
  the bias-adjustment uncertainty.
- Empirical impact on the headline: 2024 SST median shifts by
  +0.02 °C (modest modern-era polar warming uplift; matches paper's
  qualitative claim of polar-region warming attribution).

**Caveat #1 (carried over from DCENT v3.0):** each DCENT-I member's
SST-over-ocean field is paired with a specific LSAT field under the
underlying DCENT dynamical-consistency constraint. SST is the larger
reservoir and the larger uncertainty source in DCENT — its joint
constraint mainly prunes SST extremes inconsistent with the smaller
LSAT degrees of freedom, so the SST-marginal spread is plausibly
~5–15% narrower than a standalone SST product would show.

**Caveat #2 (new for DCENT-I):** the kriging-infilled modern-era
σ across 200 members is ~40% tighter than DCENT v3.0's was in
the same era (2024 σ: 0.0145 → 0.0088 °C); pre-1900 σ is comparable
(0.20 vs 0.19). This is a deliberate consequence of the infilling
model — kriging is a smooth low-variance estimator where data are
dense. The leaf still carries appropriate sparse-era uncertainty.

(The DCENT-I v1.0.0.0 release we initially planned to use covered
only 1850–2024 and required weight-based decomposition of a single
merged `ts` field via a static land/sea-fraction diagnostic. The
v1.1.0.0 release published 2026-03-30 exposes `sst`, `lsat`, and `ts`
as separate variables per member and extends the time axis through
2025-12, eliminating both prior complications.)

### 2.3 ERSSTv6 — pull 1000 native NCEI members, compute global means; frozen-offset fallback past native end

The NCEI pre-release ERSSTv6 ensemble (Huang et al., 2025; obtained from
`https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/tmp/ersstv6.ensemble/`)
exposes 1000 gridded 2°×2° monthly NetCDF/Fortran-binary fields, each a
full Monte-Carlo realization over ERSST's bias-correction parameter
space (`smult`, `stdmin`, `stdmax`, `mbias`, `adjs/b/a`, `wt1*`, etc. —
see `parameter.case.txt` in the source directory).

We download each member's `sst2d.ano.1850.2024.1991-2020.ens.NNNN.dat`
file (~135 MB each, Fortran sequential unformatted big-endian, 180×89
grid, 2100 monthly records), compute the area-weighted (cos lat) global
mean per month while masking the GrADS undef sentinel (-999.9), append
the resulting (NNNN'th, 2100) monthly series to an in-memory ensemble,
and **delete the raw .dat file before downloading the next member**.
This keeps peak disk use at ~135 MB rather than ~135 GB. The final
ensemble lands in `ERSSTv6_monthly_ensemble.csv` (2100 rows × 1002 cols:
year, month, m0001..m1000); see `pull_ersstv6_members.py`.

**Validation.** A 10-member subsample's annual median agreed with the
published ERSSTv6 best-estimate (aravg) to RMSE 0.030 °C (max abs 0.085
°C, mean bias +0.014 °C) over 1850–2024 — i.e., consistent with random
subsample noise. Inter-member σ at 2024 is ~0.025 °C, comparable to
HadSST4's 200-member σ at 2024 (~0.026 °C); pre-1900 σ grows to ~0.04
°C, capturing the expected sparse-era widening.

**Frozen-offset fallback past native end (boundary-year handling).**
The native ensemble ends Dec 2024 but the build runs through 2025.
For year y past the last native year y* (here y* = 2024):

  σ_scale(y)        = σ_HadSST4(y) / σ_HadSST4(y*)
  ensemble_median*  = nanmedian(native[y*, :])
  member_m(y)       = central(y) + (native[y*, m] − ensemble_median*) × σ_scale(y)

where `central(y)` is the published ERSSTv6 best-estimate from aravg
(re-baselined to 1981-2010 member-wise to match the rest of the build).

Each member's 2024 deviation from the ensemble median is carried
forward to 2025 (and rescaled by HadSST4's year-on-year σ ratio).
This preserves per-member trajectory continuity at the boundary —
the alternative of redrawing each member's 2025 deviation
independently would create ~0.05–0.10 °C discontinuities at 2024→2025
and inject spurious noise into per-member trend statistics ending in
2025. The σ-scaling factor captures any real coverage-driven σ growth
(e.g., the 2023–2025 drifter-fleet contraction shows up in HadSST4's
σ ratio).

The post-native end is therefore an *extrapolation*, not an
observation, and we **hard-fail** if the gap between the native
ensemble's last year and `END_YEAR` is ≥ 2. Frozen-offset is only
defensible for a single year; multi-year extrapolation requires
revisiting the design (e.g., AR(1) decay of the per-member offsets,
or restoring an explicit donor-scatter mechanism for years > y* + 1).
The threshold and the revisit-trigger message live in
`build_ocean_ensemble.py: extend_with_frozen_offset`.

### 2.4 COBE-SST2 — deterministic best estimate + HadSST4-donor uncertainty (with sparse-era σ inflation)

COBE-SST2 ships only a deterministic global-mean monthly time series
in the form we can read. For uncertainty we use a donor approach
analogous to NOAA Land in the LSAT methodology:

- Take HadSST4's native 200-member ensemble (Section 2.1),
  demean it (subtract the per-year ensemble mean), and rescale to
  match a target 1σ envelope. This gives a 200-member uncertainty
  ensemble that is anchored on HadSST4's structural uncertainty
  family.
- Add the rescaled donor to COBE-SST2's best estimate on the common
  annual axis.

**Choice of donor and rationale:** HadSST4 is the canonical in-situ
SST reconstruction with a published structural uncertainty
representation. Using HadSST4 as donor for COBE-SST2 is acceptable
because COBE-SST2 uses the same in-situ archive in the tail (it adds
satellite data in the modern era). COBE-SST2's analysis-error field
exists on the gridded NetCDF distributed by JMA/TCC but is **not
publicly distributed in the NOAA PSL mirror** we use; we therefore
cannot use it as a per-product target σ and fall back to HadSST4.
Documented limitation: COBE-SST2 contributes only its best-estimate
signal, not independent uncertainty information.

**Target σ for rescaling — sparse-era inflated.** The base target σ is
HadSST4's own annual σ. To partially compensate for the (uncaptured)
extra structural uncertainty that COBE brings from
optimal-interpolation reconstruction in the sparse pre-satellite era,
we **inflate the target σ by a year-dependent factor**:
- pre-1920: 1.3× HadSST4 σ
- 1920–1960: linearly decays from 1.3× to 1.0×
- post-1960: 1.0× HadSST4 σ

This is a deliberate compensation, not a measurement. The
sensitivity test version with inflation factor = 1.0× throughout
(identity rescaling) is reported alongside. Inflation is applied
only to COBE; the native ERSSTv6 ensemble already carries its own
sparse-era spread.

## Step 3 — Re-baseline all four ensembles to 1981–2010

For every member of every dataset (HadSST4 200, DCENT 200, ERSSTv6
1000, COBE-SST2 200):
- subtract that member's own 1981–2010 mean.

After this step all four ensembles share zero mean over 1981–2010 and
are directly comparable. The post-native-end frozen-offset fallback
for ERSSTv6 (Section 2.3) is applied *after* re-baselining, using the
re-baselined aravg central series as the 2025 anchor — so the
extrapolated 2025 members are also natively on the 1981-2010 baseline.

## Step 4 — Family tree weighting

We use an **equal-weight 4-leaf tree** (P=1/4 per leaf):

```
SST root (P=1)
├── HadSST family                              (P=1/4)
│   └── HadSST 4.2.0.0                                  (P=1/4)
├── ERSST family                               (P=1/4)
│   └── ERSSTv6                                         (P=1/4)
├── COBE family                                (P=1/4)
│   └── COBE-SST 2                                      (P=1/4)
└── DYNAMICAL-CONSTRAINT family                (P=1/4)
    └── DCENT SST                                       (P=1/4)
```

All four leaves cover both the tail (1850–1995) and head (1996–2025)
sub-periods, so there is no time-varying weight adjustment (unlike the
LSAT tree, which redistributes GloSATLAT's weight post-2021).

**Rationale.** As of v2 (May 2026), three of the four leaves carry
their own native uncertainty template: HadSST4 (200-member native bias
ensemble), DCENT SST (200-member dynamical-consistency ensemble), and
ERSSTv6 (1000-member native parameter ensemble from NCEI pre-release).
Only COBE-SST2 remains donor-imputed from HadSST4. Under the strict
"lineage of uncertainty" reading of Thorne 2026 §3.2.5 — which would
down-weight donor-imputed leaves — only COBE-SST2 would qualify for a
half-weight; the other three would each carry 1/4. We keep equal 1/4
weights at the leaf level for simplicity and because COBE-SST2's
best-estimate trajectory is still a structurally independent product
even with a borrowed σ template (the v1 stats-review recommendation to
half-weight no longer applies once two of the three previously-donored
leaves are native).

Effective independent uncertainty templates ≈ 3.5 (HadSST4, DCENT,
ERSSTv6 native; COBE shares a quarter of HadSST4's template). This
is a substantial improvement over v1 (~3.0 effective), where ~50% of
the final ensemble's uncertainty traced back to HadSST4. The
HadSST4-shaped contribution to the final 10k ensemble is now ~25%
(HadSST4's own 1/4 leaf), plus the COBE-SST2 1/4 leaf's donor-borrowed
spread = ~37.5% nominally HadSST4-shaped uncertainty, vs. ~75% in v1.

**Open issue (to revisit when additional ensembles become available):**
- COBE-SST2 analysis-error: the JMA/TCC archive exposes an analysis
  error field but only via gated forms; if extractable, switch
  COBE-SST2's target σ from "HadSST4 σ" to "COBE-derived σ".
- ERSSTv6 native ensemble: **resolved in v2** (1000-member NCEI
  pre-release acquired from Boyin Huang). The v1 plan to use the
  NOAAGlobalTempv5 500-member ensemble as a stand-in is no longer
  needed.
- Post-native-end fallback: the current frozen-offset is defensible
  for a single year (`gap == 1`). If ERSSTv6's native ensemble falls
  multiple years behind `END_YEAR`, the build hard-fails with a
  revisit message (`build_ocean_ensemble.py:
  extend_with_frozen_offset`). Re-pull the native ensemble or replace
  the fallback with a multi-year-capable design before continuing.

**Note on asymmetry with the LSAT tree.** The LSAT tree uses five
method families (HOMOG/PHA/CRU-lineage/Dynamical/CMA-homogenization)
each at P=1/5, while the SST tree uses four (HadSST/ERSST/COBE/Dynamical)
each at P=1/4. On the LSAT side the CRU-lineage family splits two ways
(CRUTEM5 + GloSATLAT) so its leaves carry P=1/10 each; all other LSAT
leaves carry P=1/5. The dynamical-constraint family is the only one
present on both sides (DCLSAT at P=1/5 on land; DCENT SST at P=1/4 on
ocean). This asymmetry is intentional — the SST catalogue lacks an
analogue of the CMA-homogenization family — and should not be
"corrected" when SST × LSAT are recombined into a GMST ensemble:
draws should still be independent across sides.

**Planned sensitivity tests (out of scope for v2):**
- **Half-weight COBE-SST2** (the only remaining donor-imputed leaf:
  alternative tree HadSST4 2/7, ERSSTv6 2/7, DCENT 2/7, COBE 1/7).
  Defensible interpretation under a strict "lineage of uncertainty"
  reading of Thorne.
- Single in-situ-only sub-family containing HadSST4, ERSSTv6,
  COBE-SST2 at P=1/2 (split 1/3, 1/3, 1/3 within), DCENT at P=1/2 —
  the most aggressive version, lumping all three in-situ-based
  products as one lineage.
- Sparse-era σ inflation factor = 1.0× for COBE (no inflation,
  identity rescaling) — to quantify how much the inflation
  contributes to tail spread.
- Frozen-offset fallback variants: AR(1) decay of the per-member
  offset across the fallback year(s), or restoring an explicit
  donor-scatter mechanism for the post-native-end year, to quantify
  the fallback-design sensitivity.
- Splice at the 1995/96 midpoint of the 1981–2010 reference period
  (matching Thorne's GMST approach).

## Step 5 — Generate 10,000-member ensemble

For each of 10,000 draws:
1. Pick a leaf dataset using the leaf weights from Step 4.
2. From that leaf's already-rebaselined ensemble (HadSST4 200, DCENT
   200, ERSSTv6 1000, COBE-SST2 200), sample one member uniformly at
   random.
3. Take that member's full 1850–2025 annual time series.
4. Record it as one ensemble member.

Heterogeneous leaf sizes (1000 vs 200) are handled by the
`build_ocean_ensemble.py` sampler — it draws independently from each
leaf with the leaf-specific member count.

The result is a 176 × 10,001 CSV (year + 10,000 anomaly columns).
Year-by-year statistics (mean, 2.5/97.5 percentiles, SD) are then
computed for plotting.

Like the LSAT case, we do **not** splice tail and head — all four
leaves cover the full record.

## Known limitations

1. **One of four products contributes no independent uncertainty in v2.**
   COBE-SST2 is donor-imputed; only its best estimate is independent.
   With equal 1/4 leaf weights, ~25% of the final ensemble's
   uncertainty originates from a HadSST4-shaped donor spread wrapped
   around COBE-SST2's best estimate. (ERSSTv6, HadSST4, DCENT all
   carry their own native ensembles in v2; this is down from ~75%
   nominally-HadSST4-shaped uncertainty in v1.)
2. **HadSST4 bias dimension is native; measurement / correlated /
   coverage dimensions are layered on as Gaussian noise** with the
   AR(1) and perfect-time-correlation approximations from the LSAT
   methodology. This is the same Morice-style shortcut applied to all
   our products with decomposed σ.
3. **DCENT SST marginal spread is conditioned on SST/LSAT joint
   consistency** and is plausibly ~10–30% narrower than a standalone
   SST product's uncertainty would be (same caveat as DCLSAT).
4. **No splicing** narrows spread relative to a Thorne-style spliced
   version; defensible here because all leaves cover the full record.
5. **Effective N is much less than 10,000** in any era where the
   ensemble draws collapse to a single product (e.g. if HadSST4 and
   DCENT SST agree closely with each other, then half of the 10k are
   members of those two ensembles and the spread is set by ~400
   distinct trajectories).
6. **HadSST4 still influences ~37.5% of the v2 ensemble** when you
   count its own 1/4 leaf plus COBE-SST2's HadSST4-donored 1/4 leaf —
   so any HadSST4-specific shape (e.g. its bias-CI assumption about
   bucket vs. engine-room mixtures) still bleeds into roughly a third
   of the ensemble. This is down from ~75% in v1 but remains material.
7. **ERSSTv6 boundary year (2025) is an extrapolation, not an
   observation.** Each 2025 member inherits its 2024 deviation from
   the 1000-member median (σ-scaled by HadSST4's year-on-year σ
   ratio). Members maintain trajectory continuity at the boundary,
   but late-year σ growth from non-coverage sources (e.g., a sudden
   bias-correction revision specific to 2025) is not captured. Hard-
   fails if the gap to `END_YEAR` reaches 2 years.
8. **Sea-ice region handling is not harmonised across products**
   (see Step 1 "Note on sea-ice region handling" for the per-product
   convention table). After the v2 switch to DCENT-I, two products
   include ice cells (COBE-SST 2 at −1.8 °C; DCENT-I via blended
   2-m air-temp through its weight-based decomposition) and two
   exclude them (HadSST4, ERSSTv6). This 2-of-4 symmetry is an
   improvement over v1's 1-of-4 configuration. Cross-product
   disagreement from sea-ice convention contributes ~0.01–0.05 °C
   to modern-era spread, absorbed by the family tree as structural
   diversity rather than reconciled.

## Safe vs unsafe downstream uses

**Safe:**
- Plotting central estimate + 90% band; comparing with the LSAT
  ensemble.
- Combining SST ensemble × LSAT ensemble × land-fraction to
  reconstruct a structurally-attributed GMST ensemble, **provided** the
  draws are independent (see GMST combination notes below).
- Reporting the ensemble mean as the SST best estimate.

**Unsafe:**
- Reporting tail-period 2.5/97.5 percentiles as if calibrated
  frequencies.
- Trend-significance p-values from member-percentile counts.
- Per-leaf attribution without running the sensitivity tests in Step 4.
- Constructing a GMST ensemble that attempts to *pair* DCENT SST and
  DCLSAT members by index. DCENT's joint dynamical consistency is
  destroyed once members are sampled independently by the family-tree
  sampler (most of the 10,000 final SST members come from non-DCENT
  leaves anyway). Do not assume the SST and LSAT files retain DCENT's
  joint pairing — they do not.

## GMST combination notes (cross-cutting)

When the SST and LSAT 10,000-member ensembles are combined to produce
a structural-uncertainty GMST ensemble, draws should be made
**independently** from each file:
- LSAT uses **DCLSAT** as donor for its imputed leaf (NOAA Land).
- SST uses **HadSST4** as donor for its one remaining imputed leaf
  (COBE-SST2 in v2).
- No donor product is shared between LSAT and SST, so cross-component
  imputation correlation is zero by construction.

Therefore: pair an LSAT member (column m_i in `land_ensemble.csv`)
with an SST member (column m_j in `ocean_ensemble.csv`) by drawing
both indices uniformly and independently from {1..10000}. Do not
attempt to align by member index — the indices have no shared
provenance.

## File outputs

- `ocean_ensemble.csv` — 176 × 10,001 (year + 10,000 members), 1850–2025,
  annual anomalies on a 1981–2010 baseline (each member centred at zero
  over 1981–2010 by construction).
- `ocean_ensemble_summary.csv` — annual mean, 2.5 / 50 / 97.5 percentiles
  and 1σ from the 10,000-member ensemble.
- `ocean_ensemble_perdataset.csv` — per-dataset annual best estimate and
  1σ across that dataset's own 200-member ensemble.
- `ocean_ensemble.png` — central estimate + 95% band, with the four
  individual best estimates overlaid for sanity-checking.
- `family_tree_ocean.png` — visual summary of the family tree.

## References

- Thorne, P. W. *et al.* (2026), ESSD preprint 2025-825 — Section 3.2.5
  on the SST-grouped family tree.
- Kennedy, J. J. *et al.* (2019), *An ensemble data set of sea-surface
  temperature change from 1850: HadSST4*, JGR-Atmos 124(14) — for the
  200-member ensemble and uncertainty decomposition.
- Huang, B. *et al.* (2025), *NOAA ERSSTv6*, J. Climate Parts I & II —
  for the ERSSTv6 methodology; the 1000-member parameter ensemble
  (pre-release) is hosted at
  `https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/tmp/ersstv6.ensemble/`
  with contact `boyin.huang@noaa.gov`.
- Hirahara, S., Ishii, M., & Fukuda, Y. (2014), *Centennial-scale SST
  analyses*, J. Climate — for COBE-SST2 methodology.
- Chan, D. *et al.* (2024), DCENT, Scientific Data 11 — for the
  dynamical-consistency joint SST/LSAT correction.
- Morice, C. P. *et al.* (2021), HadCRUT5, JGR-Atmos — for the σ
  decomposition we reuse.

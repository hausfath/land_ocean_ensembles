# Ocean-only SST ensemble — methodology

This note describes how we build a 10,000-member ensemble of annual
global-mean **sea-surface temperature (SST)** anomalies, approximately
replicating the family-tree-with-imputation approach used for GMST by
Thorne et al. (2026, ESSD discussion paper 2025-825) but applied to
ocean-only datasets and with a smaller catalogue. It mirrors the
sister methodology document `land_ensemble_methodology.md` for LSAT
and should be read alongside it.

The goal is to produce an SST analogue of the SST-grouped Thorne
ensemble — a component-domain ensemble of global SST with
structural-method uncertainty captured by the family tree. This
document concerns *only* the SST side; the sister land methodology
is described in `land_ensemble_methodology.md`.

## Input datasets (SST, qualifying products in May 2026)

| # | Dataset | Native ensemble? | Uncertainty form | Time span | Native baseline | File(s) |
|---|---------|------------------|------------------|-----------|-----------------|---------|
| 1 | HadSST 4.2.0.0 | 200 native members (gridded; reduced to 200 global means in `prepare_hadsst_global_ensemble.py`) | ensemble of global means we compute ourselves; central + decomposed σ also available in CSV | 1850-01 – 2026-03 | 1961-1990 | `Ocean Data/HadSST4/HadSST.4.2.0.0_monthly_GLOBE.csv` (central) + 4 ensemble zips (members 1-50, 51-100, 101-150, 151-200) totalling ~1.7 GB |
| 2 | ERSSTv6 (NOAA pre-release ensemble + aravg central) | 1000 native members (gridded; reduced to 1000 global means in `pull_ersstv6_members.py`) | ensemble of global means we compute ourselves; aravg deterministic retained as central anchor for years past native end | 1850-01 – 2024-12 (native); 1854-01 – 2026-04 (aravg central) | 1991-2020 | `Ocean Data/ERSSTv6/ERSSTv6_monthly_ensemble.csv` (1000 members) + `Ocean Data/ERSSTv6/aravg.mon.ocean.90S.90N.v6.1.0.202604.asc` (central) |
| 3 | COBE-SST 3 (Ishii et al. 2025) | 300 native perturbation members (gridded 1°; reduced to 300 global means in `pull_cobe_sst3_members.py`) | (i) ensemble of global means we compute from the native perturbations + (ii) HadSST4-donor pseudo-ensemble wrapping the SST3 best estimate for years where native perturbations are unavailable | 1850-01 – 2024-12 (best estimate); 1870-01 – 2024-12 (native perturbations) | 1991-2020 (cell-wise, applied per member to absolute SST) | `Ocean Data/COBE-SST3/COBE-SST3_global_monthly.csv` (central) + `Ocean Data/COBE-SST3/COBE-SST3_monthly_ensemble.csv` (300 members) |
| 4 | DCENT SST (DCENT-I `sst` field, v1.1.0.0) | 200 members (gridded; we reduce to 200 global means in `pull_dcent_i_members.py`) | ensemble of global means we compute ourselves | 1850-01 – 2025-12 | 1982-2014 | `Ocean Data/DCENT SST/DCENT_I_SST_monthly_ensemble.csv` |

### Note on the COBE family

COBE-SST3 (Ishii et al. 2025) supersedes COBE-SST2 methodologically:
it introduces a maritime air temperature (NMAT) correction in the SST
bias-adjustment step that improves joint consistency with land
records, and distributes a native 300-member perturbation ensemble
representing bias-correction and NMAT-input uncertainty. We use
COBE-SST3 as the sole COBE-family product in the family tree, hosted
as two sibling variants within a single COBE-family branch:

- **COBE-SST3 native** — the 300-member perturbation ensemble, finite
  for 1870–2024.
- **COBE-SST3 donor** — the SST3 best estimate (1850–2024) wrapped in
  a HadSST4-shaped uncertainty envelope (Section 2.5 describes the
  donor recipe). This sibling provides COBE-family uncertainty in the
  1850–1869 sliver where SST3's native perturbations don't reach.

Both siblings share SST3's best-estimate trajectory; they differ only
in *how* uncertainty is represented (native bias-correction
perturbations vs. donor-imputed HadSST4-shaped spread). Inside
1870–2024 each carries P=1/8 of the tree; in 1850–1869 the donor
sibling carries the full 1/4 of the COBE-family branch via
sibling-redirect (Step 4). In 2025+ neither sibling has SST3 data and
the COBE family contributes 0, with the three remaining method
families (HadSST, ERSST, DCENT) renormalising to ~1/3 each — see Step
5 for the implications.

This structure follows Thorne et al. (2026) in two ways. First, the
within-branch arrangement of two sibling variants mirrors Thorne's
treatment of products whose native ensemble doesn't span the full
record. Second, accepting that a method family drops out of the head
sub-period rather than borrowing from a different product matches
Thorne's qualification rule that *"for the second sub-period, data
should exist from 1981 to the last full year inclusive"* — which would
itself disqualify COBE-STEMP3 from Thorne's 2025 head tree.

A draft email to Ishii clarifying perturbation coverage and asking
about the 2025+ extension plan lives at
`Ocean Data/COBE-SST3/draft_email_to_ishii.md`. The README states
"1890–2020 with post-2020 by request"; the archive itself (as verified
2026-05-22 across m001, m099, m150, m300) carries 155 yearly files per
member spanning 1870–2024.

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
scripts (`prepare_hadsst_global_ensemble.py`, `pull_dcent_i_members.py`,
`pull_ersstv6_members.py`, `prepare_cobe_sst3_global_mean.py`,
`pull_cobe_sst3_members.py`), all of which use a `cos(lat)` weighting
over `isfinite` cells (plus an explicit `> -50 °C` undef filter for
ERSSTv6's GrADS `-999.9` sentinel).

| Product | Sea-ice handling in source | Effective handling in our build |
|---|---|---|
| **HadSST 4.2.0.0** | NaN over ice (in-situ only — no obs there) | **ice cells excluded** from global mean |
| **ERSSTv6** | `-999.9` undef sentinel over fully-ice cells; EOF reconstruction extends into partial-ice cells | **fully-ice excluded; partial-ice included as reconstructed SST** |
| **COBE-SST 3** | Analysis returns finite SST under ice with the interpolated near-freezing-point value (verified Jan 2024: 456 cells in [−1.85, −1.75] °C; 0 cells at exactly −1.8); separate `sic` field distributed alongside `sst` | **ice cells INCLUDED with native near-freezing analysis values (no fixed placeholder)** |
| **DCENT-I SST** | Full-globe kriging infill; sea-ice cells filled with 2-m air-temperature anomalies blended in (DCENT-I infilling design) | **ice cells INCLUDED via weight=0 → contribute their blended air-temp value to the sea-fraction-weighted SST mean** |

The v3 configuration (HadSST4 + ERSSTv6 exclude; COBE-SST 3 + DCENT-I
include) gives 2-of-4 ice-cell-including symmetry. Each of the two
including products uses a different convention (SST3 interpolated
near-freezing values; DCENT-I blended air-temp), and the COBE family
itself is uniformly SST3-convention throughout (the donor sibling
wraps SST3's central, so its central also follows SST3's ice
convention; the donor's spread template is HadSST4-shaped but
demeaned, so HadSST4's exclude-ice convention does not propagate
through the spread). The two systematic asymmetries discussed below
remain — they are part of the structural diversity our family tree is
meant to represent — but worth surfacing:

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
+ noise terms, DCENT 200, ERSSTv6 1000 native through 2024, COBE-SST3
300 native perturbations through 2024). The SST3-donor sibling
(Section 2.5) wraps SST3's best estimate in a HadSST4-shaped envelope
for 1850–1869 coverage and as a second within-branch uncertainty
template inside 1870–2024.

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
COBE-SST 3 takes a related approach with interpolated near-freezing
values under ice. See the Step 1 sea-ice table: DCENT-I and COBE-SST 3
include ice cells in the SST mean (with different conventions), while
HadSST4 and ERSSTv6 remain ice-excluding.

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

### 2.4 COBE-SST3 native — area-mean 300 native perturbation members from JMA archive

COBE-SST3 (Ishii et al. 2025; doi:10.2151/jmsj.2024-010) replaces the
historical bias-correction model of COBE-SST2 with one that is jointly
consistent with the COBE-LSAT3 land record via a maritime air temperature
correction (NMAT). JMA's MRI distributes a 300-member perturbation
ensemble of the gridded analysis on a 1°×1° monthly grid (the 1° fields
are derived from the native 0.25° analysis). Members differ in
bias-correction parameters and the NMAT input.

We download the per-member yearly NetCDF files at
`https://climate.mri-jma.go.jp/pub/archives/Ishii-et-al_COBE-SST3/cobe-sst3/1x1/perturb/m{NNN}/`
in a streaming manner: for each of the 300 members we fetch its 155
yearly files (1870–2024, ~1 MB each, ~155 MB transient), compute the
per-cell 1991–2020 monthly climatology from that member's own series,
subtract climatology, area-weight to a global mean over ocean cells
(`sst` is finite), and delete the member's yearly files before moving on.
The output is `Ocean Data/COBE-SST3/COBE-SST3_monthly_ensemble.csv` with
1860 rows × 302 columns (year, month, m001..m300). The 1991–2020
climatology choice is a convenience — the downstream re-baselining to
1981–2010 (Step 3) makes the within-member climatology choice
algebraically irrelevant.

Total transfer: ~46.5 GB (300 × 155 MB). The deterministic best-estimate
field (`prepare_cobe_sst3_global_mean.py`, downloads the 1850–2024
yearly best-estimate files at ~370 MB total) is retained as a diagnostic
column in `ocean_ensemble_perdataset.csv` but is **not** used to anchor
the SST3 leaf in the family tree — the SST3 leaf enters via its 300
ensemble members only. This avoids any double-weighting of SST3 (once
through its central, once through its perturbations).

**Sea-ice handling.** As tabulated in Step 1, SST3 distributes finite
`sst` under ice with the analysis's interpolated near-freezing value
(verified Jan 2024: 456 cells in [−1.85, −1.75] °C; 0 at exactly −1.8).
Our prep script masks by `isfinite(sst)`, so ice cells enter the global
mean with that interpolated value — placing SST3 with DCENT-I in the
ice-cell-including group (Step 1).

**Native-window boundaries.** SST3 perturbations are NaN before 1870
and after 2024. Inside the family tree, the native-SST3 leaf is NaN
in those years; the NaN-aware sampler (Step 5) routes the 1850–1869
share to the SST3-donor sibling (Section 2.5) via sibling-redirect.
After 2024 both COBE-family siblings are NaN — the SST3 best estimate
itself stops at 2024 — and the COBE family contributes 0 to the
ensemble in 2025+, with the three other method families renormalising
to ~1/3 each. The working perturbation window is 1870–2024 (the
README's "1890–2020 with post-2020 by request" appears outdated
relative to the 2025-10-01 archive contents; see
`Ocean Data/COBE-SST3/draft_email_to_ishii.md` for an open clarifying
question to Ishii).

**Diagnostics performed at build time.** The build script logs the
per-leaf shares at 1860 / 1900 / 2020 / 2024 / 2025 and hard-fails if
the COBE family does not carry 1/4 of the tree in any year 1850–2024
or anything other than 0 in 2025+. A soft σ-step diagnostic on either
side of 1869→1870 and 2024→2025 surfaces any wildly mis-scaled spread.

**Note on the native ensemble's spread.** The 300-member native
perturbation ensemble exhibits a *tighter* global-mean spread than the
HadSST4-donor sibling across the full record (empirically: native 1σ
≈ 0.02 °C vs. donor 1σ ≈ 0.07–0.10 °C through the sparse era; native
1σ ≈ 0.01 °C vs. donor 1σ ≈ 0.007 °C in 2024 — donor tighter only in
the modern post-satellite era). Plausible contributions:

1. **Joint NMAT-LSAT consistency constraint.** SST3's bias correction
   is anchored to a maritime air-temperature (NMAT) field that itself
   couples to COBE-LSAT3. The 300-member perturbations sample within
   this joint constraint, pruning trajectories inconsistent with land
   air temperature — narrowing the SST-marginal spread (analogous to
   DCENT-I's joint-constraint caveat, Section 2.2, though arising
   through a different mechanism).
2. **Spatial smoothing damps cell-level perturbations.** Each member
   passes through the same OI/kriging analysis structure, which
   reduces the global-mean variance from per-cell perturbations.
3. **Per-cell analysis error is distributed separately.** Each
   gridded member ships a per-cell `err` field (1σ analysis error in
   K) that represents uncertainty not spanned by the 300-member
   spread. We do not currently propagate `err` into the global-mean
   ensemble — see the open issues list below.

The HadSST4-donor sibling therefore plays a real structural-uncertainty
role inside the SST3-native window: it samples bias-correction
parameter perturbations *outside* the joint-NMAT-LSAT-consistent
subspace that the native ensemble explores. The two siblings together
cover a broader uncertainty space than either alone, and the equal
1/8 + 1/8 weighting reflects that complementarity.

### 2.5 COBE-SST3 donor — SST3 best estimate + HadSST4-donor uncertainty (with sparse-era σ inflation)

The SST3-donor sibling wraps a HadSST4-shaped uncertainty envelope
around the SST3 best-estimate central, giving COBE-family
uncertainty coverage in the 1850–1869 sliver where SST3's native
perturbations don't reach and a methodologically independent
uncertainty representation (HadSST4 bias-correction parameter space)
alongside the SST3 native sibling inside 1870–2024.

The donor recipe is analogous to NOAA Land in the LSAT methodology:

- Take HadSST4's native 200-member ensemble (Section 2.1),
  demean it (subtract the per-year ensemble mean), and rescale to
  match a target 1σ envelope. This gives a 200-member uncertainty
  ensemble anchored on HadSST4's structural uncertainty family.
- Add the rescaled donor to the SST3 best-estimate central on the
  common annual axis.
- The resulting leaf is naturally NaN in any year where the SST3
  central is NaN (i.e., 2025+), so the donor sibling falls out of
  the sampler past the SST3 best-estimate end without needing an
  explicit gap flag.

**Choice of donor and rationale.** HadSST4 is the canonical in-situ
SST reconstruction with a published structural uncertainty
representation. Using HadSST4 as donor for COBE is methodologically
defensible because both share an in-situ archive lineage (COBE adds
satellite data in the modern era). COBE-SST3 ships its own analysis
error field (`err`) alongside the gridded best estimate, but we use
the 300-member native perturbation ensemble for SST3's primary
uncertainty representation (Section 2.4) — the donor sibling adds
*structural* uncertainty diversity by sampling from HadSST4's
bias-correction parameter space rather than from SST3's own. The
two siblings therefore represent the COBE family's uncertainty
through two distinct families of perturbations.

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
only to the SST3-donor sibling; HadSST4, ERSSTv6, DCENT, and the
SST3 native sibling all carry their own native sparse-era spread.

## Step 3 — Re-baseline all ensembles to 1981–2010

For every member of every dataset (HadSST4 200, DCENT 200, ERSSTv6
1000, COBE-SST3 native 300, COBE-SST3 donor 200):
- subtract that member's own 1981–2010 mean.

After this step all five leaves share zero mean over 1981–2010 and
are directly comparable. The post-native-end frozen-offset fallback
for ERSSTv6 (Section 2.3) is applied *after* re-baselining, using the
re-baselined aravg central series as the 2025 anchor — so the
extrapolated 2025 members are also natively on the 1981-2010 baseline.

## Step 4 — Family tree weighting

We use a **4-branch tree with two SST3 siblings under the COBE branch**
(P=1/4 per branch, COBE branch split 1/8 + 1/8 inside 1870–2024):

```
SST root (P=1)
├── HadSST family                                 (P=1/4)
│   └── HadSST 4.2.0.0                                  (P=1/4)
├── ERSST family                                  (P=1/4)
│   └── ERSSTv6                                         (P=1/4)
├── COBE family                                   (P=1/4 in 1850–2024; P=0 in 2025+)
│   ├── COBE-SST 3 native (300 perturbations)        (P=1/8 in 1870–2024; P=0 elsewhere)
│   └── COBE-SST 3 donor (HadSST4-wrapped central)   (P=1/8 in 1870–2024;
│                                                     P=1/4 in 1850–1869;
│                                                     P=0 in 2025+)
└── DYNAMICAL-CONSTRAINT family                   (P=1/4)
    └── DCENT SST (DCENT-I)                             (P=1/4)
```

All non-COBE leaves cover the full 1850–2025 record. The COBE family
has two boundaries: at 1870 (native perturbations begin; donor
sibling-redirect takes over below) and at 2024 (SST3 best estimate
ends; the family contributes nothing in 2025+). Both transitions are
handled automatically by the NaN-aware sampler without an explicit
splice (Step 5).

**Rationale for the two-sibling COBE branch.** COBE-SST3's 300-member
native perturbation ensemble represents bias-correction and NMAT-input
uncertainty internal to JMA's analysis. The SST3 donor sibling
contributes a *structurally different* uncertainty template — sampling
HadSST4's bias-correction parameter space rather than SST3's — wrapped
around the same SST3 best-estimate central. Within the COBE branch
the two siblings therefore explore two different families of plausible
bias-correction perturbations applied to a shared best estimate.
Equal 1/8 weighting treats them as comparable structural-uncertainty
contributions inside the SST3-perturbation window; in 1850–1869 the
donor sibling absorbs the COBE-family branch weight (1/4) via
sibling-redirect, ensuring the COBE family retains its full
tree-level contribution wherever SST3's best estimate exists.

**Why the COBE family drops out in 2025+.** SST3's best estimate ends
in December 2024. Both siblings depend on it (native via member-level
trajectories; donor via the central anchor), so both are NaN in 2025+.
Rather than imputing an SST3 central from another product (which would
falsely advertise SST3 information we don't have), we let the COBE
family drop to zero and let the other three method families
renormalise to ~1/3 each. This matches Thorne's qualification logic
for the head sub-period (products must cover through the last full
year to qualify) and avoids manufactured 2025 SST3 values.

**Independent best estimates vs independent uncertainty templates.**
These two counts are not the same and the v3 design separates them:

- **Independent best estimates:**
  - 1850–2024: 4 (HadSST4, ERSSTv6, DCENT, SST3 — both COBE siblings
    share SST3's best-estimate central).
  - 2025+: 3 (COBE family absent).
- **Independent uncertainty templates:**
  - 1850–1869: ~3 (HadSST4, ERSSTv6 native, DCENT native; the SST3
    donor sibling's spread is HadSST4-shaped → not independent of
    HadSST4).
  - 1870–2024: ~4 (above + SST3's native perturbation template).
  - 2025+: 3 (HadSST4, ERSSTv6 native, DCENT native; COBE family
    absent).

So v3's improvement is concentrated in the 1870–2024 window, which
adds a fourth independent uncertainty template alongside four
independent best estimates. Single-number summary: ~4.0 inside the
1870–2024 window, ~3.0 in 1850–1869 and 2025+ — splitting the
difference between the two definitions is fine for prose framing.

**HadSST4-shaped contribution to the final 10k ensemble:**
- 1870–2024: 25% (HadSST4 leaf) + 12.5% (SST3-donor sibling spread)
  = **37.5% nominally HadSST4-shaped**.
- 1850–1869: 25% (HadSST4 leaf) + 25% (SST3-donor sibling absorbing
  the full COBE branch via sibling-redirect) = **50% nominally
  HadSST4-shaped**.
- 2025+: 25% (HadSST4 leaf only; COBE family absent) renormalised over
  three contributing families = **33.3% nominally HadSST4-shaped**.
- For comparison: v1 was ~75% across the whole record; v2 was ~50%
  unconditionally.

So v3 reduces HadSST4-shape concentration relative to v2 in the
1870–2024 window (the longest stretch of the record); matches v2 in
the sparse 1850–1869 sliver; and *further* reduces it in 2025+ because
the COBE family — which donored HadSST4 — drops out entirely.

**Open issue (to revisit when additional ensembles become available):**
- COBE-SST3 best estimate + perturbations 2025+: **awaiting JMA's
  operational extension** (see `Ocean Data/COBE-SST3/draft_email_to_ishii.md`).
  When 2025+ data are distributed, extend `prepare_cobe_sst3_global_mean.py`'s
  YEAR_END (and `pull_cobe_sst3_members.py`'s, if perturbations
  are simultaneously extended) and rebuild — no methodology change
  needed; the COBE family will simply contribute to more years.
- COBE-SST3 native `err` field: each gridded perturbation member also
  ships an analysis-error field. We currently use the across-member
  spread to characterise SST3 uncertainty; reading the per-cell `err`
  field would add a second within-product uncertainty term (analogous
  to HadSST4's measurement+sampling σ on top of the bias ensemble).
  Lower priority.
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
ocean). The asymmetry is intentional — the SST catalogue lacks an
analogue of the CMA-homogenization family — and reflects the distinct
catalogues qualifying in each domain.

**Planned sensitivity tests (out of scope for v3):**
- **Information-weighted COBE-branch split** (alternative to the
  current equal 1/8 + 1/8): weight the SST3 native sibling = 2/3 of
  the branch and the SST3 donor sibling = 1/3, reflecting that the
  native sibling has SST3-internal uncertainty while the donor
  sibling borrows HadSST4's. Keeping equal 1/8 + 1/8 as default.
- **Half-weight SST3-donor sibling** (the only donor-imputed sibling
  in v3): alternative tree HadSST4 2/7, ERSSTv6 2/7, DCENT 2/7, COBE
  1/7 inside 1870–2024 (with COBE = 1/14 native + 1/14 donor).
  Defensible interpretation under a strict "lineage of uncertainty"
  reading of Thorne.
- **Continue COBE family into 2025 via the SST3-donor sibling alone**
  (re-introducing a COBE-SST2-style central): would require either
  resurrecting the COBE-SST2 central as a 2025-only anchor or
  extrapolating SST3's best estimate using a model. Quantifies the
  cost of dropping COBE in 2025.
- Single in-situ-only sub-family containing HadSST4, ERSSTv6, and
  COBE-SST3 (both siblings) at P=1/2 (split 1/3 / 1/3 / 1/3 within;
  COBE sub-split further), DCENT at P=1/2 — the most aggressive
  version, lumping all three in-situ-based products as one lineage.
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
   200, ERSSTv6 1000, COBE-SST3 native 300, COBE-SST3 donor 200),
   sample one member uniformly at random.
3. Take that member's full 1850–2025 annual time series.
4. Record it as one ensemble member.

Heterogeneous leaf sizes (1000 vs 300 vs 200) are handled by the
`build_ocean_ensemble.py` sampler — it draws independently from each
leaf with the leaf-specific member count.

**Sibling-redirect handles the COBE branch's internal boundaries.**
The COBE family has time-varying leaf availability:

- In 1850–1869 the native sibling is NaN (no SST3 perturbations
  before 1870).
- In 2025+ both siblings are NaN (SST3 best estimate ends at 2024).

A generic NaN-aware sampler would redistribute missing weight
*proportionally* across all finite leaves — which for the COBE branch
in 1850–1869 would route the native sibling's 1/8 share to
HadSST/ERSST/DCENT in proportion to their weights, shrinking the COBE
family to ~22%. To preserve the COBE family's full 1/4 tree-level
share wherever any SST3-derived leaf exists, we add an explicit
*sibling-redirect* rule in `build_ocean_ensemble.py:LEAF_NAN_SIBLING`:

  - When a primary draw lands on `cobe_sst3` (native) and that leaf is
    NaN, redraw from `cobe_sst3_donor` (its sibling) rather than from
    all finite leaves.

The redirect applies whenever the sibling itself has data. In 2025+
both COBE siblings are NaN simultaneously, no sibling-redirect can
fire, and the generic renormalisation routes those draws to the
remaining three method families — collapsing the COBE family's
tree-level contribution to 0 in 2025+ and renormalising the others
to ~1/3 each. This is intentional and matches the Step 4 rationale.

**Per-member trajectory continuity.** The sampler fixes each
member's primary (leaf, member-index) once at the start, so most of
each member's 1850–2025 trajectory is a single product's coherent
record. Per-year re-draws kick in only when the primary is NaN at a
given year, and they are independent draws year-by-year. Two
boundaries are affected:

- **1869→1870.** The ~1/8 of members assigned primary=SST3 native get
  a year-by-year sibling-redirected SST3-donor draw in 1850–1869, then
  switch to their fixed SST3-native trajectory at 1870. Because both
  siblings share the SST3 best-estimate central and are baselined
  identically, the level step at the transition is bounded by the
  difference between (HadSST4 donor spread around SST3 central) and
  (SST3 native perturbation around SST3 central) at year 1870 —
  typically a small fraction of a degree.
- **2024→2025.** The ~1/4 of members assigned primary to either COBE
  sibling get a year-by-year random draw from HadSST/ERSST/DCENT for
  2025. This is a genuine per-member discontinuity: a member that
  followed a coherent SST3 trajectory through 2024 has its 2025 value
  drawn independently from a different method family.

These boundary draws do **not** affect year-by-year ensemble
statistics (mean, percentiles, σ) — those are computed from N=10,000
draws each year independently. They **do** affect per-member trajectory
statistics that span the boundary year: per-member short-window trends
ending in 2025 (e.g. a 2015–2025 trend computed from one member's
trajectory) inherit extra noise from the 2025 re-draw, and per-member
autocorrelation across the boundary is broken for the ~25% of members
affected. The "Safe vs unsafe uses" section flags trend-significance
counts that depend on per-member trajectory coherence; these were
already on the unsafe list and remain so.

The build script prints `share_<year>` columns for 1860 / 1900 / 2020
/ 2024 / 2025 and hard-fails on a family-share invariant check, so
any drift from the intended weighting surfaces at run time.

The result is a 176 × 10,001 CSV (year + 10,000 anomaly columns).
Year-by-year statistics (mean, 2.5/97.5 percentiles, SD) are then
computed for plotting.

Like the LSAT case, we do **not** splice tail and head — all four
leaves cover the full record.

## Known limitations

1. **The SST3-donor sibling's spread is HadSST4-shaped.** Within the
   COBE branch, the SST3 native sibling carries its own 300-member
   perturbation uncertainty (1870–2024) but the SST3 donor sibling
   inherits HadSST4's bias-correction-parameter spread template
   (wrapped around the SST3 best-estimate central). The donor sibling
   contributes 12.5% of the ensemble in 1870–2024 and absorbs the full
   25% of the COBE branch in 1850–1869 via sibling-redirect. Both
   COBE siblings share SST3's best-estimate trajectory; only the
   native sibling adds an SST3-internal uncertainty representation.
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
6. **HadSST4-shape influence varies by era.** 1870–2024: ~37.5% of
   the ensemble (25% HadSST4 leaf + 12.5% SST3-donor sibling spread).
   1850–1869: ~50% (HadSST4 leaf + SST3-donor sibling carrying the
   full COBE branch). 2025+: ~33% (HadSST4 leaf renormalised over the
   three remaining method families; COBE family absent). Down from
   ~75% in v1.
7. **ERSSTv6 boundary year (2025) is an extrapolation, not an
   observation.** Each 2025 member inherits its 2024 deviation from
   the 1000-member median (σ-scaled by HadSST4's year-on-year σ
   ratio). Members maintain trajectory continuity at the boundary,
   but late-year σ growth from non-coverage sources (e.g., a sudden
   bias-correction revision specific to 2025) is not captured. Hard-
   fails if the gap to `END_YEAR` reaches 2 years.
8. **COBE family does not contribute to the 2025 ensemble.** SST3's
   best estimate ends in December 2024; both COBE siblings depend on
   it and both are NaN in 2025+. The other three method families
   (HadSST, ERSST, DCENT) renormalise to ~1/3 each in 2025. The ~25%
   of members assigned primary to either COBE sibling receive an
   independent random draw from one of those three families for 2025
   — see Step 5 for the per-member trajectory implications.
9. **Sea-ice region handling is not harmonised across products**
   (see Step 1 sea-ice table). Two products include ice cells
   (COBE-SST 3 with interpolated near-freezing values; DCENT-I with
   blended 2-m air-temp) and two exclude them (HadSST4, ERSSTv6).
   Cross-product disagreement from sea-ice convention contributes
   ~0.01–0.05 °C to modern-era spread, absorbed by the family tree
   as structural diversity rather than reconciled.

## Safe vs unsafe uses

**Safe:**
- Plotting central estimate + 90% band.
- Reporting the ensemble mean as the SST best estimate.
- Comparing with the sister land LSAT ensemble for context.

**Unsafe:**
- Reporting tail-period 2.5/97.5 percentiles as if calibrated
  frequencies.
- Trend-significance p-values from member-percentile counts.
- Per-leaf attribution without running the sensitivity tests in Step 4.
- **Per-member trajectory statistics that straddle 2024→2025** —
  short-window trends ending in 2025 (e.g. a 2015–2025 trend computed
  from a single member's trajectory) carry inflated noise from the
  ~25% of members whose 2025 value is an independent draw from a
  different method family than their 2024 value (Step 5). Ensemble
  year-by-year statistics (mean, percentile, σ at year 2025) are not
  affected; the inflation is specific to trajectory-coherence-dependent
  metrics.
- Pairing DCENT SST and DCLSAT members by index across sister files.
  DCENT's joint dynamical consistency is destroyed once members are
  sampled independently by the family-tree sampler (most of the
  10,000 final SST members come from non-DCENT leaves anyway), so the
  SST and LSAT outputs do not retain DCENT's joint pairing.

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
- Ishii, M., Nishimura, A., Yasui, S., & Hirahara, S. (2025),
  *Historical High-Resolution Daily SST Analyses with Consistency to
  Monthly Land Surface Air Temperature*, J. Meteor. Soc. Japan 103(1),
  17–44 — for COBE-SST3 methodology (NMAT-anchored bias correction,
  300-member perturbation ensemble); gridded files at
  `https://climate.mri-jma.go.jp/pub/archives/Ishii-et-al_COBE-SST3/`.
- Chan, D. *et al.* (2024), DCENT, Scientific Data 11 — for the
  dynamical-consistency joint SST/LSAT correction.
- Morice, C. P. *et al.* (2021), HadCRUT5, JGR-Atmos — for the σ
  decomposition we reuse.

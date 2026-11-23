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
| 2 | CRUTEM5 v5.2.0.0 | No | uncorrelated σ + correlated σ + coverage σ + bias 2.5/97.5% limits (annual and monthly) | 1857 – 2026-09 | 1961-1990 | `Land Data/CRUTEM5/CRUTEM.5.2.0.0.component_series.global.annual.csv` (+ hemispheric annual series and `CRUTEM.5.2.0.0.anomalies.nc` for coverage masks) |
| 3 | GloSATLAT v1.0.0.0 | No | same components as CRUTEM5 (annual and monthly) | 1781-01 – 2021-12 (annual components from 1857 used) | 1961-1990 | `Land Data/GloSATLAT/GloSATLAT-1-0-0-0_component-series_global_annual.nc` |
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

### 2.2 CRUTEM5: synthesized from the published annual uncertainty components

CRUTEM5 (v5.2.0.0) publishes, for the global annual mean, the best
estimate *a(t)* and four uncertainty components (Osborn et al. 2021):
uncorrelated 1σ *u(t)* (measurement and sampling), correlated 1σ *c(t)*
(station homogenisation and climatological-normal errors, systematic for
a station but independent between stations), coverage 1σ *v(t)*
(incomplete sampling of each hemisphere's land), and bias 2.5/97.5%
limits (urbanisation and thermometer exposure). We build 200 annual
members with `lsat_synthesis.py`:

*member_k(t) = a(t) + u(t)·w_k(t) + c(t)·x_k(t) + v(t)·y_k(t) + U_k(t) + E_k(t)*

All terms are independent and each standardised process has unit
variance, so on CRUTEM5's native 1961–1990 baseline every component
reproduces its published annual σ.

- **Uncorrelated, w:** white noise across years. The published annual σ
  is the monthly σ divided by √12, i.e. independent months.
- **Correlated, x:** AR(1) with a 40-year e-folding time, from the
  homogenisation-break rate in Morice et al. (2012) ("drawn from a
  Poisson distribution with a 40 year repeat rate"). The published annual
  σ equals the monthly σ, i.e. fully correlated within a year. This term
  is small (0.01–0.02 °C after 1900).
- **Coverage, y:** correlation between years
  *corr(L) = (1 − p)·φ^L + p*, with *p* and *φ* estimated by
  `estimate_coverage_persistence.py` (Section 2.2.1).
- **Urbanisation, U:** the one-sided term of Brohan et al. (2006),
  1σ = 0.0055 °C/decade from 1900, pivoting at 1975.5 (centre of the
  1961–1990 normals): *U_k(t) = −max(0, Z_k) · 0.0055 · (max(t, 1900) −
  1975.5)/10*. As in Morice et al. (2012), negative rate draws are set
  to zero, so members lie above the best estimate before 1975 and below
  it after. This reproduces the published bias limits after 1950 to
  1e-7 °C.
- **Exposure, E:** one split-normal draw per member with per-year scales
  *s_lo(t)*, *s_hi(t)* fitted so that U + E reproduces the published
  bias 2.5/97.5% limits in every year (to ≤ 0.003 °C). E is zero after
  1949. The scales are an empirical fit, not a published model; they
  absorb an unexplained mismatch between a simple quadrature model and
  the published limits in 1905–1949.

The one-draw-per-member terms (U, E) use evenly spaced normal quantiles
in random order, so a 200-member leaf carries the intended mean offset.
Because U is one-sided, the CRUTEM5 leaf mean sits slightly below the
best estimate in recent decades (about −0.0022 °C/decade in trend after
1900). This is part of CRUTEM5's published uncertainty and is retained.
A recentred build (`--recentre`) is the sensitivity test.

The synthesis uses its own random stream (`SYNTH_SEED`). Validation with
4,000 members: each component's SD on the native baseline is within 3%
of the published annual σ (Monte Carlo noise), and the member
2.5/97.5% bias quantiles are within 0.003 °C of the published limits.

#### 2.2.1 Year-to-year persistence of the coverage error

The published files give the coverage σ for each year but not how
coverage errors correlate between years. That correlation determines how
much coverage uncertainty survives the per-member 1981–2010 rebaselining
of Step 3. We estimate it as the Met Office estimates the σ itself
(Morice et al. 2012, §5.3; Osborn et al. 2021): CRUTEM5's monthly 5°
coverage mask for year *t* is applied to every year *j* of ERA5 2 m
temperature (reduced to the 5° grid by `prepare_era5_5deg.py`). The
coverage error is the cos(lat)-weighted mean over observed boxes minus
the mean over all boxes containing land, with hemispheres combined as
(2·NH + SH)/3.

- **Check against the published σ:** using ERA5 1979–2018, the estimate
  is 0.93× CRUTEM5's published annual coverage σ (5–95% over years:
  0.86–1.03). The NH matches closely; the SH is about 20–25% low after
  1950. The ERA5 annual/monthly SD ratio is 0.41, the same as the
  published ratio.
- **Lag correlation:** correlation over *j* (ERA5 1940–2025) between
  *e(t, j)* and *e(t+L, j+L)*, averaged over mask years 1857–2010. It is
  0.29 at lag 1 and about 0.1 at lags 2–14. Detrended, it is 0.16 at
  lag 1 and about zero beyond. A bump at lags 15–16 appears with fixed
  masks and in ERA5 1979–2025 alone, so it is ERA5 sampling noise and
  is excluded from the fit (lags 1–10).
- **Fit:** *p = 0.12, φ = 0.19* (`coverage_persistence_fit.json`), used
  in the headline build. Sensitivity builds: white noise, the detrended
  fit (*p = 0, φ = 0.14*), AR(1) with a 4-year e-folding time, and fully
  persistent. Across these the 10,000-member ensemble's spread in
  long-term change (1850–1900 to 2016–2025) is 0.159–0.169 °C (headline
  0.160) and in the 1980–2025 trend 0.020–0.022 °C/decade (headline
  0.021).

### 2.3 GloSATLAT: same synthesis as CRUTEM5

GloSATLAT publishes the same four components in annual files, with the
same structure: the annual/monthly σ ratios are 0.29 (uncorrelated), 1.0
(correlated) and 0.41 (coverage), and the bias limits contain the same
urbanisation ramp. Its 200 members are built exactly as in Section 2.2
from its own components, using the CRUTEM5 coverage-persistence fit
(GloSATLAT's coverage σ is also estimated from ERA5; Morice et al. 2025).
Validation as for CRUTEM5: component SDs within 3% of published, bias
quantiles within 0.002 °C.

CRUTEM5 and GloSATLAT use the same seed on the same 1850–2025 axis, so
member *k* of each carries the same standardised draws. GloSATLAT ends
in 2021; the family tree gives its weight to CRUTEM5 for 2022–2025
(Step 4), and a mixture member drawn from GloSATLAT member *k* continues
with CRUTEM5 member *k*, keeping its error draws. A level step of about
−0.026 °C remains at the switch because GloSATLAT's 2021 best estimate
is that much below CRUTEM5's (1981–2010 baseline); it affects the 10% of
members drawn from GloSATLAT.

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

### 2.5 NOAA Land: best estimate + CRUTEM5 error model

NOAA's v6.1 aravg text has its variance columns set to −999, and the
gridded file has no land-only uncertainty field. NOAA Land's 200 members
are its own annual best estimate (aravg `anomaly` column) plus the
CRUTEM5 error model of Section 2.2: CRUTEM5's published annual σ
components (back-filled to 1850 from CRUTEM5's first year, 1857), the
same correlation structure and the same bias terms, drawn with an
independent seed.

The errors are recentred (the member-mean error is removed each year),
so the leaf mean equals NOAA's best estimate. CRUTEM5's one-sided
urbanisation term describes a product without an urban adjustment;
NOAA's pairwise homogenisation removes part of the urban signal, so
here the term adds spread but not an offset.

NOAA Land contributes no independent uncertainty information; its
uncertainty is CRUTEM5's, applied to NOAA's best estimate.

### 2.6 C-LSAT 2.1: best estimate + CRUTEM5 error model

C-LSAT 2.1 (Sun, Li, et al., Sun Yat-sen University; the v2.1 update
of the Xu et al. 2018 "C-LSAT" / Sun et al. 2021 "C-LSAT 2.0" pipeline;
hosted at `gwpu.net` and figshare
[doi:10.6084/m9.figshare.28255394](https://doi.org/10.6084/m9.figshare.28255394))
contributes a distinct method-family lineage to the LSAT side:
station-level homogenization on a merged station network that
substantially expands on GHCN-monthly's Asian and African coverage
(~25,000 stations across the updated network including non-GHCN
Chinese / Russian / Latin-American collections). The structural
distinction from NOAA PHA is real: C-LSAT uses MASH-style penalized
multiple-reference tests rather than NOAA's pairwise SNHT cascade,
and the input network differs by ~20% pre-1950.

The C-LSAT lineage is *not* new since AR6: AR6 WG1 Chapter 2 included
C-LSAT 2.0 (Sun et al., 2021) as the land component of the China-MST
product (Table 2.3, qualifying as "Land only"; LSAT trend rows
contributing to the AR6 LSAT average in Table 2.4). What is new in
this build is (a) the version update from 2.0 to 2.1, and (b)
treating C-LSAT as an explicit CMA-homogenization method family at
P=1/5 of the LSAT tree, rather than averaging it into an unweighted
multi-product mean alongside HadCRUT5 / NOAA / Berkeley / Kadow as
AR6 did.

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
- Add 200 members of the CRUTEM5 error model (Section 2.5), with its own
  independent seed, recentred on C-LSAT's best estimate.

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
- Half-weight the two leaves without their own uncertainty (NOAA Land and C-LSAT 2.1
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
2. From that leaf's already-rebaselined ensemble (Berkeley Earth 10;
   CRUTEM5, GloSATLAT, DCLSAT, NOAA Land and C-LSAT 2.1 200 each),
   sample one member uniformly at random.
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
2. **The CRU-lineage error model is partly empirical.** The
   urbanisation term follows the published model (Brohan et al. 2006),
   but the exposure scales are fitted per year to the published bias
   limits, and the year-to-year persistence of coverage error comes from
   ERA5, whose own early-period (pre-1979) inhomogeneity can look like
   persistence. ERA5 reproduces the published coverage σ to 0.93× overall
   but about 20–25% low in the Southern Hemisphere after 1950.
3. **DCLSAT marginal-LSAT spread is conditioned on SST/LSAT joint
   consistency** and is plausibly ~10–30% narrower than a standalone
   DCLSAT product's uncertainty would be.
4. **Two of six leaves contribute no independent uncertainty.** NOAA
   Land and C-LSAT 2.1 carry CRUTEM5's error model around their own best
   estimates (Sections 2.5–2.6), so with the CRU-lineage leaves about 60%
   of the final ensemble's within-leaf uncertainty follows CRUTEM5's
   published error model.
5. **Effective N in the tail is far below 10,000.** With Berkeley
   Earth's 10 native trajectories sitting under P=0.25 of the tree,
   the effective independent draws in the tail (1850–1900) is
   plausibly closer to 200–500 once all leaves are accounted for.
   Percentile-based p-values should not be used for trend
   significance — the spread is structural, not stochastic.

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
  126(3), e2019JD032361.
- Osborn, T. J. *et al.* (2021), *Land surface air temperature variations
  across the globe updated to 2019: the CRUTEM5 data set*, JGR-Atmos 126,
  [doi:10.1029/2019JD032352](https://doi.org/10.1029/2019JD032352) — the
  four published uncertainty components and the ERA5-based coverage σ.
- Brohan, P. *et al.* (2006), *Uncertainty estimates in regional and
  global observed temperature changes: a new data set from 1850*, JGR 111,
  D12106, [doi:10.1029/2005JD006548](https://doi.org/10.1029/2005JD006548)
  — one-sided urbanisation uncertainty (0.0055 °C/decade from 1900).
- Morice, C. P. *et al.* (2012), *Quantifying uncertainties in global and
  regional temperature change using an ensemble of observational
  estimates: the HadCRUT4 data set*, JGR 117, D08101,
  [doi:10.1029/2011JD017187](https://doi.org/10.1029/2011JD017187) —
  coverage σ by masked reanalysis (§5.3), the 40-year homogenisation-break
  rate, and zero-truncated urbanisation draws.
- Chan, D., Gebbie, G., Huybers, P. & Kent, E. C. (2024), *DCENT:
  Dynamically Consistent ENsemble of Temperature*, Scientific Data 11,
  article 612 — for the dynamical-consistency LSAT/SST joint correction
  underlying DCLSAT.
- Morice, C. P. *et al.* (2025), *An observational record of global
  gridded near-surface air temperature change over land and ocean from
  1781*, ESSD 17, [doi:10.5194/essd-17-7079-2025](https://doi.org/10.5194/essd-17-7079-2025)
  — GloSAT (GloSATLAT uncertainty components).
- Taylor, M. *et al.* (2025), *GloSAT land air temperature*, GDJ —
  for the LATsdb / LEK normals / exposure-bias adjustment used by
  GloSATLAT.
- Rohde, R. A. *et al.* (2013, 2020 updates), Berkeley Earth methodology
  papers — for the scalpel-based homogenization and kriging.
- Yin, X. *et al.* (2024, 2025), NOAAGlobalTemp v6/v6.1 documentation —
  for the LSAT analysis used in NOAA Land.

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
| 2 | ERSSTv6 (via NOAAGlobalTemp v6.1.0 aravg) | No | deterministic only — v6.1 aravg ocean variance columns are populated as `-999` | 1854-01 – 2026-04 | 1991-2020 (or 1971-2000 in legacy aravg framing) | `Ocean Data/ERSSTv6/aravg.mon.ocean.90S.90N.v6.1.0.202604.asc` |
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

## Step 2 — Build a per-dataset ensemble

We use native counts where available (DCENT 200) and synthesize 200
members for HadSST4 from its decomposed σ. For deterministic
products we use a donor (Section 2.3).

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

### 2.2 DCENT SST — pull 200 native DCENT members, compute SST global means

The DCENT v3.0 Harvard Dataverse release exposes 200 gridded NetCDFs
(~25 MB each, ~5 GB total) with separate `sst`, `lsat`, and
`temperature` fields. We download all 200, compute the
area-weighted (cos(lat)) global mean of the `sst` field per member per
month, store the resulting 200 × n_months matrix as
`DCENT_SST_monthly_ensemble.csv`, and delete the raw NetCDFs (the
fileId-per-member manifest in `dcent_member_fileids.json` makes the
pull reproducible).

`sst` is defined only on ocean grid cells (NaN over land). The global
mean is Σ(cos(lat) · sst) / Σ(cos(lat) · mask), so it is an
ocean-area-weighted mean restricted to ocean cells.

**Caveat (NOT a mirror of land):** each DCENT member's SST field is
*paired* with a specific LSAT field under DCENT's dynamical-consistency
constraint, so the SST and LSAT marginals are not independent samples
of standalone-product uncertainty. **However**, SST is the larger
reservoir and the larger uncertainty source in DCENT — its joint
constraint mainly prunes SST extremes that are inconsistent with the
smaller LSAT degrees of freedom, not the other way around. The
SST-marginal spread is therefore plausibly only ~5–15% narrower than
a standalone SST product would show (vs. the ~10–30% narrowing we
flag for DCLSAT). We apply no correction in v1.

### 2.3 ERSSTv6 and COBE-SST2 — deterministic best estimate + HadSST4-donor uncertainty (with sparse-era σ inflation)

Both ERSSTv6 (via NOAAGlobalTemp aravg) and COBE-SST2 ship only a
deterministic global-mean monthly time series in the form we can read.
For uncertainty we use a donor approach analogous to NOAA Land in the
LSAT methodology:

- Take HadSST4's native 200-member ensemble (Section 2.1),
  demean it (subtract the per-year ensemble mean), and rescale to
  match a target 1σ envelope. This gives a 200-member uncertainty
  ensemble that is anchored on HadSST4's structural uncertainty
  family.
- Add the rescaled donor to ERSSTv6's (or COBE-SST2's) best estimate
  on the common annual axis.

**Choice of donor and rationale:** HadSST4 is the only SST product in
this catalogue with a published structural uncertainty representation.
It is the canonical in-situ SST reconstruction. Using HadSST4 as
donor for both ERSSTv6 and COBE-SST2 is methodologically natural for
ERSSTv6 (which is also in-situ-based and shares much of the same
underlying data) and acceptable for COBE-SST2 (which adds satellite
data in the modern era but uses the same in-situ archive in the tail).
COBE-SST2's analysis-error field exists on the gridded NetCDF
distributed by JMA/TCC but is **not publicly distributed in the NOAA
PSL mirror** we use; we therefore cannot use it as a per-product
target σ and fall back to HadSST4. Documented limitation: ERSSTv6 and
COBE-SST2 contribute only their best-estimate signal, not independent
uncertainty information.

**Target σ for rescaling — sparse-era inflated.** The base target σ is
HadSST4's own annual σ. To partially compensate for the (uncaptured)
extra structural uncertainty that ERSST and COBE bring from EOF /
optimal-interpolation reconstruction in the sparse pre-satellite era,
we **inflate the target σ by a year-dependent factor**:
- pre-1920: 1.3× HadSST4 σ
- 1920–1960: linearly decays from 1.3× to 1.0×
- post-1960: 1.0× HadSST4 σ

This is a deliberate compensation, not a measurement. The
sensitivity test version with inflation factor = 1.0× throughout
(identity rescaling) is reported alongside.

## Step 3 — Re-baseline all four 200-member ensembles to 1981–2010

For every member of every dataset:
- subtract that member's own 1981–2010 mean.

After this step all four ensembles share zero mean over 1981–2010 and
are directly comparable.

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

**Rationale and a forward-looking caveat.** A strict "lineage of
uncertainty" reading of Thorne 2026 §3.2.5 would *down*-weight ERSSTv6
and COBE-SST2 because their envelopes in v1 are donor-imputed from
HadSST4 (NCEI does not currently publish a per-member ERSSTv6
ensemble, and the COBE-SST2 analysis-error field on the JMA archive is
not openly downloadable — see "Open issue" below). Half-weighting the
two donor-imputed leaves was the stats-review recommendation under
that strict reading.

We instead keep equal 1/4 weights deliberately, in anticipation that
**ERSSTv6 will acquire its own native ensemble representation** in a
future revision (likely via the NOAAGlobalTempv5/v6 500-member
ensemble that Thorne et al. 2026 use as donor for ERSST-based products,
or a future NCEI release of an ERSSTv6 ensemble proper). Once that
ensemble is wired in, ERSSTv6 will contribute its own uncertainty
template and the equal 1/4 weighting becomes the correct "structural
lineage" tree without any further re-weighting.

The transient implication is that ~50% of v1's ensemble has
HadSST4-shaped uncertainty wrapped around three different best
estimates (HadSST4's own, ERSSTv6's, COBE-SST2's). This is flagged in
the limitations section and the safe/unsafe-uses block.

**Open issue (to revisit when ensembles become available):**
- ERSSTv6 ensemble: not currently in NCEI's public `v6/access/` tree.
  When obtained (from NCEI directly, or via Boyin Huang), switch
  ERSSTv6 from "deterministic + HadSST4 donor" to its native ensemble.
- COBE-SST2 analysis-error: the JMA/TCC archive exposes an analysis
  error field but only via gated forms; if extractable, switch
  COBE-SST2's target σ from "HadSST4 σ" to "COBE-derived σ".
- NOAAGlobalTempv5 500-member ensemble (Thorne donor) — available
  locally in the parent `GMST ensemble/ManagedData/Data/NOAA_ensemble/`
  but it is the *merged* GMST product, not SST-only. Using it as the
  ERSSTv6 donor would change the uncertainty *template family* away
  from HadSST4 toward NOAA — an improvement in structural-lineage
  fidelity but with the asterisk that the ensemble is merged GMST not
  SST-only and ends in 2016.

**Note on asymmetry with the LSAT tree.** Both the LSAT tree and the
SST tree use four "method families" (HOMOG/PHA/CRU-lineage/Dynamical
on the LSAT side; HadSST/ERSST/COBE/Dynamical on the SST side). On
LSAT the CRU-lineage family splits two ways (CRUTEM5 + GloSATLAT) so
its leaves carry P=1/8 each; all other LSAT and SST leaves carry
P=1/4. The result is symmetric for DCLSAT vs DCENT SST (both P=1/4)
and for the family-level weights. This symmetry should be preserved
when SST × LSAT are recombined downstream into a GMST ensemble.

**Planned sensitivity tests (out of scope for v1):**
- **Half-weight ERSSTv6 and COBE-SST2** (the alternative tree:
  HadSST4 1/3, DCENT 1/3, ERSST 1/6, COBE 1/6). This is the
  defensible interpretation under a strict "lineage of uncertainty"
  reading of Thorne. Report as a sensitivity for transparency.
- Single in-situ-only sub-family containing HadSST4, ERSSTv6,
  COBE-SST2 at P=1/2 (split 2/4, 1/4, 1/4 within), DCENT at P=1/2 —
  the most aggressive version, lumping all three in-situ-based
  products as one lineage.
- Sparse-era σ inflation factor = 1.0× (no inflation, identity
  rescaling) — to quantify how much the inflation contributes to tail
  spread.
- Switch ERSSTv6 donor from HadSST4 to the NOAAGlobalTempv5
  500-member ensemble — to quantify the donor-choice sensitivity.
- Splice at the 1995/96 midpoint of the 1981–2010 reference period
  (matching Thorne's GMST approach).

## Step 5 — Generate 10,000-member ensemble

For each of 10,000 draws:
1. Pick a leaf dataset using the leaf weights from Step 4.
2. From that leaf's (200-member, already-rebaselined) ensemble, sample
   one member uniformly at random.
3. Take that member's full 1850–2025 annual time series.
4. Record it as one ensemble member.

The result is a 176 × 10,001 CSV (year + 10,000 anomaly columns).
Year-by-year statistics (mean, 2.5/97.5 percentiles, SD) are then
computed for plotting.

Like the LSAT case, we do **not** splice tail and head — all four
leaves cover the full record.

## Known limitations

1. **Two of four products contribute no independent uncertainty in v1.**
   ERSSTv6 and COBE-SST2 are donor-imputed; only their best estimates
   are independent. With equal 1/4 leaf weights, ~50% of the final
   ensemble's uncertainty originates from HadSST4-shaped donor spread
   wrapped around three different best-estimate trajectories. This is a
   *transitional* arrangement (see Step 4 rationale) that resolves once
   ERSSTv6 acquires its own ensemble representation; users should not
   read the v1 envelope as four-way independent structural diversity.
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
6. **HadSST4 as the universal donor** propagates HadSST4's bias and
   coverage uncertainty structure into ERSSTv6 and COBE-SST2 envelopes
   — any HadSST4-specific shape (e.g. its bias-CI assumption about
   bucket vs. engine-room mixtures) bleeds into 50% of the ensemble.

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
- SST uses **HadSST4** as donor for its imputed leaves (ERSSTv6, COBE-SST2).
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
  for the ERSSTv6 methodology.
- Hirahara, S., Ishii, M., & Fukuda, Y. (2014), *Centennial-scale SST
  analyses*, J. Climate — for COBE-SST2 methodology.
- Chan, D. *et al.* (2024), DCENT, Scientific Data 11 — for the
  dynamical-consistency joint SST/LSAT correction.
- Morice, C. P. *et al.* (2021), HadCRUT5, JGR-Atmos — for the σ
  decomposition we reuse.

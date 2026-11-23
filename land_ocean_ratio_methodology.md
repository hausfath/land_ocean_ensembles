# Land/ocean warming ratio — methodology

Companion to `compute_land_ocean_ratio.py`. Documents inputs, computation,
the statistical treatment of the ratio, caveats, and reproducibility, for the
AR7 WG1 Chapter 2 land/ocean warming-contrast paragraph.

## Inputs

- `land_ensemble.csv` — 10,000-member annual ensemble of global-mean **land**
  surface air temperature anomalies, 1850–2025, built via the hierarchical
  family-tree sampling described in `land_ensemble_methodology.md`
  (Berkeley Earth, CRUTEM5, NOAAGlobalTemp-land, C-LSAT, GloSATLAT, DCLSAT).
- `ocean_ensemble.csv` — 10,000-member annual ensemble of global-mean **ocean**
  (SST) anomalies, 1850–2025, built analogously (HadSST4, ERSSTv6, COBE-SST3,
  DCENT-SST, with COBE-SST2 in the family tree); see
  `ocean_ensemble_methodology.md`.

Each member is stored as an anomaly on a 1981–2010 baseline; the stored
baseline cancels in the warming difference (below), so it does not affect
results.

## Computation

1. **Per-member warming**: for each member, warming = mean over the recent
   period − mean over 1850–1900. Computed separately within each ensemble,
   giving 10,000 land warming values and 10,000 ocean warming values per
   period. Recent periods: 2016–2025 (primary) and 2011–2020 (AR6 period,
   for like-for-like comparison).
2. **Marginal summaries**: median and 2.5/97.5 percentiles across members
   (e.g. ocean 0.90 [0.68–1.11] °C for 2016–2025; the chapter's assessed
   ranges use 5–95%, 0.90 [0.70–1.09] °C, see SOD_TODO C3).
3. **Ratio and difference**: the land and ocean ensembles are built
   independently from different source datasets, so members are **not**
   physically paired. The ratio (and land−ocean difference) distribution is
   therefore obtained by Monte Carlo: draw a random land member and a random
   ocean member independently (with replacement; 2,000,000 draws, seed 42) and
   form land/ocean and land−ocean. Report median + 2.5/97.5 percentiles.
4. **Central ratio**: quoted as the **ratio of medians** (land median / ocean
   median), which is the cleaner, less skew-sensitive central estimate and
   matches how AR6's implied ratio is formed (1.59/0.88 = 1.81). The
   median-of-ratios (2.05) is reported alongside but not used as the headline.

## Results (1850–1900 baseline)

| Period | Land (°C) | Ocean (°C) | Ratio (of medians) | Ratio (MC, 2.5–97.5) | Land−ocean (°C) |
|---|---|---|---|---|---|
| 2016–2025 | 1.91 [1.58–2.21] | 0.90 [0.68–1.11] | **2.11** | 2.05 [1.55–2.93] | 0.97 [0.58–1.39] |
| 2011–2020 | 1.64 [1.34–1.95] | 0.78 [0.57–0.98] | **2.10** | 2.06 [1.51–3.05] | 0.84 [0.47–1.24] |

AR6 (2011–2020): land 1.59 [1.34–1.83], ocean 0.88 [0.68–1.01] → implied
ratio 1.81, implied difference 0.71 °C. Both our ratio (2.10) and difference
(0.84) exceed AR6 for the same period, driven principally by the lower
assessed ocean warming (0.78 vs 0.88, from the updated SST products with
warmer early-twentieth-century SSTs) and secondarily by slightly higher land
warming (1.64 vs 1.59).

## Statistical caveats (vetted by independent stats review)

- **Independence is conservative.** Independent pairing treats land and ocean
  warming errors as uncorrelated. They share structural components — notably a
  common DCENT ancestor (DCLSAT `lsat` / DCENT-I `sst`), the 1850–1900
  reference period, and early-record coverage choices — which induce positive
  correlation. For a land/ocean contrast (ratio or difference), positive
  correlation makes common-mode error partly cancel, so the **true intervals
  are likely narrower** than quoted. The quoted ranges are therefore a
  conservative (wide) bound, not a tight best estimate.
- **Ratio skew.** The ratio is right-skewed (upper tail ~3.0) because it can
  pair a high land member with a low ocean member — exactly the anti-correlated
  combination that shared structural error makes unlikely. The land−ocean
  **difference** is not subject to this multiplicative tail inflation and is
  the more robust uncertainty statement; it is reported alongside and
  corroborates the "larger than AR6" finding.
- **Denominator well-bounded.** Ocean warming (2.5th pct ~0.68 °C) is far from
  zero, so the ratio has no blow-up pathology.

## Reproducibility

`python3 compute_land_ocean_ratio.py` (seed 42) regenerates the table above and
writes `land_ocean_ratio_summary.csv`. Requires the two ensemble CSVs in the
same directory.

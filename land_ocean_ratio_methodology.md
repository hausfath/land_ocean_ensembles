# Land/ocean warming ratio — methodology

Companion to `compute_land_ocean_ratio.py`. Documents inputs, computation,
the statistical treatment of the ratio, caveats, and reproducibility, for the
AR7 WG1 Chapter 2 land/ocean warming-contrast paragraph.

## Inputs

- `land_ensemble.csv` — 10,000-member annual ensemble of global-mean **land**
  surface air temperature anomalies, 1850–2025, built via the hierarchical
  family-tree sampling described in `land_ensemble_methodology.md`
  (Berkeley Earth, CRUTEM5, NOAAGlobalTemp-land, C-LSAT, GloSATLAT, DCLSAT;
  ERA5 and JRA-3Q replace 1/6 of members from 1996).
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
2. **Marginal summaries**: median and 5–95% range across members (e.g.
   ocean 0.90 [0.70–1.09] °C for 2016–2025).
3. **Ratio and difference**: the land and ocean ensembles are built
   independently from different source datasets, so members are **not**
   physically paired. The ratio (and land−ocean difference) distribution is
   therefore obtained by Monte Carlo: draw a random land member and a random
   ocean member independently (with replacement; 2,000,000 draws, seed 42) and
   form land/ocean and land−ocean. Report median and 5–95% range.
4. **Trends**: the same pairing applied to per-member OLS linear trends
   (°C/decade) over the periods in the second table.
5. **Central ratio**: quoted as the **ratio of medians** (land median / ocean
   median), which is less sensitive to the skew of the ratio distribution
   and matches how AR6's implied ratio is formed (1.59/0.88 = 1.81). The
   median of the paired ratios is reported alongside.

## Results

All ranges are 5–95%.

**Warming since 1850–1900 (°C)**

| Period | Land | Ocean | Ratio of medians | Paired ratio | Land − ocean |
|---|---|---|---|---|---|
| 2016–2025 | 1.91 [1.61–2.15] | 0.90 [0.70–1.09] | **2.11** | 2.06 [1.62–2.82] | 0.98 [0.65–1.32] |
| 2011–2020 | 1.64 [1.37–1.89] | 0.78 [0.59–0.96] | **2.11** | 2.06 [1.58–2.91] | 0.84 [0.53–1.17] |

AR6 (2011–2020): land 1.59 [1.34–1.83], ocean 0.88 [0.68–1.01], implied
ratio 1.81, implied difference 0.71 °C. Both our ratio (2.11) and difference
(0.84) exceed AR6 for the same period, mainly because of the lower ocean
warming (0.78 vs 0.88, from the updated SST products with warmer
early-twentieth-century SSTs) and secondarily the slightly higher land warming
(1.64 vs 1.59).

**Linear trends (°C/decade)**

| Period | Land | Ocean | Ratio of medians | Paired ratio | Land − ocean |
|---|---|---|---|---|---|
| 1850–2025 | 0.094 [0.076–0.111] | 0.046 [0.033–0.059] | 2.03 | 1.98 [1.46–2.95] | 0.047 [0.025–0.069] |
| 1901–2025 | 0.131 [0.111–0.146] | 0.078 [0.072–0.094] | 1.68 | 1.65 [1.32–1.92] | 0.051 [0.028–0.069] |
| 1950–2025 | 0.230 [0.212–0.250] | 0.118 [0.104–0.136] | 1.95 | 1.93 [1.62–2.30] | 0.111 [0.083–0.139] |
| 1960–2025 | 0.269 [0.246–0.294] | 0.137 [0.116–0.156] | 1.97 | 1.97 [1.67–2.38] | 0.133 [0.101–0.166] |
| 1979–2025 | 0.319 [0.293–0.353] | 0.144 [0.123–0.161] | 2.21 | 2.23 [1.91–2.68] | 0.177 [0.143–0.215] |
| 1981–2025 | 0.327 [0.299–0.361] | 0.152 [0.129–0.167] | 2.16 | 2.18 [1.88–2.61] | 0.178 [0.143–0.217] |
| 1991–2025 | 0.361 [0.325–0.411] | 0.171 [0.151–0.181] | 2.11 | 2.13 [1.87–2.51] | 0.192 [0.153–0.243] |

The 1901–2025 ratio is lower because these SST products warm strongly in the
early twentieth century relative to land.

## Dependence between land and ocean errors

Independent pairing assumes land and ocean warming errors are uncorrelated.
Two kinds of dependence exist, and they push in opposite directions.

- **Shared error components** (a common DCENT ancestor for DCLSAT and DCENT-SST,
  similar early-coverage gaps) would correlate land and ocean errors
  positively. For a ratio or difference, positive correlation partly cancels,
  so the true intervals would be narrower than quoted.
- **Structural pairing across producers** goes the other way. Pairing each
  producer's own land and ocean series (per-dataset ensemble means, 2016–2025
  vs 1850–1900):

  | Producer | Land | Ocean | Ratio | Land − ocean | 1981–2025 trend ratio | 1950–2025 trend ratio |
  |---|---|---|---|---|---|---|
  | HadCRUT5 (CRUTEM5 / HadSST4) | 1.75 | 1.01 | 1.74 | 0.74 | 1.87 | 1.64 |
  | Berkeley Earth (BE land / HadSST4, modified) | 1.73 | 1.01 | 1.72 | 0.73 | 1.92 | 1.72 |
  | DCENT (DCLSAT / DCENT-SST) | 1.94 | 1.04 | 1.87 | 0.90 | 2.05 | 1.73 |
  | NOAAGlobalTemp (GHCN land / ERSSTv6) | 2.01 | 0.73 | 2.73 | 1.27 | 2.65 | 2.30 |
  | CMST3 (C-LSAT / ERSSTv6) | 2.01 | 0.73 | 2.73 | 1.27 | 2.54 | 2.24 |

  The producers with the most land warming use the SST product with the
  least warming (ERSSTv6), so across producers land and ocean warming are
  anti-correlated. That would make the independent interval too narrow. With
  five producers sharing two SST products this is weak evidence, and GloSATLAT
  and COBE-SST3 have no partner.

Rank-correlation (Gaussian copula) sensitivity, 2016–2025 (the 1950–2025 and
1981–2025 trend versions are in `land_ocean_ratio_sensitivity.csv`):

| Assumed correlation | Paired ratio | Land − ocean (°C) |
|---|---|---|
| 0 (headline) | 2.06 [1.62–2.82] | 0.98 [0.65–1.32] |
| 0.3 | 2.06 [1.68–2.72] | 0.98 [0.70–1.26] |
| 0.5 | 2.06 [1.73–2.64] | 0.98 [0.74–1.23] |
| 0.7 | 2.06 [1.79–2.55] | 0.98 [0.79–1.18] |

The central values do not depend on the assumed correlation; only the width
does. Independence is therefore kept as the headline. It is neither clearly
conservative nor clearly optimistic, and the producer-paired ratios
(1.7–2.7) are a check on the structural range.

**Other caveats**

- **Ratio skew.** The paired ratio is right-skewed (upper tail ~2.8–2.9)
  because it can pair a warm land member with a cool ocean member. The
  land−ocean difference has no multiplicative tail inflation and is the more
  robust uncertainty statement.
- **Denominator well-bounded.** Ocean warming (5th percentile 0.70 °C) is far
  from zero, so the ratio has no blow-up pathology.

## Reproducibility

`python3 compute_land_ocean_ratio.py` (seed 42) regenerates all tables above and
writes `land_ocean_ratio_summary.csv` (headline) and
`land_ocean_ratio_sensitivity.csv` (dependence checks). Requires the two ensemble CSVs in the
same directory.

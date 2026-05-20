#!/usr/bin/env bash
# Download the raw input datasets needed to reproduce the LSAT and SST
# ensembles. Total transfer ≈ 2.5 GB. DCENT (~10 GB total during DCLSAT +
# DCENT-SST extraction) is NOT downloaded here — it's pulled member-by-
# member by `pull_dcent_lsat_members.py` and `pull_dcent_sst_members.py`,
# which delete each raw NetCDF after global-mean extraction.
#
# Run from the repo root. The script is idempotent — existing files are
# kept (`curl -L --fail -sS --create-dirs -o ... -z` skips if newer on
# server). Re-run after a few months to pick up updated late-month aravg
# / HadSST4 / CRUTEM5 / GloSATLAT releases.
#
# Note: dataset names below match what the build scripts expect in
# `Land Data/` and `Ocean Data/`. Do not rename folders or files.

set -euo pipefail

curl_get() {
    local url="$1"
    local dest="$2"
    mkdir -p "$(dirname "$dest")"
    echo "  $dest"
    curl -L --fail -sS -o "$dest" "$url"
}

echo "==== LAND ===="

echo "[1/5] Berkeley Earth High-Resolution Land — TAVG ensemble"
curl_get \
    "https://storage.googleapis.com/berkeley-earth-temperature-hr/global/Land_TAVG_ensemble.txt" \
    "Land Data/Berkeley Earth Highres/Land_TAVG_ensemble.txt"

echo "[2/5] CRUTEM 5.1.0.0 — summary + component series, monthly global"
curl_get \
    "https://www.metoffice.gov.uk/hadobs/crutem5/data/CRUTEM.5.1.0.0/diagnostics/CRUTEM.5.1.0.0.summary_series.global.monthly.csv" \
    "Land Data/CRUTEM5/CRUTEM.5.1.0.0.summary_series.global.monthly.csv"
curl_get \
    "https://www.metoffice.gov.uk/hadobs/crutem5/data/CRUTEM.5.1.0.0/diagnostics/CRUTEM.5.1.0.0.component_series.global.monthly.csv" \
    "Land Data/CRUTEM5/CRUTEM.5.1.0.0.component_series.global.monthly.csv"

echo "[3/5] GloSATLAT v1.0.0.0 — summary + component series, monthly global"
curl_get \
    "https://dap.ceda.ac.uk/badc/deposited2025/GloSAT/GloSATLAT-1-0-0-0/diagnostics/summary-series/GloSATLAT-1-0-0-0_summary-series_global_monthly.nc" \
    "Land Data/GloSATLAT/GloSATLAT-1-0-0-0_summary-series_global_monthly.nc"
curl_get \
    "https://dap.ceda.ac.uk/badc/deposited2025/GloSAT/GloSATLAT-1-0-0-0/diagnostics/component-series/GloSATLAT-1-0-0-0_component-series_global_monthly.nc" \
    "Land Data/GloSATLAT/GloSATLAT-1-0-0-0_component-series_global_monthly.nc"

echo "[4/6] DCENT-I — diagnostics + members (pulled iteratively by pull_dcent_i_members.py)"
# DCENT-I v1.1.0.0 (Chan et al. 2026, GDJ; doi:10.7910/DVN/ROG38Q,
# dataset version 2.0 published 2026-03-30 with end-of-2025 extension) is the
# spatially complete kriging-infilled extension of DCENT. The pull script
# downloads each of the 200 member NetCDFs in sequence (~22 MB each, ~4.3 GB
# total), extracts both LSAT and SST global means via static land-fraction
# decomposition (see methodologies §2.4 land / §2.2 ocean), and deletes the
# raw NetCDF before pulling the next. Peak transient disk: ~22 MB.
# The diagnostics.nc file is downloaded once for the land/sea weight field.
echo "  (DCENT-I members + diagnostics pulled iteratively by pull_dcent_i_members.py)"
echo "  v3 archive of the old DCENT v3.0 derived outputs lives under"
echo "  Land Data/DCLSAT/v3_archive/ and Ocean Data/DCENT SST/v3_archive/"

echo "[5/6] NOAA Land v6.1.0 — aravg monthly + annual, land 90S-90N"
echo "  (NOAA Land aravg is used as the central; uncertainty is donor-imputed from DCENT-I LSAT)"
curl_get \
    "https://www.ncei.noaa.gov/data/noaa-global-surface-temperature/v6.1/access/timeseries/aravg.mon.land.90S.90N.v6.1.0.202604.asc" \
    "Land Data/NOAA Land/aravg.mon.land.90S.90N.v6.1.0.202604.asc"
curl_get \
    "https://www.ncei.noaa.gov/data/noaa-global-surface-temperature/v6.1/access/timeseries/aravg.ann.land.90S.90N.v6.1.0.202604.asc" \
    "Land Data/NOAA Land/aravg.ann.land.90S.90N.v6.1.0.202604.asc"

echo "[6/6] C-LSAT 2.1 (Sun et al. / CMA-homogenization)"
# The 5°×5° gridded NetCDF is hosted via the Sun Yat-sen group's portal at
# gwpu.net (linked from http://www.gwpu.net/en/h-col-103.html) and on figshare
# (10.6084/m9.figshare.28255394). The gwpu.net portal has been the timelier of
# the two — the figshare deposit was frozen at the ESSD 2025 publication with
# 1901-2023 coverage, but a separate extended NetCDF covering 1850-01 through
# 2025-12 has been distributed via gwpu.net. We pull the extended version.
# The portal currently routes downloads through a JS-rendered page that
# `curl` cannot follow directly — fetch in a browser and place at:
#     Land Data/C-LSAT/China-LSAT2.1_tavg.nc
echo "  (China-LSAT2.1_tavg.nc — manual download from gwpu.net required)"
if [ ! -f "Land Data/C-LSAT/China-LSAT2.1_tavg.nc" ]; then
    echo "  WARNING: Land Data/C-LSAT/China-LSAT2.1_tavg.nc not found."
    echo "    Download manually from http://www.gwpu.net/en/h-col-103.html"
    echo "    and place at that path before running prepare_c_lsat_global_mean.py"
fi

echo
echo "==== OCEAN ===="

echo "[1/4] HadSST 4.2.0.0 — global mean CSV + 200-member gridded ensemble (1.7 GB)"
curl_get \
    "https://www.metoffice.gov.uk/hadobs/hadsst4/data/data/HadSST.4.2.0.0_monthly_GLOBE.csv" \
    "Ocean Data/HadSST4/HadSST.4.2.0.0_monthly_GLOBE.csv"
for chunk in 1-50 51-100 101-150 151-200; do
    curl_get \
        "https://www.metoffice.gov.uk/hadobs/hadsst4/data/data/HadSST.4.2.0.0_ensemble_members_${chunk}.zip" \
        "Ocean Data/HadSST4/HadSST.4.2.0.0_ensemble_members_${chunk}.zip"
done

echo "[2/4] ERSSTv6 — aravg central + 1000-member native ensemble"
curl_get \
    "https://www.ncei.noaa.gov/data/noaa-global-surface-temperature/v6.1/access/timeseries/aravg.mon.ocean.90S.90N.v6.1.0.202604.asc" \
    "Ocean Data/ERSSTv6/aravg.mon.ocean.90S.90N.v6.1.0.202604.asc"
# The 1000-member native ensemble (~135 GB total if pulled in full) is not
# downloaded here. pull_ersstv6_members.py iterates one member at a time,
# computing the global mean from each member's 2°×2° gridded Fortran-binary
# file (~135 MB) and deleting the raw file before pulling the next. Peak
# transient disk use is ~135 MB; final CSV is ~25 MB. Wall time ~2.7 hours.
echo "  (ERSSTv6 native ensemble pulled iteratively by pull_ersstv6_members.py)"

echo "[3/4] COBE-SST 2 — gridded monthly + 1991-2020 climatology"
curl_get \
    "https://downloads.psl.noaa.gov/Datasets/COBE2/sst.mon.mean.nc" \
    "Ocean Data/COBE-SST2/sst.mon.mean.nc"
curl_get \
    "https://downloads.psl.noaa.gov/Datasets/COBE2/sst.mon.ltm.1991-2020.nc" \
    "Ocean Data/COBE-SST2/sst.mon.ltm.1991-2020.nc"

echo "[4/4] DCENT-I SST — same iterative pull as land (pull_dcent_i_members.py covers both)"
echo "  DCENT-I derives both LSAT and SST from the same member files via"
echo "  land-fraction / sea-fraction decomposition. See [4/6] in the LAND section."

echo
echo "Done. Next step:"
echo "  python3 pull_dcent_i_members.py       # ~22 MB transient peak, ~13 min — produces LSAT + SST CSVs"
echo "  python3 pull_ersstv6_members.py       # ~135 MB transient, deleted after extraction, ~2.7 hr"
echo "  python3 prepare_hadsst_global_ensemble.py"
echo "  python3 prepare_cobe_global_mean.py"
echo "  python3 prepare_c_lsat_global_mean.py"
echo "  python3 build_land_ensemble.py"
echo "  python3 build_ocean_ensemble.py"

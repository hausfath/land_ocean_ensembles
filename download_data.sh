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

echo "[4/5] DCLSAT — ensemble mean and climatology (members pulled by pull_dcent_lsat_members.py)"
curl_get \
    "https://dataverse.harvard.edu/api/access/datafile/13636717" \
    "Land Data/DCLSAT/DCENT_ensemble_mean_1850_2025.nc"
curl_get \
    "https://dataverse.harvard.edu/api/access/datafile/10393657" \
    "Land Data/DCLSAT/DCENT_monthly_climatology_1982_2014.nc"

echo "[5/5] NOAA Land v6.1.0 — aravg monthly + annual, land 90S-90N"
curl_get \
    "https://www.ncei.noaa.gov/data/noaa-global-surface-temperature/v6.1/access/timeseries/aravg.mon.land.90S.90N.v6.1.0.202604.asc" \
    "Land Data/NOAA Land/aravg.mon.land.90S.90N.v6.1.0.202604.asc"
curl_get \
    "https://www.ncei.noaa.gov/data/noaa-global-surface-temperature/v6.1/access/timeseries/aravg.ann.land.90S.90N.v6.1.0.202604.asc" \
    "Land Data/NOAA Land/aravg.ann.land.90S.90N.v6.1.0.202604.asc"

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

echo "[4/4] DCENT SST — ensemble mean (members pulled by pull_dcent_sst_members.py)"
curl_get \
    "https://dataverse.harvard.edu/api/access/datafile/13636717" \
    "Ocean Data/DCENT SST/DCENT_ensemble_mean_1850_2025.nc"

echo
echo "Done. Next step:"
echo "  python3 pull_dcent_lsat_members.py    # ~5 GB transient, deleted after extraction"
echo "  python3 pull_dcent_sst_members.py     # ~5 GB transient, deleted after extraction"
echo "  python3 pull_ersstv6_members.py       # ~135 MB transient, deleted after extraction, ~2.7 hr"
echo "  python3 prepare_hadsst_global_ensemble.py"
echo "  python3 prepare_cobe_global_mean.py"
echo "  python3 build_land_ensemble.py"
echo "  python3 build_ocean_ensemble.py"

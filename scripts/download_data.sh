#!/usr/bin/env bash
# Download the raw data used in this project into data/raw/.
# FOPH/BAG Infectious Diseases Dashboard (IDD) export API and MeteoSwiss open data (NBCN).
set -euo pipefail
mkdir -p data/raw
IDD="https://api.idd.bag.admin.ch/api/v1/export/latest"
for f in CAMPYLOBACTERIOSIS_oblig SALMONELLOSIS_oblig; do
  curl -fsSL "$IDD/$f/csv"      -o "data/raw/BAG_${f}_latest.csv"
  curl -fsSL "$IDD/$f/metadata" -o "data/raw/BAG_${f}_latest.json"
done
NBCN="https://data.geo.admin.ch/ch.meteoschweiz.ogd-nbcn"
for s in sma bas ber gve lug luz stg; do
  curl -fsSL "$NBCN/$s/ogd-nbcn_${s}_m.csv" -o "data/raw/ogd-nbcn_${s}_m.csv"
done
curl -fsSL "$NBCN/ogd-nbcn_meta_parameters.csv" -o "data/raw/ogd-nbcn_meta_parameters.csv"
echo "done; the BAG files are replaced monthly, so results can differ slightly from the archived copies"

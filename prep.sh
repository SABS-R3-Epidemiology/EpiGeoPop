#!/bin/bash
mkdir data/raw
cd data/raw
echo "Downloading global population density map..."
curl -O https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_GLOBE_R2023A/GHS_POP_E2025_GLOBE_R2023A_4326_30ss/V1-0/GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0.zip
unzip GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0.zip
echo "Downloading border files..."
curl -LO "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_0_countries_lakes.zip"
curl -LO "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_admin_1_states_provinces.zip"
curl -LO "https://naciscdn.org/naturalearth/10m/cultural/ne_10m_urban_areas_landscan.zip"
echo "Downloading population age file..."
# This server uses an outdated SSL protocol so we need to enable legacy renegotiation
OPENSSL_CONF=../../openssl.cnf curl -O "https://population.un.org/wpp/assets/Excel%20Files/1_Indicator%20(Standard)/CSV_FILES/WPP2024_PopulationByAge5GroupSex_Medium.csv.gz"
gunzip WPP2024_PopulationByAge5GroupSex_Medium.csv.gz

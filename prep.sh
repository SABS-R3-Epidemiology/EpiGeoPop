#!/bin/bash
mkdir data/raw
cd data/raw
echo "Downloading global population density map..."
curl -L "https://figshare.com/ndownloader/files/53996249" > "GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0.zip"
unzip GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0.zip
echo "Downloading border files..."
curl -L "https://figshare.com/ndownloader/files/53990156" > "ne_10m_admin_1_states_provinces.zip"
curl -L "https://figshare.com/ndownloader/files/53990159" > "ne_10m_urban_areas_landscan.zip"
curl -L "https://figshare.com/ndownloader/files/53990153" > "ne_10m_admin_0_countries_lakes.zip"
echo "Downloading population age file..."
curl -L "https://figshare.com/ndownloader/files/53996234" > "WPP2024_PopulationByAge5GroupSex_Medium.csv.gz"
gunzip WPP2024_PopulationByAge5GroupSex_Medium.csv.gz

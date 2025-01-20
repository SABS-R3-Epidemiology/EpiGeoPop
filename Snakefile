rule all:
    input:
        # "data/processed/countries/Bermuda_microcells.csv",
        # "data/processed/countries/Bermuda_pop_dist.json",
        # "data/processed/countries/Gibraltar_microcells.csv",
        # "data/processed/countries/Gibraltar_pop_dist.json",
        "data/processed/countries/Andorra_microcells.csv",
        # "data/processed/countries/Luxembourg_pop_dist.json",
        # "data/processed/countries/Netherlands_microcells.csv",
        # "data/processed/countries/Netherlands_pop_dist.json",
        # "data/processed/countries/New Zealand_microcells.csv",
        # "data/processed/countries/New Zealand_pop_dist.json",
        # "data/processed/provinces/CA-NB_microcells.csv",
        # "data/processed/provinces/CA-NB_pop_dist.json",
        # "data/processed/cities/Winnipeg_microcells.csv",
        # "data/processed/cities/Winnipeg_pop_dist.json",
        # "data/processed/cities/Dublin2_microcells.csv",
        # "data/processed/cities/Dublin2_pop_dist.json",
        # "outputs/dag.pdf"

# rule render_dag:
#     input:
#         "Snakefile"
#     output:
#         "outputs/dag.pdf"
#     shell:
#         "snakemake --dag | dot -Tpdf > outputs/dag.pdf"

rule get_country_arrays:
    input:
        "data/raw/GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0.tif",
        "data/raw/ne_10m_admin_0_countries_lakes.zip"
    output:
        "data/processed/countries/{country}.npz"
    script:
        "scripts/get_country_arrays.py"

rule get_province_arrays:
    input:
        "data/raw/GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0.tif",
        "data/raw/ne_10m_admin_1_states_provinces.zip"
    output:
        "data/processed/provinces/{province}.npz",
    script:
        "scripts/get_province_arrays.py"

rule get_city_arrays:
    input:
        "data/raw/GHS_POP_E2025_GLOBE_R2023A_4326_30ss_V1_0.tif",
        "data/raw/ne_10m_urban_areas_landscan.zip"
    output:
        "data/processed/cities/{city}.npz",
    script:
        "scripts/get_city_arrays.py"

rule plot_place:
    input:
        "data/processed/{region}/{place}.npz"
    output:
        "outputs/{region}/{place}.pdf"
    script:
        "scripts/plot_place.py"

rule country_pops:
    input:
        "data/processed/{region}/{country}.npz"
    output:
        "data/processed/{region}/{country}_population.txt"
    script:
        "scripts/country_pops.py"

rule lat_long_pops:
    input:
        "data/processed/{region}/{place}.npz",
    output:
        "data/processed/{region}/{place}_cells.csv"
    script:
        "scripts/lat_long_pops.py"

rule make_microcells:
    input:
        "data/processed/{region}/{place}_cells.csv",
        "configs/{region}/{place}_parameters.json",
        "outputs/{region}/{place}.pdf"
    output:
        "data/processed/{region}/{place}_microcells.csv"
    script:
        "scripts/microcell_conversion.py"

# rule make_pop_dist:
#     input:
#         "data/raw/WPP2022_PopulationByAge5GroupSex_Medium.csv",
#         "configs/{region}/{place}_parameters.json"
#     output:
#         "data/processed/{region}/{place}_pop_dist.json"
#     script:
#         "scripts/get_pop_dist.py"

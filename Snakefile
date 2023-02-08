rule all:
    input:
        "outputs/Bermuda.pdf",
        "data/processed/country_populations.csv",
        "outputs/dag.pdf"

rule render_dag:
    input:
        "Snakefile"
    output:
        "outputs/dag.pdf"
    shell:
        "snakemake --dag | dot -Tpdf > outputs/dag.pdf"

rule get_country_arrays:
    input:
        "data/raw/GHS_POP_E2015_GLOBE_R2019A_4326_30ss_V1_0.tif",
        "data/raw/ne_10m_admin_0_countries_lakes.zip"
    output:
        "data/processed/countries.npz"
    script:
        "scripts/get_country_arrays.py"

rule plot_countries:
    input:
        "data/processed/countries.npz"
    output:
        "outputs/Bermuda.pdf"
    script:
        "scripts/plot_countries.py"

rule country_pops:
    input:
        "data/processed/countries.npz"
    output:
        "data/processed/country_populations.csv"
    script:
        "scripts/country_pops.py"

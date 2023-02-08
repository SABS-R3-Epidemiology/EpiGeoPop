# Country population data workflow

This repository is a snakemake workflow for getting population density data for arbitrary countries.
It uses population data from the [JRC Big Data Analytics Platform](https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_MT_GLOBE_R2019A/GHS_POP_E2015_GLOBE_R2019A_4326_30ss/V1-0/), border data from [Natural Earth](https://www.naturalearthdata.com/downloads/10m-cultural-vectors/), and is largely based on Adam Symington's [excellent blog post](https://towardsdatascience.com/creating-beautiful-population-density-maps-with-python-fcdd84035e06).
This workflow is motivated by extending [epiabm](https://github.com/SABS-R3-Epidemiology/epiabm) to other countries and is very much a work in progress.

## Running

Create virtual environment (recommended)

```
python -m venv venv
source venv/bin/activate
```

Install dependencies

```
pip install -r requirements.txt
```

Downlaod the raw data (See `data/README.md` for more information)

Run the snakemake pipeline

```
snakemake --cores 1
```

# Exploring the data

Check the `outputs` directory for example population density maps.
The image `outputs/dag.svg` shows the entire workflow.
The file `data/processed/country_pops.csv` contains estimated populations for each processed country.

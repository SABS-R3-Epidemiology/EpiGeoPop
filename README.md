# Country population data workflow

This repository is a snakemake workflow for getting population density data for arbitrary countries.
It uses population data from the [JRC Big Data Analytics Platform](https://jeodpp.jrc.ec.europa.eu/ftp/jrc-opendata/GHSL/GHS_POP_MT_GLOBE_R2019A/GHS_POP_E2015_GLOBE_R2019A_4326_30ss/V1-0/), border data from [Natural Earth](https://www.naturalearthdata.com/downloads/10m-cultural-vectors/), and is largely based on Adam Symington's [excellent blog post](https://towardsdatascience.com/creating-beautiful-population-density-maps-with-python-fcdd84035e06).
This workflow is motivated by extending [epiabm](https://github.com/SABS-R3-Epidemiology/epiabm) to other countries and is very much a work in progress.

## Running

The following shows how to setup and run the Snakemake pipeline.
By default, it will create the files for running a Luxembourg simulation, but the `Snakefile` can be modified to generated files for many countries, province/states, or cities.

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

```
bash prep.sh
```

Run the snakemake pipeline

```
snakemake --cores 1
```

# Exploring the data

Check the `outputs` directory for example population density maps.
The image `outputs/dag.svg` shows the entire workflow.
The file `data/processed/countries/Luxembourg_microcells.csv` contains the generated microcells, used for input to simulations such as epiabm.
The file `data/processed/countries/Luxembourg_pop_dist.json` contains the age distribution of populations.

# Running on other regions

The Snakefile contains commented out examples of other regions to show how to generate files for other countries, provinces, and cities.
These also require a configuration file which can be copied from similar files in the `configs` directory.

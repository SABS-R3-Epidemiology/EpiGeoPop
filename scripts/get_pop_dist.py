import json

# import rasterio
# from rasterio import mask as msk
import pandas as pd
# import geopandas as gpd
import numpy as np
# from shapely.geometry import mapping


def get_pop_dist(data_path, config_path, loc_path):

    year = 2022 # TODO: year in config

    with open(config_path, 'r') as f:
        config = json.loads(f.read())
    country_name = config['country']
    df = pd.read_csv(data_path)
    df = df[df['Location'] == country_name]
    df = df[df['Time'] == year]
    # print(df)
    ages_arr = [p for p in df['PopTotal']]
    under_80 = ages_arr[:16]
    plus_80 = sum(ages_arr[16:])
    age_groups = under_80 + [plus_80]
    age_props = [pop/(sum(age_groups)) for pop in age_groups]
    # print(age_props)

    # country_df = df.loc[df['ADMIN'] == country]

    with open(loc_path, 'w') as f:
        json.dump(age_props, f)

get_pop_dist(snakemake.input[0], snakemake.input[1], snakemake.output[0])

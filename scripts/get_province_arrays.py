import rasterio
from rasterio import mask as msk
# import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import mapping


def get_province_arrays(tif_path, borders_path, province, country_path):

    tif_file = rasterio.open(tif_path)
    ghs_data = tif_file.read()

    ghs_data[0][ghs_data[0] < 0.0] = 0.0

    df = gpd.read_file(borders_path)

    country_df = df.loc[df['iso_3166_2'] == province]
    # Mask
    country_array, transform = msk.mask(tif_file, [mapping(geom) for geom in country_df.geometry.tolist()], crop=True)
    transform = transform.to_gdal()
    with open(country_path, 'wb') as f:
        np.savez_compressed(f, country_array=country_array, transform=transform)

get_province_arrays(snakemake.input[0], snakemake.input[1], snakemake.wildcards[0], snakemake.output[0])

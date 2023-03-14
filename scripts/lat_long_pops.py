import numpy as np
from affine import Affine


def lat_long_pops(country_path, out_path):

    country_data = {}
    transform_data = {}

    with open(country_path, 'rb') as f:
        npzfile = np.load(f)
        country_array = npzfile['country_array']
        transform = Affine.from_gdal(*npzfile['transform'])

    csv_lines = 'longitude,latitude,population'
    for section in country_array:
        for i in range(len(section)):
            for j in range(len(section[i])):
                pop = section[i][j]
                if pop > 0:
                    coords = transform * (j, i)
                    long, lat = coords
                    csv_lines += f'\n{long},{lat},{pop}'
    with open(out_path, 'w') as f:
        f.write(csv_lines)

lat_long_pops(snakemake.input[0], snakemake.output[0])

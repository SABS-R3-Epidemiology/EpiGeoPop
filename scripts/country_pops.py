# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap, LogNorm


def count_pop(data):
    pop = 0
    for stripe in data[0]:
        for people in stripe:
            if people > 0:
                pop += people
    return pop

def country_pops(data_path, out_path):

    country_data = {}
    with open(data_path, 'rb') as f:
        npzfile = np.load(f)
        population = int(count_pop(npzfile['country_array']))

    with open(out_path, 'w') as f:
        f.write(population)


country_pops(snakemake.input[0], snakemake.output[0])

# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap, LogNorm


def plot_countries(data_path, out_path):

    country_data = {}

    with open(data_path, 'rb') as f:
        npzfile = np.load(f)
        for country in npzfile.files:
            country_data[country] = npzfile[country]

    ourcmap = cm.get_cmap('hot_r', 460)
    newcolors = ourcmap(np.linspace(0, 1, 460))
    background_colour = np.array([0.9882352941176471, 0.9647058823529412, 0.9607843137254902, 1.0])
    newcolors[:1, :] = background_colour
    newcmp = ListedColormap(newcolors)

    for country in country_data:
        fig, ax = plt.subplots(facecolor='#FCF6F5FF')
        fig.set_size_inches(14, 7)
        ax.imshow(country_data[country][0], norm=LogNorm(), cmap=newcmp)
        # ax.imshow(europe_array[0], cmap=newcmp)
        ax.axis('off')
        plt.savefig(f'outputs/{country}.pdf', dpi=1000)

plot_countries(snakemake.input[0], snakemake.output[0])

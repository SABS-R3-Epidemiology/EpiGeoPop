# import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap, LogNorm


def plot_place(data_path, out_path):

    country_data = {}

    with open(data_path, 'rb') as f:
        npzfile = np.load(f)
        country_data = npzfile['country_array']

    ourcmap = cm.get_cmap('plasma_r', 460)
    newcolors = ourcmap(np.linspace(0, 1, 460))
    background_colour = np.array([0.9882352941176471, 0.9647058823529412, 0.9607843137254902, 1.0])
    newcolors[:1, :] = background_colour
    newcmp = ListedColormap(newcolors)

    fig, ax = plt.subplots(facecolor='#FCF6F5FF')
    fig.set_size_inches(14, 7)
    ax.imshow(country_data[0], norm=LogNorm(), cmap=newcmp)
    ax.axis('off')
    plt.savefig(out_path, dpi=1000)

plot_place(snakemake.input[0], snakemake.output[0])

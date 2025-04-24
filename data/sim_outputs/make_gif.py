import argparse
import math
import glob
import os
from PIL import Image

import pandas as pd
import matplotlib.animation
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from tqdm import tqdm


parser = argparse.ArgumentParser()
parser.add_argument("-f", "--sim_file", help="CSV file output by a simulation such as epiabm")
parser.add_argument("--duration", default=100, help="Time (ms) per frame of the GIF")
parser.add_argument("--dpi", default=300, help="Generated image resolution")
args = parser.parse_args()


def generate_colour_map(df, name, min_value=0.0, cmap=cm.Reds):
    """Generates a given color map, with the max value determined by
    the max value in a given column of the provided dataframe.
    """
    max_inf = math.log10(max(df[name])+1)
    cmap.set_under('lightgrey')
    norm = matplotlib.colors.Normalize(vmin=min_value, vmax=max_inf, clip=False)
    return cm.ScalarMappable(norm=norm, cmap=cmap)


def add_colorbar(im, width=None, pad=None, **kwargs):
    # From https://stackoverflow.com/a/76378778
    l, b, w, h = im.axes.get_position().bounds       # get boundaries
    width = width or 0.1 * w                         # get width of the colorbar
    pad = pad or width                               # get pad between im and cbar
    fig = im.axes.figure                             # get figure of image
    cax = fig.add_axes([l + w + pad, b, width, h])   # define cbar Axes
    # return fig.colorbar(im, cax=cax, **kwargs)       # draw cbar
    return fig.colorbar(cax=cax, **kwargs)       # draw cbar


def cbar_ticks(mapper):
    max_val = mapper.norm.vmax
    ticks = list(range(math.ceil(mapper.norm.vmax)))
    labels = [10**tick for tick in ticks]
    ticks = [math.log10(10**tick+1) for tick in ticks]
    ticks = [0] + ticks
    labels = [0] + labels
    return ticks, labels


def render_frame(ax, df, i, time, name, mapper, save_path='.', snaps=False):
    rows = df[df['time'] == time]
    img_grid = [[-1] * len(x_locs)] * len(y_locs)
    img_grid = np.array(img_grid)
    rows = zip(
        rows['location_x'],
        rows['location_y'],
        rows[name]
    )
    for row in rows:
        coords = f'{row[1]}-{row[0]}'
        idx = coord_map[coords]
        img_grid[idx] = math.log10(row[2]+1)
    ax.set_title(f'Time = {time}')
    im = ax.imshow(img_grid, norm=mapper.norm, cmap=mapper.get_cmap())
    if i == 0 and not snaps:
        ticks, labels = cbar_ticks(mapper)
        cbar = add_colorbar(im, mappable=mapper, ticks=ticks)
        cbar.ax.set_yticklabels(labels)
        cbar.set_label("Total infections")
    if not os.path.exists(save_path) and not snaps:
        os.makedirs(save_path)

    if not snaps:
        plt.savefig(f'{save_path}/frame-{i:03d}.png', bbox_inches='tight', dpi=im.axes.figure.dpi*4)


def make_gif(
        df,
        times,
        name='InfectionStatus.InfectMild',
        save_path = 'animation'
    ):
    file_names = []

    mapper = generate_colour_map(df, name=name)
    ax = plt.axes()

    plt.axis('off')
    pbar = tqdm(range(len(times)))
    for i in pbar:
        pbar.set_description(f'Generating frame {i}')
        t = times[i]
        render_frame(ax, df, i, t, name, mapper, save_path)

    fp_in = f"{save_path}/frame-*.png"
    fp_name = args.sim_file.split('/')[-1].replace('.csv', '')
    fp_out = f"{save_path}/{fp_name}.gif"
    img, *imgs = [Image.open(f).convert("RGB")
                  for f in sorted(glob.glob(fp_in))]
    img.save(
        fp=fp_out,
        format="GIF",
        append_images=imgs,
        save_all=True,
        duration=args.duration,
        loop=0,
        optimise=True,
    )
    for file in glob.glob(fp_in):  # Delete images after use
        os.remove(file)


def make_snaps(
        df,
        times,
        grid_dim,
        name='InfectionStatus.InfectMild',
        save_path = 'animation',
    ):
    file_names = []

    # Generate colour map to use
    mapper = generate_colour_map(df, name=name)

    # Configure subplot grid
    fig, axs = plt.subplots(grid_dim[0], grid_dim[1], sharex=True, sharey=True)
    axs = axs.ravel()

    # Determine time points to use
    plot_num = math.prod(grid_dim)

    times = [time for time in times if time % math.ceil(len(times)/plot_num) == 0]
    pbar = tqdm(range(len(times)))
    # Add each subplot
    for i in pbar:
        pbar.set_description(f'Generating frame {i}')
        t = times[i]
        render_frame(axs[i], df, i, t, name, mapper, save_path, snaps=True)
        axs[i].axis('off')

    plt.tight_layout()
    ticks, labels = cbar_ticks(mapper)
    cbar = plt.colorbar(mapper, ax=axs.tolist(), ticks=ticks)
    cbar.ax.set_yticklabels(labels)
    cbar.set_label("Total infections")
    fp_name = args.sim_file.split('/')[-1].replace('.csv', '')
    fp_out = f"{save_path}/{fp_name}_grid.png"
    plt.savefig(fp_out, bbox_inches='tight', dpi=args.dpi)


if __name__ == '__main__':
    df = pd.read_csv(args.sim_file)
    x_locs = sorted(list(set(df['location_x'])))
    y_locs = sorted(list(set(df['location_y'])))

    # Map geo coordinates to pixel indices
    coord_map = {}
    for i in range(len(x_locs)):
        x = x_locs[i]
        for j in range(len(y_locs)):
            y = y_locs[j]
            coord_map[f'{y}-{x}'] = (len(y_locs) - j - 1, i)

    names = [
        'InfectionStatus.InfectASympt',
        'InfectionStatus.InfectMild',
        'InfectionStatus.InfectGP',
        'InfectionStatus.InfectHosp',
        'InfectionStatus.InfectICU',
    ]

    # Get all infections
    df['infections'] = df[names].sum(axis=1)

    # Sum infections over all age groups, etc.
    df = df.groupby(['time', 'location_x', 'location_y'], as_index=False).sum()

    times = sorted(list(set(df['time'])))
    make_gif(df, times, name='infections')
    make_snaps(df, times, (3,3), name='infections')


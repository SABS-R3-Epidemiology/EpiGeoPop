import json
import math
import os
import numpy as np
import pandas as pd


def min_separation(x_loc, y_loc):
    """Returns minimum separation between any pair of locations"""
    assert len(x_loc) == len(y_loc), "Mismatched coordinates"
    sep = np.inf

    for x1, y1 in zip(x_loc, y_loc):
        for x2, y2 in zip(x_loc, y_loc):
            if (x1 == x2) and (y1 == y2):
                continue
            if np.linalg.norm([x2-x1, y2-y1]) < sep:
                sep = np.linalg.norm([x2-x1, y2-y1])
    return sep


def min_separation_lazy(x_loc, y_loc):
    """Returns minimum separation between the middle point and all other points
    Should be identitical to min_separation for regular grids but computes in O(n)
    instead of O(n^2)
    """
    assert len(x_loc) == len(y_loc), "Mismatched coordinates"
    sep = np.inf
    mid_idx = int(len(x_loc)/2)

    x1, y1 = x_loc[mid_idx], y_loc[mid_idx]
    for x2, y2 in zip(x_loc, y_loc):
        if (x1 == x2) and (y1 == y2):
            continue
        if np.linalg.norm([x2-x1, y2-y1]) < sep:
            sep = np.linalg.norm([x2-x1, y2-y1])
    return sep


def make_csv_row(data):
    return ','.join(str(i) for i in data) + '\n'


def make_microcells(population_path, config_path, out_path):
    with open(config_path, 'r') as f:
        config = json.loads(f.read())
    mcell_num_per_cell = config['mcell_num_per_cell']
    ave_places_per_mcell = config['ave_places_per_mcell']
    # Household count - based on average household size
    hh_freq = config['household_size_distribution']

    np.random.seed(42)

    columns = ["cell", "microcell", "location_x", "location_y",
               "household_number", "place_number", "Susceptible"]

    df = pd.read_csv(population_path)

    print('Calculating the grid size...')
    delta = min_separation_lazy(df['longitude'], df['latitude'])
    print('Done')

    f = open(out_path, 'w')
    f.write(make_csv_row(columns))
    print('Creating microcells...')
    for cell_index, row in df.iterrows():
        if cell_index % int(len(df) / 10) == 0:
            print('Creating cell {}/{}'.format(cell_index, len(df)))
        grid_len = math.ceil(math.sqrt(mcell_num_per_cell))
        m_pos = np.linspace(0, 1, grid_len)

        # Multinomial population distribution
        mcell_pop = int(row["population"] / mcell_num_per_cell)
        p = [1 / mcell_num_per_cell] * mcell_num_per_cell
        mcell_split = np.random.multinomial(row["population"], p, size=1)[0]

        ave_size = np.sum(np.multiply(np.array(range(1, len(hh_freq) + 1)),
                                      hh_freq))

        for n in range(mcell_num_per_cell):
            x = (row["longitude"]
                 + (m_pos[n % grid_len] - 0.5) * delta / grid_len)
            y = (row["latitude"]
                 + (m_pos[n // grid_len] - 0.5) * delta / grid_len)

            data_dict = {"cell": cell_index,
                         "microcell": n,
                         "location_x": x,
                         "location_y": y,
                         "Susceptible": mcell_split[n],
                         "place_number": np.random.poisson(ave_places_per_mcell),
                         "household_number": math.ceil(mcell_split[n] / ave_size)}

            f.write(make_csv_row(data_dict[c] for c in columns))

    f.close()
    print('Done')


make_microcells(snakemake.input[0], snakemake.input[1], snakemake.output[0])

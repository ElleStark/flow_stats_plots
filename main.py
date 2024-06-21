"""
Script to compute and plot flow field statistics, including mean velocities, integral scales, and turbulence intensity & isotropy.
Separate script included in repository to compute integral scales using the HPC resources (Blanca compute node at CU Boulder).
Elle Stark, May 2024
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np
import utils
plt.rcParams["mathtext.fontset"] = "dejavuserif"

def main():
    # Read in data from HDF5 file
    filename = 'D:/Re100_0_5mm_50Hz_16source_FTLE_manuscript.h5'
    with h5py.File(filename, 'r') as f:
        # x and y grids for plotting
        x_grid = f.get(f'Model Metadata/xGrid')[:].T
        y_grid = f.get(f'Model Metadata/yGrid')[:].T

        # u and v velocity field data
        u_data = f.get('Flow Data/u')[:].transpose(0, 2, 1)
        v_data = f.get('Flow Data/v')[:].transpose(0, 2, 1)

        # Spatial resolution
        dx = f.get('Model Metadata/spatialResolution')[0].item()

    # Decompose velocity fields into mean and fluctuating components of u and v
    u_mean, u_flx = utils.reynolds_decomp(u_data, time_ax=0)
    v_mean, v_flx = utils.reynolds_decomp(v_data, time_ax=0)
    
    # Calculate mean, standard deviation, RMS of field of interest




if __name__ == '__main__':
    main()


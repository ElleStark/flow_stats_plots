"""
Script to compute and plot flow field statistics, including mean velocities, integral scales, and turbulence intensity & isotropy.
Separate script included in repository to compute integral scales using the HPC resources (Blanca compute node at CU Boulder).
Elle Stark, May 2024
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np
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

    # Calculate mean, standard deviation, RMS of field of interest
    data = u_data
    mean_field = np.mean(data, axis=0, keepdims=True)
    flx_field = data - mean_field  # fluctuating component of velocity
    mean_v = np.mean(v_data, axis=0)
    mean_u = np.mean(u_data, axis=0)



if __name__ == '__main__':
    main()


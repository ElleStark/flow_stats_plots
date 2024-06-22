"""
Script to compute and plot flow field statistics, including mean velocities, integral scales, and turbulence intensity & isotropy.
Separate script included in repository to compute integral scales using the HPC resources (Blanca compute node at CU Boulder).
Elle Stark, May 2024
"""

import h5py
import matplotlib.pyplot as plt
import numpy as np
import utils
import cmasher as cmr
import scipy.io
plt.rcParams["mathtext.fontset"] = "dejavuserif"

def main():
    # Read in data from HDF5 file
    filename = 'E:/Re100_0_5mm_50Hz_16source_FTLE_manuscript.h5'
    with h5py.File(filename, 'r') as f:
        # x and y grids for plotting
        x_grid = f.get(f'Model Metadata/xGrid')[:].T
        y_grid = f.get(f'Model Metadata/yGrid')[:].T

        # u and v velocity field data
        # u_data = f.get('Flow Data/u')[:].transpose(0, 2, 1)
        # v_data = f.get('Flow Data/v')[:].transpose(0, 2, 1)

        # # Spatial resolution
        # dx = f.get('Model Metadata/spatialResolution')[0].item()

    # Decompose velocity fields into mean and fluctuating components of u and v
    # u_mean, u_flx = utils.reynolds_decomp(u_data, time_ax=0)
    # v_mean, v_flx = utils.reynolds_decomp(v_data, time_ax=0)

    # Compute turbulence intensity


    # # PLOT: mean velocity field (consistent axes version)
    # cmap = cmr.waterlily_r
    # utils.plot_field_xy(x_grid, np.flipud(y_grid), v_mean, cmap=cmap, title='mean v', filepath='ignore/v_mean.png', dpi=600)

    # # PLOTS: integral scales
    # fname = 'ignore/ILS_u_cross_stream.npy'
    # ils_ustream = np.load(fname)
    # cmap = cmr.gothic_r
    # utils.plot_field_xy(x_grid, y_grid, ils_ustream,  title='Integral Length Scale, u cross-stream', cmap=cmap, range=[0, 0.15], filepath='ignore/ils_u_cross_stream.png', dpi=600)
    its = scipy.io.loadmat('ignore/Tux_array_u_update.mat')
    its = np.array(its['Tux_array'])
    cmap = cmr.jungle_r
    utils.plot_field_xy(x_grid, y_grid, its, title='Integral Time Scale, u', cmap=cmap, range=[0, 1.5], filepath='ignore/its_u.png', dpi=600, smooth=True, trimmed=True)


if __name__ == '__main__':
    main()


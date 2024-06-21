# Utility functions for computing relavent flow statistics for turbulent flow field
# Elle Stark, University of Colorado Boulder
# June 2024
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

def reynolds_decomp(flowfield, time_ax=0):
    """decomposes time series of 2D flowfield data into its mean and fluctuating components

    Args:
        flowfield (nd-array): stack of flow data of interest.
        time_ax (int, optional): Specify time axis along which to calculate the mean. Defaults to 0.

    Returns:
        decomp_fields (list): [mean field, fluctuating field], where mean field is 2D array of flow averaged over time and 
            fluctuating field is 3D array of fluctuating component of flow at each timestep.
    """
    mean_field = np.mean(flowfield, axis=time_ax, keepdims=True)
    flx_field = flowfield - mean_field
    mean_field = mean_field[0, :, :]

    decomp_fields = [mean_field, flx_field]
    return decomp_fields

def plot_field_xy(x_grid, y_grid, field, cmap, title, filepath='plot.png', colorbar=True, save=True, dpi=300, vecs=False, field2=None):
    """plot 2D field with specified colormap and, if save=True, save at high resolution.

    Args:
        x_grid (2d-array): x locations 
        y_grid (2d-array): y locations
        field (2d-array): values for field of interest
        cmap (matlap color map): color map to use for plotting.
        title (string): plot title
        filepath (string, optional): path to desired save location. Defaults to 'plot.png'.
        colorbar (bool, optional): whether to include color bar on figure. Defaults to True.
        save (bool, optional): whether to save plot to file. Defaults to true.
        dpi (int, optional): resolution at which to save plot. Defaults to 300.
    """
    fig, ax = plt.subplots(figsize=(5.9, 4))

    if vecs:
        plt.streamplot(x_grid, y_grid, field, field2, density=1, linewidth=0.5, color='white')

    # define color map
    vmin = np.min([np.min(field), -np.max(field)])
    vmax = np.max([np.max(field), -np.min(field)])

    norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    plt.pcolormesh(x_grid, y_grid, field, cmap=cmap, norm=norm)

    plt.xlabel(r'$x$ (m)')
    plt.ylabel(r'$y$ (m)')
    plt.axis('equal')
    plt.title(title)
    if colorbar:
        ticks = [vmin, 3*vmin/4, vmin/2, vmin/4, 0, vmax/4, vmax/2, 3*vmax/4, vmax]
        plt.colorbar(ticks=ticks)
    if save:
        plt.savefig(filepath, dpi=dpi)
    else:
        plt.show()
    
# PLOTTING
#
# plt.close()
#
# # x-direction velocity
# # define u scale, with white at zero
# stat_u = mean_u
# # vmin = np.min(stat_u)
# vmin = -0.05
# # vmax = np.max(stat_u)
# vmax = 0.18
# norm = colors.TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
# # plot using pcolormesh
# fig, ax = plt.subplots(figsize=(5.9, 4))
# plt.pcolormesh(x_grid, y_grid, stat_u, cmap='PuOr', norm=norm)
# plt.xlabel(r'$x$ (m)')
# plt.ylabel(r'$y$ (m)')
# plt.axis('equal')
# plt.title(r'$\overline{u}$ (m/s)')
# # plt.title(r'$RMS_v$ (m/s)')
# plt.colorbar()
# plt.savefig('u_mean_fisherPlume.png', dpi=300)
# plt.show()


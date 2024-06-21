# Utility functions for computing relavent flow statistics for turbulent flow field
# Elle Stark, University of Colorado Boulder
# June 2024
import numpy as np
import matplotlib.pyplot as plt

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

def plot_field_xy(xgrid, ygrid, field, cmap, title, colorbar=True):
    fig, ax = plt.subplots(figsize=(5.9, 4))
    
    
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


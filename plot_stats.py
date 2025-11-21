# Script to quicklly generate plots for plume statistics
# Elle Stark October 2024

import h5py
import matplotlib.pyplot as plt
import numpy as np
import utils
import cmasher as cmr

filename = 'D:/singlesource_2d_extended/Re100_0_5mm_50Hz_singlesource_2d.h5'
with h5py.File(filename, 'r') as f:
    # x and y grids for plotting
    x_grid = f.get(f'Model Metadata/xGrid')[:].T
    y_grid = f.get(f'Model Metadata/yGrid')[:].T

# mean_u = np.load('D:/singlesource_2d_extended/mean_u_0to180s.npy')
# mean_v = np.load('D:/singlesource_2d_extended/mean_v_0to180s.npy')
# tke = np.load('ignore/extended_sim/tke_extendedsim.npy')
its = np.load('ignore/extended_sim/ITS_extendedSim_u_gauss_cut0.4.npy')
# print(f'ITS dims: {its.shape}')


# cmap = cmr.waterlily_r
# utils.plot_field_xy(x_grid, y_grid, mean_u[:-1, :-1].T, cmap=cmap, range=[-0.15, 0.15], title='mean u', filepath='ignore/extended_sim/u_mean_v2.png', dpi=600)
# utils.plot_field_xy(x_grid, y_grid, mean_v[:-1, :-1].T, cmap=cmap, range=[-0.15, 0.15], title='mean v', filepath='ignore/extended_sim/v_mean_v2.png', dpi=600)

cmap = cmr.torch
utils.plot_field_xy(x_grid, y_grid, its[:, :], cmap=cmap, range=[0, 1.0], title='Integral Time Scale, x-direction', filepath='ignore/extended_sim/its_u_colors.png', dpi=600, shading='gouraud')

# cmap = cmr.eclipse
# utils.plot_field_xy(x_grid, y_grid, tke, cmap=cmap, range=[0, 0.005], title='', filepath='ignore/extended_sim/tke.png', dpi=600, shading='gouraud')




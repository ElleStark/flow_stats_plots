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
        u_data = f.get('Flow Data/u')[:].transpose(0, 2, 1)
        v_data = f.get('Flow Data/v')[:].transpose(0, 2, 1)

        # Spatial resolution
        dx = f.get('Model Metadata/spatialResolution')[0].item()

    # Decompose velocity fields into mean and fluctuating components of u and v
    u_mean, u_flx = utils.reynolds_decomp(u_data, time_ax=0)
    v_mean, v_flx = utils.reynolds_decomp(v_data, time_ax=0)


    # Compute turbulent kinetic energy
    tke = 0.5 * (np.mean(u_flx**2, axis=0) + np.mean(v_flx**2, axis=0))
    # t_intensity = np.sqrt(tke) / np.sqrt(u_mean**2 + v_mean**2)

    # PLOT: turbulent intensity and turbulent kinetic energy
    cmap = cmr.ember
    utils.plot_field_xy(x_grid, y_grid, tke, title='turbulent kinetic energy', cmap=cmap, filepath='ignore/tke.png', dpi=600, trimmed=False)
    # utils.plot_field_xy(x_grid, y_grid, t_intensity, title='turbulence intensity', cmap=cmap, range=[0, 0.8], filepath='ignore/t_intensity_trimmed.png', dpi=600, trimmed=True)

    # PLOT: mean velocity field (consistent axes version)
    cmap = cmr.waterlily_r
    utils.plot_field_xy(x_grid, y_grid, u_mean, cmap=cmap, title='mean v', range=[-0.1785, 0.1785], filepath='ignore/v_mean_flipped.png', dpi=600)
    utils.plot_field_xy(x_grid, y_grid, v_mean, cmap=cmap, title='mean u', filepath='ignore/u_mean_flipped.png', dpi=600)

    # # PLOTS: integral scales
    # ils_uxstream = np.load('ignore/ILS_u_cross_stream.npy')
    # ils_vstream = np.load('ignore/ILS_vstreamwise.npy')
    ils_ustream = np.load('ignore/ILS_ustreamwise.npy')
    ils_vxstream = np.load('ignore/ILS_v_cross_stream.npy')

    cmap = cmr.gothic_r
    # utils.plot_field_xy(x_grid, y_grid, np.flipud(ils_uxstream),  title='Integral Length Scale, u cross-stream', cmap=cmap, range=[0, 0.15], filepath='ignore/ils_u_cross_stream_trimmed_alt.png', dpi=600, trimmed=True)
    # utils.plot_field_xy(x_grid, y_grid, np.flipud(ils_vstream),  title='Integral Length Scale, v streamwise', cmap=cmap, range=[0, 0.15], filepath='ignore/ils_v_streamwise_trimmed_alt.png', dpi=600, trimmed=True)
    utils.plot_field_xy(x_grid, y_grid, np.flipud(ils_ustream),  title='Integral Length Scale, u streamwise', cmap=cmap, range=[0, 0.15], filepath='ignore/ils_u_streamwise_flipped.png', dpi=600, trimmed=False)
    utils.plot_field_xy(x_grid, y_grid, np.flipud(ils_vxstream),  title='Integral Length Scale, v cross-stream', cmap=cmap, range=[0, 0.15], filepath='ignore/ils_v_cross_stream_flipped.png', dpi=600, trimmed=False)
 
    its_u = scipy.io.loadmat('ignore/Tux_array_u_update.mat')
    its_u = np.array(its_u['Tux_array'])
    its_v = scipy.io.loadmat('ignore/Tux_array_v.mat')
    its_v = np.array(its_v['Tux_array'])

    cmap = cmr.jungle_r
    utils.plot_field_xy(x_grid, y_grid, its_u, title='Integral Time Scale, u', cmap=cmap, range=[0, 1.5], filepath='ignore/its_u_trimmed.png', dpi=600, smooth=True, trimmed=True)
    utils.plot_field_xy(x_grid, y_grid, its_v, title='Integral Time Scale, v', cmap=cmap, range=[0, 1.5], filepath='ignore/its_v_trimmed_alt.png', dpi=600, smooth=True, trimmed=True)


    # Compute Reynolds Stress Tensor
    u_prime_sq = np.mean(u_flx**2, axis=0)
    v_prime_sq = np.mean(v_flx**2, axis=0)
    uv_prime = np.mean(u_flx*v_flx, axis=0)
    aspect_ratio = np.minimum(u_prime_sq, v_prime_sq)/np.maximum(u_prime_sq, v_prime_sq)
    print(aspect_ratio.shape)
    mean_ar = np.mean(aspect_ratio[123:723, :])
    print(f"mean aspect ratio: {mean_ar}")

    # PLOTS: components of Reynolds stress
    cmap = cmr.dusk_r
    utils.plot_field_xy(x_grid, y_grid, u_prime_sq, title='Reynolds Stress: u normal', cmap=cmap, filepath='ignore/rs_uprimesq.png', dpi=600)
    utils.plot_field_xy(x_grid, y_grid, v_prime_sq, title='Reynolds Stress: v normal', cmap=cmap, filepath='ignore/rs_vprimesq.png', dpi=600)
    utils.plot_field_xy(x_grid, y_grid, uv_prime, title='Reynolds Stress: uv shear', cmap=cmap, filepath='ignore/rs_uvprime.png', dpi=600)
    utils.plot_field_xy(x_grid, y_grid, u_prime_sq/v_prime_sq, title='Reynolds Stress Anisotropy: u_stress/v_stress', cmap=cmap, range=[1, 5], filepath='ignore/rs_uvanisotropy.png', dpi=600)
    utils.plot_field_xy(x_grid, y_grid, aspect_ratio, title='Reynolds Stress Anisotropy: Aspect Ratio', cmap=cmap, filepath='ignore/rs_normalaspectratio.png', dpi=600)

if __name__ == '__main__':
    main()


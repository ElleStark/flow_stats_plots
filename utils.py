# Utility functions for computing relavent flow statistics for turbulent flow field
# Elle Stark, University of Colorado Boulder
# June 2024
import numpy as np

def reynolds_decomp(flowfield, time_ax=0):
    mean_field = np.mean(flowfield, axis=time_ax, keepdims=True)
    flx_field = flowfield - mean_field
    mean_field = mean_field[0, :, :]

    decomp_fields = [mean_field, flx_field]
    return decomp_fields


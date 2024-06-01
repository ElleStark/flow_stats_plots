# Source code to run on CU Boulder HPC resources: Blanca cluster
# Computes integral length scales for multisource plume dataset for u and v in both x-stream and streamwise directions
# Elle Stark May 2024
#!/usr/bin/env python3

import h5py
import logging
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpi4py import MPI
import numpy as np
import time

# Set up logging for convenient messages
# Note, Blanca cluster runs Python 2.7 as of May 2024, which doesn't allow f-strings
logger = logging.getLogger('ilsmpipy')
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setFormatter(logging.Formatter("%(asctime)s:%(name)s:%(levelname)s:%(message)s"))
logger.addHandler(handler)
INFO = logger.info
WARN = logger.warn
DEBUG = logger.debug


# MPI setup and related data
comm_world = MPI.COMM_WORLD
rank = comm_world.Get_rank()
# number of processes will be determined from Slurm job script (.sh file)
# ntasks for this script ideally inteded to be ~1000 (using pre-emptable Blanca nodes) 
num_procs = comm_world.Get_size()
print('RUNNING ON {num_procs} CORES.'.format(num_procs=num_procs))


if (rank==0):
    # process 0 reads in data from PetaLibrary
    DEBUG('Reading data from PetaLibrary')

    start_t = time.time()

    file_path = '/pl/active/odor2action/Data/odorPlumeData/multisource-plume-simulations/hdf5/Re100_0_5mm_50Hz_16source_FTLE_manuscript.h5'

    with h5py.File(file_path, 'r') as f:
        u_data = f.get('Flow Data/u')[:-1, :, :].transpose(2, 1, 0)  # original: [time, x, y] array with dims [3001, 1001, 846]

    total_t = time.time() - start_t
    # INFO('u array read from PetaLibrary in {total_t} sec'.format(total_t=total_t))

    u_slice = u_data[0, 0, 1:10]
    # DEBUG('First 10 vals of u array: {u_slice}'.format(u_slice=u_slice))

    # Store output in rc_scratch
    # out_path = '/rc_scratch/elst4602/tests/uStatTests.npy'
    # np.save(out_path, u_slice)

    # FOR CASE WITH 1 TASK:

    # QC plot: original u_data at time 0
    plt.pcolormesh(u_data[:, :, 0])
    plt.savefig('/rc_scratch/elst4602/tests/u_orig.png')

    # reshape u_data into 1D array
    u_data = np.reshape(u_data, -1) 

    # QC plot: reshape and plot one timestep of u data
    u_check = np.reshape(u_data, (1001, 846, 3000))
    plt.pcolormesh(u_data[:, :, 0])
    plt.savefig('/rc_scratch/elst4602/tests/u_reshaped.png')


# Compute number of elements per spatial grid and per processor
x, y, z = u_data.size
tsteps_per_proc = z // num_procs
remainder = z % num_procs

# Define start and end indices for chunk for each process
start_idx = rank * tsteps_per_proc + min(rank, remainder)
end_idx = start_idx + tsteps_per_proc + (1 if rank < remainder else 0)

# Buffer array to fill in each process
u_chunk = np.empty((x, y, end_idx-start_idx), dtype='d')  # dtype='d' specifies double-precision floating point numbers, i.e. np.float64
 
# scatter u data from process 0 into variable 'u_data' variable in other processes
comm_world.Scatterv([u_data, (tsteps_per_proc + (1 if i < remainder else 0)) * x * y for i in range(num_procs)], u_chunk, root=0)

# INTEGRAL LENGTH SCALE COMPUTATION FOR EACH PROCESS:



# each process returns computed ILS matrix for the selected time step
#if (my_rank != 0):
 #   comm_world.gather(ils_col, dest=0, tag=1)
#else:
    # average all timesteps to get ILS at each point


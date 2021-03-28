import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm as colormaps
import scipy
import scipy.stats
import scipy.integrate
from numba import jit
from sklearn import mixture
from matplotlib.colors import LogNorm
import pickle

import mylib as my
import mdtraj as mdt

sqrt_2pi = np.sqrt(2 * np.pi)

def stdize(x):
    return (x - np.mean(x)) / np.std(x)

def unitize(x):
	return (x - np.min(x)) / (np.max(x) - np.min(x))

def files_exist(filenames):
    for f_name in filenames:
        if(not os.path.isfile(f_name)):
            return False
    return True

def get_present_elements(top, select='protein'):
    els = []
    for i in top.select(select):
        el = top.atom(i).element.symbol
        if(el not in els):
            els.append(el)
    return els

# =============== global paths ================
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'K')
res_path = os.path.join(root_path, 'res')

# ================ params =================
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns
Dt = dt * 5000     # ns
readframes_timestep = 1
supercell = np.array([1, 1, 2], dtype=np.intc)
traj_filename = 'npt_nojump.xtc'
pbctraj_filename = 'npt.xtc'
topol_filename = 'topol.top'
initial_pdb_filename = 'initial_npt.pdb'
gmx_exe_name = 'gmx_ser_20'

# ============== arg parse =================
supercell_str = ''.join([str(x) for x in supercell])

[Tmp, recomp, Zstep], _ = \
    my.parse_args(sys.argv[1:], ['-temp', '-recomp', '-Zstep'], \
                  possible_values=[None, ['0', '1'], None], \
                  possible_arg_numbers=[[1], [0, 1], [0, 1]], \
                  default_values=[None, ['0'], ['0.5']])
Tmp = float(Tmp[0])
recomp = (recomp[0] == '1')
Zstep = float(Zstep[0])

# =================== per-model paths ======================
model_name = 'flucts_temp' + my.f2str(T) + '_0'
model_path = os.path.join(run_path, model_name)
traj_filepath = os.path.join(model_path, traj_filename)
pbctraj_filepath = os.path.join(model_path, pbctraj_filename)
topol_filepath = os.path.join(model_path, topol_filename)
init_pdb_filepath = os.path.join(model_path, initial_pdb_filename)
wcrd_filepath = os.path.join(model_path, 'water_coord.npy')
data_files = [wcrd_filepath]
recomp = recomp or not files_exist(data_files)

# ======================== load data & initial process ==============================

init_traj = mdt.load(traj_filepath, stride=100000, top=init_pdb_filepath)
top = init_traj.topology
waterO_atoms_inds = top.select('water and name O')
if(verbose):
    print('top: ', top)

if(recomp):
#   ========================== load & preproc data =========================
    traj = mdt.load(traj_filepath, stride=readframes_timestep, top=init_pdb_filepath)
    N_frames = traj.xyz.shape[0]
    N_atoms = traj.xyz.shape[1]
    time = np.arange(0, N_frames) * Dt
    #timecut_ind = (time >= time_cut)
    #time = time[timecut_ind]
    w_crd = traj.xyz[:, waterO_atoms_inds, :]
    N_water = w_crd.shape[1]
    if(verbose):
        print('traj shape: ', traj.xyz.shape)

    with open(wcrd_filepath, 'wb') as f:
        np.save(f, w_crd)
else:
    with open(wcrd_filepath, 'rb') as f:
        w_crd = np.load(f)

    N_frames = w_crd.shape[0]
    N_water = w_crd.shape[1]
    #start_time_ind = np.ceil(time_cut / Dt)
    #time = np.arange(start_time_ind, start_time_ind + N_frames) * Dt

plt.show()
print('ALL DONE')

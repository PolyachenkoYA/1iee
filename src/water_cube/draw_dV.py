import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing
import matplotlib.pyplot as plt 

import mylib as my

import gromacs.formats as gmx

# =============== paths ================
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'water_cube')
res_path = os.path.join(root_path, 'res')
nvt0_filename = 'nvt'
nvt1_filename = 'nvt_bigbox'

def run_it(cmd, shell=False):
    print(cmd)
    sp.run(cmd, shell=shell)

def read_last_line(filename):
    with open(filename, 'rb') as f:
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        last_line = f.readline().decode()
        return last_line
    
def load_xvg(filepath, failsafe=True):
    if(failsafe):
        if(not os.path.isfile(filepath)):
            print('file "' + filepath + '" was not found. Aborting.')
            exit(1)
    xvg_file = gmx.XVG()
    xvg_file.read(filepath)
    return xvg_file
    
def proc_half(filepath, cut_time):
    xvg_filepath = filepath + '.xvg'
    xvg_file = load_xvg(xvg_filepath)  # timestep, T, P
    N_fields = len(xvg_file.names)
    time = xvg_file.array[0]
    stab_time_ind = (time < cut_time)
    cut_time_P = xvg_file.array[2][stab_time_ind]
    N_cut = len(cut_time_P)
    P_mean = np.mean(cut_time_P)
    P_std = np.std(cut_time_P)
    d_P = P_std / np.sqrt(N_cut)
    
    gro_filepath = filepath + '.gro'
    box_sizes_line = read_last_line(gro_filepath)
    sizes = np.float_(box_sizes_line.split())
    V = np.prod(sizes)
    
    return P_mean, P_std, d_P, V

def proc_model(time, id, cut_time_0=0, cut_time_1=0):
    model_name = 'dV_time' + str(time) + '_' + str(id)
    model_path = os.path.join(run_path, model_name)
    
    #print(time, id)
    # stab time is ~20-30
    P0, P0_std, d_P0, V0 = proc_half(os.path.join(model_path, nvt0_filename), 100)
    P1, P1_std, d_P1, V1 = proc_half(os.path.join(model_path, nvt1_filename), 100)
    
    dV_rel = 1 / (V1 / V0 - 1)  # V1 > V0
    K = (P0 - P1) * dV_rel
    d_K = np.sqrt(d_P1**2 + d_P0**2) * dV_rel
    
    return K * 1e5, d_K * 1e5  # atm -> Pa
    
    
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns

# ================ params =================

times = [0.25, 0.5, 1.0, 2.0]
ids = [[0, 1, 2, 3],
       [0, 1, 2, 3],
       [0, 1, 2, 3],
       [0, 1, 2]]

#times = [0.25, 0.5, 1.0]
#ids = [0, 1, 2, 3]

# ============== arg parse ====================
N = len(ids)
N_t = len(times)
K_arr = np.zeros((N_t, N))
d_K_arr = np.zeros((N_t, N))
K = np.zeros(N_t)
K_std = np.zeros(N_t)
d_K = np.zeros(N_t)
for t_i, time in enumerate(times):
    for i, id in enumerate(ids[t_i]):
        K_arr[t_i, i], d_K_arr[t_i, i] = proc_model(time, id)
        #print(K_arr[t_i, i], end=' ')
        #sp.run(['cp', K_filepath, './K_' + suff + '.dat'])
        #print()
#    print(time)
    K[t_i] = np.mean(K_arr[t_i, :])
    K_std[t_i] = np.std(K_arr[t_i, :])
    d_K[t_i] = K_std[t_i] / np.sqrt(len(K_arr[t_i, :]))

print(K)
print(d_K / K)
print(K_std / K)

#fig, ax = my.get_fig('traj length (ns)', 'K', title='K(t)')
#ax.plot(times, K_arr)
#ax.errorbar(times, K, yerr=d_K)

#plt.show()

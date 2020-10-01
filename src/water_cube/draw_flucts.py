import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing
import matplotlib.pyplot as plt 

import mylib as my

def run_it(cmd, shell=False):
    print(cmd)
    sp.run(cmd, shell=shell)
    
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns

# =============== paths ================
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'water_cube')
res_path = os.path.join(root_path, 'res')
main_mdp_filename_base = 'npt'
K_filename = 'K.dat'

# ================ params =================

times = [0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5]
times = [0.2500, 0.3053, 0.3727, 0.4551, 0.5558, 0.6786, 0.8286, 1.012, 1.235, 1.5085, 1.842, 2.249, 2.746, 3.3535, 4.095, 5.000]
times = [0.25, 0.5, 1.0, 2.0]
ids = [0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12]

# ============== arg parse ====================
N = len(ids)
N_t = len(times)
K_arr = np.zeros((N_t, N))
K = np.zeros(N_t)
d_K = np.zeros(N_t)
for t_i, time in enumerate(times):
    for i, id in enumerate(ids):
        suff = my.f2str(time) + '_' + str(id)
        model_path = os.path.join(run_path, 'time' + suff)
        K_filepath = os.path.join(model_path, K_filename)
        K_arr[t_i, i] = 1 / np.loadtxt(K_filepath)
        print(K_arr[t_i, i], end=' ')
        #sp.run(['cp', K_filepath, './K_' + suff + '.dat'])
        #print()
    print(time)
    #K[t_i] = np.mean(K_arr[t_i, :])
    #d_K[t_i] = np.std(K_arr[t_i, :])

print(K)
print(d_K)
#fig, ax = my.get_fig('traj length (ns)', 'K', title='K(t)')
#ax.plot(times, K_arr)
#ax.errorbar(times, K, yerr=d_K)

#plt.show()
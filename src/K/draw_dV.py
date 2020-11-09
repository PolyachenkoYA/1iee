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
exe_path = os.path.join(root_path, 'src', 'K')
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
    
def proc_half(filepath, cut_time, to_draw=False, title=None):
    cptsave_filepath = filepath + '_prev.cpt'
    if(not os.path.isfile(cptsave_filepath)):
        print('WARNING\nNo "' + cptsave_filepath + '" found. Computation was probably corrupted.')
    xvg_filepath = filepath + '.xvg'
    xvg_file = load_xvg(xvg_filepath)  # timestep, T, P
    N_fields = len(xvg_file.names)
    time = xvg_file.array[0] * 1e-3
    stab_time_ind = (time > cut_time)
    time_cut = time[stab_time_ind]
    P = xvg_file.array[2][:]
    Pcut = P[stab_time_ind]
    N_cut = len(Pcut)
    Pcut_mean = np.mean(Pcut)
    Pcut_std = np.std(Pcut)
    d_Pcut = Pcut_std / np.sqrt(N_cut)

    if(to_draw):
        if(title is None):
            title = filepath + ' | P(t)'
        fig, ax = my.get_fig('time (ns)', 'P (atm)', title=title)
        
        ax.plot(time, P, label='data')
        ax.plot([cut_time] * 2, [min(P), max(P)], label='$t_{stab} = ' + str(cut_time) + '$')
        ax.legend()
    
    gro_filepath = filepath + '.gro'
    box_sizes_line = read_last_line(gro_filepath)
    sizes = np.float_(box_sizes_line.split())
    V = np.prod(sizes)
    
    return Pcut_mean, Pcut_std, d_Pcut, V

def proc_model(Ttau, dV, cut_time_0=5, cut_time_1=5, to_draw=False):
    model_name = 'dV_Ttau' + my.f2str(Ttau) + '_dVmult' + my.f2str(dV)
    model_path = os.path.join(run_path, model_name)
    
    #print(temp, id)
    # stab time is ~20-30
    P0, P0_std, d_P0, V0 = \
        proc_half(os.path.join(model_path, nvt0_filename), \
                  cut_time_0, to_draw=to_draw, \
                  title=model_name + '/initial_box')

    P1, P1_std, d_P1, V1 = \
        proc_half(os.path.join(model_path, nvt1_filename), \
                  cut_time_1, to_draw=to_draw, \
                  title=model_name + '/big_box')
    
    dV_rel = (V1 - V0) / (V1 + V0) * 2  # V1 > V0
    K = - (P1 - P0) / dV_rel
    d_K = np.sqrt(d_P1**2 + d_P0**2) / dV_rel
    
    return K * 1e5, d_K * 1e5  # atm -> Pa
    

# ================ params =================
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns

'''
# ============================== single plot =======================
temp = 35.0
id = 0
to_draw = True
cut_time = 0.5

K, d_K = proc_model(temp, id, cut_time_0=cut_time, cut_time_1=cut_time, to_draw=to_draw)
print(K * 1e-9, d_K * 1e-9)
plt.show()
'''

# =============================== K(tau_T, dV) ========================
dVs = np.array([2.5, 3, 3.5, 4])
#dVs = [3]
Ttaus = np.array([16, 32, 64, 128])
tau_to_draw = (np.array([0, 0, 0, 0]) == 1)
dV_to_draw = (np.array([1, 0, 0, 0]) == 1)

to_draw_P_time = False

# ======== proc ========
N_tau = len(Ttaus)
N_dV = len(dVs)
K = np.zeros((N_dV, N_tau))
d_K = np.zeros((N_dV, N_tau))
for dV_i, dV in enumerate(dVs):
    for Ttau_i, Ttau in enumerate(Ttaus):
        K[dV_i, Ttau_i], d_K[dV_i, Ttau_i] = proc_model(Ttau, dV, to_draw=to_draw_P_time)

#    K[t_i] = np.mean(K_arr[t_i])
#    K_std[t_i] = np.std(K_arr[t_i])
#    d_K[t_i] = K_std[t_i] / np.sqrt(len(K_arr[t_i]))

print(K)
print(d_K / K)

# ======== draw ========
if(np.any(dV_to_draw)):
    fig_tau, ax_tau = my.get_fig(r'$\tau_T (ps)$', '$K$ (Pa)', title=r'$K(\tau_T, dV)$', xscl='log', yscl='log')
    for dV_i, dV in enumerate(dVs):
        if(dV_to_draw[dV_i]):
            ax_tau.errorbar(Ttaus, K[dV_i, :], yerr=d_K[dV_i, :], label=('dV/V = 1e-' + str(dV)), fmt='o', mfc='none')
    ax_tau.legend()

if(np.any(tau_to_draw)):
    fig_dV, ax_dV = my.get_fig(r'$log_{10}(dV/V)$', '$K$ (Pa)', title=r'$K(\tau_T, dV)$')
    for Ttau_i, Ttau in enumerate(Ttaus):
        if(tau_to_draw[Ttau_i]):
            ax_dV.errorbar(-dVs, K[:, Ttau_i], yerr=d_K[:, Ttau_i], label=(r'$\tau_T$ = ' + str(Ttau) + ' (ps)'), fmt='o', mfc='none')
    ax_dV.legend()

plt.show()


'''
# ========================== K(T) ================================
# ======== local params ========
temps = [0.1, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0]
#temps = [0.1, 15.0, 30.0, 45.0]
#temps = [5.0, 10.0, 20.0, 25.0, 35.0, 40.0, 50.0, 55.0]

ids =  [[0, 1, 2, 3],
        [0, 4],
        [0, 4],
        [0, 1, 2, 3],
        [0, 4],
        [0, 4],
        [0, 1, 2, 3],
        [0, 4],
        [0, 4],
        [0, 1, 2, 3],
        [0, 4],
        [0, 4]]

ids =  [range(8)]*8

#ids =  [[0, 1, 2, 3],
#        [0, 1, 2, 3],
#        [0, 1, 2, 3],
#        [0, 1, 2, 3]]

# ======== proc ========
N = len(ids)
N_t = len(temps)
K_arr = [[]] * N_t
d_K_arr = [[]] * N_t
K = np.zeros(N_t)
K_std = np.zeros(N_t)
d_K = np.zeros(N_t)
for t_i, temp in enumerate(temps):
    K_arr[t_i] = np.zeros(len(ids[t_i]))
    d_K_arr[t_i] = np.zeros(len(ids[t_i]))
    for i, id in enumerate(ids[t_i]):
        K_arr[t_i][i], d_K_arr[t_i][i] = proc_model(temp, id)
        #print(K_arr[t_i, i], end=' ')
        #sp.run(['cp', K_filepath, './K_' + suff + '.dat'])
        #print()

#    print(temp)
#    print(K_arr[t_i][i], d_K_arr[t_i][i])
    K[t_i] = np.mean(K_arr[t_i])
    K_std[t_i] = np.std(K_arr[t_i])
    d_K[t_i] = K_std[t_i] / np.sqrt(len(K_arr[t_i]))

print(K)
print(d_K / K)
print(K_std / K)

# ======== draw ========
fig, ax = my.get_fig('T (C)', 'K (Pa)', title='K(T)')
#ax.errorbar(temps, K, yerr=d_K, label='averaged', fmt='o', mfc='none')

K_buf = np.zeros(N_t)
d_K_buf = np.zeros(N_t)
for i in range(8):
    for t_i in range(N_t):
        K_buf[t_i] = K_arr[t_i][i]
        d_K_buf[t_i] = d_K_arr[t_i][i]
    ax.errorbar(temps, K_buf, yerr=d_K_buf, label=str(i), fmt='o', mfc='none')

ax.legend()
plt.show()
'''
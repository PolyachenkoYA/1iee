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
run_path = os.path.join(root_path, 'run', 'flucts_conv_smallCompr')
exe_path = os.path.join(root_path, 'src', 'K')
res_path = os.path.join(root_path, 'res')
npt_filename = 'npt'

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
            return None
    xvg_file = gmx.XVG()
    xvg_file.read(filepath)
    return xvg_file

def proc_series(data, inds):
    filt_data = data[inds]
    N = len(filt_data)
    avrg = np.mean(filt_data)
    std = np.std(filt_data)
    dlt = std / np.sqrt(N)
    return filt_data, avrg, std, dlt, N
    
def proc_model(P_tau, compr, time, model_id=0, temp=35.0, cut_time=1000):
    model_name = 'flucts_Ptau' + my.f2str(P_tau) + '_compr' + my.f2str(compr) + '_time' + my.f2str(time) + '_' + str(model_id)
    model_path = os.path.join(run_path, model_name)
    filepath = os.path.join(model_path, npt_filename)
    cptsave_filepath = filepath + '_prev.cpt'
    if(not os.path.isfile(cptsave_filepath)):
        print('WARNING\nNo "' + cptsave_filepath + '" found. Computation was probably corrupted.')
    xvg_filepath = filepath + '.xvg'
    xvg_file = load_xvg(xvg_filepath)  # timestep, P, V, T, rho
    if(xvg_file is None):
        return None, None, None
    N_fields = len(xvg_file.names)
    time = xvg_file.array[0]
    stab_time_ind = (time < cut_time)

    T_cut, T_mean, T_std, d_T, N_cut = proc_series(xvg_file.array[1], stab_time_ind)
    P_cut, P_mean, P_std, d_P, _ = proc_series(xvg_file.array[2] * 1e5, stab_time_ind)
    V_cut, V_mean, V_std, d_V, _ = proc_series(xvg_file.array[3] * 1e-27, stab_time_ind)
    rho_cut, rho_mean, rho_std, d_rho, _ = proc_series(xvg_file.array[4], stab_time_ind)
    
    gro_filepath = filepath + '.gro'
    box_sizes_line = read_last_line(gro_filepath)
    sizes = np.float_(box_sizes_line.split()) * 1e-9
    V0 = np.prod(sizes)

    K0 = None
    if(gromacs_provided):
        Kdat_filepath = os.path.join(model_path, 'K.dat')
        if(os.path.isfile(Kdat_filepath)):
            if(os.stat(Kdat_filepath).st_size > 0):
                K0 = np.loadtxt(Kdat_filepath)
            else:
                print('file "' + Kdat_filepath + '" is empty')
                K0 = 0
        else:
            print('no "' + Kdat_filepath + '" found')
    

    if(abs(V_std / V_mean) < sys.float_info.epsilon):
        K = sys.float_info.epsilon
        V_m4 = 0
        d_V_std2 = 0
        d_K = sys.float_info.epsilon
    else:
        K = V_mean * T_mean * 1.38e-23 / V_std**2
        if(K < 0):
            print('K<0: ',V_mean, T_mean)
        V_m4 = np.mean((V_cut - V_mean)**4)
        d_V_std2 = np.sqrt((V_m4 - (N_cut - 3) / (N_cut - 1) * V_std**4)/N_cut)
        d_K = K * np.sqrt((d_V/V_mean)**2 + (d_T/T_mean)**2 + (d_V_std2/V_std**2)**2)

    return K, d_K, K0 # atm -> Pa
    
# ================ params =================
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns
equil_maxsol_poly = [-2.9516, 1117.2]   # maxsol = np.polyval(equil_maxsol_poly, T), [T] = C (not K)
temps = np.array([0.1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55])
P_taus = np.array([10, 20, 50, 100, 200, 400, 800, 1600])
comprs = np.array([1e-9, 2e-9, 4e-9, 1e-8])
times = np.array([5.0, 10.0, 20.0, 40.0])

filt_peaks = False
gromacs_provided = True

# ============== arg parse ====================
N_temps = len(temps)
N_times = len(times)
N_comprs = len(comprs)
N_Ptaus = len(P_taus)
N_temps = 1
K_data = np.empty((N_times, N_comprs, N_Ptaus))
d_K_data = np.empty((N_times, N_comprs, N_Ptaus))
K0_data = np.empty((N_times, N_comprs, N_Ptaus))
for time_i, time in enumerate(times):
    for Ptau_i, P_tau in enumerate(P_taus):
        for compr_i, compr in enumerate(comprs):
            K_data[time_i, compr_i, Ptau_i], d_K_data[time_i, compr_i, Ptau_i], K0_data[time_i, compr_i, Ptau_i] = \
                proc_model(P_tau, compr, time, model_id=0, temp=35.0, cut_time=1000)
            #print(K_arr[t_i, i], end=' ')
            #sp.run(['cp', K_filepath, './K_' + suff + '.dat'])
            #print()

# ===== mean over identical IDs ====
#K[t_i] = np.mean(K_arr[t_i])
#K_std[t_i] = np.std(K_arr[t_i])
#d_K[t_i] = K_std[t_i] / np.sqrt(len(K_arr[t_i]))

print('all data:')
print(K_data)
print(d_K_data / K_data)
print(K0_data)

times = [5.0]
for time_i, time in enumerate(times):
        for compr_i, compr in enumerate(comprs):
            K_draw = K_data[time_i, compr_i, :]
            draw_ind = ~(np.equal(K_draw, None) | np.equal(K_draw, sys.float_info.epsilon))
            if(np.any(draw_ind)):
                if(filt_peaks):
                    K_mean = np.mean(K_draw[draw_ind])
                    for i, K in enumerate(K_draw):
                        if(draw_ind[i]):
                            #draw_ind[i] = draw_ind[i] and (max(K/K_mean, 1 / (K / K_mean + 0.001)) < 5.0)
                            draw_ind[i] = draw_ind[i] and (K < 2e9)

                fic, ax = my.get_fig(r'$\tau_P$', r'K (Pa)', r'K($\tau_P$)| time (ns) = ' + str(time) + ', compr = ' + str(compr) + ', T=$35C^\circ$', xscl='log', yscl='log')
                ax.errorbar(P_taus[draw_ind], K_draw[draw_ind], yerr=d_K_data[time_i, compr_i, draw_ind], fmt='o', mfc='none', label='mine')

                if(gromacs_provided):
                    K0_draw = K0_data[time_i, compr_i, :]
                    draw0_ind = ~np.equal(K0_draw, None)
                    ax.scatter(P_taus[draw_ind], K0_draw[draw_ind], marker='o', facecolor='none', edgecolor='red', label='gromacs')
                ax.legend()

plt.show()

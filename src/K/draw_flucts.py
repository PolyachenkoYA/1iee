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
#run_path = os.path.join(root_path, 'run', 'flucts_conv_smallCompr')
#run_path = os.path.join(root_path, 'run', 'flucts_conv_InitTemp&immediateMeasure')
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'K')
res_path = os.path.join(root_path, 'res')
npt_filename = 'npt'
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
            return None
    xvg_file = gmx.XVG()
    xvg_file.read(filepath)
    return xvg_file

def proc_dV_half(filepath, cut_time, to_draw=False, title=None):
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

#def proc_model(Ttau, dV, cut_time_0=5, cut_time_1=5, to_draw=False):
def proc_dV_model(temp, model_id=0, cut_time_0=5, cut_time_1=5, to_draw=False):
    #model_name = 'dV_Ttau' + my.f2str(Ttau) + '_dVmult' + my.f2str(dV)
    model_name = 'dV_temp' + my.f2str(temp)
    model_path = os.path.join(run_path, model_name)
    
    #print(temp, id)
    # stab time is ~20-30
    P0, P0_std, d_P0, V0 = \
        proc_dV_half(os.path.join(model_path, nvt0_filename), \
                  cut_time_0, to_draw=to_draw, \
                  title=model_name + '/initial_box')

    P1, P1_std, d_P1, V1 = \
        proc_dV_half(os.path.join(model_path, nvt1_filename), \
                  cut_time_1, to_draw=to_draw, \
                  title=model_name + '/big_box')
    
    dV_rel = (V1 - V0) / (V1 + V0) * 2  # V1 > V0
    K = - (P1 - P0) / dV_rel
    d_K = np.sqrt(d_P1**2 + d_P0**2) / dV_rel
    
    return K * 1e5, d_K * 1e5  # atm -> Pa


def proc_series(data, cut_time, time, xlbl='time (ns)', ylbl='', yk=1, title='', to_draw=False):
    tk = 1e3  # ps -> ns
    time = time / tk
    cut_time /= tk
    data = data / yk
    inds = (time > cut_time)
    filt_data = data[inds]
    N = len(filt_data)
    avrg = np.mean(filt_data)
    std = np.std(filt_data)
    dlt = std / np.sqrt(N)

    if(to_draw):
        if(not title):
            title = ylbl + '(t)'
        fig, ax = my.get_fig(xlbl, ylbl, title=title)
        ax.plot([min(time), max(time)], [avrg] * 2, '--', label=(r'mean = $' + my.f2str(avrg) + ' \pm ' + my.f2str(std) + r'$'), color='green')
        ax.plot([cut_time] * 2, [min(data), max(data)], label=('stab time = ' + my.f2str(cut_time)), color='red')
        ax.plot(time[inds], data[inds], '.', label='used data', color='black', markersize=3)
        ax.plot(time[~inds], data[~inds], '.', label='relax', color='black', markersize=1)
        ax.legend()

    return filt_data * yk, avrg * yk, std * yk, dlt * yk, N

#def proc_fluct_model(P_tau, compr, time, model_id=0, temp=35.0, cut_time=5, draw_T=False, draw_P=False, draw_V=False, draw_rho=False):
def proc_fluct_model(temp, model_id=0, cut_time=5, draw_T=False, draw_P=False, draw_V=False, draw_rho=False, gromacs_provided=False):
    #model_name = 'flucts_Ptau' + my.f2str(P_tau) + '_compr' + my.f2str(compr) + '_time' + my.f2str(time) + '_' + str(model_id)
    #model_name = 'flucts_temp' + my.f2str(temp) + '_' + str(model_id)
    model_name = 'watercube_T' + my.f2str(temp)
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
    time = xvg_file.array[0] * 1e-3
    stab_time_ind = (time < cut_time)

    #title=r'$\tau_P = ' + str(P_tau) + r'$, $\kappa_T = ' + str(compr) + r'$, 
    title = 'Temp = ' + my.f2str(temp)
    Ptau_to_draw = 512
    compr_to_draw = 3e-4
    T_cut, T_mean, T_std, d_T, N_cut = proc_series(xvg_file.array[2], cut_time, time, ylbl='T (K)', title=title, to_draw=(draw_T))
    P_cut, P_mean, P_std, d_P, _ = proc_series(xvg_file.array[3] * 1e5, cut_time, time, yk=1e5, ylbl='P (atm)', title=title, to_draw=(draw_P))
    V_cut, V_mean, V_std, d_V, _ = proc_series(xvg_file.array[4] * 1e-27, cut_time, time, yk=1e-27, ylbl=r'V ($nm^3$)', title=title, to_draw=(draw_V))
    rho_cut, rho_mean, rho_std, d_rho, _ = proc_series(xvg_file.array[5], cut_time, time, yk=1000, ylbl=r'$\rho (g/cm^3)$', title=title, to_draw=(draw_rho))
    
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
                K0 = sys.float_info.epsilon
        else:
            print('no "' + Kdat_filepath + '" found')
    

    if(abs(V_std / V_mean) < sys.float_info.epsilon):
        K = sys.float_info.epsilon
        V_m4 = 0
        d_V_std2 = 0
        d_K = sys.float_info.epsilon
    else:
        K = V_mean * T_mean * 1.38e-23 / V_std**2
        V_m4 = np.mean((V_cut - V_mean)**4)
        d_V_std2 = np.sqrt((V_m4 - (N_cut - 3) / (N_cut - 1) * V_std**4)/N_cut)
        d_K = K * np.sqrt((d_V/V_mean)**2 + (d_T/T_mean)**2 + (d_V_std2/V_std**2)**2)

    return K, d_K, K0 # atm -> Pa
    
# ================ general params =================
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns
equil_maxsol_poly = [-2.9516, 1117.2]   # maxsol = np.polyval(equil_maxsol_poly, T), [T] = C (not K)
temps = np.array([0.1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55])
temps = np.array([30.0])
P_taus = np.array([200, 400, 800, 1600, 3200, 6400])
P_taus = np.array([200, 400, 800, 1600, 3200, 6400, 12800, 25000, 50000, 100000, 200000, 400000, 800000, 1600000, 3200000, 6400000])
P_taus = np.array([4, 8, 16, 32, 64, 128, 256, 512])
comprs = np.array([2e-4, 3e-4, 4e-4])
comprs = np.array([3e-4])
times = np.array([20.0, 40.0])
stab_time = 0.5

# ================ K(T) ===================
draw_all = False
draw_T = False or draw_all
draw_P = False or draw_all
draw_V = True or draw_all
draw_rho = False or draw_all

do_flucts = True
do_dV = False

use_gmx_K = True
N_temps = len(temps)

units = 1e6
fig, ax = my.get_fig('T (C)', 'K (MPa)', title='K(T)')

if(do_flucts):
    fl_K = np.empty(N_temps)
    fl_d_K = np.empty(N_temps)
    fl_K_gmx = np.empty(N_temps)

    for temp_i, temp in enumerate(temps):
        fl_K[temp_i], fl_d_K[temp_i], fl_K_gmx[temp_i] = \
            proc_fluct_model(temp, model_id=0, cut_time=stab_time, gromacs_provided=use_gmx_K, \
                      draw_T=draw_T, draw_P=draw_P, draw_V=draw_V, draw_rho=draw_rho)

    print('flucts:')
    print(fl_K_gmx)
    print(fl_K)
    print(fl_d_K / fl_K)
    
    ax.errorbar(temps, fl_K / units, yerr=fl_d_K / units, label='flucts', fmt='o', mfc='none')

if(do_dV):
    dV_K = np.empty(N_temps)
    dV_d_K = np.empty(N_temps)
    
    for temp_i, temp in enumerate(temps):
        dV_K[temp_i], dV_d_K[temp_i] = \
            proc_dV_model(temp, model_id=0, cut_time_0=stab_time, cut_time_1=stab_time, to_draw=draw_all)

    print('\ndV:')
    print(dV_K)
    print(dV_K / dV_d_K)

    ax.errorbar(temps, dV_K / units, yerr=dV_d_K / units, label=r'$dV/V \sim 1e-3$', fmt='o', mfc='none')

#ax.set_yticks(np.arange(3500, 5001, step=500))
ax.legend()

plt.show()


'''
# ================ compr & tau_P search ==============
filt_peaks = True
gromacs_provided = True
draw_all = True
draw_T = False or draw_all
draw_P = False or draw_all
draw_V = False or draw_all
draw_rho = False or draw_all

N_temps = len(temps)
N_times = len(times)
N_comprs = len(comprs)
N_Ptaus = len(P_taus)
N_temps = 1
for time_i, time in enumerate(times):
    for Ptau_i, P_tau in enumerate(P_taus):
        for compr_i, compr in enumerate(comprs):
            K_data[time_i, compr_i, Ptau_i], d_K_data[time_i, compr_i, Ptau_i], K0_data[time_i, compr_i, Ptau_i] = \
                proc_model(P_tau, compr, time, model_id=0, temp=35.0, cut_time=stab_time, \
                          draw_T=draw_T, draw_P=draw_P, draw_V=draw_V, draw_rho=draw_rho)
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

for time_i, time in enumerate(times):
        for compr_i, compr in enumerate(comprs):
            K_draw = K_data[time_i, compr_i, :]
            draw_ind = ~(np.equal(K_draw, None) | np.less_equal(abs(K_draw), sys.float_info.epsilon * 5))
            if(np.any(draw_ind)):
                if(filt_peaks):
                    K_mean = np.mean(K_draw[draw_ind])
                    for i, K in enumerate(K_draw):
                        if(draw_ind[i]):
                            #draw_ind[i] = draw_ind[i] and (max(K/K_mean, 1 / (K / K_mean + 0.001)) < 5.0)
                            draw_ind[i] = draw_ind[i] and (K < 100e8)

                fig, ax = my.get_fig(r'$\tau_P$', r'K (Pa)', r'K($\tau_P$)| time (ns) = ' + str(time) + ', compr = ' + str(compr) + ', T=$35C^\circ$', xscl='log', yscl='log')
                ax.errorbar(P_taus[draw_ind], K_draw[draw_ind], yerr=d_K_data[time_i, compr_i, draw_ind], fmt='o', mfc='none', label='mine')

                if(gromacs_provided):
                    K0_draw = K0_data[time_i, compr_i, :]
                    draw0_ind = ~(np.equal(K0_draw, None) | np.less_equal(abs(K0_draw), sys.float_info.epsilon * 5))
                    ax.scatter(P_taus[draw0_ind], K0_draw[draw0_ind], marker='o', facecolor='none', edgecolor='red', label='gromacs')
                ax.legend()

plt.show()
'''

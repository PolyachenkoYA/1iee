import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing
import matplotlib.pyplot as plt
import scipy
import scipy.stats
import scipy.integrate
from numba import jit
from sklearn import mixture
from matplotlib.colors import LogNorm

import mylib as my
import mdtraj as mdt

sqrt_2pi = np.sqrt(2 * np.pi)

def files_exist(filenames):
    for f_name in filenames:
        if(not os.path.isfile(f_name)):
            return False
    return True

# =============== paths ================
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'water_cube')
res_path = os.path.join(root_path, 'res')

# ================ params =================

T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns
Dt = 0.5     # ns
traj_filename = 'npt_nojump.xtc'
topol_filename = 'topol.top'
initial_pdb_filename = 'initial_npt.pdb'

N_r2_bins = 50
r2_min = 7e0
r2_max = 1e3
R_cut = 0.75
S_cut = 0.5

# ============== arg parse ==========================
supercell_str = ''.join([str(x) for x in supercell])

[time_cut, verbose], _ = \
    my.parse_args(sys.argv[1:], ['-time_cut', '-verbose'], \
                  possible_values=[None, ['0', '1']], \
                  possible_arg_numbers=[[0, 1], [0, 1]], \
                  default_values=[['0'], ['1']])
time_cut = float(time_cut[0])
verbose = (verbose[0] == '1')

model_name = 'flucts_temp' + my.f2str(T) + '_0'
model_path = os.path.join(run_path, model_name)
traj_filepath = os.path.join(model_path, traj_filename)
topol_filepath = os.path.join(model_path, topol_filename)
init_pdb_filepath = os.path.join(model_path, initial_pdb_filename)
R_filename = os.path.join(model_path, 'R_time' + my.f2str(time_cut) + '.npy')
D_filename = os.path.join(model_path, 'D_time' + my.f2str(time_cut) + '.npy')
S_filename = os.path.join(model_path, 'S_time' + my.f2str(time_cut) + '_sgm' + my.f2str(sgm) + '.npy')
displs_filename = os.path.join(model_path, 'displs_time' + my.f2str(time_cut) + '.npy')
maxdispl_filename = os.path.join(model_path, 'maxdispl_time' + my.f2str(time_cut) + '.npy')
data_files = [R_filename, S_filename, D_filename, maxdispl_filename, displs_filename]
recomp = recomp or not files_exist(data_files)
if(recomp):
    if(not os.path.isfile(init_pdb_filepath)):
        os.chdir(model_path)
        my.run_it('gmx_mpi trjconv -s npt.gro -f npt.xtc -skip 1000000000 -o ' + initial_pdb_filename + ' < output_whole_sys0.in')
    if(not os.path.isfile(traj_filepath)):
        os.chdir(model_path)
        my.run_it('gmx_mpi trjconv -s npt.gro -f npt.xtc -pbc nojump -o ' + traj_filename + ' < output_whole_sys0.in')
    
    traj = mdt.load(traj_filepath, top=init_pdb_filepath)
    top = traj.topology
    N_frames = traj.xyz.shape[0]
    N_atoms = traj.xyz.shape[1]
    time = np.arange(0, N_frames) * Dt
    w_crd = traj.xyz[:, top.select('water'), :]
    N_water = w_crd.shape[1]
    if(verbose):
        print('traj shape: ', traj.xyz.shape)
        print('top: ', top)

    timecut_ind = time > time_cut
    all_displs = np.sum((w_crd[timecut_ind, :, :] - w_crd[0, :, :])**2, axis=2)
    time = time[timecut_ind]
    N_frames = len(time)
    max_displ = np.max(all_displs, axis=0)
    min_displ = np.min(all_displs, axis=0)
    r2_range = max_displ - min_displ
    displ_ind = np.argsort(max_displ)
    w_crd = w_crd[:, displ_ind, :]
    all_displs = all_displs[:, displ_ind]
    max_displ = max_displ[displ_ind]
#    if(sgm is None):
#        sgm = max(max_displ) / N_frames / 2
#    else:
#        sgm = float(sgm[0])
    
    D = np.empty((N_water,))
    R = np.empty((N_water,))
    R_pvalue = np.empty((N_water,))
    I2 = np.empty((N_water,))
    I4 = np.empty((N_water,))
    S = np.empty((N_water,))
    for i in range(N_water):
        linfit = np.polyfit(time, all_displs[:, i], 1)
        D[i] = linfit[0] / 6
        #R[i], R_pvalue[i] = scipy.stats.pearsonr(time, all_displs[:, i])
        R[i] = my_R(time, all_displs[:, i])
    
        #I2[i] = np.log(np.sum(all_displs[:, i]**2) / (N_frames * r2_range[i]**2 / 3))
        #I4[i] = np.log(np.sum(all_displs[:, i]**4) / (N_frames * r2_range[i]**4 / 5))
        S[i] = get_S(all_displs[:, i], sgm)
    
        if(verbose):
            if(not i%100):
                print('done: ' + my.f2str((i+1)/N_water * 100) + ' %      \r', end='')

    with open(S_filename, 'wb') as f:
        np.save(f, S)
    with open(R_filename, 'wb') as f:
        np.save(f, R)
    with open(D_filename, 'wb') as f:
        np.save(f, D)
    with open(maxdispl_filename, 'wb') as f:
        np.save(f, max_displ)
    with open(displs_filename, 'wb') as f:
        np.save(f, all_displs)
else:
    with open(R_filename, 'rb') as f:
        R = np.load(f)
    with open(S_filename, 'rb') as f:
        S = np.load(f)
    with open(D_filename, 'rb') as f:
        D = np.load(f)
    with open(maxdispl_filename, 'rb') as f:
        max_displ = np.load(f)
    with open(displs_filename, 'rb') as f:
        all_displs = np.load(f)
    N_water = len(R)
    N_frames = all_displs.shape[0]
    start_time_ind = np.ceil(time_cut / Dt)
    time = np.arange(start_time_ind, start_time_ind + N_frames) * Dt
    
R_min_ind = np.argmin(R)
R_max_ind = np.argmax(R)
S_min_ind = np.argmin(S)
S_max_ind = np.argmax(S)
#I2_max_ind = np.argmax(I2)

SR_inds, [S_mean, R_mean], ax_SR = clustering_2D(S, R, '$S$', '$R$', gauss_cut=gauss_cut, verbose=verbose)
ax_SD = clustering_2D(S, D, '$S$', '$D$', gauss_cut=gauss_cut, verbose=verbose, n_comp=1)
ax_Sr = clustering_2D(S, np.sqrt(max_displ), '$S$', '$r_{max}$', gauss_cut=gauss_cut, verbose=verbose, n_comp=1)

N_SR_cut = np.sum((R > R_cut) & (S > S_cut))
N_SR_clust = np.sum(SR_inds)
#N_SD = np.sum(SD_inds)
if(verbose):
    print('N_cut: ', N_SR_cut, '; fraction: ', N_SR_cut / N_water)
    print('N_clust: ', N_SR_clust, '; fraction: ', N_SR_clust / N_water)
else:
    print(N_SR_cut, N_SR_clust, N_water, S_mean, R_mean)
    #print(N_SR_clust / N_water)

if(verbose):
    if(draw_hists):
        #r2_hist, r2_bins = make_hist(max_displ, r2_min, r2_max, N_r2_bins, scl='log', title='$r_{final}^2$ ($nm^2$)')
        #D_hist, D_bins = make_hist(D, -1e10, 1e10, N_r2_bins, scl='linear', title='$D$ ($nm^2/ns$)')
        R_hist, R_bins = make_hist(R, -2, 2, N_r2_bins, scl='linear', title='R')
        #R_hist_cut, R_bins_cut = make_hist(R, R_cut, 2, N_r2_bins, scl='linear', title='R')
    
        #fig_I, ax_I = my.get_fig('$I_2$', '$I_4$')
        #ax_I.scatter(I2, I4, s=4)
        #I2_hist, I2_bins = make_hist(I2, -20, 20, 100, scl='linear', title='$I_2$')
        #I4_hist, I4_bins = make_hist(I4, -20, 20, 100, scl='linear', title='$I_4$')
    
        #S_hist, S_bins = make_hist(S, -20, 20, 100, scl='linear', title='$S$')

        ax_SR.plot([S_cut, S_cut], [min(R), max(R)], c='red')
        ax_SR.plot([min(S), max(S)], [R_cut, R_cut], c='red')
    
    if(draw_extremes):
        y_lbl = '$\sqrt{r^2}$ ($nm$)'
        x_lbl = '$t$ (ns)'
    
        #fig, ax = my.get_fig(x_lbl, y_lbl, title=('$I_2 = ' + my.f2str(I2[I2_max_ind]) + '$'))
        #ax.plot(time, all_displs[:, I2_max_ind])
    
        #fig_most_mobile, ax_most_mobile = my.get_fig(x_lbl, y_lbl, title='the most mobile water')
        #ax_most_mobile.plot(time, np.sqrt(all_displs[:, -1]))
    
        #fig_least_mobile, ax_least_mobile = my.get_fig(x_lbl, y_lbl, title='the least mobile water')
        #ax_least_mobile.plot(time, np.sqrt(all_displs[:, 0]))
    
        #fig_minR, ax_minR = my.get_fig(x_lbl, y_lbl, title='min R water')
        #ax_minR.plot(time, np.sqrt(all_displs[:, R_min_ind]))
    
        #fig_minS, ax_minS = my.get_fig(x_lbl, y_lbl, title=('max S water; $S = ' + my.f2str(S[S_max_ind]) + '$'))
        #ax_minS.plot(time, np.sqrt(all_displs[:, S_max_ind]))

        draw_Rcut(all_displs, time, S, S[S_max_ind], sgm)
        draw_Rcut(all_displs, time, R, R_cut, sgm)
        draw_Rcut(all_displs, time, S, S_cut, sgm)
    
    if(draw_Rcuts):
        y_lbl = '$\sqrt{r^2}$ ($nm$)'
        y_lbl = '$r^2$ ($nm^2$)'
        x_lbl = '$t$ (ns)'
        for r_cut in draw_Rcuts:
            draw_Rcut(all_displs, time, R, r_cut, sgm)

    plt.show()
'''
    T     Rcut   d
    0.1   0.82   0.05
    10    0.46   0.02
    15    0.75   0.5
'''

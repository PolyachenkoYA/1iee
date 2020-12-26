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

import mylib as my
import mdtraj as mdt

def make_hist(data, mn, mx, N_bins, title='data', scl='linear'):
    min_data = data.min()
    max_data = data.max()
    grid_min = max(mn, min_data)
    grid_max = min(mx, max_data)
    lin_scl = (scl == 'linear')
    if(lin_scl):
        bins = np.concatenate(([min_data], np.linspace(grid_min, grid_max, N_bins - 1), [max_data]))
    else:
        bins = np.concatenate(([min_data], np.logspace(np.log10(grid_min), np.log10(grid_max), N_bins - 1), [max_data]))

    hist = np.histogram(data, bins=bins)

    fig, ax = my.get_fig('index', title, yscl=scl)
    ax.plot(np.arange(len(data)), data, '.')

    fig_hist, ax_hist = my.get_fig(title, 'n', xscl=scl)
    ax_hist.hist(data, bins=bins, density=lin_scl)

    return hist, bins

#@jit
def diff_S(x, data, sgm):
    #p = np.mean(scipy.stats.norm.pdf(x, data[np.newaxis].T, sgm), axis=0)
    x_grid, data_grid = np.meshgrid(x, data)
    p = np.mean(np.exp(-((x_grid - data_grid) / sgm)**2 / 2), axis=0)
    return -np.multiply(p, np.log(p))

#@jit
def get_S(data, Nint=100):
    N = len(data)
    rng = max(data) - min(data)
    sgm = rng / N
    x1 = min(data) - 5 * sgm
    x2 = max(data) + 5 * sgm
#    S, S_err = scipy.integrate.quad(diff_S, x1, x2, args=(data, sgm))

    x_draw = np.linspace(x1, x2, Nint)
    #p = np.mean(scipy.stats.norm.pdf(x_draw, data[np.newaxis].T, sgm), axis=0)
    S = np.mean(diff_S(x_draw, data, sgm)) / np.sqrt(2 * np.pi)
    return S

def draw_Rcut(all_displs, time, R, R_cut, Ndgt=5):
    ind = np.argmin(np.abs(R - R_cut))
    data = all_displs[:, ind]
    N = len(data)
    rng = max(data) - min(data)

    fig, ax = my.get_fig('$t$ (ns)', '$r^2$ ($nm^2$)', title=r'$R \approx ' + my.f2str(R[ind], Ndgt) + '$')
    #ax.plot(time, np.sqrt(data))
    ax.plot(time, data)

    sgm = rng / N
    N_draw = 1000
    x_draw = np.linspace(min(data) - 5*sgm, max(data) + 5*sgm, N_draw)
    distr = np.mean(scipy.stats.norm.pdf(x_draw, data[np.newaxis].T, sgm), axis=0)
    #S, S_err = scipy.integrate.quad(diff_S, min(x_draw), max(x_draw), args=(data, sgm))
    S = get_S(data)

    fig, ax = my.get_fig('$r^2$', 'n', title=r'$S \approx ' + my.f2str(S, Ndgt) + '$')
    ax.plot(x_draw, distr)


# =============== paths ================
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'K')
res_path = os.path.join(root_path, 'res')

# ================ params =================

T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns
Dt = 0.5     # ns
supercell = np.array([1, 1, 2], dtype=np.intc)
traj_filename = 'npt_nojump.xtc'
topol_filename = 'topol.top'
initial_pdb_filename = 'initial_npt.pdb'

N_r2_bins = 50
r2_min = 3e0
r2_max = 3e2
R_cut = 0.7
S_cut = 0.0325
verbose = False

# ============== arg parse ==========================
supercell_str = ''.join([str(x) for x in supercell])

[T, draw_extremes, draw_hists, time_cut, draw_Rcuts], _ = \
    my.parse_args(sys.argv[1:], ['-temp', '-extremes', '-hists', '-time_cut', '-Rcut'], \
                  possible_values=[None, ['0', '1'], ['0', '1'], None, None], \
                  possible_arg_numbers=[[1], [0, 1], [0, 1], [0, 1], None], \
                  default_values=[None, ['0'], ['0'], ['10'], []])
T = float(T)
time_cut = float(time_cut[0])
draw_extremes = (draw_extremes[0] == '1')
draw_hists = (draw_hists[0] == '1')
draw_Rcuts = [float(r) for r in draw_Rcuts]

model_name = 'flucts_temp' + my.f2str(T) + '_0'
model_path = os.path.join(run_path, model_name)
traj_filepath = os.path.join(model_path, traj_filename)
topol_filepath = os.path.join(model_path, topol_filename)
init_pdb_filepath = os.path.join(model_path, initial_pdb_filename)

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

D = np.empty((N_water,))
R = np.empty((N_water,))
R_pvalue = np.empty((N_water,))
I2 = np.empty((N_water,))
I4 = np.empty((N_water,))
S = np.empty((N_water,))
for i in range(N_water):
    linfit = np.polyfit(time, all_displs[:, i], 1)
    D[i] = linfit[0] / 6
    R[i], R_pvalue[i] = scipy.stats.pearsonr(time, all_displs[:, i])

    I2[i] = np.log(np.sum(all_displs[:, i]**2) / (N_frames * r2_range[i]**2 / 3))
    I4[i] = np.log(np.sum(all_displs[:, i]**4) / (N_frames * r2_range[i]**4 / 5))
    S[i] = get_S(all_displs[:, i])

    if(verbose or True):
        if(not i%100):
            print('done: ' + my.f2str((i+1)/N_water * 100) + ' %      \r', end='')
R_min_ind = np.argmin(R)
I2_max_ind = np.argmax(I2)
Su = (S - np.mean(S)) / np.std(S)
Ru = (R - np.mean(R)) / np.std(R)

print('\n', np.sum((R > R_cut) & (S > S_cut)))
#print(np.sum(R > R_cut))

if(verbose):
    if(draw_hists):
    #    r2_hist, r2_bins = make_hist(max_displ, r2_min, r2_max, N_r2_bins, scl='log', title='$r_{final}^2$ ($nm^2$)')
    #    D_hist, D_bins = make_hist(D, -1e10, 1e10, 100, scl='linear', title='$D$ ($nm^2/ns$)')
    #    R_hist, R_bins = make_hist(R, -2, 2, 100, scl='linear', title='R')
    #    R_hist_cut, R_bins_cut = make_hist(R, R_cut, 2, 50, scl='linear', title='R')
    
        #fig_I, ax_I = my.get_fig('$I_2$', '$I_4$')
        #ax_I.scatter(I2, I4, s=4)
        #I2_hist, I2_bins = make_hist(I2, -20, 20, 100, scl='linear', title='$I_2$')
        #I4_hist, I4_bins = make_hist(I4, -20, 20, 100, scl='linear', title='$I_4$')
    
        #S_hist, S_bins = make_hist(S, -20, 20, 100, scl='linear', title='$S$')
        fig_SR, ax_SR = my.get_fig('$S$', '$R$')
        ax_SR.scatter(S, R, s=1)
    
    
    if(draw_extremes):
        y_lbl = '$\sqrt{r^2}$ ($nm$)'
        x_lbl = '$t$ (ns)'
    
        fig, ax = my.get_fig(x_lbl, y_lbl, title=('$I_2 = ' + my.f2str(I2[I2_max_ind]) + '$'))
        ax.plot(time, all_displs[:, I2_max_ind])
    
        fig_most_mobile, ax_most_mobile = my.get_fig(x_lbl, y_lbl, title='the most mobile water')
        ax_most_mobile.plot(time, np.sqrt(all_displs[:, -1]))
    
        fig_least_mobile, ax_least_mobile = my.get_fig(x_lbl, y_lbl, title='the least mobile water')
        ax_least_mobile.plot(time, np.sqrt(all_displs[:, 0]))
    
        fig_minR, ax_minR = my.get_fig(x_lbl, y_lbl, title='min R water')
        ax_minR.plot(time, np.sqrt(all_displs[:, R_min_ind]))
    
        draw_Rcut(all_displs, time, R, R_cut)
    
    if(draw_Rcuts):
        y_lbl = '$\sqrt{r^2}$ ($nm$)'
        y_lbl = '$r^2$ ($nm^2$)'
        x_lbl = '$t$ (ns)'
        for r_cut in draw_Rcuts:
            #draw_Rcut(all_displs, R, r_cut, x_lbl, y_lbl)
            draw_Rcut(all_displs, time, S, r_cut)

    plt.show()
'''
    T     Rcut   d
    0.1   0.82   0.05
    10    0.46   0.02
    15    0.75   0.5
'''

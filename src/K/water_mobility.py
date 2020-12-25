import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing
import matplotlib.pyplot as plt
import scipy

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

def draw_Rcut(all_displs, R, R_cut, x_lbl, y_lbl, Ndgt=5):
    ind = np.argmin(np.abs(R - R_cut))
    fig, ax = my.get_fig(x_lbl, y_lbl, title=r'$R \approx ' + my.f2str(R[ind], Ndgt) + '$')
    ax.plot(time, np.sqrt(all_displs[:, ind]))

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
R_cut = 0.46

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
print('traj shape: ', traj.xyz.shape)
print('top: ', top)

timecut_ind = time > time_cut
all_displs = np.sum((w_crd[timecut_ind, :, :] - w_crd[0, :, :])**2, axis=2)
time = time[timecut_ind]
max_displ = np.max(all_displs, axis=0)
displ_ind = np.argsort(max_displ)
w_crd = w_crd[:, displ_ind, :]
all_displs = all_displs[:, displ_ind]
max_displ = max_displ[displ_ind]

D = np.empty((N_water, 1))
R = np.empty((N_water, 1))
R_pvalue = np.empty((N_water, 1))
for i in range(N_water):
    linfit = np.polyfit(time, all_displs[:, i], 1)
    D[i] = linfit[0] / 6
    R[i], R_pvalue[i] = scipy.stats.pearsonr(time, all_displs[:, i])
R_min_ind = np.argmin(R)

if(draw_hists):
    r2_hist, r2_bins = make_hist(max_displ, r2_min, r2_max, N_r2_bins, scl='log', title='$r_{final}^2$ ($nm^2$)')
    D_hist, D_bins = make_hist(D, -1e10, 1e10, 100, scl='linear', title='D')
    R_hist, R_bins = make_hist(R, -2, 2, 100, scl='linear', title='R')
    R_hist_cut, R_bins_cut = make_hist(R, R_cut, 2, 50, scl='linear', title='R')

if(draw_extremes):
    y_lbl = '$\sqrt{r^2}$ ($nm$)'
    x_lbl = '$t$ (ns)'
    fig_most_mobile, ax_most_mobile = my.get_fig(x_lbl, y_lbl, title='the most mobile water')
    ax_most_mobile.plot(time, np.sqrt(all_displs[:, -1]))

    fig_least_mobile, ax_least_mobile = my.get_fig(x_lbl, y_lbl, title='the least mobile water')
    ax_least_mobile.plot(time, np.sqrt(all_displs[:, 0]))

    fig_minR, ax_minR = my.get_fig(x_lbl, y_lbl, title='min R water')
    ax_minR.plot(time, np.sqrt(all_displs[:, R_min_ind]))

    draw_Rcut(all_displs, R, R_cut, x_lbl, y_lbl)

if(draw_Rcuts):
    y_lbl = '$\sqrt{r^2}$ ($nm$)'
    x_lbl = '$t$ (ns)'
    for r_cut in draw_Rcuts:
        draw_Rcut(all_displs, R, r_cut, x_lbl, y_lbl)

'''
    T     Rcut   d
    0.1   0.82   0.05
    10    0.46   0.02
'''

plt.show()



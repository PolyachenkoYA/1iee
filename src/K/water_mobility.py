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

def comp_hist(data, mn=None, mx=None, N_bins=100, lin_scl=True, dens=False):
    min_data = data.min()
    max_data = data.max()
    grid_min = min_data if(mn is None) else max(mn, min_data)
    grid_max = max_data if(mx is None) else min(mx, max_data)
    bins = np.concatenate(([] if min_data == grid_min else [min_data], \
                           np.linspace(grid_min, grid_max, N_bins + 1) if lin_scl else np.logspace(np.log10(grid_min), np.log10(grid_max), N_bins + 1), \
                           [] if max_data == grid_max else [max_data]))

    hist, _ = np.histogram(data, bins=bins, density=dens)

    return hist, bins, grid_min, grid_max

def make_hist(data, mn=None, mx=None, N_bins=100, y_lbl='y', scl='linear', title=None, plot_vs_index=False, dens=False, draw=True):
    lin_scl = (scl=='linear')
    hist, bins, grid_min, grid_max = comp_hist(data, mn=mn, mx=mx, N_bins=N_bins, lin_scl=lin_scl, dens=dens)
    max_dens = np.max(hist)

    if(draw):
        if(plot_vs_index):
            fig, ax = my.get_fig('index', y_lbl, yscl=scl)
            ax.plot(np.arange(len(data)), data, '.')

        if(title is None):
            title = '$p(' + y_lbl + ')$'
        p_lbl = '$p_{dens}$' if(dens) else 'N'
        fig_hist, ax_hist = my.get_fig(y_lbl, p_lbl, xscl=scl, title=title)
        plt_hist, plt_bins, _ = ax_hist.hist(data, bins=bins, density=dens)

        if(mn is not None):
            if(grid_min < mn):
                ax_hist.plot([mn, mn], [0, max_dens], c='red')
        if(mx is not None):
            if(grid_min > mx):
                ax_hist.plot([mx, mx], [0, max_dens], c='red')
    else:
        ax_hist = None

    return hist, bins, ax_hist

def my_R(x, y):
    x = x - np.mean(x)
    y = y - np.mean(y)
    return np.sum(np.multiply(x,y)) / np.sqrt(np.sum(x**2) * np.sum(y**2))

def diff_p(x, data, sgm):
    dx = abs(x[1] - x[0])
    #p = np.mean(scipy.stats.norm.pdf(x, data[np.newaxis].T, sgm), axis=0)

    x_grid, data_grid = np.meshgrid(x, data)
    p = np.mean(np.exp(-((x_grid - data_grid) / sgm)**2 / 2), axis=0) / (sqrt_2pi * sgm)       # density
    #p = np.mean(np.exp(-(x_grid - data_grid)**2 / 2), axis=0) * (dx / sqrt_2pi)   # prob

    return p

#@jit
def diff_S(x, data, sgm):
    p = diff_p(x, data, sgm)
    p[p < 2 * np.finfo(float).eps] = 1   # p*ln(p) = 0 for such p, so we can as well put p=1 to get the same 0
    return -np.multiply(p, np.log(p))

#@jit
#def get_S(data, Nint=100):
def get_S(data, sgm, Nint=500, S_margin=3):
    N = len(data)
    x1 = min(data) - S_margin * sgm
    x2 = max(data) + S_margin * sgm
    rng = x2 - x1
    sgm = rng / N
    #sgm = min(rng / 2, sgm)
#    S, S_err = scipy.integrate.quad(diff_S, x1, x2, args=(data, sgm))

    x_draw = np.linspace(x1, x2, Nint)
    #p = np.mean(scipy.stats.norm.pdf(x_draw, data[np.newaxis].T, sgm), axis=0)
    S = np.trapz(diff_S(x_draw, data, sgm), x=x_draw)  # density
    #S = np.sum(diff_S(x_draw, data, sgm))          # prob

    return S - (1 + np.log(2 * np.pi * sgm**2)) / 2

def draw_Rcut(all_displs, time, x, x_cut, sgm, Ndgt=5, S_margin=3):
    ind = np.argmin(np.abs(x - x_cut))
    data = all_displs[:, ind]
    N = len(data)
    rng = max(data) - min(data)
    R, R_pvalue = scipy.stats.pearsonr(time, data)
    D, _, _, _ = scipy.linalg.lstsq(time[:, np.newaxis], data[:, np.newaxis])
    D /= 6
    log_D = np.log10(D / D_bulk[T])

    #fig, ax = my.get_fig('$t$ (ns)', '$\sqrt{r^2}$ ($nm$)', title=r'$i = ' + str(ind) + r'; R_{r^2} \approx ' + my.f2str(R, Ndgt) + '$')
    fig, ax = my.get_fig('$t$ (ns)', '$r^2$ ($nm^2$)', title=r'$i = ' + str(ind) + r'; logD = ' + my.f2str(log_D, Ndgt) + '$')
    #ax.plot(time, np.sqrt(data))
    ax.plot(time, data)
    ax.plot([0, max(time)], [0, 6*D*max(time)])
    #ax.plot(time, data)

    #sgm = rng / N
    N_draw = 1000
    x_draw = np.linspace(min(data) - S_margin*sgm, max(data) + S_margin*sgm, N_draw)
    #distr = np.mean(scipy.stats.norm.pdf(x_draw, data[np.newaxis].T, sgm), axis=0)
    #S, S_err = scipy.integrate.quad(diff_S, min(x_draw), max(x_draw), args=(data, sgm))
    distr = diff_p(x_draw, data, sgm)
    S = get_S(data, sgm, S_margin=S_margin, Nint=N_draw)

    fig, ax = my.get_fig('$r^2$', 'n', title=r'$i = ' + str(ind) + r'; S \approx ' + my.f2str(S, Ndgt) + '$')
    ax.plot(x_draw, distr)

def clustering_2D(x, y, x_lbl, y_lbl, n_comp=2, gauss_cut=1, N_X_grid=300, N_Y_grid=300, verbose=False, title=None):
    if(verbose):
        fig, ax = my.get_fig(x_lbl, y_lbl, title=title)
    else:
        ax = None

    if(n_comp == 1):
        if(verbose):
            ax.scatter(x, y, s=1)
        return ax
    else:
        N_water = len(x)
        clasifier = mixture.GaussianMixture(n_components=n_comp, covariance_type='full')
        X_train = np.array([x, y]).T
        clasifier.fit(X_train)

        labels = clasifier.predict(X_train).astype(bool)
        gauss_means = clasifier.means_
        gauss_covs = clasifier.covariances_
        free_water_mean_ind = np.argmax(gauss_means[:, 1])
        free_water_labels = (clasifier.predict(X_train) == free_water_mean_ind)
        free_water_mean = gauss_means[free_water_mean_ind, :]
        free_water_cov = gauss_covs[free_water_mean_ind, :, :]
        #free_probs = scipy.stats.multivariate_normal(mean=free_water_mean, cov=free_water_cov).pdf(X_train)
        inv_cov = np.linalg.inv(free_water_cov)
        gauss_dists = np.empty((N_water,))
        for i in range(N_water):
            v = X_train[i, :] - free_water_mean
            gauss_dists[i] = np.dot(v, np.dot(inv_cov, v))
        #free_water_inds = (gauss_dists < gauss_cut) & free_water_labels
        free_water_inds = (gauss_dists < gauss_cut)

        if(verbose):
            X, Y = np.meshgrid(np.linspace(min(x), max(x), N_X_grid), np.linspace(min(y), max(y), N_Y_grid))
            XX = np.array([X.ravel(), Y.ravel()]).T
            Z = -clasifier.score_samples(XX)
            Z = Z.reshape(X.shape)
            CS = ax.contour(X, Y, Z, norm=LogNorm(vmin=1.0, vmax=100.0),levels=np.logspace(-1, 2, 15))
            CB = plt.colorbar(CS, extend='both')

            ax.scatter(x, y, s=1, c='green', label='all')
            ax.scatter(x[free_water_inds], y[free_water_inds], s=1, c='blue', label='free')
            ax.legend(loc='best')

        return free_water_inds, free_water_mean, ax
#@jit
def get_D_dyn(w_crd, N_D, step, win):
    N_frames = w_crd.shape[0]
    N_water = w_crd.shape[1]
    D_dyn = np.empty((N_D, N_water))
    win_ind = np.arange(win) + 1
    win_time = win_ind[:, np.newaxis] * Dt
    for i in range(N_water):
        for j in np.arange(N_D):
            displ = (w_crd[j * step + win_ind, i, :] - w_crd[j * step, i, :])**2
            y = np.sum(displ, axis=1)[:, np.newaxis]
            D_dyn[j, i], _, _, _ = scipy.linalg.lstsq(win_time, y)

    return D_dyn / 6

def color_from_mob(filepath, inds, x, y, N_bins=100, max_vdw=2.3, min_vdw=0.8, vdw_k=10, to_draw=False):
    N = len(x);
    mob = unitize(stdize(x) + stdize(y))
    mob_sort_inds = np.argsort(mob)
    mob = mob[mob_sort_inds]
    mob_hist, mob_bins, mob_ax = make_hist(mob, -2, 2, N_bins, y_lbl='$mob$', dens=False, title=os.path.split(filepath)[1], draw=to_draw)

    mob_bins = (mob_bins[1:] + mob_bins[:-1]) / 2
    #vdw_rads = np.minimum(min_vdw + vdw_k * np.abs(mob - np.mean(mob))**2, max_vdw)
    vdw_rads = (np.max(mob_hist) / mob_hist[np.minimum(np.intc(np.floor(mob / abs(mob_bins[2] - mob_bins[1]))), len(mob_hist) - 1)]) * min_vdw
    vdw_rads = np.minimum(vdw_rads, max_vdw)
    inds = inds[mob_sort_inds]
    rgb_colors = colormaps.cool(mob)
    hex_colors = [matplotlib.colors.to_hex(c) for c in rgb_colors]

    with open(filepath, 'w') as f:
        #f.write('select :SOL@OW;\ncolor red sel;\nrepr sphere sel;\nvdwdefine 0.8 sel;\n~show sel;\n~select sel;\n\n')
        for i in range(N):
            ser_str = '@/serialNumber=' + str(inds[i])
            f.write('repr sphere %s;\nvdwdefine %f %s;\ncolor %s %s;\ndisplay %s;\n\n' % (ser_str, vdw_rads[i], ser_str, hex_colors[i], ser_str, ser_str))

    return mob, hex_colors

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
D_bulk = { 0.1: 2.0833,
           5.0: 2.3548,
          10.0: 2.6759,
          15.0: 2.9911,
          20.0: 3.3278,
          25.0: 3.7376,
          30.0: 4.1991,
          35.0: 4.4784,
          40.0: 4.8772,
          45.0: 5.3572,
          50.0: 5.7396,
          55.0: 6.21449}  # C : nm^2/ns
d_D_bulk = { 0.1: 0.03,
             5.0: 0.03,
            10.0: 0.016,
            15.0: 0.015,
            20.0: 0.04,
            25.0: 0.03,
            30.0: 0.04,
            35.0: 0.03,
            40.0: 0.01,
            45.0: 0.08,
            50.0: 0.016,
            55.0: 0.11}  # C : nm^2/ns

dt = 2e-6    # 1 fs = 1e-6 ns
Dt = 0.5     # ns
supercell = np.array([1, 1, 2], dtype=np.intc)
traj_filename = 'npt_nojump.xtc'
pbctraj_filename = 'npt.xtc'
topol_filename = 'topol.top'
initial_pdb_filename = 'initial_npt.pdb'
gmx_exe_name = 'gmx_mpi'

r2_min = 7e0
r2_max = 1e3
R_cut = 0.75
S_cut = 0.5
S_margin = 3

do_r2range = False
do_I = False
do_R = False
do_S = False
do_D = False
do_Ddyn = False
do_waterDist = True
readframes_timestep = 1

# ============== arg parse =================
supercell_str = ''.join([str(x) for x in supercell])

[T, draw_extremes, draw_hists, time_cut, draw_Rcuts, verbose, recomp, sgm, gauss_cut, D_timewindow, D_timestep, N_bins, draw_clustering], _ = \
    my.parse_args(sys.argv[1:], ['-temp', '-extremes', '-hists', '-time_cut', '-Rcut', '-verbose', '-recomp', '-sgm', '-gauss_cut', '-D_timewindow', '-D_timestep', '-N_bins' ,'-clustering'], \
                  possible_values=[None, ['0', '1'], ['0', '1'], None, None, ['0', '1'], ['0', '1'], None, None, None, None, None, ['0', '1']], \
                  possible_arg_numbers=[[1], [0, 1], [0, 1], [0, 1], None, [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1]], \
                  default_values=[None, ['0'], ['0'], ['0'], [], ['0'], ['0'], ['0.5'], ['5.0'], ['2'], ['1'], ['100'], ['0']])
T = float(T)
sgm = float(sgm[0])
N_bins = int(N_bins[0])
time_cut = float(time_cut[0])
gauss_cut = float(gauss_cut[0])
draw_extremes = (draw_extremes[0] == '1')
draw_hists = (draw_hists[0] == '1')
draw_clustering = (draw_clustering[0] == '1')
verbose = (verbose[0] == '1')
recomp = (recomp[0] == '1')
D_timewindow = int(D_timewindow[0])
D_timestep = int(D_timestep[0])
draw_Rcuts = [float(r) for r in draw_Rcuts]

# =================== per-model paths ======================
model_name = 'flucts_temp' + my.f2str(T) + '_0'
model_path = os.path.join(run_path, model_name)
traj_filepath = os.path.join(model_path, traj_filename)
pbctraj_filepath = os.path.join(model_path, pbctraj_filename)
topol_filepath = os.path.join(model_path, topol_filename)
init_pdb_filepath = os.path.join(model_path, initial_pdb_filename)
mobility_S_R_filepath = os.path.join(model_path, 'mob_S_R_time' + my.f2str(time_cut) + '_sgm' + my.f2str(sgm) + '.chc')
mobility_S_D_filepath = os.path.join(model_path, 'mob_S_D_time' + my.f2str(time_cut) + '_sgm' + my.f2str(sgm) + '.chc')
mobility_S_logD_filepath = os.path.join(model_path, 'mob_S_logD_time' + my.f2str(time_cut) + '_sgm' + my.f2str(sgm) + '.chc')
select_water_filepath = os.path.join(model_path, 'water_inds.chc')
R_filename = os.path.join(model_path, 'R_time' + my.f2str(time_cut) + '.npy')
D_filename = os.path.join(model_path, 'D_time' + my.f2str(time_cut) + '.npy')
Ddyn_filename = os.path.join(model_path, 'Ddyn_W' + my.f2str(D_timewindow) + '_s' + my.f2str(D_timestep) + '.npy')
S_filename = os.path.join(model_path, 'S_time' + my.f2str(time_cut) + '_sgm' + my.f2str(sgm) + '.npy')
displs_filename = os.path.join(model_path, 'displs_time' + my.f2str(time_cut) + '.npy')
p_w_dists_filename = os.path.join(model_path, 'p_w_dist.npy')
p_w_elementHist_filename = os.path.join(model_path, 'p_w_elementHist.npy')
p_w_mindistInds_filename = os.path.join(model_path, 'p_w_mindistInds.npy')
pNoH_w_dists_filename = os.path.join(model_path, 'pNoH_w_dist.npy')
pNoH_w_elementHist_filename = os.path.join(model_path, 'pNoH_w_elementHist.npy')
pNoH_w_mindistInds_filename = os.path.join(model_path, 'pNoH_w_mindistInds.npy')
wcrd_filename = os.path.join(model_path, 'water_coord.npy')
data_files = [R_filename, S_filename, D_filename, displs_filename, wcrd_filename, Ddyn_filename, p_w_dists_filename, pNoH_w_dists_filename, p_w_elementHist_filename, pNoH_w_elementHist_filename, p_w_mindistInds_filename, pNoH_w_mindistInds_filename]
recomp = recomp or not files_exist(data_files)

# ======================== load data & initial process ==============================

init_traj = mdt.load(traj_filepath, stride=100000, top=init_pdb_filepath)
init_traj_pbc = mdt.load(pbctraj_filepath, stride=100000, top=init_pdb_filepath)
top = init_traj.topology
top_pbc = init_traj_pbc.topology
waterO_atoms_inds = top.select('water and name O')
pbctraj_protein_atoms_inds = top_pbc.select('protein')
pbctraj_proteinNoH_atoms_inds = top_pbc.select("protein and (element != H)")
pbctraj_waterO_atoms_inds = top_pbc.select('water and name O')
present_elements = get_present_elements(top_pbc, select='protein')
present_elements_inds = {}
for i in range(len(present_elements)):
    present_elements_inds[present_elements[i]] = i
if(verbose):
    print('top: ', top)

if(recomp):
    if(not os.path.isfile(init_pdb_filepath)):
        os.chdir(model_path)
        my.run_it(gmx_exe_name + ' trjconv -s npt.gro -f npt.xtc -skip 1000000000 -o ' + initial_pdb_filename + ' < output_whole_sys0.in')
    if(not os.path.isfile(traj_filepath)):
        os.chdir(model_path)
        my.run_it(gmx_exe_name + ' trjconv -s npt.gro -f npt.xtc -pbc nojump -o ' + traj_filename + ' < output_whole_sys0.in')

#   ========================== load & preproc data =========================
    traj = mdt.load(traj_filepath, stride=readframes_timestep, top=init_pdb_filepath)
    traj_pbc = mdt.load(pbctraj_filepath, stride=readframes_timestep, top=init_pdb_filepath)
    N_frames = traj.xyz.shape[0]
    N_atoms = traj.xyz.shape[1]
    time = np.arange(0, N_frames) * Dt
    timecut_ind = (time >= time_cut)
    time = time[timecut_ind]
    w_crd = traj.xyz[:, waterO_atoms_inds, :]
    protein_PBCcrd = traj_pbc.xyz[:, pbctraj_protein_atoms_inds, :]
    proteinNoH_PBCcrd = traj_pbc.xyz[:, pbctraj_proteinNoH_atoms_inds, :]
    water_PBCcrd = traj_pbc.xyz[:, pbctraj_waterO_atoms_inds, :]
    N_water = w_crd.shape[1]
    N_protein = protein_PBCcrd.shape[1]
    N_proteinNoH = proteinNoH_PBCcrd.shape[1]
    if(verbose):
        print('traj shape: ', traj.xyz.shape)
        print('protein coords shape: ', protein_PBCcrd.shape)
        print('proteinNoH coords shape: ', proteinNoH_PBCcrd.shape)

    if(do_waterDist):
        protein_PBCcrd_big = np.empty((N_frames, N_protein * 27, 3))
        proteinNoH_PBCcrd_big = np.empty((N_frames, N_proteinNoH * 27, 3))
        for ti in range(N_frames):
            for xi in range(3):
                for yi in range(3):
                    for zi in range(3):
                        glob_i = xi*9 + yi*3 + zi
                        cell_shift = np.array([(xi - 1) * traj_pbc.unitcell_lengths[ti, 0], \
                                               (yi - 1) * traj_pbc.unitcell_lengths[ti, 1], \
                                               (zi - 1) * traj_pbc.unitcell_lengths[ti, 2]])
                        protein_PBCcrd_big[ti, glob_i * N_protein : (glob_i + 1) * N_protein, :] = \
                            protein_PBCcrd[ti, :                                            , :] + cell_shift
                        proteinNoH_PBCcrd_big[ti, glob_i * N_proteinNoH : (glob_i + 1) * N_proteinNoH, :] = \
                            proteinNoH_PBCcrd[ti, :                                                  , :] + cell_shift
    
    
            if(verbose):
                print('bigPBC done: ' + my.f2str((ti+1)/N_frames * 100) + ' %      \r', end='')

#   =================================== heavy process ======================================
    if(do_r2range or do_D or do_Ddyn or do_R or do_S):
        all_displs = np.sum((w_crd[timecut_ind, :, :] - w_crd[0, :, :])**2, axis=2)
        with open(displs_filename, 'wb') as f:
            np.save(f, all_displs)

    if(do_D):
        D = np.empty((N_water,))
        for i in range(N_water):
            #D[i] = np.polyfit(time, all_displs[:, i], 1)
            D[i], _, _, _ = scipy.linalg.lstsq(time[:, np.newaxis], all_displs[:, i][:, np.newaxis])
        D /= 6
        with open(D_filename, 'wb') as f:
            np.save(f, D)

    if(do_Ddyn):
        N_D_dyn = (N_frames - D_timewindow - 1) // D_timestep + 1
        timewindow_ind = np.arange(D_timewindow)
        timewindow_time = timewindow_ind[:, np.newaxis] * Dt
        D_dyn = np.empty((N_D_dyn, N_water))
        for i in range(N_water):
            for j in np.arange(N_D_dyn):
                displ = (w_crd[j * D_timestep + timewindow_ind, i, :] - w_crd[j * D_timestep, i, :])**2
                D_dyn[j, i], _, _, _ = \
                    scipy.linalg.lstsq(timewindow_time, np.sum(displ, axis=1)[:, np.newaxis])

            if(verbose or True):
                print('D_dyn done: ' + my.f2str((i+1)/N_water * 100) + ' %      \r', end='')

        D_dyn /= 6
        with open(Ddyn_filename, 'wb') as f:
            np.save(f, D_dyn)

    if(do_R):
        R = np.empty((N_water,))
        R_pvalue = np.empty((N_water,))
        for i in range(N_water):
            R[i], R_pvalue[i] = scipy.stats.pearsonr(time, all_displs[:, i])
            #R[i] = my_R(time, all_displs[:, i])
        with open(R_filename, 'wb') as f:
            np.save(f, R)

    if(do_I):
        I2 = np.empty((N_water,))
        I4 = np.empty((N_water,))
        for i in range(N_water):
            I2[i] = np.log(np.sum(all_displs[:, i]**2) / (N_frames * r2_range[i]**2 / 3))
            I4[i] = np.log(np.sum(all_displs[:, i]**4) / (N_frames * r2_range[i]**4 / 5))

    if(do_S):
        S = np.empty((N_water,))
        for i in range(N_water):
            S[i] = get_S(all_displs[:, i], sgm)
        with open(S_filename, 'wb') as f:
            np.save(f, S)

    if(do_waterDist):
        p_w_dist = np.empty((N_water, N_frames))
        p_w_mindistInds = np.empty((N_water, N_frames), dtype=np.intc)
        p_w_mindistElementHist = {}
        for el in present_elements:
            p_w_mindistElementHist[el] = np.zeros(N_frames, dtype=np.intc)
        pNoH_w_dist = np.empty((N_water, N_frames))
        pNoH_w_mindistInds = np.empty((N_water, N_frames), dtype=np.intc)
        pNoH_w_mindistElementHist = {}
        for el in present_elements:
            pNoH_w_mindistElementHist[el] = np.zeros(N_frames, dtype=np.intc)

        for ti in range(N_frames):
            prot_crd = protein_PBCcrd_big[ti, :, :]
            protNoH_crd = proteinNoH_PBCcrd_big[ti, :, :]
            for i in range(N_water):
                w_crd_local = water_PBCcrd[ti, i, :]

                dists = np.sum((prot_crd - w_crd_local)**2, axis=1)
                p_w_mindistInds[i, ti] = np.argmin(dists)
                p_w_dist[i, ti] = dists[p_w_mindistInds[i, ti]]
                mindist_glob_ind = pbctraj_protein_atoms_inds[p_w_mindistInds[i, ti] % N_protein]
                p_w_mindistElementHist[top_pbc.atom(mindist_glob_ind).element.symbol][ti] += 1
        
                dists = np.sum((protNoH_crd - w_crd_local)**2, axis=1)
                pNoH_w_mindistInds[i, ti] = np.argmin(dists)
                pNoH_w_dist[i, ti] = dists[pNoH_w_mindistInds[i, ti]]
                mindist_glob_ind = pbctraj_proteinNoH_atoms_inds[pNoH_w_mindistInds[i, ti] % N_proteinNoH]
                pNoH_w_mindistElementHist[top_pbc.atom(mindist_glob_ind).element.symbol][ti] += 1

                if(verbose or True):
                    glob_i = i + ti * N_water
                    if(glob_i % 1000):
                        print('P_W_dist done: ' + my.f2str((glob_i + 1) / (N_frames * N_water) * 100) + ' %      \r', end='')

        p_w_dist = np.sqrt(p_w_dist)
        pNoH_w_dist = np.sqrt(pNoH_w_dist)
        for el in present_elements:
            p_w_mindistElementHist[el] = p_w_mindistElementHist[el] / N_water
            pNoH_w_mindistElementHist[el] = pNoH_w_mindistElementHist[el] / N_water

        with open(p_w_dists_filename, 'wb') as f:
            np.save(f, p_w_dist)
        with open(pNoH_w_dists_filename, 'wb') as f:
            np.save(f, pNoH_w_dist)
        with open(p_w_mindistInds_filename, 'wb') as f:
            np.save(f, p_w_mindistInds)
        with open(pNoH_w_mindistInds_filename, 'wb') as f:
            np.save(f, pNoH_w_mindistInds)
        with open(p_w_elementHist_filename, 'wb') as f:
            np.save(f, p_w_mindistElementHist)
        with open(pNoH_w_elementHist_filename, 'wb') as f:
            np.save(f, pNoH_w_mindistElementHist)

    with open(wcrd_filename, 'wb') as f:
        np.save(f, w_crd)
else:
    with open(wcrd_filename, 'rb') as f:
        w_crd = np.load(f)
    if(do_R):
        with open(R_filename, 'rb') as f:
            R = np.load(f)
    if(do_S):
        with open(S_filename, 'rb') as f:
            S = np.load(f)
    if(do_D):
        with open(D_filename, 'rb') as f:
            D = np.load(f)
    if(do_r2range):
        with open(displs_filename, 'rb') as f:
            all_displs = np.load(f)
    if(do_Ddyn):
        with open(Ddyn_filename, 'rb') as f:
            D_dyn = np.load(f)
    if(do_waterDist):
        with open(p_w_dists_filename, 'rb') as f:
            p_w_dist = np.load(f)
        with open(p_w_mindistInds_filename, 'rb') as f:
            p_w_mindistInds = np.load(f)
        with open(p_w_elementHist_filename, 'rb') as f:
            p_w_mindistElementHist = np.load(f, allow_pickle=True).item()
        with open(pNoH_w_dists_filename, 'rb') as f:
            pNoH_w_dist = np.load(f)
        with open(pNoH_w_mindistInds_filename, 'rb') as f:
            pNoH_w_mindistInds = np.load(f)
        with open(pNoH_w_elementHist_filename, 'rb') as f:
            pNoH_w_mindistElementHist = np.load(f, allow_pickle=True).item()

    N_frames = w_crd.shape[0]
    N_water = w_crd.shape[1]
    start_time_ind = np.ceil(time_cut / Dt)
    time = np.arange(start_time_ind, start_time_ind + N_frames) * Dt
    if(do_Ddyn):
        N_D_dyn = D_dyn.shape[0]

if(do_r2range):
    max_displ = np.max(all_displs, axis=0)
    min_displ = np.min(all_displs, axis=0)
    r2_range = max_displ - min_displ
#    if(sgm is None):
#        sgm = max(max_displ) / N_frames / 2
#    else:
#        sgm = float(sgm[0])

if(do_Ddyn):
    D_grid = np.empty((N_D_dyn, N_bins))
    time_grid = np.empty((N_D_dyn, N_bins))
    hist_grid = np.empty((N_D_dyn, N_bins))
    for i in range(N_D_dyn):
        time_grid[i, :] = time[i]
        hist_grid[i, :], bins, grid_min, grid_max = comp_hist(D_dyn[i, :], -1, 1e10, N_bins)
        D_grid[i, :] = (bins[1:] + bins[:-1]) / 2
        hist_grid[i, :] /= np.sum(hist_grid[i, :]) * abs(bins[2] - bins[1])
    
        if(verbose):
            print('postproc done: ' + my.f2str((i+1)/N_D_dyn * 100) + ' %      \r', end='')

if(do_waterDist):
    av_prot_dist_norm = np.empty((N_water,))
    filter_sigma = 2.5
    av_prot_dist = np.mean(p_w_dist, axis=1)
    for i in range(N_water):
        av_prot_dist_norm[i] = np.mean(p_w_dist[i, np.abs(p_w_dist[i, :] - np.mean(p_w_dist[i, :])) / np.std(p_w_dist[i, :]) < filter_sigma])
    is_close_to_protein = (av_prot_dist_norm < 0.225)
    pNoH_w_mindistElementHistAvg = {}
    p_w_mindistElementHistAvg = {}

    for el in present_elements:
        p_w_mindistElementHistAvg[el] = np.mean(p_w_mindistElementHist[el])
        pNoH_w_mindistElementHistAvg[el] = np.mean(pNoH_w_mindistElementHist[el])

    if(verbose):
        protDist_distrElsTime_fig, protDist_distrElsTime_ax = my.get_fig('time (ns)', '$prob$', title='dist distr($t$) in elements', yscl='log')
        for el in present_elements:
            protDist_distrElsTime_ax.plot(time, p_w_mindistElementHist[el], label=el)
        protDist_distrElsTime_ax.legend()

        protDistNoH_distrElsTime_fig, protDistNoH_distrElsTime_ax = my.get_fig('time (ns)', '$prob$', title='$dist_{no H}$ distr($t$) in elements', yscl='log')
        for el in present_elements:
            protDistNoH_distrElsTime_ax.plot(time, pNoH_w_mindistElementHist[el], label=el)
        protDistNoH_distrElsTime_ax.legend()

        protDist_distrEls_fig, protDist_distrEls_ax = my.get_fig('Element', '$prob$', title='dist distr by elements')
        protDist_distrEls_ax.bar(list(p_w_mindistElementHistAvg.keys()), p_w_mindistElementHistAvg.values(), label='whole protein', alpha=0.5)
        protDist_distrEls_ax.bar(list(pNoH_w_mindistElementHistAvg.keys()), pNoH_w_mindistElementHistAvg.values(), label='no H', alpha=0.5)
        protDist_distrEls_ax.legend()

        protDist_distrElsLog_fig, protDist_distrElsLog_ax = my.get_fig('Element', '$prob$', title='dist distr by elements', yscl='log')
        protDist_distrElsLog_ax.bar(list(p_w_mindistElementHistAvg.keys()), p_w_mindistElementHistAvg.values(), label='whole protein', alpha=0.5)
        protDist_distrElsLog_ax.bar(list(pNoH_w_mindistElementHistAvg.keys()), pNoH_w_mindistElementHistAvg.values(), label='no H', alpha=0.5)
        protDist_distrElsLog_ax.legend()
    
        protDist_distr_fig, protDist_distr_ax = my.get_fig('dist (nm)', '$p_{dens} (1/nm)$', title='$O_{water} - protein$ dist')
        protDist_distr_ax.hist(p_w_dist.flatten(), N_bins, label='whole protein', density=True, alpha=0.5)
        protDist_distr_ax.hist(pNoH_w_dist.flatten(), N_bins, label='no H', density=True, alpha=0.5)
        protDist_distr_ax.legend()
    
        #_, _, protDist_distr_ax = make_hist(pNoH_w_dist.flatten(), N_bins=N_bins, y_lbl=u'$d$ (nm)', title='$O_{water} - protein_{~H}$ distances', dens=True, scl='log')
        #protDist_distr_ax.hist(pNoH_w_dist.flatten(), N_bins, alpha=0.5, label=)

        #make_hist(av_prot_dist_norm, N_bins=N_bins, y_lbl=u'$<d>$ (nm)', title='$<d_{prot}>$ distr', dens=True)
        #make_hist(av_prot_dist, N_bins=N_bins, y_lbl=u'$<d_{reg}>$ (nm)', title='$<d_{prot, reg}>$; filter: $'+str(filter_sigma)+'\sigma$', dens=True)
        #make_hist(p_w_dist.flatten(), N_bins=N_bins, y_lbl=u'$d$ (nm)', title='$d$ distr', dens=True)
        #make_hist(p_w_dist.flatten(), N_bins=N_bins, y_lbl=u'$d$ (nm)', title='$O_{water}$ - protein distances', dens=True, scl='log')
        #_, _, Dprot_distr_loglog_ax = make_hist(p_w_dist.flatten(), N_bins=N_bins, y_lbl=u'$d$ (nm)', title='$d$ distr', dens=True, scl='log')
        #Dprot_distr_loglog_ax.set_yscale('log')

        #make_hist(protein_PBCcrd[:, :, 0].flatten(), N_bins=N_bins, y_lbl='$x$', dens=True)
        #make_hist(protein_PBCcrd[:, :, 1].flatten(), N_bins=N_bins, y_lbl='$y$', dens=True)
        #make_hist(protein_PBCcrd[:, :, 2].flatten(), N_bins=N_bins, y_lbl='$z$', dens=True)

if(do_R):
    R_min_ind = np.argmin(R)
    R_max_ind = np.argmax(R)
if(do_S):
    S_min_ind = np.argmin(S)
    S_max_ind = np.argmax(S)
if(do_I):
    I2_max_ind = np.argmax(I2)
if(do_D):
    log_D = np.log10(D / D_bulk[T])

if(do_S and do_R):
    color_from_mob(mobility_S_R_filepath, waterO_atoms_inds, S, R, N_bins)
    N_SR_cut = np.sum((R > R_cut) & (S > S_cut))
    N_SR_clust = np.sum(SR_inds)

    if(verbose):
        print('N_cut: ', N_SR_cut, '; fraction: ', N_SR_cut / N_water)
        print('N_clust: ', N_SR_clust, '; fraction: ', N_SR_clust / N_water)
    else:
        print(N_SR_cut, N_SR_clust, N_water, S_mean, R_mean)
        print(N_SR_clust / N_water)

if(do_S and do_D):
    color_from_mob(mobility_S_D_filepath, waterO_atoms_inds, S, D, N_bins)
    color_from_mob(mobility_S_logD_filepath, waterO_atoms_inds, S, log_D, N_bins)
    N_SD = np.sum(SD_inds)

if(do_D):
    D_hist, D_bins, grid_min, grid_max = comp_hist(D, -1, D_bulk[T], N_bins, dens=False)
    dens = D_hist
    max_ind = np.argmax(dens)
    D_hist_x = (D_bins[1:] + D_bins[:-1]) / 2
    D_cr = D_hist_x[max_ind + np.argmin(np.abs(dens[max_ind : ] - 0.5))]
    D_cr = np.sort(D)[int(len(D) * 0.9)]

    if(verbose):
        print(D_cr)

if(verbose):
    if(do_waterDist):
        pass
        #fig_Dprot_log, ax_Dprot_log = my.get_fig('$time$ (ns)', '$d$ (nm)', yscl='log', title='$d_{protein}(t)$')
        #fig_Dprot, ax_Dprot = my.get_fig('$time$ (ns)', '$d$ (nm)', yscl='linear', title='$d_{protein}(t)$')
        #for i in np.random.randint(N_water, size=N_bins):
        #    ax_Dprot_log.plot(time, p_w_dist[i, :])
        #    ax_Dprot.plot(time, p_w_dist[i, :])
    
    if(do_Ddyn):
        fig_relD, ax_relD = my.get_fig('$time (ns)$', '$D/D_0$', yscl='linear')
        fig_relD, ax_relD = my.get_fig('$time (ns)$', '$log_{10}(D/D_0)$')
        #for i in range(N_bins, N_bins+1):
        for i in np.random.randint(N_water, size=N_bins):
            ax_relD.plot(time[range(N_D_dyn)], D_dyn[:, i] / D_dyn[0, i])
            ax_relD.plot(time[range(N_D_dyn)], np.log10(D_dyn[:, i] / D_dyn[0, i]))
    
        fig_Ddyn, ax_Ddyn = my.get_fig('$log_{10}(D)$ ($nm^2/ns$)', 'time (ns)', projection='3d', zlbl='$p_{dens}$', title='D(t)')
        D_hist_surf = ax_Ddyn.plot_surface(np.log10(D_grid[:, 2:]), time_grid[:, 2:], hist_grid[:, 2:], cmap=matplotlib.cm.coolwarm)
        fig_Ddyn.colorbar(D_hist_surf)

    if(draw_clustering):
        if(do_S and do_R):
            SR_inds, [S_mean, R_mean], ax_SR = clustering_2D(S, R, '$S$', '$R$', gauss_cut=gauss_cut, verbose=verbose, title='$R(S)$')
            clustering_2D(S, R, '$S$', '$R$', gauss_cut=gauss_cut, verbose=verbose, title='$R(S)$', n_comp=1)

        if(do_S and do_D):
            clustering_2D(S, D, '$S$', '$D$', gauss_cut=gauss_cut, verbose=verbose, n_comp=1, title='$D(S)$')
            clustering_2D(S, log_D, '$S$', '$log_{10}(D/D_{bulk})$', gauss_cut=gauss_cut, verbose=verbose, n_comp=1, title='$D(S)$')
            clustering_2D(np.log10(S), log_D, '$log_{10}(S)$', '$log_{10}(D/D_{bulk})$', gauss_cut=gauss_cut, verbose=verbose, n_comp=1, title='$D(S)$')

        if(do_S and do_r2range):
            clustering_2D(S, np.sqrt(max_displ), '$S$', '$r_{max}$', gauss_cut=gauss_cut, verbose=verbose, n_comp=1)
        
        if(do_S and do_R and do_D):
            fig1, ax1 = my.get_fig('$SR$', '$SD$')
            ax1.scatter(S_R_mob, S_D_mob)
    
            fig2, ax2 = my.get_fig('$SR$', '$S-lgD$')
            ax2.scatter(S_R_mob, S_logD_mob)
    
            fig3, ax3 = my.get_fig('$S-lgD$', '$SD$')
            ax3.scatter(S_logD_mob, S_D_mob)

    if(draw_hists):
        if(do_r2range):
            r2_hist, r2_bins = make_hist(max_displ, r2_min, r2_max, N_bins, scl='log', y_lbl='$r_{final}^2$ ($nm^2$)')

        if(do_D):
            D_hist, D_bins, D_ax = make_hist(D, -1, D_bulk[T], N_bins, scl='linear', dens=False, \
                                       y_lbl='$D$ ($nm^2/ns$)', \
                                       title='$T = ' + my.f2str(T) + '$ $(C^{\circ}); D_{bulk} = ' + my.f2str(D_bulk[T]) + '$ $(nm^2/ns$)')
            Dlog_hist, Dlog_bins, Dlog_ax = make_hist(log_D, -6, 0, N_bins, y_lbl='$log_{10}(D/D_{bulk})$', dens=True)
            D_ax.plot([D_cr, D_cr], [0, max(D_hist)], c='red')

        if(do_R):
            R_hist, R_bins = make_hist(R, -2, 2, N_bins, scl='linear', y_lbl='R')
            R_hist_cut, R_bins_cut = make_hist(R, R_cut, 2, N_bins, scl='linear', y_lbl='R')
    
        if(do_I):
            fig_I, ax_I = my.get_fig('$I_2$', '$I_4$')
            ax_I.scatter(I2, I4, s=4)
            I2_hist, I2_bins = make_hist(I2, -20, 20, 100, scl='linear', y_lbl='$I_2$')
            I4_hist, I4_bins = make_hist(I4, -20, 20, 100, scl='linear', y_lbl='$I_4$')
        
        if(do_S):
            S_hist, S_bins, S_ax = make_hist(S, -20, 20, N_bins, y_lbl='$S$', dens=True)

            ax_SR.plot([S_cut, S_cut], [min(R), max(R)], c='red')
            ax_SR.plot([min(S), max(S)], [R_cut, R_cut], c='red')

    if(draw_extremes):
        y_lbl = '$\sqrt{r^2}$ ($nm$)'
        x_lbl = '$t$ (ns)'

        if(do_r2range):
            fig, ax = my.get_fig(x_lbl, y_lbl, title=('$I_2 = ' + my.f2str(I2[I2_max_ind]) + '$'))
            ax.plot(time, all_displs[:, I2_max_ind])

            fig_most_mobile, ax_most_mobile = my.get_fig(x_lbl, y_lbl, title='the most mobile water')
            ax_most_mobile.plot(time, np.sqrt(all_displs[:, -1]))
    
            fig_least_mobile, ax_least_mobile = my.get_fig(x_lbl, y_lbl, title='the least mobile water')
            ax_least_mobile.plot(time, np.sqrt(all_displs[:, 0]))
    
            fig_minR, ax_minR = my.get_fig(x_lbl, y_lbl, title='min R water')
            ax_minR.plot(time, np.sqrt(all_displs[:, R_min_ind]))
    
            fig_minS, ax_minS = my.get_fig(x_lbl, y_lbl, title=('max S water; $S = ' + my.f2str(S[S_max_ind]) + '$'))
            ax_minS.plot(time, np.sqrt(all_displs[:, S_max_ind]))

        if(do_S):
            draw_Rcut(all_displs, time, S, S[S_max_ind], sgm)
            draw_Rcut(all_displs, time, S, S_cut, sgm)

        if(do_D):
            draw_Rcut(all_displs, time, D, D_cut, sgm)

    if(draw_Rcuts):
        y_lbl = '$\sqrt{r^2}$ ($nm$)'
        y_lbl = '$r^2$ ($nm^2$)'
        x_lbl = '$t$ (ns)'
        for r_cut in draw_Rcuts:
            #draw_Rcut(all_displs, time, log_D, r_cut, sgm)
            draw_Rcut(all_displs, time, S, r_cut, sgm)

plt.show()
print('ALL DONE')

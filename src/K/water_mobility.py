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

import mylib as my
import mdtraj as mdt

sqrt_2pi = np.sqrt(2 * np.pi)

def files_exist(filenames):
    for f_name in filenames:
        if(not os.path.isfile(f_name)):
            return False
    return True

def comp_hist(data, mn, mx, N_bins, lin_scl=True, dens=False):
    min_data = data.min()
    max_data = data.max()
    grid_min = max(mn, min_data)
    grid_max = min(mx, max_data)
    bins = np.concatenate(([] if min_data == grid_min else [min_data], \
                           np.linspace(grid_min, grid_max, N_bins + 1) if lin_scl else np.logspace(np.log10(grid_min), np.log10(grid_max), N_bins + 1), \
                           [] if max_data == grid_max else [max_data]))

    hist, _ = np.histogram(data, bins=bins, density=dens)

    return hist, bins, grid_min, grid_max

def make_hist(data, mn, mx, N_bins, y_lbl='y', scl='linear', title=None, plot_vs_index=False, dens=False):
    lin_scl = (scl=='linear')
    hist, bins, grid_min, grid_max = comp_hist(data, mn, mx, N_bins, lin_scl=lin_scl, dens=dens)
    max_dens = (np.max(dens) if dens else np.max(hist) / np.sum(hist) / abs(bins[2] - bins[1]))

    if(plot_vs_index):
        fig, ax = my.get_fig('index', y_lbl, yscl=scl)
        ax.plot(np.arange(len(data)), data, '.')

    if(title is None):
        title = '$p(' + y_lbl + ')$'
    fig_hist, ax_hist = my.get_fig(y_lbl, '$p_{dens}$', xscl=scl, title=title)
    plt_hist, plt_bins, _ = ax_hist.hist(data, bins=bins, density=dens)

    if(grid_min < mn):
        ax_hist.plot([mn, mn], [0, max_dens], c='red')
    if(grid_min > mx):
        ax_hist.plot([mx, mx], [0, max_dens], c='red')

    return hist, bins, ax_hist

def my_R(x, y):
    x = x - np.mean(x)
    y = y - np.mean(y)
    return np.sum(np.multiply(x,y)) / np.sqrt(np.sum(x**2) * np.sum(y**2))

#@jit
def diff_S(x, data, sgm):
    x /= sgm
    dx = abs(x[1] - x[0])

    #p = np.mean(scipy.stats.norm.pdf(x, data[np.newaxis].T, sgm), axis=0)

    x_grid, data_grid = np.meshgrid(x, data)
    p = np.mean(np.exp(-(x_grid - data_grid)**2 / 2), axis=0) / (sqrt_2pi * sgm)       # density
    #p = np.mean(np.exp(-(x_grid - data_grid)**2 / 2), axis=0) * (dx / sqrt_2pi)   # prob

    p[p < 2 * np.finfo(float).eps] = 1   # p*ln(p) = 0 for such p, so we can as well put p=1 to get the same 0

    return -np.multiply(p, np.log(p))

#@jit
#def get_S(data, Nint=100):
def get_S(data, sgm, Nint=200):
    N = len(data)
    rng = max(data) - min(data)
    #sgm = rng / N
    sgm = min(rng / 2, sgm)
    x1 = min(data) - 3 * sgm
    x2 = max(data) + 3 * sgm
#    S, S_err = scipy.integrate.quad(diff_S, x1, x2, args=(data, sgm))

    x_draw = np.linspace(x1, x2, Nint)
    #p = np.mean(scipy.stats.norm.pdf(x_draw, data[np.newaxis].T, sgm), axis=0)
    S = np.mean(diff_S(x_draw, data, sgm)) * rng  # density
    #S = np.sum(diff_S(x_draw, data, sgm))          # prob

    return S

def draw_Rcut(all_displs, time, x, x_cut, sgm, Ndgt=5):
    ind = np.argmin(np.abs(x - x_cut))
    data = all_displs[:, ind]
    N = len(data)
    rng = max(data) - min(data)
    R, R_pvalue = scipy.stats.pearsonr(time, data)

    fig, ax = my.get_fig('$t$ (ns)', '$\sqrt{r^2}$ ($nm$)', title=r'$i = ' + str(ind) + r'; R_{r^2} \approx ' + my.f2str(R, Ndgt) + '$')
    ax.plot(time, np.sqrt(data))
    #ax.plot(time, data)

    #sgm = rng / N
    N_draw = 1000
    x_draw = np.linspace(min(data) - 3*sgm, max(data) + 3*sgm, N_draw)
    distr = np.mean(scipy.stats.norm.pdf(x_draw, data[np.newaxis].T, sgm), axis=0)
    #S, S_err = scipy.integrate.quad(diff_S, min(x_draw), max(x_draw), args=(data, sgm))
    S = get_S(data, sgm)

    fig, ax = my.get_fig('$r^2$', 'n', title=r'$i = ' + str(ind) + r'; S \approx ' + my.f2str(S, Ndgt) + '$')
    ax.plot(x_draw, distr)

def clustering_2D(x, y, x_lbl, y_lbl, n_comp=2, gauss_cut=1, N_X_grid=300, N_Y_grid=300, verbose=False):
    ax = None
    if(n_comp == 1):
        if(verbose):
            fig, ax = my.get_fig(x_lbl, y_lbl)
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
        free_water_inds = (gauss_dists < gauss_cut) & free_water_labels
        
        if(verbose):
            fig, ax = my.get_fig(x_lbl, y_lbl)
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

# =============== paths ================
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
topol_filename = 'topol.top'
initial_pdb_filename = 'initial_npt.pdb'

r2_min = 7e0
r2_max = 1e3
R_cut = 0.75
S_cut = 0.5

# ============== arg parse ==========================
supercell_str = ''.join([str(x) for x in supercell])

[T, draw_extremes, draw_hists, time_cut, draw_Rcuts, verbose, recomp, sgm, gauss_cut, D_timewindow, D_timestep, N_bins], _ = \
    my.parse_args(sys.argv[1:], ['-temp', '-extremes', '-hists', '-time_cut', '-Rcut', '-verbose', '-recomp', '-sgm', '-gauss_cut', '-D_timewindow', '-D_timestep', '-N_bins'], \
                  possible_values=[None, ['0', '1'], ['0', '1'], None, None, ['0', '1'], ['0', '1'], None, None, None, None, None], \
                  possible_arg_numbers=[[1], [0, 1], [0, 1], [0, 1], None, [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1]], \
                  default_values=[None, ['0'], ['0'], ['0'], [], ['0'], ['0'], ['6.0'], ['5.0'], ['2'], ['1'], ['100']])
T = float(T)
sgm = float(sgm[0])
N_bins = int(N_bins[0])
time_cut = float(time_cut[0])
gauss_cut = float(gauss_cut[0])
draw_extremes = (draw_extremes[0] == '1')
draw_hists = (draw_hists[0] == '1')
verbose = (verbose[0] == '1')
recomp = (recomp[0] == '1')
D_timewindow = int(D_timewindow[0])
D_timestep = int(D_timestep[0])
draw_Rcuts = [float(r) for r in draw_Rcuts]

model_name = 'flucts_temp' + my.f2str(T) + '_0'
model_path = os.path.join(run_path, model_name)
traj_filepath = os.path.join(model_path, traj_filename)
topol_filepath = os.path.join(model_path, topol_filename)
init_pdb_filepath = os.path.join(model_path, initial_pdb_filename)
R_filename = os.path.join(model_path, 'R_time' + my.f2str(time_cut) + '.npy')
D_filename = os.path.join(model_path, 'D_time' + my.f2str(time_cut) + '.npy')
Ddyn_filename = os.path.join(model_path, 'Ddyn_W' + my.f2str(D_timewindow) + '_s' + my.f2str(D_timestep) + '.npy')
S_filename = os.path.join(model_path, 'S_time' + my.f2str(time_cut) + '_sgm' + my.f2str(sgm) + '.npy')
displs_filename = os.path.join(model_path, 'displs_time' + my.f2str(time_cut) + '.npy')
maxdispl_filename = os.path.join(model_path, 'maxdispl_time' + my.f2str(time_cut) + '.npy')
wcrd_filename = os.path.join(model_path, 'water_coord.npy')
data_files = [R_filename, S_filename, D_filename, maxdispl_filename, displs_filename, wcrd_filename, Ddyn_filename]
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
    N_D_dyn = (N_frames - D_timewindow - 1) // D_timestep + 1
    timewindow_ind = np.arange(D_timewindow)
    timewindow_time = timewindow_ind[:, np.newaxis] * Dt
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
    D_dyn = np.empty((N_D_dyn, N_water))
    R = np.empty((N_water,))
    R_pvalue = np.empty((N_water,))
    I2 = np.empty((N_water,))
    I4 = np.empty((N_water,))
    S = np.empty((N_water,))
    for i in range(N_water):
        #linfit = np.polyfit(time, all_displs[:, i], 1)
        linfit, _, _, _ = scipy.linalg.lstsq(time[:, np.newaxis], all_displs[:, i][:, np.newaxis])
        D[i] = linfit / 6
        #R[i], R_pvalue[i] = scipy.stats.pearsonr(time, all_displs[:, i])
        #R[i] = my_R(time, all_displs[:, i])
        #I2[i] = np.log(np.sum(all_displs[:, i]**2) / (N_frames * r2_range[i]**2 / 3))
        #I4[i] = np.log(np.sum(all_displs[:, i]**4) / (N_frames * r2_range[i]**4 / 5))
        #S[i] = get_S(all_displs[:, i], sgm)

        #for j in np.arange(N_D_dyn):
        #    displ = (w_crd[j * D_timestep + timewindow_ind, i, :] - w_crd[j * D_timestep, i, :])**2
        #    D_dyn[j, i], _, _, _ = \
        #        scipy.linalg.lstsq(timewindow_time, np.sum(displ, axis=1)[:, np.newaxis])

        if(verbose):
            if(not i%100):
                print('main done: ' + my.f2str((i+1)/N_water * 100) + ' %      \r', end='')
    D_dyn /= 6

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
    with open(wcrd_filename, 'wb') as f:
        np.save(f, w_crd)
#    with open(Ddyn_filename, 'wb') as f:
#        np.save(f, D_dyn)

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
    with open(wcrd_filename, 'rb') as f:
        w_crd = np.load(f)
    with open(Ddyn_filename, 'rb') as f:
        D_dyn = np.load(f)
    N_water = len(R)
    N_frames = all_displs.shape[0]
    start_time_ind = np.ceil(time_cut / Dt)
    time = np.arange(start_time_ind, start_time_ind + N_frames) * Dt
    N_D_dyn = D_dyn.shape[0]
    #N_water = D_dyn.shape[1]

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
    
#R_min_ind = np.argmin(R)
#R_max_ind = np.argmax(R)
#S_min_ind = np.argmin(S)
#S_max_ind = np.argmax(S)
#I2_max_ind = np.argmax(I2)

#SR_inds, [S_mean, R_mean], ax_SR = clustering_2D(S, R, '$S$', '$R$', gauss_cut=gauss_cut, verbose=verbose)
#ax_SD = clustering_2D(S, D, '$S$', '$D$', gauss_cut=gauss_cut, verbose=verbose, n_comp=1)
#ax_Sr = clustering_2D(S, np.sqrt(max_displ), '$S$', '$r_{max}$', gauss_cut=gauss_cut, verbose=verbose, n_comp=1)

#N_SR_cut = np.sum((R > R_cut) & (S > S_cut))
#N_SR_clust = np.sum(SR_inds)
#N_SD = np.sum(SD_inds)
#if(verbose):
#    print('N_cut: ', N_SR_cut, '; fraction: ', N_SR_cut / N_water)
#    print('N_clust: ', N_SR_clust, '; fraction: ', N_SR_clust / N_water)
#else:
#    print(N_SR_cut, N_SR_clust, N_water, S_mean, R_mean)
    #print(N_SR_clust / N_water)

D_hist, D_bins, grid_min, grid_max = comp_hist(D, -1, D_bulk[T], N_bins, dens=True)
dens = D_hist
max_ind = np.argmax(dens)
D_hist_x = (D_bins[1:] + D_bins[:-1]) / 2
D_cr = D_hist_x[max_ind + np.argmin(np.abs(dens[max_ind : ] - 0.5))]
#print(T, ' : ', D_cr)
print(D_cr)

if(verbose):
    #fig_relD, ax_relD = my.get_fig('$time (ns)$', '$D/D_0$', yscl='linear')
    #fig_relD, ax_relD = my.get_fig('$time (ns)$', '$log_{10}(D/D_0)$')
    #for i in range(N_bins, N_bins+1):
    #for i in np.random.randint(N_water, size=N_bins):
        #ax_relD.plot(time[range(N_D_dyn)], D_dyn[:, i] / D_dyn[0, i])
        #ax_relD.plot(time[range(N_D_dyn)], np.log10(D_dyn[:, i] / D_dyn[0, i]))

#    fig_Ddyn, ax_Ddyn = my.get_fig('$log_{10}(D)$ ($nm^2/ns$)', 'time (ns)', projection='3d', zlbl='$p_{dens}$', title='D(t)')
#    D_hist_surf = ax_Ddyn.plot_surface(np.log10(D_grid[:, 2:]), time_grid[:, 2:], hist_grid[:, 2:], cmap=matplotlib.cm.coolwarm)
#    fig_Ddyn.colorbar(D_hist_surf)

    if(draw_hists):
        #r2_hist, r2_bins = make_hist(max_displ, r2_min, r2_max, N_bins, scl='log', y_lbl='$r_{final}^2$ ($nm^2$)')
        D_hist, D_bins, D_ax = make_hist(D, -1, D_bulk[T], N_bins, scl='linear', dens=True, \
                                   y_lbl='$D$ ($nm^2/ns$)', \
                                   title='$T = ' + my.f2str(T) + '$ $(C^{\circ}); D_{bulk} = ' + my.f2str(D_bulk[T]) + '$ $(nm^2/ns$)')
        #R_hist, R_bins = make_hist(R, -2, 2, N_bins, scl='linear', y_lbl='R')
        #R_hist_cut, R_bins_cut = make_hist(R, R_cut, 2, N_bins, scl='linear', y_lbl='R')

        D_ax.plot([D_cr, D_cr], [0, max(dens)], c='red')
    
        #fig_I, ax_I = my.get_fig('$I_2$', '$I_4$')
        #ax_I.scatter(I2, I4, s=4)
        #I2_hist, I2_bins = make_hist(I2, -20, 20, 100, scl='linear', y_lbl='$I_2$')
        #I4_hist, I4_bins = make_hist(I4, -20, 20, 100, scl='linear', y_lbl='$I_4$')
    
        #S_hist, S_bins = make_hist(S, -20, 20, 100, scl='linear', y_lbl='$S$')

#        ax_SR.plot([S_cut, S_cut], [min(R), max(R)], c='red')
#        ax_SR.plot([min(S), max(S)], [R_cut, R_cut], c='red')
    
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

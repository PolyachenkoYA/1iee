import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm as colormaps
from numba import jit
from matplotlib.colors import LogNorm
import dill
import scipy as sci

import mylib as my
import mdtraj as mdt
import gromacs.formats as gmx

def rollavg_convolve_edges(a, n):
    assert n % 2 == 1
    #return sci.convolve(a, np.ones(n, dtype='float'), 'same') / sci.convolve( np.ones(len(a)), np.ones(n), 'same')
    return np.convolve(a, np.ones(n, dtype='float'), 'same') / np.convolve( np.ones(len(a)), np.ones(n), 'same')

def files_exist(filenames):
    for f_name in filenames:
        if(not os.path.isfile(f_name)):
            return False
    return True

def load_xvg(filepath, failsafe=2):
    if(failsafe in [1, 2]):
        if(not os.path.isfile(filepath)):
            print('file "' + filepath + '" was not found. Aborting.')
            if(failsafe == 1):
                return None
            elif(failsafe == 2):
                sys.exit(1)
    xvg_file = gmx.XVG()
    xvg_file.read(filepath)
    return xvg_file
    
def parse_xvg(xvg_data, fields):
    N_frames = len(xvg_data.array[0, :])
    N_f = len(fields)
    res = np.zeros((N_f, N_frames))
    for fi in range(N_f):
        ind = np.where([name == fields[fi] for name in xvg_data.names])[0]
        if(ind.size == 1):
            res[fi, :] = xvg_data.array[ind[0] + 1, :]   # +1 because 0th column is time and it's not in fields
        else:
            print('field "' + fields[fi] + '" not found of found more than once (ind = ' + str(ind) + ')')

    return res


# =============== global paths ================
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'K')
res_path = os.path.join(root_path, 'res')

# ============== arg parse =================

[Tmp, recomp, Zstep, avg_hist_time, extra_water, to_save_pics, to_draw_pics, comprZ], _ = \
    my.parse_args(sys.argv[1:], ['-temp', '-recomp', '-Zstep', '-avg_time', '-extra_water', '-save_pics', '-draw_pics', '-comprZ'], \
                  possible_values=[None, ['0', '1'], None, None, None, ['0', '1'], ['0', '1'], None], \
                  possible_arg_numbers=[[1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1]], \
                  default_values=[None, ['0'], ['0.5'], ['2.5'], ['0'], ['0'], ['1'], ['10']])
Tmp = float(Tmp)
recomp = (recomp[0] == '1')
Zstep = float(Zstep[0])
avg_hist_time = float(avg_hist_time[0])
extra_water = int(extra_water[0])
to_save_pics = (to_save_pics[0] == '1')
to_draw_pics = (to_draw_pics[0] == '1') or to_save_pics
comprZ = float(comprZ[0])

# ================ params =================
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns,      dt of the numerica scheme
Dt = dt * 5000     #                Dt of dumps
m_H2O = 0.018 / (6.02 * 10**23)
gas_max_rho = 0.1
readframes_timestep = 1
supercell = np.array([1, 1, 2], dtype=np.intc)
verbose = True
to_draw_dist = to_draw_pics and False
to_draw_P = to_draw_pics and True
to_draw_L = to_draw_pics and True

step_hist_av = int(round(avg_hist_time / Dt))
rho_atm = 101000 * 0.018 / (8.31 * (Tmp + T_C2K))
supercell_str = ''.join([str(x) for x in supercell])

# =================== per-model paths ======================
model_name = 'flucts_t4p2005_temp' + my.f2str(Tmp) + '_extW' + str(extra_water) + '_comprZ' + str(comprZ)
model_path = os.path.join(run_path, model_name)
traj_filepath = os.path.join(model_path, 'npt.xtc')
xvg_filepath = os.path.join(model_path, 'npt.xvg')
topol_filepath = os.path.join(model_path, 'topol.top')
init_pdb_filepath = os.path.join(model_path, 'em_nojump.pdb')
wcrd_filepath = os.path.join(model_path, 'water_coord.npy')
xvg_bin_filepath = os.path.join(model_path, 'xvg.pkl')
init_pdb_bin_filepath = os.path.join(model_path, 'initPDB.pkl')
_filepath = os.path.join(model_path, 'npt.xtc')

data_files = [xvg_bin_filepath, init_pdb_bin_filepath, wcrd_filepath]
recomp = recomp or not files_exist(data_files)

#   ========================== load & preproc data =========================
if(os.path.isfile(xvg_bin_filepath) and not True):  # pickle breaks xvg files
    print('loading "' + xvg_bin_filepath + '"')
    with open(xvg_bin_filepath, 'rb') as f:
        xvg_file = dill.load(f)
else:
    print('reading "' + xvg_filepath + '"')
    xvg_file = load_xvg(xvg_filepath)
    with open(xvg_bin_filepath, 'wb') as f:
        dill.dump(xvg_file, f)        
# ===================================
if(os.path.isfile(init_pdb_bin_filepath) and not recomp):
    print('loading "' + init_pdb_bin_filepath + '"')
    with open(init_pdb_bin_filepath, 'rb') as f:
        init_traj = dill.load(f)        
else:
    print('reading "' + init_pdb_filepath + '"')
    init_traj = mdt.load(traj_filepath, stride=100000, top=init_pdb_filepath)
    with open(init_pdb_bin_filepath, 'wb') as f:
        dill.dump(init_traj, f)        
top = init_traj.topology
# =================================
if(os.path.isfile(wcrd_filepath) and not recomp):
    print('loading "' + wcrd_filepath + '"')
    with open(wcrd_filepath, 'rb') as f:
        w_crd = np.load(f)
else:
    print('reading "' + traj_filepath + '"')
    traj = mdt.load(traj_filepath, stride=readframes_timestep, top=init_pdb_filepath)
    
    waterO_atoms_inds = init_traj.topology.select('water and name O')
    w_crd = traj.xyz[:, waterO_atoms_inds, :]
    
    with open(wcrd_filepath, 'wb') as f:
        np.save(f, w_crd)        

    if(verbose):
        print('traj shape: ', traj.xyz.shape)

#start_time_ind = np.ceil(time_cut / Dt)
#time = np.arange(start_time_ind, start_time_ind + N_frames) * Dt

N_frames = w_crd.shape[0]
N_water = w_crd.shape[1]
time = np.arange(0, N_frames) * Dt
if(verbose):
    print('top: ', top)
    print('w_crd shape: ', w_crd.shape)

# =================== process ======================
Lbox = parse_xvg(xvg_file, ['Box-X', 'Box-Y', 'Box-Z'])
Pbox = parse_xvg(xvg_file, ['Pres-XX', 'Pres-YY', 'Pres-ZZ'])
time_xvg = xvg_file.array[0] * 1e-3
dt_xvg = abs(time_xvg[1] - time_xvg[0])

w_z = w_crd[:, :, 2]
Zmin = np.min(w_z[:])
Zmax = np.max(w_z[:])
L = Zmax - Zmin
N_Zbins = int(L / Zstep)
Zstep = L / N_Zbins
Zbins = np.arange(N_Zbins) * Zstep + Zmin
Zbins = np.append(Zbins, Zmax)
Zcentrs = (Zbins[1:] + Zbins[:-1]) / 2
water_hist = np.zeros((N_frames, N_Zbins))
d_water_hist = np.zeros((N_frames, N_Zbins))
for ti in range(N_frames):
    xvg_closest_time_ind = np.argmin(np.abs(time_xvg - time[ti]))
    dV = (Zbins[1:] - Zbins[:-1]) * Lbox[0, xvg_closest_time_ind] * Lbox[1, xvg_closest_time_ind] * 10**(-27)
    water_hist[ti, :] = np.histogram(w_z[ti, :], bins=Zbins)[0] * (m_H2O / dV)
    d_water_hist[ti, :] = np.sqrt(water_hist[ti, :] * (1 - water_hist[ti, :] / N_water)) * (m_H2O / dV)

N_av = int(N_frames / step_hist_av)
water_hist_average = np.zeros((N_av, N_Zbins))
d_water_hist_average = np.zeros((N_av, N_Zbins))
hist_av_inds = np.zeros((N_av, 2), dtype=np.intc)
for ti in range(N_av - 1):
    hist_av_inds[ti, 0] = step_hist_av * ti
    hist_av_inds[ti, 1] = step_hist_av * (ti + 1)
hist_av_inds[N_av - 1, 0] = step_hist_av * (N_av - 1)
hist_av_inds[N_av - 1, 1] = N_frames
gas_rho = np.zeros(N_av)
d_gas_rho = np.zeros(N_av)
for ti in range(N_av):
    water_hist_average[ti, :] = np.mean(water_hist[hist_av_inds[ti, 0] : hist_av_inds[ti, 1], :], axis=0)
    #d_water_hist_average[ti, :] = np.std(water_hist[hist_av_inds[ti, 0] : hist_av_inds[ti, 1], :], axis=0) / np.sqrt(hist_av_inds[ti, 1] - hist_av_inds[ti, 0])
    d_water_hist_average[ti, :] = \
        np.sum(d_water_hist[hist_av_inds[ti, 0] : hist_av_inds[ti, 1], :] ** 2, axis=0) / (hist_av_inds[ti, 1] - hist_av_inds[ti, 0])
    gas_Z_ind = water_hist_average[ti, :] < gas_max_rho
    gas_rho[ti] = np.mean(water_hist_average[ti, gas_Z_ind])
    d_gas_rho[ti] = np.sum(d_water_hist_average[ti, gas_Z_ind] ** 2) / np.sum(gas_Z_ind)

if(to_draw_pics):
    fig, ax = my.get_fig(r'$t$ (ns)', r'$\rho_w$ (kg / $m^3$)', title=r'$\rho_w(t)$; $T$ = ' + my.f2str(Tmp) + r' $C^\circ$; $\rho_{sat} = ' + my.f2str(rho_atm) + '$ kg / $m^3$')
    t_draw = (time[hist_av_inds[:, 0]] + time[hist_av_inds[:, 1] - 1]) / 2
    ax.errorbar(t_draw, gas_rho, d_gas_rho, fmt='o', markerfacecolor='None')
    
    if(to_save_pics):
        pic_rho_t_filepath = os.path.join(model_path, 'rho(t).png')
        fig.savefig(pic_rho_t_filepath)

if(to_draw_L):
    fig_z, ax_z = my.get_fig(r'$t$ (ns)', r'$L_z$ (nm)', title=r'$L_z(t)$; $\beta = ' + my.f2str(comprZ) + ' (1/bar)$')
    ax_z.plot(time_xvg, Lbox[2, :])
    fig_xy, ax_xy = my.get_fig(r'$t$ (ns)', r'$L$ (nm)', title=r'$L(t)$; $T = ' + my.f2str(comprZ) + ' (1/bar)$')
    ax_xy.plot(time_xvg, Lbox[0, :], label=r'$L_x$')
    ax_xy.plot(time_xvg, Lbox[1, :], label=r'$L_y$')
    ax_xy.legend()
    
    if(to_save_pics):
        pic_Lz_t_filepath = os.path.join(model_path, 'Lz(t).png')
        pic_Lxy_t_filepath = os.path.join(model_path, 'Lxy(t).png')
        fig_z.savefig(pic_Lz_t_filepath)
        fig_xy.savefig(pic_Lxy_t_filepath)

if(to_draw_P):
    fig_Pz, ax_Pz = my.get_fig(r'$t$ (ns)', r'$P_z$ (nm)', title=r'$P_z(t)$; $T = ' + my.f2str(Tmp) + ' C^\circ$')
    ax_Pz.plot(time_xvg, Pbox[2, :], color = my.colors[2], alpha=0.1, label='no avg')
    ax_Pz.plot(time_xvg, rollavg_convolve_edges(Pbox[2, :], 51), color = my.colors[2], alpha=0.5, label='avg = ' + my.f2str(dt_xvg * 51) + ' ns')
    ax_Pz.plot(time_xvg, rollavg_convolve_edges(Pbox[2, :], 251), color = my.colors[2], label='avg = ' + my.f2str(dt_xvg * 251) + ' ns')

    fig_Pxy, ax_Pxy = my.get_fig(r'$t$ (ns)', r'$P$ (bar)', title=r'$P(t)$; $T = ' + my.f2str(Tmp) + ' C^\circ$')
    ax_Pxy.plot(time_xvg, Pbox[1, :], color = my.colors[1], alpha=0.1, label='no avg')
    ax_Pxy.plot(time_xvg, rollavg_convolve_edges(Pbox[1, :], 51), color = my.colors[1], alpha=0.5, label='avg = ' + my.f2str(dt_xvg * 51) + ' ns')
    ax_Pxy.plot(time_xvg, rollavg_convolve_edges(Pbox[1, :], 251), color = my.colors[1], label=r'$P_y$; avg = ' + my.f2str(dt_xvg * 251) + ' ns')
    ax_Pxy.plot(time_xvg, Pbox[0, :], color = my.colors[0], alpha=0.1, label='no avg')
    ax_Pxy.plot(time_xvg, rollavg_convolve_edges(Pbox[0, :], 51), color = my.colors[0], alpha=0.5, label='avg = ' + my.f2str(dt_xvg * 51) + ' ns')
    ax_Pxy.plot(time_xvg, rollavg_convolve_edges(Pbox[0, :], 251), color = my.colors[0], label=r'$P_y$; avg = ' + my.f2str(dt_xvg * 251) + ' ns')
    ax_Pxy.legend()
    
    fig_Pav, ax_Pav = my.get_fig(r'$t$ (ns)', r'$P$ (bar)', title=r'$P(t)$; $T = ' + my.f2str(Tmp) + ' C^\circ$')
    ax_Pav.plot(time_xvg, rollavg_convolve_edges(Pbox[0, :], 251), color = my.colors[0], label=r'$P_x$; avg = ' + my.f2str(dt_xvg * 251) + ' ns')
    ax_Pav.plot(time_xvg, rollavg_convolve_edges(Pbox[1, :], 251), color = my.colors[1], label=r'$P_y$; avg = ' + my.f2str(dt_xvg * 251) + ' ns')
    ax_Pav.plot(time_xvg, rollavg_convolve_edges(Pbox[2, :], 251), color = my.colors[2], label=r'$P_z$; avg = ' + my.f2str(dt_xvg * 251) + ' ns')
    ax_Pav.legend()

    if(to_save_pics):
        pic_Pz_t_filepath = os.path.join(model_path, 'Pz(t).png')
        pic_Pxy_t_filepath = os.path.join(model_path, 'Pxy(t).png')
        pic_Pav_t_filepath = os.path.join(model_path, 'Pav(t).png')
        fig_Pz.savefig(pic_Pz_t_filepath)
        fig_Pxy.savefig(pic_Pxy_t_filepath)
        fig_Pav.savefig(pic_Pav_t_filepath)

if(to_draw_dist):
    fig_dist, ax_dist = my.get_fig(r'z (nm)', r'$\rho_w$ (kg / $m^3$)', title=r'$\rho_w(z)$', yscl='log')
    plt.ion()
    plt.show()
    for ti in range(N_av):
        ax_dist.clear()
        plt.yscale('log')
        ax_dist.bar(Zcentrs, water_hist_average[ti, :], width=Zstep, \
               facecolor=my.colors[4], edgecolor=my.colors[4])
        ax_dist.plot([min(Zcentrs), max(Zcentrs)], [rho_atm] * 2, color=my.colors[1], label=r'$x$')
        #ax.legend()
        fig_dist.suptitle(r'$t$ \in $[' + my.f2str(time[hist_av_inds[ti, 0]]) + ' - ' + my.f2str(time[hist_av_inds[ti, 1] - 1]) + ']$ (ns)')
        plt.draw()
        plt.pause(0.001)
        input("Press Enter to continue...")
else:
    plt.show()

print('ALL DONE')

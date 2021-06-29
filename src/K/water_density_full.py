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
        for fj in range(len(fields[fi])):
            ind = np.where([name == fields[fi][fj] for name in xvg_data.names])[0]
            if(ind.size == 1):
                break
            
        if(ind.size == 1):
            res[fi, :] = xvg_data.array[ind[0] + 1, :]   # +1 because 0th column is time and it's not in fields
        else:
            print('field "' + fields[fi] + '" not found of found more than once (ind = ' + str(ind) + ')')

    return res


# =============== global paths ================
#root_path = my.git_root_path()
with open('git_root_path', 'r') as f:
    root_path = f.readline()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'K')
res_path = os.path.join(root_path, 'res')

# ============== arg parse =================

[Tmp, recomp, Zstep, avg_hist_time, extra_water, to_save_pics, to_draw_pics, comprZ, model_id, do_dist, do_P, do_L], _ = \
    my.parse_args(sys.argv[1:], ['-temp', '-recomp', '-Zstep', '-avg_time', '-extra_water', '-save_pics', '-draw_pics', '-comprZ', '-id', '-do_dist', '-do_P', '-do_L'], \
                  possible_values=[None, ['0', '1'], None, None, None, ['0', '1'], ['0', '1'], None, None, ['0', '1'], ['0', '1'], ['0', '1']], \
                  possible_arg_numbers=[[1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1]], \
                  default_values=[None, ['0'], ['0.5'], ['5'], ['0'], ['0'], ['1'], ['0'], ['0'], ['1'], ['1'], ['1']])
Tmp = float(Tmp)
recomp = (recomp[0] == '1')
Zstep = float(Zstep[0])
avg_hist_time = float(avg_hist_time[0])
extra_water = int(extra_water[0])
to_save_pics = (to_save_pics[0] == '1')
to_draw_pics = (to_draw_pics[0] == '1') or to_save_pics
comprZ = float(comprZ[0])
model_id = int(model_id[0])
do_dist = (do_dist[0] == '1')
do_P = (do_P[0] == '1')
do_L = (do_L[0] == '1')

# ================ params =================
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns,      dt of the numerica scheme
Dt = dt * 5000     #                Dt of dumps
m_H2O = 0.018 / (6.02 * 10**23)
gas_max_rho = 0.1
readframes_timestep = 1
supercell = np.array([1, 1, 2], dtype=np.intc)
verbose = True
clr_x = my.colors[0]
clr_y = my.colors[3]
clr_z = my.colors[2]

step_hist_av = int(round(avg_hist_time / Dt / 2)) * 2 + 1
rho_atm = 101000 * 0.02 * 0.018 / (8.31 * (Tmp + T_C2K))
sat_params = [27.373, 5973.2, 1.002]
rho_atm = np.exp(sat_params[0] - sat_params[1] / (Tmp + T_C2K)) / (Tmp + T_C2K)**sat_params[2] / 1000
supercell_str = ''.join([str(x) for x in supercell])

# =================== per-model paths ======================
#model_name = 'flucts_t4p2005_temp' + my.f2s(Tmp) + '_extW' + str(extra_water) + '_AnisoXYComprZ' + str(comprZ) + '_id' + str(model_id)
#model_name = 'flucts_t4p2005_temp' + my.f2s(Tmp) + '_extW' + str(extra_water) + '_comprZ' + str(comprZ)
model_name = 'flucts_t4p2005_temp' + my.f2s(Tmp) + '_extW' + str(extra_water)
model_path = os.path.join(run_path, model_name)
traj_filepath = os.path.join(model_path, 'npt.xtc')
xvg_filepath = os.path.join(model_path, 'npt.xvg')
topol_filepath = os.path.join(model_path, 'topol.top')
init_pdb_filepath = os.path.join(model_path, 'em_nojump.pdb')
wcrd_filepath = os.path.join(model_path, 'water_coord.npy')
xvg_bin_filepath = os.path.join(model_path, 'xvg.pkl')
init_pdb_bin_filepath = os.path.join(model_path, 'initPDB.pkl')
#_filepath = os.path.join(model_path, 'npt.xtc')

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
Lbox = parse_xvg(xvg_file, [['Box-X', 'Box-XX'], ['Box-Y', 'Box-YY'], ['Box-Z', 'Box-ZZ']])
Pbox = parse_xvg(xvg_file, [['Pres-XX'], ['Pres-YY'], ['Pres-ZZ']])
Ptbox = parse_xvg(xvg_file, [['Pres-XY', 'Pres-YX'], ['Pres-XZ', 'Pres-ZX'], ['Pres-YZ', 'Pres-ZY']])
time_xvg = xvg_file.array[0] * 1e-3
dt_xvg = abs(time_xvg[1] - time_xvg[0])
N_xvg = Lbox.shape[1]

w_z = w_crd[:, :, 2]
Zmin = np.min(w_z[:])
#Zmin = 0
Zmax = np.max(w_z[:])
#Zmax = np.max(Lbox[2, :])
L = Zmax - Zmin
N_Zbins = int(np.ceil(L / Zstep))
Zstep = L / N_Zbins
Zbins = np.arange(N_Zbins) * Zstep + Zmin
Zbins = np.append(Zbins, Zmax)
Zcenters = (Zbins[1:] + Zbins[:-1]) / 2
water_hist = np.zeros((N_frames, N_Zbins))
d_water_hist = np.zeros((N_frames, N_Zbins))
for ti in range(N_frames):
    xvg_closest_time_ind = np.argmin(np.abs(time_xvg - time[ti]))
    dV = (Zbins[1:] - Zbins[:-1]) * Lbox[0, xvg_closest_time_ind] * Lbox[1, xvg_closest_time_ind] * 10**(-27)
    water_hist[ti, :] = np.histogram(w_z[ti, :], bins=Zbins)[0] * (m_H2O / dV)
    d_water_hist[ti, :] = np.sqrt(water_hist[ti, :] * (1 - water_hist[ti, :] / N_water)) * (m_H2O / dV)
    
z1_gas = 8   # left side of the layer
z2_gas = 25   # right side of the layer
t_rho_stab = 30
gas_z_ind = np.logical_or(Zcenters < z1_gas - Zstep, Zcenters > z2_gas + Zstep)
gas_z_grid_ind, t_stab_grid_ind = np.meshgrid(gas_z_ind, time > t_rho_stab)
w_stab_gas_ind = np.logical_and(t_stab_grid_ind, gas_z_grid_ind)
rho_mean = np.mean(water_hist[w_stab_gas_ind])
d_rho_mean = np.sqrt(np.sum(d_water_hist[w_stab_gas_ind] ** 2)) / np.sum(w_stab_gas_ind)

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
        np.sqrt(np.sum(d_water_hist[hist_av_inds[ti, 0] : hist_av_inds[ti, 1], :] ** 2, axis=0)) / (hist_av_inds[ti, 1] - hist_av_inds[ti, 0])
    #gas_rho_ind = water_hist_average[ti, :] < gas_max_rho
    gas_rho_ind = np.logical_or(Zcenters < z1_gas - Zstep, Zcenters > z2_gas + Zstep)
    gas_rho[ti] = np.mean(water_hist_average[ti, gas_rho_ind])
    d_gas_rho[ti] = np.sqrt(np.sum(d_water_hist_average[ti, gas_rho_ind] ** 2)) / np.sum(gas_rho_ind)

title_suff = 'T = ' + my.f2s(Tmp) + r' $C^{\circ}$; $H_2 O_{ext}$ = ' + str(extra_water) + ' ($H_2 O$/cell)'

if(to_draw_pics):
    def reg_mean(y, n=1, sgm=1):
        for i in range(n):
            y = y[np.abs(y - np.mean(y)) < np.std(y) * sgm]
        return np.mean(y)

    draw_units = 1e3   # kg to g
    fig, ax = my.get_fig(r'$t$ (ns)', r'$\rho_w$ (g / $m^3$)', \
                        title=r'$\rho_w(t)$; ' + title_suff)
    t_draw = (time[hist_av_inds[:, 0]] + time[hist_av_inds[:, 1] - 1]) / 2
    gas_rho_stab = my.regularize(gas_rho[t_draw > t_rho_stab], n=0, sgm=1)
    rho_mean_stab = np.mean(gas_rho_stab)
    d_rho_mean_stab = np.std(gas_rho_stab) / np.sqrt(len(gas_rho_stab))
    ax.errorbar(t_draw, gas_rho * draw_units, d_gas_rho * draw_units, fmt='o', markerfacecolor='None', label='data')
    ax.plot([min(t_draw), max(t_draw)], [rho_atm * draw_units] * 2, '--', label=r'$\rho_{sat} = ' + my.f2s(rho_atm * draw_units) + '$')
    #ax.plot([min(t_draw), max(t_draw)], [rho_mean_stab * draw_units] * 2, '--', label=r'$\rho_{reg} = ' + my.f2s(rho_mean_stab * draw_units) + ' \pm ' + my.f2s(d_rho_mean_stab * draw_units) + '$')
    ax.plot([min(t_draw), max(t_draw)], [rho_mean * draw_units] * 2, '--', label=r'$\rho_{mean} = ' + my.f2s(rho_mean * draw_units) + ' \pm ' + my.f2s(d_rho_mean * draw_units) + '$')
    ax.legend()
    
    if(to_save_pics):
        pic_rho_t_filepath = os.path.join(model_path, 'rho(t).png')
        fig.savefig(pic_rho_t_filepath)

if(do_L):
    def plot_L(ax, L, t, lbl, clr, t_cut=75):
        L_mean = np.mean(L[t > t_cut])
        d_L_mean = np.std(L[t > t_cut])
        ax.plot(t, L, '.', label=lbl + r' = $' + my.f2s(L_mean) + r' \pm ' + my.f2s(d_L_mean) + r'$', color=clr, alpha=0.03)
        ax.plot([max([min(t), t_cut]), max(t)], [L_mean] * 2, color=clr, label='_nolegend_')
        
        return L_mean, d_L_mean

    fig_z, ax_z = my.get_fig(r'$t$ (ns)', r'$L_z$ (nm)', title=r'$L_z(t)$; ' + title_suff)
    Lz_mean, d_Lz_mean = plot_L(ax_z, Lbox[2, :], time_xvg, r'$L_z$', clr_z)
    ax_z.legend()
    
    fig_xy, ax_xy = my.get_fig(r'$t$ (ns)', r'$L$ (nm)', title=r'$L(t)$; ' + title_suff)
    Lx_mean, d_Lx_mean = plot_L(ax_xy, Lbox[0, :], time_xvg, r'$L_x$', clr_x)
    Ly_mean, d_Ly_mean = plot_L(ax_xy, Lbox[1, :], time_xvg, r'$L_y$', clr_y)
    ax_xy.legend()
    
    print('Lx, Ly: ', [Lx_mean, d_Lx_mean, Ly_mean, d_Ly_mean])
    
    if(to_save_pics):
        pic_Lz_t_filepath = os.path.join(model_path, 'Lz(t).png')
        pic_Lxy_t_filepath = os.path.join(model_path, 'Lxy(t).png')
        fig_z.savefig(pic_Lz_t_filepath)
        fig_xy.savefig(pic_Lxy_t_filepath)

if(do_P):
    def plot_pressure(ax, t, P, dt, title, clr, av1=51, av2=251, lw=2, t_cut=30):
        P_av1 = rollavg_convolve_edges(P, av1)
        P_av2 = rollavg_convolve_edges(P, av2)
        P_mean = np.mean(P[t > t_cut])
        d_P_mean = np.std(P[t > t_cut])
        ax.plot(t, P, color = clr, alpha=0.1, label = title + ' = $' + my.f2s(P_mean) + r' \pm ' + my.f2s(d_P_mean) + '$', linewidth=lw)
        ax.plot(t, P_av1, color = clr, alpha=0.5, label='avg = ' + my.f2s(dt * av1) + ' ns', linewidth=lw)
        ax.plot(t, P_av2, color = clr, label='avg = ' + my.f2s(dt * av2) + ' ns', linewidth=lw)
        return P_mean, d_P_mean
    
    fig_Pz, ax_Pz = my.get_fig(r'$t$ (ns)', r'$P_z$ (nm)', title=r'$P_z(t)$; ' + title_suff)
    Pz_mean, d_Pz_mean = plot_pressure(ax_Pz, time_xvg, Pbox[2, :], dt_xvg, '$P_z$; ', clr_z, av2=step_hist_av)
    ax_Pz.legend()

    fig_Pxy, ax_Pxy = my.get_fig(r'$t$ (ns)', r'$P$ (bar)', title=r'$P(t)$; ' + title_suff)
    Py_mean, d_Py_mean = plot_pressure(ax_Pxy, time_xvg, Pbox[1, :], dt_xvg, '$P_y$', clr_y, av2=step_hist_av)
    Px_mean, d_Px_mean = plot_pressure(ax_Pxy, time_xvg, Pbox[0, :], dt_xvg, '$P_x$; ', clr_x, av2=step_hist_av)
    Pxz_mean, d_Pxz_mean = plot_pressure(ax_Pxy, time_xvg, Ptbox[1, :], dt_xvg, '$P_{xz}$; ', tuple([c * 0.8 for c in clr_x]), av2=step_hist_av, lw=1)
    Pyz_mean, d_Pyz_mean = plot_pressure(ax_Pxy, time_xvg, Ptbox[2, :], dt_xvg, '$P_{yz}$; ', tuple([c * 0.8 for c in clr_y]), av2=step_hist_av, lw=1)
    ax_Pxy.legend()
        
    fig_Pav, ax_Pav = my.get_fig(r'$t$ (ns)', r'$P$ (bar)', title=r'$P(t)$; avg = ' + my.f2s(dt_xvg * step_hist_av) + ' (ns); ' + title_suff)
    ax_Pav.plot(time_xvg, rollavg_convolve_edges(Pbox[0, :], step_hist_av), \
                color=clr_x, label=r'$P_x$')
    ax_Pav.plot(time_xvg, rollavg_convolve_edges(Pbox[1, :], step_hist_av), \
                color=clr_y, label=r'$P_y$')
    ax_Pav.plot(time_xvg, rollavg_convolve_edges(Pbox[2, :], step_hist_av), \
                color=clr_z, label=r'$P_z$')
    ax_Pav.plot(time_xvg, rollavg_convolve_edges(Ptbox[1, :], step_hist_av), \
                color=tuple([c * 0.8 for c in clr_x]), label=r'$P_{xz}$', linewidth=1)
    ax_Pav.plot(time_xvg, rollavg_convolve_edges(Ptbox[2, :], step_hist_av), \
                color=tuple([c * 0.8 for c in clr_y]), label=r'$P_{yz}$', linewidth=1)
    ax_Pav.legend()

    if(to_save_pics):
        pic_Pz_t_filepath = os.path.join(model_path, 'Pz(t).png')
        pic_Pxy_t_filepath = os.path.join(model_path, 'Pxy(t).png')
        pic_Pav_t_filepath = os.path.join(model_path, 'Pav(t).png')
        fig_Pz.savefig(pic_Pz_t_filepath)
        fig_Pxy.savefig(pic_Pxy_t_filepath)
        fig_Pav.savefig(pic_Pav_t_filepath)

if(do_dist):
    fig_dist, ax_dist = my.get_fig(r'z (nm)', r'$\rho_w$ (kg / $m^3$)', title=r'$\rho_w(z)$; ' + title_suff, yscl='log')
    plt.ion()
    plt.show()
    y_lims = [np.amin(water_hist_average[water_hist_average > 0]) * 0.9, np.amax(water_hist_average) * 1.1]
    #print(Zmin)
    for ti in range(N_av):
        ax_dist.clear()
        plt.yscale('log')
        #print(np.transpose(np.array([Zcenters, water_hist_average[ti, :]])))
        ax_dist.bar(Zcenters, water_hist_average[ti, :], width=Zstep, \
               facecolor=my.colors[4], edgecolor=my.colors[4])
        ax_dist.plot([min(Zcenters), max(Zcenters)], [rho_atm] * 2, color=my.colors[1], label=r'$\rho_{sat}$')
        ax_dist.set_ylim(y_lims)            

        ax_dist.legend()
        fig_dist.suptitle(r'$\rho_w(z)$' + \
                          r'; $t \in [' + my.f2s(time[hist_av_inds[ti, 0]]) + r' - ' + my.f2s(time[hist_av_inds[ti, 1] - 1]) + r']$ (ns)' + \
                          r'; $H_2 O_{ext}$ = ' + str(extra_water) + \
                          r'; T = ' + my.f2s(Tmp))   # + 1014 - 180 for old w_extra
        plt.draw()
        if(to_save_pics):
            pic_rho_dist_filepath = os.path.join(model_path, 'rho_dist_t' + my.f2s(time[hist_av_inds[ti, 0]]) + '_' + my.f2s(time[hist_av_inds[ti, 1] - 1]) + '.png')
            fig_dist.savefig(pic_rho_dist_filepath)
        else:
            plt.pause(0.001)
            input("Press Enter to continue...")
else:
    plt.show()

print('ALL DONE')

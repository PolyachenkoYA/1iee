import sys
import os
import numpy as np
import subprocess as sp
import shutil
import matplotlib.pyplot as plt

import gromacs.formats as gmx

# ========== paths ============
root_path = sp.run(['git', 'rev-parse', '--show-toplevel'], stdout=sp.PIPE, text=True).stdout[:-1]  # -1 to cut the '\n'
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'exe')
xvgres_path = os.path.join(root_path, 'res', 'xvg')
gmx_exe = 'gmx_mpi'
save_mode = 'save'
draw_mode = 'draw'
process_mode = 'proc'
short_mode = 'short'
temperature_feat_str = 'T'
pressure_feat_str = 'P'
modes_flag = '-m'
features_flag = '-f'
all_modes = [process_mode, draw_mode, save_mode, short_mode]
all_features = [temperature_feat_str, pressure_feat_str]
energy_features_ids = {temperature_feat_str : '15',
                       pressure_feat_str    : '17'}

# ========= path-independent defines =========

def f2str(x, n=3):
    return '%s' % float(('%.' + str(n) + 'g') % x)

def get_fig(xlbl, ylbl, title=None, xscl='linear', yscl='linear', projection=None, zscl='linear', zlbl='z'):
    if(title is None):
        title = ylbl + '(' + xlbl + ')'
    fig = plt.figure()
    ax = fig.gca(projection=projection)
    plt.title(title)
    plt.xlabel(xlbl)
    plt.ylabel(ylbl)
    plt.xscale(xscl)
    plt.yscale(yscl)
    if(projection == '3d'):
        ax.set_zscale(zscl)
        ax.set_zlabel(zlbl)
    return fig, ax


def safe_copy(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(src, dst)

def index_safe(list, element, default=-1):
    return list.index(element) if element in list else default


# ========= path-dependent defines =========

def print_usage_and_exit(exe_name, code=1):
    print('usage:\n' + exe_name + '   [' + modes_flag + '    (' + '/'.join(all_modes) + ')]   [' + features_flag + '   (' + '/'.join(all_features) + ')]')
    exit(code)

def process_model(model_name, trajectory, modes=[process_mode], features=[pressure_feat_str]):
    model_path = os.path.join(run_path, model_name)
    os.chdir(model_path)
    
    base_filename = model_name + '_' + trajectory
    xvg_filename = base_filename + '_' + '_'.join(features) + '.xvg'
    xvg_filepath = os.path.join(xvgres_path, xvg_filename)

    if(not os.path.isfile(xvg_filepath)):
        if(not process_mode in modes):
            print('file "' + xvg_filepath + '" was not found.\nGenerating from the trajectory "' + model_path + '_' + trajectory + '"')
        modes.append(process_mode)

    if(process_mode in modes):
        input_line = ' '.join([energy_features_ids[f] for f in features] + ['0'])
        sp.run([gmx_exe, 'energy', '-f', trajectory + '.edr', '-o', xvg_filename], input=input_line, text=True)
        safe_copy(xvg_filename, xvg_filepath)
        print('"' + xvg_filepath + '" saved')

    xvg_file = gmx.XVG()
    xvg_file.read(xvg_filepath)
    N_fields = len(xvg_file.names)
    time = xvg_file.array[0]
    cut_time = cut_times[model_name]
    stab_time_ind = (time < cut_time)
    data = xvg_file.array[1:]
    mean_vals = np.zeros(N_fields)
    std_vals = np.zeros(N_fields)
    for i, field in enumerate(xvg_file.names):
        cut_data = data[i][stab_time_ind]
        mean_val = np.mean(cut_data)
        val_std = np.std(cut_data)
        mean_vals[i] = mean_val
        std_vals[i] = val_std
        if(((draw_mode in modes) or (save_mode in modes)) and not (short_mode in modes)):
            plot_title = '$t_{cut} = ' + f2str(cut_time) + '$; ' + field + ' = $' +  f2str(mean_val) + '\pm' + f2str(val_std) + '$'
            fig, ax = get_fig('time (ps)', field, title=plot_title)
            xvg_file.plot(maxpoints=None, columns=[0, i+1])
            
            if(save_mode in modes):
                pic_path = os.path.join(xvgres_path, base_filename + '_' + field + '.jpg')
                plt.savefig(pic_path)
                print('"' + pic_path + '" saved')

    return (mean_vals, std_vals)


# ========== arg parse =========
def pick_args(list, separators, ind):
    return list[ind+1 : min(separators[separators > ind])]

args = sys.argv[1:]
argc = len(args)

mode_args_ind = index_safe(args, modes_flag)
features_args_ind = index_safe(args, features_flag)
flags_inds = np.array([mode_args_ind, features_args_ind, len(args)])
modes = pick_args(args, flags_inds, mode_args_ind)
features = pick_args(args, flags_inds, features_args_ind)
features.sort(key = lambda f: int(energy_features_ids[f]))
for m in modes:
    if(not (m in all_modes)):
        print_usage_and_exit(sys.argv[0])
for f in features:
    if(not (f in all_features)):
        print_usage_and_exit(sys.argv[0])

# ========== process ===========
target_pressure = 1
maxsol = [1024, 1088, 1152, 1216, 1274]
#maxsol = [1024, 1274]

model_names = [str(n) for n in maxsol]
cut_times = {}
for m in model_names:
    cut_times[m] =  1000
N_models = len(model_names)
print(model_names)

means = [[]] * N_models
stds = [[]] * N_models
for i, model_name in enumerate(model_names):
    means[i], stds[i] = process_model(model_name, 'nvt', modes, features)
means = np.array(means).T
stds = np.array(stds).T
print('values:\n', means)
print('stds:\n', stds)

pressure_feat_ind = index_safe(features, pressure_feat_str)
pressure = means[pressure_feat_ind]
d_pressure = stds[pressure_feat_ind]
pfit2 = np.polyfit(maxsol, pressure, 2, w=1/d_pressure)
x_fit = np.linspace(min(maxsol), max(maxsol), 1000)
pressure_fit = np.polyval(pfit2, x_fit)
pfit2[2] -= target_pressure
maxsol_roots = np.roots(pfit2)
target_maxsol = maxsol_roots[np.argmin(abs(maxsol_roots - maxsol[np.argmin(pressure)]))]
slope = 2 * pfit2[0] * target_maxsol + pfit2[1]
print('target maxsol = ', target_maxsol)
print('slope = ', slope)
if((draw_mode in modes) or (save_mode in modes)):
    fig, ax = get_fig('maxsol', 'Pressure (bar)', title='Pressure(maxsol)')
    ax.errorbar(maxsol, pressure, yerr=d_pressure,
                marker='o', markerfacecolor='None', linestyle='', color='black')
    ax.plot(x_fit, pressure_fit, color='black')
    if(save_mode in modes):
        pic_path = os.path.join(xvgres_path, 'Pressure.jpg')
        plt.savefig(pic_path)
        print('"' + pic_path + '" saved')

    if(draw_mode in modes):
        plt.show()


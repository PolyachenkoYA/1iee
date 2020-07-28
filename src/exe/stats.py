import sys
import os
import numpy as np
import subprocess as sp
import shutil
import matplotlib.pyplot as plt
import mylib as my

import gromacs.formats as gmx

# ========== paths ============
root_path = sp.run(['git', 'rev-parse', '--show-toplevel'], stdout=sp.PIPE, text=True).stdout[:-1]  # -1 to cut the '\n'
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'exe')
default_output_dir = 'xvg'
gmx_exe = 'gmx_mpi'
save_mode = 'save'
draw_mode = 'draw'
process_mode = 'proc'
short_mode = 'short'
temperature_feat_str = 'T'
pressure_feat_str = 'P'
modes_flag = '-mode'
features_flag = '-feat'
def_output_dir_flag = '-dir'
all_modes = [process_mode, draw_mode, save_mode, short_mode]
all_features = [temperature_feat_str, pressure_feat_str]
energy_features_ids = {temperature_feat_str : '15',
                       pressure_feat_str    : '17'}

# ========= path-independent defines =========
def y_range(y, margin=0.05):
    mn = min(y)
    mx = max(y)
    r = mx - mn
    return (mn - r * margin, mx + r * margin)

def print_usage_and_exit(str='', exe_name=sys.argv[0], code=1):
    print('usage:\n' + exe_name + '   [ ' + def_output_dir_flag + ' output_dir (' + default_output_dir + ')]   [' + modes_flag + '    (' + '/'.join(all_modes) + ')]   [' + features_flag + '   (' + '/'.join(all_features) + ')]')
    if(str != ''):
        print(str)
    exit(code)

# ========= path-dependent defines =========

def process_model(model_path, trajectory, cut_time, xvgres_path=default_output_dir, modes=[process_mode], features=[pressure_feat_str]):
    model_name = model_path.replace('/', '_')
    model_path = os.path.join(run_path, model_path)
    base_filename = '_'.join([model_name, trajectory])
    xvg_filename = '_'.join([base_filename] + features) + '.xvg'
    xvg_filepath = os.path.join(xvgres_path, xvg_filename)

    os.chdir(model_path)
    if(not os.path.isfile(xvg_filepath)):
        if(not process_mode in modes):
            print('file "' + xvg_filepath + '" was not found.\nGenerating from the trajectory "' + model_path + '_' + trajectory + '"')
        modes.append(process_mode)

    if(process_mode in modes):
        input_line = ' '.join([energy_features_ids[f] for f in features] + ['0'])
        sp.run([gmx_exe, 'energy', '-f', trajectory + '.edr', '-o', xvg_filename], input=input_line, text=True)
        my.safe_copy(xvg_filename, xvg_filepath)
        print('"' + xvg_filepath + '" saved')

    xvg_file = gmx.XVG()
    xvg_file.read(xvg_filename)
    N_fields = len(xvg_file.names)
    time = xvg_file.array[0]
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
            plot_title = '$t_{cut} = ' + my.f2str(cut_time) + '$; ' + field + ' = $' +  my.f2str(mean_val) + '\pm' + my.f2str(val_std) + '$'
            fig, ax = my.get_fig('time (ps)', field, title=plot_title)
            xvg_file.plot(maxpoints=None, columns=[0, i+1])
            
            if(save_mode in modes):
                pic_path = os.path.join(xvgres_path, base_filename + '_' + field + '.jpg')
                plt.savefig(pic_path)
                print('"' + pic_path + '" saved')

    return (mean_vals, std_vals)


# ========== arg parse =========

args = sys.argv[1:]
argc = len(args)
if(argc < 2):
    print_usage_and_exit()

[output_dir, modes, features] = my.parse_args(args, [def_output_dir_flag, modes_flag, features_flag])
if(isinstance(features, str)):
    features = [features]
features.sort(key = lambda f: int(energy_features_ids[f]))

if(len(output_dir) == 0):
    output_dir = default_output_dir
elif(not isinstance(output_dir, str)):
    print_usage_and_exit(str='output dir "' + output_dir + '" is invalid')

for m in modes:
    if(not (m in all_modes)):
        print_usage_and_exit(str='mode "' + m + '" is invalid')
for f in features:
    if(not (f in all_features)):
        print_usage_and_exit(str='feature "' + f + '" is invalid')
res_path = os.path.join(root_path, 'res', output_dir)
model_path = os.path.join(run_path, output_dir)

# ========== process ===========
target_pressure = 1
stab_time = 1000  # ps
maxsol = [1024, 1088, 1152, 1216, 1274]
#maxsol = [1024, 1274]
maxsol = [1053, 1054, 1055, 1056]
maxsol = [1053, 1054, 1055, 1056, 1024, 1088, 1152, 1216, 1274]
maxsol = [1053, 1054, 1055, 1056, 1024, 1088]

model_names = [str(n) for n in maxsol]
N_models = len(model_names)
print('models:\n', model_names)
print('cut time = ', stab_time)

means = [[]] * N_models
stds = [[]] * N_models
for i, model_name in enumerate(model_names):
    means[i], stds[i] = process_model(os.path.join(output_dir, model_name), 'nvt', stab_time, modes=modes, features=features, xvgres_path=res_path)
means = np.array(means).T
stds = np.array(stds).T
print('values:\n', means)
print('stds:\n', stds)

pressure_feat_ind = my.index_safe(features, pressure_feat_str)
pressure = means[pressure_feat_ind]
d_pressure = stds[pressure_feat_ind]
pfit2 = np.polyfit(maxsol, pressure, 2, w=1/d_pressure)
pfit1 = np.polyfit(maxsol, pressure, 1, w=1/d_pressure)
x_fit = np.linspace(min(maxsol), max(maxsol), 1000)
pressure_fit1 = np.polyval(pfit1, x_fit)
pfit1[1] -= target_pressure
target_maxsol1 = np.roots(pfit1)[0]
slope1 = pfit1[0]

pressure_fit2 = np.polyval(pfit2, x_fit)
pfit2[2] -= target_pressure
maxsol_roots2 = np.roots(pfit2)
target_maxsol2 = maxsol_roots2[np.argmin(abs(maxsol_roots2 - maxsol[np.argmin(pressure)]))]
slope2 = 2 * pfit2[0] * target_maxsol2 + pfit2[1]

print('target maxsol1 = ', target_maxsol1)
print('slope1 = ', slope1)
print('target maxsol2 = ', target_maxsol2)
print('slope2 = ', slope2)
if((draw_mode in modes) or (save_mode in modes)):
    fig, ax = my.get_fig('maxsol', 'Pressure (bar)', title='Pressure(maxsol)')
    ax.errorbar(maxsol, pressure, yerr=d_pressure,
                marker='o', markerfacecolor='None', linestyle='', color='black')
    ax.plot(x_fit, pressure_fit1, color='black', label='lin fit; ms = ' + my.f2str(target_maxsol1, 6))
    ax.plot(x_fit, pressure_fit2, '--', color='black', label='q fit; ms = ' + my.f2str(target_maxsol2, 6))
    ax.set_ylim(y_range(pressure))
    ax.legend()
    if(save_mode in modes):
        pic_path = os.path.join(xvgres_path, 'Pressure_' + output_dir.replace('/', '_') + '.jpg')
        plt.savefig(pic_path)
        print('"' + pic_path + '" saved')

    if(draw_mode in modes):
        plt.show()


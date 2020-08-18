import sys
import os
import re
import numpy as np
import subprocess as sp
import shutil
import matplotlib.pyplot as plt
import mylib as my

import gromacs.formats as gmx

# ========== paths ============
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'maxsol')
res_path = os.path.join(root_path, 'res')
default_output_dir = 'xvg'
gmx_exe = 'gmx_mpi'
#gmx_exe = 'gmx_angara'
gmx_exe = 'gmx_serial'
save_mode = 'save'
draw_mode = 'draw'
process_mode = 'proc'
short_mode = 'short'
temperature_feat_str = 'T'
pressure_feat_str = 'P'
modes_flag = '-mode'
features_flag = '-feat'
output_dir_flag = '-dir'
stab_time_flag = '-stab_time'
maxsol_prefix = 'maxsol'
all_modes = [process_mode, draw_mode, save_mode, short_mode]
all_features = [temperature_feat_str, pressure_feat_str]
energy_features_ids = {temperature_feat_str : '15',
                       pressure_feat_str    : '17'}
stab_time_default = 1000

# ========= path-independent defines =========
def y_range(y, margin=0.05):
    mn = min(y)
    mx = max(y)
    r = mx - mn
    return (mn - r * margin, mx + r * margin)

def print_usage_and_exit(custom_str='', exe_name=sys.argv[0], code=1):
    print('usage:\n' + exe_name + '   ' + \
                       output_dir_flag + ' output_dir   ' + \
                 '[' + modes_flag + '    (' + '/'.join(all_modes) + ')]   ' + \
                 '[' + features_flag + '   (' + '/'.join(all_features) + ')]   ' + \
                 '[' + stab_time_flag + '   stab_time (' + str(stab_time_default) + ')]')
    if(custom_str != ''):
        print(custom_str)
    exit(code)

# ========= path-dependent defines =========

def process_model(model_path, trajectory, cut_time, xvgres_path=default_output_dir, modes=[process_mode], features=[pressure_feat_str]):
    model_name = model_path.replace('/', '_')
    model_path = os.path.join(run_path, model_path)
    base_filename = '_'.join([model_name, trajectory])
    xvg_filename = '_'.join([base_filename] + features) + '.xvg'
    xvg_filepath = os.path.join(xvgres_path, xvg_filename)

    if(not os.path.isfile(xvg_filepath)):
        if(not process_mode in modes):
            print('file "' + xvg_filepath + '" was not found.\nGenerating from the trajectory "' + model_path + '_' + trajectory + '"')
        modes.append(process_mode)

    if(process_mode in modes):
        os.chdir(model_path)
        input_line = ' '.join([energy_features_ids[f] for f in features] + ['0'])
        cmd = [gmx_exe, 'energy', '-f', trajectory + '.edr', '-o', xvg_filename]
        print(cmd)
        sp.run(cmd, input=input_line, text=True)
        my.safe_copy(xvg_filename, xvg_filepath)
        print('"' + xvg_filepath + '" saved')

    xvg_file = gmx.XVG()
    xvg_file.read(xvg_filepath)
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
            ax.plot(np.array([1, 1]) * cut_time, [min(data[i]), max(data[i])], color='red')
            
            if(save_mode in modes):
                pic_path = os.path.join(xvgres_path, base_filename + '_' + field + '.jpg')
                plt.savefig(pic_path)
                print('"' + pic_path + '" saved')

    return (mean_vals, std_vals)


# ========== arg parse =========

args = sys.argv[1:]
argc = len(args)
if(argc < 6):
    print_usage_and_exit()

[output_dir, modes, features, stab_time], correct_input = \
    my.parse_args(args, \
                  [output_dir_flag, modes_flag, features_flag, stab_time_flag], \
                  possible_values=[None, all_modes, all_features, None], \
                  possible_arg_numbers=[[1], ['+'], [1, 2], [0, 1]])
features.sort(key = lambda f: int(energy_features_ids[f]))
stab_time = int(stab_time[0]) if len(stab_time) > 0 else stab_time_default

model_res_path = os.path.join(res_path, output_dir)
model_path = os.path.join(run_path, output_dir)

# ========== process ===========
target_pressure = 1
maxsol = [1024, 1088, 1152, 1216, 1274]
#maxsol = [1024, 1274]
maxsol = [1053, 1054, 1055, 1056]
maxsol = [1053, 1054, 1055, 1056, 1024, 1088, 1152, 1216, 1274]
maxsol = [1053, 1054, 1055, 1056, 1024, 1088]
maxsol = [1040, 1050, 1055, 1060, 1070]
maxsol = []

if(len(maxsol) == 0):
    dir_for_list = os.path.join(res_path, output_dir)
    if(not os.path.isdir(dir_for_list)):
        dir_for_list = os.path.join(run_path, output_dir)
    #maxsol = [int(maxsol_dir[len(maxsol_prefix):]) for maxsol_dir in os.listdir(dir_for_list)]
    maxsol = [int(re.search(r'.*' + maxsol_prefix + r'(\d+).*', maxsol_dir).group(1)) for maxsol_dir in os.listdir(dir_for_list)]
    
model_names =  [(maxsol_prefix + str(n)) for n in maxsol]
N_models = len(model_names)
print('models:\n', model_names)
print('cut time = ', stab_time)

means = [[]] * N_models
stds = [[]] * N_models
for i, model_name in enumerate(model_names):
    means[i], stds[i] = process_model(os.path.join(output_dir, model_name), 'nvt', stab_time, modes=modes, features=features, xvgres_path=model_res_path)
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
        pic_path = os.path.join(os.path.join(res_path, os.path.dirname(output_dir)), 'Pressure_' + output_dir.replace('/', '_') + '.jpg')
        plt.savefig(pic_path)
        print('"' + pic_path + '" saved')
        maxsol_filename = os.path.join(res_path, os.path.dirname(output_dir) + '_' + str(stab_time) + '.txt')
        with open(maxsol_filename, 'a') as outfile:
            print(target_maxsol1, file=outfile, end=' ')
            print(maxsol_filename, 'appended')

    if(draw_mode in modes):
        plt.show()


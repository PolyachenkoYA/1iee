import sys
import os
import numpy as np
import subprocess as sp
import shutil
import matplotlib.pyplot as plt

import gromacs.formats as gmx

# ==============================================

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

def print_usage_and_exit(exe_name, code=1):
    print('usage:\n' + exe_name + '   [modes (' + process_mode + '/' + draw_mode + '/' + save_mode + ')]')
    exit(code)

def safe_copy(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(src, dst)

def process_model(model_name, trajectory, modes=[None]):
    model_path = os.path.join(run_path, model_name)
    os.chdir(model_path)
    
    base_filename = model_name + '_' + trajectory
    xvg_filename = base_filename + '_T_P.xvg'
    xvg_filepath = os.path.join(xvgres_path, xvg_filename)

    if(process_mode in modes):
        sp.run([gmx_exe, 'energy', '-f', trajectory + '.edr', '-o', xvg_filename], input='15 17 0', text=True)
        safe_copy(xvg_filename, xvg_filepath)

    if(draw_mode in modes):
        xvg_file = gmx.XVG()
        xvg_file.read(xvg_filepath)
        time = xvg_file.array[0]
        cut_time = cut_times[model_name]
        stab_time_ind = (time < cut_time)
        data = xvg_file.array[1:]
        for i, field in enumerate(xvg_file.names):
            cut_data = data[i][stab_time_ind]
            mean_val = np.mean(cut_data)
            val_std = np.std(cut_data)
            plot_title = '$t_{cut} = ' + f2str(cut_time) + '$; ' + field + ' = $' +  f2str(mean_val) + '\pm' + f2str(val_std) + '$'
            fig, ax = get_fig('time (ps)', field, title=plot_title)
            xvg_file.plot(maxpoints=None, columns=[0, i+1])
            
            if(save_mode in modes):
                pic_path = os.path.join(xvgres_path, base_filename + '_' + field + '.jpg')
                plt.savefig(pic_path)
                print('"' + pic_path + '" saved')

# ==================================================

# ========== paths ============
root_path = sp.run(['git', 'rev-parse', '--show-toplevel'], stdout=sp.PIPE, text=True).stdout[:-1]  # -1 to cut the '\n'
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'exe')
xvgres_path = os.path.join(root_path, 'res', 'xvg')
gmx_exe = 'gmx_mpi'
save_mode = 'save'
draw_mode = 'draw'
process_mode = 'proc'

# ========== arg parse =========
args = sys.argv[1:]
argc = len(args)

modes = [process_mode] if argc == 0 else args
for mode in modes:
    if(not (mode in [process_mode, draw_mode, save_mode])):
        print_usage_and_exit(sys.argv[0])

# ========== process ===========
model_names = [str(n) for n in [1, 1024, 1088, 1152, 1216, 1274]]
model_names = ['1024']
cut_times = {}
for m in model_names:
    cut_times[m] =  1000

for model_name in model_names:
    process_model(model_name, 'nvt', modes)

if(draw_mode in modes):
    plt.show()

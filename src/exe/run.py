import os
import sys
import numpy as np
import subprocess as sp

root_path = sp.run(['git', 'rev-parse', '--show-toplevel'], stdout=sp.PIPE, text=True).stdout[:-1]  # -1 to cut the '\n'
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'exe')
nvtmdp_filename = 'nvt.mdp'

maxsol_arr = [1040, 1050, 1055, 1060, 1070]
temperatue_arr = [1, 10, 15, 20, 25, 30, 40] + 273
jobs = ['112', '113', '114', '122']

maxsol_arr = [1050, 1060]
temperature_arr = [10, 20] + 273
temperature = [20] + 273
jobs = ['112', '113']

for j in jobs:
    for t in temperature_arr:
        for ms in maxsol_arr:
            path = os.path.join(j, str(t), str(ms))
            mdp_filepath = os.path.join(run_path, path, nvtmdp_filename)
            sp.run('./clear_restore.sh ' + path, shell=True)
            sp.run('python change_mdp.py -in ' + mdp_filepath + ' -out ' + mdp_filepath + ' -flds gen_temp ' + str(t) + ' ref_t' + str(t))
            sp.run('./preproc.sh ' + path + ' 4 ' + str(ms), shell=True)
            sp.run('./mainrun.sh ' + path + ' 4', shell=True)

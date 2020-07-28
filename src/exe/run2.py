import os
import sys
import numpy as np
import subprocess as sp

def run_it(cmd, shell=False):
    print(cmd)
    sp.run(cmd, shell=shell)

root_path = sp.run(['git', 'rev-parse', '--show-toplevel'], stdout=sp.PIPE, text=True).stdout[:-1]  # -1 to cut the '\n'
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'exe')
nvtmdp_filename = 'nvt.mdp'
omp_cores = 4

maxsol_arr = np.array([1040, 1050, 1055, 1060, 1070])
temperatue_arr = np.array([1, 10, 15, 20, 25, 30, 35, 40]) + 273
jobs = np.array(['112', '113', '114', '122'])

#maxsol_arr = np.array([1050])
temperature_arr = np.array([25, 30]) + 273
jobs = np.array(['111'])

for j in jobs:
    for t in temperature_arr:
        for ms in maxsol_arr:
            path = os.path.join('job' + j, 't' + str(t), 'maxsol' + str(ms))
            start_pdb_file = '1iee' + j + '_prot4gmx.pdb'
            mdp_filepath = os.path.join(run_path, path, nvtmdp_filename)
            run_it('./clear_restore.sh ' + path, shell=True)
            run_it(['python', 'change_mdp.py', '-in', mdp_filepath, '-out', mdp_filepath, '-flds', 'gen_temp', str(t), 'ref_t', str(t)])
            run_it('./preproc.sh ' + path + ' ' + str(ms) + ' ' + start_pdb_file, shell=True)
            run_it('./mainrun.sh ' + path + ' ' + str(omp_cores), shell=True)

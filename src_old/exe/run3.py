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
omp_cores = 1
gpu_id = 0

# appromixate maxsol : P(maxsol) = 1 bar
# 1127, 1096, 1080, 1065, 1050, 1036, 1022, 1010
jobs = np.array([[1, 1, 2], [1, 1, 3], [1, 1, 4], [1, 2, 2]])
temperature_arr = np.array([1, 10]) + 273
maxsol_arr = np.array([[1112, 1122, 1127, 1132, 1142],
                       [1081, 1091, 1096, 1101, 1111],
                       [1065, 1075, 1080, 1085, 1075],
                       [1050, 1060, 1065, 1070, 1080],
                       [1035, 1045, 1050, 1055, 1065],
                       [1021, 1031, 1036, 1041, 1051],
                       [1007, 1017, 1022, 1027, 1032],
                       [995 , 1005, 1010, 1015, 1025]])

jobs = np.array([[1, 1, 4]])
#temperature_arr = np.array([20]) + 273
#maxsol_arr = np.array([1050])

for j in jobs:
    for t_i, t in enumerate(temperature_arr):
        for ms in maxsol_arr[t_i]:
            j_strs = [str(x) for x in j]
            j_str = ''.join(j_strs)
            path = os.path.join('job' + j_str, 't' + str(t), 'maxsol' + str(ms))
            start_pdb_file = '1iee' + j_str + '_prot4gmx.pdb'
            mdp_filepath = os.path.join(run_path, path, nvtmdp_filename)
            
            run_it('./clear_restore.sh ' + path, shell=True)
            run_it(['python', 'change_mdp.py', '-in', mdp_filepath, '-out', mdp_filepath, '-flds', 'gen_temp', str(t), 'ref_t', str(t)])
            run_it(' '.join(['./preproc.sh', path] + j_strs + [str(ms * np.prod(j)), start_pdb_file]), shell=True)
            run_it(' '.join(['./mainrun.sh', path, str(omp_cores), str(gpu_id)]), shell=True)

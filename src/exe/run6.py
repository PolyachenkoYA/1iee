import os
import sys
import numpy as np
import subprocess as sp

import mylib as my

def run_it(cmd, shell=False):
    print(cmd)
    sp.run(cmd, shell=shell)

root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'exe')
nvtmdp_filename = 'nvt.mdp'
omp_cores = 6
gpu_id = 2
T_C2K = 273.15

# appromixate maxsol : P(maxsol) = 1 bar
# 1127, 1096, 1080, 1065, 1050, 1036, 1022, 1010
jobs = np.array([[1, 1, 2], [1, 1, 3], [1, 1, 4], [1, 2, 2]])
temperature_arr = np.array([1, 10, 15, 20, 25, 30, 35, 40])
maxsol_arr = {1:  [1112, 1122, 1127, 1132, 1142],
              10: [1081, 1091, 1096, 1101, 1111],
              15: [1065, 1075, 1080, 1085, 1075],
              20: [1050, 1060, 1065, 1070, 1080],
              25: [1035, 1045, 1050, 1055, 1065],
              30: [1021, 1031, 1036, 1041, 1051],
              35: [1007, 1017, 1022, 1027, 1032],
              40: [995 , 1005, 1010, 1015, 1025]}

jobs = np.array([[1, 1, 1]])
temperature_arr = np.array([30])
#maxsol_arr = {20:[1050]}

for j in jobs:
    for t in temperature_arr:
        for ms in maxsol_arr[t]:
            j_strs = [str(x) for x in j]
            j_str = ''.join(j_strs)
            path = os.path.join('job' + j_str, 't' + str(t), 'maxsol' + str(ms))
            start_pdb_file = '1iee' + j_str + '_prot4gmx.pdb'
            mdp_filepath = os.path.join(run_path, path, nvtmdp_filename)
            
            run_it('./clear_restore.sh ' + path, shell=True)
            run_it(['python', 'change_mdp.py', '-in', mdp_filepath, '-out', mdp_filepath, '-flds', 'gen_temp', str(t + T_C2K), 'ref_t', str(t + T_C2K)])
            run_it(' '.join(['./preproc.sh', path] + j_strs + [str(ms * np.prod(j)), start_pdb_file]), shell=True)
            run_it(' '.join(['./mainrun.sh', path, str(omp_cores), str(gpu_id)]), shell=True)

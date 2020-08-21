import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing

import mylib as my

def run_it(cmd, shell=False):
    print(cmd)
    sp.run(cmd, shell=shell)
    
omp_cores = 6
gpu_id = -1

[output_dir, modes, features, stab_time], correct_input = \
    my.parse_args(args, ['-omp', '-gpu_id'], \
                  possible_values=[range(multiprocessing.cpu_count()), None], \
                  possible_arg_numbers=[])

# =============== paths ================
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'maxsol')
res_path = os.path.join(root_path, 'res')
nvtmdp_filename = 'nvt.mdp'

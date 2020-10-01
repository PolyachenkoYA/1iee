import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing

import mylib as my

def run_it(cmd, shell=False):
    print(cmd)
    sp.run(cmd, shell=shell)
    
cpu_N = multiprocessing.cpu_count()

[omp_N, gpu_id], correct_input = \
    my.parse_args(args, ['-omp', '-gpu_id'], \
                  possible_values=[range(cpu_N), None], \
                  possible_arg_numbers=[[0, 1], [0, 1]], \
                  default_values=[cpu_N, -1])

# =============== paths ================
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'maxsol')
res_path = os.path.join(root_path, 'res')
nvtmdp_filename = 'nvt.mdp'

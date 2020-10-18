import os
import sys
import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt

import mylib as my

# ============================= params =================================
omp = 3  # for serial, omp for mpi is 1
mpi = 4  # for MPI, mpi for omp>1 is 1
N_gpu = 4
dV_mult = 1.001
#host_ids = [27, 28, 29, 30]
temp_i = 7
id = 0
Ptau_i = 0
gpu_id = Ptau_i % N_gpu

root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'K')
res_path = os.path.join(root_path, 'res')

# =========================== input parse =======================================
args = sys.argv[1:]
argc = len(args)
if(argc not in [0]):
    print('usage:\n' + sys.argv[0])
    exit(1)
    
#mainrun_mode = int(args[0])
#mainrun_mode = 1
cmd = 'python run_K_flucts.py -omp ' + str(omp) + ' -mpi ' + str(mpi) + ' -gpu_id ' + str(gpu_id) + ' -param_ids ' + str(temp_i) + ' ' + str(Ptau_i) + ' -id ' + str(id) + ' -preproc_mode force'  # flucts
my.run_it(cmd)
# =============== launch tasks =========
#for id_i, id in enumerate(range(0, 1)):
#    for temp_i in range(7, 8):
        #cmd = 'python run_K_dV.py -omp ' + str(omp) + ' -mpi ' + str(mpi) + ' -gpu_id ' + str(gpu_id) + ' -mainrun_mode ' + str(mainrun_mode) + ' -preproc_mode if_needed -temp_ids ' + str(temp_i) + ' -id ' + str(id) + ' -dV_mult ' + str(dV_mult)  # dV 
        #cmd = 'python run_K_flucts.py -omp ' + str(omp) + ' -mpi ' + str(mpi) + ' -gpu_id ' + str(gpu_id) + ' -param_ids ' + str(temp_i) + ' ' + str(compr_i) + ' -id ' + str(id) + ' -preproc_mode force'  # flucts

        #global_i = id_i + temp_i * 1
        #host_i = host_ids[global_i % len(host_ids)]
        #job_name = 'gromacs' + 'F_' + str(temp_i) + '_' + str(id)
        #my.run_it('sbatch --nodelist=host' + str(host_i) + ' --reservation=test -J ' + job_name  + ' -p max1n -N 1 --ntasks-per-node=' + str(mpi) + ' --gres=gpu:1 --wrap="' + cmd + '"')
        #my.run_it("screen -d -m -S " + job_name + " bash -c '" + cmd + "'")
    
        #my.run_it(cmd)
# sbatch --reservation=test -J gromacs1 -p max1n -N 1 --ntasks-per-node=1 --gres=gpu:1 --wrap="python run.py -omp 6 -mpi 4 -do_mainrun 0 -time_ids 0"

# ================ extract  maxsol_0 for the matlab plotting ===============
#if(not len(args) in [1, 2]):
#    print('usage:\n' + sys.argv[0] + '   job_dir   [stab_time(1000)]')
#    exit(1)
#job_name = args[0]
#stab_time = 1000 if(len(args) < 2) else float(args[1])
#
#t_arr = [1, 10, 15, 20, 25, 30, 35, 40]
#txt_datafile_path = os.path.join(res_path, job_name + '_' + str(350) + '.txt')
#if(os.path.isfile(txt_datafile_path)):
#    os.remove(txt_datafile_path)
#for t in t_arr:
#    run_it('python stats.py -dir ' + os.path.join(job_name, 't' + str(t)) + ' -mode short save -feat P -stab_time ' + str(stab_time), shell=True)

# ================================== cell size convergence ============================

#for job_name in job_names:
#    for stab_time in stab_times:
#        txt_datafile_path = os.path.join(res_path, job_name + '_' + str(stab_time) + '.txt')
#        if(os.path.isfile(txt_datafile_path)):
#            os.remove(txt_datafile_path)
#        for t in t_arr:
#            run_it('python stats.py -dir ' + os.path.join(job_name, 't' + str(t)) + ' -mode short save -feat P -stab_time ' + str(stab_time), shell=True)

import os
import sys
import numpy as np
import subprocess as sp
import matplotlib.pyplot as plt

import mylib as my

# ============================= params =================================
omp = 12  # for serial, omp for mpi is 1
mpi = 4  # for MPI, mpi for omp>1 is 1
N_gpu = 4
dV_mult = 1.001
#host_ids = [27, 28, 29, 30]
temp_i = 7
compr_i = 1
time_i = 0
Ttau_i = 4
id = 0

root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'K')
res_path = os.path.join(root_path, 'res')

# =========================== input parse =======================================
args = sys.argv[1:]
argc = len(args)
if(argc not in [0]):
    print('usage:\n' + sys.argv[0] + '')
    exit(1)
    

# ======= flucts ===========
#time_i = int(args[0])
#Ptau_i = int(args[1])
#gpu_id = Ptau_i % N_gpu
#cmd = 'python run_K_flucts.py -omp ' + str(omp) \
#                          + ' -mpi ' + str(mpi) \
#                          + ' -gpu_id ' +  str(gpu_id) \
#                          + ' -param_ids ' + str(temp_i) \
#                                     + ' ' + str(compr_i) \
#                                     + ' ' + str(time_i) \
#                                     + ' ' + str(Ptau_i * 2 + 1) \
#                          + ' -id ' + str(id) \
#  
                        + ' -preproc_mode force'  # flucts
'''
# ===== dV ===============
gpu_id = Ttau_i % N_gpu
#mainrun_mode = int(args[0])
mainrun_mode = 1
cmd = 'python run_K_dV.py -omp ' + str(omp) \
                      + ' -mpi ' + str(mpi) \
                      + ' -gpu_id ' +  str(gpu_id) \
                      + ' -param_ids ' + str(temp_i) \
                                 + ' ' + str(time_i) \
                                 + ' ' + str(Ttau_i) \
                      + ' -id ' + str(id) \
                      + ' -preproc_mode force' \
                      + ' -dV_mult ' + str(dV_mult) \
                      + ' -mainrun_mode ' + str(mainrun_mode)
my.run_it(cmd)
'''

# =============== launch tasks =========
for dV_i, dV_mult in enumerate([3, 3.5, 4]):
    for Ttau_i in range(4):
        gpu_id = (Ttau_i + dV_i * 4) % N_gpu
        cmd = 'python run_K_dV.py -omp ' + str(omp) \
                              + ' -mpi ' + str(mpi) \
                              + ' -gpu_id ' +  str(gpu_id) \
                              + ' -param_ids '+ str(Ttau_i) \
                              + ' -id ' + str(id) \
                              + ' -preproc_mode force' \
                              + ' -dV_mult ' + str(dV_mult) \
                              + ' -mainrun_mode ' + str(mainrun_mode)

        #cmd = 'python run_K_flucts.py -omp ' + str(omp) + ' -mpi ' + str(mpi) + ' -gpu_id ' + str(gpu_id) + ' -param_ids ' + str(temp_i) + ' ' + str(compr_i) + ' -id ' + str(id) + ' -preproc_mode force'  # flucts

        #global_i = id_i + temp_i * 1
        #host_i = host_ids[global_i % len(host_ids)]
        #job_name = 'gromacs' + 'F_' + str(temp_i) + '_' + str(id)
        #my.run_it('sbatch --nodelist=host' + str(host_i) + ' --reservation=test -J ' + job_name  + ' -p max1n -N 1 --ntasks-per-node=' + str(mpi) + ' --gres=gpu:1 --wrap="' + cmd + '"')
        #my.run_it("screen -d -m -S " + job_name + " bash -c '" + cmd + "'")
    
        my.run_it(cmd)
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

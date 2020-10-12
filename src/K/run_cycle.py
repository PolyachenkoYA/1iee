import os
import sys
import numpy as np
import subprocess as sp

import mylib as my

omp = 6  # for serial, omp for mpi is 1
mpi = 4  # for MPI, mpi for omp>1 is 1
dV_mult = 1.001
host_ids = [27, 28, 29, 30]

root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'water_cube')
res_path = os.path.join(root_path, 'res')

args = sys.argv[1:]
argc = len(args)
if(argc not in [1]):
    print('usage:\n' + sys.argv[0] + '   marinrun_mode (1/2)')
    exit(1)
    
mainrun_mode = int(args[0])

# =============== launch tasks =========
for time_i in range(1):
    for id_i, id in enumerate(range(0, 1)):
        # cmd = 'python run.py -omp ' + str(omp) + ' -mpi ' + str(mpi) + ' -do_mainrun yes -preproc_mode force -time_ids ' + str(time_i) + ' -id ' + str(id)  # fluct water_cube K
        cmd = 'python run_K_dV.py -omp ' + str(omp) + ' -mpi ' + str(mpi) + ' -mainrun_mode ' + str(mainrun_mode) + ' -preproc_mode if_needed -time_ids ' + str(time_i) + ' -id ' + str(id) + ' -dV_mult ' + str(dV_mult)  # dV water_cube K
        job_name = 'gromacs' + str(time_i) + '_' + str(id) + '_' + str(mainrun_mode)

        #host_i = host_ids[(id_i + time_i * 4) % len(host_ids)]
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

# ================ post proc all ===============
#if(not len(args) in [0]):
#    print('usage:\n' + sys.argv[0])
#    exit(1)

#job_names = ['job111', 'job112', 'job113', 'job114']
#stab_times = [200, 400, 600, 800, 1000]
#t_arr = [1, 10, 15, 20, 25, 30, 35, 40]

#job_names = ['job111', 'job112']
#stab_times = [600, 800]

#for job_name in job_names:
#    for stab_time in stab_times:
#        txt_datafile_path = os.path.join(res_path, job_name + '_' + str(stab_time) + '.txt')
#        if(os.path.isfile(txt_datafile_path)):
#            os.remove(txt_datafile_path)
#        for t in t_arr:
#            run_it('python stats.py -dir ' + os.path.join(job_name, 't' + str(t)) + ' -mode short save -feat P -stab_time ' + str(stab_time), shell=True)

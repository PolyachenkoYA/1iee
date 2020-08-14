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
res_path = os.path.join(root_path, 'res')
args = sys.argv[1:]

# =============== launch tasks =========
#for job_i in range(3):
#    for temp_i in range(2,3):
#        run_it('sbatch --reservation=test -J gromacs' + str(i) + ' -p max1n -N 1 --ntasks-per-node=1 --gres=gpu:1 --wrap="python run.py -omp 6 -continue 1 -temp_ids ' + str(temp_i) + ' -job_ids ' + str(job_i) + '"', shell=True)
        #run_it("screen -d -m -S " + run_str + " bash -c 'python " + run_str + "'", shell=True)
# sbatch --reservation=test -J gromacs1 -p max1n -N 1 --ntasks-per-node=1 --gres=gpu:1 --wrap="python run.py"

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
if(not len(args) in [0]):
    print('usage:\n' + sys.argv[0])
    exit(1)

job_names = ['job111', 'job112', 'job113', 'job114']
stab_times = [200, 400, 600, 800, 1000]
t_arr = [1, 10, 15, 20, 25, 30, 35, 40]

#job_names = ['job111', 'job112']
#stab_times = [600, 800]

for job_name in job_names:
    for stab_time in stab_times:
        txt_datafile_path = os.path.join(res_path, job_name + '_' + str(stab_time) + '.txt')
        if(os.path.isfile(txt_datafile_path)):
            os.remove(txt_datafile_path)
        for t in t_arr:
            run_it('python stats.py -dir ' + os.path.join(job_name, 't' + str(t)) + ' -mode short save -feat P -stab_time ' + str(stab_time), shell=True)

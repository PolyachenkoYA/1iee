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

#for i in range(1, 3):
#    run_str = 'run' + str(i) + '.py'
    #run_it('sbatch -J gmx_' + str(i) + ' -p max8n -N 8 --ntasks-per-node=1 --gres=gpu:1 --wrap="python ' + run_str + '"', shell=True)
#    run_it("screen -d -m -S " + run_str + " bash -c 'python " + run_str + "'", shell=True)

t_arr = [274, 283, 288, 293, 298, 303, 308, 313]
#t_arr = [274, 283]
job_name = args[0]
txt_datafile_path = os.path.join(res_path, job_name + '.txt')
if(os.path.isfile(txt_datafile_path)):
    os.remove(txt_datafile_path)
for t in t_arr:
    run_it('python stats.py -dir ' + os.path.join(job_name, 't' + str(t)) + ' -mode short save -feat P', shell=True)

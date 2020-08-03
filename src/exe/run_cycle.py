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

for i in range(1, 3):
    run_str = 'run' + str(i) + '.py'
    #run_it('sbatch -J gmx_' + str(i) + ' -p max8n -N 8 --ntasks-per-node=1 --gres=gpu:1 --wrap="python ' + run_str + '"', shell=True)
    run_it("screen -d -m -S " + run_str + " bash -c 'python " + run_str + "'", shell=True)

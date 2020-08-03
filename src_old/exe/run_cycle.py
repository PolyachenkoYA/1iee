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

for i in range(1, 9):
    run_str = 'run' + str(i) + '.py'
    run_it("screen -d -m -S " + run_str + " bash -c 'python " + run_str + "'", shell=True)

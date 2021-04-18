import os
import sys

import mylib as my

comprs = [3e-6, 3e-8, 3e-10, 3e-12, 3e-14, 0]
Tmps = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
for i in range(len(comprs)):
#for i in range(len(Tmps)):
    my.run_it('python water_density_full.py -temp 35 -extra_water 30 -comprZ ' +  str(comprs[i]) + ' -save_pics 1')
    #my.run_it('sbatch -J gmx' + str(i) + ' -p max1n -N 1 --ntasks-per-node=1 --wrap="python run_K_flucts.py -param_ids 7 -omp 12 -mpi 1 -extra_water 30 -comprZ ' +  str(comprs[i]) + '"')

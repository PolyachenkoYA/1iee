import os
import sys

import mylib as my

comprs = [0, 3e-14]
Tmps = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
ids = [7, 0, 1, 2, 3, 4, 5]

ids = [1, 2, 3, 4]
comprs = [0]

for id in ids:
    for i in range(len(comprs)):
#for i in range(len(Tmps)):
        #my.run_it('python water_density_full.py -temp 35 -extra_water 30 -comprZ ' +  str(comprs[i]) + ' -save_pics 1 -id ' + str(id))
        my.run_it('sbatch -J gmx1_' + str(i) + '_' + str(id) + ' -p max1n -N 1 --ntasks-per-node=1 --wrap="python run_K_flucts.py -param_ids 7 ' + str(id) + ' -omp 12 -mpi 1 -extra_water 30 -comprZ ' +  str(comprs[i]) + '"')

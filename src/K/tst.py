import os
import sys

import mylib as my

comprs = [0, 3e-14]
Tmps = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
ids = [7, 0, 1, 2, 3, 4, 5]
extra_waters = [-40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110]
#extra_waters = [-40]

ids = [1, 2, 3, 4]
comprs = [0]

#for id in ids:
#    for i in range(len(comprs)):
#for i in range(len(Tmps)):
for i in range(len(extra_waters)):
        #my.run_it('python water_density_full.py -temp 35 -extra_water 30 -comprZ ' +  str(comprs[i]) + ' -save_pics 1 -id ' + str(id))
        my.run_it('sbatch -J gmx' + str(i) + ' -p max1n -N 1 --ntasks-per-node=1 --wrap="python run_K_flucts.py -param_ids 7 -omp 12 -mpi 1 -extra_water ' + str(extra_waters[i]) + ' -comprZ 0"')

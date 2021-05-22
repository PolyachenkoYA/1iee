import os
import sys

import mylib as my

comprs = [0, 3e-14]
Tmps = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
ids = [7, 0, 1, 2, 3, 4, 5]
extra_waters = [-40, 110, 0, -30, -20, -10, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
extra_waters = [-500,  -420,  -330,  -250,  -160,   -80,   180,   260,   350,   440,   520,   610]
extra_waters = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450]
extra_waters = [150, 200, 250]

ids = [1, 2, 3, 4]
comprs = [0]

#for id in ids:
#    for i in range(len(comprs)):
#for i in range(len(Tmps)):
for i in range(len(extra_waters)):
        my.run_it('python water_density_full.py -temp 35 -extra_water ' + str(extra_waters[i]) + ' -comprZ 0 -save_pics 1 -avg_time 10 -Zstep 1 -recomp 1')
        #my.run_it('sbatch -J gmx' + str(i) + ' -p max1n -N 1 --ntasks-per-node=1 --wrap="python run_K_flucts.py -param_ids 7 -omp 12 -mpi 1 -extra_water ' + str(extra_waters[i]) + ' -comprZ 0"')
        

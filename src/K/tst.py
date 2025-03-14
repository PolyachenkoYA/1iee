import os
import sys

import mylib as my

comprs = [0, 3e-14]
Tmps = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
ids = [7, 0, 1, 2, 3, 4, 5]
ids = [1, 2, 3, 4]
comprs = [0]

extra_waters = [-40, 110, 0, -30, -20, -10, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
extra_waters = [-500,  -420,  -330,  -250,  -160,   -80,   180,   260,   350,   440,   520,   610]
extra_waters = [0, 50, 100, 150, 200, 250, 300, 350, 400, 450]
extra_waters = [0, 100, 200, 300, 400]
extra_waters = [0, 150, 300, 450]
extra_waters = [600, 750]
extra_waters = [0, 150, 300, 450, 600, 750]
extra_waters = [450, 600, 750]
#extra_waters = [0]

Tmps = [8]  # 26,2
#Tmps = [2, 4, 8, 12, 18, 22, 27]

mpi = 2
omp = 6
partition = 'max1n'

#for id in ids:
#    for i in range(len(comprs)):
for ti in range(len(Tmps)):
    for i in range(len(extra_waters)):
        #my.run_it('python water_density_full.py -temp ' + str(Tmps[ti]) + ' -extra_water ' + str(extra_waters[i]) + ' -do_dist 1 -save_pics 1 -avg_time 10 -Zstep 1 -recomp 1')
        
        my.run_it('sbatch -J gmx' + str(Tmps[ti]) + '_' + str(extra_waters[i]) + ' -p ' + partition + ' -N ' + str(mpi) + ' --ntasks-per-node=' + str(omp) + ' --wrap="python run_K_flucts.py -Tmp ' + str(Tmps[ti])  + ' -omp ' + str(omp) + ' -mpi ' + str(mpi) + ' -extra_water ' + str(extra_waters[i]) + '"')
        
        #my.run_it('sbatch -J gmx'   + str(extra_waters[i]) + ' -p ' + partition + ' -N ' + str(mpi) + ' --ntasks-per-node=1                --wrap="python run_K_flucts.py -param_ids 8 -omp ' + str(omp) + ' -mpi ' + str(mpi) + ' -extra_water ' + str(extra_waters[i]) + '"')
        

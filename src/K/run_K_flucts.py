import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing

import mylib as my

# =============== paths ================
#root_path = my.git_root_path()
with open('git_root_path', 'r') as f:
    root_path = f.readline()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'K')
res_path = os.path.join(root_path, 'res')
main_mdp_filename_base = 'npt'
eql_mdp_filename_base = 'eql'
no_preproc = 'no'
preproc_if_permitted = 'ask'
preproc_if_needed = 'if_needed'
preproc_force = 'force'
minimE_filename_base = 'em'
init_pdb_filename = '1iee112_prot4gmx.pdb'

# ================ params =================

N_gpus = 1
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns
compr = 0.0003
time = 65
omp_default = multiprocessing.cpu_count()
equil_maxsol_poly = [-2.9516, 1117.2]   # maxsol = np.polyval(equil_maxsol_poly, T), [T] = C (not K)
temps = np.array([0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55])
P_taus = np.array([200, 400, 800, 1600, 3200, 6400, 12800, 25000, 50000, 100000, 200000, 400000, 800000, 1600000, 3200000, 6400000])
P_taus = np.array([4, 8, 16, 32, 64, 128, 256, 512])
comprs = np.array([2e-4, 3e-4, 4e-4])
times = np.array([20.0, 40.0])
#mdrun_mode = 'serial'
mdrun_mode = 'slurm'

# ============== arg parse ====================
N_omp_max = multiprocessing.cpu_count()
possible_gpu_ids = [str(i) for i in range(N_gpus)] + ['-1']
possible_omps = [str(i) for i in range(1, N_omp_max + 1)]
possible_preproc = [no_preproc, preproc_if_permitted, preproc_force, preproc_if_needed]
[omp_cores, mpi_cores, gpu_id, do_mainrun, preproc_mode, param_ids, extra_water, compressibility_Z], _ = \
    my.parse_args(sys.argv[1:], ['-omp', '-mpi', '-gpu_id', '-do_mainrun', '-preproc_mode', '-param_ids', '-extra_water', '-comprZ'], \
                  possible_values=[possible_omps, None, possible_gpu_ids, ['no', 'yes'], possible_preproc, None, None, None], \
                  possible_arg_numbers=[[0, 1], [0, 1], [0, 1], [0, 1], [0, 1], None, [0, 1], [0, 1]], \
                  default_values=[[str(omp_default)], ['0'], ['-1'], ['yes'], [preproc_if_needed], ['0'], ['0'], ['0']])
param_ids = [int(i) for i in param_ids]
omp_cores = int(omp_cores[0])
mpi_cores = int(mpi_cores[0])
gpu_id = int(gpu_id[0])
do_mainrun = (do_mainrun[0] == 'yes')
preproc_mode = preproc_mode[0]
extra_water = int(extra_water[0])
compressibility_Z = float(compressibility_Z[0])
# ===================== cycle ===================
# flucts K

temp = temps[param_ids[0]]
#time = 40.0
#gpu_id = param_ids[0] % N_gpus
#compr = comprs[param_ids[1]]
#time = times[param_ids[2]]
#for Ptau_i, P_tau in enumerate(P_taus[param_ids[3:]]):
for _ in range(1):
    maxsol = int((round(np.polyval(equil_maxsol_poly, temp)) + extra_water) * 2)
    nsteps = int(round(time / dt))
    
    model_name = 'flucts_t4p_temp' + my.f2str(temp) + '_extW' + str(extra_water) + '_AnisoXYComprZ' + str(compressibility_Z) + '_id' + str(param_ids[1])
    mdp_filepath = os.path.join(run_path, model_name, main_mdp_filename_base + '.mdp')
    eql_filepath = os.path.join(run_path, model_name, eql_mdp_filename_base + '.mdp')
    checkpoint_filepath = os.path.join(run_path, model_name, main_mdp_filename_base + '.cpt')
    continue_comp = do_mainrun and os.path.isfile(checkpoint_filepath)
        
    if(preproc_mode in [preproc_if_permitted, preproc_force, preproc_if_needed]):
        if(os.path.isfile(checkpoint_filepath)):
            print('WARNING:\npreproc_mode is ' + preproc_mode + ', the checkpoint file "' + checkpoint_filepath + '" was found.')
            if(preproc_mode == preproc_if_permitted):
                skip = (input('Proceed? (y/N)') != 'y')
            elif(preproc_mode == preproc_force):
                skip = False
            elif(preproc_mode == preproc_if_needed):
                skip = True
        else:
            skip = False
    
    else:    # no_preproc
        skip = True 
                        
    if(skip):
        print('Skipping preproc')
    else:
        print('Running the model from the initial "' + init_pdb_filename + '"')
        #sys.exit(1)
        my.run_it('./clear_restore.sh ' + model_name)
        my.run_it(['python', 'change_mdp.py', '-in', mdp_filepath, '-out', mdp_filepath, '-flds', 'ref-t', str(temp + T_C2K), \
                                                                                                  'nsteps', str(nsteps), \
                                                                                                  'gen-temp', str(temp + T_C2K), \
                                                                                                  'gen-seed', str(param_ids[1]),
                                                                                                  'compressibility', '3e-4 3e-4 ' + str(compressibility_Z) + ' 3e-4 0 0'])
#                                                                                                  'compressibility', '3e-4 ' + str(compressibility_Z)])
        my.run_it(' '.join(['./preproc.sh', model_name, str(omp_cores), str(mpi_cores), str(gpu_id), '1', '1', '2', str(maxsol), init_pdb_filename]))
        my.run_it(' '.join(['./mainrun.sh', model_name, str(omp_cores), str(mpi_cores), str(gpu_id), '1iee_wions', minimE_filename_base, mdrun_mode, '1']))
        my.run_it(' '.join(['./save_initial_nojump.sh', model_name, minimE_filename_base]))

    if(do_mainrun):
        my.run_it(' '.join(['./mainrun.sh', model_name, str(omp_cores), str(mpi_cores), str(gpu_id), minimE_filename_base, main_mdp_filename_base, mdrun_mode, '0' if skip else '1']))
        my.run_it(' '.join(['./postproc_flucts.sh', model_name]))

#sbatch -J gmx -p max1n -N 1 --ntasks-per-node=1 --wrap="python run_K_flucts.py -param_ids 7 -omp 12 -mpi 1 -extra_water 30"

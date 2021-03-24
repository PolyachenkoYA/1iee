import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing

import mylib as my

def run_it(cmd, shell=False):
    print(cmd)
    sp.run(cmd, shell=shell)

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

# ================ params =================

N_gpus = 4
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns
compr = 0.0003
time = 10
omp_default = 3
equil_maxsol_poly = [-2.9516, 1117.2]   # maxsol = np.polyval(equil_maxsol_poly, T), [T] = C (not K)
temps = np.array([0.1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55])
P_taus = np.array([200, 400, 800, 1600, 3200, 6400, 12800, 25000, 50000, 100000, 200000, 400000, 800000, 1600000, 3200000, 6400000])
P_taus = np.array([4, 8, 16, 32, 64, 128, 256, 512])
comprs = np.array([2e-4, 3e-4, 4e-4])
times = np.array([20.0, 40.0])
mdrun_mode = 'serial'
mdrun_mode = 'slurm'

# ============== arg parse ====================
N_omp_max = multiprocessing.cpu_count()
possible_gpu_ids = [str(i) for i in range(N_gpus)] + ['-1']
possible_omps = [str(i) for i in range(1, N_omp_max + 1)]
possible_preproc = [no_preproc, preproc_if_permitted, preproc_force, preproc_if_needed]
[omp_cores, mpi_cores, gpu_id, do_mainrun, preproc_mode, param_ids, model_id], _ = \
    my.parse_args(sys.argv[1:], ['-omp', '-mpi', '-gpu_id', '-do_mainrun', '-preproc_mode', '-param_ids', '-id'], \
                  possible_values=[possible_omps, None, possible_gpu_ids, ['no', 'yes'], possible_preproc, None, None], \
                  possible_arg_numbers=[[0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1]], \
                  default_values=[[str(omp_default)], ['0'], ['-1'], ['yes'], [preproc_force], ['0'], ['0']])
param_ids = [int(i) for i in param_ids]
omp_cores = int(omp_cores[0])
mpi_cores = int(mpi_cores[0])
gpu_id = int(gpu_id[0])
do_mainrun = (do_mainrun[0] == 'yes')
preproc_mode = preproc_mode[0]
model_id = model_id[0]
# ===================== cycle ===================
# flucts K

temp = temps[param_ids[0]]
#time = 40.0
#gpu_id = param_ids[0] % N_gpus
#compr = comprs[param_ids[1]]
#time = times[param_ids[2]]
#for Ptau_i, P_tau in enumerate(P_taus[param_ids[3:]]):
for _ in range(1):
    maxsol = int(round(np.polyval(equil_maxsol_poly, temp) * 2))
    nsteps = int(round(time / dt))
    
    model_name = 'flucts_temp' + my.f2str(temp) + '_' + model_id
    mdp_filepath = os.path.join(run_path, model_name, main_mdp_filename_base + '.mdp')
    eql_filepath = os.path.join(run_path, model_name, eql_mdp_filename_base + '.mdp')
    checkpoint_filepath = os.path.join(run_path, model_name, main_mdp_filename_base + '.cpt')
    continue_comp = do_mainrun and os.path.isfile(checkpoint_filepath)
        
    if(preproc_mode in [preproc_if_permitted, preproc_force, preproc_if_needed]):
        if(os.path.isfile(checkpoint_filepath)):
            print('WARNING:\npreproc_mode is ' + preproc_mode + ', but the checkpoint file "' + checkpoint_filepath + '" was found.\nRunning the model from the initial .pdb.')
            if(preproc_mode == preproc_if_permitted):
                skip = (input('Proceed? (y/N)') != 'y')
                print(1)
            elif(preproc_mode == preproc_force):
                skip = False
                print(2)
            elif(preproc_mode == preproc_if_needed):
                skip = True
                print(3)
        else:
            skip = False
    
    else:    # no_preproc
        skip = True 
                        
    if(skip):
        print('Skipping')
    else:
        my.run_it('./clear_restore.sh ' + model_name)
        my.run_it(['python', 'change_mdp.py', '-in', mdp_filepath, '-out', mdp_filepath, '-flds', 'ref-t', str(temp + T_C2K), \
                                                                                                  'nsteps', str(nsteps), \
                                                                                                  'gen-temp', str(temp + T_C2K), \
                                                                                                  'gen-seed', str(model_id)])
        my.run_it(' '.join(['./preproc.sh', model_name, str(omp_cores), str(mpi_cores), str(gpu_id), '1', '1', '2', str(maxsol), '1iee112_prot4gmx.pdb']))
        my.run_it(' '.join(['./mainrun.sh', model_name, str(omp_cores), str(mpi_cores), str(gpu_id), '1iee_wions', minimE_filename_base, mdrun_mode]))
        my.run_it(' '.join(['./save_initial_nojump.sh', model_name, minimE_filename_base]))

    if(do_mainrun):
        my.run_it(' '.join(['./mainrun.sh', model_name, str(omp_cores), str(mpi_cores), str(gpu_id), minimE_filename_base, main_mdp_filename_base, mdrun_mode]))
        my.run_it(' '.join(['./postproc_flucts.sh', model_name]))

import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing

import mylib as my

T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns

# =============== paths ================
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'water_cube')
res_path = os.path.join(root_path, 'res')
main_mdp_filename_base = 'npt'
no_preproc = 'no'
preproc_if_permitted = 'ask'
preproc_force = 'force'

# ================ params =================

times = [0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5]
times = [0.2500, 0.3053, 0.3727, 0.4551, 0.5558, 0.6786, 0.8286, 1.012, 1.235, 1.5085, 1.842, 2.249, 2.746, 3.3535, 4.095, 5.000]
times = [0.25, 0.5, 1.0, 2.0]
temps = [0.1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]
temps = [0, 27, 77]

# ============== arg parse ====================
params = temps
N_gpus = 4
N_omp_max = multiprocessing.cpu_count()
possible_gpu_ids = [str(i) for i in range(N_gpus)] + ['-1']
possible_omps = [str(i) for i in range(1, N_omp_max + 1)]
possible_param_ids = [str(i) for i in range(len(params))] + ['all']
possible_paramids_numbers = range(1, len(params) + 1)
possible_preproc = [no_preproc, preproc_if_permitted, preproc_force]
[omp_cores, mpi_cores, gpu_id, do_mainrun, preproc_mode, param_ids, model_id], _ = \
    my.parse_args(sys.argv[1:], ['-omp', '-mpi', '-gpu_id', '-do_mainrun', '-preproc_mode', '-param_ids', '-id'], \
                  possible_values=[possible_omps, None, possible_gpu_ids, ['no', 'yes'], possible_preproc, possible_param_ids, None], \
                  possible_arg_numbers=[[0, 1], [0, 1], [0, 1], [0, 1], [0, 1], possible_paramids_numbers, [1, 0]], \
                  default_values=[[str(N_omp_max)], ['1'], ['0'], ['yes'], [preproc_force], None, ['0']])
if(param_ids[0] == 'all'):
    param_ids = [str(i) for i in range(1, len(params) + 1)]
param_ids = [int(i) for i in param_ids]
omp_cores = int(omp_cores[0])
mpi_cores = int(mpi_cores[0])
gpu_id = int(gpu_id[0])
do_mainrun = (do_mainrun[0] == 'yes')
preproc_mode = preproc_mode[0]
# ===================== cycle ===================
# flucts K
for param_i in param_ids:
    time = 2.0
    gpu_id = param_i % N_gpus
    Tmp = temps[param_i]
    nsteps = np.intc(time / dt)
    model_name = os.path.join('watercube_T' + my.f2str(Tmp))
    mdp_filepath = os.path.join(run_path, model_name, main_mdp_filename_base + '.mdp')
    checkpoint_filepath = os.path.join(run_path, model_name, main_mdp_filename_base + '.cpt')
    continue_comp = do_mainrun and os.path.isfile(checkpoint_filepath)
    
    if(preproc_mode in [preproc_if_permitted, preproc_force]):
        if(os.path.isfile(checkpoint_filepath)):
            print('WARNING:\npreproc_mode is ' + preproc_mode + ', but the checkpoint file "' + checkpoint_filepath + '" was found.\nRunning the model from the initial .pdb.')
            if(preproc_mode == preproc_if_permitted):
                skip = (input('Proceed? (y/N)') != 'y')
            elif(preproc_mode == preproc_force):
                skip = False
                
            if(skip):
                print('Skipping')
                continue
        my.run_it('./clear_restore.sh ' + model_name)
        my.run_it(['python', 'change_mdp.py', '-in', mdp_filepath, '-out', mdp_filepath, '-flds', 'ref-t', str(Tmp + T_C2K), 'nsteps', str(nsteps)])
        my.run_it(' '.join(['./preproc.sh', model_name, str(omp_cores), '1', str(gpu_id)]))
    
    if(do_mainrun):
        #my.run_it(' '.join(['./mainrun_slurm.sh', model_name, '1', str(mpi_cores), str(gpu_id)]))
        my.run_it(' '.join(['./mainrun_serial.sh', model_name, str(omp_cores), str(gpu_id), main_mdp_filename_base]))
        my.run_it(' '.join(['./postproc_fluct.sh', model_name]))

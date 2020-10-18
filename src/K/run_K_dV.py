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
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'K')
res_path = os.path.join(root_path, 'res')
no_preproc = 'no'
preproc_if_permitted = 'ask'
preproc_if_needed = 'if_needed'
preproc_force = 'force'
preproc_simple = 'simple'
mdp_0_filename_base = 'nvt'
mdp_1_filename_base = 'nvt_bigbox'

# ================ params =================

N_gpus = 4
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns
equil_maxsol_poly = [-2.9516, 1117.2]   # maxsol = np.polyval(equil_maxsol_poly, T), [T] = C (not K)
temps = [0.1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55]

# ============== arg parse ====================
N_omp_max = multiprocessing.cpu_count()
possible_gpu_ids = [str(i) for i in range(N_gpus)] + ['-1']
possible_omps = [str(i) for i in range(1, N_omp_max + 1)]
possible_temp_ids = [str(i) for i in range(len(temps))] + ['all']
possible_tempids_numbers = range(0, len(temps) + 1)
possible_preproc = [no_preproc, preproc_if_permitted, preproc_force, preproc_if_needed, preproc_simple]
[omp_cores, mpi_cores, gpu_id, mainrun_mode, preproc_mode, temp_ids, model_id, dV_mult], _ = \
    my.parse_args(sys.argv[1:], ['-omp', '-mpi', '-gpu_id', '-mainrun_mode', '-preproc_mode', '-temp_ids', '-id', '-dV_mult'], \
                  possible_values=[possible_omps, None, possible_gpu_ids, ['0', '1', '2'], possible_preproc, possible_temp_ids, None, None], \
                  possible_arg_numbers=[[0, 1], [0, 1], [0, 1], [0, 1], [0, 1], possible_tempids_numbers, [1], [0, 1]], \
                  default_values=[[str(N_omp_max)], ['1'], ['0'], ['0'], [preproc_simple], ['all'], None, ['1.001']])
if(temp_ids[0] == 'all'):
    temp_ids = [str(i) for i in range(1, len(temps) + 1)]
temp_ids = [int(i) for i in temp_ids]
omp_cores = int(omp_cores[0])
mpi_cores = int(mpi_cores[0])
gpu_id = int(gpu_id[0])
mainrun_mode = int(mainrun_mode[0])
preproc_mode = preproc_mode[0]
dV_mult = float(dV_mult[0])

# ===================== cycle ===================
for temp_i in temp_ids:
    temp = temps[temp_i]
    maxsol = int(round(np.polyval(equil_maxsol_poly, temp) * 2))

    model_name = os.path.join('dV_temp' + my.f2str(temp) + '_' + model_id)
    mdp_filepath_0 = os.path.join(run_path, model_name, mdp_0_filename_base + '.mdp')
    mdp_filepath_1 = os.path.join(run_path, model_name, mdp_1_filename_base + '.mdp')
    checkpoint_filepath_0 = os.path.join(run_path, model_name, mdp_0_filename_base + '.cpt')
    checkpoint_filepath_1 = os.path.join(run_path, model_name, mdp_1_filename_base + '.cpt')
    
    if(preproc_mode in [preproc_if_permitted, preproc_force, preproc_if_needed]):
        if(mainrun_mode == 1):
            ask_if_preproc = os.path.isfile(checkpoint_filepath_0)
            found_checkpoint_filepath = checkpoint_filepath_0
        elif(mainrun_mode == 2):
            found_checkpoint_filepath = []
            for cp_fp in [checkpoint_filepath_0, checkpoint_filepath_1]:
                if(os.path.isfile(cp_fp)):
                    found_checkpoint_filepath.append(cp_fp)
            ask_if_preproc = bool(found_checkpoint_filepath)
        else:  # mainrun_mode == 0
            ask_if_preproc = False
            
        if(ask_if_preproc):
            print('WARNING:\npreproc_mode is "' + preproc_mode + '", but the checkpoint file "' + str(found_checkpoint_filepath) + '" was found.')
            if(preproc_mode == preproc_if_permitted):
                skip = (input('Proceed from the initial .pdb? (y/N)') != 'y')
            elif(preproc_mode == preproc_force):
                skip = False
            elif(preproc_mode == preproc_if_needed):
                skip = True
        else:
            skip = False
    elif(preproc_mode in [preproc_simple]):
        if(mainrun_mode == 1):
            skip = False
        elif(mainrun_mode == 2):
            skip = True
    else:
        skip = True
                
    if(skip):
        print('Skipping preproc')
    else:
        print('Running the model from the initial .pdb')
        my.run_it('./clear_restore.sh ' + model_name)
        my.run_it(['python', 'change_mdp.py', '-in', mdp_filepath_0, '-out', mdp_filepath_0, '-flds', 'ref-t', str(temp + T_C2K)])
        my.run_it(['python', 'change_mdp.py', '-in', mdp_filepath_1, '-out', mdp_filepath_1, '-flds', 'ref-t', str(temp + T_C2K)])
        my.run_it(' '.join(['./preproc.sh', model_name, str(omp_cores), '1', str(gpu_id), '1', '1', '2', str(maxsol), '1iee112_prot4gmx.pdb']))
    
    if(mainrun_mode == 1):
        #my.run_it(' '.join(['./mainrun_slurm.sh', model_name, '1', str(mpi_cores), str(gpu_id),  mdp_0_filename_base]))
        my.run_it(' '.join(['./mainrun_serial.sh', model_name, str(omp_cores), '1', str(gpu_id),  mdp_0_filename_base]))
        my.run_it(' '.join(['./change_box_size.sh', model_name, mdp_0_filename_base + '.gro', mdp_1_filename_base + '.gro', str(dV_mult)]))
        
#    if(mainrun_mode == 2):
        #my.run_it(' '.join(['./mainrun_slurm.sh', model_name, '1', str(mpi_cores), str(gpu_id),  mdp_1_filename_base]))
        my.run_it(' '.join(['./mainrun_serial.sh', model_name, str(omp_cores), '1', str(gpu_id),  mdp_1_filename_base]))
        my.run_it(' '.join(['./postproc_dV.sh', model_name]))

/
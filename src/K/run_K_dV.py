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

# dV: dL/L = 1e-2.5
N_gpus = 4
omp_default = 3
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns
equil_maxsol_poly = [-2.9516, 1117.2]   # maxsol = np.polyval(equil_maxsol_poly, T), [T] = C (not K)
supercell = np.array([1, 1, 2], dtype=np.intc)
temps = np.array([0.1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55])
times = np.array([20, 40])
Ttaus = np.array([4, 8, 16, 32, 64, 128, 256])

temp = 35.0
time = 20.0
Ttau = 128.0

# ============== arg parse ==========================
supercell_str = ''.join([str(x) for x in supercell])
initial_pdb_filename = '1iee' + supercell_str + '_prot4gmx.pdb'
N_omp_max = multiprocessing.cpu_count()

possible_gpu_ids = [str(i) for i in range(N_gpus)] + ['-1']
possible_omps = [str(i) for i in range(1, N_omp_max + 1)]
possible_preproc = [no_preproc, preproc_if_permitted, preproc_force, preproc_if_needed, preproc_simple]
[omp_cores, mpi_cores, gpu_id, mainrun_mode, preproc_mode, param_ids, model_id, dV_mult], _ = \
    my.parse_args(sys.argv[1:], ['-omp', '-mpi', '-gpu_id', '-mainrun_mode', '-preproc_mode', '-param_ids', '-id', '-dV_mult'], \
                  possible_values=[possible_omps, None, possible_gpu_ids, ['0', '1', '2'], possible_preproc, None, None, None], \
                  possible_arg_numbers=[[0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1], [0, 1]], \
                  default_values=[[str(omp_default)], ['1'], ['0'], ['1'], [preproc_force], ['0'], ['0'], ['2.5']])
param_ids = [int(i) for i  in param_ids]
omp_cores = int(omp_cores[0])
mpi_cores = int(mpi_cores[0])
gpu_id = int(gpu_id[0])
mainrun_mode = int(mainrun_mode[0])
preproc_mode = preproc_mode[0]
model_id = int(model_id[0])
dV_mult = float(dV_mult[0])   # log10(V_new/V_old - 1)

temp = temps[param_ids[0]]
#time = times[param_ids[1]]
#Ttau = Ttaus[param_ids[0] + 2]
gpu_id = param_ids[0] % N_gpus
# ===================== cycle ===================
#for temp_i in temp_ids:
#for dV_mult in [4, 2.5, 3.5]:
for _ in range(1):
    maxsol = int(round(np.polyval(equil_maxsol_poly, temp) * np.prod(supercell)))
    nsteps = int(round(time / dt))

    model_name = 'dV_temp' + my.f2str(temp)
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

        my.run_it(['python', 'change_mdp.py', '-in', mdp_filepath_0, '-out', mdp_filepath_0, '-flds', 'ref-t', str(temp + T_C2K), \
                                                                                                      'nsteps', str(nsteps), \
                                                                                                      'gen-temp', str(temp + T_C2K), \
                                                                                                      'gen-seed', str(model_id), \
                                                                                                      'tau-t', str(Ttau)])
        my.run_it(['python', 'change_mdp.py', '-in', mdp_filepath_1, '-out', mdp_filepath_1, '-flds', 'ref-t', str(temp + T_C2K), \
                                                                                                      'nsteps', str(nsteps), \
                                                                                                      'tau-t', str(Ttau)])

        my.run_it(' '.join(['./preproc.sh', model_name, str(omp_cores), '1', str(gpu_id), supercell_str[0], supercell_str[1], supercell_str[2], str(maxsol), initial_pdb_filename]))
    
    if(mainrun_mode == 1):
        #my.run_it(' '.join(['./mainrun_slurm.sh', model_name, '1', str(mpi_cores), str(gpu_id),  mdp_0_filename_base]))
        my.run_it(' '.join(['./mainrun_serial.sh', model_name, str(omp_cores), '1', str(gpu_id),  mdp_0_filename_base]))
        V_new = 1 + 10**(-dV_mult)
        my.run_it(' '.join(['./change_box_size.sh', model_name, mdp_0_filename_base + '.gro', mdp_1_filename_base + '.gro', str(V_new)]))
        
#    if(mainrun_mode == 2):
        #my.run_it(' '.join(['./mainrun_slurm.sh', model_name, '1', str(mpi_cores), str(gpu_id),  mdp_1_filename_base]))
        my.run_it(' '.join(['./mainrun_serial.sh', model_name, str(omp_cores), '1', str(gpu_id),  mdp_1_filename_base]))
        my.run_it(' '.join(['./postproc_dV.sh', model_name]))



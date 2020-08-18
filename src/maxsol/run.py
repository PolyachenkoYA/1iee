import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing

import mylib as my

def run_it(cmd, shell=False):
    print(cmd)
    sp.run(cmd, shell=shell)
    
#omp_cores = 6
#gpu_id = -1
T_C2K = 273.15
    

# =============== paths ================
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'maxsol')
res_path = os.path.join(root_path, 'res')
nvtmdp_filename = 'nvt.mdp'

# ================ params =================
# appromixate maxsol : P(maxsol) = 1 bar
# 1127, 1096, 1080, 1065, 1050, 1036, 1022, 1010
jobs = np.array([[1, 1, 1], [1, 1, 2], [1, 1, 3], [1, 1, 4], [1, 2, 2]])
temperature_arr = np.array([1, 10, 15, 20, 25, 30, 35, 40])

maxsol_arr = {1:  [1112, 1122, 1127, 1132, 1142],
              10: [1081, 1091, 1096, 1101, 1111],
              15: [1065, 1075, 1080, 1085, 1095],
              20: [1050, 1060, 1065, 1070, 1080],
              25: [1035, 1045, 1050, 1055, 1065],
              30: [1021, 1031, 1036, 1041, 1051],
              35: [1007, 1017, 1022, 1027, 1037],
              40: [995 , 1005, 1010, 1015, 1025]}

maxsol_arr = { 1: [1132, 1142],
              10: [1101, 1111],
              15: [1095],
              20: [1070, 1080],
              25: [1055, 1065],
              30: [1041, 1051],
              35: [1027, 1032],
              40: [1015, 1025]}

#maxsol_arr = {20:[1050]}

# ============== arg parse ====================
N_gpus = 1
N_omp_max = 12
possible_gpu_ids = [str(i) for i in range(N_gpus)]
possible_omps = [str(i) for i in range(1, N_omp_max + 1)]
possible_temp_ids = [str(i) for i in range(len(temperature_arr))]
possible_temp_nums = range(len(temperature_arr) + 1)
possible_jobs_ids = [str(i) for i in range(len(jobs))]
possible_jobs_nums = range(len(jobs) + 1)
[jobs_ids, temperature_ids, omp_cores, gpu_id, mainrun_only], _ = \
    my.parse_args(sys.argv[1:], ['-job_ids', '-temp_ids', '-omp', '-gpu_id', '-continue'], \
                  possible_values=[possible_jobs_ids, possible_temp_ids, possible_omps, possible_gpu_ids, ['0', '1']], \
                  possible_arg_numbers=[possible_jobs_nums, possible_temp_nums, [0, 1], [0, 1], [0, 1]], \
                  default_values=[possible_jobs_ids, possible_temp_ids, [str(multiprocessing.cpu_count())], ['-1'], ['0']])
jobs_ids = [int(i) for i in jobs_ids]
temperature_ids = [int(i) for i in temperature_ids]
omp_cores = int(omp_cores[0])
gpu_id = int(gpu_id[0])
mainrun_only = (mainrun_only[0] == '1')

# ===================== cycle ===================
for job_i in jobs_ids:
    job = jobs[job_i]
    job_strs = [str(x) for x in job]
    job_str = ''.join(job_strs)
    start_pdb_file = '1iee' + job_str + '_prot4gmx.pdb'
    
    for t_i in temperature_ids:
        t = temperature_arr[t_i]
        
        for ms in maxsol_arr[t]:
            path = os.path.join('job' + job_str, 't' + str(t), 'maxsol' + str(ms))
            mdp_filepath = os.path.join(run_path, path, nvtmdp_filename)
            checkpoint_filepath = os.path.join(run_path, path, 'nvt.cpt')
            continue_comp = mainrun_only and os.path.isfile(checkpoint_filepath)
            
            if(not continue_comp):
                if(mainrun_only):
                    print('WARNING:\nmainrun_only is set to be True, but the checkpoint file "' + checkpoint_filepath + '" was not found.\nRunning the model from the initial .pdb.')
                    #goon = input('WARNING:\nmainrun_only is set to be True, but the checkpoint file "' + checkpoint_filepath + '" was not found.\nRunning the model from the initial .pdb. Continue (y/n)?')
                    #if(goon != 'y'):
                    #    print('Skipping')
                    #    continue
                run_it('./clear_restore.sh ' + path, shell=True)
                run_it(['python', 'change_mdp.py', '-in', mdp_filepath, '-out', mdp_filepath, '-flds', 'gen_temp', str(t + T_C2K), 'ref_t', str(t + T_C2K)])
                run_it(' '.join(['./preproc.sh', path] + job_strs + [str(ms * np.prod(job)), start_pdb_file]), shell=True)
                run_it(' '.join(['./minim.sh', path, str(omp_cores), str(gpu_id)]), shell=True)
                
            run_it(' '.join(['./nvt.sh', path, str(omp_cores), str(gpu_id), '1' if(continue_comp) else '0']), shell=True)

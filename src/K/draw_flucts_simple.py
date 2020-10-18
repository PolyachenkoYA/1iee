import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing
import matplotlib.pyplot as plt 

import mylib as my

def run_it(cmd, shell=False):
    print(cmd)
    sp.run(cmd, shell=shell)
    

# =============== paths ================
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'K')
res_path = os.path.join(root_path, 'res')
main_mdp_filename_base = 'npt'
K_filename = 'K.dat'

# ============================================ Ptau, compr, traj_len =================================
# ============================== arg parse ============================
if(not len(args) in [0]):
    print('usage:\n' + sys.argv[0])
    exit(1)

# ============================= params =================================
dV_mult = 1.001
temp_i = 7
times = [5.0, 10.0, 40.0]
times = [5.0, 10.0]
id = 0
P_taus = np.array([10, 20, 50, 100, 200, 400, 800, 1600], dtype=np.float)
comprs = np.array([1e-9, 2e-9, 4e-9, 1e-8])
N_P = len(P_taus)

# ========================= post proc all ==========================

N_compr = len(comprs)
N_t = len(times)
K_arr = np.empty((N_t, N_P, N_compr))
for time_i, time in enumerate(times):
    for P_i, P_tau in enumerate(P_taus):
        for compr_i, compr in enumerate(comprs):
            model_path = os.path.join(run_path, 'flucts_Ptau' + str(P_tau) + '_compr' + str(compr) + '_time' + str(time) + '_' + str(id))
            Kdat_filepath = os.path.join(model_path, 'K.dat')
            if(os.path.isfile(Kdat_filepath)):
                if(os.stat(Kdat_filepath).st_size > 0):
                    K = np.loadtxt(Kdat_filepath)
                else:
                    print('file "' + Kdat_filepath + '" is empty')
                    K = 0
            else:
                K = None
                print('no "' + Kdat_filepath + '" found')
            K_arr[time_i, P_i, compr_i] = K
print(K_arr)

for time_i, time in enumerate(times):
        for compr_i, compr in enumerate(comprs):
            fic, ax = my.get_fig(r'$\tau_P$', r'K (Pa)', r'K($\tau_P$)| time (ns) = ' + str(time) + ', compr = ' + str(compr) + ', T=$35C^\circ$', xscl='log')
            K_draw = K_arr[time_i, :, compr_i]
            draw_ind = ~np.equal(K_draw, None)
            ax.plot(P_taus[draw_ind], K_draw[draw_ind])
plt.show()
exit(0)

# ============================================== ?? cell size convergence ?? =========================================
#for job_name in job_names:
#    for stab_time in stab_times:
#        txt_datafile_path = os.path.join(res_path, job_name + '_' + str(stab_time) + '.txt')
#        if(os.path.isfile(txt_datafile_path)):
#            os.remove(txt_datafile_path)
#        for t in t_arr:
#            run_it('python stats.py -dir ' + os.path.join(job_name, 't' + str(t)) + ' -mode short save -feat P -stab_time ' + str(stab_time), shell=True)


# ============================================ traj len analysis =================================
# ================ params =================

T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns
ids = [0, 1, 2, 3, 5, 6, 7, 8, 9, 10, 11, 12]

# ============== process ====================
N = len(ids)
N_t = len(times)
K_arr = np.zeros((N_t, N))
K = np.zeros(N_t)
d_K = np.zeros(N_t)
for t_i, time in enumerate(times):
    for i, id in enumerate(ids):
        suff = my.f2str(time) + '_' + str(id)
        model_path = os.path.join(run_path, 'time' + suff)
        K_filepath = os.path.join(model_path, K_filename)
        K_arr[t_i, i] = np.loadtxt(K_filepath)
        print(K_arr[t_i, i], end=' ')
        #sp.run(['cp', K_filepath, './K_' + suff + '.dat'])
        #print()
    print(time)
    #K[t_i] = np.mean(K_arr[t_i, :])
    #d_K[t_i] = np.std(K_arr[t_i, :])

print(K)
print(d_K)
#fig, ax = my.get_fig('traj length (ns)', 'K', title='K(t)')
#ax.plot(times, K_arr)
#ax.errorbar(times, K, yerr=d_K)

#plt.show()

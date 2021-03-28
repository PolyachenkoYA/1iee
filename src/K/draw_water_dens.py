import os
import sys
import numpy as np
import subprocess as sp
import multiprocessing
import matplotlib.pyplot as plt 
from struct import unpack

import mylib as my

import gromacs.formats as gmx

# =============== paths ================
#root_path = my.git_root_path()
with open('git_root_path') as f:
    root_path = f.readline()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'K')
res_path = os.path.join(root_path, 'res')
npt_filename = 'npt'
nvt0_filename = 'nvt'
nvt1_filename = 'nvt_bigbox'

def read_last_line(filename):
    with open(filename, 'rb') as f:
        f.seek(-2, os.SEEK_END)
        while f.read(1) != b'\n':
            f.seek(-2, os.SEEK_CUR)
        last_line = f.readline().decode()
        return last_line
    
def load_xvg(filepath, failsafe=True):
    if(failsafe):
        if(not os.path.isfile(filepath)):
            print('file "' + filepath + '" was not found. Aborting.')
            return None
    xvg_file = gmx.XVG()
    xvg_file.read(filepath)
    return xvg_file

def get_frames_number(filepath):
    with open(filepath, 'rb') as f:
        tag = f.read(8) # Tag: magic number and number of atoms
        frame_size = 92 + magic_for_xtc_size(f.read(84)[-4:]) # Size of frame in bytes
        file_size = os.stat(filepath).st_size
        return file_size / frame_size
    
# ================ general params =================
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns
equil_maxsol_poly = [-2.9516, 1117.2]   # maxsol = np.polyval(equil_maxsol_poly, T), [T] = C (not K)

temp = 30
model_id = 0

model_name = 'flucts_temp' + my.f2str(temp) + '_' + str(model_id)
model_path = os.path.join(run_path, model_name)
filepath = os.path.join(model_path, 'water_dens.xvg')

xvg_file = load_xvg(filepath)
#print(xvg_file)
z = xvg_file.array[0]
rho = xvg_file.array[1]
rho_0 = 10**5 * 0.018 / 8.31 / 300

#fig, ax = my.get_fig('$z$ (nm)', '$\rho$ (kg/$m^3$)', title='$\rho(z)$')
fig, ax = my.get_fig('z (nm)', 'rho_water (kg/m3)', title='rho_w(z)', yscl='log')

ax.plot(z, rho, label='water')
ax.plot([min(z), max(z)], [rho_0] * 2, label='saturated vapor')

ax.legend()

plt.savefig('water_dens.png')
plt.show()

print('DONE')

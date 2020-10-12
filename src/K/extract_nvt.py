import os
import sys
import re
import numpy as np
import subprocess
import multiprocessing
import matplotlib.pyplot as plt 

import mylib as my
    
T_C2K = 273.15
dt = 2e-6    # 1 fs = 1e-6 ns

# =============== paths ================
root_path = my.git_root_path()
run_path = os.path.join(root_path, 'run')
exe_path = os.path.join(root_path, 'src', 'water_cube')
res_path = os.path.join(root_path, 'res')

# ============ arg parse ===============
args = sys.argv[1:]
argc = len(args)
if(argc not in [2]):
	print('usage:\n' + sys.argv[0] + '   input_filename   output_filename')
	exit(1)
nvt_filename = args[0]
nvt_pressure_filename = args[1]

# =============== process ==============
num_re = r'-?[0-9]+.?[0-9]*'
pressure_file = open(nvt_pressure_filename, 'w')
with open(nvt_filename) as nvt_file:
	for line in nvt_file:
		match_P = re.search(r'Pressure\ +(' + num_re + '\ +' + num_re + '\ +' + num_re + '\ +' + num_re + ')', line)
		if(match_P):
			print(match_P.group(1), file=pressure_file)
		match_T = re.search(r'Temperature\ +(' + num_re + '\ +' + num_re + '\ +' + num_re + '\ +' + num_re + ')', line)
		if(match_T):
			print(match_T.group(1), file=pressure_file)
pressure_file.close()


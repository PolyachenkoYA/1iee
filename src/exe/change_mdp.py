import os
import sys
import subprocess as sp
import numpy as np
import mylib as my
import re

def print_usage_and_exit(exe_name=sys.argv[0], exit_code=1):
    print('usage:\n' + exe_name + '   ' + infile_flag + ' infile   ' + outfile_flag + ' outfile   ' + fields_flag + ' field1 value1 [f2 v2 ...]')
    exit(exit_code)

# ============= consts ==============
infile_flag = '-in'
outfile_flag = '-out'
fields_flag = '-flds'

# ============= arg parse ============
args = sys.argv[1:]
argc = len(args)
if(np.mod(argc, 2) != 1):
    print_usage_and_exit()

[infile_name, outfile_name, fields] = my.parse_args(args, [infile_flag, outfile_flag, fields_flag])
if(np.any([(len(i) == 0) for i in [infile_name, outfile_name, fields]])):
    print_usage_and_exit()
if(np.mod(len(fields), 2) != 0):
    print_usage_and_exit()

fields = np.array(fields).reshape((len(fields) // 2), 2)

# ========== proc ============
infile = open(infile_name, 'r')
infile_textdata = infile.readlines()
infile.close()
outfile = open(outfile_name, 'w')
fields_found = [False] * len(fields)
for line in infile_textdata:
    new_line = line
    for field_i, field in enumerate(fields):
        ptrn = r'(' + field[0] + r'\s+=\s).*(\s+;.*)'
        match = re.search(ptrn, new_line)
        if(match):
            new_line = match.group(1) + field[1] + match.group(2) + '\n'
            fields_found[field_i] = True
    print(new_line, end='', file=outfile)
outfile.close()

for i, f in enumerate(fields_found):
    if(not f):
        print('field "' + fields[i][0] + '" not found in "' + infile_name + '"', file=sys.stderr)

print('file "' + outfile_name + '" written')
import os
import sys
import re

import mylib as my

[topol_filepath, N], _ = \
     my.parse_args(sys.argv[1:], ['-file', '-N'], \
                   possible_values=[None, None], \
                   possible_arg_numbers=[[1], [1]], \
                   default_values=[None, None])
N = int(N)

with open(topol_filepath, 'r') as f:
	topol_strs = f.readlines()

molecules_partition_header = '[ molecules ]'
molecule_line_pattern = r'([a-zA-Z_0-9]+)(\s+)([1-9][0-9]*)'

molecules_partition_ind = topol_strs.index(molecules_partition_header + '\n') + 2
molecules_strs = topol_strs[molecules_partition_ind:]
topol_strs = topol_strs + molecules_strs * (N - 1)
open(topol_filepath, 'w').writelines(topol_strs)

#open(topol_filepath, 'w').writelines(
#    topol_strs[: molecules_partition_ind] + 
#    [
#        re.sub(molecule_line_pattern, r'\g<1>\g<2>' + str(int(re.match(molecule_line_pattern, s).group(3)) * N), s)
#        for s in topol_strs[molecules_partition_ind : ]
#    ]
#)

#for str_i in range(molecules_partition_ind, len(topol_strs)):
#	m = re.match(molecule_line_pattern, topol_strs[str_i])
#	topol_strs[str_i] = re.sub(molecule_line_pattern, \
#	       r'\g<1>\g<2>' + str(int(m.group(3)) * N), \
#	       topol_strs[str_i])
#open(topol_filepath, 'w').writelines(topol_strs)

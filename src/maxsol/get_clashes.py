import sys
import os
import re

min_dist = 0.01

def get_res_id(word):
    m = re.search(r'#(\d+):(\d+)\.', word)
    return (int(m.group(1)), int(m.group(2)))
    
def add_id(clashes, word, max_n):
    (mdl, num) = get_res_id(word)
    if(mdl < max_n):
        if(not (num in clashes[mdl])):
            clashes[mdl].append(num)
    return clashes
    
def print_clash(mdl, clashes):
    print('delete #' + str(mdl) + ':' + ",".join([str(cl) for cl in clashes]))

# ================================================

narg = len(sys.argv)
if(not (narg in [3])):
    print('usage:\n' + sys.argv[0] + '   filename   max_model_num')
    exit(1)
filename = sys.argv[1]
max_model_n = int(sys.argv[2])

# ==========

file = open(filename, "r")

clashes = [[]] * max_model_n   # max_model_n empty lists
for line in file.readlines():
    if(line[0] == '#'):
        words = line.split()
        dist = float(words[3])
        if(dist < min_dist):
            clashes = add_id(clashes, words[0], max_model_n)
            clashes = add_id(clashes, words[1], max_model_n)
            
for i in range(len(clashes)):
    clashes[i].sort()
    print_clash(i, clashes[i])


file.close()

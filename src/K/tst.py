import gromacs.formats as gmx
import os
import pickle

def load_xvg(filepath, failsafe=2):
    xvg_file = gmx.XVG()
    xvg_file.read(filepath)
    return xvg_file

def test_fnc(xvg_file):
    print('names: ', xvg_file.names)
    print('shape: ', xvg_file.array.shape)
    print('names after: ', xvg_file.names)    
    print('shape after: ', xvg_file.array.shape)
    
xvg_filepath = 'test.txt'
pickle_filepath = 'test.pkl'

xvg_file = load_xvg(xvg_filepath)
print('names out: ', xvg_file.names)
print('shape out: ', xvg_file.array.shape)
test_fnc(xvg_file)

pickle.dump(xvg_file, open(pickle_filepath, 'wb'))
xvg_file = pickle.load(open(pickle_filepath, 'rb'))
print('\n2 out: ', xvg_file.names)
test_fnc(xvg_file)

os.remove(pickle_filepath)

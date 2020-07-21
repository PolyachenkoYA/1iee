# conda activate moleculekit
### it has py3.6 and htmd installed

from moleculekit.molecule import Molecule
from moleculekit.tools.preparation import proteinPrepare
import sys

args = sys.argv[1:]
if(len(args) != 1):
	print('usage:\n' + sys.argv[0] + '   infile.pdb')

filename = args[0]

mol = Molecule(filename)

mol.filter('protein or water')

molPrep, prepData = proteinPrepare(mol, pH=4.5, returnDetails=True)

print(molPrep)

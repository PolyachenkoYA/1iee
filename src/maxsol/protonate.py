#conda create -n moleculekit python=3.6
#conda activate moleculekit
#
#conda install moleculekit -c acellera
#conda install htmd-pdb2pqr -c acellera
#conda install pandas=0.25.3

from moleculekit.molecule import Molecule
from moleculekit.tools.preparation import proteinPrepare
import sys

args = sys.argv[1:]
if(not(len(args) in [1, 2])):
	print('usage:\n' + sys.argv[0] + '   infile.pdb   [outfile.pdb]')

input_filename = args[0]
output_filename = input_filename if(len(args) == 1) else args[1]
if(input_filename == output_filename):
	if(input(input_filename + ' will be overwritten. Continue? (y/[n])\n') != 'y'):
		print('Aborted')
		exit(1)

molecule = Molecule(input_filename)
molecule.filter('protein or water')
molecule_prepared, prep_data = proteinPrepare(molecule, pH=4.5, returnDetails=True)

molecule_prepared.write(output_filename)
print(output_filename + ' written')

### '###' - human actions/comments
### '#'   - commands to execute for gmx/chimera
### =================================================

set -e
gmx_exeutable=gmx_mpi

argc=$#
if [ $argc -ne 2 ] && [ $argc -ne 3 ]
then
        printf "usage:\n$0   job_name   ompN   [maxsol]\n"
        exit 1
fi
if [ $argc == 3 ]
then
	maxsol=$3
else
	maxsol=0
fi
job_id=$1
root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/exe
omp=$2
cd $root_path
cd $run_path

### Download 1iee.pdb from database

# chimera --start "Unit Cell" --start "Command Line" --start "Reply Log" pdbs/1iee.pdb

### Tools -> HigherOrderStructure -> UnitCell -> outline

### Tools -> HigherOrderStructure -> UnitCell -> make copies

# findclash # test self colorClashes true setAttrs true log true namingStyle command overlapCutoff 2 selectClashes true saveFile clashesALL.txt
### Options-> Number of sells 3 3 3, offset -1 -1 -1
### there will be a number of lines in Reply Log. They have 'dist' field. Some of them have dist >~0.2. Some have <0.01. Use those with <0.01 in the next step.
# python exe/get_chashes.py clashesALL.txt 8
### get all water-residues (molecules) to delete duw to overlap

### these are residues numbers (HOH molecule's numbers in case of solvent) given by the findclash
# delete #0-7:1092,1103,1139,1169

# combine # newchainids true name 1iee_cryst close true

### write the whole system (all 0-7 groups were encorporated into a single group 0 by 'combine' command above) relative to the outline (group 8)
# write relative 0 8 pdbs/1iee_cryst.pdb


# python ../src/exe/protonate.py 1iee_cryst.pdb 1iee_prot.pdb
### this can be done on playmolecule/ProteinPrepare website

# ../src/exe/playmol2gmx.sh 1iee_prot.pdb 1iee_prot4gmx.pdb

$gmx_exeutable pdb2gmx -f 1iee"$job_id"_prot4gmx.pdb -o 1iee_init.gro -water tip4p -missing < protonation_gromacs.in
$gmx_exeutable editconf -f 1iee_init.gro -o 1iee_newbox.gro -c -box 7.7061   7.7061   3.7223
$gmx_exeutable grompp -f ions.mdp -c 1iee_newbox.gro -p topol.top -o ions.tpr -maxwarn 5
$gmx_exeutable genion -s ions.tpr -o 1iee_wions.gro -p topol.top -pname NA -nname CL -neutral -rmin 0.28  < genion_gromacs.in
$gmx_exeutable solvate -cp 1iee_wions.gro -cs tip4p.gro -o 1iee_solv.gro -p topol.top -maxsol $maxsol

cd $root_path
cd $exe_path
./mainrun.sh $job_id $omp

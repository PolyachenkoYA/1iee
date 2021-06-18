#!/bin/bash
### '###' - human actions/comments
### '#'   - commands to execute for gmx/chimera
### =================================================

set -e
gmx_serial=gmx_mpi
gmx_serial=$HOME/gromacs-2020/build_lin/bin/gmx
gmx_serial=gmx_ser_20

gmx_mdrun=gmx_mpi
gmx_mdrun=$HOME/gromacs-2020/build/bin/gmx_mpi
#gmx_mdrun=$gmx_serial

argc=$#
if [ $argc -ne 8 ] && [ $argc -ne 9 ]
then
        printf "usage:\n$0   job_name   ompN   mpiN   gpu_id   [nx ny nz   [maxsol   [start_pdbf_file]]]\n"
        exit 1
fi

job_id=$1
omp=$2
mpiN=$3
gpu_id=$4
nx=$5
ny=$6
nz=$7
lx=$(echo 7.7061*$nx | bc)
ly=$(echo 7.7061*$ny | bc)
lz=$(echo 3.7223*$nz | bc)
maxsol=$8
if [ $argc -gt 8 ]
then
        start_pdb_file=$9
else
        start_pdb_file=1iee_prot4gmx.pdb
fi
#root_path=$(git rev-parse --show-toplevel)
root_path=`cat git_root_path`
run_path=run/$job_id
exe_path=src/K
topol_filename=topol.top
cd $root_path
cd $run_path

### Download 1iee.pdb from database

# chimera --start "Unit Cell" --start "Command Line" --start "Reply Log" pdbs/1iee.pdb

### Tools -> HigherOrderStructure -> UnitCell -> outline

### Tools -> HigherOrderStructure -> UnitCell -> make copies

# findclash # test self colorClashes true setAttrs true log true namingStyle command selectClashes true saveFile clashesALL.txt
### Options-> Number of sells 3 3 3, offset -1 -1 -1
### there will be a number of lines in Reply Log. They have 'dist' field. Some of them have dist >~0.2. Some have <0.01. Use those with <0.01 in the next step.
# python exe/get_chashes.py clashesALL.txt 8
### get all water-residues (molecules) deleted due to overlap

### these are residues numbers (HOH molecule's numbers in case of solvent) given by the findclash
# delete #0-7:1092,1103,1139,1169
# for 112, 113, 114
# delete #0-7:1092,1093,1103,1139,1169

# combine # newchainids true name 1iee_cryst close true

### write the whole system (all 0-7 groups were encorporated into a single group 0 by 'combine' command above) relative to the outline (group 8)
# write relative 0 8 1iee_cryst.pdb

# python ../src/exe/protonate.py 1iee_cryst.pdb 1iee_prot.pdb
### this can be done on playmolecule/ProteinPrepare website

# ../src/exe/playmol2gmx.sh 1iee_prot.pdb 1iee_prot4gmx.pdb

$gmx_serial pdb2gmx -f $start_pdb_file -o 1iee_init.gro -missing -p $topol_filename < protonation_gromacs.in
$gmx_serial editconf -f 1iee_init.gro -o 1iee_newbox.gro -c -box $lx $ly $lz
$gmx_serial solvate -cp 1iee_newbox.gro -cs amber03w.ff/tip4p2005.gro -o 1iee_solv.gro -p $topol_filename -maxsol $maxsol
#$gmx_serial solvate -cp 1iee_newbox.gro -cs tip4p.gro                  -o 1iee_solv.gro -p $topol_filename -maxsol $maxsol
$gmx_serial grompp -f ions.mdp -c 1iee_solv.gro -p $topol_filename -o 1iee_wions.tpr -maxwarn 5
$gmx_serial genion -s 1iee_wions.tpr -o 1iee_wions.gro -p $topol_filename -pname NA -nname CL -neutral -rmin 0.2  < genion_gromacs.in

cd $root_path
cd $exe_path

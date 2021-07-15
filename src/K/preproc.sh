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

nx="1"
ny="1"
nz="2"
Nx="1"
Ny="1"
Nz="1"
maxsol="0"
start_pdb_file=1iee_prot4gmx.pdb
do_1phase="0"

argc=$#
if [ $argc -ne 11 ] && [ $argc -ne 12 ] && [ $argc -ne 13 ] && [ $argc -ne 8 ] && [ $argc -ne 5 ]
then
	printf "usage:\n$0   job_name   ompN   mpiN   gpu_id   [Nx Ny Nz($Nx,$Ny,$Nz)   [nx ny nz($nx,$ny,$nz)   [maxsol($maxsol)   [start_pdbf_file($start_pdb_file)   [1phase($do_1phase)]]]]]\n"
	exit 1
fi

job_id=$1
omp=$2
mpiN=$3
gpu_id=$4
if [ $argc -gt 4 ]
then
	Nx=$5
	Ny=$6
	Nz=$7
fi
if [ $argc -gt 7 ]
then
	nx=$8
	ny=$9
	nz=${10}
fi
Nbox=$(echo $Nx*$Ny*$Nz | bc)
lx=$(echo 7.7061*$nx | bc)
ly=$(echo 7.7061*$ny | bc)
lz=$(echo 3.7223*$nz | bc)
Lx=$(echo $lx*$Nx | bc)
Ly=$(echo $ly*$Ny | bc)
Lz=$(echo $lz*$Nz | bc)
Lz_big=$(echo $Lz+25 | bc)
if [ $argc -gt 10 ]
then
	maxsol=${11}
fi
if [ $argc -gt 11 ]
then
	start_pdb_file=${12}
fi
if [ $argc -gt 12 ]
then
	do_1phase=${13}
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

init_gro=1iee_init.gro
smallbox_gro=1iee_smallbox.gro
bigbox_gro=1iee_bigbox.gro
gro_to_ionize=1iee_solv.gro
wions_tpr=1iee_wions.tpr
wions_gro=1iee_wions.gro

$gmx_serial pdb2gmx -f $start_pdb_file -o $init_gro -missing -p $topol_filename < protonation_gromacs.in
$gmx_serial editconf -f $init_gro -o $smallbox_gro -c -box $lx $ly $lz
$gmx_serial solvate -cp $smallbox_gro -cs amber03w.ff/tip4p2005.gro -o $gro_to_ionize -p $topol_filename -maxsol $maxsol   # tip4p/2005
#$gmx_serial solvate -cp $smallbox_gro -cs tip4p.gro                 -o $gro_to_ionize -p $topol_filename -maxsol $maxsol   # tip4p
$gmx_serial grompp -f ions.mdp -c $gro_to_ionize -p $topol_filename -o $wions_tpr -maxwarn 5
$gmx_serial genion -s $wions_tpr -o $wions_gro -p $topol_filename -pname NA -nname CL -neutral -rmin 0.2  < genion_gromacs.in
$gmx_serial genconf -f $wions_gro -nbox $Nx $Ny $Nz -o $bigbox_gro
python $root_path/$exe_path/multiply_topology.py -file $topol_filename -N $Nbox
if [[ "$do_1phase" == "0" ]]
then
	gro_2phase=1iee_2phase.gro
	$gmx_serial editconf -f $bigbox_gro -o $gro_2phase -c -box $Lx $Ly $Lz_big
	mv $gro_2phase $bigbox_gro
fi

cd $root_path
cd $exe_path

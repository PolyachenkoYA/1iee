#!/bin/bash
### '###' - human actions/comments
### '#'   - commands to execute for gmx/chimera
### =================================================

set -e
gmx_serial=gmx_mpi
gmx_serial=gmx_serial
gmx_serial=gmx_ser_gpu

gmx_mdrun=gmx_mpi
gmx_mdrun=gmx_angara
gmx_mdrun="$gmx_serial"

argc=$#
if [ $argc -ne 3 ]
then
        printf "usage:\n$0   job_name   ompN   mpiN\n"
        exit 1
fi
job_id=$1
ompN=$2
mpiN=$3
root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/water_cube

cd $root_path
cd $run_path

$gmx_serial solvate -cs tip4p -o solvated.gro -box 7.70610 7.70610 7.44460 -p topol.top
#$gmx_serial solvate -cs tip4p -o solvated.gro -box 3 3 3 -p topol.top

$gmx_serial grompp -f min.mdp -c solvated.gro -p topol.top -o min.tpr
$gmx_mdrun mdrun -v -deffnm min -ntmpi $mpiN -ntomp $ompN -pin on

$gmx_serial grompp -f min2.mdp -c min.gro -p topol.top -o min2.tpr
$gmx_mdrun mdrun -v -deffnm min2 -ntmpi $mpiN -ntomp $ompN -pin on

$gmx_serial grompp -f eql.mdp -c min2.gro -p topol.top -o eql.tpr
$gmx_mdrun mdrun -deffnm eql -ntmpi $mpiN -ntomp $ompN -pin on

$gmx_serial grompp -f eql2.mdp -c eql.gro -p topol.top -o eql2.tpr
$gmx_mdrun mdrun -deffnm eql2 -ntmpi $mpiN -ntomp $ompN -pin on

cp eql2.gro nvt.gro
#$gmx_serial grompp -f npt.mdp -c eql2.gro -p topol.top -o npt.tpr

cd $root_path
cd $exe_path

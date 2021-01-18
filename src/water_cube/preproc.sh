#!/bin/bash
### '###' - human actions/comments
### '#'   - commands to execute for gmx/chimera
### =================================================

set -e
gmx_serial=gmx_mpi
#gmx_serial=gmx_serial
#gmx_serial=gmx_ser_gpu

gmx_mdrun=gmx_mpi
gmx_mdrun=gmx_angara
gmx_mdrun="$gmx_serial"

argc=$#
if [ $argc -ne 4 ]
then
        printf "usage:\n$0   job_name   ompN   mpiN   gpu_id\n"
        exit 1
fi
job_id=$1
ompN=$2
mpiN=$3
gpu_id=$4
root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/water_cube

cd $root_path
cd $run_path

$gmx_serial solvate -cs tip4p -o solvated.gro -box 7.70610 7.70610 7.44460 -p topol.top
#$gmx_serial solvate -cs tip4p -o solvated.gro -box 3 3 3 -p topol.top

$gmx_serial grompp -f min.mdp -c solvated.gro -p topol.top -o min.tpr
if [ $mpiN -gt 1 ]
then
    $gmx_mdrun mdrun -v -deffnm min -ntmpi $mpiN -ntomp $ompN  -gpu_id $gpu_id
else
    $gmx_mdrun mdrun -v -deffnm min -ntomp $ompN  -gpu_id $gpu_id
fi

$gmx_serial grompp -f min2.mdp -c min.gro -p topol.top -o min2.tpr
if [ $mpiN -gt 1 ]
then
    $gmx_mdrun mdrun -v -deffnm min2 -ntmpi $mpiN -ntomp $ompN  -gpu_id $gpu_id
else
    $gmx_mdrun mdrun -v -deffnm min2 -ntomp $ompN  -gpu_id $gpu_id
fi

$gmx_serial grompp -f eql.mdp -c min2.gro -p topol.top -o eql.tpr
if [ $mpiN -gt 1 ]
then
    $gmx_mdrun mdrun -v -deffnm eql -ntmpi $mpiN -ntomp $ompN  -gpu_id $gpu_id
else
    $gmx_mdrun mdrun -v -deffnm eql -ntomp $ompN  -gpu_id $gpu_id
fi

$gmx_serial grompp -f eql2.mdp -c eql.gro -p topol.top -o eql2.tpr
if [ $mpiN -gt 1 ]
then
    $gmx_mdrun mdrun -v -deffnm eql2 -ntmpi $mpiN -ntomp $ompN  -gpu_id $gpu_id
else
    $gmx_mdrun mdrun -v -deffnm eql2 -ntomp $ompN  -gpu_id $gpu_id
fi

cp eql2.gro nvt.gro
cp eql2.gro npt.gro
#$gmx_serial grompp -f npt.mdp -c eql2.gro -p topol.top -o npt.tpr

cd $root_path
cd $exe_path

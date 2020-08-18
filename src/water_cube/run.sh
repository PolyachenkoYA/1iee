#!/bin/bash
### '###' - human actions/comments
### '#'   - commands to execute for gmx/chimera
### =================================================

set -e
gmx_serial=gmx_ser_gpu

gmx_mdrun=gmx_mpi
gmx_mdrun=gmx_angara
gmx_mdrun=$gmx_serial

argc=$#
if [ $argc -ne 1 ]
then
        printf "usage:\n$0   job_name\n"
        exit 1
fi
job_id=$1
root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/maxsol

./clear_restore.sh $job_id

cd $root_path
cd $run_path

$gmx_serial solvate -cs tip4p -o solvated.gro -box 7.70610 7.70610 7.44460 -p topol.top

$gmx_serial grompp -f min.mdp -c solvated.gro -p topol.top -o min.tpr
$gmx_mdrun mdrun -v -deffnm min

$gmx_serial grompp -f min2.mdp -c min.gro -p topol.top -o min2.tpr
$gmx_mdrun mdrun -v -deffnm min2

$gmx_serial grompp -f eql.mdp -c min2.gro -p topol.top -o eql.tpr
$gmx_mdrun mdrun -deffnm eql

$gmx_serial grompp -f eql2.mdp -c eql.gro -p topol.top -o eql2.tpr
$gmx_mdrun mdrun -deffnm eql2

#$gmx_serial grompp -f prd.mdp -c eql2.gro -p topol.top -o prd.tpr
#$gmx_mdrun mdrun -deffnm prd

$gmx_serial energy -f eql.edr -o eql-temp.xvg < temp.in
$gmx_serial energy -f min.edr -o min-energy.xvg < en1.in
$gmx_serial energy -f min2.edr -o min2-energy.xvg < en2.in

cd $root_path
cd $exe_path

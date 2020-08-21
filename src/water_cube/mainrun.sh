#!/bin/bash

set -e
gmx_serial=gmx_mpi

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
exe_path=src/water_cube

cd $root_path
cd $run_path

$gmx_serial grompp -f prd.mdp -c eql2.gro -p topol.top -o prd.tpr
$gmx_mdrun mdrun -deffnm prd -ntomp 6 -pin on

cd $root_path
cd $exe_path

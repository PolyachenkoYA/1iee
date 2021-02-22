#!/bin/bash

set -e
gmx_serial=gmx_mpi
gmx_serial=gmx_ser_newhead

gmx_mdrun=../../..$HOME/.local/gromacs/bin/gmx_mpi
#gmx_mdrun=gmx_angara
#gmx_mdrun="$gmx_serial"

argc=$#
if [ $argc -ne 5 ]
then
        printf "usage:\n$0   job_name   ompN   mpiN   gpu_id   name\n"
        exit 1
fi
job_id=$1
ompN=$2
mpiN=$3
gpu_id=$4
name=$5
root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/K

cd $root_path
cd $run_path

$gmx_serial grompp -f $name.mdp -c $name.gro -p topol.top -o $name.tpr

if [ $gpu_id -eq -1 ]
then
    trun -m 1 -ppn=$mpiN $gmx_mdrun mdrun -v -deffnm $name -ntomp $ompN -cpi $name.cpt -cpo $name.cpt
else
    trun -m 1 -ppn=$mpiN $gmx_mdrun mdrun -v -deffnm $name -ntomp $ompN -cpi $name.cpt -cpo $name.cpt -gpu_id $gpu_id
fi

cd $root_path
cd $exe_path

#!/bin/bash

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
        printf "usage:\n$0   job_name   ompN   gpu_id   name\n"
        exit 1
fi
job_id=$1
ompN=$2
gpu_id=$3
name=$4
root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/water_cube

cd $root_path
cd $run_path

$gmx_serial grompp -f $name.mdp -c $name.gro -p topol.top -o $name.tpr
#$gmx_serial mdrun -deffnm $name -ntomp $ompN -ntmpi $mpiN -dlb no -pin on -cpi $name.cpt -cpo $name.cpt -maxh 24 -nb gpu -gpu_id $gpu_id
$gmx_mdrun mdrun -v -deffnm $name -ntomp $ompN -cpi $name.cpt -cpo $name.cpt -gpu_id $gpu_id
#sbatch -J gromacs -p max1n -N 1 --reservation=test --ntasks-per-node=$mpiN --gres=gpu:1 --wrap="$cmd"

cd $root_path
cd $exe_path

#!/bin/bash

set -e
gmx_serial=gmx_ser_newhead
#gmx_serial=gmx_serial
#gmx_serial=gmx_ser_gpu

gmx_mdrun=gmx_mpi
#gmx_mdrun=gmx_angara
#gmx_mdrun="$gmx_serial"

argc=$#
if [ $argc -ne 7 ]
then
        printf "usage:\n$0   job_name   ompN   mpiN   gpu_id   in_name   out_name   run_mode\n"
        exit 1
fi
job_id=$1
ompN=$2
mpiN=$3
gpu_id=$4
in_name=$5
name=$6
run_mode=$7
root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/K
timeout=0

cmd="mdrun -v -deffnm $name -ntomp $ompN -cpi $name.cpt -cpo $name.cpt -pin on"
if [ $gpu_id -ge 0 ]
then
    cmd="$cmd -gpu_id $gpu_id"
fi
if [ $timeout -gt 0 ]
then
    cmd="$cmd -maxh $timeout"
fi
if [ $mpiN -gt 0 ]
then
    cmd="$cmd -ntmpi $mpiN"
fi

cd $root_path
cd $run_path

$gmx_serial grompp -f $name.mdp -c $in_name.gro -p topol.top -o $name.tpr -maxwarn 2

if [[ "$run_mode" == "serial" ]]
then
    $gmx_serial $cmd
elif [[ "$run_mode" == "slurm" ]]
then
    srun --ntasks-per-node=$mpiN $gmx_mdrun $cmd
elif [[ "$run_mode" == "trun" ]]
then
    trun -m 1 -ppn=$mpiN $gmx_mdrun $cmd
else
    echo "wrong run_mode: $run_mode"
fi

#sbatch -J gromacs -p max1n -N 1 --reservation=test --ntasks-per-node=$mpiN --gres=gpu:1 --wrap="$cmd"

cd $root_path
cd $exe_path

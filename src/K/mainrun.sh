#!/bin/bash

set -e
gmx_serial=gmx_mpi
gmx_serial=gmx_ser_20

gmx_mdrun=gmx_mpi
gmx_mdrunonly=mdrun_mpionly_20
gmx_mdrun=gmx_mpi2
#gmx_mdrun=gmx_sermpi_20
#gmx_mdrun="$gmx_serial"

argc=$#
if [ $argc -ne 9 ]
then
        printf "usage:\n$0   job_name   ompN   mpiN   gpu_id   in_filename   out_name   run_mode   do_grompp(1/0)   do_dlb(1/0)\n"
        exit 1
fi
job_id=$1
ompN=$2
mpiN=$3
gpu_id=$4
in_filename=$5
name=$6
run_mode=$7
do_grompp=$8
do_dlb=$9
#root_path=$(git rev-parse --show-toplevel)
root_path=`cat git_root_path`
run_path=run/$job_id
exe_path=src/K
timeout=24

cmd="srun --ntasks-per-node=$ompN $gmx_mdrun mdrun -v -deffnm $name -ntomp 1 -cpi $name.cpt -cpo $name.cpt -pin off"
#cmd="srun --ntasks-per-node=1 $gmx_mdrunonly -v -deffnm $name -ntomp $ompN -cpi $name.cpt -cpo $name.cpt -pin off"
if [ $gpu_id -ge -1 ]
then
    cmd="$cmd -nb gpu"
    if [ $gpu_id -ge 0 ]
    then
        cmd="$cmd -gpu_id $gpu_id"
    fi
fi
if [ $do_dlb -eq 1 ]
then
    cmd="$cmd -dlb yes"
else
    cmd="$cmd -dlb no"
fi
if [ $timeout -gt 0 ]
then
    cmd="$cmd -maxh $timeout"
fi
#if [ $mpiN -gt 0 ]
#then
#    cmd="$cmd -ntmpi $mpiN"
#fi

cd $root_path
cd $run_path

if [ $do_grompp -eq 1 ]
then
    $gmx_serial grompp -f $name.mdp -c $in_filename -p topol.top -o $name.tpr -maxwarn 2
fi

if [[ "$run_mode" == "serial" ]]
then
    echo "serial run not supported"
    #$gmx_serial mdrun $cmd
elif [[ "$run_mode" == "slurm" ]]
then
    echo $cmd
    $cmd
elif [[ "$run_mode" == "trun" ]]
then
    echo "trun not supported"
    #trun -m 1 -ppn=$mpiN $gmx_mdrun mdrun $cmd
else
    echo "wrong run_mode: $run_mode"
fi

#sbatch -J gromacs -p max1n -N 1 --reservation=test --ntasks-per-node=$mpiN --wrap="$cmd"

cd $root_path
cd $exe_path

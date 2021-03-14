#!/bin/bash
### '###' - human actions/comments
### '#'   - commands to execute for gmx/chimera
### =================================================

set -e
gmx_serial=gmx_mpi
gmx_serial=gmx_ser_gpu
gmx_serial=gmx_ser_newhead

argc=$#
if [ $argc -ne 2 ]
then
    printf "usage:\n$0   job_name   name\n"
    exit 1
fi

job_id=$1
name=$2

root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/K
topol_filename=topol.top
cd $root_path
cd $run_path

$gmx_serial trjconv -s $name.tpr -f $name.gro -pbc nojump -o $name\_nojump.gro < output_whole_sys0.in

#$gmx_serial grompp -f eql.mdp -c eql.gro -p $topol_filename -o eql.tpr
#if [ $gpu_id -eq -1 ]
#then
#        $gmx_mdrun mdrun -v -deffnm eql -ntomp $omp
#else
#        $gmx_mdrun mdrun -v -deffnm eql -ntomp $omp -gpu_id $gpu_id
#fi

cd $root_path
cd $exe_path

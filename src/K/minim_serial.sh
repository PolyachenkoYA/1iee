#!/bin/bash
### '###' - human actions/comments
### '#'   - commands to execute for gmx/chimera
### =================================================

set -e
gmx_serial=gmx_mpi
gmx_serial=gmx_ser_gpu
gmx_serial=gmx_ser_newhead

argc=$#
if [ $argc -ne 5 ]
then
    printf "usage:\n$0   job_name   ompN   mpiN   gpu_id   name\n"
    exit 1
fi

job_id=$1
omp=$2
mpiN=$3
gpu_id=$4
name=$5

root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/K
topol_filename=topol.top
cd $root_path
cd $run_path

$gmx_serial grompp -f minim.mdp -c 1iee_wions.gro -p $topol_filename -o $name.tpr

if [ $gpu_id -eq -1 ]
then
    $gmx_serial mdrun -v -deffnm $name -ntomp $omp
else
    $gmx_serial mdrun -v -deffnm $name -ntomp $omp -gpu_id $gpu_id
fi

$gmx_serial trjconv -s $name.tpr -f $name.gro -pbc nojump -o $name\_nojump.gro < output_whole_sys0.in
cp $name.gro eql.gro

#$gmx_serial grompp -f eql.mdp -c eql.gro -p $topol_filename -o eql.tpr
#if [ $gpu_id -eq -1 ]
#then
#        $gmx_mdrun mdrun -v -deffnm eql -ntomp $omp
#else
#        $gmx_mdrun mdrun -v -deffnm eql -ntomp $omp -gpu_id $gpu_id
#fi

cp eql.gro npt.gro
cp eql.gro nvt.gro

cd $root_path
cd $exe_path
#./mainrun.sh $job_id $omp

#!/bin/bash

set -e
gmx_serial=gmx_ser_gpu
#gmx_serial=gmx_serial

gmx_mdrun=gmx_angara
gmx_mdrun="$gmx_serial"

argc=$#
if [ $argc -ne 2 ] && [ $argc -ne 3 ] && [ $argc -ne 4 ] && [ $argc -ne 5 ]
then
        printf "usage:\n$0   job_name   ompN   [gpu_id (-1=no_gpu)   [continue(0)   [mpiN(1)]]]\n"
        exit 1
fi
if [ $argc -gt 2 ]
then
        gpu_id=$3
else
        gpu_id=-1
fi
if [ $argc -gt 3 ]
then
        cont=$4
else
        cont=0
fi
if [ $argc -gt 4 ]
then
        mpiN=$5
else
        mpiN=1
fi
job_id=$1
root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/maxsol
ompN=$2
cd $root_path
cd $run_path

# ====================================

if [ $cont -eq 0 ]
then
        $gmx_serial grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
fi
if [ $gpu_id -eq -1 ]
then
        $gmx_mdrun mdrun -v -deffnm nvt -ntomp $ompN -ntmpi $mpiN -nsteps 1000000 -cpi nvt.cpt -maxh 24 -cpo nvt.cpt
else
        $gmx_mdrun mdrun -v -deffnm nvt -ntomp $ompN -ntmpi $mpiN -gpu_id $gpu_id -nsteps 1000000 -cpi nvt.cpt -maxh 24 -cpo nvt.cpt
fi
$gmx_serial trjconv -s em.tpr -f nvt.gro -pbc nojump -o nvt_nojump.gro < output_whole_sys0.in

#rm nvt.trr

#$gmx_serial grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
#$gmx_mdrun mdrun -v -deffnm npt -ntomp $omp -gpu_id $gpu_id
#$gmx_serial trjconv -s em.tpr -f npt.gro -pbc nojump -o npt_nojump.gro < output_whole_sys0.in

#$gmx_serial grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_mainrun.tpr
#$gmx_mdrun mdrun -v -deffnm md_mainrun -ntomp $omp -gpu_id $gpu_id
#$gmx_serial trjconv -s em.tpr -f md_mainrun.gro -pbc nojump -o md_mainrun_nojump.gro < output_whole_sys0.in

cd $root_path
cd $exe_path


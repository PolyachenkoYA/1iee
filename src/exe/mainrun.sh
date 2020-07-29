#!/bin/bash

set -e
gmx_executable=gmx

argc=$#
if [ $argc -ne 2 ]; then
        printf "usage:\n$0   job_name   ompN\n"
        exit 1
fi
job_id=$1
root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/exe
omp=$2
cd $root_path
cd $run_path

# ====================================

$gmx_executable grompp -f minim.mdp -c 1iee_wions.gro -p topol.top -o em.tpr
$gmx_executable mdrun -v -deffnm em -ntomp $omp
$gmx_executable trjconv -s em.tpr -f em.gro -pbc nojump -o em_nojump.gro < output_whole_sys0.in

exit 0

$gmx_executable grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
$gmx_executable mdrun -v -deffnm nvt -ntomp $omp
$gmx_executable trjconv -s em.tpr -f nvt.gro -pbc nojump -o nvt_nojump.gro < output_whole_sys0.in

rm nvt.trr

#$gmx_executable grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
#$gmx_executable mdrun -v -deffnm npt -ntomp $omp
#$gmx_executable trjconv -s em.tpr -f npt.gro -pbc nojump -o npt_nojump.gro < output_whole_sys0.in

#$gmx_executable grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_mainrun.tpr
#$gmx_executable mdrun -v -deffnm md_mainrun -ntomp $omp
#$gmx_executable trjconv -s em.tpr -f md_mainrun.gro -pbc nojump -o md_mainrun_nojump.gro < output_whole_sys0.in

cd $root_path
cd $exe_path

#!/bin/bash

set -e

omp=$1

cd ../../run

# ====================================

#gmx_mpi grompp -f minim.mdp -c 1iee_solv.gro -p topol.top -o em.tpr
#gmx_mpi mdrun -v -deffnm em -ntomp $omp -gpu_id 0 -ntmpi 1
#gmx_mpi trjconv -s em.tpr -f em.gro -pbc nojump -o em_nojump.gro < output_whole_sys0.in

gmx_mpi grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx_mpi mdrun -v -deffnm nvt -ntomp $omp -gpu_id 0 -ntmpi 1
gmx_mpi trjconv -s nvt.tpr -f nvt.gro -pbc nojump -o nvt_nojump.gro < output_whole_sys0.in

#gmx_mpi grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
#gmx_mpi mdrun -v -deffnm npt -ntomp $omp -gpu_id 0 -ntmpi 1
#gmx_mpi trjconv -s npt.tpr -f npt.gro -pbc nojump -o npt_nojump.gro < output_whole_sys0.in

#gmx_mpi grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_mainrun.tpr
#gmx_mpi mdrun -v -deffnm md_mainrun -ntomp $omp -gpu_id 0 -ntmpi 1
#gmx_mpi trjconv -s md_mainrun.tpr -f md_mainrun.gro -pbc nojump -o md_mainrun_nojump.gro < output_whole_sys0.in


cd ../src/exe

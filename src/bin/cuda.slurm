#!/bin/bash

#SBATCH -p max1n
#SBATCH -J gromacs
#SBATCH --gres=gpu:1

mpi=$1
omp=$2

gmx grompp -f minim.mdp -c 1iee_solv.gro -p topol.top -o em.tpr
srun --ntasks-per-node=$mpi gmx_mpi mdrun -v -deffnm em -ntomp $omp

gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
srun --ntasks-per-node=$mpi gmx_mpi mdrun -v -deffnm nvt -ntomp $omp

gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
srun --ntasks-per-node=$mpi gmx_mpi mdrun -v -deffnm npt -ntomp $omp

gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_mainrun.tpr
srun --ntasks-per-node=$mpi gmx_mpi mdrun -v -deffnm md_mainrun -ntomp $omp

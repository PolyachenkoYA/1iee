#!/bin/bash

set -e

omp=$1

cd ../run

$HOME/gromacs-2018.8/build_serial/bin/gmx grompp -f minim.mdp -c 1iee_solv.gro -p topol.top -o em.tpr
$HOME/gromacs-2018.8/build_serial/bin/gmx mdrun -v -deffnm em -ntomp $omp -gpu_id 0 -ntmpi 1

$HOME/gromacs-2018.8/build_serial/bin/gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
$HOME/gromacs-2018.8/build_serial/bin/gmx mdrun -v -deffnm nvt -ntomp $omp -gpu_id 0 -ntmpi 1

$HOME/gromacs-2018.8/build_serial/bin/gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
$HOME/gromacs-2018.8/build_serial/bin/gmx mdrun -v -deffnm npt -ntomp $omp -gpu_id 0 -ntmpi 1

$HOME/gromacs-2018.8/build_serial/bin/gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_mainrun.tpr
$HOME/gromacs-2018.8/build_serial/bin/gmx mdrun -v -deffnm md_mainrun -ntomp $omp -gpu_id 0 -ntmpi 1

cd ../init

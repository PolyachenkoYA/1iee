#!/bin/bash
### '###' - human actions/comments
### '#'   - commands to execute for gmx/chimera
### =================================================

set -e
gmx_serial=gmx_ser_newhead
#gmx_serial=gmx_serial
#gmx_serial=gmx_ser_gpu

gmx_mdrun=gmx_mpi
#gmx_mdrun=gmx_angara
#gmx_mdrun="$gmx_serial"

argc=$#
if [ $argc -ne 2 ]
then
        printf "usage:\n$0   job_name   out_name\n"
        exit 1
fi
job_id=$1
name=$2
root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/water_cube

cd $root_path
cd $run_path

$gmx_serial solvate -cs amber03w.ff/tip4p2005.gro -o $name.gro -box 7.70610 7.70610 7.44460 -p topol.top
#$gmx_serial solvate -cs tip4p -o solvated.gro -box 3 3 3 -p topol.top

cd $root_path
cd $exe_path

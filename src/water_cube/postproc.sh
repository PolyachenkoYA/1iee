#!/bin/bash

set -e
gmx_serial=gmx_mpi

gmx_mdrun=gmx_mpi
gmx_mdrun=gmx_angara
gmx_mdrun=$gmx_serial

argc=$#
if [ $argc -ne 1 ]
then
        printf "usage:\n$0   job_name\n"
        exit 1
fi
job_id=$1
root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/water_cube

./clear_restore.sh $job_id

cd $root_path
cd $run_path

$gmx_serial energy -f eql.edr -o eql-temp.xvg < temp.in
$gmx_serial energy -f min.edr -o min-energy.xvg < en1.in
$gmx_serial energy -f min2.edr -o min2-energy.xvg < en2.in

cd $root_path
cd $exe_path

#!/bin/bash

set -e
gmx_serial=gmx_mpi
gmx_serial=gmx_ser_gpu

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
nvt_summary_filename=nvt_summary.txt
nvt_bigbox_summary_filename=nvt_bigbox_summary.txt

cd $root_path
cd $run_path

$gmx_serial energy -f nvt.edr -o nvt.xvg < nvt.in > $nvt_summary_filename
$gmx_serial energy -f nvt_bigbox.edr -o nvt_bigbox.xvg < nvt.in > $nvt_bigbox_summary_filename

python $root_path/$exe_path/extract_nvt.py $nvt_summary_filename nvt_pressure.dat
python $root_path/$exe_path/extract_nvt.py $nvt_bigbox_summary_filename nvt_bigbox_pressure.dat

cd $root_path
cd $exe_path

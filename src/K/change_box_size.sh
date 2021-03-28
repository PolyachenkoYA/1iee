#!/bin/bash -f

set -e

gmx_serial=gmx_mpi
gmx_serial=$HOME/gromacs-2020/build_lin/bin/gmx
gmx_serial=gmx_ser_20

argc=$#
if [ $argc -ne 4 ]
then
        printf "usage:\n$0   job_id   old_name   new_name   multiplier\n"
        exit 1
fi

job_id=$1
old_name=$2
new_name=$3
mult=$4

#root_path=$(git rev-parse --show-toplevel)
root_path=`cat git_root_path`
run_path=run/$job_id
exe_path=src/water_cube

cd $root_path
cd $run_path

new_box_sizes=$(gawk -v m=$mult '
END{
	print $1*m, $2*m, $3*m
}
' < $old_name)

$gmx_serial editconf -f $old_name -o $new_name -box $new_box_sizes

#pwd
cd $root_path
#pwd
cd $exe_path

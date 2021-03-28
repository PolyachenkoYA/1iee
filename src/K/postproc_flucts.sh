#!/bin/bash

set -e
gmx_serial=gmx_mpi
gmx_serial=$HOME/gromacs-2020/build_lin/bin/gmx
gmx_serial=gmx_ser_20

#gmx_mdrun=gmx_mpi
#gmx_mdrun=$HOME/gromacs-2020/build/bin/gmx_mpi
#gmx_mdrun=gmx_angara
#gmx_mdrun=$gmx_serial

npt_summary_filename=npt_summary_measures.txt

argc=$#
if [ $argc -ne 1 ]
then
        printf "usage:\n$0   job_name\n"
        exit 1
fi
job_id=$1
#root_path=$(git rev-parse --show-toplevel)
root_path=`cat git_root_path`
run_path=run/$job_id
exe_path=src/K

cd $root_path
cd $run_path

$gmx_serial energy -f npt.edr -o npt.xvg -fluct_props -driftcorr < npt.in > $npt_summary_filename
$root_path/$exe_path/extract_K.sh $npt_summary_filename > ./K.dat
$gmx_serial trjconv -f npt.xtc -o npt_nojump.xtc -s em_nojump.gro -pbc nojump <  output_whole_sys0.in

cd $root_path
cd $exe_path

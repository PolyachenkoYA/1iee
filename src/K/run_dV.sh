job_id=$1
ompN=$2

mode=nvt
mpiN=1
gpu_id=-1

initial_base=$mode
new_base=$mode\_bigbox
initial_filename=$initial_base.gro
new_filename=$new_base.gro

./clear_restore.sh $job_id
#cp ../../run/save/* ../../run/$job_id/
./preproc.sh $job_id $ompN $mpiN $gpu_id 1 1 2 2118 1iee112_prot4gmx.pdb
./mainrun_serial.sh $job_id $ompN $mpiN $gpu_id $initial_base
./change_box_size.sh $job_id $initial_filename $new_filename 1.001
./mainrun_serial.sh $job_id $ompN $mpiN $gpu_id $new_base
./postproc_dV.sh $job_id

# sbatch --nodelist=host27 -J gromacs -p max1n -N 1 --reservation=test --ntasks-per-node=1 --gres=gpu:1 --wrap="./run_dV.sh tst 6 1 nvt"

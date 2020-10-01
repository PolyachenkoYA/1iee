job_id=$1
ompN=$2
mpiN=$3
mode=$4

initial_base=$mode
new_base=$mode\_bigbox
initial_filename=$initial_base.gro
new_filename=$new_base.gro

./clear_restore.sh $job_id
./preproc.sh $job_id $ompN $mpiN
./mainrun_serial2.sh $job_id $ompN $mpiN 0 eql2 $initial_base
./change_box_size.sh $job_id $initial_filename $new_filename 1.01
./mainrun_serial2.sh $job_id $ompN $mpiN 0 $initial_base $new_base
./postproc_dV.sh $job_id

# sbatch --nodelist=host27 -J gromacs -p max1n -N 1 --reservation=test --ntasks-per-node=1 --gres=gpu:1 --wrap="./run_dV.sh tst 6 1 nvt"

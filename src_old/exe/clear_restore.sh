argc=$#
if [ $argc -ne 1 ]; then
	printf "usage:\n$0   job_name\n"
        exit 1
fi
job_id=$1
root_path=$(git rev-parse --show-toplevel)
exe_path=src/exe
src_path=src
run_path=run/$job_id

cd $root_path

rm -rf $run_path
mkdir -p $run_path
cp $src_path/*.mdp $run_path
cp $src_path/*.pdb $run_path
cp $src_path/*.in $run_path

cd $exe_path

echo "$run_path cleared"

set -e
gmx_exeutable=gmx_mpi

argc=$#
if [ $argc -ne 2 ] && [ $argc -ne 3 ]
then
        printf "usage:\n$0   job_name\n"
        exit 1
fi
job_id=$1
root_path=$(git rev-parse --show-toplevel)
run_path=run/$job_id
exe_path=src/maxsol
cd $root_path
cd $run_path

$gmx_exeutable trjconv -s md_mainrun.tpr -f md_mainrun.xtc -o md_main_noPBC.xtc -pbc mol -center
# we want to center Protein (1)
# and return the whole system (0)

$gmx_exeutable rms -s md_mainrun.tpr -f md_main_noPBC.xtc -o rmsd.xvg -tu ns
# rmsd of the structure is interesting, so fit to the Backbone (4)
# and analyse Backbone (4)

$gmx_exeutable rms -s em.tpr -f md_main_noPBC.xtc -o rmsd_xtal.xvg -tu ns
# this is rmsd relative to the crystal state. Choices are the same - fit Backbone (4)
# and analyse Backbone (4)

$gmx_exeutable gyrate -s md_mainrun.tpr -f md_main_noPBC.xtc -o gyrate.xvg
# giration radius of a Protein (1) represent the stability

cd $root_path
cd $exe_path

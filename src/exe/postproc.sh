set -e
cd ../run

$HOME/gromacs-2018.8/build_serial/bin/gmx trjconv -s md_mainrun.tpr -f md_mainrun.xtc -o md_main_noPBC.xtc -pbc mol -center
# we want to center Protein (1)
# and return the whole system (0)

$HOME/gromacs-2018.8/build_serial/bin/gmx rms -s md_mainrun.tpr -f md_main_noPBC.xtc -o rmsd.xvg -tu ns
# rmsd of the structure is interesting, so fit to the Backbone (4)
# and analyse Backbone (4)

$HOME/gromacs-2018.8/build_serial/bin/gmx rms -s em.tpr -f md_main_noPBC.xtc -o rmsd_xtal.xvg -tu ns
# this is rmsd relative to the crystal state. Choices are the same - fit Backbone (4)
# and analyse Backbone (4)

$HOME/gromacs-2018.8/build_serial/bin/gmx gyrate -s md_mainrun.tpr -f md_main_noPBC.xtc -o gyrate.xvg
# giration radius of a Protein (1) represent the stability

cd ../init

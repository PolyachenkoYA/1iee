### '###' - human actions/comments
### '#'   - commands to execute for gmx/chimera
### =================================================

set -e
cd ../../run

### Download 1iee.pdb from database

# chimera --start "Unit Cell" --start "Command Line" --start "Reply Log" pdbs/1iee.pdb

### Tools -> HigherOrderStructure -> UnitCell -> outline

### Tools -> HigherOrderStructure -> UnitCell -> make copies

# findclash # test self colorClashes true setAttrs true log true namingStyle command overlapCutoff 2 selectClashes true saveFile clashesALL.txt
### Options-> Number of sells 3 3 3, offset -1 -1 -1
### there will be a number of lines in Reply Log. They have 'dist' field. Some of them have dist >~0.2. Some have <0.01. Use those with <0.01 in the next step.
# python exe/get_chashes.py clashesALL.txt 8
### get all water-residues (molecules) to delete duw to overlap

### these are residues numbers (HOH molecule's numbers in case of solvent) given by the findclash
# delete #0-7:1092,1103,1139,1169

# combine # newchainids true name 1iee_cryst close true

### write the whole system (all 0-7 groups were encorporated into a single group 0 by 'combine' command above) relative to the outline (group 8)
# write relative 0 8 pdb/1iee_cryst.pdb


### go to website and protonate (prepare protein) -> 1iee_prot.pdb
### this potentially might by done using pythong code given by the website.


../src/exe/playmol2gmx.sh 1iee_prot.pdb 1iee_prot4gmx.pdb

gmx pdb2gmx -f 1iee_prot4gmx.pdb -o 1iee_init.gro -water spce < protonation_gromacs.in
gmx editconf -f 1iee_init.gro -o 1iee_newbox.gro -bt triclinic -box 7.7061 7.7061 3.7223 -noc
gmx grompp -f ions.mdp -c 1iee_newbox.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o 1iee_wions.gro -p topol.top -pname NA -nname CL -neutral -rmin 0.28  < genion_gromacs.in
#gmx editconf -f 1iee_wions.gro -o 1iee_wions_cent.gro -c
gmx solvate -cp 1iee_wions.gro -cs spc216.gro -o 1iee_solv.gro -p topol.top

#################################################

gmx grompp -f minim.mdp -c 1iee_solv.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em -ntomp 6


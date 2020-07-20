### '###' - human actions/comments
### '#'   - commands to execute for gmx/chimera
### =================================================

set -e
cd ../run

### Download 1iee.pdb from database

# chimera --start "Unit Cell" --start "Command Line" --start "Reply Log" pdbs/1iee.pdb

### Tools -> HigherOrderStructure -> UnitCell -> outline

### Tools -> HigherOrderStructure -> UnitCell -> make copies

# combine # newchainids true name 1iee_cryst close true

### solvent only
# findclash :SOL test self colorClashes true setAttrs true log true namingStyle command overlapCutoff 2.5 selectClashes true saveFile clashesSOL.txt
###########  OR  ###########
### all atoms
# findclash # test self colorClashes true setAttrs true log true namingStyle command overlapCutoff 2.5 selectClashes true saveFile clashesALL.txt
### Options-> Number of sells 3 3 3, offset -1 -1 -1
### there will be a number of lines in Reply Log. They have 'dist' field. Some of them have dist >~0.2. Some have <0.01. Use those with <0.01 in the next step.
# ./get_clashes.sh clashesALL.txt
### get all water-residues (molecules) to delete duw to overlap

### these are residues numbers (HOH molecule's numbers in case of solvent) given by the findclash
# delete #7:1103,#6:1103,#3:1103,#1:1103,#0:1103,#4:1103,#5:1103,#2:1103,#6:1139,#5:1139,#4:1139,#3:1139,#2:1139,#1:1139,#7:1139,#0:1139,#6:1169,#1:1169,#4:1169,#2:1169,#4:1092,#2:1092,#1:1092,#7:1169,#7:1092,#6:1092,#0:1092,#0:1169,#5:1169,#3:1169,#5:1092,#3:1092

### write the whole system (all 0-7 groups were encorporated into a single group 0 by 'combine' command above) relative to the outline (group 8)
# write relative 0 8 pdb/1iee_cryst.pdb


./clear_restore.sh

./playmol2gmx.sh 1iee_cryst.pdb 1iee_cryst_prot.pdb

#grep -v HOH 1iee.pdb | grep -v -w CL | grep -v -w NA > 1iee_clean.pdb

#gmx pdb2gmx -f 1iee_cryst_prot.pdb -o 1iee_init.gro -water tip4p -asp -his -glu
###                                         asp(52) -1 ->  0
###                                         his(15)  0 -> +1
### possibly (pKa = 6.89 vs 7.0 in gromacs) glu(35) -1 ->  0
gmx pdb2gmx -f 1iee_cryst_prot.pdb -o 1iee_init.gro -water spce < protonation_gromacs.in

##############gmx editconf -f 1iee_init.gro -o 1iee_newbox.gro -bt triclinic -box 7.7061 7.7061 3.7223 -noc
##############gmx grompp -f ions.mdp -c 1iee_newbox.gro -p topol.top -o ions.tpr
##############gmx genion -s ions.tpr -o 1iee_wions.gro -p topol.top -pname NA -nname CL -neutral -rmin 0.4
##############gmx trjconv -f 1iee_wions.gro -s 1iee_wions.gro -o 1iee_wions.pdb

### add header to the pdb

#gmx pdb2gmx -f 1iee_wions_cr.pdb -o 1iee_wions_cr.gro -water spce
#gmx editconf -f 1iee_wions_cr.gro -o 1iee_wions_cr_cent.gro -c
#gmx solvate -cp 1iee_wions_cr_cent.gro -cs spc216.gro -o 1iee_solv.gro -p topol.top

#################################################

#gmx grompp -f minim.mdp -c 1iee_solv.gro -p topol.top -o em.tpr
#gmx mdrun -v -deffnm em -ntomp 6



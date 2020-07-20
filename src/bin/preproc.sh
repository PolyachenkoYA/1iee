### '###' - human actions/comments
### '#'   - commands to execute for gmx/chimera
### =================================================

set -e
cd ../../run

### Download 1iee.pdb from database

# chimera --start "Unit Cell" --start "Command Line" --start "Reply Log" 1iee.pdb

### Tools -> HigherOrderStructure -> UnitCell -> outline

### Tools -> HigherOrderStructure -> UnitCell -> make copies

# combine # newchainids false name 1iee_wions close true

### solvent only
# findclash :SOL test self colorClashes true setAttrs true log true namingStyle simple overlapCutoff 1.95 selectClashes true saveFile clashesSOL.txt
###########  OR  ###########
### all atoms
# findclash # test self colorClashes true setAttrs true log true namingStyle simple overlapCutoff 1.95 selectClashes true saveFile clashesALL.txt
### there will be a number of lines in Reply Log. They have 'dist' field. Some of them have dist >~0.2. Some have <0.01. Use those with <0.01 in the next step.

### these are residues numbers (HOH molecule's numbers in case of solvent) given by the findclash
# delete :1173:1248:1611:1656:1686:2049:2094:2124:2487:2562

### write the whole system (all 0-7 groups were encorporated into a single group 0 by 'combine' command above) relative to the outline (group 8)
# write relative 0 8 1iee_wions_cr.pdb



#grep -v HOH 1iee.pdb | grep -v -w CL | grep -v -w NA > 1iee_clean.pdb
gmx pdb2gmx -f 1iee.pdb -o 1iee_init.gro -water tip4p
###                                         asp(52) -1 ->  0
###                                         his(15)  0 -> +1
### possibly (pKa = 6.89 vs 7.0 in gromacs) glu(35) -1 ->  0
gmx editconf -f 1iee_init.gro -o 1iee_newbox.gro -bt triclinic -box 7.7061 7.7061 3.7223 -noc
gmx grompp -f ions.mdp -c 1iee_newbox.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o 1iee_wions.gro -p topol.top -pname NA -nname CL -neutral -rmin 0.4
gmx trjconv -f 1iee_wions.gro -s 1iee_wions.gro -o 1iee_wions.pdb

### add header to the pdb

#gmx pdb2gmx -f 1iee_wions_cr.pdb -o 1iee_wions_cr.gro -water spce
#gmx editconf -f 1iee_wions_cr.gro -o 1iee_wions_cr_cent.gro -c
#gmx solvate -cp 1iee_wions_cr_cent.gro -cs spc216.gro -o 1iee_solv.gro -p topol.top

#################################################

#gmx grompp -f minim.mdp -c 1iee_solv.gro -p topol.top -o em.tpr
#gmx mdrun -v -deffnm em -ntomp 6



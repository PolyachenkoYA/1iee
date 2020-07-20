set -e
cd ../run

gmx pdb2gmx -f 1iee_wions_cr.pdb -o 1iee_wions_cr.gro -water spce
gmx editconf -f 1iee_wions_cr.gro -o 1iee_wions_cr_cent.gro -c
gmx solvate -cp 1iee_wions_cr_cent.gro -cs spc216.gro -o 1iee_solv.gro -p topol.top

cd ../init

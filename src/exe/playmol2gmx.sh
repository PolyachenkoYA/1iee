#!/bin/bash -f

tmpfile=$(mktemp /tmp/abc-script.XXXXXX)

awk '
function sub_pair(letter, residue){
    sub(sprintf("H%s2 %s", letter, residue), sprintf("H%s1 %s", letter, residue));
    sub(sprintf("H%s3 %s", letter, residue), sprintf("H%s2 %s", letter, residue));
}
BEGIN{}
{
# LYS
    sub_pair("B", "LYS");
    sub_pair("G", "LYS");
    sub_pair("D", "LYS");
    sub_pair("E", "LYS");

    sub_pair("B", "PHE");

    sub_pair("A", "GLY");

    sub_pair("B", "ARG");
    sub_pair("G", "ARG");
    sub_pair("D", "ARG");

    sub_pair("B", "CYX");

    sub_pair("B", "GLU");
    sub_pair("G", "GLU");

    sub_pair("B", "LEU");

    sub_pair("B", "MET");
    sub_pair("G", "MET");

    sub_pair("B", "HIP");

    sub_pair("B", "ASP");

    sub_pair("B", "ASN");

    sub_pair("B", "TYR");

    sub_pair("B", "SER");

    sub_pair("B", "TRP");

    sub_pair("B", "GLH");
    sub_pair("G", "GLH");

    sub_pair("B", "GLN");
    sub_pair("G", "GLN");

    sub_pair("B", "ASH");

    sub_pair("G1", "ILE");

    sub_pair("B", "PRO");
    sub_pair("G", "PRO");
    sub_pair("D", "PRO");

    print
}
END{}' < $1 > $tmpfile

rm $1
mv $tmpfile $1

gmx pdb2gmx -f output.pdb -o 1iee_init.gro -water spce < input_protonation

dir_id=$1

cd ..

rm "../run$1/*"
cp *.mdp ../run$1
cp *.pdb ../run$1
cp *.in ../run$1

cd exe

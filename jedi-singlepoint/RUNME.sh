#!/bin/bash
# Clean up
rm -rf output
mkdir output
cd output
# *** Copy input
#cp ../input/system.pdb .
cp ../input/system_allh.pdb .
cp ../input/ligand.pdb . 
cp ../input/*.mdp .
cp ../input/jedi.params .
cp ../input/plumed.dat .
# *** Generate paraneters for ligand or other cofactors if needed
# THIS PART IS MISSING 
# *** Generate input files for JEDI ***
../jedi-setup.py -i system_allh.pdb -l ligand.pdb -c 0.4 -s 0.15 -a apolar.pdb -p polar.pdb -g grid.pdb
# *** Build system ***
# may need optional keywords to line below if ligands are present 
pdb2gmx -f system_allh.pdb -o system_processed.gro -p system.top -ff amber99sb-ildn -water tip3p -ignh 
# *** Single point equilibration ***
grompp -f sp.mdp -c system_processed.gro -p system.top -o sp.tpr
mdrun -deffnm sp -plumed plumed.dat -v

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
../jedi-setup.py -i system_allh.pdb -l ligand.pdb -c 0.6 -s 0.15 -a apolar.pdb -p polar.pdb -g grid.pdb
#python ../../jedi-setup.py -i system_allh.pdb -s 0.15 -a apolar.pdb -p polar.pdb -g grid.pdb -l ligand.pdb -f no 
# *** Build system ***
# may need optional keywords to line below if ligands are present 
gmx pdb2gmx -f system_allh.pdb -o system_processed.gro -p system.top -ff amber99sb-ildn -water tip3p -ignh 
# *** Single point equilibration ***
gmx grompp -f useless.mdp -c system_processed.gro -p system.top -o sp.tpr
gmx mdrun -deffnm sp -plumed plumed.dat -v -rerun system_processed.gro -ntomp $1
#gmx mdrun -ntomp $1 -deffnm sp -plumed plumed.dat -v -pin on -nsteps $2 
#gmx mdrun -deffnm sp -plumed plumed.dat -v -rerun system_processed.gro

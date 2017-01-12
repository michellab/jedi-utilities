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
gmx pdb2gmx -f system_allh.pdb -o system_processed.gro -p system.top -ff amber99sb-ildn -water tip3p -ignh 
gmx editconf -f system_processed.gro -o system_newbox.gro -c -d 1.0 -bt cubic
gmx solvate -cp system_newbox.gro -cs spc216.gro -o system_solv.gro -p system.top
gmx grompp -f min.mdp -c system_solv.gro -p system.top -o ions.tpr
gmx genion -s ions.tpr -o system_solv_ions.gro -p system.top -neutral #CHECK WHY DOES NOT WORK -conc 0.15
# *** Minimization ***
gmx grompp -f min.mdp -c system_solv_ions.gro -p system.top -o min.tpr
gmx mdrun -deffnm min 
# *** NVT equilibration ***
gmx grompp -f nvt.mdp -c min.gro -p system.top -o nvt.tpr
gmx mdrun -deffnm nvt -plumed plumed.dat -v
# *** NPT equilibration ***
gmx grompp -f npt.mdp -c nvt.gro -p system.top -o npt.tpr
gmx mdrun -deffnm npt -plumed plumed.dat -v
# *** Production NPT ***
gmx grompp -f md.mdp -c npt.gro -p system.top -o md.tpr -maxwarn 1
gmx mdrun -deffnm md -plumed plumed.dat -v


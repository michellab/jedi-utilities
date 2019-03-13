#!/bin/bash
# Gromacs executable
gmx=gmxserial
# Clean up
rm -rf output-SITE
mkdir output-SITE
cd output-SITE
# *** Copy input
#cp ../input/system.pdb .
cp ../input/protein.pdb .
cp ../input/cofactor.pdb .
cp ../input/waterbox.pdb .
cp ../input/reference_ligand.pdb .
cp ../input/*.mdp .
cp ../input/jedi.params .
cp ../input/plumed-SITE.dat .
cp ../input/probe_2a_nobond.pdb .
# *** Copy paraneters for ligand or other cofactors if needed
cp ../parameters/cofactor_GMX.* .
# *** Build system ***
cat waterbox.pdb >> protein.pdb
$gmx pdb2gmx -f protein.pdb -o protein.gro -p protein.top -ff amber99sb-ildn -water tip3p -ignh
# The box line has to be deleted and reinserted
tail -n 1 protein.gro > box.txt
sed -i '$ d' protein.gro
# Merge Coordinates of the cofactor
rm -r cofactor_GMX.gro
$gmx editconf -f cofactor.pdb -o cofactor_GMX.gro
python ../merge_coordinates_gro.py protein.gro cofactor_GMX.gro 14800_prep.gro
cat box.txt >> 14800_prep.gro
rm -rf temp.gro
#Merge atom types
python ../merge_atomtypes.py cofactor_GMX.itp
sed  "/forcefield.itp/a ; Include Ligand Atomtypes" protein.top > 14800_prep.top
sed -i  "/; Include Ligand Atomtypes/a #include \"LigandAtomTypes.itp\"" 14800_prep.top
sed -i  "/#include \"LigandAtomTypes.itp\"/a ; Include Ligand Topologies" 14800_prep.top

for topology in $(ls *GMX.itp); do
      sed -i "/; Include Ligand Topologies/a  #include \"$topology\"\n" 14800_prep.top
      molecule=`echo $topology | awk 'BEGIN {FS="_GMX"} {print $1}'`
      echo "SAH    1" >>  14800_prep.top
   done
# *** Generate input files for JEDI ***
$gmx editconf -f 14800_prep.gro -o 14800_prep.pdb
rm -rf temp.pdb
../jedi-setup.py -i 14800_prep.pdb -l reference_ligand.pdb -c 0.6 -s 0.15 -a apolar.pdb -p polar.pdb -g grid.pdb

# *** Single point equilibration ***
$gmx grompp -f useless.mdp -c 14800_prep.gro -p 14800_prep.top -o sp.tpr
$gmx mdrun -deffnm sp -plumed plumed-SITE.dat -v -rerun 14800_prep.gro 
# *** Visualize the grid ***
python ../grid_analysis.py -i grid-step-0.xyz -a acti-step-0.txt -o activity-grid.pdb


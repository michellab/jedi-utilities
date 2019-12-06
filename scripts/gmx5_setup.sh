# 0. Define functions and variables

# functions
function globexist # To check if ligands exist ([-e *GMX*] crashes if there is more than 1)
  {
   test -e "$1" -o -L "$1"
  }
#paths
utilities=~/github/jedi-utilities
mdp=$utilities/mdp
gmx=/home/joan/local/gromacs-5.1.4-plumed2_julien/bin/gmx_mpi
merge_gro=$utilities/scripts/merge_coordinates_gro.py
merge_atomtypes=$utilities/scripts/merge_atomtypes.py

#variables
system=$1

# MD parameters
forcefield=amber99sb-ildn
solvent=tip3p
boxtype=dodecahedron
boxsize=1.0

$gmx pdb2gmx -f $system'_protein.pdb' -o $system'_protein.gro' -p $system'_protein.top' -ff $forcefield -water $solvent -ignh

if globexist  *GMX* ; then
   # The box line  nd the waters have to be deleted and reinserted
   tail -n 1 $system'_protein.gro' > box.txt
   sed -i '$ d' $system'_protein.gro'
   grep 'SOL' $system'_protein.gro' >  water.txt
   sed -i '/SOL/d' $system'_protein.gro'
   grep 'HOH' $system'_protein.gro' >> water.txt
   sed -i '/HOH/d' $system'_protein.gro'
   grep 'WAT' $system'_protein.gro' >> water.txt
   sed -i '/WAT/d' $system'_protein.gro'
   # Merge Coordinates
   if globexist ligand*GMX.gro; then
      cof=$(ls LIG*GMX.gro)
      lig=ligand_GMX.gro
      gmxs=$(echo $cof $lig)
   else
      gmxs=$(ls *GMX.gro)
   fi
   python $merge_gro $system'_protein.gro' $gmxs $system'_prep.gro'
   cat water.txt  >> $system'_prep.gro'
   cat box.txt    >> $system'_prep.gro'

   grep 'SOL' $system'_protein.top' >  water.top
   sed -i '/SOL/d' $system'_protein.top'

   if globexist ligand*GMX.itp; then
      cof=$(ls LIG*GMX.itp)
      lig=ligand_GMX.itp
      itp=$(echo $cof $lig)
   else
      itp=$(ls *GMX.itp)
   fi
   python $merge_atomtypes $itp
   sed  "/forcefield.itp/a ; Include Ligand Atomtypes" $system'_protein.top' > $system'.top'
   sed -i  "/; Include Ligand Atomtypes/a #include \"LigandAtomTypes.itp\"" $system'.top'
   sed -i  "/#include \"LigandAtomTypes.itp\"/a ; Include Ligand Topologies" $system'.top'
   for topology in $itp; do
      sed -i "/; Include Ligand Topologies/a  #include \"$topology\"\n" $system'.top'
      molecule=`echo $topology | awk 'BEGIN {FS="_GMX"} {print $1}'`
      #sed -i "/SOL/i $molecule    1" $system'.top' # Doesn't work when there is no water in the structure
      echo "trying to sed"
      echo "$molecule    1" >>  $system'.top'
   done
   cat water.top  >> $system'.top'
else
   mv $system'_protein.top' $system'.top'
   mv $system'_protein.gro' $system'_prep.gro'
fi

#generate box
$gmx editconf -f $system'_prep.gro' -o $system'_box.gro' -bt $boxtype -d $boxsize
# Solvate the system
$gmx solvate -cp $system'_box.gro' -cs spc216.gro -p $system'.top' -o $system'_water.gro'

# Neutralize
$gmx grompp -f $mdp'/useless.mdp' -c $system'_water.gro' -p $system'.top' -o $system'_ions.tpr'
echo SOL | $gmx genion -s $system'_ions.tpr' -o $system'_ions.gro' -p $system'.top' -pname NA -nname CL -neutral



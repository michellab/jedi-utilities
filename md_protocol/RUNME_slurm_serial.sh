function globexist # To check if ligands exist ([-e *GMX*] crashes if there is more than 1)
  {
   test -e "$1" -o -L "$1"
  }

gmx=/home/joan/local/gromacs-5.1.4.serial/bin/gmx_mpi
plumed=/home/joan/local/bin/plumed
ntomp=16

md="
    function globexist # To check if ligands exist ([-e *GMX*] crashes if there is more than 1)
  {
   test -e "$1" -o -L "$1"
  }
  
    mkdir jedi_files
    mkdir complex apoWet apoDry
    cd jedi_files
      cp ../in/protein_apo.pdb ../in/ligand.pdb ../../../jedi_files/* . 
      python ../../../../../scripts/jedi-setup.py -i protein_apo.pdb -l ligand.pdb -c 0.6 -s 0.15 -a apolar.pdb -p polar.pdb -g grid.pdb
      cp jedi.params jedi_monitor.dat ../complex
      cp jedi.params jedi_monitor.dat ../apoWet
      cp jedi.params jedi_monitor.dat ../apoDry
      $plumed driver --mf_pdb protein_apo.pdb --plumed jedi_monitor.dat
    cd ..

    cd in/
     python ../../../../../scripts/water_replace.py ligand.pdb # Generate the waters that replace the ligand
     cat protein_sol.pdb HOH_ligand.pdb > apoWet_protein.pdb
     cp protein_sol.pdb apoDry_protein.pdb
     cp protein_sol.pdb complex_protein.pdb
     cp complex_protein.pdb *acpype/*GMX.gro *acpype/*GMX.itp ../complex/
     cp apoWet_protein.pdb ../apoWet/
     cp apoDry_protein.pdb ../apoDry/
     if  globexist LIG*acpype; then
      cp  LIG*acpype/*GMX.gro LIG*acpype/*GMX.itp ../apoWet/
      cp  LIG*acpype/*GMX.gro LIG*acpype/*GMX.itp ../apoDry/
     fi
    cd ..

    cd complex
     bash ../../../../../scripts/gmx5_setup.sh complex
     if [ ! -f complex_ions.gro ]; then continue; fi
     $gmx editconf -f complex_ions.gro -o complex_ions.pdb
     grep 'LIG' complex_ions.pdb > ligand_centered.pdb
     python ../../../../../scripts/jedi-setup.py -i complex_ions.pdb -l ligand_centered.pdb -ln LIG -c 0.6 -s 0.15 # the atom indices change from pdb to gro so we have to do it again
     cp apolar.pdb polar.pdb grid.pdb ../apoWet
     cp apolar.pdb polar.pdb grid.pdb ../apoDry
     $gmx grompp -f ../../../../../mdp/em_10000_sd.mdp -c complex_ions.gro -p complex.top -o complex_em_sd.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm complex_em_sd -plumed jedi_monitor.dat -nsteps 10000
     if [ ! -f complex_em_sd.gro ]; then continue; fi
     mv COLVAR sd_COLVAR
     $gmx grompp -f ../../../../../mdp/em_10000_cg.mdp -c complex_em_sd.gro -p complex.top -o complex_em_cg.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm complex_em_cg -plumed jedi_monitor.dat -nsteps 10000
     if [ ! -f complex_em_cg.gro ]; then continue; fi
     mv COLVAR cg_COLVAR
     $gmx grompp -f ../../../../../mdp/nvt.mdp -c complex_em_cg.gro -p complex.top -o complex_nvt.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm complex_nvt -nsteps 100000 -plumed jedi_monitor.dat
     if [ ! -f complex_nvt.gro ]; then continue; fi
     mv COLVAR nvt_COLVAR
     $gmx grompp -f ../../../../../mdp/nvt.mdp -c complex_nvt.gro -p complex.top -o complex_npt.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm complex_npt -nsteps 100000 -plumed jedi_monitor.dat
     if [ ! -f complex_npt.gro ]; then continue; fi
     mv COLVAR npt_COLVAR
     $gmx grompp -f ../../../../../mdp/md.mdp -c complex_npt.gro -p complex.top -o complex_md.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm complex_md -nsteps 25000000 -plumed jedi_monitor.dat
     if [ ! -f complex_md.gro ]; then continue; fi
    cd ..

    cd apoWet
     bash ../../../../../scripts/gmx5_setup.sh apoWet
     if [ ! -f apoWet_ions.gro ]; then continue; fi
     $gmx grompp -f ../../../../../mdp/em_10000_sd.mdp -c apoWet_ions.gro -p apoWet.top -o apoWet_em_sd.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm apoWet_em_sd -plumed jedi_monitor.dat -nsteps 10000
     if [ ! -f apoWet_em_sd.gro ]; then continue; fi
     mv COLVAR sd_COLVAR
     $gmx grompp -f ../../../../../mdp/em_10000_cg.mdp -c apoWet_em_sd.gro -p apoWet.top -o apoWet_em_cg.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm apoWet_em_cg -plumed jedi_monitor.dat -nsteps 10000
     if [ ! -f apoWet_em_cg.gro ]; then continue; fi
     mv COLVAR cg_COLVAR
     $gmx grompp -f ../../../../../mdp/nvt.mdp -c apoWet_em_cg.gro -p apoWet.top -o apoWet_nvt.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm apoWet_nvt -nsteps 100000 -plumed jedi_monitor.dat
     if [ ! -f apoWet_nvt.gro ]; then continue; fi
     mv COLVAR nvt_COLVAR
     $gmx grompp -f ../../../../../mdp/nvt.mdp -c apoWet_nvt.gro -p apoWet.top -o apoWet_npt.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm apoWet_npt -nsteps 100000 -plumed jedi_monitor.dat
     if [ ! -f apoWet_npt.gro ]; then continue; fi
     mv COLVAR npt_COLVAR
     $gmx grompp -f ../../../../../mdp/md.mdp -c apoWet_npt.gro -p apoWet.top -o apoWet_md.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm apoWet_md -nsteps 25000000 -plumed jedi_monitor.dat
     if [ ! -f apoWet_md.gro ]; then continue; fi
    cd ..

    cd apoDry
     bash ../../../../../scripts/gmx5_setup.sh apoDry
     if [ ! -f apoDry_ions.gro ]; then continue; fi
     $gmx grompp -f ../../../../../mdp/em_10000_sd.mdp -c apoDry_ions.gro -p apoDry.top -o apoDry_em_sd.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm apoDry_em_sd -plumed jedi_monitor.dat -nsteps 10000
     if [ ! -f apoDry_em_sd.gro ]; then continue; fi
     mv COLVAR sd_COLVAR
     $gmx grompp -f ../../../../../mdp/em_10000_cg.mdp -c apoDry_em_sd.gro -p apoDry.top -o apoDry_em_cg.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm apoDry_em_cg -plumed jedi_monitor.dat -nsteps 10000
     if [ ! -f apoDry_em_cg.gro ]; then continue; fi
     mv COLVAR cg_COLVAR
     $gmx grompp -f ../../../../../mdp/nvt.mdp -c apoDry_em_cg.gro -p apoDry.top -o apoDry_nvt.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm apoDry_nvt -nsteps 100000 -plumed jedi_monitor.dat
     if [ ! -f apoDry_nvt.gro ]; then continue; fi
     mv COLVAR nvt_COLVAR
     $gmx grompp -f ../../../../../mdp/nvt.mdp -c apoDry_nvt.gro -p apoDry.top -o apoDry_npt.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm apoDry_npt -nsteps 100000 -plumed jedi_monitor.dat
     if [ ! -f apoDry_npt.gro ]; then continue; fi
     mv COLVAR npt_COLVAR
     $gmx grompp -f ../../../../../mdp/md.mdp -c apoDry_npt.gro -p apoDry.top -o apoDry_md.tpr
     $gmx mdrun -ntomp $ntomp -v -deffnm apoDry_md -nsteps 25000000 -plumed jedi_monitor.dat
     if [ ! -f apoDry_md.gro ]; then continue; fi
    cd ..
  "


slurm="#!/bin/bash
#SBATCH --job-name=$pdb
#SBATCH -o JEDI2_MDtest.out
#SBATCH -e JEDI2_MDtest.err
#SBATCH -p serial 
#SBATCH -n 16 
#SBATCH -N 1
##SBATCH --gres=gpu:1 #this line is commented out in the cpu scripts
#SBATCH --time 48:00:00

"

rm -rfI output
cp -r input output
cd output

cd benchmark
for pdb in $(ls -d */); do
 cd $pdb

slurm="#!/bin/bash
#SBATCH --job-name=$pdb
#SBATCH -o JEDI2_MDtest.out
#SBATCH -e JEDI2_MDtest.err
#SBATCH -p serial 
#SBATCH -n 16 
#SBATCH -N 1
##SBATCH --gres=gpu:1 #this line is commented out in the cpu scripts
#SBATCH --time 48:00:00
module load gromacs/5.1.4
"

   echo "Processing system $pdb"
   echo "$slurm" > mdsub.sh
   echo "$md" >> mdsub.sh
 cd ..
done


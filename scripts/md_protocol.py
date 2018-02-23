import argparse, os,sys, copy, math
import mdtraj, numpy, gromacs

def split_pdb(frame,lignames,cofnames):
    '''
    Separate the structure frame of the complex
    in protein(s), ligand(s), water(s), and cofactor(s).
    A PDB file will be generated for the protein and
    the waters and ligands and cofactors.
    '''
    protein_water_list=[]
    ligands_list={}
    cofactors_list={}

    top=frame.topology

    # Separate atoms in protein, cofactor, ligand, water
    for i in range(0,frame.n_atoms):
        if i in top.select('protein'):
           protein_water_list.append(i)

        elif str(top.atom(i).residue)[0:3] in lignames:
           #FIXME This is going to put ligands that have same residue 
           #FIXME name and residue number (ligands in different chains 
           #FIXME of a dimer, for instance) in the same pdb file.
           if str(top.atom(i).residue) not in ligands_list.keys(): 
              ligands_list[str(top.atom(i).residue)]=[i]
           else:
              ligands_list[str(top.atom(i).residue)].append(i)

        elif str(top.atom(i).residue)[0:3] in cofnames:
           if str(top.atom(i).residue)[3:] not in ligands_list.keys():
              cofactors_list[str(top.atom(i).residue)]=[i]
           else:
              cofactors_list[str(top.atom(i).residue)].append(i)

        elif i in top.select('water'):
           protein_water_list.append(i)

  
    # Generate a new complex pdb file with the order: protein-cofactor-ligand-water (for JEDI purposes)
    all_atoms=[]
    for atom in protein_water_list:
        all_atoms.append(atom)
    for molecule in cofactors_list.keys():
        for atom in cofactors_list[molecule]:
            all_atoms.append(atom)
    for molecule in ligands_list.keys():
        for atom in ligands_list[molecule]:
            all_atoms.append(atom)
    system_renum=frame.atom_slice(all_atoms)
    system_renum.save_pdb('system_renum.pdb')

    system_std=frame.atom_slice(protein_water_list)
    system_std.save_pdb('system_std.pdb')    

    small_molec=[]
    for cof in cofactors_list.keys():
        namepdb=cof+'.pdb'
        small_molec.append(cof)
        cof_frame=frame.atom_slice(cofactors_list[cof])
        cof_frame.save_pdb()
        
    for lig in ligands_list.keys():
        namepdb=lig+'.pdb'
        small_molec.append(lig)
        lig_frame=frame.atom_slice(ligands_list[lig])
        lig_frame.save_pdb()

    return small_molec

def ligand_topologies(small_molec):

    for element in small_molec:
        cmd='babel -ipdb '+element+'.pdb -omol2 '+element+'.mol2' #FIXME There is a babel python wrapper but I don't know how to use it
        os.system(cmd)
        cmd='acpype -i '+element+'mol2 -c user'
        os.system(cmd) #FIXME Same as babel
        filecheck=element+'.acpype/'+element+'_GMX.gro'
        if not os.path.isfile(filecheck):
           print "Something went wrong with acpype. Exiting."
           sys.exit()


if __name__=='__main__':
     
     
     small_molec=split_pdb(frame,lignames,cofnames)



















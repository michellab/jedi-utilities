import argparse, os,sys, copy, math
import mdtraj, numpy, gromacs

def merge_molecules(bonds):
    out=[bonds[0].tolist()]
    for bond in bonds:
        blist=bond.tolist()
        for o in out:
            if set(blist).intersection(set(o)):
                o[:] = list(set(blist)) + list(set(o))
                break
            else:
                out.append(blist)
    return out 

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
    
    
    # Find separate molecules
    table,bonds=top.to_dataframe()
    '''
    molecules={0:[bonds[0][0],bonds[0][1]]}
    key_init=0
    for bond in bonds[1:]:
        new_molecule=True
        for key in molecules.keys():
            if bond[0] in molecules[key] and bond[1] not in molecules[key]:
               molecules[key].append(bond[1])
               new_molecule=False
               break
            elif bond[0] not in molecules[key] and bond[1] in molecules[key]:
               molecules[key].append(bond[0])
               new_molecule=False
               break
            elif bond[0] in molecules[key] and bond[1] in molecules[key]:
               new_molecule=False
               break
        if new_molecule==False:
           continue
        else:
           key_init=key_init+1
           molecules[key_init]=[bond[0],bond[1]]
    '''
    molecules=merge_molecules(bonds)
    print molecules
    sys.exit()

    # Separate atoms in protein, cofactor, ligand, water
    for i in range(0,frame.n_atoms):
        if i in top.select('protein'):
           protein_water_list.append(i)

        elif str(top.atom(i).residue)[0:3] in lignames:
           #FIXME This is going to put ligands that have same residue 
           #FIXME name and residue number (ligands in different chains 
           #FIXME of a dimer, for instance) in the same pdb file.
           #nameres=str(top.atom(i).residue)+str(top.atom(i).chain)
           #print nameres
           #sys.exit()
           if str(top.atom(i).residue) not in ligands_list.keys(): 
              ligands_list[str(top.atom(i).residue)]=[i]
           else:
              ligands_list[str(top.atom(i).residue)].append(i)
           
        elif str(top.atom(i).residue)[0:3] in cofnames:
           if str(top.atom(i).residue) not in cofactors_list.keys():
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
        cof_frame.save_pdb(namepdb)
        
    for lig in ligands_list.keys():
        namepdb=lig+'.pdb'
        small_molec.append(lig)
        lig_frame=frame.atom_slice(ligands_list[lig])
        lig_frame.save_pdb(namepdb)

    return small_molec

def ligand_topologies(small_molec):

    for element in small_molec:
        cmd='babel -ipdb '+element+'.pdb -omol2 '+element+'.mol2' #FIXME There is a babel python wrapper but I don't know how to use it
        os.system(cmd)
        cmd='/home/joan/github/acpype/acpype.py -i '+element+'.mol2'
        os.system(cmd) #FIXME Same as babel
        filecheck=element+'.acpype/'+element+'_GMX.gro'
        if not os.path.isfile(filecheck):
           print "Something went wrong with acpype. Exiting."
           sys.exit()

def find_oxygens(ligand):
    oxygen_coordinates=[]
    top_oxygen=mdtraj.Topology()
    chain=top_oxygen.add_chain()
    oxygen = mdtraj.element.oxygen
    top=ligand.topology
    for i in range(0,ligand.n_atoms):
        if top.atom(i).element==mdtraj.element.oxygen:
           residue = top_oxygen.add_residue("HOH",chain)
           top_oxygen.add_atom("O",oxygen,residue)
           oxygen_coordinates.append(ligand.xyz[0][i])
    xyz=numpy.array(oxygen_coordinates)
    oxygens=mdtraj.Trajectory(xyz,top_oxygen)
    oxygens.save("original_oxy.pdb")
    n_atoms=oxygens.n_atoms
    return oxygens,oxygen_coordinates,top_oxygen

def add_oxygens(ligand,oxygens,oxygen_coordinates,top_oxygen):
    noreplace=[]
    oxygen = mdtraj.element.oxygen
    chain=top_oxygen.add_chain()
    top=ligand.topology
    for i in range(0,ligand.n_atoms):
        if top.atom(i).element==mdtraj.element.oxygen:
           print "Found an old oxygen. Skipping"
           continue
        else:
           i_coords=ligand.xyz[0][i]
           for j in range(0,oxygens.n_atoms):
               j_coords=oxygens.xyz[0][j]
               dist2=((j_coords[0]-i_coords[0])**2+(j_coords[1]-i_coords[1])**2+(j_coords[2]-i_coords[2])**2)
               if dist2 < 0.0729:
                  noreplace.append(i)
                  break
           if i not in noreplace:
               residue = top_oxygen.add_residue("HOH",chain)
               top_oxygen.add_atom("O",oxygen,residue)
               oxygen_coordinates.append(ligand.xyz[0][i])
               xyz=numpy.array(oxygen_coordinates)
               oxygens=mdtraj.Trajectory(xyz,top_oxygen)
               oxygens=add_oxygens(ligand,oxygens,oxygen_coordinates,top_oxygen)

    return oxygens

def output_oxygens(oxygens,name):
    savename='HOH_'+name+'.pdb'
    oxygens.save(savename)
    return 0



if __name__=='__main__':
     
     frame=mdtraj.load_pdb('../input/1qzr.pdb')
     lignames=['CDX']
     cofnames=['ANP','MG'] 
     small_molec=split_pdb(frame,lignames,cofnames)
     ligand_topologies(small_molec)
     for name in small_molec:
         if name[0:3] in lignames:
            ligname=name+'.pdb'
            ligand=mdtraj.load_pdb(ligname)
            oxygens,oxygen_coordinates,top_oxygen=find_oxygens(ligand)
            oxygens=add_oxygens(ligand,oxygens,oxygen_coordinates,top_oxygen)
            z=output_oxygens(oxygens,name)


















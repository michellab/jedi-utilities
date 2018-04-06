import sys, numpy
import mdtraj

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





if __name__ == '__main__':
    
    ligand=mdtraj.load(sys.argv[1])
    ligand_name=sys.argv[1].split('.')[0]
    oxygens,oxygen_coordinates,top_oxygen=find_oxygens(ligand)
    print "Original oxygens: ", oxygens
    new_oxygens=add_oxygens(ligand,oxygens,oxygen_coordinates,top_oxygen)
    print "New oxygens: ", new_oxygens
    output_oxygens(new_oxygens,ligand_name)
    

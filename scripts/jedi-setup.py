#!/usr/bin/env python
#
""" 
 @AUTHOR Julien Michel. Sep 2015.
 @USAGE
  jedi-setup.py [-h] [-i [INPUT]] [-v] [-l [LIGAND]] [-r [REGION]]
                     [-a [APOLAR]] [-p [POLAR]] [-g [GRID]]
                     [-c [CUTOFF]] [-s [SPACING]]

Setup input files for a Gromacs/Plumed JEDI calculation.

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT], --input [INPUT]
                        Name of input PDB(GRO) file.
  -v, --version         Get version information about this script.
  -l [LIGAND], --ligand [LIGAND]
                        Optional ligand PDB to position the grid
  -r [REGION], --region [REGION]
                        Optional region file to define grid position
  -a [APOLAR], --apolar [APOLAR]
                        name of output apolar atoms file. Default apolar.pdb
  -p [POLAR], --polar [POLAR]
                        name of output polar atoms file. Default polar.pdb
  -g [GRID], --grid [GRID]
                        name of output grid file. Default grid.pdb
  -c [CUTOFF], --cutoff [CUTOFF]
                        the maximum distance (in nm) to extend the ligand's
                        bounding box
  -s [SPACING], --spacing [SPACING]
                        the spacing between neighboring grid points (in nm)
  -w [yes/no], --ignore_waters [yes/no]
                        specify if the water molecules have to be ignored
                        or not. Default is yes.
  -k [yes/no], --crop [yes/no]
                        delete (or not) the grid points that are further 
                        from the ligand than the cutoff. Default is yes.
  -f [yes/no], --full [yes/no]
                        generate a grid that covers the whole protein.
                        Default is no.
  -fp [POCKET] --fpocket [POCKET]
                        Generate the grid and the binding site from a
                        PDB file containing fpocket alpha spheres.
  -ln [LIGANDNAME] --ligandname [LIGANDNAME]
                        In liganded simulations, skip the residues with
                        the name of the ligand

jedi-setup.py is distributed under the GPL.

 @DEPENDENCIES:
 mdtraj, numpy

"""
#
# 02.11.15 FIXME) Add code to customise the created PDB files
# Need to:
# - retain original PDB indices. 
# - replace occupancy column with atomic masses
#
import argparse, os,sys, copy, math
import mdtraj, numpy 

parser = argparse.ArgumentParser(description="Setup input files for a Gromacs/Plumed JEDI calculation.",
                                 epilog="jedi-setup.py is distributed under the GPL.",
                                 prog="jedi-setup.py")
parser.add_argument('-i','--input', nargs="?",
                    help='Name of input PDB(GRO) file.')
parser.add_argument('-v','--version',action="store_true",
                    help='Get version information about this script.')
parser.add_argument('-l','--ligand', nargs="?",
                    help="Optional ligand PDB to position the grid")
parser.add_argument('-r','--region', nargs="?",
                    help="Optional region file to define grid position")
parser.add_argument('-a','--apolar', nargs="?",
                    help="name of output apolar atoms file. Default apolar.pdb",
                    default="apolar.pdb")
parser.add_argument('-p','--polar', nargs="?",
                    help="name of output polar atoms file. Default polar.pdb",
                    default="polar.pdb")
parser.add_argument('-g','--grid', nargs="?",
                    help="name of output grid file. Default grid.pdb",
                    default="grid.pdb")
parser.add_argument('-c','--cutoff',nargs="?",
                     default="0.50",
                     help="the maximum distance (in nm) to extend the ligand's bounding box")
parser.add_argument('-s','--spacing', nargs="?",
                    default="0.15",
                    help="the spacing between neighboring grid points (in nm)")
parser.add_argument('-w', '--ignore_waters',nargs="?",
                    default="yes",
                    help="omit water oxygens in apolar atoms file")
parser.add_argument('-k', '--crop',nargs="?",
                    default="yes",
                    help="delete grid points that are further from the ligand than the cutoff")
parser.add_argument('-f', '--full',nargs="?",
                    default="no",
                    help="Generate a grid that overlaps the whole protein")
parser.add_argument('-fp', '--fpocket',nargs="?",
                    help="Generate the grid and the binding site from a PDB file containing fpocket alpha spheres.")
parser.add_argument('-ln', '--ligandname',nargs="?",
                    help="Skip the residues that have the name of the ligand")

def parse(parser):
    args = parser.parse_args()

    if args.input is None:
        parser.print_help()
        print ("""@@@
ERROR ! You must specificy an input PDB. See usage above. 
@@@""")
        sys.exit(-1)

    if not os.path.isfile(args.input):
        parser.print_help()
        print ("""@@@
ERROR ! The input PDB cannot be found. See usage above.
@@@""")
        sys.exit(-1)

    # Catch lack of ligand and region file if the -f/--full flag is not present
    if ( (args.ligand is None and args.region is None and args.fpocket is None) or 
         (
          (args.ligand is not None and args.region is not None) or 
          (args.fpocket is not None and args.region is not None) or 
          (args.ligand is not None and args.fpocket is not None) 
         ) 
       ) and (args.full=="no"):
        parser.print_help()
        print ("""@@@
ERROR ! You must specify a ligand, a region file or fpocket alpha spheres as input.
@@@""")
        sys.exit(-1)

    # Catch specified but absent ligand file
    if not args.ligand is None and not os.path.isfile(args.ligand):
        parser.print_help()
        print ("""@@@
ERROR ! The specified input ligand file cannot be found. See usage above.
@@@""")
        sys.exit(-1)

    # Catch specified but absent region file
    if not args.region is None and not os.path.isfile(args.region):
        parser.print_help()
        print ("""@@@
ERROR ! The specified input region file cannot be found. See usage above.
@@@""")
        sys.exit(-1)

    if args.fpocket is not None:
       if not os.path.isfile(args.fpocket):
          print "File "+args.fpocket+" could not be found. Exiting"
          sys.exit()
       print "Grid and binding site are going to be generated from fpocket alpha spheres from file "+args.fpocket


    # Inform the user of whether water oxygens are omitted or not
    if (args.ignore_waters is None) or (args.ignore_waters == "yes"):
       print "Water oxygens are going to be ignored"
    elif args.ignore_waters == "no":
       print "Water oxygens are NOT going to be ignored"
    else:
       print """@@@
ERROR ! The -w flag can only have "yes" or "no". Default is "yes".
@@@"""
       sys.exit(-1)

    # Inform the user of whether the grid is going to be cropped or not
    if not (args.crop is None or args.crop == "yes" or args.crop=="no"):
       print """@@@
ERROR ! The -k flag can only have "yes" or "no". Default is "yes".
@@@"""
       sys.exit(-1)
       
    elif ((args.crop is None) or (args.crop == "yes")) and (args.full=="no") and args.ligand is not None:
       print "The grid points further than ", args.cutoff, " nm from the ligand are going to be deleted."
    elif (args.crop == "no") and (args.full=="no"):
       print "All grid points are going to be kept."
    elif args.full=="yes":
       print "Overlapping a grid onto the whole protein"

    if args.full=="yes":
       print "The grid is going to overlap the whole protein"
   
    print (args)
    return args.input, args.ligand, float(args.cutoff), args.region, float(args.spacing),\
        args.apolar,args.polar,args.grid,args.ignore_waters,args.crop,args.full,args.fpocket,\
        args.ligandname


def loadStructure(pdbfile,origin=None):
    """Input: pdbfile: A pdb file name and the origin from which the grid will be defined (optional)
    Output: A datastructure that holds a pdb frame and (optional) a list with necessary parameters as bfactor
    """
    struc = mdtraj.load(pdbfile)
    # Remove the water molecules unless the specifies otherwise
    if wat == 'yes': 
       struc=struc.remove_solvent()
    if origin=='fpocket':
       filein=open(pdbfile,'r')
       bfact=[]
       for line in filein:
           if line.startswith('ATOM') or line.startswith('HETATM'):
              bfact.append(float(line.split()[9])/10)
       filein.close()
    else:
       bfact=None
    return struc,bfact

def loadRegion(regionfile):
    """
    Input: regionfile: a file that contains the definiton of a rectangular box
    Output: a datastructure defining the bounding box
    """
    region={'min':[],'max':[]}
    filein=open(regionfile,'r')
    for line in filein:
        if line.startswith('#') or len(line.split())!=3:
           continue
        line=line.split()    
        region['min'].append(float(line[1]))
        region['max'].append(float(line[2]))
    filein.close()
    return region 

#def selectAlignment(frame, rule="backbone"):
#    """Input: frame: a mdtraj frame
#              rule (optional): How to select atoms to retain from frame.
#              rule can be 'backbone' == name C, CA, N, O, HA, H
#              rule can be 'Calpha' == name CA
#              Default is backbone
#    """
#    subset = frame.topology.select(rule)
#    alignment = frame.atom_slice(subset)
#    return alignment, subset

def defineGrid(frame, full, ligand=None, lig_cutoff=5.0, region=None, spacing=0.15, fpocket=None, bfact=None):
    """Input: frame: a mdtraj frame
              ligand (optional): a mdtraj frame with ligand coordinates
              lig_cutoff (optional): the maximal distance between ligand
              atoms and grid points. In angstroms.
              region (optional): a txt file with min/max coordinates of a
              paralleliped.
              fpocket (optional): an mdtraj frame containing the coordinates of the alpha spheres
       Output:
              grid: a mdtraj frame
       Throw exception if ligand and region and fpocket are all None
    """
    # Build bounding box around ligand. If ligand is none, use region file
    #
    # !!!NOPBC CHECK !!!
    #
    mincoords = [99999.0, 99999.0, 99999.0]
    maxcoords = [-99999.0, -99999.0, -99999.0]
    #min_crd=99999.0
    #max_crd=-99999.0
#    print "ligand ", ligand
#    print "region " ,region
#    print "fpocket ", fpocket
#    sys.exit()
    if full=="no":
        if (ligand is not None):
            for i in range(0,ligand.n_atoms):
               # min_crd_at=min(ligand.xyz[0][i])
               # if min_crd_at<min_crd:
               #    min_crd=min_crd_at
               # max_crd_at=max(ligand.xyz[0][i])
               # if max_crd_at>max_crd:
               #    max_crd=max_crd_at
            #mincoords=[min_crd]*3
            #maxcoords=[max_crd]*3
                atcoord = ligand.xyz[0][i]
                #print (atcoord)
                if atcoord[0] < mincoords[0]:
                    mincoords[0] = atcoord[0]
                if atcoord[1] < mincoords[1]:
                    mincoords[1] = atcoord[1]
                if atcoord[2] < mincoords[2]:
                    mincoords[2] = atcoord[2]
                if atcoord[0] > maxcoords[0]:
                    maxcoords[0] = atcoord[0]
                if atcoord[1] > maxcoords[1]:
                    maxcoords[1] = atcoord[1]
                if atcoord[2] > maxcoords[2]:
                    maxcoords[2] = atcoord[2]
            # Now extend bounding box by lig_cutoff
            for i in range(0,3):
                mincoords[i] = mincoords[i] - lig_cutoff
                maxcoords[i] = maxcoords[i] + lig_cutoff
        elif region is not None:
            mincoords = region['min']
            maxcoords = region['max']
        elif (fpocket is not None):
            for i in range(0,fpocket.n_atoms):
                atcoord = fpocket.xyz[0][i]
                #print (atcoord)
                if atcoord[0] < mincoords[0]:
                    mincoords[0] = atcoord[0]
                if atcoord[1] < mincoords[1]:
                    mincoords[1] = atcoord[1]
                if atcoord[2] < mincoords[2]:
                    mincoords[2] = atcoord[2]
                if atcoord[0] > maxcoords[0]:
                    maxcoords[0] = atcoord[0]
                if atcoord[1] > maxcoords[1]:
                    maxcoords[1] = atcoord[1]
                if atcoord[2] > maxcoords[2]:
                    maxcoords[2] = atcoord[2]
            # Now extend bounding box by lig_cutoff
            for i in range(0,3):
                mincoords[i] = mincoords[i] - max(bfact) 
                maxcoords[i] = maxcoords[i] + min(bfact)

    elif full=="yes":
        for i in range(0,frame.n_atoms):
            atcoord=frame.xyz[0][i]
            if atcoord[0] < mincoords[0]:
                mincoords[0] = atcoord[0]
            if atcoord[1] < mincoords[1]:
                mincoords[1] = atcoord[1]
            if atcoord[2] < mincoords[2]:
                mincoords[2] = atcoord[2]
            if atcoord[0] > maxcoords[0]:
                maxcoords[0] = atcoord[0]
            if atcoord[1] > maxcoords[1]:
                maxcoords[1] = atcoord[1]
            if atcoord[2] > maxcoords[2]:
                maxcoords[2] = atcoord[2]
        # extending the grid by 0.3 nm so it extends a bit from the protein
        for i in range(0,3): 
            mincoords[i] = mincoords[i] - 0.3 
            maxcoords[i] = maxcoords[i] + 0.3

    #print ("Mininum grid coordinates: %s " % mincoords)
    #print ("Maximum grid coordinates: %s " % maxcoords)
    # Populate paralleliped with evenly spaced grid points
    # work out how many N grid points to add
    xstep = int ( (maxcoords[0] - mincoords[0] )/ spacing ) + 1
    ystep = int ( (maxcoords[1] - mincoords[1] )/ spacing ) + 1
    zstep = int ( (maxcoords[2] - mincoords[2] )/ spacing ) + 1
    print ("Number of grid points along each axis (before cropping): %s %s %s " % (xstep,ystep,zstep))
    Ngrid = xstep*ystep*zstep
    print ("Total number of grid points (before cropping): %s " % Ngrid)
    
    # Generate a list of coordinates of grid points
    coords = []
    for i in range(0,xstep):
        for j in range(0,ystep):
            for k in range(0,zstep):
                c_x = mincoords[0] + i * spacing
                c_y = mincoords[1] + j * spacing
                c_z = mincoords[2] + k * spacing
                coords.append( [ c_x, c_y, c_z ] )
    
    # Delete (if required) the grid points that are too far from the ligand
    if (ligand is not None) and (crop=='yes'):
       cropcoords=[]
       for gridpoint in coords:
           for atom in ligand.xyz[0]:
               distance=math.sqrt((atom[0]-gridpoint[0])**2+\
                                  (atom[1]-gridpoint[1])**2+\
                                  (atom[2]-gridpoint[2])**2)
               if distance<=lig_cutoff:
                  cropcoords.append(gridpoint)
                  break
       coords=cropcoords    
       Ngrid=len(coords)
       print ("Total number of grid points (after cropping): %s " % Ngrid)
    
    elif fpocket is not None and crop=='yes':
       cropcoords=[]
       for gridpoint in coords:
           for sphere in fpocket.xyz[0]:
               i=int(numpy.where(fpocket.xyz[0]==sphere)[0][0])
               radius=bfact[i]
               distance=math.sqrt((sphere[0]-gridpoint[0])**2+\
                                  (sphere[1]-gridpoint[1])**2+\
                                  (sphere[2]-gridpoint[2])**2)
               if distance<=radius:
                  cropcoords.append(gridpoint)
                  break
       coords=cropcoords
       Ngrid=len(coords)
       print ("Total number of grid points (after cropping): %s " % Ngrid)

    elif (full=="yes"):
       cropcoords=[]
       for gridpoint in coords:
           for atom in frame.xyz[0]:
               distance=math.sqrt((atom[0]-gridpoint[0])**2+\
                                  (atom[1]-gridpoint[1])**2+\
                                  (atom[2]-gridpoint[2])**2)
               if distance<=1.0:
                  cropcoords.append(gridpoint)
                  break
       coords=cropcoords
       Ngrid=len(coords)
       print ("Total number of grid points (after cropping): %s " % Ngrid)


    top = mdtraj.Topology()
    # build one mdtraj.topology.chain
    chain = top.add_chain()
    # build one mdtraj.topology.residue
    residue = top.add_residue("GRI", chain)
    # build N mdtraj.topology.atoms
    hydrogen = mdtraj.element.hydrogen
    for i in range(0,Ngrid):
        top.add_atom("GRI", hydrogen, residue)
    # Loop over every atom in the topology and create a xyz np array
    xyz = numpy.array(coords)#, dfloat=float32)
    # Find a way to build a mdtraj.traj from top + xyz
    grid = mdtraj.Trajectory(xyz, top)
    return grid, mincoords, maxcoords

def selectPolarApolar(frame, grid_min, grid_max, ligand=None,ligname):
    """Input: frame: a mdtraj frame
              grid_min: minimum grid coordinates
              grid_max: maximim grid coordinates
       Output:
              polar: a mdtraj frame containing polar atoms from system
              apolar: a mdtraj frame containing apolar atoms from system
    """
    #
    #!!! NO PBC CHECK
    #
    polar_list = []
    apolar_list = []
    for i in range(0,frame.n_atoms):
        if frame.topology.atom(i).residue.name==ligname:
           continue
        atcoord = frame.xyz[0][i]
        if (atcoord[0] > grid_min[0] and
            atcoord[0] < grid_max[0] and
            atcoord[1] > grid_min[1] and
            atcoord[1] < grid_max[1] and
            atcoord[2] > grid_min[2] and 
            atcoord[2] < grid_max[2] ):
            if (ligand is not None) and (crop=='yes'): # Crop the apolar and polar files if requested
               for ligatom in ligand.xyz[0]:
                   distance=math.sqrt((ligatom[0]-atcoord[0])**2+\
                                      (ligatom[1]-atcoord[1])**2+\
                                      (ligatom[2]-atcoord[2])**2)
                   if distance<=lig_cutoff:
                       if frame.topology.atom(i).element.symbol in ['N','O']:
                          polar_list.append(i)
                       elif frame.topology.atom(i).element.symbol in ['C','S']:
                          apolar_list.append(i)
                       break
            else:
               if (frame.topology.atom(i).element.symbol in ['N','O']):
                   polar_list.append(i)
               elif (frame.topology.atom(i).element.symbol in ['C','S']):
                   apolar_list.append(i)
    #print polar_list
    #print apolar_list
    polar = system.atom_slice(polar_list)
    apolar = system.atom_slice(apolar_list)
    return polar, polar_list, apolar, apolar_list

def outputPdb(frame, outfile="output.pdb", indices=None):
    """
    Input: frame: a mdtraj frame
           outfile (optional): the name of the pdb file to write frame to
    """
    frame.save(outfile)
    # Not pretty ! Ideally would modify mdtraj's API
    if indices is not None:
        rstream = open(outfile)
        buffer = rstream.readlines()
        rstream.close()
        wstream = open('temp.pdb','w')
        pdb_idx = 0
        for line in buffer:
            if line.startswith("ATOM"):
                new_idx = indices[pdb_idx] + 1 # Because starts at 0
                mass = frame.topology.atom(pdb_idx).element.mass
                radius = frame.topology.atom(pdb_idx).element.radius
                #print pdb_idx, new_idx, mass
                if new_idx < 10:
                    str_idx = '    %s' % new_idx
                elif new_idx < 100:
                    str_idx= '   %s' % new_idx
                elif new_idx < 1000:
                    str_idx = '  %s' % new_idx
                elif new_idx < 10000:
                    str_idx = ' %s' % new_idx
                else:
                    str_idx = '%s' % new_idx
                str_mass = '%-5.2f' % mass
                str_radius = '%-5.2f' % radius
                new_line = line[0:6] + str_idx + line[11:55] + str_mass + '  ' + str_radius + line[66:] 
                #import pdb; pdb.set_trace()
                wstream.write(new_line)
                pdb_idx +=1
                #print "replaced index ",pdb_idx," with ", str_idx
            else:
                wstream.write(line)
        wstream.close()
        cmd = "mv temp.pdb %s" % outfile
        os.system(cmd)
        #sys.exit(-1)
    # If indices is not None
    # load file
    # for every ATOM line
    # update index
    # also update Occupancy column to have atomic mass

def centerGrid(grid_data, polar, apolar):
    # Compute com polar/apolar
    com = [0.0, 0.0, 0.0]
    tot_mass = 0.0
    idx = 0
    for coords in polar.xyz[0]:
        mass = polar.topology.atom(idx).element.mass
        for i in range(0,3):
            com[i] += coords[i]*mass
        tot_mass += mass
        idx += 1
    idx = 0
    for coords in apolar.xyz[0]:
        mass = apolar.topology.atom(idx).element.mass
        for i in range(0,3):
            com[i] += coords[i]*mass
        tot_mass += mass
        idx += 1
    for i in range(0,3):
        com[i] /= tot_mass
    #print (com)
    # Compute cog grid
    cog = [0.0,0.0,0.0]
    idx = 0
    for coords in grid_data[0].xyz[0]:
        for i in range(0,3):
            cog[i] += coords[i]
        idx += 1
    for i in range(0,3):
        cog[i] /= idx
    #print (cog)
    #print (grid_data[0].xyz[0][0])
    # Update grid coordinates by delta_cog
    delta_cog = [com[0]-cog[0], com[1]-cog[1], com[2]-cog[2]]
    #print (delta_cog)
    #delta_cog = [0.0, 0.0, 0.0]

    idx=0
    for coords in grid_data[0].xyz[0]:
        newcoord = coords
        for i in range(0,3):
            newcoord[i] = coords[i]+delta_cog[i]
        grid_data[0].xyz[0][idx] = newcoord
        idx += 1
    #print (grid_data[0].xyz[0][0])
    #import pdb; pdb.set_trace()

if __name__ == '__main__':
    print ("*** jedi setup beginning *** ")
    # Parse command line arguments
    system_pdb, ligand_pdb, lig_cutoff, region_dim, spacing,\
        apolar_pdb, polar_pdb, grid_pdb, wat, crop, full, fpocket_pdb, ligname = parse(parser)
    # Load protein coordinates
    system,bfact = loadStructure(system_pdb)

    # Load ligand coordinates/region definition if necessary
    if full=="no":
        if ligand_pdb is not None:
            ligand,bfact = loadStructure(ligand_pdb)
        else:
            ligand = None

        if region_dim is not None:
            region = loadRegion(region_dim)
        else:
            region = None
        
        if fpocket_pdb is not None:
            fpocket,bfact=loadStructure(fpocket_pdb,"fpocket")
        else:
            fpocket=None
    elif full=="yes":
        ligand=None
        region=None
        fpocket=None
    else:
        print "Something went wrong. Check the code."
        sys.exit()
    # construct alignment
    #alignment, alignment_indices = selectAlignment(system, rule="backbone")
    # construct grid, polar and apolar
    print "System ", system
    if full=="no":
       print "ligand ", ligand
       print "lig_cutoff ", lig_cutoff
       print "region: ", region, " Type ", type(region)
    print "spacing ", spacing
    #sys.exit("This is just a test, modifications done between line 367 and this point")
    grid_data = defineGrid(system,full,ligand=ligand, lig_cutoff=lig_cutoff,\
                                   region=region, spacing=spacing,fpocket=fpocket,bfact=bfact)
    polar, polar_indices, apolar, apolar_indices =\
        selectPolarApolar(system, grid_data[1], grid_data[2],ligand,ligname)
    # Now center grid on COM of polar+apolar region
    centerGrid(grid_data, polar, apolar)
    #sys.exit(-1)
    # Write output
    outputPdb(grid_data[0], outfile=grid_pdb, indices=None)
    outputPdb(polar, outfile=polar_pdb, indices=polar_indices)
    outputPdb(apolar, outfile=apolar_pdb, indices=apolar_indices)
    #outputPdb(alignment, outfile=alignment_pdb, indices=alignment_indices)
    print ("*** jedi setup complete *** ")

#!/usr/bin/env python
#
""" 
 @AUTHOR Julien Michel. Sep 2015.
 @USAGE
  jedi-setup.py [-h] [-i [INPUT]] [-v] [-l [LIGAND]] [-r [REGION]]
                     [-a [APOLAR]] [-p [POLAR]] [-n [ALIGNMENT]] [-g [GRID]]
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
  -n [ALIGNMENT], --alignment [ALIGNMENT]
                        name of output alignment atoms file. Default
                        alignment.pdb
  -g [GRID], --grid [GRID]
                        name of output grid file. Default grid.pdb
  -c [CUTOFF], --cutoff [CUTOFF]
                        the maximum distance (in nm) to extend the ligand's
                        bounding box
  -s [SPACING], --spacing [SPACING]
                        the spacing between neighboring grid points (in nm)

jedi-setup.py is distributed under the GPL.

 @DEPENDENCIES:
 mdtraj, numpy

"""
import argparse, os,sys, copy
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
parser.add_argument('-n','--alignment', nargs="?",
                    help="name of output alignment atoms file. Default alignment.pdb",
                    default="alignment.pdb")
parser.add_argument('-g','--grid', nargs="?",
                    help="name of output grid file. Default grid.pdb",
                    default="grid.pdb")
parser.add_argument('-c','--cutoff',nargs="?",
                     default="0.50",
                     help="the maximum distance (in nm) to extend the ligand's bounding box")
parser.add_argument('-s','--spacing', nargs="?",
                    default="0.15",
                    help="the spacing between neighboring grid points (in nm)")

def parse(parser):
    args = parser.parse_args()

    if args.input is None:
        parser.print_help()
        print """@@@
ERROR ! You must specificy an input PDB. See usage above. 
@@@"""
        sys.exit(-1)

    if not os.path.isfile(args.input):
        parser.print_help()
        print """@@@
ERROR ! The input PDB cannot be found. See usage above.
@@@"""
        sys.exit(-1)

    # Catch lack of ligand and region file
    if ( (args.ligand is None and args.region is None) or
         (not args.ligand is None and not args.region is None) ):
        parser.print_help()
        print """@@@
ERROR ! You must specify a ligand OR a region file as input.
@@@"""
        sys.exit(-1)

    # Catch specified but absent ligand file
    if not args.ligand is None and not os.path.isfile(args.ligand):
        parser.print_help()
        print """@@@
ERROR ! The specified input ligand file cannot be found. See usage above.
@@@"""
        sys.exit(-1)

    # Catch specified but absent region file
    if not args.region is None and not os.path.isfile(args.region):
        parser.print_help()
        print """@@@
ERROR ! The specified input region file cannot be found. See usage above.
@@@"""
        sys.exit(-1)


    print args
    return args.input, args.ligand, float(args.cutoff), args.region, float(args.spacing),\
        args.apolar,args.polar, args.alignment, args.grid

def loadStructure(pdbfile):
    """Input: pdbfile: A pdb file name
    Output: A datastructure that holds a pdb frame
    """
    struc = mdtraj.load(pdbfile)

    return struc

def loadRegion(regionfile):
    """
    Input: regionfile: a file that contains the definiton of a rectangular box
    Output: a datastructure defining the bounding box
    """
    return 0

def selectAlignment(frame, rule="backbone"):
    """Input: frame: a mdtraj frame
              rule (optional): How to select atoms to retain from frame.
              rule can be 'backbone' == name C, CA, N, O, HA, H
              rule can be 'Calpha' == name CA
              Default is backbone
    """
    subset = frame.topology.select(rule)
    alignment = frame.atom_slice(subset)
    return alignment

def defineGrid( frame, ligand=None, lig_cutoff=5.0, region=None, spacing=0.15):
    """Input: frame: a mdtraj frame
              ligand (optional): a mdtraj frame with ligand coordinates
              lig_cutoff (optional): the maximal distance between ligand
              atoms and grid points. In angstroms.
              region (optional): a txt file with min/max coordinates of a
              paralleliped.
       Output:
              grid: a mdtraj frame
       Throw exception if ligand and region are both None
    """
    # Build bounding box around ligand. If ligand is none, use region file
    #
    # !!!NOPBC CHECK !!!
    #
    mincoords = [99999.0, 99999.0, 99999.0]
    maxcoords = [-99999.0, -99999.0, -99999.0]
    if ligand is not None:
        for i in range(0,ligand.n_atoms):
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
        for i in range(0,2):
            mincoords[i] = mincoords[i] - lig_cutoff
            maxcoords[i] = maxcoords[i] + lig_cutoff
    else:
        pass
        #mincoords = region['min']
        #maxcoords = region['max']
    print ("Mininum grid coordinates: %s " % mincoords)
    print ("Maximum grid coordinates: %s " % maxcoords)
    # Populate paralleliped with evenly spaced grid points
    # work out how many N grid points to add
    xstep = int ( (maxcoords[0] - mincoords[0] )/ spacing ) + 1
    ystep = int ( (maxcoords[1] - mincoords[1] )/ spacing ) + 1
    zstep = int ( (maxcoords[2] - mincoords[2] )/ spacing ) + 1
    print ("Number of grid points along each axis : %s %s %s " % (xstep,ystep,zstep))
    Ngrid = xstep*ystep*zstep
    print ("Total number of grid points: %s " % Ngrid)
    top = mdtraj.Topology()
    # build one mdtraj.topology.chain
    chain = top.add_chain()
    # build one mdtraj.topology.residue
    residue = top.add_residue("GRI", chain)
    # build N mdtraj.topology.atoms
    germanium = mdtraj.element.germanium
    for i in range(0,Ngrid):
        top.add_atom("GRI", germanium, residue)
    # Loop over every atom in the topology and create a xyz np array
    coords = []
    for i in range(0,xstep):
        for j in range(0,ystep):
            for k in range(0,zstep):
                c_x = mincoords[0] + i * spacing
                c_y = mincoords[1] + j * spacing
                c_z = mincoords[2] + k * spacing
                coords.append( [ c_x, c_y, c_z ] )
    xyz = numpy.array(coords)#, dfloat=float32)
    # Find a way to build a mdtraj.traj from top + xyz
    grid = mdtraj.Trajectory(xyz, top)
    return grid, mincoords, maxcoords

def selectPolarApolar(frame, grid_min, grid_max):
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
        atcoord = frame.xyz[0][i]
        if (atcoord[0] > grid_min[0] and
            atcoord[0] < grid_max[0] and
            atcoord[1] > grid_min[1] and
            atcoord[1] < grid_max[1] and
            atcoord[2] > grid_min[2] and 
            atcoord[2] < grid_max[2] ):
            if frame.topology.atom(i).element.symbol in ['N','O']:
                polar_list.append(i)
            elif frame.topology.atom(i).element.symbol in ['C','S']:
                apolar_list.append(i)
    #print polar_list
    #print apolar_list
    polar = system.atom_slice(polar_list)
    apolar = system.atom_slice(apolar_list)
    return polar, apolar

def outputPdb(frame, outfile="output.pdb"):
    """
    Input: frame: a mdtraj frame
           outfile (optional): the name of the pdb file to write frame to
    """
    frame.save(outfile)


if __name__ == '__main__':
    print ("*** jedi setup beginning *** ")
    # Parse command line arguments
    system_pdb, ligand_pdb, lig_cutoff, region_dim, spacing,\
        apolar_pdb, polar_pdb, alignment_pdb, grid_pdb = parse(parser)
    # Load protein coordinates
    system = loadStructure(system_pdb)

    # Load ligand coordinates/region definition
    if ligand_pdb is not None:
        ligand = loadStructure(ligand_pdb)
    else:
        ligand = None

    if region_dim is not None:
        region = loadRegion(region_dim)
    else:
        region = None
    # construct alignment
    # FIXME HOW TO RETAIN ORIGINAL PDB INDICES ?
    # FIXME HOW TO PUT ATOMIC MASSES IN OCCUPANCY COLUMN
    alignment = selectAlignment(system, rule="backbone")
    # construct grid, polar and apolar
    grid_data = defineGrid(system, ligand=ligand, lig_cutoff=lig_cutoff,\
                               region=region, spacing=spacing)
    # FIXME HOW TO RETAIN ORIGINAL PDB INDICES ?
    polar, apolar = selectPolarApolar(system, grid_data[1], grid_data[2])
    #sys.exit(-1)
    # Write output
    outputPdb(grid_data[0], outfile=grid_pdb)
    outputPdb(polar, outfile=polar_pdb)
    outputPdb(apolar, outfile=apolar_pdb)
    outputPdb(alignment, outfile=alignment_pdb)
    print ("*** jedi setup complete *** ")

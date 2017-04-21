import argparse, os,sys, copy, math
import mdtraj, numpy

parser = argparse.ArgumentParser(description="",
                                 epilog="",
                                 prog="")
parser.add_argument('-i','--input_grid', nargs="?",
                    help='Name of the grid xyz file')
parser.add_argument('-a','--input_activities',nargs="?",
                    help='Name of the input_activities file')
parser.add_argument('-o','--output_grid',nargs="?",
                    help='Output pdb file with the input_activities as the occupancy')
parser.add_argument('-k', '--crop',nargs="?",
                    help="Delete grid points with activity 0")

def parse(parser):
     args = parser.parse_args()
     if args.input_grid is None:
        parser.print_help()
        print ("""@@@
ERROR ! You must specificy an xyz grid file and an input_activities file. See usage above. 
@@@""")
        sys.exit(0)

     if not os.path.isfile(args.input_grid):
        parser.print_help()
        print ("""@@@
ERROR ! The input grid xyz file cannot be found. See usage above.
@@@""")
        sys.exit(-1)

     if not os.path.isfile(args.input_activities):
        parser.print_help()
        print ("""@@@
ERROR ! The input_activities file cannot be found. See usage above.
@@@""")
        sys.exit(-1)

     if args.output_grid is None:
        parser.print_help()
        print ("""@@@
ERROR ! Grid output file must be specified. See usage above.
@@@""")
        sys.exit(-1)

     if args.crop is None or args.crop=="yes":
        args.crop="yes"
        print "grid points with activity 0 are going to be ignored" 
     elif args.crop=="no":
        print "Keeping all grid points"
     else:
        print 'the -k/--crop flag only accepts "yes" or "no. Default is "yes".'

     return args.input_grid, args.input_activities, args.output_grid, args.crop

def getGridCoordinates(input_grid):
    coords=[]
    grid_file=open(input_grid,'r')
    numline=0
    num_points=0
    for line in grid_file:
        numline=numline+1
        line=line.split()
        if numline==1:
           num_points=int(line[0])
        elif numline==2:
           continue
        else:
           coords.append([float(line[1]),float(line[2]),float(line[3])])
    grid_file.close()
    return coords,num_points

def getActivities(input_activities):
    activities=[]
    activities_file=open(input_activities,'r')
    for line in activities_file:
        activities.append(float(line.split()[0]))
    activities_file.close()
    return activities

def getGridObject(coords,activities,num_points):
    top = mdtraj.Topology()
    chain = top.add_chain()
    residue = top.add_residue("GRI", chain)
    germanium = mdtraj.element.germanium
    for i in range(0,num_points):
        top.add_atom("GRI", germanium, residue)
    xyz = numpy.array(coords)
    grid_object = mdtraj.Trajectory(xyz, top)
    return grid_object
    
def outputPdb(grid_object,activities,output_grid):
    grid_object.save('temp.pdb')
    tempfile=open('temp.pdb','r')
    outfile=open(output_grid,'w')
    for line in tempfile:
        if line.startswith('ATOM'):
           pdb_index=int(line.split()[1])
           activities_index=pdb_index-1 # indices in pdb files start at 1
           activity_pdb='{0:.2f}'.format(activities[activities_index]) #Occuppancies in pdb files only have 2 decimals
           new_line=line[0:56]+str(activity_pdb)+line[60:]
           outfile.write(new_line)
    tempfile.close()
#    os.system('rm temp.pdb')
    outfile.close()
    return 0

if __name__ == '__main__':
   
   input_grid,input_activities,output_grid,crop=parse(parser)
   coords,num_points=getGridCoordinates(input_grid)
   activities=getActivities(input_activities)
   grid_object=getGridObject(coords,activities,num_points)
   outputPdb(grid_object,activities,output_grid) 
   

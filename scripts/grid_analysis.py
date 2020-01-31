#!/usr/bin/env python
#
help_message=\
""" 
 @AUTHOR Joan Clark-Nicolas. May 2017.
 @USAGE  python grid_analysis.py jedi-setup.py [-h] [-i [grid xyz input]] [-a [activities txt input]] 
                                               [-o [grid pdb output]]
       
Generate a PDB file of the JEDI grid with the activities of the grid points as occupancy.

optional arguments:
  -h, --help            Show this help message and exit. 
  -i [INPUT], --input_grid [INPUT]
                        Name of input xyz grid file.
  -a [INPUT], --input_activities [INPUT]
                        Name of input txt activities file
  -o [OUTPUT], '--output_grid [OUTPUT]
                        Name of output pdb grid file. The occupancy column has been replaced with the activities.
  -k [yes/no], --crop [yes/no]
                        delete (or not) the grid points that have null activity. Default is yes.
grid_analysis.py is distributed under the GPL.

 @DEPENDENCIES:
 mdtraj, numpy

"""

import argparse, os,sys, copy, math
import mdtraj, numpy

parser = argparse.ArgumentParser(description="",
                                 epilog="",
                                 prog="")
parser.add_argument('-i','--input_grid', nargs="?",
                    help='Name of the grid xyz file')
parser.add_argument('-a','--input_activities',nargs="?",
                    help='Name of the input_activities file')
parser.add_argument('-p','--property',nargs="?", default="activity",
                    help='Property to print as occupancy (activity, lig_i, Contact_i, Exposure_i)')
parser.add_argument('-o','--output_grid',nargs="?",
                    help='Output pdb file with the input_activities as the occupancy')
parser.add_argument('-k', '--crop',nargs="?",
                    help="Delete grid points with activity 0")

def parse(parser):
     args = parser.parse_args()

     if args.input_grid is None:
        print ("""@@@
ERROR ! You must specificy an xyz grid file and an input_activities file. See usage above. 
@@@""")
        sys.exit(0)

     if not os.path.isfile(args.input_grid):
        print (help_message)
        print ("""@@@
ERROR ! The input grid xyz file cannot be found. See usage above.
@@@""")
        sys.exit(-1)

     if not os.path.isfile(args.input_activities):
        print (help_message)
        print ("""@@@
ERROR ! The input_activities file cannot be found. See usage above.
@@@""")
        sys.exit(-1)

     if args.output_grid is None:
        print (help_message)
        print ("""@@@
ERROR ! Grid output file must be specified. See usage above.
@@@""")
        sys.exit(-1)

     if args.crop is None or args.crop=="yes":
        args.crop="yes"
        print ("grid points with activity 0 are going to be ignored") 
     elif args.crop=="no":
        print ("Keeping all grid points")
     else:
        print ('the -k/--crop flag only accepts "yes" or "no. Default is "yes".')

     properties=["activity","lig","contact","exposure","all"]
     if (args.property) not in properties:
        print ("@@@\n\
         ERROR ! -p/--property must have one of the following values: "+",".join(properties)+". Default is activity.\n"+\
         "@@@")
        sys.exit(-1)


     return args.input_grid, args.input_activities, args.property, args.output_grid, args.crop

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
           coords.append([float(line[1])/10,float(line[2])/10,float(line[3])/10])
    grid_file.close()
    return coords,num_points

def getActivities(input_activities):
    activities=[]
    Contact_i=[]
    Exposure_i=[]
    activities_file=open(input_activities,'r')
    for line in activities_file:
        if (line.startswith("Point")):
           continue
        activities.append(float(line.split()[1]))
        Contact_i.append(float(line.split()[2]))
        Exposure_i.append(float(line.split()[3]))
    activities_file.close()
    return activities,Contact_i,Exposure_i

def getGridObject(coords,activities,num_points):
    top = mdtraj.Topology()
    chain = top.add_chain()
    residue = top.add_residue("GRI", chain)
    #germanium = mdtraj.element.germanium
    hydrogen = mdtraj.element.hydrogen
    for i in range(0,num_points):
        #top.add_atom("GRI", germanium, residue)
       top.add_atom("GRI", hydrogen, residue) # JJJ modification to visualize the grid as non bonded atoms
    xyz = numpy.array(coords)
    grid_object = mdtraj.Trajectory(xyz, top)
    return grid_object
    
def outputPdb(grid_object,activities,output_grid):
    grid_object.save('temp.pdb')
    tempfile=open('temp.pdb','r')
    outfile=open(output_grid,'w')
    i = 1 # JJJ modification to visualiize the grid as non bonded atoms
    for line in tempfile:
        if line.startswith('ATOM'):
           pdb_index=int(line.split()[1])
           activities_index=pdb_index-1 # indices in pdb files start at 1
           if activities[activities_index]==0 and crop=='yes':
              continue
           activity_pdb='{0:.2f}'.format(activities[activities_index]) #Occuppancies in pdb files only have 2 decimals
           myline=line.replace('A   0','A%4d'%i)# JJJ modification to visualize the grid as non bonded atoms
           new_line=myline[0:56]+str(activity_pdb)+myline[60:]
           outfile.write(new_line)
           outfile.write("TER\n") # JJJ modification to visualize the grid as non bonded atoms
           i = i + 1 #JJJ modification to visualize the grid as non bonded atoms
    tempfile.close()
#    os.system('rm temp.pdb')
    outfile.close()
    return 0

if __name__ == '__main__':
   print (help_message) 
   input_grid,input_activities,prop,output_grid,crop=parse(parser)
   coords,num_points=getGridCoordinates(input_grid)
   activities,Contact_i,Exposure_i=getActivities(input_activities)
   grid_object=getGridObject(coords,activities,num_points)
   output_grid_base=output_grid.split('.pdb')[0]
   if (prop=="activity"):
      output_grid =output_grid_base+".activity.pdb"
      outputPdb(grid_object,activities,output_grid)
   elif(prop=="contact"):
      output_grid =output_grid_base+".contact.pdb"
      outputPdb(grid_object,Contact_i,output_grid)
   elif (prop=="exposure"):
      output_grid =output_grid_base+".exposure.pdb"
      outputPdb(grid_object,Exposure_i,output_grid)
   elif (prop=="all"):
      output_grid =output_grid_base+".activity.pdb"
      outputPdb(grid_object,activities,output_grid)
      output_grid =output_grid_base+".contact.pdb"
      outputPdb(grid_object,Contact_i,output_grid)
      output_grid =output_grid_base+".exposure.pdb"
      outputPdb(grid_object,Exposure_i,output_grid)

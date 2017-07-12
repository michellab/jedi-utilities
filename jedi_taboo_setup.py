import argparse, os,sys, copy, math
import mdtraj, numpy

parser = argparse.ArgumentParser(description="Analyse a GROMACS trajectory \
                                              find the snapshot with the best JEDI score \
                                              and generate a structure to bias away from",
                                 epilog="jedi_taboo_setup.py is distributed under the GPL.",
                                 prog="jedi-setup.py")
parser.add_argument('-g','--gro', nargs="?",
                    help='Name of input GRO file.')
parser.add_argument('-t','--trr', nargs="?",
                    help="Name of the input trr file")
parser.add_argument('-c','--colvar', nargs="?",
                    help="PLUMED COLVAR file containing time in ps and JEDI values. Default 'COLVAR'",
                    default='COLVAR')
parser.add_argument('-j','--jedi', nargs="?",
                    help="Name the JEDI CV has in your PLUMED input file. default 'je'",
                    default="je")
parser.add_argument('-s','--selection', nargs="?",
                    help="Atom set that will be used to calculate the RMSD bias. Values are 'CA', 'backbone' or 'site'")
parser.add_argument('-a','--apolar', nargs="?",
                    help="PDB file that contains the apolar atoms of the binding site. Only to be used with '-s site'",
                    default='apolar.pdb')
parser.add_argument('-p','--polar', nargs="?",
                    help="PDB file that contains the polar atoms of the binding site. Only to be used with '-s site'",
                    default='polar.pdb')
parser.add_argument('-o','--output', nargs="?",
                    help="PDB file that will be used as a reference for RMSD calculations.'",
                    default='ref.pdb')


def parse(parser):
    args = parser.parse_args()

    if args.gro is None:
        parser.print_help()
        print ("""@@@
ERROR ! You must specificy an input GRO. See usage above. 
@@@""")
        sys.exit(-1)

    if not os.path.isfile(args.gro):
        parser.print_help()
        print ("""@@@
ERROR ! The input GRO cannot be found. See usage above.
@@@""")
        sys.exit(-1)

    if args.trr is None:
        parser.print_help()
        print ("""@@@
ERROR ! You must specificy an input TRR. See usage above. 
@@@""")
        sys.exit(-1)

    if not os.path.isfile(args.trr):
        parser.print_help()
        print ("""@@@
ERROR ! The input TRR cannot be found. See usage above.
@@@""")
        sys.exit(-1)

    if not os.path.isfile(args.colvar):
        parser.print_help()
        print ("""@@@
ERROR ! The input COLVAR file cannot be found. See usage above.
@@@""")
        sys.exit(-1)

    if (args.selection != 'CA') and (args.selection != 'backbone')  and (args.selection != 'site'):
        print ("""@@@
ERROR ! The -s/--selection option only accepts "CA", "backbone" or "site".
@@@""")
        sys.exit(-1)
  
    return args.gro, args.trr, args.colvar, args.selection, args.jedi, args.apolar, args.polar, args.output

    



def findMaxJEDI(colvar,name):
    filein=open(colvar,'r')
    jedi=[]
    time=[]
    for line in filein:
        line=line.split()
        if line[0].startswith('#'):
           if name not in line:
              print "Name "+name+" not found in file "+colvar+". EXITING."
              sys.exit(-1)
           ind=line.index(name)-2
        else:
           time.append(int(float(line[0])))
           jedi.append(float(line[ind]))

    maxJEDI=max(jedi)
    maxJEDItime=time[jedi.index(maxJEDI)]

    return maxJEDI, maxJEDItime

def selectFrame(gro,trr,maxJEDItime):
    traj=mdtraj.load_trr(trr, top=gro)
    time=traj.time.tolist()
   # print time
   # sys.exit()
    framenumber=time.index(float(maxJEDItime))
    coords=traj.xyz[framenumber]
    top=traj.topology
    frame=mdtraj.Trajectory(coords,top)
#    frame.save_pdb('frame.pdb')
    print frame
    return frame, framenumber

def genRef(selection,apolar,polar,frame,outfile):
    reflist=[]
    if selection=='site':
       for pdb in [apolar,polar]:
           filein=open(pdb,'r')  
           for line in filein:
               if line.startswith('ATOM') or line.startswith('HETATM'):
                  reflist.append(int(line.split()[1])-1)
                  
           filein.close()
       reflist=sorted(reflist)

    elif selection=='CA':
       frame.save_pdb('frame.pdb')
       filein=open('frame.pdb','r')
       for line in filein:
           if line.startswith('ATOM') or line.startswith('HETATM'):
              line=line.split()
              if line[2]=='CA':
                 reflist.append(int(line[1])-1)
       filein.close()

    elif selection=='backbone':
       frame.save_pdb('frame.pdb')
       filein=open('frame.pdb','r')
       for line in filein:
           if line.startswith('ATOM') or line.startswith('HETATM'):
              line=line.split()
              if (line[2]=='CA') or (line[2]=='C') or (line[2]=='O') or (line[2]=='N'):
                 reflist.append(int(line[1])-1)
       filein.close()

    #This is to keep atom masses and indices because mdtraj doesn't seem to do it 
    ref=frame.atom_slice(reflist)
    ref.save_pdb('temp.pdb')
#    sys.exit()
    filein=open('temp.pdb','r')
    fileout=open(outfile,'w')
    idx=0
    for line in filein:
        if line.startswith('ATOM') or line.startswith('HETATM'):
           pdb_idx=reflist[idx]+1
           mass = frame.topology.atom(reflist[idx]).element.mass
           if pdb_idx < 10:
                str_idx = '    %s' % pdb_idx
           elif pdb_idx < 100:
                str_idx= '   %s' % pdb_idx
           elif pdb_idx < 1000:
                str_idx = '  %s' % pdb_idx
           elif pdb_idx < 10000:
                str_idx = ' %s' % pdb_idx
           else:
                str_idx = '%s' % pdb_idx
           
           str_mass = '%-5.2f' % mass
 
           line = line[0:6] + str_idx + line[11:55] + str_mass + line[60:61] + '1.00' + line[66:] # modifying b-factor so the rmsd is properly calculated
           idx=idx+1
        fileout.write(line)
#    os.system('rm frame.pdb temp.pdb')
                 

if __name__ == '__main__':
    
    gro, trr, colvar, selection, jedi, apolar, polar, output = parse(parser)
    
    maxJEDI, maxJEDItime = findMaxJEDI(colvar,jedi)
    print "max JEDI value is ", maxJEDI, " at ", maxJEDItime, "ps."
    
    frame, framenumber = selectFrame(gro,trr,maxJEDItime)
    print "you got max jedi at snapshot ",framenumber
    
    genRef(selection,apolar,polar,frame,output)

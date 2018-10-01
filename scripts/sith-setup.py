#!/usr/bin/env python
#
"""
 @AUTHOR Joan Clark-Nicolas. Oct 2017.
 @USAGE
  ///FILL IN WITH USAGE///

Structure Investigation by Taboo SearcH
Perform an iterative taboo search using CVs of choice.

optional arguments

///FILL IN WITH ARGUMENTS
  -h, --help            show this help message and exit

sith.py is distributed unfer a GPL license.

 @DEPENDENCIES:
 mdtraj, numpy, scipy, pymp

"""
import argparse
import mdtraj
import sys


parser = argparse.ArgumentParser(description="Perform an iterative taboo search using the JEDI collective variable",
                                 epilog="sith.py is distributed under a GPL license.",
                                 prog="sith.py")
parser.add_argument('-a','--apolar', nargs="?",
                help="name of output apolar atoms file. Default apolar.pdb",
                default="apolar.pdb")
parser.add_argument('-p','--polar', nargs="?",
                help="name of output polar atoms file. Default polar.pdb",
                default="polar.pdb")
parser.add_argument('-c','--cv', nargs="?",
                help="Type of CV to use with SITH. Default is COM",
                default="COM")
parser.add_argument('-o','--plumed', nargs="?",
                help="Output plumed file",
                default="COM")

def parse(parser):
    args = parser.parse_args()
    return args.apolar, args.polar, args.cv, args.plumed

def getAtoms(apolar,polar,cv):
    if  cv=="COM":
        comlist={}
        for pdb in [apolar,polar]:
            site=mdtraj.load_pdb(pdb)
            for atom in site.topology.atoms:
                if atom.residue not in comlist.keys():
                   comlist[atom.residue]=[atom.serial]
                else:
                   comlist[atom.residue].append(atom.serial)

            return comlist
    else:
        print "The selected CV is not implemented yet."
        sys.exit()

    
def genCV(cv,plumed,cvlist):
    fileout=open(plumed,'w')
    if cv=="COM":
       labels=[]
       for key in  cvlist.keys():
          line='COM_'+str(key)+": COM ATOMS="+','.join(str(num) for num in cvlist[key])+'\n'
          fileout.write(line)
          labels.append('COM_'+str(key))
       fileout.write('\n')
       distance_labels=[]
       for i in range(0,len(labels)):
           label1=labels[i]
           for j in range(i+1,len(labels)):
               label2=labels[j] 
               if label1==label2:
                  continue
               dislab='d_'+label1+'_'+label2
               distance_labels.append(dislab)
               line=dislab+': DISTANCE ATOMS='+label1+','+label2+'\n'
               fileout.write(line)

       fileout.write('\n')
               
       line='SITH ARG='+','.join(distance_labels)+\
            " SITHSTRIDE=50000 SITHFILE=clusters.dat DC=0.01 DELTA0=0.05 CVSTRIDE=2500 CVFILE=cvs.dat TYPOT=LOWER_WALLS\n\n"
       fileout.write(line)
       line="PRINT ARG="+','.join(distance_labels)+" FILE=COLVAR STRIDE=2500"
       fileout.write(line)
       
       fileout.close()
    
          

if __name__=="__main__":

     apolar,polar,cv,plumed=parse(parser)
     cvlist=getAtoms(apolar,polar,cv)
     genCV(cv,plumed,cvlist)
     print "SITH plumed input file has been generated. You might need to change some parameters"
     

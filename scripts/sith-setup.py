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
import numpy as np


parser = argparse.ArgumentParser(description="Perform an iterative taboo search using the JEDI collective variable",
                                 epilog="sith.py is distributed under a GPL license.",
                                 prog="sith.py")
parser.add_argument('-i','--input_pdb', nargs="?",
                help="Name of the pdb reference file, only needed for certain types of CVs",
                default=None)
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
                default="plumed.dat")

def parse(parser):
    args = parser.parse_args()
    return args.input_pdb, args.apolar, args.polar, args.cv, args.plumed

def getAtoms(apolar,polar,cv,refStruct):
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
    elif cv=="TORSION":
        torlist={}
        if refStruct is None:
            sys.exit("You need the reference structure if you use TORSION. Exiting.")
        structure=mdtraj.load_pdb(refStruct)
        atRes=[[],[]]
        for atom in structure.topology.atoms:
            atRes[0].append(int(atom.serial))
            atRes[1].append(str(atom.residue))
        #atRes=np.array(atRes)

        reslist=[]
        for pdb in [apolar,polar]:
            site=mdtraj.load_pdb(pdb)
            for res in site.topology.residues:
                if res not in reslist:
                   reslist.append(str(res))

        for i in range(0,len(mdtraj.compute_chi1(structure)[0])):
            atomNumber=mdtraj.compute_chi1(structure)[0][i][0]
            resPos=atRes[0].index(atomNumber+1)
            residue=atRes[1][resPos]
            if residue not in reslist: 
               continue
            label="chi1_"+residue
            atomlist=[]
            for atom in mdtraj.compute_chi1(structure)[0][i]:
                atomlist.append(str(atom))
            torlist[label]=atomlist

        for i in range(0,len(mdtraj.compute_chi2(structure)[0])):
            atomNumber=mdtraj.compute_chi2(structure)[0][i][0]
            resPos=atRes[0].index(atomNumber+1)
            residue=atRes[1][resPos]
            if residue not in reslist:
               continue
            label="chi2_"+residue
            atomlist=[]
            for atom in mdtraj.compute_chi2(structure)[0][i]:
                atomlist.append(str(atom))
            torlist[label]=atomlist

        for i in range(0,len(mdtraj.compute_chi3(structure)[0])):
            atomNumber=mdtraj.compute_chi3(structure)[0][i][0]
            resPos=atRes[0].index(atomNumber+1)
            residue=atRes[1][resPos]
            if residue not in reslist:
               continue
            label="chi3_"+residue
            atomlist=[]
            for atom in mdtraj.compute_chi3(structure)[0][i]:
                atomlist.append(str(atom))
            torlist[label]=atomlist

        for i in range(0,len(mdtraj.compute_chi4(structure)[0])):
            atomNumber=mdtraj.compute_chi4(structure)[0][i][0]
            resPos=atRes[0].index(atomNumber+1)
            residue=atRes[1][resPos]
            if residue not in reslist:
               continue
            label="chi4_"+residue
            atomlist=[]
            for atom in mdtraj.compute_chi4(structure)[0][i]:
                atomlist.append(str(atom))
            torlist[label]=atomlist
        
        return torlist
        #print torlist
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
       cv_labels=[]
       for i in range(0,len(labels)):
           label1=labels[i]
           for j in range(i+1,len(labels)):
               label2=labels[j] 
               if label1==label2:
                  continue
               dislab='d_'+label1+'_'+label2
               cv_labels.append(dislab)
               line=dislab+': DISTANCE ATOMS='+label1+','+label2+'\n'
               fileout.write(line)

       fileout.write('\n')

    elif cv=="TORSION":
         cv_labels=[]
         for key in  cvlist.keys():
             line=str(key)+": TORSION ATOMS="+','.join(str(num) for num in cvlist[key])+'\n'
             fileout.write(line)
             label_sin='sin_'+key
             line=label_sin+': MATHEVAL ARG='+key+' FUNC=sin(x) PERIODIC=NO'+'\n'
             fileout.write(line)
             label_cos='cos_'+key
             line=label_cos+': MATHEVAL ARG='+key+' FUNC=cos(x) PERIODIC=NO'+'\n'
             fileout.write(line)
             cv_labels.append(label_sin)
             cv_labels.append(label_cos)

    line='SITH ARG='+','.join(cv_labels)+\
         " SITHSTRIDE=50000 SITHFILE=clusters.dat DC=0.01 DELTA0=0.05 CVSTRIDE=2500 CVFILE=cvs.dat TYPOT=LOWER_WALLS\n\n"
    fileout.write(line)
    line="PRINT ARG="+','.join(cv_labels)+" FILE=COLVAR STRIDE=2500"
    fileout.write(line)
     
    fileout.close()
    
          

if __name__=="__main__":

     input_pdb,apolar,polar,cv,plumed=parse(parser)
     cvlist=getAtoms(apolar,polar,cv,input_pdb)
     genCV(cv,plumed,cvlist)
     print "SITH plumed input file has been generated. You might need to change some parameters"
     

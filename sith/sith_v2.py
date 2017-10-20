#!/usr/bin/env python
#
"""
 @AUTHOR Julien Michel. Sep 2015.
 @USAGE
  ///FILL IN WITH USAGE///

Perform an iterative taboo search using the JEDI collective variable

optional arguments

///FILL IN WITH ARGUMENTS
  -h, --help            show this help message and exit

sith.py is distributed unfer a GPL license.

 @DEPENDENCIES:
 mdtraj, numpy, scipy

"""

import argparse
import mdtraj
import numpy, scipy, math
import glob, os, sys

parser = argparse.ArgumentParser(description="Perform an iterative taboo search using the JEDI collective variable",
                                 epilog="sith.py is distributed under a GPL license.",
                                 prog="jedi-setup.py")
#
parser.add_argument('-i','--input', nargs="?",
                    help='input file')


# Get all the parameters from the input file
def parse(parser):

    args = parser.parse_args()
    if args.input is None:
        parser.print_help()
        print ("""@@@
ERROR ! You must specificy an input file. See usage above. 
@@@""")
        sys.exit(-1)
    elif not os.path.isfile(args.input):
        parser.print_help()
        print ("""@@@
ERROR ! The input cannot be found. See usage above.
@@@""")
        sys.exit(-1)
    else:
        parameters={}
        filein=open(args.input,'r')
        for line in filein:
            if line.startswith('#') or line=="\n":
                continue
            if '#' in line:
                line=line.split('#')[0] # allows you to put comments after the parameters
            if '=' not in line:
                print "ERROR: line '"+line+"' is not properly formatted. See usage."
                sys.exit(-1)

            line=line.replace(' ','').strip('\n').split('=')

            if len(line)!=2:
                print "ERROR: line '"+line+"' is not properly formatted. See usage."
                sys.exit()
            if line[0] in parameters.keys():
                print "ERROR: parameter '"+label[0]+"' is already used. Check input file. Exiting."
                sys.exit()
            
            parameters[line[0]]=line[1]

    return parameters


# Check that the chosen CV is supported and get the necessary parameters
def getCV(parameters):

    supported_cvs=['JEDI']

    if 'cv' not in parameters.keys():
       print "ERROR: Main CV must be defined with the option 'cv='."
       print "Supported CVs so far are: "+','.join(supported_cvs)+'. Exiting'
       sys.exit()
    elif parameters['cv'] not in supported_cvs:
       print "ERROR: CV '"+parameters['cv']+"' is not supported yet."
       print "Supported CVs so far are: "+','.join(supported_cvs)+'. Exiting'
       sys.exit()
    else:
       cv=parameters['cv']
    
    if cv=='JEDI':

       if 'apolar' not in parameters.keys():
          parameters['apolar']='apolar.pdb'
       if not os.path.isfile(parameters['apolar']):
          print "ERROR: File '"+parameters['apolar']+"' not found. Exiting."
          sys.exit()

       if 'polar' not in parameters.keys():
          parameters['polar']='polar.pdb'
       if not os.path.isfile(parameters['polar']):
          print "ERROR: File '"+parameters['polar']+"' not found. Exiting."
          sys.exit()

       if 'grid' not in parameters.keys():
          parameters['grid']='grid.pdb'
       if not os.path.isfile(parameters['grid']):
          print "ERROR: File '"+parameters['grid']+"' not found. Exiting."
          sys.exit()

       if 'ligand' not in parameters.keys():
          parameters['ligand']=None
       elif not os.path.sdfile(parameters['ligand']):
          print "ERROR: File '"+parameters['ligand']+"' not found. Exiting."
          sys.exit()

       if 'jedi_params' not in parameters.keys():
          parameters['jedi_params']='jedi.params'
       if not os.path.isfile(parameters['jedi_params']):
          print "ERROR: File '"+parameters['jedi_params']+"' not found. Exiting."
          sys.exit()

       if 'stride' not in parameters.keys():
          print "ERROR: Plumed CV calculation stride must be defined with option 'stride='. Exiting"
          sys.exit()

       if 'summary' not in parameters.keys():
          parameters['summary']='jedi_stats.dat'
 
       if 'gridstride' not in parameters.keys():
          parameters['gridstride']='0'

       if 'metaD_sigma' not in parameters.keys():
          parameters['metaD_sigma']=None

       if 'dumpderivatives' not in parameters.keys():
          parameters['dumpderivatives']='0'

       return parameters   


# Check that the metric is supported and it works wit the CV.
def getMetric(parameters):

    supported_metrics=['SC_TORSION']
    supported_cvs_metric={'SC_TORSION':['JEDI']}
    
    if 'metric' not in parameters.keys():
       print "ERROR: The metric that drives the taboo search must be defined with the option 'metric='. Exiting."
       print "Supported metrics so far are: "+','.join(supported_metrics)+'. Exiting'
       sys.exit()
    elif parameters['metric'] not in supported_metrics:
       print "ERROR: Metric '"+parameters['metric']+"' is not supported yet."
       print "Supported metrics so far are: "+','.join(supported_metrics)+'. Exiting'
       sys.exit()
    else:
       metric=parameters['metric']
       cv=parameters['cv']
       if cv not in supported_cvs_metric[metric]:
          print "ERROR: Collective variable '"+cv+"' not supported with metric '"+metric+"'."
          print "This metric supports the following CVs: "+','.join(supported_cvs_metric[metric])
          sys.exit()
    
def setupMetric(parameters):
    
    metric=parameters['metric']
    cv=parameters['cv']
    
    if metric=='SC_TORSION': #FIXME: need to specify that it has to be pdb and raise an error if it's not
       if 'reference_structure' not in parameters.keys():
          print "ERROR: you need to provide a reference structure to measure side chain torsions." 
          print "Use option 'reference_structure=' to provide a file. Exiting."
          sys.exit()
       elif not os.path.isfile(parameters['reference_structure']):
          print "ERROR: file '"+parameters['reference_structure']+"' not found. Exiting."
          sys.exit()
       else:
          structure=parameters['reference_structure']

       if cv=='JEDI':
          apolar=parameters['apolar']
          polar=parameters['polar']

          residues_str=[]
          residues=[]
          for pdb in [apolar,polar]:
              structure=mdtraj.load_pdb(pdb)
              for res in structure.topology.residues:
                  if str(res) not in residues_str and "GLY" not in str(res):
                     residues_str.append(str(res))
                     residues.append(res)
          print "There are", len(residues_str), "residues in the defined binding site (which could be the whole protein!)"
          
          struct=mdtraj.load_pdb(structure)
          
          

          chi={}
          
          for i in range(0,len(mdtraj.compute_chi1(struct)[0])):
              atoms=[]
              for atom in mdtraj.compute_chi1(struct)[0][i]:
                  for at2 in struct.topology.atoms:
                      if int(at2.index)==atom and str(at2.residue) in residues_str:
                         atoms.append(str(atom+1)) # +1 because they start at 0 in mdtraj but at 1 in GROMACS
                         name='chi1_'+str(at2.residue)
              if len(atoms)==4:
                 chi[name]=atoms
              elif len(atoms)!=0:
                 print "something went wrong when assigning dihedrals. A dihedral can't have a number of atoms different than 4"
                 break

          for i in range(0,len(mdtraj.compute_chi2(struct)[0])):
              atoms=[]
              for atom in mdtraj.compute_chi2(struct)[0][i]:
                  for at2 in struct.topology.atoms:
                      if int(at2.index)==atom and str(at2.residue) in residues_str:
                         atoms.append(str(atom+1)) # +1 because they start at 0 in mdtraj but at 1 in GROMACS
                         name='chi2_'+str(at2.residue)
              if len(atoms)==4:
                 chi[name]=atoms
              elif len(atoms)!=0:
                 print "something went wrong when assigning dihedrals. A dihedral can't have a number of atoms different than 4"
                 break

          for i in range(0,len(mdtraj.compute_chi3(struct)[0])):
              atoms=[]
              for atom in mdtraj.compute_chi3(struct)[0][i]:
                  for at2 in struct.topology.atoms:
                      if int(at2.index)==atom and str(at2.residue) in residues_str:
                         atoms.append(str(atom+1)) # +1 because they start at 0 in mdtraj but at 1 in GROMACS
                         name='chi3_'+str(at2.residue)
              if len(atoms)==4:
                 chi[name]=atoms
              elif len(atoms)!=0:
                 print "something went wrong when assigning dihedrals. A dihedral can't have a number of atoms different than 4"
                 break

          for i in range(0,len(mdtraj.compute_chi4(struct)[0])):
              atoms=[]
              for atom in mdtraj.compute_chi4(struct)[0][i]:
                  for at2 in struct.topology.atoms:
                      if int(at2.index)==atom and str(at2.residue) in residues_str:
                         atoms.append(str(atom+1)) # +1 because they start at 0 in mdtraj but at 1 in GROMACS
                         name='chi4_'+str(at2.residue)
              if len(atoms)==4:
                 chi[name]=atoms
              elif len(atoms)!=0:
                 print "something went wrong when assigning dihedrals. A dihedral can't have a number of atoms different than 4"
                 break

          fileout=open('torsions_sidechains.dat','w')
          for key in chi.keys():
              line=key+': TORSION ATOMS='+','.join(chi[key])+'\n'
              fileout.write(line)
              line='sin_'+key+': MATHEVAL ARG='+key+' FUNC=sin(x) PERIODIC=NO'+'\n'
              fileout.write(line)
              line='cos_'+key+': MATHEVAL ARG='+key+' FUNC=cos(x) PERIODIC=NO'+'\n'
              fileout.write(line)
          fileout.close()
          print "The sines and cosines of "+str(len(chi.keys()))+ " torsions are going to be used as a metric."
          print "This is a total of "+str(3*len(chi.keys()))+" variables"
          return chi

          

 
if __name__ == '__main__':

    ######## SET UP PARAMETERS AND RUN SANITY CHECKS #########
    
    parameters=parse(parser)
    
    parameters=getCV(parameters) #Set up the CV of use

    getMetric(parameters)


   ######### DO STUFF ######################
   
    metric_input=setupMetric(parameters)

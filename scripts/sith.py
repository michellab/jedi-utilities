#!/usr/bin/env python
#
"""
 @AUTHOR Joan Clark-Nicolas. Oct 2017.
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
from scipy import spatial
import glob, os, sys
import time as tiempo
import psutil

parser = argparse.ArgumentParser(description="Perform an iterative taboo search using the JEDI collective variable",
                                 epilog="sith.py is distributed under a GPL license.",
                                 prog="jedi-setup.py")
#
parser.add_argument('-i','--input', nargs="?",
                    help='input file')
parser.add_argument('-d','--debug',action='store_true',
                    help='Do not remove any files (rename them if necessary)')

# Funtions that make your life easier

def isfloat(value):
  try:
    float(value)
    return True
  except ValueError:
    return False

# Get all the parameters from the input file
def parse(parser):

    args = parser.parse_args()
    print args.debug
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
        if args.debug is True:
           parameters['debug']=True
        else:
           parameters['debug']=False
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

    # Getting Target extension           
    supported_target_extensions=['pdb','gro','trr','dcd','crd','crdbox','g96','trr','trj','xtc']
    if 'target' not in parameters.keys():
        parameters['target']=None
    else:
        target_ext=parameters['target'].split('.')[-1]
        if target_ext not in supported_target_extensions:
           print "ERROR: Target extension '"+target_ext+"' not supported."
           print "Supported target extensions are: "+','.join(supported_target_extensions)+". Exiting."
           sys.exit()
        else:
           parameters['target_extension']=target_ext

    # Getting MD engine executables
    supported_mden=['GROMACS']
    if 'md_engine' not in parameters.keys():
       parameters['md_engine']='GROMACS'
       print "MD Engine not specified. Going to use "+parameters['md_engine']+'.'
    elif parameters['md_engine'] not in supported_mden:
       print "ERROR: MD engine '"+parameters['md_engine']+"is not supported (yet?)."
       print "Supported MD engines are: "+','.join(supported_mden)
       sys.exit()
    if 'simtime' not in parameters.keys():
        print "Please specify the simulation time using option 'simtime='. Exiting"
        sys.exit()
    if 'dt' not in parameters.keys():
        print "Please specify MD time step using option 'dt=. Exiting'"
        sys.exit

    nsteps=str(int(float(parameters['simtime'])/float(parameters['dt'])))
    parameters['nsteps']=nsteps

    if parameters['md_engine']=='GROMACS':
       if 'gmx_path' not in parameters.keys():
          parameters['gmx_path']='/usr/bin/gmx'
          print "looking for gromacs in /usr/bin/gmx"
       if not os.path.isfile(parameters['gmx_path']):
          print "ERROR: GROMACS executable cannot be found. Exiting"
          sys.exit()
       # Check that we have all the files that gromacs needs
       if 'mdp' not in parameters.keys():
          print "Gromacs mdp file was not specified. Will use 'md.mdp'"
          parameters['mdp']='md.mdp'
       elif not os.path.isfile(parameters['mdp']):
          print "Error: file '"+parameters['mdp']+"'. Could not be found. Exiting."
          sys.exit()
       
       if 'tpr' not in parameters.keys():
          print "You must supply the tpr file for the fist iteration with option 'tpr='. Exiting."
          sys.exit()
       elif not os.path.isfile(parameters['tpr']):
          print "file "+parameters['tpr']+" not found. Exiting"
          sys.exit()

       if 'top' not in parameters.keys():
          print "You must supply the top file for grompp with option 'top='. Exiting."
          sys.exit()
       elif not os.path.isfile(parameters['top']):
          print "file "+parameters['top']+" not found. Exiting"
          sys.exit()
       else:
          filein=open(parameters['top'],'r')
          for line in filein:
              if line.startswith(';'):
                 continue
              elif "[ moleculetype ]" in line:
                 break
              else:
                 if "#include" in line and ".ff/" not in line:
                    line=line.split()
                    file_to_search=line[1].strip("\"")
                    if not os.path.isfile(file_to_search):
                       print "Error: file '"+file_to_search+"' was not found but is included in the GROMACS top file. Exiting."
                       sys.exit()

       if 'ndx' not in parameters.keys():
          print "You must supply the index file for grompp with 'ndx='. Exiting."
          sys.exit()
       elif not os.path.isfile(parameters['ndx']):
          print "file "+parameters['tpr']+" not found. Exiting"
          sys.exit()



    # Getting PLUMED executable
    if 'plumed_exec' not in parameters.keys():
        parameters['plumed_exec']='/usr/bin/plumed'
        print "looking for PLUMED in /usr/bin/plumed"
    if not os.path.isfile(parameters['plumed_exec']):
        print "WARNING! PLUMED executable not found." 
        print "You might be able to continue if you didn't specify a target AND you are sure that your MD code is patched."

    
    # MPI related stuff
    if 'mpi_exec' not in parameters.keys():
        print "looking for mpirun in /usr/bin/mpirun." 
        print "WARNING: THIS PROGRAM WILL CRASH IF YOU DON'T HAVE MPI."
        parameters['mpi_exec']='/usr/bin/mpirun'
 
    # Getting info to restart in case of crash:
    if 'st_iter' not in parameters.keys():
        parameters['st_iter']=0
    elif int(parameters['st_iter'])!=0:
        print "WARNING. Using restarting function still to be debugged!!!!!"
        #sys.exit()
    if 'clusters_hist' not in parameters.keys(): #history of sampled clusters
        parameters['clusters_hist']={}
    else:
        print "Restarting function has still to be debugged. Do NOT use. Exiting"
        sys.exit()
        parameters['clusters_hist']=parameters['clusters_hist'].split(',') 

    if 'include' not in parameters.keys():
        parameters['include']=[]
    else:
        print "Restarting function has still to be debugged. Do NOT use. Exiting"
        sys.exit()
        parameters['include']=parameters['include'].split(',')
        for bias in parameters['include']:
            if not os.path.isfile(bias):
               print "Biasing file '"+bias+"' needed and not found. Exiting"
               sys.exit

    #Getting clustering info    
    supported_clusterings=['density_Laio']
    if 'clustering' not in parameters.keys():
        print "Clustering method not specified. Will use the density method by Rodriguez and Laio (Science (2014), 344(6191) 1492-1496)"
        parameters['clustering']='density_Laio'
    elif parameters['clustering'] not in supported_clusterings:
        print "Clustering method must be specified with option 'clustering='."
        print "Supported clustering methods are: "+','.join(supported_clusterings)+". Exiting."
        sys.exit()
    if parameters['clustering']=='density_Laio':
       print "Clustering method selected is the density method by Rodriguez and Laio (Science (2014), 344(6191) 1492-1496)"
       d0_methods=['avg-2sd','avg-sd','avg']
       if 'cluster_dc' not in parameters.keys():
          print "The value for dc or a way to calculate it has not been specified. Will use d0=euclidMat_avg-2*euclidMat_sd"
          parameters['cluster_dc']='avg-2sd'
       elif isfloat(parameters['cluster_dc'])==True:
          print "The value ", parameters['cluster_dc'], " is going to be used as dc for the clustering"
          parameters['dc_value']=float(parameters['cluster_dc'])
          parameters['cluster_dc']='value'
       elif parameters['cluster_dc'] not in d0_methods:
          print "The method to calculate dc is not supported."
          print "Currently supported methods are: "+','.join(d0_methods)+". exiting"
          sys.exit()
    
    if 'sampltime' not in parameters.keys():
        print "The time in ps for which a populated cluster has to be sampled must be defined with option 'sampltime='. Exiting."
        sys.exit()
    
    return parameters


def setup_queue(parameters):
    supported_queues=['slurm']
    
    if 'q_system' not in parameters.keys():
       parameters['q_system']='slurm'
    elif parameters['q_system'] not in supported_queues:
        print "ERROR: Queue system '"+parameters['q_system']+"' is not supported."
        print "Supported queue systems are: "+','.join(supported_queues)
        sys.exit()
    
    if 'q_name' not in parameters.keys():
       print "Error: The queue name must be specified with option 'q_name='. Exiting."
       sys.exit()
    
    if 'q_time' not in parameters.keys():
        print "Maximum time for the job was not specified. Setting it to 24 hours."
        parameters['q_time']="24:00:00"
    
    if parameters['q_system']=='slurm':
       print "Setting up que jobs for slurm"
       if 'sbatch' not in parameters.keys(): 
          parameters['sbatch']='/usr/local/bin/sbatch'
          print "looking for sbatch in /usr/local/bin/sbatch"
       if not os.path.isfile(parameters['sbatch']):
          print "Error: "+parameters['sbatch']+" not found. Exiting." 
          sys.exit()

    if parameters['q_system']=='archer':
       print "Setting up que jobs for ARCHER."
       print "If specified, the variable nthreads will be overriden woith 24 times \
              the number of nodes to be used"
       if 'num_nodes' not in parameters.keys():
          print "Number of nodes was not selected. Setting it to 1"
          parameters['num_nodes']=1
       num_nodes=int(parameters['num_nodes'])
       nthreads=num_nodes*24
       parameters['nthreads']=nthreads
       if 'budget' not in parameters.keys():
           print "Budget code needs to be specified with option 'budget='. Exiting."
           sys.exit()
       
       if 'qsub' not in parameters.keys():
          print "looking for qsub in /opt/pbs/12.2.401.141761/bin/qsub"
          parameters['qsub']='/opt/pbs/12.2.401.141761/bin/qsub'
          
           

       

    return parameters

def genEnVar(parameters):

    if 'max_iter' not in parameters.keys():
        print "Maximum number of iterations not specified. Setting it to 1000"    
        parameters['max_iter']=1000

    if 'nthreads' not in parameters.keys():
        print "The number of available cores was not specified. Setting it to 1."
        parameters['nthreads']='1'
    if 'ngpu' not in parameters.keys():
        print "The number of available gpus has to be specified with option 'ngpu='. Exiting"
        sys.exit()
    if 'ntomp' not in parameters.keys():
        print "The number of cores to be used in each replica was not specified. Setting it to 1."
        parameters['ntomp']='1'
    

    nthreads=parameters['nthreads']
    ntomp=parameters['ntomp']
    ngpu=parameters['ngpu']
    nsim=int(int(nthreads)/int(ntomp))
    parameters['nsim']=nsim
    if int(nthreads)<int(ntomp) or int(nthreads)%int(ntomp)!=0:
       print "The number of available cores has to be equal or bigger than and a multiple of the number"
       print "of cores to be used in each replica. Exiting"
       sys.exit()
    if int(ngpu)!=0:
       if nsim%int(ngpu)!=0:
          print "The number of simulations has to be a multiple of the number of GPUS for MPI related reasons." 
          print "Please fix that by editing ntomp, nthreads or ngpu [nsim=nthreads/ntomp]. Exiting"
          sys.exit()
    print "Each iteration will run "+str(nsim)+" replicas with "+str(ntomp)+" cores each"
    
    return parameters

# Check that the chosen CV is supported and get the necessary parameters
def getCV(parameters):
    
    supported_cvs=['JEDI']

    if 'cv' not in parameters.keys(): #FIXME: this assumes that only one CV is taken as the main one
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
    return parameters

def getBias(prameters):
    supported_biases=['RESTRAINT','LOWER_WALLS','UPPER_WALLS','MOVINGRESTRAINT_L']
    if 'bias' not in parameters.keys():
       print "CV biasing method was not chosen. It is not going to be biased."
       parameters['bias']=None
    elif parameters['bias'] not in supported_biases:
       print "ERROR: biasing method '"+parameters['bias']+"' is not supported."
       print "Supported biasing methods are: "+','.join(supported_biases)+". Exiting."
       sys.exit()
    bias=parameters['bias']
    
    if 'mts_cv' not in parameters.keys():
       print "Stride for multiple time step in biasing force calculation for CV not specified. Setting it to 1."
       parameters['mts_cv']=1
    if 'mts_metric' not in parameters.keys():
       print "Stride for multiple time step in biasing force calculation in metric not specified. Setting it to 1."
       parameters['mts_metric']=1

    if bias=='RESTRAINT' or  bias=='LOWER_WALLS' or bias=='UPPER_WALLS' or bias=='MOVINGRESTRAINT_L':
       if 'at_cv' not in parameters.keys():
           print "ERROR: You must define an equilibrium value for your CV with option 'at_cv='. Exiting"
       if 'kappa_cv' not in parameters.keys():
           print "KAPPA for CV not specified. setting it to 100 kJ/(mol*(nm**2))"
           parameters['kappa_cv']=100
       if 'at_metric' not in parameters.keys():
           print "ERROR: You must define an equilibrium value for your metric with option 'at_metric='. Exiting"
       if 'kappa_metric' not in parameters.keys():
           print "KAPPA for metric not specified. setting it to 100 kJ/(mol*(nm**2))"
           parameters['kappa_metric']=100
    
    return parameters
    

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
 
       torsion_specs=['SC_list','TOR_list']
       for spec in torsion_specs:
           atLeastOne=False
           if spec not in parameters.keys():
              parameters[spec]=None
           else:
              atLeastOne=True
       if atLeastOne==False:
          print "You need to specify how the torsions will be chosen (whole SC (SC_list) or just some TOR (TOR_list)). eXITING."
          sys.exit()
       if parameters['SC_list'] is not None:
          residues_str=parameters['SC_list'].split(',')
       elif parameters['TOR_list'] is not None:
          torsions={}
          residues_str=[]
          for residue in parameters['TOR_list'].split(','):
              res=residue.split(':')[0]
              residues_str.append(res)
              torsions[res]=[]
              for torsion in residue.split(':')[1].split(';'):
                  torsions[res].append(torsion)
          print residues_str
          print torsions
          #sys.exit()
       elif cv=='JEDI':
          apolar=parameters['apolar']
          polar=parameters['polar']

          residues_str=[]
          residues=[]
          for pdb in [apolar,polar]:
              site=mdtraj.load_pdb(pdb)
              for res in site.topology.residues:
                  if str(res) not in residues_str and "GLY" not in str(res):
                     residues_str.append(str(res))
                     residues.append(res)

       print "The residues that are going to be taken into account are:"
       print residues_str
       print "There are", len(residues_str), "residues in the defined binding site (which could be the whole protein!)"
          
       struct=mdtraj.load_pdb(structure) 
          
       chi={}
          
       for i in range(0,len(mdtraj.compute_chi1(struct)[0])):
           atoms=[]
           for atom in mdtraj.compute_chi1(struct)[0][i]:
               for at2 in struct.topology.atoms:
                   if int(at2.index)==atom and str(at2.residue) in residues_str:
                      if parameters['TOR_list']==None or (str(at2.residue) in torsions.keys() and 'chi1' in torsions[str(at2.residue)]):
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
                      if parameters['TOR_list']==None or (str(at2.residue) in torsions.keys() and 'chi2' in torsions[str(at2.residue)]):
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
                      if parameters['TOR_list']==None or (str(at2.residue) in torsions.keys() and 'chi3' in torsions[str(at2.residue)]):
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
                      if parameters['TOR_list']==None or (str(at2.residue) in torsions.keys() and 'chi4' in torsions[str(at2.residue)]):
                         atoms.append(str(atom+1)) # +1 because they start at 0 in mdtraj but at 1 in GROMACS
                         name='chi4_'+str(at2.residue)
           if len(atoms)==4:
              chi[name]=atoms
           elif len(atoms)!=0:
              print "something went wrong when assigning dihedrals. A dihedral can't have a number of atoms different than 4"
              break

       fileout=open('metric.dat','w')
       for key in chi.keys():
           line=key+': TORSION ATOMS='+','.join(chi[key])+'\n'
           fileout.write(line)
           line='sin_'+key+': MATHEVAL ARG='+key+' FUNC=sin(x) PERIODIC=NO'+'\n'
           fileout.write(line)
           line='cos_'+key+': MATHEVAL ARG='+key+' FUNC=cos(x) PERIODIC=NO'+'\n'
           fileout.write(line)
       fileout.close()
       print "The sines and cosines of "+str(len(chi.keys()))+ " torsions are going to be used as a metric."
       print "This is a total of "+str(2*len(chi.keys()))+" variables"
       return chi

          
def analyse_target(parameters,metric_input):
    
    metric=parameters['metric']
    cv=parameters['cv']
    target=parameters['target']
    
    fileout=open('metric_reference.dat','w')
    line='INCLUDE FILE=metric.dat\n'
    fileout.write(line)
    
    if cv=='JEDI':
       line='cv: JEDI'
       line=line+' APOLAR='+parameters['apolar']
       line=line+' POLAR='+parameters['polar']
       line=line+' GRID='+parameters['grid']
       if parameters['ligand'] is not None:
          line=line+' SITE='+parameters['ligand']
       line=line+' PARAMETERS='+parameters['jedi_params']
       line=line+' STRIDE='+parameters['stride']
       line=line+' SUMMARY='+parameters['summary']
       line=line+' GRIDSTRIDE='+parameters['gridstride']
       line=line+' DUMPDERIVATIVES='+parameters['dumpderivatives']
       if parameters['metaD_sigma'] is not None:
          line=line+' SIGMA='+parameters['metaD_sigma']
       line=line+'\n'
       fileout.write(line)

    if metric=='SC_TORSION':
       metric_lab=[]
       for key in metric_input.keys():
           sin='sin_'+key
           metric_lab.append(sin)
           cos='cos_'+key
           metric_lab.append(cos)

    line='PRINT ARG=cv,'+','.join(metric_lab)+' STRIDE=1 FILE=TARGET.colvar' # stride is 1 because we'll use driver
    fileout.write(line)

    fileout.close()
           
    if parameters['target_extension']=='xyz' or parameters['target_extension']=='gro':
       inputflag='--i'
    else:
       inputflag='--mf_'
    driver=parameters['plumed_exec']+' driver '+inputflag+parameters['target_extension']+' '+parameters['target']+' --plumed metric_reference.dat'
    os.system(driver)

    time,metric_list_val,cv_list_val,metric_avg,metric_sd,cv_avg,cv_sd=analyse_plumed_output(parameters,'TARGET.colvar',-1)

    print "Target CV is: "+str(cv_avg)+" +/- "+str(cv_sd)
    print "Target metric is: "+str(metric_avg)+" +/- "+str(metric_sd)

    return metric_avg,metric_sd,cv_avg,cv_sd

def analyse_plumed_output(parameters, nameIn,iteration):
     
    #FIXME: this is assuming that the main CV is only one AND it is the first one after 'time'. Should not be very difficult to modify.
    filein=open(nameIn,'r')
    cv_list_val=None
    metric_list_val=None
    time=numpy.array([])
    dt=float(parameters['dt'])
    stride=int(parameters['stride'])
    strideps=dt*stride
    lenline=None
    linum=0
    skip_first=False
    for line in filein:
        if line.startswith('#'):
           continue
        else:
           line=line.split()
           if lenline==None:
              lenline=len(line)
           else:
              if len(line)!=lenline:
                 message="File '"+nameIn+"' appears to be malformed at time "+str(time)+" ps, probably due to a crash. Skipping line"
                 continue
           if iteration==0 and float(line[0])==0 and skip_first==False:
              skip_first=True
           elif iteration==0 and float(line[0])==0 and skip_first==True:
              continue
           elif iteration>0 and float(line[0])==0:
              continue
           time_ps=linum*strideps
           time=numpy.append(time,time_ps)
           cv_lst_line=[]
           cv_lst_line.append(float(line[1]))
           if cv_list_val is None:
              cv_list_val=numpy.array(cv_lst_line)
           else:
              cv_list_val=numpy.append(cv_list_val,cv_lst_line)
           val_lst_line=numpy.array([])
           for val in line[2:]:
               val_lst_line=numpy.append(val_lst_line,float(val))
           if metric_list_val is None:
              metric_list_val=numpy.array([val_lst_line])
           else:
              metric_list_val=numpy.append(metric_list_val,[val_lst_line],axis=0)
        linum=linum+1

    #time=numpy.asarray(time)
    #cv_list_val=numpy.asarray(cv_list_val)
    #metric_list_val=numpy.asarray(metric_list_val)
    #print  type(metric_list_val)

    cv_avg=numpy.average(cv_list_val,axis=0)
    cv_sd=numpy.std(cv_list_val,axis=0)
 
    metric_avg=numpy.average(metric_list_val,axis=0)
    metric_sd=numpy.std(metric_list_val,axis=0)      

    return  time, metric_list_val, cv_list_val, metric_avg, metric_sd, cv_avg, cv_sd

def gen_plumed_input(parameters,include,clusters):

    metric=parameters['metric']
    cv=parameters['cv']
    bias=parameters['bias']
    stride=parameters['stride']

    fileout=open('taboo_bias.dat','w')
    
    fileout.write('INCLUDE FILE=metric.dat\n')    

    if cv=="JEDI":
       line='cv: JEDI'
       line=line+' APOLAR='+parameters['apolar']
       line=line+' POLAR='+parameters['polar']
       line=line+' GRID='+parameters['grid']
       if parameters['ligand'] is not None:
          line=line+' SITE='+parameters['ligand']
       line=line+' PARAMETERS='+parameters['jedi_params']
       line=line+' STRIDE='+parameters['stride']
       line=line+' SUMMARY='+parameters['summary']
       line=line+' GRIDSTRIDE='+parameters['gridstride']
       line=line+' DUMPDERIVATIVES='+parameters['dumpderivatives']
       if parameters['metaD_sigma'] is not None:
          line=line+' SIGMA='+parameters['metaD_sigma']
       line=line+'\n'
       fileout.write(line)
    
    if parameters['debug']==False and bias is not None:
       if bias=="LOWER_WALLS" or bias=="UPPER_WALLS" or bias=="RESTRAINT" or bias=='MOVINGRESTRAINT_L':
          at_cv=parameters['at_cv']
          kappa_cv=parameters['kappa_cv']
          mts_cv=parameters['mts_cv']
          at_metric=parameters['at_metric']
          kappa_metric=parameters['kappa_metric']
          mts_metric=parameters['mts_metric']
   
          line='res_cv: LOWER_WALLS ARG=cv AT='+str(at_cv)+' KAPPA='+str(kappa_cv)+' STRIDE='+str(mts_cv)+'\n'
          fileout.write(line)

    for center in clusters.keys():
        time_str=str(int(center))
        line='INCLUDE FILE=dist_metric_'+time_str+'.dat\n'
        fileout.write(line)


    for bias in include:
        line='INCLUDE FILE='+bias+'\n'
        fileout.write(line)
        
 
    if metric=='SC_TORSION':
       metric_lab=[]
       for key in metric_input.keys():
           sin='sin_'+key
           metric_lab.append(sin)
           cos='cos_'+key
           metric_lab.append(cos)
    else:
       metric_lab=metric_input.keys()

    line='PRINT ARG=cv,'+','.join(metric_lab)+' STRIDE='+stride+' FILE=COLVAR'
    fileout.write(line)
    fileout.close()

def build_metric_bias(parameters,clusters,metric_arr,iteration,metric_input,time):

    include=[]
    # If there are no populated clusters, there is no metric bias to be built.
    if len(clusters.keys())==0:
       return include
    else:
       clusdist_avg=parameters['at_metric']
       dt=float(parameters['dt'])
       stride=int(parameters['stride'])
       strideps=dt*stride
       sampltime=float(parameters['sampltime'])
       kappa_metric=float(parameters['kappa_metric'])
       kappas={}
       for center in clusters.keys():
           kappas[center]=(len(clusters[center])*strideps/sampltime)*kappa_metric
           time_str=str(int(center))
           nameOut="dist_metric_"+time_str+".dat"
           fileout=open(nameOut,'w')
           time_index=numpy.where(time==center)[0][0]
           labelstemp=[]
           loc=0
           for key in metric_input.keys():
               line='sin_'+key+'_'+time_str+': MATHEVAL ARG='+key+' FUNC=sin(x)-'+\
                     metric_arr[time_index][loc].astype(str)+' PERIODIC=NO'+'\n'
               #print line
               fileout.write(line)
               loc += 1
               labelstemp.append('sin_'+key+'_'+time_str)
               line='cos_'+key+'_'+time_str+': MATHEVAL ARG='+key+' FUNC=cos(x)-'+\
                     metric_arr[time_index][loc].astype(str)+' PERIODIC=NO'+'\n'
               fileout.write(line)
               #print line
               loc += 1
               labelstemp.append('cos_'+key+'_'+time_str)
           line='dist2_'+time_str+': COMBINE ARG='+','.join(labelstemp)+\
               ' POWERS='+','.join(['2']*len(metric_arr[time_index]))+' PERIODIC=NO\n'
           fileout.write(line)
           line='dist_'+time_str+': MATHEVAL ARG=dist2_'+time_str+' FUNC=sqrt(x) PERIODIC=NO'
           fileout.write(line)
           fileout.close()
      
    nameOut='metric_bias.'+str(iteration)+'.dat'
    fileout=open(nameOut,'w')
    for center in clusters.keys():
        time_str=str(int(center))
        if parameters['bias']=='LOWER_WALLS' or parameters['bias']=='UPPER_WALLS' or parameters['bias']=='RESTRAINT':
           line= line='rest_'+time_str+'_'+str(iteration)+': '+parameters['bias']+' ARG=dist_'+time_str+\
                      ' AT='+str(clusdist_avg)+' KAPPA='+str(kappas[center])+'\n'
        elif parameters['bias']=='MOVINGRESTRAINT_L':
           first=int(int(parameters['nsteps'])*0.2)
           line='MOVINGRESTRAINT ...\n'+\
                'ARG=dist_'+time_str+'\n'+\
                'VERSE=L\n'+\
                'STEP0=0 AT0=0 KAPPA0='+str(parameters['kappa_metric'])+'\n'+\
                'STEP1='+str(parameters['nsteps'])+' AT1='+str(clusdist_avg)+' KAPPA1='+str(kappas[center])+'\n'+\
                '... MOVINGRESTRAINT\n'
        fileout.write(line)
    fileout.close()
    include.append(nameOut)

    return include
    

def submit_calc(parameters,iteration):
    nameOut='iteration'+str(iteration)+'.sh'
    fileout=open(nameOut,'w')
    
    q_system=parameters['q_system']
    q_name=parameters['q_name']
  
    # Writing input file

    if q_system=='slurm':
      line='#!/bin/bash\n'
      line=line+'#SBATCH --job-name=iter'+str(iteration)+'\n'
      line=line+'#SBATCH -o iter'+str(iteration)+'.out\n'
      line=line+'#SBATCH -e iter'+str(iteration)+'.err\n'
      line=line+'#SBATCH -p '+q_name+'\n'
      line=line+'#SBATCH -n '+str(parameters['nthreads'])+'\n'
      line=line+'#SBATCH -N 1\n' # request entire node
      if int(parameters['ngpu'])>0:
         line=line+'#SBATCH --gres=gpu:'+str(parameters['ngpu'])+'\n'
      line=line+'#SBATCH --time '+parameters['q_time']+'\n'
      fileout.write(line)
    elif q_system=='archer':
      line='#!/bin/bash --login\n'
      line=line+'#PBS -l select='+str(parameters['num_nodes'])+'\n'
      line=line+'#PBS -N iter'+str(iteration)+'\n'
      line=line+'#PBS -A '+parameters['budget']+'\n'
      line=line+'#PBS -l walltime='+parameters['q_time']+'\n'
      line=line+'#PBS -q '+q_name+'\n\n'
      line=line+'export OMP_NUM_THREADS='+str(parameters['ntomp'])+'\n'
      line=line+'cd $PBS_O_WORKDIR\n'

    if parameters['md_engine']=='GROMACS':
       if iteration==0:
          for i in range(0,int(parameters['nsim'])):
             cmd='cp '+parameters['tpr']+' iteration0'+str(i)+'.tpr'
             os.system(cmd)
       deffnm='iteration'+str(iteration)
       line=parameters['mpi_exec']+' -n '+str(parameters['nsim'])+' '+\
            parameters['gmx_path']+' mdrun -v -ntomp '+str(parameters['ntomp'])+\
            ' -deffnm '+deffnm+' -nsteps '+parameters['nsteps']+' -plumed taboo_bias.dat -multi '+str(parameters['nsim'])+'\n'
       fileout.write(line)

       gettime=tiempo.strftime("%c")
       name_time=gettime.split()[2]+gettime.split()[1]+gettime.split()[4]+'_'+gettime.split()[3].replace(':','_')
       line='touch '+name_time # This is just to check whether the calculation finished or not
       fileout.write(line)
    fileout.close()

    # Sumbitting the job
    if q_system=='slurm':
       cmd=parameters['sbatch']+' '+nameOut
       os.system(cmd)
   
    isfi=os.path.isfile(name_time)    
    while isfi==False: # the program will loop here until it finds te calculation_finished.ok file
      tiempo.sleep(1) 
      isfi=os.path.isfile(name_time)
    
    cmd='rm '+name_time
    os.system(cmd)

     
def combine_trajectories(iteration,parameters,restart):
    if restart==False:
       if parameters['md_engine']=='GROMACS':
     
          # Combine MD files
          totaltraj="totaltraj.trr"
          trajs=glob.glob('iteration'+str(iteration)+'*.trr')
          #print glob.glob('iteration'+str(iteration)+'*.trr')
          #print trajs
          fileout=open('c.txt','w')
          c_len=len(trajs)
          if os.path.isfile(totaltraj):
             c_len=c_len+1
          line='c\n'*(c_len)
          fileout.write(line)
          fileout.close()
          if not os.path.isfile(totaltraj):
             cmd='gmx trjcat -f '+' '.join(trajs)+' -o iteration'+str(iteration)+'.trr -keeplast -cat -settime < c.txt'
             os.system(cmd)
          else:
             cmd='gmx trjcat -f totaltraj.trr '+' '.join(trajs)+' -o iteration'+str(iteration)+'.trr -keeplast -cat -settime < c.txt'
             os.system(cmd)
          cmd='mv iteration'+str(iteration)+'.trr totaltraj.trr'
          os.system(cmd)
          if parameters['debug']==False:
             cmd='rm iteration'+str(iteration)+'*'
             os.system(cmd)
       
    # Combine PLUMED COLVAR files
    totalplumed='COLVAR'
    colvars=glob.glob('COLVAR.*')
    if not os.path.isfile(totalplumed):
       cmd='cat '+' '.join(colvars)+' > '+totalplumed
       os.system(cmd)
    else:
       cmd='cat '+totalplumed+' '+' '.join(colvars)+' > tmp.cv' 
       os.system(cmd)
       cmd='mv tmp.cv '+totalplumed
       os.system(cmd)
    if parameters['debug']==False:
       cmd='rm COLVAR.*'
       os.system(cmd)
    
    time,metric_list_val,cv_list_val,metric_avg,metric_sd,cv_avg,cv_sd=analyse_plumed_output(parameters,'COLVAR',iteration)
    
    #print time
   
    return cv_list_val, metric_list_val,time

def clustering(time,values,iteration,parameters):


    if parameters['clustering']=='density_Laio':
        #calculate euclidean distances:
       # print "Memory usage before calculating euclidMat"
       # print(psutil.virtual_memory()) #only for debug purposes
 
        print "Calculating euclidean distances"
        #print values
        #print values.shape
        euclidMat=spatial.distance.cdist(values,values,'euclidean')
        #print euclidMat
        #print euclidMat.shape
        euclidMat_stats=spatial.distance.pdist(values,'euclidean')
        #print euclidMat_stats
        euclidMat_avg=numpy.mean(euclidMat_stats)
        euclidMat_sd=numpy.std(euclidMat_stats)
        print "euclidMat_avg = ", euclidMat_avg
        print "euclidMat_sd = ", euclidMat_sd
        #sys.exit()

        #calculte rho for each data point
        if parameters['cluster_dc']=='value':
           d0=parameters['dc_value']
        if parameters['cluster_dc']=='avg-2sd':
           d0=euclidMat_avg-2*euclidMat_sd 
        elif parameters['cluster_dc']=='avg-sd':
           d0=euclidMat_avg-euclidMat_sd
        elif parameters['cluster_dc']=='avg':
           d0=euclidMat_avg
        
        print "Calculating rho"
        rho=[]
        for i in range(0,len(values)):
            rhoi=0.
            for j in range(0,len(values)):
                if j==i:
                   continue
                elif j<i:
                   #print j, i
                   if euclidMat[i][j]<d0:
                      rhoi=rhoi+1
                elif j>i:
                   if euclidMat[j][i]<d0:
                      rhoi=rhoi+1
            rho.append(rhoi)
#        print rho
#        sys.exit()
        
        print "calculating delta"
        #calculate delta for each data point:
        delta=[]
        for i in range(0,len(values)):
            deltai=999999999.
            for j in range(0,len(values)):
                if j==i:
                   continue
                elif j<i:
                   if (rho[j]>rho[i]) and (euclidMat[i][j]<deltai):
                      deltai=euclidMat[i][j]
                elif j>i:
                   if (rho[j]>rho[i]) and (euclidMat[j][i]<deltai):
                      deltai=euclidMat[j][i]
            delta.append(deltai)
#        print delta
#        sys.exit()
        # correct the value for the point with maximum density. 
        #This is done different than in the paper. In the paper
        # the delta value for this point is the distance with
        # the furthest point.
        deltai=0.
        for i in range(0,len(values)):
            if rho[i]==max(rho):
               print "max value ", rho[i], " found in position ", i
               disti=0.
               for j in range(0,len(values)):
                   if j==i:
                      continue
                   elif j<i:
                      distj=euclidMat[i][j]
                   elif j>i:
                      distj=euclidMat[j][i]
                   if distj>disti:
                      disti=distj
               delta[i]=disti    
#        print delta
#        sys.exit()
      
        fileout=open('rhodelta.txt','w')
        for i in range(1,len(values)):
            line=str(time[i])+' '+str(rho[i])+' '+str(delta[i])+'\n'
            fileout.write(line)
        fileout.close()

        rhodelta=numpy.array([time,rho,delta]).transpose().astype(float)
        #print rhodelta

        maxdelta=max(rhodelta[:,2])
        #print "Maximum density is:", maxdelta
        
        print "getting cluster centers"
        clusterCenters=[]
        for i in range(0,len(values)):
            if delta[i] >= maxdelta*0.95:
               if len(clusterCenters)==0: # The paper doesn't consider the case where 2 points of the dataset are the same. Added this trick.
                  clusterCenters.append(i)
#                  print i, time[i]
               else:
                  dis=[]
                  for j in clusterCenters:
                      dis.append(euclidMat[j][i])
                  if 0. not in dis:
#                     print i, j, time[i], dis
                     clusterCenters.append(i)
       
        #print "Found",len(clusterCenters), "cluster centers"

        print "Assigning snapshots to each cluster center"
        clusters={}
        for center in clusterCenters:
            #print "Found center", center, "at time=", time[center], "picoseconds."
            clusters[time[center]]=[]

        # Assign each point to the cluster corresponding to the closest cluster center
        for i in range(0,len(values)):
            disti=[]
            for j in clusterCenters:
                if j==i:
                    disti.append(0.)
                elif j>i:
                    disti.append(euclidMat[i][j])
                elif j<i:
                    disti.append(euclidMat[j][i])
            mindisti=min(disti)
            jmindist=disti.index(mindisti)
            center=clusterCenters[jmindist]
            clusters[time[center]].append(time[i])

    #check how many elements are in each cluster
    print "Deciding if a cluster has been well sampled or it's an outlier"
    totalElements=0
    clusters_forward={}
    outliers_forward={}
    clusters_file='clusters.'+str(iteration)+'.txt'
    fileout=open(clusters_file,'w')
    for center in clusters.keys():
        sampltime=len(clusters[center])*float(parameters['stride'])*float(parameters['dt'])
        #print fraction
        if sampltime >= float(parameters['sampltime']):
            line="Cluster with center at "+str(center)+ " picoseconds has been sampled for "+str(sampltime)+ " picoseconds. --CLUSTER\n"
            fileout.write(line)
            clusters_forward[center]=clusters[center]
        else:
            line="Cluster with center at "+str(center)+ " picoseconds has been sampled for "+str(sampltime)+ " picoseconds. --OUTLIER\n"
            fileout.write(line)
            outliers_forward[center]=clusters[center]

        totalElements=totalElements+len(clusters[center])
    fileout.close()
    NClust=len(clusters_forward.keys())
    Noutli=len(outliers_forward.keys())
#    print "Total elements clustered:", totalElements, ". Number of observations: ",len(values)," (MUST BE THE SAME)"
#    print "Number of clusters: ", NClust
#    print "Number of outliers: ", Noutli
#    print "---------------------------------------------------------"


#    print clusters_forward    
#    sys.exit()
    print "deleting big variables"
    del euclidMat, euclidMat_stats
    return clusters_forward,outliers_forward


def save_clusters(parameters,clusters,clustype,iteration):
    if parameters['md_engine']=='GROMACS':
       trr='totaltraj.trr'
       tpr=parameters['tpr']
       ndx=parameters['ndx']
       gmx=parameters['gmx_path']
       namedir='iter'+str(iteration)
       cmd='mkdir '+namedir
       os.system(cmd)
       for time in clusters.keys():
           time=str(int(float(time)))
           nameOut=namedir+'/'+clustype+'_'+time+'.pdb'
           if not os.path.isfile(nameOut): # the time of each center should be the same since we are combining iterations
              cmd=gmx+' trjconv -f '+trr+' -s '+tpr+' -n '+ndx+' -b '+str(int(time))+' -e '+str(int(time))+\
                  ' -pbc mol -ur compact -center -o '+nameOut+'<<OUT\n3\n0\n'
              os.system(cmd)
           

def generate_restarts(clusters,outliers,iteration,parameters):

    # Generate an inverted distribution for taboo search.
    # Every cluster/outlier will be given a normalised weight 
    # inverse to the weight in the real distribution (cite paper)
    all_clust={}

    # Combine clusters and outliers and calculate the total number of snapshots
    num_snaps=0

    for key in clusters.keys():
        all_clust[key]=clusters[key]
        num_snaps=num_snaps+len(all_clust[key])
    for key in outliers.keys():
        all_clust[key]=outliers[key]
        num_snaps=num_snaps+len(all_clust[key])

    # Assign to each cluster/outlier weight equal to the inverse of their statistical distribution
    inv_weights={}
    sum_inv_weights=0
    for key in all_clust.keys():
        inv_weights[key]=float(num_snaps)/float(len(all_clust[key]))
        sum_inv_weights=sum_inv_weights+inv_weights[key]
#        print "cluster "+str(key)+" has an inverse weight of "+str(inv_weights[key])
#    print  "the sum of inverse weights is: "+str(sum_inv_weights)

    # Normalise the inverse weights and draw times
    norm_weights={}
    for key in inv_weights.keys():
        norm_weights[key]=inv_weights[key]/float(sum_inv_weights)
#        print "cluster "+str(key)+" has a normalised weight of "+str(norm_weights[key])

    if len(inv_weights.keys())>=int(parameters['nsim']):
       times_restarts=numpy.random.choice(a=norm_weights.keys(),p=norm_weights.values(),size=int(parameters['nsim']),replace=False)
    else:
       rest=int(parameters['nsim'])-len(inv_weights.keys())
       times_restarts=numpy.random.choice(a=norm_weights.keys(),p=norm_weights.values(),size=len(inv_weights.keys()),replace=False)
       times_restarts=numpy.append(times_restarts, numpy.random.choice(a=norm_weights.keys(),p=norm_weights.values(),size=rest,replace=True))
    
    if parameters['md_engine']=='GROMACS':

       gmx=parameters['gmx_path']
       trr='totaltraj.trr'
       tpr=parameters['tpr']
       ndx=parameters['ndx']
       top=parameters['top']
       mdp=parameters['mdp']

       rep=0
       nameOut="restarts_iteration"+str(iteration+1)+".txt"
       fileout=open(nameOut,'w')

       for instant in times_restarts:
           time_str=str(int(float(instant)))
           gro='restart'+str(iteration)+'_'+time_str+'.gro'
           line="Replica "+str(rep)+" of Iteration "+str(iteration+1)+" will be started from the snapshot at "+time_str+" ps.\n"
           fileout.write(line)
           if not os.path.isfile(gro):
              cmd=gmx+' trjconv -f '+trr+' -s '+tpr+' -n '+ndx+' -b '+str(int(time_str))+' -e '+str(int(time_str))+' -o '+gro+'  <<EOF\n0\n'
              os.system(cmd)
           cmd=gmx+' grompp -f '+mdp+' -c '+gro+' -p '+top+' -n '+ndx+' -o '+'iteration'+str(iteration+1)+str(rep)+'.tpr'
           os.system(cmd)
           rep=rep+1
       fileout.close()
       if parameters['debug']==False:
          os.system('rm restart*gro *mdout*')
    
def calc_avg(iteration,time,clusters,cv_arr,metric_arr,metricAvgTarget,metricSDTarget,cvAvgTarget,cvSDTarget):
    
    nameOut='distance_target.'+str(iteration)
    fileout=open(nameOut,'w')

    fileout.write('Center Time CV_avg CV_sd Superposition_CV Dist_metric_avg Dist_metric_sd Superposition_metric')

    cv_avgs=[]
    cv_sds=[]
    metric_avgs=[]
    metric_sds=[]   

    for time_center in clusters.keys():
        fileout.write('\n')
        cv_vals=[]
        dist_vals=[]
        for time_element in clusters[time_center]:
            index=numpy.where(time==time_element)[0][0]
            cv_vals.append(cv_arr[index])
            dist=numpy.linalg.norm(metric_arr[index]-metricAvgTarget)
            dist_vals.append(dist)
        cv_vals=numpy.array(cv_vals).astype(float)
        dist_vals=numpy.array(dist_vals).astype(float)
        cv_avg=cv_vals.mean()
        cv_sd=cv_vals.std()
        dist_avg=dist_vals.mean()
        dist_sd=dist_vals.std()
        super_cv="NA" #Need to decide how to calculate it
        super_metric="NA"
        center=clusters.keys().index(time_center)
        line=str(center)+' '+str(time_center)+' '+\
             str(cv_avg)+' '+str(cv_sd)+' '+super_cv+' '+\
             str(dist_avg)+' '+str(dist_sd)+' '+super_metric
        fileout.write(line)
    fileout.close()
    
    

if __name__ == '__main__':

    ######## SET UP PARAMETERS AND RUN SANITY CHECKS #########
    
    parameters=parse(parser) # Get all parameters from input file

    parameters=setup_queue(parameters)

    parameters=genEnVar(parameters)
    
    parameters=getCV(parameters) #Set up the CV of use

    parameters=getMetric(parameters) # Check that the chosen metric is supported and can work with the chosen CV

    parameters=getBias(parameters)


   ######### PREPARE FILES TO BE USED  ######################
   
    metric_input=setupMetric(parameters) # Get the necessary parameters for the metric of use

    if parameters['target'] is not None:
       metricAvgTarget,metricSDTarget,cvAvgTarget,cvSDTarget=analyse_target(parameters,metric_input)
    else:
       metricAvgTarget,metricSDTarget,cvAvgTarget,cvSDTarget=None,None,None,None
    
    ##### RESTART CALCULATION IF NECESSARY ###################
    st_iter=int(parameters['st_iter'])
    if st_iter>0:
#       print "The restarting function is not complete yet. Exiting."
#       sys.exit()
       cv_arr,metric_arr,time=combine_trajectories(st_iter-1,parameters,True) 
       clusters,outliers=clustering(time,metric_arr,st_iter-1,parameters)
       save_clusters(parameters,clusters,'cluster',st_iter-1)
       save_clusters(parameters,outliers,'outlier',st_iter-1)
       generate_restarts(clusters,outliers,st_iter-1,parameters)
       if parameters['target'] is not None:
           calc_avg(st_iter,time,clusters,cv_arr,metric_arr,metricAvgTarget,metricSDTarget,cvAvgTarget,cvSDTarget) 
       nametpr='iteration'+str(st_iter)+'0.tpr'
       totaltraj="totaltraj.trr"
       colvarfile="COLVAR"
       names=[nametpr,totaltraj,colvarfile]
       for name in names:
          if not os.path.isfile(name):
             print "File "+name+" not found. Exiting"
             sys.exit()
#       print "The restarting function is not complete yet. Exiting."
#       sys.exit()

    ######## DO INTERESTING STUFF (WHEN IT IS NOT CRASHING) ####################

    max_iter=int(parameters['max_iter'])

    for iteration in range(st_iter,max_iter):
        if iteration>0:
           include=build_metric_bias(parameters,clusters,metric_arr,\
                                     iteration,metric_input,time)
        else:
           include=[]
           clusters={}

        taboo_plumedat=gen_plumed_input(parameters,include,clusters)
        print "Iteration "+str(iteration)+" is going to be submitted."
        qsub_md=submit_calc(parameters,iteration)  
        cv_arr,metric_arr,time=combine_trajectories(iteration,parameters,False)
        clusters,outliers=clustering(time,metric_arr,iteration,parameters)
        save_clusters(parameters,clusters,'cluster',iteration)
        save_clusters(parameters,outliers,'outlier',iteration)
        if parameters['target'] is not None:
           calc_avg(iteration,time,clusters,cv_arr,metric_arr,metricAvgTarget,metricSDTarget,cvAvgTarget,cvSDTarget)
        generate_restarts(clusters,outliers,iteration,parameters)

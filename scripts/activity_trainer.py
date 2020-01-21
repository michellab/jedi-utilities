import argparse, glob, os, sys, copy, math
import mdtraj
import numpy as np
import pandas as pd
import scipy as scp
from scipy import stats
import matplotlib
matplotlib.use('Agg') # set the backend before importing pyplot
import matplotlib.pyplot as plt
import math
    
def parser():
    """
    This should probably only take the name of a text file which has the PDB
    codes of the training set and their druggability value. This should be
    inside a folder where every subfolder has the protein, the ligand and the
    grid. maybe also the binding site files from jedi-setup.py?
    """
    parser = argparse.ArgumentParser(description="Optimise the parameters used to calculate the JEDI descriptors",
                                 epilog="jedi_trainer.py is distributed under the GPL.",
                                 prog="jedi-setup.py")
    parser.add_argument('-i','--input', nargs="?",
                        help='File containing PDB codes and druggability measures of the trained structures.\
                              Must be in a directory that contains a subdirectory for each PDB.')
    parser.add_argument('-c','--cutoff', nargs="?",type=str,default=0.6,
                        help='Maximum distance from the protein for a ligand atom to be considered\
                              Default will be 0.6.')

    args=parser.parse_args()

    return args

def get_atoms(system,lig):
    lig_indices=[]
    for atom in system.topology.atoms:
        if (atom.residue.name!=lig or atom.element==mdtraj.element.hydrogen):
           continue
        lig_indices.append(atom.index)

    neighbours= mdtraj.compute_neighbors(system, 0.6, lig_indices, haystack_indices=None, periodic=True)
    bsite_indices=[]
    for neighbour in neighbours[0]:
        if neighbour not in lig_indices:
           bsite_indices.append(neighbour)
    return lig_indices, bsite_indices

def get_distances(system,lig_indices,bsite_indices,thetas,pdb,cutoff):
    minimum_contact=[]
    mind={}
    for theta in thetas:
        mind[theta]=[]

    for lig_idx in lig_indices:
        atom_pairs=[]
        for bsite_idx in bsite_indices:
            if system.topology.atom(bsite_idx).element==mdtraj.element.hydrogen:
               continue
            atom_pairs.append([lig_idx,bsite_idx])
        distances = mdtraj.compute_distances(system, atom_pairs, periodic=True, opt=True)
        r_min=min(distances[0])
        if (r_min<0.15):
           print("R_min of ",r_min, "found for pdb ",pdb)
        elif (r_min > cutoff):
           continue
        minimum_contact.append(min(distances[0]))

        for theta in thetas:
            sum_exp=0
            for r in distances[0]:
                if (r>cutoff):
                   continue
                sum_exp += np.exp(theta/r)
            mind_k=theta/np.log(sum_exp)
            mind[theta].append(mind_k)

    return minimum_contact, mind

def S_off(v,v0,delta):
    m=(v-v0)/delta-(1-math.sqrt(2))/2
    if m<=0:
       s_off=1
    elif 0<m and m<1:
       s_off=((1-m**2)**2)*(1+2*m**2)
    else:
       s_off=0
    return s_off


def find_pcents(df,col):
    pcents=[]
    for pcent in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]:
        idx=int(len(df[col])*pcent)
        val=sorted(df[col])[idx]
        pcents.append(val)
        print (pcent, idx, val)
    return pcents

def get_contacts(system,lig_indices,bsite_indices,v0,delta,cutoff):
    contacts=[]
    for lig_idx in lig_indices:
        depth_j=0
        atom_pairs=[]
        for bsite_idx in bsite_indices:
            if system.topology.atom(bsite_idx).element==mdtraj.element.hydrogen:
               continue
            atom_pairs.append([lig_idx,bsite_idx])
        distances = mdtraj.compute_distances(system, atom_pairs, periodic=True, opt=True)
        for distance in distances[0]:
            if (distance > cutoff):
               continue
            depth_j+=S_off(distance,v0,delta)
        contacts.append(depth_j)

    return contacts

if __name__=="__main__":
    # Load structures as mdtraj objects
    args=parser()
    dataset=pd.read_csv(args.input,sep=" ")
    cutoff=args.cutoff

    r_min=[]

    thetas=xrange(1,101)
    mindist={}
    for theta in thetas:
        mindist[theta]=[]

    for i in xrange(0,len(dataset)):
        system=mdtraj.load(dataset["PDB"][i])
        lig_indices, bsite_indices=get_atoms(system,dataset["LIGAND"][i])
        
        minimum_contact, mind=get_distances(system,lig_indices,bsite_indices,thetas,dataset["PDB"][i],cutoff)
        r_min+=minimum_contact
        for theta in thetas:
            mindist[theta]+=mind[theta]
   

    dists=pd.DataFrame()
    dists["r_min"]=r_min
    r2=[]
    for theta in thetas:
        slope,intercept,rvalue,pvalue,stderr = scp.stats.linregress(r_min,mindist[theta])
        dists[str(theta)]=mindist[theta]
        r2.append(rvalue**2)
        print(rvalue**2) 
 
    corr=pd.DataFrame()
    corr["theta"]=thetas
    corr["r2"]=r2

    find_pcents(dists,'r_min')
    find_pcents(dists,'50')
         
    pcents=find_pcents(dists,'r_min')
    rmin_max=max(dists['r_min'])
    contacts_pcents=pd.DataFrame()
    for val in pcents:
        contacts=[]
        r_hydro=val
        deltar_hydro=rmin_max-val
        for i in xrange(0,len(dataset)):
            system=mdtraj.load(dataset["PDB"][i])
            lig_indices, bsite_indices=get_atoms(system,dataset["LIGAND"][i])
            contacts+=get_contacts(system,lig_indices,bsite_indices,r_hydro,deltar_hydro,cutoff)
        label="r_hydro="+str(r_hydro)+" deltar="+str(deltar_hydro)
        contacts_pcents[label]=contacts
 
    dists.to_pickle("dists.pkl")
    corr.to_pickle("corr.pkl")
    contacts_pcents.to_pickle("contacts.pkl")

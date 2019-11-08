import argparse, glob, os, sys, copy, math
import mdtraj
import numpy as np
import pandas as pd
import scipy as scp
from scipy import stats
    
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
    parser.add_argument('-d','--druggability', nargs="?",
                    help='Name of the value used as druggability measure. Must be present in --input as a column')
    parser.add_argument('-dt','--d_theta', nargs="?", default=0.1,
                    help='increment for theta optimisation')
    parser.add_argument('-tmi','--theta_max_iter', nargs="?", default=100,
                    help='Maximum number of iterations for the optimisation of theta')
    parser.add_argument('-tt','--theta_tolerance', nargs="?", default=0.0001,
                    help='Maximum number of iterations for the optimisation of theta')
    args=parser.parse_args()

    return args


def load_structures(dataset):
    grids=[]
    ligands=[]
    proteins=[]
    for pdb in dataset['PDB']:
        proteins.append(mdtraj.load(pdb+'/system_apo.pdb'))
        ligands.append(mdtraj.load(pdb+'/ligand.pdb'))
        grids.append(mdtraj.load(pdb+'/grid.pdb'))
    dataset['grid']=grids
    dataset['ligand']=ligands
    dataset['protein']=proteins
    return dataset
        


def get_r_min(point,atoms):
    r_min=np.inf
    for atom in atoms:
        r=np.linalg.norm(point-atom)
        if r<r_min:
            r_min=r
    return r_min

def calc_mindist(point, atoms, theta):
    sum_exp=0
    for atom in atoms:
        sum_exp += np.exp(theta/np.linalg.norm(point-atom))
    mindist=theta/np.log(sum_exp)
    return mindist
    

def theta_trainer(dataset,r2_target,tolerance,max_iter,d_theta):
    """
    This will get the best set of r2, slope and intercept between mindist
    and r_min. We want the line to be as close as possible to y=x. It returns
    the optimum value of theta and it should pickle a file with all theta, r2,
    slope and intercept so that a nice plot can be made.
    """

    thetas=pd.DataFrame()
    distances=pd.DataFrame()
    
    #Get real minimum distances gridpoint-protein
    r_min_lst=[]
    for i in range(0,len(dataset['PDB'])):
        protein=dataset['protein'][i].xyz[0]
        grid=dataset['grid'][i].xyz[0]
        for point in grid:
            r_min=get_r_min(point,protein)
            r_min_lst.append(r_min)
    distances['r_min']=r_min_lst

    #Calculate mindist according to Equation 6 JEDI paper for every trial theta value
    #and check the correlation with the actual values
    r2_lst=[]
    slope_lst=[]
    intercept_lst=[]
    theta_range=[]
    r2_max=0.0
    theta_selected=None
    
    theta=0.0
    n_iter=0
    while (abs(r2_max-r2_target)>tolerance):
        theta_range.append(theta)
        mindist_lst=[]
        for i in range(0,len(dataset['PDB'])):
            protein=dataset['protein'][i].xyz[0]
            grid=dataset['grid'][i].xyz[0]
            for point in grid:
                mindist=calc_mindist(point,protein,theta)
                mindist_lst.append(mindist)
        distances[str(theta)]=mindist_lst
        slope,intercept,rvalue,pvalue,stderr = scp.stats.linregress(mindist_lst,r_min_lst)
        slope_lst.append(slope)
        intercept_lst.append(intercept)
        r2=rvalue**2
        r2_lst.append(rvalue**2)
        print(n_iter, theta, r2)
        if (r2>r2_max):
            r2_max=r2
        theta=theta+d_theta
        n_iter=n_iter+1
        if n_iter>max_iter:
            print ("Theta optimisation reached ",max_iter, "iterations without convergence.")
            break
        
    thetas['theta']=theta_range
    thetas['r2']=r2_lst
    thetas['slope']=slope_lst
    thetas['intercept']=intercept_lst

    thetas.to_pickle("thetas.pkl")
    distances.to_pickle("distances.pkl")
    
    return theta_selected

def cc_min_trainer(dataset, theta):
    """
    Input: protein, grid, ligand, all as MDTraj objects
    """
    return 0

if __name__=="__main__":
    args=parser()
    dataset=pd.read_csv(args.input,sep=" ")
    dataset=load_structures(dataset)
    i=theta_trainer(dataset,r2_target=1.0,tolerance=args.theta_tolerance,max_iter=args.theta_max_iter,d_theta=args.d_theta)
    print(i)
    
 
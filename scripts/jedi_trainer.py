import argparse, os,sys, copy, math
import mdtraj
import numpy as np
import pandas as pd
    
def parser():
    """
    This should probably only take the name of a text file which has the PDB
    codes of the training set and their druggability value. This should be
    inside a folder where every subfolder has the protein, the ligand and the
    grid. maybe also the binding site files from jedi-setup.py?
    """
    args=0
    return args


def get_r_min(point,atoms):
    r_min=numpy.inf
    for atom in atoms:
        r=distance(point,atom)
        if r<r_min:
            r_min=r
    return r_min

def calc_mindist(point, atoms, theta):
    sum_exp=0
    for atom in atoms:
        sum_exp += np.exp(theta/distance(point,atom))
    mindist=theta/np.log(sum_exp)
    return mindist
    

def theta_trainer(dataset,r2_target=1.00,tolerance=0.001,maxIter=100000,theta_max=100.0):
    """
    This will get the best set of r2, slope and intercept between mindist
    and r_min. We want the line to be as close as possible to y=x. It returns
    the optimum value of theta and it should pickle a file with all theta, r2,
    slope and intercept so that a nice plot can be made.
    """
    #Get real minimum distances gridpoint-protein
    thetas=pd.DataFrame()
    r_min_lst=[]
    for system in dataset:
        for point in grid:
            r_min=get_r_min(point,site)
            r_min_lst.append(r_min)
            
    theta_range=np.linspace(0,theta_max,maxIter)
    
    theta_lst=[]
    r2_lst=[]
    slope_lst=[]
    intercept_lst=[]
    r2_max=0.0
    for theta in theta_range:
        mindist_lst=[]
        for system in dataset:
           for point in grid:
               mindist=calc_mindist(point,site,theta)
               mindist_lst.append(mindist)
        r2, slope, intercept = lin_reg(mindist_lst,r_min_lst)
        theta_lst.append(theta)
        r2_lst.append(r2)
        slope_lst.append(slope)
        intercept_lst.append(intercept)
        if (r2-r2_target)<=tolerance:
            theta_selected=theta
            break
        elif (r2>r2_max):
            r2_max=r2
            theta_selected=theta
            
    thetas.to_pickle("thetas.pkl")
    return theta_selected

def cc_min_trainer(dataset, theta):
    """
    Input: protein, grid, ligand, all as MDTraj objects
    
    """
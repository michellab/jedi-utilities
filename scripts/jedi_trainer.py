import argparse, glob, os, sys, copy, math
import mdtraj
import numpy as np
import pandas as pd
import scipy as scp
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
    
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
    parser.add_argument('-dt','--d_theta', nargs="?", default=0.1, type=float,
                    help='increment for theta optimisation')
    parser.add_argument('-tmi','--theta_max_iter', nargs="?", default=100, type=int,
                    help='Maximum number of iterations for the optimisation of theta')
    parser.add_argument('-tt','--theta_tolerance', nargs="?", default=0.0001, type=float,
                    help='Maximum number of iterations for the optimisation of theta')
    parser.add_argument('-lmd','--lig_mindist', nargs="?", default=0.2, type=float,
                    help='maximum distance gridpoint-ligand to consider them overlapping')
    parser.add_argument('-plc','--prot_lig_cutoff', nargs="?", default=0.5, type=float,
                    help='maximum distance protein-ligand to consider them interacting')
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

def get_points_ligand(dataset, lig_mindist, prot_lig_cutoff):
    dataset_points_ligand=[]
    for i in range(0,len(dataset['PDB'])):
        
        ligand_all=dataset['ligand'][i].xyz[0]
        protein=dataset['protein'][i].xyz[0]
        ligand=[]
        # Select ligand atoms that are close to the protein
        for atom in ligand_all:
            r_min=get_r_min(atom,protein)
            if r_min<prot_lig_cutoff:
               ligand.append(atom)
        
        # Select grid points that are close to those atoms
        grid=dataset['grid'][i].xyz[0]
        points_ligand=[]
        for j in range(0,len(grid)):
            r_min=get_r_min(grid[j],ligand)
            if r_min<=lig_mindist:
               points_ligand.append(j)
        dataset_points_ligand.append(points_ligand)
    dataset['points_ligand']=dataset_points_ligand
    return dataset

def cc_min_trainer(dataset, theta, percentage):
    mindist_active=[]
    for i in range(0,len(dataset['PDB'])):
        protein=dataset['protein'][i].xyz[0]
        grid=dataset['grid'][i].xyz[0]
        for j in range(0,len(grid)):
            if (j not in dataset['points_ligand'][i]):
                continue
            mindist=calc_mindist(grid[j],protein,theta) # This is used in the calculation of JEDI, so need to use Eq 6 from paper
            if (mindist<0.15):
                print(dataset['PDB'][i], j)
            mindist_active.append(mindist)
    mindist_active.sort(reverse=True) # This is a Soff, so a bigger value will mean less active points
    print("------------ CCmin+deltaCC --------------")
    for pcent in [0.2,0.4,0.6,0.8,1.0]:
        i=int(len(mindist_active)*pcent)-1
        cc_i=mindist_active[i]
        print(pcent, cc_i)
    pcent=percentage/100
    i=int(len(mindist_active)*pcent)-1
    cc=mindist_active[i]
    print("CCmin+deltaCC value that covers ", percentage, "\% of points overlapping with a ligand atom: ", cc)
    print("------------ CCmin+deltaCC --------------")

    #Plot a histogram for Thesis / paper purposes
    fig=plt.figure(figsize=(20,10))
    plt.hist(mindist_active, density=True, stacked=True, bins = 10)
    plt.xlabel("Minimum distance grid point - protein atom (nm)",fontsize=20)
    plt.ylabel("Probability", fontsize=20)
    plt.savefig("Mindist_hist.png",dpi=300)
    #FIXME add something here to select the CCmin+deltaCC we want to cover a certain percentage of the histogram
    return cc

def get_neighbours(dataset, GPmin=0.25, GPmax=0.35,max_allowed_neighbours=38):
    neighbours=[]
    max_neighbours=0
    for i in range(0,len(dataset['PDB'])):
        grid=dataset['grid'][i].xyz[0] #set of vectors of dimension 3
        points_ligand=dataset['points_ligand'][i] # set of integers (indices of gridpoints colliding with ligand atoms)
        neighbours_i=[] #list of lists that has the neighbours of each grid point that collides with a ligand
        for j in range(0,len(points_ligand)):
            neighbours_j_count=0
            neighbours_j=[]
            point_lig=grid[j]
            for k in range(0,len(grid)):
                point=grid[k]
                r=np.linalg.norm(point-point_lig)
                if (GPmin <= r) and (r<=GPmax):
                   neighbours_j.append(k)
                   neighbours_j_count += 1
                   if neighbours_j_count>max_allowed_neighbours:
                       sys.exit("Too many neighbours detected. You are doing something wrong.")
                   if neighbours_j_count>max_neighbours:
                      max_neighbours=neighbours_j_count   
            neighbours_i.append(neighbours_j)
        neighbours.append(neighbours_i)
    print("max neighbours: ", max_neighbours)
    dataset['Neighbours']=neighbours

    return dataset

def cc2_min_trainer(dataset,theta,percentage):
    mindist_active=[]
    for i in range(0,len(dataset['PDB'])):
        protein=dataset['protein'][i].xyz[0]
        grid=dataset['grid'][i].xyz[0]
        neighbours=dataset['Neighbours'][i]
        for neighbour_list in neighbours:
            for j in neighbour_list: # Notice that we are not caring who they are neighbors with
                mindist=calc_mindist(grid[j],protein,theta)
                mindist_active.append(mindist)

    mindist_active.sort()
    print("------------ CC2min+deltaCC2 --------------")
    for pcent in [0.2,0.4,0.6,0.8,1.0]:
        i=int(len(mindist_active)*pcent)-1
        cc_i=mindist_active[i]
        print(pcent, cc_i)
    pcent=percentage/100
    i=int(len(mindist_active)*pcent)-1
    cc2=mindist_active[i]
    print("CC2min (not deltaCC2 because that's a S_off) value that covers ", percentage, "\% of points overlapping with a ligand atom: ", cc2)
    print("------------ CC2min+deltaCC2 --------------")
    
    #Plot a histogram for Thesis / paper purposes
    fig=plt.figure(figsize=(20,10))
    plt.hist(mindist_active, density=True,stacked=True, bins = 10)
    plt.xlabel("Minimum distance grid point neighbour - protein atom (nm)",fontsize=20)
    plt.ylabel("Probability", fontsize=20)
    plt.savefig("Mindist2_hist.png",dpi=300)
    return cc2

def Emin_trainer(dataset,theta,cc2,percentage):
    num_neighbours=[]
    for i in range(0,len(dataset['PDB'])):
        protein=dataset['protein'][i].xyz[0]
        grid=dataset['grid'][i].xyz[0]
        neighbours=dataset['Neighbours'][i] #notice these are only the neighbours of gps overlapping the ligand
        for j in range(0,len(neighbours)): # Notice that we are not caring who they are neighbors with
            count_neighbours=0
            for k in neighbours[j]:
                mindist=calc_mindist(grid[k],protein,theta)
                if mindist<=cc2:
                    count_neighbours += 1
            num_neighbours.append(count_neighbours)
    
    num_neighbours.sort(reverse=True)
    print("------------ Emin+deltaE --------------")
    for pcent in [0.2,0.4,0.6,0.8,1.0]:
        i=int(len(num_neighbours)*pcent)-1
        E_i=num_neighbours[i]
        print(pcent, E_i)
    pcent=percentage/100
    i=int(len(num_neighbours)*pcent)-1
    E=num_neighbours[i]
    print("Emin+deltaE value that covers ", percentage, "\% of points overlapping with a ligand atom: ", cc2)
    print("------------ Emin+deltaE --------------")
    
    #Plot a histogram for Thesis / paper purposes
    fig=plt.figure(figsize=(20,10))
    plt.hist(num_neighbours, density=True,stacked=True, bins = 10)
    plt.xlabel("Exposure in number of grid points",fontsize=20)
    plt.ylabel("Probability", fontsize=20)
    plt.savefig("Exposure_hist.png",dpi=300)


if __name__=="__main__":
    # Load structures as mdtraj objects
    args=parser()
    dataset=pd.read_csv(args.input,sep=" ")
    dataset=load_structures(dataset)

    # Optimise theta so that the calculated mindist is as well correlated as possible with the actual minimum distance gridpoint-atom
    #theta=theta_trainer(dataset,r2_target=1.0,tolerance=args.theta_tolerance,max_iter=args.theta_max_iter,d_theta=args.d_theta)
    
    # Optimise CCmin+deltaCC so that Son_mind is equal to 1 for the points within 0.2 nm of any ligand heavy atom
    # 1) Get a list of grid points within 0.2 nm of any ligand heavy atoms.
    dataset=get_points_ligand(dataset, args.lig_mindist,args.prot_lig_cutoff)
    #CC=cc_min_trainer(dataset, 5.0,80.0)
    dataset=get_neighbours(dataset)
    #CC2=cc2_min_trainer(dataset,5.0,80.0)
    E=Emin_trainer(dataset,5.0,0.21833799852702668,80.0)
    print(dataset)

    
 
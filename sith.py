
# coding: utf-8

# In[1]:

import argparse
import mdtraj
import numpy, scipy, math
import glob, os, sys


# In[2]:

os.chdir('/home/joan/PhD/year_2/hPNMT/MD/1HNN_apo/taboo/')


# In[3]:

def findResidues(apolar,polar):
    residues_str=[]
    for pdb in [apolar,polar]:
        structure=mdtraj.load_pdb(pdb)
        for res in structure.topology.residues:
            if str(res) not in residues_str and "GLY" not in str(res):
               residues_str.append(str(res))
    print "There are", len(residues_str), "residues"
    return residues_str


# In[4]:

def findAtomGroups(structure,residues):
    struct=mdtraj.load_pdb(structure)

    if metric == "COM":
        sideChains={}
        for res in struct.topology.residues:
            if str(res) in residues:
               sideChains[str(res)]=[]
               for atom in res.atoms:
                   if atom.is_sidechain==True and str(atom.element) is not "hydrogen":
                      sideChains[str(res)].append(str(atom.serial))
        return sideChains

    elif metric == "TORSION":
        chi={}
        for i in range(0,len(mdtraj.compute_chi1(struct)[0])):
            name='chi1_'+str(i)
            atoms=[]
            for atom in mdtraj.compute_chi1(struct)[0][i]:
                atoms.append(str(atom))
            chi[name]=atoms

        for i in range(0,len(mdtraj.compute_chi2(struct)[0])):
            name='chi2_'+str(i)
            atoms=[]
            for atom in mdtraj.compute_chi2(struct)[0][i]:
                atoms.append(str(atom))
            chi[name]=atoms

        for i in range(0,len(mdtraj.compute_chi3(struct)[0])):
            name='chi3_'+str(i)
            atoms=[]
            for atom in mdtraj.compute_chi3(struct)[0][i]:
                atoms.append(str(atom))
            chi[name]=atoms

        for i in range(0,len(mdtraj.compute_chi4(struct)[0])):
            name='chi4_'+str(i)
            atoms=[]
            for atom in mdtraj.compute_chi4(struct)[0][i]:
                atoms.append(str(atom))
            chi[name]=atoms

        return chi


# In[5]:

def monitor_metric_plumedat(atomGroups,target,iteration,replica,torsions_target,labels):
    if metric=='COM':
        sys.exit("COM bias is not implemented yet")
        return 0
    elif metric=='TORSION':
        plumedPrintFile='colvar_torsions.dat.'+str(iteration)#+'.'+str(replica)
        jePrintFile='JEDI.'+str(iteration) 
        nameOut='torsion_monitor.dat.'+str(iteration)#+'.'+str(replica)
        fileout=open(nameOut,'w')
        line='je: JEDI APOLAR='+apolar+' POLAR='+polar+' GRID='+grid+' PARAMETERS=jedi.params'+' STRIDE=1 SUMMARY=jedi_stats.dat.'+str(iteration)+'.'+str(replica)+' GRIDSTRIDE='+str(stride)+'\n'
        fileout.write(line)
        printarg=[]
        for key in atomGroups.keys():
            line=key+': TORSION ATOMS='+','.join(atomGroups[key])+'\n'
            fileout.write(line)
            line='sin_'+key+': MATHEVAL ARG='+key+' FUNC=sin(x) PERIODIC=NO'+'\n'
            printarg.append('sin_'+key)
            fileout.write(line)
            line='cos_'+key+': MATHEVAL ARG='+key+' FUNC=cos(x) PERIODIC=NO'+'\n'
            fileout.write(line)
            printarg.append('cos_'+key)
        line='PRINT ARG='+','.join(printarg)+' STRIDE=1 FILE='+plumedPrintFile+'\n'
        fileout.write(line)
        line='PRINT ARG=je STRIDE=1 FILE=JEDI.'+str(iteration)#+'.'+str(replica)+'\n'
        fileout.write(line)
        fileout.close()
        if (target is not None) and (torsions_target is None) and (labels is None):
            if target.split('.')[1]=='pdb':
               torsions_target_cmd='plumed driver --mf_pdb '+target+' --plumed '+nameOut
            elif target.split('.')[1]=='gro':
               torsions_target_cmd='plumed driver --igro '+target+' --plumed '+nameOut
            else:
                sys.exit('The target has to be in either pdb or gro format')
            print(torsions_target_cmd)
            os.system(torsions_target_cmd)
            cmd='mv '+plumedPrintFile+' torsions_target.dat && mv '+jePrintFile+' je_target.dat'
            os.system(cmd)
            filein=open('torsions_target.dat','r')
            for line in filein:
                if line.startswith('#'):
                    line=line.split()
                    labels=line[3:]
                else:
                    line=line.split()
                    torsions_target=numpy.array(line[1:]).astype(float)
        elif target is None:
            torsions_target=None
            
    return nameOut,plumedPrintFile,torsions_target,labels


# In[6]:

def jedi_taboo_plumedat(monitor_cv,iteration,replica,apolar,polar,grid,stride,kappa,at,mts,include):
    if iteration == 0:
        fileout=open('jedi_taboo.dat','w')
        line='je: JEDI APOLAR='+apolar+' POLAR='+polar+' GRID='+grid+' PARAMETERS=jedi.params'+            ' STRIDE='+str(stride)+' SUMMARY=jedi_stats.dat.'+str(iteration)+'.'+str(replica)+' GRIDSTRIDE='+str(stride)+'\n'
        fileout.write(line)
        line='resje: LOWER_WALLS ARG=je AT='+str(at)+' KAPPA='+str(kappa)+' STRIDE='+str(mts)+'\n'
        fileout.write(line)
        #line='INCLUDE FILE='+monitor_cv+'\n'
        #fileout.write(line)
        #line='PRINT ARG=je STRIDE='+str(stride)+' FILE=JEDI.'+str(iteration)+'.'+str(replica)
        #fileout.write(line)
        fileout.close()
    else:
        nameOut='jedi_taboo.'+str(iteration)+'.'+str(replica)+'.dat'
        fileout=open(nameOut,'w')
        line='je: JEDI APOLAR='+apolar+' POLAR='+polar+' GRID='+grid+' PARAMETERS=jedi.params'+            ' STRIDE='+str(stride)+' SUMMARY=jedi_stats.dat.'+str(iteration)+'.'+str(replica)+' GRIDSTRIDE='+str(stride)+'\n'
        fileout.write(line)
        line='resje: LOWER_WALLS ARG=je AT='+str(at)+' KAPPA='+str(kappa)+' STRIDE='+str(mts)+'\n'
        fileout.write(line)
        #line='INCLUDE FILE='+monitor_cv+'\n'
        #fileout.write(line)
        #line='PRINT ARG=je STRIDE='+str(stride)+' FILE=JEDI.'+str(iteration)+'.'+str(replica)+'\n'
        #fileout.write(line)
        for bias in include:
            line='INCLUDE FILE='+bias+'\n'
            fileout.write(line)
        fileout.close()
    return include
        


# In[7]:

def analysis_output(colvar_je,colvar_tor,residues,torsions_target,strideps):
    #get JEDI and time values
    time=[]
    time_line=0.
    jedi=[]
    filein=open(colvar_je,'r')
    for line in filein:
        if line.startswith('#'):
            continue
        else:
            time.append(time_line*strideps)
            time_line += 1
            jedi.append(line.split()[1])
    time=numpy.array(time)
    jedi=numpy.array(jedi)
    
    #get the torsions and torsional distances with the target if any
    if torsions_target is not None:
        torsion_dist=[]
    else:
        torsion_dist=None
    torsions_arr=[]
    filein=open(colvar_tor,'r')
    for line in filein:
        if line.startswith('#'):
            continue
        else:
            line=line.split()
            torsions=numpy.array(line[1:]).astype(float)
            torsions_arr.append(torsions)
            if torsions_target is not None:
                print type(torsions[0]), len(torsions)
                print type(torsions_target[0]), len(torsions_target)
                dist=numpy.linalg.norm(torsions-torsions_target)
                torsion_dist.append(dist)
    torsions_arr=numpy.array(torsions_arr)
    return time, jedi, torsions_arr, torsion_dist    


# In[8]:

def Clustering(time,values,iteration):
    #calculate euclidean distances:
    print "Calculating euclidean distances"
    euclidMat=[]
    euclidMat_stats=[]
    for i in range (0,len(values)):
        euclid=[]
        for j in range(0,len(values)):
            if i==j:
               euclid.append(99999999.)
            elif j>i:
               euclid.append(numpy.linalg.norm(values[i]-values[j]))
               euclidMat_stats.append(numpy.linalg.norm(values[i]-values[j]))
            elif j<i:
               euclid.append(99999999.)
        euclidMat.append(euclid)
    euclidMat=numpy.array(euclidMat)
    #print euclidMat
    euclidMat_avg=numpy.mean(euclidMat_stats)
    euclidMat_sd=numpy.std(euclidMat_stats)


    print "euclidMat_avg = ", euclidMat_avg
    print "euclidMat_sd = ", euclidMat_sd
    
    #calculte rho for each data point
    if metric == "COM":
       d0=euclidMat_avg/2
    elif metric == "TORSION":
       d0=euclidMat_avg-2*euclidMat_sd #FIXME: depends on euclidmat. Should find a standard.
    rho=[]
    for i in range(0,len(values)):
        rhoi=0.
        for j in range(0,len(values)):
            if j==i:
               continue
            elif j>i:
               #print j, i
               if euclidMat[i][j]<d0:
                  rhoi=rhoi+1
            elif j<i:
               if euclidMat[j][i]<d0:
                  rhoi=rhoi+1
        rho.append(rhoi)
    #print rho

    #calculate delta for each data point:
    delta=[]
    for i in range(0,len(values)):
        deltai=999999999.
        for j in range(0,len(values)):
            if j==i:
               continue
            elif j>i:
               if (rho[j]>rho[i]) and (euclidMat[i][j]<deltai):
                  deltai=euclidMat[i][j]
            elif j<i:
               if (rho[j]>rho[i]) and (euclidMat[j][i]<deltai):
                  deltai=euclidMat[j][i]
        delta.append(deltai)
        
    # correct the value for the point with maximum density. 
    #This is done different than in the paper. In the paper
    # the delta value for this point is the distance with
    # the furthest point.
    deltai=0.
    for i in range(0,len(delta)):
        if delta[i]!=999999999. and delta[i]>deltai:
            deltai=delta[i]
    for i in range(0,len(delta)):
        if delta[i]==999999999.:
            delta[i]=deltai

    fileout=open('rhodelta.txt','w')
    for i in range(1,len(values)):
        line=str(time[i])+' '+str(rho[i])+' '+str(delta[i])+'\n'
        fileout.write(line)
    fileout.close()
    
    rhodelta=numpy.array([time,rho,delta]).transpose().astype(float)
    print rhodelta
    
    maxdelta=max(rhodelta[:,2])
    print "Maximum density is:", maxdelta

    clusterCenters=[]
    for i in range(0,len(values)):
        if delta[i] >= maxdelta*0.75:
           clusterCenters.append(i)
    print "Found",len(clusterCenters), "cluster centers"

    clusters={}
    for center in clusterCenters:
        print "Found center", center, "at time=", time[center], "picoseconds."
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
        jmindist=disti.index(min(disti))
        center=clusterCenters[jmindist]
        clusters[time[center]].append(time[i])

    #check how many elements are in each cluster
    totalElements=0
    clusters_forward={}
    outliers_forward={}
    clusters_file='clusters.'+str(iteration)+'.txt'
    fileout=open(clusters_file,'w')
    for center in clusters.keys():
        ps_snap=stride*dt
        fraction=len(clusters[center])*ps_snap
        #print fraction
        if fraction >= 5000:
            line="Cluster with center at "+str(center)+ "picoseconds has been sampled for "+str(fraction)+ " picoseconds. --CLUSTER\n"
            fileout.write(line)
            clusters_forward[center]=clusters[center]
        else:
            line="Cluster with center at "+str(center)+ "picoseconds has been sampled for "+str(fraction)+ " picoseconds. --OUTLIER\n"
            fileout.write(line)
            outliers_forward[center]=clusters[center]

        totalElements=totalElements+len(clusters[center])
    fileout.close()
    NClust=len(clusters_forward.keys())
    Noutli=len(outliers_forward.keys())
    print "Total elements clustered:", totalElements, ". Number of observations: ",len(values)," (MUST BE THE SAME)"
    print "Number of clusters: ", NClust
    print "Number of outliers: ", Noutli
    print "---------------------------------------------------------"

   
#    print clusters_forward    
#    sys.exit()
    return clusters_forward,outliers_forward
       


# In[9]:

def calc_avg(clusters,time,jedi,torsion_dist,iteration):
    nameOut='distance_target.'+str(iteration)
    fileout=open(nameOut,'w')
    fileout.write('Center Time jedi_avg jedi_sd dist_avg dist_sd')
    for time_center in clusters.keys():
        fileout.write('\n')
        jedi_vals=[]
        dist_vals=[]
        for time_element in clusters[time_center]:
            index=numpy.where(time==time_element)[0][0]
            print index
            jedi_vals.append(jedi[index])
            dist_vals.append(torsion_dist[index])
        jedi_vals=numpy.array(jedi_vals).astype(float)
        dist_vals=numpy.array(dist_vals).astype(float)
        jedi_avg=jedi_vals.mean()
        jedi_sd=jedi_vals.std()
        dist_avg=dist_vals.mean()
        dist_sd=dist_vals.std()
        center=clusters.keys().index(time_center)
        line=str(center)+' '+str(time_center)+' '+str(jedi_avg)+' '+str(jedi_sd)+' '+str(dist_avg)+' '+str(dist_sd)
        fileout.write(line)
    fileout.close()
            


# In[10]:

def save_clusters(clusters,iteration,ndx):
    trr='iteration'+str(iteration)+'.trr'
    tpr='iteration'+str(iteration)+'.0.tpr'
    for time in clusters.keys():
        time=str(int(float(time)))
        nameOut='iteration'+str(iteration)+'_'+time+'.pdb'
        cmd='gmx trjconv -f '+trr+' -s '+tpr+' -n '+ndx+' -b '+str(int(time))+' -e '+str(int(time))+            ' -pbc mol -ur compact -center -o '+nameOut+'<<OUT\n3\n0\n'
        os.system(cmd)


# In[20]:

def build_biasfile(atomGroups,labels,clusters,time,torsions,iteration,include):
    # find average distance between clusters. We'll push away with a lower wall at this distance
    if len(clusters.keys())==0:
        return include
    clusdist=[]
    if len(clusters.keys())==1:
        clusdist_avg=0
    else:
        for time_center1 in clusters.keys():
            time_index1=numpy.where(time==time_center1)[0][0]
            #print time_center, time_index, time[time_index]
            for time_center2 in clusters.keys():
                if time_center2==time_center1:
                    continue
                time_index2=numpy.where(time==time_center2)[0][0]
                dist=numpy.linalg.norm(torsions[time_index1]-torsions[time_index2])
                clusdist.append(dist)
        clusdist_avg=numpy.mean(numpy.array(clusdist).astype(float))
    
    #build a bias file from this set of torsions
    
    #compute torsions, sines and cosines
    nameOut='torsion_bias.'+str(iteration)+'.dat'
    fileout=open(nameOut,'w')
    #for key in atomGroups.keys(): # all this is already in torsion monitor
        #line=key+': TORSION ATOMS='+','.join(atomGroups[key])+'\n'
        #fileout.write(line)
        #line='sin_'+key+': MATHEVAL ARG='+key+' FUNC=sin(x) PERIODIC=NO'+'\n'
        #fileout.write(line)
        #line='cos_'+key+': MATHEVAL ARG='+key+' FUNC=cos(x) PERIODIC=NO'+'\n'
        #fileout.write(line)
        
    #compute distances between torsion vector and references    
    dist_labs=[]
    
    for time_center in clusters.keys():
        time_index=numpy.where(time==time_center)[0][0]
        labelstemp=[]
        loc=0 # this is to locate the reference values for sin and cos
        for key in atomGroups.keys():
            line='sin_'+key+'_'+str(time_index)+': MATHEVAL ARG='+key+' FUNC=sin(x)-'+torsions[time_index][loc].astype(str)+' PERIODIC=NO'+'\n'
            #print line
            fileout.write(line)
            loc += 1
            labelstemp.append('sin_'+key+'_'+str(time_index))
            line='cos_'+key+'_'+str(time_index)+': MATHEVAL ARG='+key+' FUNC=cos(x)-'+torsions[time_index][loc].astype(str)+' PERIODIC=NO'+'\n'
            fileout.write(line)
            #print line
            loc += 1
            labelstemp.append('cos_'+key+'_'+str(time_index))
        line='dist2_'+str(iteration)+'_'+'_'+str(time_index)+': COMBINE ARG='+','.join(labelstemp)+             ' POWERS='+','.join(['2']*len(torsions[time_index]))+' PERIODIC=NO\n'
        '''
        if using plumed 2.3 or later, put the following line between dist2 and POWERS:
        ' PARAMETERS='+','.join(torsions[time_index].astype(str))+\
        and remove the subtracted elements from sin and cos and readjust the shape of the loop
        '''
        fileout.write(line)
        line='dist_'+str(iteration)+'_'+str(time_index)+\
             ': MATHEVAL ARG=dist2_'+str(iteration)+'_'+'_'+str(time_index)+\
             ' FUNC=sqrt(x) PERIODIC=NO\n'
        fileout.write(line)
        dist_labs.append('dist_'+str(iteration)+'_'+str(time_index))
    
    #Build the bias for all distances in the same line
    print clusdist_avg, clusdist_avg
    line='rest_'+str(iteration)+': LOWER_WALLS ARG='+','.join(dist_labs)+\
        ' AT='+','.join([str(clusdist_avg)]*len(dist_labs))+' KAPPA='+','.join([str(kappa)]*len(dist_labs))+'\n'
    fileout.write(line)
    fileout.close
    
    include.append(nameOut)
    return include
    


# In[12]:

def combine_trajectories(iteration,monitor_file):
    TotalTraj="totaltraj.trr"
    trajs=glob.glob('iteration'+str(iteration)+'*.trr')
    fileout=open('c.txt','w')
    line='c\n'*len(trajs+1)
    fileout.write(line)
    fileout.close()
    if iteration==0:
       cmd='mv iteration0.0.trr totaltraj.trr'
       os.system(cmd)
    else:
       cmd='gmx trjcat -f totaltraj.trr '+' '.join(trajs)+' -o iteration'+str(iteration)+'.trr -cat -settime < c.txt'
       os.system(cmd)
       cmd="mv iteration"+str(iteration)+" totaltraj.trr"
    #z=raw_input('continue? ')
    colvar_je='JEDI.'+str(iteration)
    #cmd='cat JEDI.'+str(iteration)+'.* > JEDI.'+str(iteration)
    cmd='plumed driver --mf_trr totaltraj.trr --plumed '+monitor_file
    os.system(cmd)
    colvar_tor='colvar_torsions.dat.'+str(iteration)
    #cmd='cat colvar_torsions.dat.'+str(iteration)+'.* > colvar_torsions.dat.'+str(iteration)
    #os.system(cmd)
    return colvar_je, colvar_tor
    
    
    


# In[13]:

def generate_restart(outliers,iteration,mdp,top,ndx):
    trr='iteration'+str(iteration)+'.trr'
    tpr='iteration'+str(iteration)+'.0.tpr' #tpr file is not important
    restarts=[]
    replica_index=[]
    for time in outliers.keys():
        time_str=str(int(float(time)))
        nameOut='restart'+str(iteration)+'_'+time_str+'.gro'
        cmd='gmx trjconv -f '+trr+' -s '+tpr+' -n '+ndx+' -b '+str(int(time_str))+' -e '+str(int(time_str))+' -o '+nameOut+'<<EOF\n0\n'
        os.system(cmd)
        cmd='gmx grompp -f '+mdp+' -c '+nameOut+' -p '+top+' -n '+ndx+' -o '+'iteration'+str(iteration+1)+            '.'+str(outliers.keys().index(time))+'.tpr'
        os.system(cmd)
        restarts.append('iteration'+str(iteration+1)+'.'+str(outliers.keys().index(time))+'.tpr')
        replica_index.append(str(outliers.keys().index(time)))
    return restarts,replica_index
        


# In[23]:

apolar='apolar.pdb'
polar='polar.pdb'
grid='grid.pdb'
stride=1250
dt=0.004
strideps=stride*dt
kappa=500
mts=4
at=7.3
metric='TORSION'
structure='1HNN.pdb'
target='2G8N_equil.pdb'
stride=1250
maxIter=100
ntomp=4
nthreads=20
gpu=4
ndx='1HNN_index.ndx'
mdp='production_50ns_4fs.mdp'
top='1HNN.top'


# In[15]:

if __name__=='__main__':
    residues = findResidues(apolar,polar)
    atomGroups = findAtomGroups(structure,residues)    


# In[ ]:

for iteration in range(0,maxIter):
    if iteration==0:
        #continue
        monitor_cv,monitor_out,torsions_target,labels=monitor_metric_plumedat(atomGroups,target,iteration,0,None,None)
        #sys.exit()
        include=jedi_taboo_plumedat(monitor_cv,iteration,0,apolar,polar,grid,stride,kappa,at,mts,[])
        tpr='iteration0.0.tpr'
        if not os.path.isfile(tpr):
            sys.exit("gromacs tpr file not found. ABORTING.")
        if nthreads > 16:
           ntomp_iter0=16
        else:
           ntomp_iter0=nthreads
        cmd='gmx mdrun -v -deffnm '+tpr.split('.')[0]+'.'+tpr.split('.')[1]+' -plumed jedi_taboo.dat -ntomp '+str(ntomp_iter0)+' -gpu_id 0 -nsteps 1250000'
        os.system(cmd)
        #cmd='mv iteration0.0.trr iteration0.trr'
        #os.system(cmd)
        #colvar_je='JEDI.0'
        #colvar_tor='colvar_torsions.dat.0'
        colvar_je,colvar_tor=combine_trajectories(iteration,monitor_cv)
        print colvar_je, colvar_tor
        #wait=raw_input("continue")
        time,jedi,torsions,torsion_dist=analysis_output(colvar_je,colvar_tor,residues,torsions_target,strideps)
        clusters,outliers=Clustering(time,torsions,iteration)
        save_clusters(clusters,iteration,ndx)
        calc_avg(clusters,time,jedi,torsion_dist,iteration)
        cmd='rm \#* *-step-*'
        os.system(cmd)
    else:
        restarts,replica_index=generate_restart(outliers,iteration-1,mdp,top,ndx)
        include=build_biasfile(atomGroups,labels,clusters,time,torsions,iteration-1,include)
        print "Iteration "+str(iteration)+': running '+str((len(replica_index)))+' replicas with '+str(ntomp)+' cores each.'
        nameOut='iteration'+str(iteration)+'.sh'
        
        fileout=open(nameOut,'w')
        nsim=0
        gpu_id=0
        for replica in replica_index: # so for each outlier we restart from
            if nsim*ntomp >= nthreads:
                fileout.write('wait\n')
                nsim=0
            nsim += 1
            monitor_cv,monitor_out,torsions_target,labels=monitor_metric_plumedat(atomGroups,target,iteration,replica,torsions_target,labels)
            jedi_taboo_plumedat(monitor_cv,iteration,replica,apolar,polar,grid,stride,kappa,at,mts,include)
            nametpr='iteration'+str(iteration)+'.'+str(replica)
            nameplumed='jedi_taboo.'+str(iteration)+'.'+str(replica)+'.dat'
            line='gmx mdrun -v -deffnm '+nametpr+' -plumed '+nameplumed+' -ntomp '+str(ntomp)+' -gpu_id '+str(gpu_id)+' -nsteps 1250000 &\n'
            fileout.write(line)
            gpu_id=gpu_id+1
            if gpu_id > gpu-1:
               gpu_id=0
        fileout.write('wait')
        fileout.close()
        cmd='bash '+nameOut
        os.system(cmd)
        colvar_je,colvar_tor=combine_trajectories(iteration,monitor_cv)
        cmd='rm iteration'+str(iteration)+'.*.trr'
        os.system(cmd)
        time,jedi,torsions,torsion_dist=analysis_output(colvar_je,colvar_tor,residues,torsions_target,strideps)
        clusters,outliers=Clustering(time,torsions,iteration)
        save_clusters(clusters,iteration,ndx)
        calc_avg(clusters,time,jedi,torsion_dist,iteration)
        cmd='rm \#* *-step-*'
        os.system(cmd)


# In[ ]:




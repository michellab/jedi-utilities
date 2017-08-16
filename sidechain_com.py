import argparse,numpy,os,mdtraj,scipy.stats,sys

parser = argparse.ArgumentParser(description="Takes a protein binding site and generates a \
                                              PLUMED input to calculate the COM of each side\
                                              chain and the distances between COMs",
                                 epilog="something.py is distributed under the GPL.",
                                 prog="something.py")
parser.add_argument('-i','--input', nargs="?",
                    help='Name of the pdb file containing the whole protein structure')
parser.add_argument('-a','--apolar', nargs="?",
                    help='Name of the pdb file containing the apolar atoms of the binding site',
                    default='apolar.pdb')
parser.add_argument('-p','--polar', nargs="?",
                    help='Name of the pdb file containing the polar atoms of the binding site',
                    default='polar.pdb')
parser.add_argument('-o','--output', nargs="?",
                    help='Name of the input file to supply to PLUMED',
                    default='com.dat')
parser.add_argument('-s','--stride', nargs="?",
                    help='PLUMED printing stride in ps (not timesteps!!)',
                    default='1')
parser.add_argument('-t','--trr', nargs="?",
                    help='gromacs trr trajectory to process with PLUMED driver')
parser.add_argument('-g','--gro', nargs="?",
                    help='Name of input GRO file.')
parser.add_argument('-d','--dat', nargs="?",
                    help='Name of PLUMED input DAT file.')
parser.add_argument('-r','--rmsd', nargs="?",
                    help='target RMSD value')
parser.add_argument('-k','--kappa', nargs="?",
                    help='force constant for the harmonic bias')
parser.add_argument('-m','--metric',nargs='?',
                    help="metric used to cluster points")
def parse(parser):
    args = parser.parse_args()
    return args.apolar,args.polar,args.input,args.output,args.stride,args.trr,args.gro,args.dat,args.rmsd,args.kappa,args.metric



def findResidues(apolar,polar): 
    residues_str=[]
    for pdb in [apolar,polar]:
        structure=mdtraj.load_pdb(pdb)
        for res in structure.topology.residues:
            if str(res) not in residues_str and "GLY" not in str(res):
               residues_str.append(str(res))
    print "There are", len(residues_str), "residues" 
    return residues_str

def findAtomGroups(structure,residues): #change this function to do it with mdtraj
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
 
def genPlumedinput(atomGroups,output,stride):
    fileout=open(output,'w')
    if metric == "COM":
        labels_com=[]
        for key in atomGroups.keys():
            label=key
            line=label+': COM ATOMS='+','.join(atomGroups[key])+'\n'
            fileout.write(line)
            labels_com.append(label)
    
        labels_dist=[]
        for i in range(0,len(labels_com)):
            for j in range(i,len(labels_com)):
                if j!=i:
                   label=labels_com[i]+'_'+labels_com[j]
                   line=label+': DISTANCE ATOMS='+labels_com[i]+','+labels_com[j]+'\n'
                   fileout.write(line)
                   labels_dist.append(label)
        print "there are ",len(labels_dist), "distances to print"

    elif metric == "TORSION":
        labels_tor=[]
        for key in atomGroups.keys():
            label=key
            line=label+': TORSION ATOMS='+','.join(atomGroups[key])+'\n'
            fileout.write(line)
            labels_tor.append(label)

        labels_dist=labels_tor
        print "there are ",len(labels_dist), "distances to print"

         # WORK OUT IF THERE ARE TORSIONS TO PLAY WITH AND HOW TO PRINT THEM

    line='PRINT ARG='+','.join(labels_dist)+' STRIDE=1 FILE=COMCOLVAR'
    fileout.write(line)    
    
    fileout.close()

def Plumed(trr,output,stride):
    command='plumed driver --mf_trr '+trr+' --plumed '+output
    os.system(command)
   
    # Get values from PLUMED output
    time=[]
    values=[]
    print "Putting values into vectors"
    filein=open('COMCOLVAR','r')
    for line in filein:
        line=line.split()
        if (line[0]=='#!'): 
           if 'SET' in line:
               continue
           labels=line[3:]
        else:
           time.append(float(line[0])*float(stride))
           values.append(map(float,line[1:]))
    values=numpy.array(values)
    print values
    #sys.exit()
    return labels, time, values



def Clustering(labels,time,values): # Done as in Rodriguez & Laio, Science(2014),344,6191,1492-1496 :
    
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
    print euclidMat
    euclidMat_avg=numpy.mean(euclidMat_stats)
    euclidMat_sd=numpy.std(euclidMat_stats)


    print "euclidMat_avg = ", euclidMat_avg
    print "euclidMat_sd = ", euclidMat_sd
    
    #calculte rho for each data point
    if metric == "COM":
       d0=euclidMat_avg/2
    elif metric == "TORSION":
       d0=euclidMat_avg-2*euclidMat_sd
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

    Nclust=0
    delta_avg=numpy.mean(delta)
    delta_sd=numpy.std(delta)
    if metric == "COM":
       dnorm=0.999
    elif metric == "TORSION":
       dnorm=0.999

    for deltai in delta:
        delta_norm=scipy.stats.norm(loc=delta_avg,scale=delta_sd).cdf(deltai)
        if delta_norm > dnorm:
           Nclust+=1

    # Select number of clusters by visual inspection and assign centers
    #Nclust=int(raw_input("Plot the rhodelta file and tell me how many cluster centers you observe:\n"))

    
   
    deltaCenters=sorted(delta,reverse=True)[0:Nclust]

    clusterCenters=[]
    for i in range(0,len(values)):
        if delta[i] in deltaCenters:
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
    for center in clusters.keys():
        fraction=float(len(clusters[center]))/float(len(values))
        print fraction
        if fraction >= 0.05:
           print "Cluster with center at ", center, "picoseconds has ", len(clusters[center]), "elements. --SAVED"
           clusters_forward[center]=clusters[center]
        else:
           print "Cluster with center at ", center, "picoseconds has ", len(clusters[center]), "elements. --DISCARDED"
        
        totalElements=totalElements+len(clusters[center])
    clusters=clusters_forward
    Nclust=len(clusters.keys())
    print "Total elements clustered:", totalElements, ". Number of observations: ",len(values)," (MUST BE THE SAME)"
    print "Number of clusters: ", Nclust
    return clusters

def genRefs(gro,trr,apolar,polar,clusters):

    refs=[]
    name=trr.split('.')[0]

    traj=mdtraj.load_trr(trr, top=gro)
    time=traj.time.tolist()

    #select atoms in the bsite
    reflist=[]
    for pdb in [apolar,polar]:
       filein=open(pdb,'r')
       for line in filein:
           if line.startswith('ATOM') or line.startswith('HETATM'):
              reflist.append(int(line.split()[1])-1)
       filein.close()
    reflist=sorted(reflist)
        
        
    for key in clusters.keys():
        outfile='ref_'+name+'_'+str(int(key))+'.pdb'
        snapfile='snap_'+name+'_'+str(int(key))+'.pdb'
        framenumber=time.index(key)
        coords=traj.xyz[framenumber]
        top=traj.topology
        frame=mdtraj.Trajectory(coords,top)
        frame.save_pdb(snapfile)
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
        filein.close()
        fileout.close()
        os.system('rm temp.pdb')
        refs.append(outfile)
    return refs

def plumeDatGenerator(dat,refs,rmsd,kappa):
    iteration=refs[0].split('_')[1] 
    labels=[]
    nameOut=iteration+'.dat'
    
    fileout=open(nameOut,'w')
    for namepdb in refs:
        time_str=namepdb.split('_')[2].split('.')[0]
        label='rmsd_'+iteration+'_'+time_str
        labels.append(label)
        line=label+': RMSD REFERENCE='+namepdb+' TYPE=OPTIMAL\n LOWER_WALLS ARG='+label+' AT='+rmsd+' KAPPA='+kappa+'\n'
        fileout.write(line)
    fileout.close()

    filein=open(dat,'r')
    fileout=open('temp.dat','w')
    for line in filein:
        if line.startswith('PRINT'):
           lineprint=line
        else:
           fileout.write(line)
 
    include_str='INCLUDE FILE='+nameOut+'\n'
    fileout.write(include_str)

    newLinePrint=''
    for element in lineprint.split():
        if element=='PRINT':
           newLinePrint=newLinePrint+element
        elif element.startswith('STRIDE') or element.startswith('FILE'):
           newLinePrint=newLinePrint+' '+element
        elif element.startswith('ARG'):
           for label in labels:
               element=element+','+label
           newLinePrint=newLinePrint+' '+element
    fileout.write(newLinePrint)
    fileout.close()
    command = 'mv temp.dat '+dat
    os.system(command)

   
            




if __name__=='__main__':
    
    apolar,polar,structure,output,stride,trr,gro,dat,rmsd,kappa,metric = parse(parser)
    
    residues = findResidues(apolar,polar)
    
    atomGroups = findAtomGroups(structure,residues)

    genPlumedinput(atomGroups,output,stride)

    labels,time,values=Plumed(trr,output,stride)

    clusters=Clustering(labels,time,values)
    
    refs=genRefs(gro,trr,apolar,polar,clusters)
    
    plumeDatGenerator(dat,refs,rmsd,kappa)

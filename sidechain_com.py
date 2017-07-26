import argparse,numpy,os,mdtraj

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


def parse(parser):
    args = parser.parse_args()
    return args.apolar,args.polar,args.input,args.output,args.stride,args.trr,args.gro



def findResidues(apolar,polar): # Change this function to do it with mdtraj
    residues=[]
    for pdb in [apolar,polar]:
        filein=open(pdb,'r')
        for line in filein:
            if line.startswith("ATOM") or line.startswith ("HETATM"):
               line=line.split()
               if line[5] not in residues:
                  residues.append(line[5])
        filein.close()
    print "There are", len(residues), "residues" 
    return residues

def findSideChains(structure,residues): #change this function to do it with mdtraj
    filein=open(structure,'r')
    backbone=['CA','C','O','N']
    SOLIONS=['SOL','WAT','HOH','NA','CL']
    sideChains={}
    currentResidue='NONE'
    for line in filein:
        if line.startswith("ATOM") or line.startswith ("HETATM"):
           line=line.split()
           if line[3] in SOLIONS:
              break
           if line[5] not in residues:
              continue
           else:
               if line[5]!=currentResidue:
                  currentResidue=line[5]
                  sideChains[currentResidue]=[]
               if line[2] not in backbone: # Add something to keep hydrogens out
                  sideChains[currentResidue].append(line[1])
    filein.close()
    return sideChains

def genPlumedinput(sideChains,output,stride):
    fileout=open(output,'w')
    labels_com=[]
    for key in sideChains.keys():
        label='com'+key
        line=label+': COM ATOMS='+','.join(sideChains[key])+'\n'
        fileout.write(line)
        labels_com.append(label)
    
    labels_dist=[]
    for i in range(0,len(labels_com)):
        for j in range(i,len(labels_com)):
            if j!=i:
               label='dist'+labels_com[i]+labels_com[j]
               line=label+': DISTANCE ATOMS='+labels_com[i]+','+labels_com[j]+'\n'
               fileout.write(line)
               labels_dist.append(label)
    print "there are ",len(labels_dist), "distances to print"

    line='PRINT ARG='+','.join(labels_dist)+' STRIDE=1 FILE=COMCOLVAR'
    fileout.write(line)    
    
    fileout.close()

def Plumed(trr,output,stride):
    command='plumed driver --mf_trr '+trr+' --plumed '+output
    os.system(command)
   
    # Get values from PLUMED output
    time=[]
    values=[]
    filein=open('COMCOLVAR','r')
    for line in filein:
        line=line.split()
        if line[0]=='#!':
           labels=line[3:]
        else:
           time.append(float(line[0])*float(stride))
           values.append(map(float,line[1:]))
    values=numpy.array(values)
    return labels, time, values



def Clustering(labels,time,values): # Done as in Rodriguez & Laio, Science(2014),344,6191,1492-1496 :
    
   #calculate euclidean distances:
    euclidMat=[]
    euclidMat_avg=[]
    for i in range (0,len(values)):
        euclid=[]
        for j in range(0,len(values)):
            if i==j:
               euclid.append(99999999.)
            elif j>i:
               euclid.append(numpy.linalg.norm(values[i]-values[j]))
               euclidMat_avg.append(numpy.linalg.norm(values[i]-values[j]))
            elif j<i:
               euclid.append(99999999.)
        euclidMat.append(euclid)
    euclidMat=numpy.array(euclidMat)
    #print euclidMat
    euclidMat_avg=numpy.mean(euclidMat_avg)
    
    #calculte rho for each data point
    d0=euclidMat_avg/2 #This parameter needs to be optimised for each system and CV  
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


    # Select number of clusters by visual inspection and assign centers
    Nclust=int(raw_input("Plot the rhodelta file and tell me how many cluster centers you observe:\n"))
   
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
    for center in clusters.keys():
        print "cluster with center at ", center, "picoseconds has ", len(clusters[center]), "elements"
        totalElements=totalElements+len(clusters[center])
    print "Total elements clustered:", totalElements, ". Number of observations: ",len(values)," (MUST BE THE SAME)"
   
    return clusters

def genRefs(gro,trr,apolar,polar,clusters):

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
        framenumber=time.index(key)
        coords=traj.xyz[framenumber]
        top=traj.topology
        frame=mdtraj.Trajectory(coords,top)
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
            




if __name__=='__main__':
    
    apolar,polar,structure,output,stride,trr,gro = parse(parser)
    
    residues = findResidues(apolar,polar)
    
    sideChains = findSideChains(structure,residues)

    genPlumedinput(sideChains,output,stride)

    labels,time,values=Plumed(trr,output,stride)

    clusters=Clustering(labels,time,values)
    
    genRefs(gro,trr,apolar,polar,clusters)
    
    

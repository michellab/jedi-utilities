import argparse

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
                    default='apolar.pdb')
parser.add_argument('-o','--output', nargs="?",
                    help='Name of the input file to supply to PLUMED',
                    default='com.dat')
parser.add_argument('-s','--stride', nargs="?",
                    help='PLUMED printing stride',
                    default='1')


def parse(parser):
    args = parser.parse_args()
    return args.apolar,args.polar,args.input,args.output,args.stride



def findResidues(apolar,polar):
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

def findSideChains(structure,residues):
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

    line='PRINT ARG='+','.join(labels_dist)+' STRIDE='+stride+' FILE=COMCOLVAR'
    fileout.write(line)    
    
    fileout.close()
        

if __name__=='__main__':
    
    apolar,polar,structure,output,stride = parse(parser)
    
    residues = findResidues(apolar,polar)
    
    sideChains = findSideChains(structure,residues)

    genPlumedinput(sideChains,output,stride)

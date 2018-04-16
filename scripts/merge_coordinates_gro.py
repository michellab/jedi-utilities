import glob,os,sys

#protein_in=sys.argv[1]
#complex_out=sys.argv[-1]

prot_in=open(sys.argv[1],'r') # Protein input *.gro file
complex_out=open(sys.argv[-1],'w') # Complex output *.gro file


prot_tmp=open('prot.tmp','w')
water_tmp=open('water.tmp','w')

i=0
for line in prot_in:
    i=i+1
    if i==1:
       complex_out.write(line) # write *.gro header
    elif i==2:
       nAtoms=int(line.split()[0]) # store the initial number of atoms in a veriable
    #elif 'HOH' not in line and len(line.split())!=3: # Not necessary if we delete the last line is previously removed
    elif 'HOH' not in line: # 
       prot_tmp.write(line)
    else:
       water_tmp.write(line)

prot_in.close()
prot_tmp.close()
water_tmp.close()

ligTempList=[]
for index in range(2,len(sys.argv)):
    i=0
    ligand=sys.argv[index]
    lig_in=open(ligand,'r')
    ligTempName=ligand+'.tmp'
    ligTempList.append(ligTempName)
    lig_tmp=open(ligTempName,'w')
    for line in lig_in:
        i=i+1
        if i==2:
           nAtoms=nAtoms+int(line.split()[0]) # update number of atoms
        elif i!=1 and len(line.split())!=3:
           lig_tmp.write(line)
    lig_in.close()

nAtoms=' '+str(nAtoms)+'\n'
complex_out.write(nAtoms) # write updated number of atoms

prot_tmp=open('prot.tmp','r')
for line in prot_tmp:
    complex_out.write(line) # write protein coordinates
prot_tmp.close()
os.system('rm prot.tmp')

water_tmp=open('water.tmp','r')
for line in water_tmp:
    complex_out.write(line) # write water coordinates
water_tmp.close()
os.system('rm water.tmp')

for ligand in ligTempList:
    lig_tmp=open(ligand,'r')
    for line in lig_tmp:
        complex_out.write(line) # write ligand coordinates
    lig_tmp.close()
    rmCommand='rm '+ligand
    os.system(rmCommand)
'''
water_tmp=open('water.tmp','r')
for line in water_tmp:
    complex_out.write(line) # write water coordinates
water_tmp.close()
os.system('rm water.tmp')
'''

complex_out.close()

print "Coordinates merged"


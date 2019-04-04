import os,glob,sys,random

def atomnames(oldname,newname,replacement,used_news): # check if the atomnames are repeated. Assuming that 'name' and 'bondtype' are always the same
    #print "Entering atomnames functions"
    if newname in replacement.keys() or newname in replacement.values():
       newname, used_news=gnn(used_news)
       atomnames(oldname,newname,replacement,used_news) # this will generate a "recursive loop" that will keep verifying that newname is not in atomlist
    else:
       replacement[oldname]=newname
    return replacement, used_news

def linereplace(line,oldname,replacement):
    line=line.replace(' '+oldname+' ',' '+replacement[oldname]+' ') # replace the old atom name by the new one
    return line
def gnn(used_news):
    number='0123456789'
    alpha='abcdefghijklmnopqrstuvwxyz'
    newname=random.choice(alpha)+random.choice(number)
    while newname in used_news:
      newname=random.choice(alpha)+random.choice(number)
    used_news.append(newname)
    return newname, used_news

ligands_atomtypes=open('LigandAtomTypes.itp','w')
ligands_atomtypes.write('[ atomtypes ]\n')

replacement={}
used_news=[]
for index in range(1,len(sys.argv)):
    ligand=sys.argv[index]
    comment='; Atom types for file '+ligand+'\n'
    ligands_atomtypes.write(comment)
    filein=open(ligand,'r')
    ligand_tmp=open(ligand+'.tmp','w')
    switch_atomtypes=0
    switch_atoms=0
    for line in filein:
        if line=='\n':
            #print 'washere nothing'
            switch_atomtypes=0
            switch_atoms=0
            ligand_tmp.write(line)
        
        elif '[ atomtypes ]' in line:
            #print 'washere atomtypes'
            switch_atomtypes=1
        elif switch_atomtypes==1 and ';' not in line.split()[0]:
        #elif switch_atoms==1 and ';' not in line.split()[0]:
            #print replacement 
            oldname=line.split()[0]
            newname=line.split()[1]
            replacement,used_news=atomnames(oldname,newname,replacement,used_news)
            line=linereplace(line,oldname,replacement)
            ligands_atomtypes.write(line)
  
        elif '[ atoms ]' in line:
             #print 'washere atom'
             switch_atoms=1
             ligand_tmp.write(line)  
        elif switch_atoms==1 and line.split()[0]!=';':
             #print 'washere switch2'
             oldname=line.split()[1]
             line=linereplace(line,oldname,replacement)
             ligand_tmp.write(line)

        else:
             ligand_tmp.write(line)


    filein.close()
    ligand_tmp.close()
    command='mv '+ligand+'.tmp '+sys.argv[index]        
    os.system(command)
ligands_atomtypes.close()


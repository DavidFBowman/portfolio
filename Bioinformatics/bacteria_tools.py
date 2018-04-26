#Bacteria Tools for use in week 3 of course 2 of Bioinformatics


def make_RNA_table():
    #converts codon to single letter peptide. The stop codon is represented as x
    fid='C:\Users\David\Dropbox\Projects\Courses\Bioinformatics\GenomeSequencing\W3\RNA_codon_table_1.txt'
    f=open(fid,'r')
    lines=f.readlines()
    rna={}
    for i in lines:
        line=i.split()
        if len(line)==2:
            rna.update({line[0]:line[1]})
        else:
            rna.update({line[0]:'X'})
    return rna
    
def make_pep_table():
    #dict giving peptide letters with corresponding codons
    fid='C:\Users\David\Dropbox\Projects\Courses\Bioinformatics\GenomeSequencing\W3\RNA_codon_table_1.txt'
    f=open(fid,'r')
    lines=f.readlines()
    rna=[]
    pep=[]
    pep_dict={}
    for i in lines:
        line=i.split()
        if len(line)==2:
            rna.append(line[0])
            pep.append(line[1])
        else:
            rna.append(line[0])
            pep.append('X')
    uniques=set(pep)
    for i in uniques:
        a=[j for j, x in enumerate(pep) if x == i]
        b=[]
        for k in a: b.append(rna[k])
        pep_dict.update({i:b})
    return pep_dict
    
def DNA_to_RNA(string):
    rna=''
    for i in range(0,len(string)):
        if string[i]=='T':
            rna+='U'
        else:
            rna+=string[i]
    return rna

def RNA_to_DNA(string):
    rna=''
    for i in range(0,len(string)):
        if string[i]=='U':
            rna+='T'
        else:
            rna+=string[i]
    return rna
    
def reverse_complement(dna):
    dna=dna[::-1]
    dnaout=''
    for i in range(0,len(dna)):
        if dna[i]=='A':
            dnaout+='T'
        if dna[i]=='T':
            dnaout+='A'
        if dna[i]=='C':
            dnaout+='G'
        if dna[i]=='G':
            dnaout+='C'
    return dnaout
    
def all_peps(rna,codon_to_pep):
    #given a string of rna with len%3==0, will return each pep from the sequencewise 3mers
    peps=''
    for i in range(0,len(rna)/3):
        peps+=codon_to_pep[rna[i*3:(i+1)*3]]
    return peps
        
       
def encode_peptide(dna,pep,codon_to_pep,pep_to_codon):
    #given a string of rna and a string of peps will return all possible instances of the pep within the rna
    dna2=reverse_complement(dna)
    rna=DNA_to_RNA(dna)
    rna2=DNA_to_RNA(dna2)
    matches=[]
    for i in range(0,len(rna)-(len(pep)*3)):
        check=all_peps(rna[i:i+(len(pep)*3)],codon_to_pep)
        check2=all_peps(rna2[i:i+(len(pep)*3)],codon_to_pep)
        if check==pep:
            matches.append(RNA_to_DNA(rna[i:i+(len(pep)*3)]))
        if check2==pep:
            matches.append(reverse_complement(RNA_to_DNA(rna2[i:i+(len(pep)*3)])))            
    return matches
    
def peptide_mass(peps,mass_table):
    masses=[0]
    for pep in peps:
        mass=0
        for i in pep:
            mass+=mass_table[i]
        masses.append(mass)
    return sorted(masses)
    
def cyclic_peptide_masses(string):
    #given a peptide string, return masses of all cyclic permutations through length 1:len(string)
    fid='C:\Users\David\Dropbox\Projects\Courses\Bioinformatics\GenomeSequencing\W3\peptide_masses.txt'
    mass_table = {}
    with open(fid,'r') as f:
        for line in f:
            (key, val) = line.split()
            mass_table[key] = int(val)
    doublestring=string+string
    perms=[string]
    for i in range(1,1+len(doublestring)/2-1):
        for j in range(0,len(doublestring)/2):
            perms.append(doublestring[j:j+i])
    masses=peptide_mass(perms,mass_table)            
    return masses
    
def linear_peptide_masses(string):
    #given a peptide string, return masses of all cyclic permutations through length 1:len(string)
    #this creates way too many strings, don't use it for anything massive without fixing it first
    fid='C:\Users\David\Dropbox\Projects\Courses\Bioinformatics\GenomeSequencing\W3\peptide_masses.txt'
    mass_table = {}
    with open(fid,'r') as f:
        for line in f:
            (key, val) = line.split()
            mass_table[key] = int(val)
    perms=[]
    for i in range(1,1+len(string)):
        for j in range(0,len(string)):
            perms.append(string[j:j+i])
    masses=peptide_mass(set(perms),mass_table) 
    return masses            
            

def RNA_to_amino(string,codon_to_pep):
    amino=''
    for i in range(0,len(string)/3):
        amino+=codon_to_pep[string[i*3:(i+1)*3]]

    return amino
            
        
    
        


codon_to_pep=make_RNA_table()
pep_to_codon=make_pep_table()

#outstring=''
#for i in range(0,len(RNA_string)/3):
#    rna_r=RNA_string[(i*3):(i*3)+3]
#    outstring+=(rna[rna_r])

    
    
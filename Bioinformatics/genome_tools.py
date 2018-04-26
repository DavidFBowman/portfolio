""" Simple tools for reconstructing a genome.
"""


def kmer_composition(full,k):
    """ Return the kmer composition of a string.
    i.e all of the constituent components of length k.
    """
    kmers=[]
    for i in range(0,len(full)-k+1):
        kmers.append(full[i:i+k])
    return sorted(kmers)
    
def write_kmers(kmers,fid):
    """ Join all of a list of kmers into a single string and write to file.
    Use this for testing on the kmer composition of a known genome sequence.
    """
    #fid='C:\Users\David\Dropbox\Projects\Courses\Bioinformatics\GenomeSequencing\W1\data.txt''
    txt = '\n'.join(kmers)
    f = open(fid,'w')
    f.write(txt)
    f.close()
    
def read_kmers():
    """ Return all of the kmers found in a file. """
    fid='C:\Users\David\Dropbox\Projects\Courses\Bioinformatics\GenomeSequencing\W1\kmers.txt'
    f=open(fid,'r')
    lines=f.readlines()
    kmers=[]
    for i in range(0,(len(lines)/2)):
        kmers.append(lines[i*2][:-1])
    #kmers.append('TGCA')
    return kmers
        
    
def remove_duplicate(seq):
    """ Remove duplicate elements."""
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]
    
def duplicate_index(l):
    """Returns index of duplicate 1-mers in string
    """
    d=[]
    for i in range(0,len(l)):
        for j in range(0,len(l)):
            if l[i]==l[j] and i!=j:
                d.append(j)
    j=0
    for i in range(0,len(d)/2):
        del d[(i*2)-j]
        j+=1
    return d
                
    
    
def de_brujin_graph(kmers,k):
    """ Returns k-1 edge overlap kmers as a list of [k-1 prefixes,[valid k-1 suffixes]].
    Takes a string, converts to kmers, finds all overlaps for each prefix and then merges the duplicates.
    """
    
    #kmers=kmer_composition(string,k)
    suffix=[]
    prefix=[]
    brujin=[]
    for i in kmers:
        suffix.append(i[1:])
        prefix.append(i[:-1])
    for i in range(0,len(prefix)):
        match=[]
        for j in range(0,len(prefix)):
            if prefix[i]==prefix[j]:
                match.append(suffix[j])
        brujin.append(match)
    d=duplicate_index(prefix)
    prefix[:] = [ item for i,item in enumerate(prefix) if i not in d ]
    brujin[:] = [ item for i,item in enumerate(brujin) if i not in d ]
    return prefix,brujin
    
def number_edges(into,outof):
    """ Converts the values of outof to the corresponding indices of into.
    """
    outof_num=[]
    into_num=[]
    unvisited=[]
    for i in outof:
        a=[]
        b=[]
        for j in i:
            try:
                a.append(into.index((j)))
                b.append(1)
            except ValueError:
                a.append(False)
                b.append(0)
        outof_num.append(a)
        unvisited.append(b)
    for i in range(0,len(into)):
        into_num.append(i)
    return into_num,outof_num,unvisited
    
def edge_counter(outof):
    out_count=[]
    for i in outof:
        out_count.append(len(i))
    return out_count
    
        
def walk(into,outof,edges,unvisited,visited):
    """Performs a walk through the kmers attempting to combine univisted combinations."""
    valid_route=True
    pos=0
    cycle=[]
    while valid_route==True:
        print pos
        cycle.append(pos)
        valid_route=False
        for i in range(0, len(outof[pos])):
            if visited[pos][i]==1:
                valid_route=True
                #print outof[pos]
                unvisited[pos][i]=0
                pos=outof[pos][i]
                
                break
            else:
                print 'no valid route'
                print pos
                print outof[pos]
                print i
                print visited[pos]
    return cycle,unvisited
        

def reconstrcut_genome():
    """Returns indices of valid walk through kmers to form a genome."""
    kmers=read_kmers()
    into,outof=de_brujin_graph(sorted(kmers),4)
    into_num,outof_num,visited=number_edges(into,outof)
    cycle,visited=walk(into_num,outof_num,edge_counter,visited)

 
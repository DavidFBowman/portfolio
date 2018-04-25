"""Uses a Monte Carlo balls algorithm to create a cell of clustered vacancies.

These can be read into any other calculation to set which atomic sites are occupied
and which are vacant.

The aim of this is usually to check affect or correlations on higher order calculations.

"""

import numpy as np
import random
import matplotlib.pyplot as plt
import time
import os


outfile=(os.path.dirname(__file__)+'/'+"pos_correlations_32_0p036_10loop.dat")
NN_dist_file=os.path.dirname(__file__)+'/'+"NN_array_32cubed_dist.dat"
NN_index_file=os.path.dirname(__file__)+'/'+"NN_array_32cubed_index.dat"

path='C:\\Users\\David\\Dropbox\\PhD\\BallsAndSprings\\Correlate\\Plots\\'

vac=0.0357 # Percentage of vacancy sites in system
T=300   # Starting temperature
k1=1e-1 # Scaling parameter used as a product with the temperature to control convergence
max_nn=1.5 # Distance within which a pair of atoms are considered nearest neighbours
loops=10 # How many iterations of the algorithm should be run per temperature

# The correlation variable (along with K) determines the nature of the correlations
# The first value gives the energy scaling for an occupied site and the second value for a vacancy
# [5e1,1] with k=0.1 gives small positive clusters
# [1, 5e1] with k=10 gives weak negative clustering but no super structure forming

correlation=[5e1,1] 

random.seed(312312)
nx,ny,nz=8,8,8 # Number of unit cells
crys_coord=np.zeros((48,3))
crys=np.zeros((48*nx*ny*nz,3))

#Setting oxygen(2) atomic positions for pyrochlore, u variables can be obtained from the
#literature as the displacement paramters

u=0.20485
uu=0.04516
uuu=0.29515

crys_coord[0,:]=[0,1-u,0]
crys_coord[1,:]=[0,u,0]
crys_coord[2,:]=[u,0,0]
crys_coord[3,:]=[0.25,0.25,uu]
crys_coord[4,:]=[uuu,0.5,0]
crys_coord[5,:]=[0.5,uuu,0]
crys_coord[6,:]=[0.5,1-uuu,0]
crys_coord[7,:]=[1-uuu,0.5,0]

crys_coord[8,:]=[1-u,0,0 ]
crys_coord[9,:]=[0.75,0.75,uu ]
crys_coord[10,:]=[0,0,u ]
crys_coord[11,:]=[0,0.5,uuu  ]
crys_coord[12,:]=[uu,0.25,0.25]
crys_coord[13,:]=[0.25,uu,0.25]
crys_coord[14,:]=[0.25,0.5-uu,0.25 ]
crys_coord[15,:]=[0.5,0.5,u]

crys_coord[16,:]=[0.5-uu,0.25,0.25]
crys_coord[17,:]=[0.5,0,uuu]
crys_coord[18,:]=[0.5+uu,0.75,0.25]
crys_coord[19,:]=[0.75,0.5+uu,0.25]
crys_coord[20,:]=[0.75,1-uu,0.25]
crys_coord[21,:]=[1-uu,0.75,0.25]
crys_coord[22,:]=[0,uuu,0.5]
crys_coord[23,:]=[0,1-uuu,0.5]

crys_coord[24,:]=[u,0.5,0.5]
crys_coord[25,:]=[0.25,0.25,1-uu]
crys_coord[26,:]=[uuu,0,0.5]
crys_coord[27,:]=[0.25,0.75,0.5+uu]
crys_coord[28,:]=[0.5,u,0.5]
crys_coord[29,:]=[0.5,1-u,0.5]
crys_coord[30,:]=[0.5+u,0,0.5]
crys_coord[31,:]=[0.75,0.75,0.5-uu]

crys_coord[32,:]=[0.75,0.25,0.5+uu]
crys_coord[33,:]=[1-u,0.5,0.5]
crys_coord[34,:]=[0,0.5,1-uuu]
crys_coord[35,:]=[0,0,1-u]
crys_coord[36,:]=[uu,0.75,0.75]
crys_coord[37,:]=[0.25,0.5+uu,0.75]
crys_coord[38,:]=[0.25,1-uu,0.75]
crys_coord[39,:]=[0.5,0,1-uuu]

crys_coord[40,:]=[0.5-uu,0.75,0.75]
crys_coord[41,:]=[0.5,0.5,1-u]
crys_coord[42,:]=[0.5+uu,0.25,0.75]
crys_coord[43,:]=[0.75,uu,0.75]
crys_coord[44,:]=[0.75,0.5-uu,0.75]
crys_coord[45,:]=[1-uu,0.25,0.75]
crys_coord[46,:]=[0.25,0.75,1-uu]
crys_coord[47,:]=[0.75,0.25,1-uu]

xtal_vac=np.ones(len(crys)).astype(int)

    
def energy2(NN,dist,vac_list,correlation=correlation):
    """Returns the cumulative nearest neighbour bond energy for the system."""
    cumul_E=0
    for i in range(0,len(NN[0])):
        cumul_E+=correlation[vac_list[i]]/(dist[i])
    return cumul_E
    
def metro(E1,E2,T,k=k1):
    """Runs the metropolis algorithm with success probability scaled by temperature"""
    accept=np.exp( (E2-E1)/(k*T))
    if (accept>(random.random())):
        return 1
    else:
        return 0
    
def distance(ptA,ptB):
    """Calculates the displacement between two 3D vectors"""
    dist=ptA-ptB
    return(dist[0]**2+dist[1]**2+dist[2]**2)
    
def writeNN(NN_dist, NN_index,NN_dist_file, NN_index_file):
    """Write the nearest neighbour data to text files"""
    fdist=open(NN_dist_file,'w')
    findex=open(NN_index_file,'w')
    for i in range(0,len(NN_dist)):
        for j in range(0,len(NN_dist[i])):
            fdist.write(str(NN_dist[i][j]))
            findex.write(str(NN_index[i][0][j]))
            fdist.write(' ')
            findex.write(' ')
        fdist.write('\n')
        findex.write('\n')
    fdist.close()
    findex.close()

def write_correlations(xtal_vac,outfile):
    """Write vacancy status of an atom to file for reading into the balls and springs model"""
    fcorr=open(outfile,'w')
    for i in range(0,len(xtal_vac)):
        fcorr.write(str(xtal_vac[i]))
    fcorr.close()
    
def find_NN(crys, max_nn=max_nn):
    """Returns the max_nn number of nearest neighbour atoms for each atom in the crystal"""
    NN_dist=[]
    NN_index=[]
    for i in range(0,len(crys)):
        NN_i=np.sum(np.abs(crys[i]-crys )**2,axis=-1)**(1./2)
        index=np.where((NN_i<max_nn) & (NN_i>0.01))
        dist=NN_i[index]
        NN_dist.append(dist)
        NN_index.append(index)
    return NN_dist,NN_index        

def read_NN(str1,str2):
    """Read nearest neighbour distances from a file to save recalcuating it each 
    time the code is run. This format is dictated by the writeNN method.
    """
    NN_dist=np.load(str1)
    NN_index=np.load(str2)   
    return NN_dist, NN_index  
    
        
# Building xtal
cell=-1            
for i in range(0,nx):
    for j in range(0,ny):
        for k in range(0,nz):
            cell+=1
            for a in range(0,len(crys_coord)):
                site=(cell*len(crys_coord))+a
                crys[site,0]=crys_coord[a,0]+i
                crys[site,1]=crys_coord[a,1]+j 
                crys[site,2]=crys_coord[a,2]+k
start_time=time.time()               
# Setting nearest neighbours
NN=np.zeros((len(crys)*2,max_nn))
NN.fill(100)                
 
# Run the first below line to calculate nearest neighbours
# Run the second alternatively to read this data from a file

NN_dist,NN_index=find_NN(crys)             
#NN_dist,NN_index=read_NN(os.path.dirname(__file__)+'/NN_array_32cubed_dist.npy',os.path.dirname(__file__)+'/NN_array_32cubed_index.npy')  

# Add vacancies to system            
vac_count=0
max_vac=len(xtal_vac)*vac
for j in range(0,100):
    for i in range(0,len(crys)):
        if (vac_count>max_vac):
            break
        if (random.random()<vac):
            xtal_vac[i]=0
            vac_count+=1

# Set the desired temperatures to run the algorithm on
# If you have issues with convergence try using a logarithmic descent

temp_list=np.linspace(T,1,T)
#temp_list=np.logspace(10,1e-6,500,base=1.1)

accept_rate=np.zeros([len(temp_list)])

# Run the Simulated Annealing algorithm
T_index=-1            
for j in temp_list:
    T_index+=1
    for k in range(0,loops):
        count=0        
        for i in range(0,len(crys)):
            v_site=int(random.random()*len(crys))
            # Pick a vacancy
            if (xtal_vac[v_site]==0):
                new_site=NN_index[v_site][0][int(random.random()*len(NN_index[v_site][0]))]
                # Is new neighbour a valid site for a swap?
                if (xtal_vac[new_site]==1):
                    E1=energy2(NN_index[v_site],NN_dist[v_site],xtal_vac[NN_index[v_site]])
                    E2=energy2(NN_index[new_site],NN_dist[new_site],xtal_vac[NN_index[new_site]])
                    accept=metro(E1,E2,j)
                    count+=1
                    if (accept==1):
                        xtal_vac[new_site]=0
                        xtal_vac[v_site]=1
                        accept_rate[T_index]+=1
        accept_rate[T_index]=accept_rate[T_index]/count                 
print("--- %s seconds ---" % (time.time() - start_time))

write_correlations(xtal_vac,outfile)  

# Plot a 3D map of the clustering. Only run this for a small test cell.                                
fig=plt.figure()
ax=fig.add_subplot(111, projection='3d')
for i in range(0,len(crys)):
    if (xtal_vac[i]==0):
        ax.scatter(crys[i,0], crys[i,1], crys[i,2],  c='r', marker='o')
    if (xtal_vac[i]==1):
        ax.scatter(crys[i,0], crys[i,1], crys[i,2],  c='b', marker='o',s=0.1)
plt.show()


tstr= 'size'+str(nx)+'x'+str(ny)+' pars='+str(correlation)+' k='+str(k1) \
      +' iter='+str(len(temp_list))+'x'+str(loops)+' nn_dist='+str(max_nn)+' vacpercent='+str(vac)
#plt.title(tstr)
tstr=tstr.replace('.','p')

tstr=tstr.replace(' ','_')
outfile=('C:\\Users\\David\\Dropbox\\PhD\\BallsAndSprings\\Correlate\\ClusterSizeVacs'+'/'+tstr+"possmall.dat")
write_correlations(xtal_vac,outfile)
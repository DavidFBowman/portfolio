from scipy.signal import convolve
import scipy.interpolate as inter
import numpy as np
import matplotlib.pyplot as plt

class Phonon:

	"""Phonon class for importing, manipulating and plotting data from .castep files.
	
	For a quick result just run initialise_from_file with the calculated data file path as an argument and it will run the below functions
	
	Having initialised the object, data can be generated using the grabData function which will read from a file of path fid
	Experimental data from MERLIN processed via HORACE of form x,y,E can be loaded in via loadExpData
	Make allQ will work out the range of data to be considered given the above input files are loaded
	fill_intensities will calcualte the intensities as S(Q,omega) and trim_intensities will set all values that are no present in the experimental
	data to zero for easy comparison.
	
	
	"""

    def __init__(self):
		"""Sets data types for initial parameters"""
        self.unitCell=[]
        self.atPos=[]
        self.atType=[]
        self.atMass=[]
        self.hkl=[]
        self.branchEnergy=[]
        self.eigens=[]
        self.formParams=[]
        self.noOfEigFiles=0
		
		
	def extractPhonons(fileName):
		"""Populate phonon parameters for analysis"""
		if fileName != "None" or "":
			file1 = open(fileName,'r')
			noOfIons, noOfBranches, noOfWavVec, units, unitCell, atPos, atType, atMass = readHeader(file1)

			hkl, branchEnergy, eigens = readEigens(noOfIons, noOfBranches, noOfWavVec, file1)
			
			formParams = []
			for x in atType:
				formParams.append(diffplotmod.readFormFactorParameters(x))
        
    return unitCell, atPos, atType,atMass, hkl, branchEnergy, eigens, formParams
        
    def grab_data(self,fid):
        """takes data from .phonon CASTEP file"""
        self.calcFid=fid      
        unitCellTemp, atPosTemp, atTypeTemp,atMassTemp, hklTemp, branchEnergyTemp, eigensTemp, formParamsTemp =  extractPhonons(self.calcFid)
        branchEnergyTemp = branchEnergyTemp/8.06554            
        self.unitCell.append(unitCellTemp)
        self.unitCell=np.squeeze(self.unitCell)
        self.atPos.append(atPosTemp)
        self.atPos=np.squeeze(self.atPos)
        self.atType.append(atTypeTemp)
        self.atMass.append(atMassTemp)
        self.atMass=np.squeeze(self.atMass)
        self.hkl.append(hklTemp)
        self.hkl=np.squeeze(self.hkl)
        self.branchEnergy.append(branchEnergyTemp)
        self.branchEnergy=np.squeeze(self.branchEnergy)
        self.eigens.append(eigensTemp)
        self.eigens=np.squeeze(self.eigens)
        self.formParams.append(formParamsTemp)
        self.formParams=np.squeeze(self.formParams)
        self.noOfEigFiles += 1
        
    def loadExpData(self,fid0):
        """Takes data from experiment output from HORACE of type x y E""" 
        self.expFid=fid0
        data= np.loadtxt(self.expFid)
        datQ=data[:,0]
        datE=data[:,1]
        datInt=data[:,2]
        datE=datE[datQ >= 0]
        self.datE=datE
        datInt=datInt[datQ >= 0]
        self.datInt=datInt
        datQ=datQ[datQ >= 0]
        self.datQ=datQ
        self.entMaxE=max(datE)
        self.entMinE=-max(datE)
        self.QRangeActual=int(max(datQ))
        self.QRange=int(np.ceil(self.QRangeActual/2.))  

	def calcDiffPat(unitCell,atPos, atMass, formParams, eigens, branchEnergy, energyRange, hkl,temp, mix, width,elastic,scale):
    """Calculates the diffraction pattern using the equation for Q
	   also convolves the results with a voigt function to give an output consistent with the experimental resolution 
	"""
	#hBar = 1.05457148e-34
		def voigt(x,mix,back,sigma,height): #this is a pseudo voigt function, good for diffraction peaks
			def gauss(x,sigma):
				return exp(-((x/(0.600561*sigma))**2))
			def lor(x,sigma):
				return 1+((x)/(0.5*sigma))**2
			ans = height*(mix/lor(x,sigma) +((1-mix)*gauss(x,sigma)))+back
			return ans
		size = 5000    
		atNo = len(atPos)
		recipLattice = diffplotmod.RecipLatt(unitCell)
		Q = dot(recipLattice, hkl)
		xaxis = linspace(energyRange[0],energyRange[1], size)
		yaxis = zeros(size)
		formFactor = []
		for x in formParams:
			formFactor.append(x[9])
		for j in range(len(branchEnergy)):
			scatPow = complex(0.,0.)
			for i in range(atNo):
				working = complex(0.,0.)
				working += exp(1J*dot(Q,atPos[i,:]))
				working *= dot(Q,array([complex(eigens[j*atNo+i,0],-eigens[j*atNo+i,1]),
										complex(eigens[j*atNo+i,2],-eigens[j*atNo+i,3]),
										complex(eigens[j*atNo+i,4],-eigens[j*atNo+i,5])]))
				working *= formFactor[i]/sqrt(atMass[i])
				scatPow += working
			intensity = abs(scatPow * scatPow.conjugate())/abs(branchEnergy[j])
			if branchEnergy[j]<=energyRange[1] and branchEnergy[j] >=energyRange[0]:
				yaxis[(abs(xaxis-branchEnergy[j])).argmin()]+=real(intensity* boseFactor(branchEnergy[j], temp))
			if -branchEnergy[j]<=energyRange[1] and -branchEnergy[j] >=energyRange[0]:
				yaxis[(abs(xaxis+branchEnergy[j])).argmin()]+=real(intensity * boseFactor(-branchEnergy[j], temp))
		yaxis = yaxis*scale
		yaxis[(abs(xaxis-0.)).argmin()]+=elastic
		yVoigt = zeros(size)
		i = 0
		for x in linspace(energyRange[0],energyRange[1], size):
			yVoigt[i] = real(voigt(x, mix,0.,width, 1.))
			i+=1
			
		
		yaxis = convolve(yaxis,yVoigt, 'same')
            
    return xaxis, yaxis
	
	def findWavVec(meashkl, calchkl):
	"""Given a reciprocal space value find the nearest calculated values index"""

		def calcDist(meashkl, currhkl):
			dist = (meashkl[0]-currhkl[0])**2
			dist +=(meashkl[1]-currhkl[1])**2
			dist +=(meashkl[2]-currhkl[2])**2
			dist = sqrt(dist)
			return real(dist)
		calchkl=calchkl+100#this is so negative numbers don't cause discontinuities
		meashkl=meashkl+100
		meashkl = array([math.fmod(meashkl[0],1.),math.fmod(meashkl[1],1.),math.fmod(meashkl[2],1.)])
		shortestDist = calcDist(meashkl, calchkl[0,:])
		shortestIndex = 0
		
		for i in range(1,len(calchkl[:,0])):
			currhkl = array([math.fmod(calchkl[i,0],1.),math.fmod(calchkl[i,1],1.),math.fmod(calchkl[i,2],1.)])
			newDist = calcDist(meashkl, currhkl)
			if newDist < shortestDist:
				shortestDist = newDist
				shortestIndex = i
		#print calchkl-100
		#print calchkl[shortestIndex,:]-100
		return shortestIndex
     
    def calcIntensity(self,Kpt,minE=-25.0,maxE=25.0):
        """Returns calculated intensities within specified energy range at a given KPoint"""
        meashkl = Kpt
        temp = 30.0
        mix = 0.0
        width = 2.5
        energyRange = np.array([minE,maxE])
        elastic = 1.0
        idx = findWavVec(meashkl, self.hkl)
        xaxis, yaxis = calcDiffPat(self.unitCell,self.atPos, self.atMass, self.formParams, self.eigens[idx,:,:], self.branchEnergy[idx,:], energyRange, meashkl,temp, mix, width,elastic,1.)
        return yaxis
        
    def whichFile(self):
        """ Reads title of input file and returns values to help construct allQ
            1- nnH
            2- HHHmn
            3- nHH
        """
        slashPos=[pos for pos, char in enumerate(self.expFid) if char == '/'][-1] #find last instance of '/'
        title=self.expFid[slashPos+1:len(self.expFid)]
        print title
        self.title=title
        if 'oo' in title:
            self.kType=1
            self.kOffset=int(title[2])
            self.calcFid='Y2Ti2O7_00l.phonon'
        if 'HHH' in title:
            self.kType=2
            self.kOffset=int(title[4])
            self.calcFid='Y2Ti2O7_hhh.phonon'  
        if 'HH' in title and 'HHH' not in title and 'HH2H' not in title:
            self.kType=3
            self.kOffset=int(title[2])
            self.calcFid='Y2Ti2O7_hh0.phonon'
        if 'HH2H' in title:
            self.kType=4
            self.kOffset=int(title[5])
            self.calcFid='Y2Ti2O7_112.phonon' 
        else:
            print 'unknown file type'   
    
    def makeAllQ(self):
        """Generates all Q points based on input file type"""
        allQ=np.zeros(((len(self.hkl)-1)*self.QRange,3))
        QLen=len(self.hkl)-1
        if self.kType==1:
            for i in range(0,self.QRange):
                allQ[i*QLen:(i+1)*QLen,0]=self.kOffset
                allQ[i*QLen:(i+1)*QLen,1:3]=self.hkl[0:QLen,1:3]+i+(float(self.kOffset)/2.) 
        if self.kType==2:
            
            for i in range(0,self.QRange*2):
                allQ[i*QLen:(i+1)*QLen,:]=self.hkl[0:QLen,:]+i
                allQ[i*QLen:(i+1)*QLen,1:3]=allQ[i*QLen:(i+1)*QLen,1:3]-(self.kOffset/2.)         
        if self.kType==3:            
            for i in range(0,self.QRange):
                allQ[i*QLen:(i+1)*QLen,0]=self.hkl[0:QLen,0]+(2*i)
                allQ[i*QLen:(i+1)*QLen,1:3]=self.hkl[0:QLen,1:3]+i+(float(self.kOffset)/2.)
        if self.kType==4:
            for i in range(0,self.QRange):
                allQ[i*QLen:(i+1)*QLen,0]=self.hkl[0:QLen,0]+(2*i)
                allQ[i*QLen:(i+1)*QLen,1:3]=self.hkl[0:QLen,1:3]+(3*i)                    
        self.allQ=allQ
        
    def fill_intensities(self,energyBins=5000):
        """Creates array of calculated intensities, removes negative energy"""
        intensities=np.zeros((len(self.allQ),energyBins))
        j=0
        for i in self.allQ:
            intensities[j,:]=self.calcIntensity(i)
            j+=1
        intensities=np.transpose(intensities)
        self.intensities=intensities[energyBins/2:energyBins,:]
        
    def trim_intensities(self):
        """Reshapes intensities to same shape as data in calculated range"""
        intensities=self.intensities
        if (self.QRange*2 != self.QRangeActual):
        #if the QRange is odd this causes issues with the cell conversion that stretches the x-axis. This clause fixes that.    
            endPoint=int((float(self.QRangeActual)/float(self.QRange*2))*np.shape(self.intensities)[1])
            intensities=self.intensities[:,0:endPoint]
        
        ##make dummy array to trim sim array
        EArray=np.transpose(np.tile(np.linspace(0,25,2500),[np.shape(intensities)[1],1]))
        QArray=np.tile(np.linspace(0,self.QRange,np.shape(intensities)[1]),[2500,1])        
        ##remove sections outside of data range linearly (before getting to the NANs)
        trimIntensity=intensities[(EArray>=0) & (EArray<=max(self.datE)) & (QArray<=self.QRangeActual)]        
        self.cols=len(EArray[(EArray>=0) & (EArray<=max(self.datE))])/np.shape(EArray)[1]
        self.rows=len(trimIntensity)/self.cols
        colsj=np.complex(0,self.cols)
        rowsj=np.complex(0,self.rows)
        trimIntensity=np.reshape(trimIntensity,[self.cols,self.rows])
        
        ##smooth data
        datIntBlur=blur_image(np.reshape(self.datInt,(np.size(np.unique(self.datE)),np.size(np.unique(self.datQ)))),4)
        ##make nan mask from the exp data array (ravel is just conving to 1D)
        datMask=(np.ravel(datIntBlur)==np.ravel(datIntBlur))
        vals = datMask     
        #interpolation array
        pts = np.array([[i,j] for i in np.linspace(0,1,np.shape(datIntBlur)[0]) for j in np.linspace(0,1,np.shape(datIntBlur)[1])] )
        grid_x, grid_y = np.mgrid[0:1:colsj, 0:1:rowsj]
        grid_z = inter.griddata(pts, vals, (grid_x, grid_y), method='linear')        
        ##set areas outside exp range to 0s
        trimIntensity=np.multiply(trimIntensity,grid_z)
        self.trimIntensity=trimIntensity
        self.grid_z=grid_z
        self.dataConvolve=datIntBlur
        
    def plot_calculation(self,cutOff=10):
        """Plots calculated data"""
        plt.contourf(np.linspace(0,self.QRangeActual,600),np.linspace(0,int(max(self.datE)),2500),self.intensities,np.linspace(1e-16,cutOff,255))
        plt.show()
        
    def plot_data(self,cutOff=200):
        """plots full calculation range"""
        plt.contourf(np.linspace(0,max(self.datQ),np.shape(self.dataConvolve)[1]),np.linspace(0,int(max(self.datE)),np.shape(self.dataConvolve)[0]),self.dataConvolve,np.linspace(0,cutOff,255))
        plt.show()
        
    def plot_all(self,cutOffCalc=10,cutOffData=200,save=False):
        """plots calculation alongside MERLIN data"""
        f, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=True)
        ax1.contourf(np.linspace(0,self.QRangeActual,self.rows),np.linspace(0,int(max(self.datE)),self.cols),self.trimIntensity,np.linspace(1e-16,cutOffCalc,255))
        ax1.set_title(self.title)        
        ax2.contourf(np.linspace(0,max(self.datQ),np.shape(self.dataConvolve)[1]),np.linspace(0,int(max(self.datE)),np.shape(self.dataConvolve)[0]),self.dataConvolve,np.linspace(0,cutOffData,255))
        ax2.get_yticklabels()[-1].set_visible(False)
        ## Fine-tune figure; make subplots close to each other and hide x ticks for
        ## all but bottom plot.
        f.subplots_adjust(hspace=0)
        plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
        yticks = ax2.yaxis.get_major_ticks()
        yticks[-1].label1.set_visible(False)        
        yyl=plt.ylabel('')
		
        yyl.set_position((yyl.get_position()[0],1)) # This sets the top of the bottom axis as the reference point.
        yyl.set_verticalalignment('center')
        plt.xlabel('Q (r.l.u)')
        plt.ylabel('Energy (meV)')
        plt.xlim(xmax=findXLim(self.grid_z,self.QRangeActual) )
        plt.show()
        if save==True:
            plt.savefig('images/'+self.title+'.png')
            
    def gridGen(self,central=[4.11,4.11,1],width=1.5,res=3):
		""" Generates a grid of reciprocal space positions corresponding to the input file"""
        ax1=np.array([1,1,0])
        ax2=np.array([0,0,1])
        ax3=np.array([1,-1,0])        
        kpts=[]
        kpts=np.zeros((res**3,3))
        start=central-(width*(ax1+ax2+ax3))
        l=0
        for i in range(0,res):
            for j in range(0,res):
                for k in range(0,res):
                    kpts[l,:]=(start+width*((2*ax1/float(i+1))+(2*ax2/float(j+1))+(2*ax3/float(k+1))))
                    l+=1        
        rhomb=np.array(([4.980299, 4.980299, 0.0],[4.980299,0.0,4.980299],[0.0,4.980299,4.980299]))
        cubic=np.array(([9.960598, 0.0, 0.0],[0.0,9.960598,0.0],[0.0,0.0,9.960598]))
        rRhomb=reciplatt(rhomb)
        rCubic=reciplatt(cubic)        
        kptsRhomb=transformAll(kpts,rCubic,rRhomb)
        self.kpts=kpts
        self.kptsRhomb=kptsRhomb 

                                                                                                                            
        
    def initialise_from_file(self,fid):
        """ Runs all the above functions required to create a comparison between exp and calc data"""
        self.loadExpData(fid)
        self.whichFile()
        self.grab_data(self.calcFid)
        self.makeAllQ()
        self.fill_intensities()
        self.trim_intensities()
              

def blur_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction. 
        (from http://scipy.github.io/old-wiki/pages/Cookbook/SignalSmooth)
    """
    g = gauss_kern(n, sizey=ny)
    improc = convolve(im,g, mode='full')
    return(improc)    

def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions
    (from http://scipy.github.io/old-wiki/pages/Cookbook/SignalSmooth)
    """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()
    
def findXLim(arrIn, QRangeActual):
    """finds the first zero instance along x to find plot range"""
    for i in range(1,np.shape(arrIn)[1]):
        if arrIn[0,-i]!=0:
            return int(np.ceil(QRangeActual*(np.shape(arrIn)[1]-i)/float(np.shape(arrIn)[1])))
            
def transformAll(Q,rlatt1,rlatt2):
    """transform from latt1 to latt2 system assuming reciprocal lattices.
		Use this to transform a rhombohedral system to cubic or vice versa.
	"""
    q=np.zeros((len(Q),3))
    for i in range(0,len(Q)):
        q[i,:]=np.dot(np.linalg.inv(rlatt2),np.dot(Q[i,:],rlatt1))
    return q

def reciplatt(a):
	""" Convert from cartesian to reciprocal space"""
    a1=a[:,0]
    a2=a[:,1]
    a3=a[:,2]
    b1=2*np.pi*(np.cross(a2,a3))/(np.dot(a1,np.cross(a2,a3)))
    b2=2*np.pi*(np.cross(a3,a1))/(np.dot(a2,np.cross(a3,a1)))
    b3=2*np.pi*(np.cross(a1,a2))/(np.dot(a3,np.cross(a1,a2)))
    b=np.array([b1,b2,b3])
    return b;
            
			
def readHeader(file1):
	""" Takes header information from phonon file"""
	file1.readline()
	line = file1.readline()
	noOfIons= int(line.split()[3])
	line = file1.readline()
	noOfBranches = int(line.split()[3])
	line = file1.readline()
	noOfWavVec = int(line.split()[3])
	line = file1.readline()
	units = line.split()[2]
	file1.readline()#skip the IR units
	file1.readline()#skip ramen units
	file1.readline()#skip unit cell vectors title
	unitCell = zeros((3,3))
	line=file1.readline()#a
	unitCell[0,:] = array([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
	line=file1.readline()#b
	unitCell[1,:] = array([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
	line=file1.readline()#c
	unitCell[2,:] = array([float(line.split()[0]),float(line.split()[1]),float(line.split()[2])])
	unitCell = transpose(unitCell)
	file1.readline()#skip the title for ionic positions
	atPos = zeros((noOfIons,3))
	atType = []
	atMass = zeros(noOfIons)
	for i in range(noOfIons):
		line = file1.readline()
		lineCont = line.split()
		atPos[i,0] = float(lineCont[1]) 
		atPos[i,1] = float(lineCont[2]) 
		atPos[i,2] = float(lineCont[3]) 
		atPos[i,:] = dot(unitCell,atPos[i,:])#converts atpos into cartesian
		atType.append(lineCont[4])
		atMass[i] = float(lineCont[5])

return noOfIons, noOfBranches, noOfWavVec, units, unitCell, atPos, atType, atMass

def readEigens(noOfIons, noOfBranches, noOfWavVec, file1):
	"""Takes eigenvector data from .castep file"""
	hkl = zeros((noOfWavVec, 3))
	eigens = zeros((noOfWavVec,noOfBranches*noOfIons, 6))
	branchEnergy = zeros((noOfWavVec,noOfBranches))
	file1.readline()#skip the END Header comment
	for i in range(noOfWavVec):
		line = file1.readline()#read the Wave Vector
		lineCont = line.split()
		hkl[i,0] = float(lineCont[2])#h
		hkl[i,1] = float(lineCont[3])#k
		hkl[i,2] = float(lineCont[4])#l
		
		for j in range(noOfBranches):
			line = file1.readline()
			lineCont = line.split()
			branchEnergy[i,j] = float(lineCont[1]) # either in meV or cm^-1
			
		file1.readline() #skips the phonon eigenvectors title
		file1.readline() #skips the column titles
			
		for j in range(noOfBranches*noOfIons):
			line = file1.readline()
			lineCont = line.split()
			eigens[i,j,0] = float(lineCont[2])
			eigens[i,j,1] = float(lineCont[3])
			eigens[i,j,2] = float(lineCont[4])
			eigens[i,j,3] = float(lineCont[5])
			eigens[i,j,4] = float(lineCont[6])
			eigens[i,j,5] = float(lineCont[7])
return hkl, branchEnergy, eigens
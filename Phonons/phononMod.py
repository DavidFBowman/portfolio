from scipy.signal import convolve
import scipy.interpolate as inter
import numpy as np
import matplotlib.pyplot as plt
import VononCalcMods as VCM

class Phonon:
    def __init__(self):
        self.unitCell=[]
        self.atPos=[]
        self.atType=[]
        self.atMass=[]
        self.hkl=[]
        self.branchEnergy=[]
        self.eigens=[]
        self.formParams=[]
        self.noOfEigFiles=0
        
    def grab_data(self,fid):
        """takes data from .phonon CASTEP file"""
        self.calcFid=fid      
        unitCellTemp, atPosTemp, atTypeTemp,atMassTemp, hklTemp, branchEnergyTemp, eigensTemp, formParamsTemp =  VCM.extractPhonons(self.calcFid)
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
        """takes data from experiment output from HORACE of type x y E""" 
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
     
    def calcIntensity(self,Kpt,minE=-25.0,maxE=25.0):
        """returns calculated intensities within specified energy range at a given KPoint"""
        meashkl = Kpt
        temp = 30.0
        mix = 0.0
        width = 2.5
        energyRange = np.array([minE,maxE])
        elastic = 1.0
        idx = VCM.findWavVec(meashkl, self.hkl)
        xaxis, yaxis = VCM.calcDiffPat(self.unitCell,self.atPos, self.atMass, self.formParams, self.eigens[idx,:,:], self.branchEnergy[idx,:], energyRange, meashkl,temp, mix, width,elastic,1.)
        return yaxis
        
    def whichFile(self):
        """ reads title of input file and returns values to help construct allQ
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
        #for the hhh you'll need to change the length to QRange*2 in the line below
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
            #allQ=allQ             
        if self.kType==3:
            #QLen=30#len(self.hkl[0:30,:]
            #allQ=np.zeros(((30)*self.QRange,3))            
            for i in range(0,self.QRange):
                allQ[i*QLen:(i+1)*QLen,0]=self.hkl[0:QLen,0]+(2*i)
                allQ[i*QLen:(i+1)*QLen,1:3]=self.hkl[0:QLen,1:3]+i+(float(self.kOffset)/2.)
        if self.kType==4:
            for i in range(0,self.QRange):
                allQ[i*QLen:(i+1)*QLen,0]=self.hkl[0:QLen,0]+(2*i)
                allQ[i*QLen:(i+1)*QLen,1:3]=self.hkl[0:QLen,1:3]+(3*i)                    
        self.allQ=allQ
        
    def fill_intensities(self,energyBins=5000):
        """creates array of calculated intensities, removes negative energy"""
        intensities=np.zeros((len(self.allQ),energyBins))
        j=0
        for i in self.allQ:
            intensities[j,:]=self.calcIntensity(i)
            j+=1
        intensities=np.transpose(intensities)
        self.intensities=intensities[energyBins/2:energyBins,:]
        
    def trim_intensities(self):
        """reshapes intensities to same shape as data in calculated range"""
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
        #pts = np.array([[i,j] for i in np.linspace(0,1,np.size(np.unique(self.datE))) for j in np.linspace(0,1,np.size(np.unique(self.datQ)))] )
        pts = np.array([[i,j] for i in np.linspace(0,1,np.shape(datIntBlur)[0]) for j in np.linspace(0,1,np.shape(datIntBlur)[1])] )
        grid_x, grid_y = np.mgrid[0:1:colsj, 0:1:rowsj]
        grid_z = inter.griddata(pts, vals, (grid_x, grid_y), method='linear')        
        ##set areas outside exp range to 0s
        trimIntensity=np.multiply(trimIntensity,grid_z)
        self.trimIntensity=trimIntensity
        self.grid_z=grid_z
        self.dataConvolve=datIntBlur
        
    def plot_calculation(self,cutOff=10):
        """plots calculated data"""
        plt.contourf(np.linspace(0,self.QRangeActual,600),np.linspace(0,int(max(self.datE)),2500),self.intensities,np.linspace(1e-16,cutOff,255))
        plt.show()
        
    def plot_data(self,cutOff=200):
        """plots full calculation range"""
        plt.contourf(np.linspace(0,max(self.datQ),np.shape(self.dataConvolve)[1]),np.linspace(0,int(max(self.datE)),np.shape(self.dataConvolve)[0]),self.dataConvolve,np.linspace(0,cutOff,255))
        plt.show()
        
    def plot_all(self,cutOffCalc=10,cutOffData=200,save=False):
        """plots calculation with data"""
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
        """initialises based on experimental data input, calcuates and reshapes intensity"""
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
    #transform from latt1 to latt2 system assuming reciprocal lattices
    q=np.zeros((len(Q),3))
    for i in range(0,len(Q)):
        q[i,:]=np.dot(np.linalg.inv(rlatt2),np.dot(Q[i,:],rlatt1))
    return q

def reciplatt(a):
    a1=a[:,0]
    a2=a[:,1]
    a3=a[:,2]
    b1=2*np.pi*(np.cross(a2,a3))/(np.dot(a1,np.cross(a2,a3)))
    b2=2*np.pi*(np.cross(a3,a1))/(np.dot(a2,np.cross(a3,a1)))
    b3=2*np.pi*(np.cross(a1,a2))/(np.dot(a3,np.cross(a1,a2)))
    b=np.array([b1,b2,b3])
    return b;
            
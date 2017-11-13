from pipe import *
import astropy.io.fits as pf
import os.path
from shutil import copyfile, move


night_idx=0
snap_per_night=1
snap_idx=0

nssp=10
obs_idx=0
alphaX=0.0
alphaY=0.0

couplage=0.2

lambdaIR=1.65e-06
condi=30.

class Class_sys:
    def set(self, tomo,dm,fourier,its):
        self.tomo=tomo
        self.dm=dm
        self.fourier=fourier
        self.its=its

    def get_actuPos(self):
        x=self.dm.csX*self.tomo.diam()/2.
        y=self.dm.csY*self.tomo.diam()/2.
        return (x,y)


def init_sys(FILEPATH):
    tomo=Tomo_tiled(FILEPATH,"./")
    nMeas  =tomo.nMeas()
    nMeasTS=tomo.nMeasTS()
    print "nMeas:", nMeas,"\tnMeasTS:",nMeasTS

    dm=cmd.makeDM(tomo.get_Nssp()[-1],couplage,tomo.obs())
    fourier=FourierSpace(dm.Ndiam,tomo.diam(),lambdaIR)
    its=tomo.get_Nw()-1 #truth sensor index
    sys=Class_sys()
    sys.set(tomo,dm,fourier,its)
    return sys
    

def generate_otf(tomo,fourier):
    """
    tomo    : (Tomo_tiled)          : 
    fourier : (FourierSpace)        :
    """
    return OTF_telescope(fourier,tomo)
    


def generate_idx(dm,fourier,ncov):
    # indices de la carte 2D carree de taille nsspXnssp, ou il y a 'vraiment' des sous-pupilles
    #  creation du tableau des decalages
    Np = dm.Ndiam;
    print "NP:",Np
    # dm.csI: valid subapertures
    xx=np.tile(np.arange(Np),(Np,1)).flatten('C')[dm.csI]
    xx=-np.tile(xx,(xx.size,1))
    dx=xx-xx.T

    yy=np.tile(np.arange(Np),(Np,1)).flatten('F')[dm.csI]
    yy=-np.tile(yy,(yy.size,1))
    dy=yy-yy.T
    #  // transformation des decalages en indice de tableau
    dx += Np; 
    dy += Np; 

    ind = dx.flatten("F")+(Np*2-1)*(dy.flatten("F")-1)

    return ind.flatten()


def generate_abs2fi( N, ud, D, dactupix, coupling):

    y=np.tile((np.arange(N)+1-(N/2+1)) * ud /(D/2.), (N,1)).T
    x=np.copy(y.T)
    fi=cmd.funcInflu(x,y,coupling)
    abs2fi=np.abs(np.fft.fft2(fi.T))**2

    return abs2fi.T


def generate_Dx(tomo,dm,fourier,its):
    #creating Dx
    D = cmd.intermat(dm,tomo.obs(),its,tomo.get_Nssp(),tomo.get_Nsubap(),tomo.diam())
    Dx = cmd.computeCommandMat(D,condi=condi)


    return Dx


def write_fits(data,fileName):
    hdu = fits.PrimaryHDU(data)
    hdu.writeto(fileName)


#print "generate_files return: otf,idx,abs2fi,Dx, psf"
def generate_files(sys,cvvFile=None):
    otf    =generate_otf(sys.tomo,sys.fourier).T
    abs2fi =generate_abs2fi(sys.fourier.N,sys.fourier.ud,sys.tomo.diam(),sys.fourier.dactupix, sys.dm.x0 ).T
    Dx     =generate_Dx(sys.tomo,sys.dm,sys.fourier,sys.its).T
    print "shape Dx",Dx.shape
    idx    =generate_idx(sys.dm,sys.fourier,Dx.shape[1])
    if(cvvFile is not None):
        ff=pf.open(cvvFile);
        Cvv=ff[0].data.astype(np.float64)
        ff.close()
        dphi_tomo=intersample(getMap(Cvv,1,sys.dm),sys.fourier.N,sys.fourier.ud,sys.tomo.diam(),sys.fourier.dactupix, sys.dm.x0, sys.fourier.lambdaIR )
        print sys.fourier.N,sys.fourier.ud,sys.tomo.diam(),sys.fourier.dactupix, sys.dm.x0, sys.fourier.lambdaIR
        psf=addBwFitting2Psf(sys.tomo,sys.fourier,dphi_tomo)
    else:
         psf=None
         dphi_tomo=None

    return otf,idx,abs2fi,Dx,psf,dphi_tomo


def saveGeneratedFiles(otf,idx,abs2fi,Dx,psf=None,PATH="./"):
    pf.writeto(PATH+"/otftel.fits",otf.astype(dtype=np.float64))
    pf.writeto(PATH+"/idx.fits",idx.astype(dtype=np.int64))
    pf.writeto(PATH+"/abs2fi.fits",abs2fi.astype(dtype=np.float64))
    pf.writeto(PATH+"/Dx.fits",Dx.astype(dtype=np.float64))
    if(psf is not None):
        pf.writeto("/psf.fits",psf.astype(dtype=np.float64))



def printdif(a,b):
    tmp=np.copy(a)
    tmp[np.where(tmp==0)]=1
    print "min,max:",a.min(),a.max()
    ind=np.abs(a-b).argmax()
    print "max diff  ", a.item(ind)-b.item(ind),"at",ind, ": ",a.item(ind)," / ",b.item(ind)
    ind=np.abs((a-b)/tmp).argmax()
    print "normalized", a.item(ind)-b.item(ind),"at",ind, ": ",a.item(ind)," / ",b.item(ind)






def create_sysParam(diamtel,nNgs,nLgs,lgs_cst,ncpu=1,nPsf=1,seed=1234,lgsFlux=7.e6,ngsMag=None):
    """create a system parameter

    diamtel : double : telescope diameter
    nNgs    : int    : number of ngs
    nLgs    : int    : number of lgs
    lgs_cst : double : 
    ncpu    : int    : (default 1)number of threads (usefull only in lapack case)
    nPsf    : int    : (default 1)side of PSFs map
    seed    : int    : (default 1234)random seed 
    lgsFlux : float  : photon return at M1 for lgs (photon/m^2/s)
    ngsMag  : float or array of float : magnitude of the star (ngs)

    """

    np.random.seed(seed)

    Nw=nNgs+nLgs+1

    nssp=long(diamtel*2)
    obs=0.2
    itTime=0.004
    lgs_alt     =100000
    GsAlt       =np.zeros(Nw)   ; GsAlt[:nLgs]+=1./lgs_alt
    type_wfs    =np.ones(Nw)    ; type_wfs[:nLgs]+=1
    XPup        =np.zeros(Nw)
    YPup        =np.zeros(Nw)
    thetaML     =np.zeros(Nw)
    thetaCam    =np.zeros(Nw)
    sensibilite =np.ones(Nw)
    tracking    =np.ones(3)
    pasDPHI     =0.0001
    pixSizeNGS  =0.3
    pixSizeLGS  =1
    lambdaNGS   =650e-9
    lambdaLGS   =589e-9
    bdw         =3300e-10
    TNGS        =0.425
    TLGS        =0.382
    Tatmo       =0.84
    RON         =3
    spot_width  =1
    lgs_depth   =5000
    nTarget     =1
    targetX     =0
    targetY     =0

    # launch angle diameter for lgs (arcmin)
    radius=7.4 #for EELT
    #convert to arcsec
    radius=radius/2.*60
    #scale to telescope diameter
    radius*=(diamtel/38.542)

    # lgs positions (regular polygone)
    lgs_xpos=np.arange(nLgs)*2*np.pi/nLgs
    lgs_ypos=np.sin(lgs_xpos)*radius
    lgs_xpos=np.cos(lgs_xpos)*radius

    #ngs position (random)
    ngs_xpos=np.random.random(nNgs+1)*radius*2-radius
    ngs_ypos=np.random.random(nNgs+1)*radius*2-radius

    #concatenate lgs and ngs positions
    xpos=np.concatenate((lgs_xpos,ngs_xpos),axis=0)
    ypos=np.concatenate((lgs_ypos,ngs_ypos),axis=0)

    #reset truth sensor to axes
    xpos[-1]=0#xTarget
    ypos[-1]=0#yTarget


    nTarget=nPsf*nPsf
    targetX=np.arange(nPsf)-nPsf/2      #create list of index for x coordinates, center on 0
    targetX=np.tile(targetX,(nPsf,1))   #extend into array
    targetX=targetX*radius/2            #scale fov
    targetY=targetX.T                   #y coordinates
    
    targetX=targetX.flatten()           #flatten to 1d array
    targetY=targetY.flatten()           #flatten to 1d array


    if(ngsMag is None):
        mr=np.zeros(Nw-nLgs)+13
    elif(isinstance(ngsMag,float)):
	mr=np.zeros(Nw-nLgs)+ngsMag
    elif(isinstance(ngsMag,np.ndarray)):
        if(ngsMag.size != nNgs+1):
            raise ValueError("the ngsMag array must contain nNgs+1 elements")
        mr=np.copy(ngsMag)
    else:
        raise ValueError("ngsMag must be a float or a numpy array of float")


    #str(array)[1:-1].replace('\n','')): print array and remove new line and brakets

    f=open("sys-params.txt","w")

    f.write("tel diam in meters\n"              +str(diamtel))
    f.write("\ncent obs\n"                      +str(obs))
    f.write("\nTframe (s)\n"                    +str(itTime))
    f.write("\nN WFS\n"                         +str(Nw))
    f.write("\nnlgs\n"                          +str(nLgs))
    f.write("\nnTarget\n"                       +str(nPsf*nPsf))
    f.write("\nNssp n subaps per wfs\n"         +str(nssp))
    f.write("\nGsAlt\n"                         +str(GsAlt)[1:-1].replace('\n',''))
    f.write("\ntype\n"                          +str(type_wfs)[1:-1].replace('\n',''))
    f.write("\nalphaX in arcsec //the last value is for the truth sensor, overwrite by the function matcov_init_tomo_tiled\n"+str(xpos)[2:-1].replace('\n',''))
    f.write("\nalphaY in arcsec //the last value is for the truth sensor, overwrite by the function matcov_init_tomo_tiled\n"+str(ypos)[2:-1].replace('\n',''))
    f.write("\nXPup\n"                          +str(XPup)[1:-1].replace('\n',''))
    f.write("\nYPup\n"                          +str(YPup)[1:-1].replace('\n',''))
    f.write("\nthetaML\n"                       +str(thetaML)[1:-1].replace('\n',''))
    f.write("\nthetaCam\n"                      +str(thetaCam)[1:-1].replace('\n',''))
    f.write("\nsensibilite\n"                   +str(sensibilite)[1:-1].replace('\n',''))
    f.write("\ntracking\n"                      +str(tracking)[1:-1].replace('\n',''))
    f.write("\npasDphi\n"                       +str(pasDPHI))
    f.write("\nncpu\n"                          +str(ncpu))
    f.write("\nmagnitude of NGS\n"              +str(mr)[1:-1].replace('\n',''))
    f.write("\nFlux for LGS\n"                  +str(lgsFlux))
    f.write("\npixel size for NGS\n"            +str(pixSizeNGS))
    f.write("\npixel size for LGS\n"            +str(pixSizeLGS))
    f.write("\nwave length for NGS (meters)\n"  +str(lambdaNGS))
    f.write("\nwave length for LGS (meters)\n"  +str(lambdaLGS))
    f.write("\nbandwidth (meters)\n"            +str(bdw))
    f.write("\ntransmission for NGS\n"          +str(TNGS))
    f.write("\ntransmission for LGS\n"          +str(TLGS))
    f.write("\natmosphere transmission\n"       +str(Tatmo))
    f.write("\nRead Out Noise (nb of e-)\n"     +str(RON))
    f.write("\nlgs_cst\n"                       +str(lgs_cst))
    f.write("\nspot_width (arcsec)\n"           +str(spot_width))
    f.write("\nlgs_alt (meters)\n"              +str(lgs_alt))
    f.write("\nlgs_depth (meters)\n"            +str(lgs_depth))
    f.write("\ntargetX  (arcsec)\n"             +str(targetX)[1:-1].replace('\n',''))
    f.write("\ntargetY  (arcsec)\n"             +str(targetY)[1:-1].replace('\n',''))

    f.close()

def createData(FILEPATH):
    sys=init_sys(FILEPATH)
    otf,idx,abs2fi,Dx,psf,dphi_tomo=generate_files(sys,cvvFile)
    saveGeneratedFiles(otf,idx,abs2fi,Dx,psf,PATH=FILEPATH)



def createFiles(diamtel,nNgs,nLgs,lgs_cst,ncpu=1,nPsf=1,seed=1234,cvvFile=None,FILEPATH="./",lgsFlux=7.e6,ngsMag=None):
    """create and save MOAO pipeline inputs files

    createFiles(diamtel,nNgs,nLgs,lgs_cst,ncpu=1,nPsf=1,seed=1234,cvvFile=None,FILEPATH="./")
    diamtel : double : telescope diameter
    nNgs    : int    : number of ngs
    nLgs    : int    : number of lgs
    lgs_cst : double : constant added to lgs covariances (simulate the fact that lgs cannot measure low order modes)
    ncpu    : int    : (default 1)number of threads (usefull only in lapack case)
    nPsf    : int    : (default 1)side of PSFs map
    seed    : int    : (default 1234)random seed (for NGS position generation)
    cvvFile : str : (optional) fits file containig a cvv matrix that will be used to generate the psf
    lgsFlux : float  : (default 7.e6) photon return at M1 for lgs (photon/m^2/s)
    ngsMag  : float or array of float : magnitude of the star (ngs)


    write the files:
        Dx.fits
	otftel.fits
	abs2fi.fits
	idx.fits
	psf.fits (if a cvvFile is provided)

    the matrices are written as fits files (needed by the intersample of the MOAO pipeline)
    the current directory must not contain files with the same name as the written files
    """
#FILEPATH: str : path to the directory containing the system parameters file (sys-params.txt)

    create_sysParam(diamtel,nNgs,nLgs,lgs_cst,ncpu,nPsf,seed,lgsFlux,ngsMag)
    if(FILEPATH!="./"):
        move("./sys-params.txt",FILEPATH+"/sys-params.txt")

    err=0
    if os.path.isfile(FILEPATH+"Dx.fits"):
        print "Dx.fits     already exists"
	err=1
    if os.path.isfile(FILEPATH+"otftel.fits"):
        print "otftel.fits already exists"
	err=1
    if os.path.isfile(FILEPATH+"abs2fi.fits"):
        print "abs2fi.fits already exists"
	err=1
    if os.path.isfile(FILEPATH+"idx.fits"):
        print "idx.fits    already exists"
	err=1
    if (os.path.isfile(FILEPATH+"psf.fits") and cvvFile is not None):
        print "psf.fits    already exists"
	err=1

    if err>0 :
        raise ValueError("output files already exists")

    sys=init_sys(FILEPATH)
    otf,idx,abs2fi,Dx,psf,dphi_tomo=generate_files(sys,cvvFile)
    saveGeneratedFiles(otf,idx,abs2fi,Dx,psf,PATH=FILEPATH)

    return sys

print "\nTo create input files for MOAO pipeline, use the createFiles function\n"
print createFiles.__doc__

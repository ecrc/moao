from pipe import *
from generate_file import *
from astropy.io import fits

#generation of the parameter file and input data files
#createFiles(8,3,4,0.1,ncpu=1,nPsf=1,seed=1234,cvvFile=None,FILEPATH="./")

filePath="../../datafile/check/nx456_nLayers2_wfs6_Nssp10/"

#create and initialize tomo struct
tomo=Tomo_tiled(filePath)
tomo.update_atm_params(0,8,0,0)
#number of (actual) measurements
NM=tomo.nMeas()-tomo.nMeasTS()
#nuber of measurements of the truth sensor
NTS=tomo.nMeasTS()

#create and initialize structure for 
isample=Intersample(tomo.get_Nssp()[tomo.get_Nw()-1],tomo.diam(),filePath)

#create array to work with (covriance matrices, reconstructor...
#python call routines working in column major (fortran layout)
#TODO: hide layout management in wrapper layer
Cmm=np.asfortranarray(np.zeros((NM,NM),dtype=np.float64))
Ctm=np.asfortranarray(np.zeros((NTS,NM),dtype=np.float64))
R  =np.asfortranarray(np.zeros((NTS,NM),dtype=np.float64))
Ctt=np.asfortranarray(np.zeros((NTS,NTS),dtype=np.float64))
Dx =np.asfortranarray(fits.getdata(filePath+"Dx.fits").astype(np.float64).T)

print "TOR"
tomo.matcov(Cmm,0,0,NM,1)
tomo.matcov(R,0,0,NTS,3)
py_reconstructor(Cmm, R, NM, NTS)
print "CEE/CVV"
tomo.matcov(Cmm,0,0,NM,1)
tomo.matcov(Ctm,0,0,NTS,3)
tomo.matcov(Ctt,0,0,NTS,4)
Cee,Cvv=py_computeCee_Cvv(Cmm,Ctt,Ctm,R,Dx)
print "psf"
psf=isample.process(Cvv)


import numpy as np
from numpy.linalg import svd as svd

class deformableMirror:
    def __init__(self):

        self.Nactu=0
        self.NactuDm=0
        self.Ntilt=0
        self.Ndiam=0
        self.x0=0

        self.interDm=np.zeros(0)
        self.csX=np.zeros(0)
        self.csY=np.zeros(0)
        self.csI=np.zeros(0)
        self.pertuVoltage=np.zeros(0)
        self.gain=np.zeros(0)
        self.alpha=np.zeros(0)

    def set_Nactu(self,n):
        self.Nactu=n

    def set_NactuDm(self,n):
        self.NactuDm=n

    def set_Ntilt(self,n):
        self.Ntilt=n

    def set_Ndiam(self,n):
        self.Ndiam=n

    def set_x0(self,x):
        self.x0=x

    def set_interDm(self,i):
        self.interDm=i

    def set_csX(self,c):
        self.csX=c

    def set_csY(self,c):
        self.csY=c

    def set_csI(self,c):
        self.csI=c

    def set_pertuVlotage(self,p):
        self.pertuVoltage=p

    def set_gain(self,g):
        self.gain=g

    def set_alpha(self,a):
        self.alpha=a

def makeDM(nssp, couplage, obs):
    """
    TODO struct dm
    """
    dm=deformableMirror()
    x=np.tile(np.linspace(-1,1,nssp+1),(nssp+1,1))
    y=np.copy(x.T)
    r=np.sqrt(x*x+y*y).flatten()
    nn=np.where((r<(1.0+1.4/nssp))*1 * (r>(obs-1./nssp))*1)[0]
    dm.set_Nactu(nn.size)
    dm.set_NactuDm(nn.size)
    dm.set_Ntilt(0)
    dm.set_Ndiam(nssp+1)
    dm.set_csX(x.flatten()[nn])
    dm.set_csY(y.flatten()[nn])
    # Indices of actuators = [3,4,5,6, 10,11,12,13,14,15, 17 ... ] in a Ndiam X Ndiam image
    dm.set_csI(nn)
    #  // valeur de x0 a passer a la fonction funcInflu(x,y,x0) pour obtenir un
    #  // couplage defini pour l'actu voisin, avec x et y exprimes en "rayon pupille"
    dm.set_x0(np.sqrt(-2/np.log(couplage))/nssp)
    print "Nber of actuators through diam :",dm.Ndiam
    print "Total number of actuators      :",dm.Nactu
    return dm

def intermat(dm, obs,its, Nssp ,Nsubap,diam):
    """
    nssp=Nssp[its] // A tester avec une structure tomo
    """
    nssp = Nssp[its]
    xy=vectmesFab(nssp, 2.0, obs, 0, 0, nssp, 0)
    #  matrix of distances (in x and y) from subaperture centre to DM
    csSize=dm.csX.size
    dx=np.tile(xy[0],(csSize,1)) - np.tile(dm.csX,(xy[0].shape[0],1)).T
    dy=np.tile(xy[1],(csSize,1)) - np.tile(dm.csY,(xy[1].shape[0],1)).T
    ssize = 2./nssp
    #  creation of MI
    """
    nsubap = Nsubap[its] // A tester avec une structure tomo
    """
    nsubap = Nsubap[its]
    #mia=np.zeros((dm.NactuDm,nsubap * 2))
    mia=np.zeros((nsubap * 2,dm.NactuDm))
    #  // here <mia> has no unit, it is just a difference between the edges of a subap
    #  // The result is expressed in meters.
    mia[:nsubap,:]= d_funcInflu(dx, dy, ssize, dm.x0).T
    mia[nsubap:,:]= d_funcInflu(dy, dx, ssize, dm.x0).T
    #  here we convert mia into an angle (in rad) per volt
    #  We need to divide the OPD by the size of subap to get an angle (radians)
    mia/=diam/nssp
    # And for grata-compatibility, we switch from radians to arcseconds
    mia *= 206265.
    return mia

def computeCommandMat(mia, thresh=0.0, condi=0.0, nmf=-1,start=0):
    # SVD
    U, w, V = svd(mia)
    # Check eigen values
    nmf = check_eigen(w[start:],thresh,np.sqrt(condi))
    # Inverting eigen values
    eigen_inv = invert_eigen(w, nmf)
    # Computing the pseudo-inverse with output which has transposed shape wrt input
    if(U.shape[0]<V.shape[0]):
        diag=np.eye(V.shape[0],U.shape[1])
        diag=diag*eigen_inv
    else:
        diag=np.eye(U.shape[1],V.shape[0])
        diag=(diag*eigen_inv).T
#    diag=np.eye(V.shape[0],U.shape[1])
#    diag=diag*eigen_inv
    mca=np.dot(np.dot(V.T,diag),U.T)
    return mca

def vectmesFab(nssp, diamtel,obs,XPup,YPup,diamPup,thetaML):
    x=np.linspace(-diamtel/2., diamtel/2., nssp+1).T
    x=np.tile(zcen(x),(nssp,1))
    y=np.copy(x.T)
    nn=np.where(getValidSubapArray(nssp, 1.0, obs))
    # taking magnification factor into account
    G = diamPup / float(nssp)
    x *= G
    y *= G
    # taking rotation into account
    th = thetaML * np.pi/180.;
    x_ = x*np.cos(th)-y*np.sin(th);
    y_ = x*np.sin(th)+y*np.cos(th);
    # taking pupil offsets (misregistration deviations) into account
    x_ += XPup
    y_ += YPup
    return x_.flatten()[nn],y_.flatten()[nn]

def d_funcInflu(x,y,ssize,x0):
#/* DOCUMENT d_funcInflu(x,y,ssize,x0)
#
#   Returns the difference at the edge of the subaperture. The
#   variables x, y, ssize and x0 are all in the same units. The
#   returned value does not depend on the chosen unit, as it is just a
#   difference between the 2 edges.
#
#   This difference is an optical path difference, expressed in METERS.
#   
    ss = ssize / 2
    return funcInflu(x + ss,y,x0) - funcInflu(x - ss,y,x0);

def funcInflu(x,y,x0):
#/* DOCUMENT opd_metres = funcInflu(x,y,x0)
#
#   The arguments <x>, <y> and <x0> must all be in the same units.
#   
#   In the particular case where <x0> is set to x0=dm.x0, <x> and <y>
#   must be expressed with no unit, they are variables defined over the
#   unitary circle of R=1.
#
#   In any case, the function returns the influence function (=OPD,
#   optical path difference) of the actuators expressed in METERS.
#     
#
#  // allez on va dire que ce sont des metres !
    return 1.e-6*np.exp( -(x*x+y*y)/(2*x0*x0) )

def check_eigen(eigen_val, thresh, condi):
    if condi != 0.0:
        thresh = np.max(eigen_val)/condi
    if thresh == 0.0:
        thresh = 1e-4
    nn = eigen_val[np.where(eigen_val<thresh)]
    if nn.size == 0:
        nmf = 0
    else:
        nmf = nn.size
    print('Number of filtered modes = ' + str(nmf))
    return nmf

def invert_eigen(eigen_val,nmf):
    if nmf>0:
        W = np.copy(eigen_val)
        for i in range(nmf):
            W[np.size(W)-1-i] = 1.0
        W = 1.0/W
        for i in range(nmf):
            W[np.size(W)-1-i] = 0.0
    else:
        W = 1.0/eigen_val
    return W

#from HRAA/trunk/codes/lib/python/tools.py
def getValidSubapArray(nssp, rext, rint, return2d=False):
    # The Grata case, tip-tilt sensor only.
    if(nssp==1):
        return [1];
    # to avoid some bug that eliminates useful central subapertures when obs=0.286
    if((nssp==7) and (rint>0.285 and rint<0.29)):
        rint=0.285;
        print "cas particulier"
    x = zcen(np.linspace(-1, 1, num=nssp+1))
    xx = []
    for i in range(nssp):
        xx = np.hstack((xx,x))
    x = np.reshape(xx,(nssp,nssp))
    y = np.transpose(x)
    r = np.sqrt(x*x+y*y)
    valid2dext = ((r<rext))*1
    valid2dint = ((r>=rint))*1
    valid2d = valid2dint*valid2dext
    if(return2d):
        return valid2d
    else:
        valid = np.reshape(valid2d, [nssp*nssp])
    return valid.tolist()
    
#from HRAA/trunk/codes/lib/python/tools.py
def zcen(data):
    data = np.array(data)
    if(len(data.shape)>1):
        print "oups zcen with dims > 1 not coded yet..."
        return 0
    tmp = tmp2 = []
    for i in range(len(data)-1):
        tmp = (float(data[i]) + float(data[i+1]))/2.
        tmp2 = np.append(tmp2, tmp)
    return tmp2

#"""
#MAIN TEST WITH RUN.I CONDITIONS
#"""
#dm  = makeDM(10,0.2,0.2)
#mia = intermat(dm,0.2,7,10,76,4.2)
#mca = computeCommandMat(mia,condi=30.0)
#print mca

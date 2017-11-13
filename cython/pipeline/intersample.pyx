from astropy.io import fits

cdef class FourierSpace:
    def __cinit__(self, int Nactu, double Dtel,lambdaIR):
        """
        :parameters:
            Nactu       : int    :number of actuators across diameter
            Dtel        : double :tel dima in meters
            lambdaIR    :        :IR wavelength in meterd
        """
        self.dactu=Dtel/(Nactu-1)
        # one will choose the number of pixels so that the area covered by
        # the DM in the phase-spectrum space is contained in the image with
        # enough margin at the edges
        self.N=8*2**long(np.log2(Nactu)+0.5)
        #number of pixels in the images of phase spectrum (Wfit, Waniso, ...)
        self.N=max(self.N,512)
        # one choose the total field of the Dphi image so that its contains a
        # 2*D area with some margin of D/2 at the edges, i.e. a total of 3*D ...
        self.champ_Dphi = 3*Dtel #taille de l'image de Dphi
        #size of the pixel of Dphi image
        self.ud = self.champ_Dphi/float(self.N)

        # we now want dactu to be a multiple of ud
        km = long(self.dactu / self.ud + 1)
        print "Number of samples in dactu : %d\n", km
        self.ud = self.dactu / km
        self.champ_Dphi = self.N * self.ud
        self.dactupix = km

        print "ud                         : ", self.ud;
        print "champ_Dphi                 : ", self.champ_Dphi;

        self.uk = 1./self.champ_Dphi # taille pixel dans espace Wiener (metres^-1)

        self.lambdaIR = lambdaIR
        RASC = 206264.8062471
        self.uz = self.uk * lambdaIR * RASC
        self.uld = self.uk * Dtel



def getMap( np.ndarray[ndim=2,dtype=cython_real_t] Caa, double average, dm ):
    """  map = getMap( Caa, average,dm );

        Re-organize coefficients of matrix Caa and returns a 2D map of the
        spatial covariance
    """     
    # indices de la carte 2D carree de taille nsspXnssp, ou il y a 'vraiment' des sous-pupilles
    #  creation du tableau des decalages
    cdef int Np = dm.Ndiam;
    xx=np.tile(np.arange(Np),(Np,1)).flatten('C')[dm.csI]
    xx=-np.tile(xx,(xx.size,1))
    dx=xx-xx.T

    yy=np.tile(np.arange(Np),(Np,1)).flatten('F')[dm.csI]
    yy=-np.tile(yy,(yy.size,1))
    dy=yy-yy.T
    #  // transformation des decalages en indice de tableau
    dx += Np; 
    dy += Np; 

    # transformation d'un couple de decalages (dx,dy) en un indice du tableau 'Map'
    cdef np.ndarray[ndim=2,dtype=cython_real_t] Map =np.zeros((Np*2-1,Np*2-1),cython_real)
    cdef np.ndarray[ndim=2,dtype=cython_real_t] div =np.zeros((Np*2-1,Np*2-1),cython_real)
    ind = dx.flatten("F")+(Np*2-1)*(dy.flatten("F")-1)-1


    cdef int i
    cdef int x,y,s
    s=Map.shape[0]
    for i in range(ind.size):
        Map.itemset(ind[i],Map.item(ind[i])+Caa.item(i))
        div.itemset(ind[i],div.item(ind[i])+1)

    div[np.where(div==0)]=1
    if( average==1 ):
        Map /= div

    return Map;
    
    
def intersample(np.ndarray[ndim=2,dtype=cython_real_t] Cvvmap, int N, double ud, cython_real_t D, int dactupix,double coupling,lambdaIR):
    """res = intersample( Cvvmap, N, ud, D, dactupix, couplingfact, lambdaIR)

        Cvvmap is the 'map' of the Cvv matrix (cov matrix of tomo error
        expressed on volts). The "volts" unit must be used together with
        the influence function funcInflu(x,y,dm.x0) expressed in meters.

        Then, the result of intersample is in meter^2.
        :parameters:
            Cvvmap  : np.ndarray[ndim=2,dtype=np.float64_t] : 'map' of the Cvv matrix
            N       : int   :
            ud      : double: size of pixel of Dphi space (meters)
            D       : double: telescope diameter
            dactupix: int   : inter-actuator pitch expressed in pixels of size ud
            coupling: double: coupling factor of one influ func to the neighbour

    """

    print "Interpolating Dphi map"
    # computation of influence function
    y=np.tile((np.arange(N)+1-(N/2+1)) * ud /(D/2.), (N,1)).T
    x=np.copy(y.T)
    fi=cmd.funcInflu(x,y,coupling)

    Map=np.zeros((N,N))
    ncmap = N/2+1;
 
    #  size of the side of Cvvmap (always odd number)
    #  ncov = dimsof(Cvvmap)(0);
    #  nber of elements on each side of the center of Cvvmap
    nelem=(Cvvmap.shape[0]-1)/2
    start=ncmap-nelem*dactupix-1
    end=ncmap+nelem*dactupix
    Map[start:end:dactupix,start:end:dactupix]=Cvvmap.T

    #  Computing the phase correlation function
    #  One should have corr(0) = phase_variance.
    #  Computing corr(0) is done using the <v_i^2> (diagonal of Cvv).
    #  We decided that <v^2> is the average value for 1 single actuator (i.e. it's the average
    #  of the diagonal of Cvv).
    #  Then the average phase variance over the pupil equals to
    #  (1/S_pupil) * $_pupil(fi^2) * Nactu * <v^2>
    #  with S_pupil the surface, and $ is an integral. Nactu needs to be here
    #  because if it wasn't, we'd have computed the phase variance over the pupil
    #  with only 1 single actu moving.
    #  So, in our formula, we have replaced the value of (S_pupil/Nactu) by (dactu^2).
    #  The (dactu^2) needs to be expressed in pixels because our integral $(fi^2) is not
    #  a real integral : it's just summing pixels instead.
    #
    #  corr = fft(abs(fft(fi))^2 * fft(map), -1).re / (1.*numberof(fi) * dactupix^2) ;

    corr=np.fft.ifft2(np.abs(np.fft.fft2(fi))**2*np.fft.fft2(Map)).real /(1.*dactupix**2)
    #  From correlation to Dphi
    #  Dphi(r) = 2*C(0) - 2*C(r)
    #  We take advantage we need to do a multiplication to multiply by another factor
    #  in the same line. This is to translate dphi from m^2 into rd^2
    fact = 2 * (2*np.pi/lambdaIR)**2
    dphi = fact*corr[ncmap-1,ncmap-1] - fact*corr
    #dphi = np.roll(np.roll(dphi,ncmap-1,axis=0),ncmap-1,axis=1)
    dphi=np.fft.fftshift(dphi)
    #
    #  and ... we don't do that any more
    #   dphi1  = dphi / (abs(FTOtel)+1);
    #  
    return dphi;



def autocorrelation(a):
    """
   comutes the autocorrelation so that
   max(aa) == sum(a^2)
    """
    if(a.ndim==2):
        b = np.abs(np.fft.fft2(a))
        b = np.fft.ifft2(b*b).real*a.size
    elif(a.ndim==1):
        b = np.abs(np.fft.fft(a))
        b = np.fft.ifft(b*b).real
    else:
        print "error: autocorrelation: expect dim 1 or 2"
        return
    n2 = a.size # N*N
    b /= n2
    return b






def OTF_telescope(fourier,tomo):
    """otf = OTF_telescope(fourier)

   Computes the OTF of the telescope, so that
   > fft(OTF_telescope()).re
   produces a PSF normalized with max(psf)=SR=1.0

    """
    N = fourier.N
    # computation of pupil
    x = fourier.ud / (tomo.diam()/2.) * (np.arange(N)+1-(N/2+1))   # x exprime en rayon pupille
    x2=np.tile(x*x,(x.size,1))
    r=np.sqrt(x2+x2.T)

    pup=(r<=1.0)*1 * (r>tomo.obs())*1
    #  factor that will normalize the psf
    #  with PSF(0)=1.00 when diffraction-limited
    #  surface_pup_m2 = tomo.tel.diam^2*(1-tomo.tel.obs^2)*pi/4;
    surface_pup_m2=tomo.diam()**2*(1-tomo.obs()**2)*np.pi/4.
    surface_pup_pix=surface_pup_m2/fourier.ud**2
    factnorm=surface_pup_pix**2
    #  compute FTO of telescope. Computing the psf using
    #  just fft(FTO).re produces a psf with max(psf)=SR
    #  In fact, FTOtel is normalized so that sum(FTOtel)=1.
    #  FTOtel = autocorrelation(pup) / factnorm;
    FTOtel=autocorrelation(pup)/factnorm
    return FTOtel;


def addBwFitting2Psf(tomo,fourier,dphi_tomo):
    psf=np.fft.ifft2(OTF_telescope(fourier,tomo)*np.exp(-0.5*dphi_tomo)).real
    psf=np.fft.fftshift(psf*(psf.size))
    return psf



cdef class Intersample:
    def __cinit__(self,Nssp,diam,files_path):
        cdef isample_struct isample_str
        self.isample=isample_str
        Dx=fits.getdata(files_path+"/Dx.fits").T
        nact=Dx.shape[0]
        intersample_prepare(&self.isample,nact*nact,Nssp+1,diam,files_path)
        if(intersample_init(&self.isample)!=1):
            raise ValueError("intersample_init failed")

    def process(self,np.ndarray[ndim=2,dtype=np.float64_t]Cvv):
        intersample_process(&self.isample,<double*>Cvv.data)
        psf=np.asfortranarray(np.zeros((self.isample.N,self.isample.N)))
        for i in range(self.isample.N):
            for j in range(self.isample.N):
                psf[i,j]=self.isample.dphi[j*self.isample.N+i]
        return psf

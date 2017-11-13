cdef class Tomo_tiled:

    def __cinit__(self,bytes sys_path, bytes atm_path):
        cdef tomo_struct tomo_str
        self.tomo=tomo_str
        cdef char *sp=sys_path
        cdef char *ap=atm_path
        self.tomo.sys_path=sp
        self.tomo.atm_path=ap
        init_tomo_sys(&self.tomo)

    def __dealloc__(self):
        matcov_free_tomo_tiled(&self.tomo)


    def init_tomo_tiled(self,long nssp,bytes sys_path,bytes atm_path,int night_idx, int snapshots_per_night, 
        int snapshot_idx, int obs_idx, real_t alphaX, real_t alphaY):
            return matcov_init_tomo_tiled( &self.tomo, sys_path,atm_path,night_idx, snapshots_per_night,
                                    snapshot_idx, obs_idx,alphaX, alphaY)


    def update_tomo_tiled(self,t=-1):
        matcov_update_tomo_tiled(&self.tomo,t)

    def update_atm_params(self,int night_idx, int snapshots_per_night, int snapshot_idx, int obs_idx):
        return  matcov_update_atm_params(&self.tomo, night_idx,snapshots_per_night,snapshot_idx, obs_idx)

    def free_tiled(self):
        matcov_free_tomo_tiled(&self.tomo)

    def nMeas(self):
        return matcov_getNumMeasurements(&self.tomo)

    def nMeasTS(self):
        return matcov_getNumMeasurementsTS(&self.tomo)

    def cn2(self):
        c=np.zeros(self.tomo.Nlayer)
        for i in range(self.tomo.Nlayer):
            c[i]=self.tomo.cn2[i]
        return c

    def get_Nssp(self):
        c=np.zeros(self.tomo.Nw,dtype=np.int64)
        for i in range(self.tomo.Nw):
            c[i]=self.tomo.Nssp[i]
        return c

    def get_GsAlt(self):
        c=np.zeros(self.tomo.Nw,dtype=cython_real)
        for i in range(self.tomo.Nw):
            c[i]=self.tomo.GsAlt[i]
        return c

    def matcov(self,np.ndarray[ndim=2,dtype=cython_real_t] A,m0,n0,lda,part):
        matcov_comp_tile(<real_t*>A.data,A.shape[0],A.shape[1],m0,n0,lda,&self.tomo,part)



    def get_Nw(self):
        return self.tomo.Nw

    def get_Nsubap(self):
        cdef int i
        cdef np.ndarray[ndim=1,dtype=np.int64_t] Nsubap=np.zeros(self.tomo.Nw,dtype=np.int64)
        for i in range(self.tomo.Nw):
            Nsubap[i]=self.tomo.Nsubap[i]
        return Nsubap


    def get_u(self):
        U=np.zeros(self.tomo.Nx*self.tomo.Nlayer)
        for i in range(U.size):
            U[i]=self.tomo.u[i]
        return U

    def get_v(self):
        V=np.zeros(self.tomo.Nx*self.tomo.Nlayer)
        for i in range(V.size):
            V[i]=self.tomo.v[i]
        return V


    def get_X(self):
        U=np.zeros(self.tomo.Nx)
        for i in range(U.size):
            U[i]=self.tomo.X[i]
        return U

    def get_Y(self):
        V=np.zeros(self.tomo.Nx)
        for i in range(V.size):
            V[i]=self.tomo.Y[i]
        return V

    def obs(self, real_t o=-1.):
        if(o<0):
            return self.tomo.obs
        else:
            self.tomo.obs=o

    def diam(self, real_t d=-1.):
        if(d<0):
            return self.tomo.DiamTel
        else:
            self.tomo.DiamTel=d


    def get_alphaX(self):
        cdef np.ndarray[ndim=1,dtype=cython_real_t] a=np.zeros(self.tomo.Nw)
        cdef int i
        for i in range(self.tomo.Nw):
            a[i]=self.tomo.alphaX[i]
        return a

    def get_alphaY(self):
        cdef np.ndarray[ndim=1,dtype=cython_real_t] a=np.zeros(self.tomo.Nw)
        cdef int i
        for i in range(self.tomo.Nw):
            a[i]=self.tomo.alphaY[i]
        return a


    def get_XPup(self):
        cdef np.ndarray[ndim=1,dtype=cython_real_t] p=np.zeros(self.tomo.Nw)
        cdef int i
        for i in range(self.tomo.Nw):
            p[i]=self.tomo.XPup[i]
        return p

    def get_YPup(self):
        cdef np.ndarray[ndim=1,dtype=cython_real_t] p=np.zeros(self.tomo.Nw)
        cdef int i
        for i in range(self.tomo.Nw):
            p[i]=self.tomo.YPup[i]
        return p

    def get_thetaML(self):
        cdef np.ndarray[ndim=1,dtype=cython_real_t] t=np.zeros(self.tomo.Nw)
        cdef int i
        for i in range(self.tomo.Nw):
            t[i]=self.tomo.thetaML[i]
        return t

    def get_diamPup(self):
        cdef np.ndarray[ndim=1,dtype=cython_real_t] t=np.zeros(self.tomo.Nw)
        cdef int i
        for i in range(self.tomo.Nw):
            t[i]=self.tomo.diamPup[i]
        return t

    def get_sspSize(self):
        cdef np.ndarray[ndim=1,dtype=cython_real_t] s=np.zeros(self.tomo.Nw)
        cdef int i
        for i in range(self.tomo.Nw):
            s[i]=self.tomo.sspSize[i]
        return s

    def get_h(self):
        cdef np.ndarray[ndim=1,dtype=cython_real_t] h=np.zeros(self.tomo.Nlayer)
        cdef int i
        for i in range(self.tomo.Nlayer):
            h[i]=self.tomo.h[i]
        return h

    def get_L0(self):
        cdef np.ndarray[ndim=1,dtype=cython_real_t] l=np.zeros(self.tomo.Nlayer)
        cdef int i
        for i in range(self.tomo.Nlayer):
            l[i]=self.tomo.L0[i]
        return l


    def extract(self, i=-1,j=-1):
        """extract the data of a the Tomo as a list
        if i>=0 and j>=0 the list will contain the data of the wfs number i to j
        """
        l=[]

        if(i<0 or j<0):
            l.append(self.tomo.Nw)
            l.append(self.tomo.Nx)
            l+=self.get_X().tolist()
            l+=self.get_Y().tolist()
            l.append(self.tomo.DiamTel)
            l.append(self.tomo.obs)
            l+=self.get_Nsubap().tolist()
            l+=self.get_Nssp().tolist()
            l+=self.get_GsAlt().tolist()
            l+=self.get_alphaX().tolist()
            l+=self.get_alphaY().tolist()
            l+=self.get_XPup().tolist()
            l+=self.get_YPup().tolist()
            l+=self.get_thetaML().tolist()
            l+=self.get_diamPup().tolist()
            l+=self.get_sspSize().tolist()

        l.append(self.tomo.Nlayer)
        l+=self.cn2().tolist()
        l+=self.get_h().tolist()
        l+=self.get_L0().tolist()

        return l

def save_tomo(bytes filename, Tomo_tiled t):
    """save a tomo and the corresponding covariance matrix
    
    :parameters:
        filename: (string) : name of the file to write (replace existing file)
        t       : (Tomo)   : tomo to be saved
    """
    tomo_struct=t.extract()
    Cmm=t.get_data()

    np.savez(filename,tomo_struct=tomo_struct,Cmm=0)




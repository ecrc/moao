cimport numpy as np
cimport cython

include "par.pxi"

IF USE_SINGLE:
    ctypedef float real_t
    ctypedef np.float32_t cython_real_t

ELSE:
    ctypedef double real_t
    ctypedef np.float64_t cython_real_t




cdef extern from "moao_defs.h":
    pass


#cdef extern from "moao_lapack.h":
#    void reconstructor(long nmeas, long nmeasts, real_t *Cmm, long ldCmm, real_t *Ctm, long ldCtm, long ngal);
#    void compute_Cee_Cvv(long nmeas, long nmeasts, long nact, real_t *Cmm, long ldCmm, real_t *Cpp, long ldCpp, real_t *Cpm, long ldCpm, real_t *R, long ldR, real_t *Dx, long ldDx, real_t *Cee, long ldCee, real_t *Cvv, long ldCvv, real_t *Tmp, long ldTmp);
#    void Xgeadd(char trans, long M, long  N, real_t alpha, real_t *A, long ldA, real_t beta, real_t *C, long ldC);


#/usr/include/cblas.h
cdef extern from "cblas.h":
    cdef enum CBLAS_ORDER:
        CblasRowMajor=101
        CblasColMajor=102
    cdef enum CBLAS_TRANSPOSE:
        CblasNoTrans=111
        CblasTrans=112
        CblasConjTrans=113
    cdef enum CBLAS_UPLO:
        CblasUpper=121
        CblasLower=122
    cdef enum CBLAS_DIAG:
        CblasNonUnit=131
        CblasUnit=132
    cdef enum CBLAS_SIDE:
        CblasLeft=141
        CblasRight=142

#    void cblas_Xtrsm(CBLAS_ORDER Order, CBLAS_SIDE Side,
#                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
#                 CBLAS_DIAG Diag, int M, int N,
#                 real_t alpha, real_t *A, int lda,
#                 real_t *B, int ldb);
#
#    void cblas_Xsyr2k(CBLAS_ORDER Order, CBLAS_UPLO Uplo,
#                  CBLAS_TRANSPOSE Trans, int N, int K,
#                  real_t alpha, real_t *A, int lda,
#                  real_t *B, int ldb, real_t beta,
#                  real_t *C, int ldc);
#
#    void cblas_Xsymm(CBLAS_ORDER Order, CBLAS_SIDE Side,
#                 CBLAS_UPLO Uplo, int M, int N,
#                 real_t alpha, real_t *A, int lda,
#                 real_t *B, int ldb, real_t beta,
#                 real_t *C, int ldc);
#
#    void cblas_Xgemm(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
#                 CBLAS_TRANSPOSE TransB, int M, int N,
#                 int K, real_t alpha, real_t *A,
#                 int lda, real_t *B, int ldb,
#                 real_t beta, real_t *C, int ldc);


#/home/ndoucet/build/lapack-3.6.0/LAPACKE/include/
cdef extern from "lapacke.h":

    cdef int LAPACK_ROW_MAJOR
    cdef int LAPACK_COL_MAJOR

#    int LAPACKE_Xpotrf( int matrix_layout, char uplo, int n, real_t* a,int lda )
#    int LAPACKE_Xlacpy( int matrix_layout, char uplo, int m, int n, const real_t* a, int lda, real_t* b, int ldb )

    cdef void matcov_styc(tomo_struct tomo, real_t *data)
    cdef void matcov_cpp_styc(tomo_struct tomo, real_t *data)
    cdef void printptr(real_t *ptr)



cdef extern from "tomo_struct.h":
    cdef struct tomo_struct:
        long Nw;  # number of wavefront sensors
        long nTarget; #number of targets

        # pointers on arrays containing corrdinates of sub-apertures
        # X and Y are biiiiig arrays from yorick containing all the subap
        # coordinates of all the WFSs, one after the other.
        real_t *X;
        real_t *Y;

        real_t DiamTel; # telescope Diameter

        real_t obs;#telescope obscuration

        # array of the number of subap of each WFS, contains Nw elements
        long *Nsubap;

        # array of the number of subap of each WFS along the telescop diameter, contains Nw elements
        long *Nssp;

        # array of the inverse of the guide star altitude (in 1/meters), contains Nw elements
        real_t *GsAlt;

        # type of WFS, 0, 1, 2 or 3. 0 is unused, 1=NGS, 2=LGS, 3=TipTilt-guide star
        int *type;

        # Pointing directions of WFS
        real_t *alphaX;         # pointing direction in X, arcseconds
        real_t *alphaY;         # pointing direction in Y, arcseconds
        real_t *targetX;        # pointing direction in X, radian
        real_t *targetY;        # pointing direction in Y, radian

        # Deviations of WFSs
        real_t *XPup;           # pupil shift of the WFS, in meters
        real_t *YPup;           # pupil shift of the WFS, in meters
        real_t *thetaML;        # rotation of microlenses
        real_t *thetaCam;       # rotation of camera
        real_t *sensibilite;    # sensitivity coeff of this WFS
        real_t *diamPup;        # magnification factor of this WFS
        real_t *sspSize;        # subaperture size of this WFS

        # PROFILE
        long Nlayer;            # number of layers in the profile
        real_t *cn2;            # profile strengh, units TBD ..
        real_t *h;              # altitude of layers (meters)
        real_t *L0;             # outer scale (meters)

        real_t rmax;            # maximum distance between subapertures (computed with yorick)
        real_t *tracking;       # telescope tracking error parameters (x^2, y^2 and xy), units : arcsec^2

        real_t pasDPHI;         # Precision of DPHI precomputation.
        int ncpu;               #Number of CPU used (only with openMP)
        int part;               #Computed part of the cov. matrix. 0: complete 1: cmm 2: cpp 3: cpm
        int Nx;
        int Nslopes;
        int nlgs;

        real_t lgs_cst;
        real_t noise_var;
        real_t spot_width;
        real_t lgs_alt;
        real_t lgs_depth;

        real_t lambdaNGS;	    # lambda for the NGS in meter
        real_t lambdaLGS;       # lambda for the LGS in meter
        real_t sNGS2;           # square of the seeing at NGS lambda
        real_t sLGS2;           # square of the seeing at LGS lambda

        real_t *pixSize;

        real_t Tatmo;           # atmosphere transmission (in visible at z=30°
        real_t *throughput;     # throughput of NGS

        int RON;                # read out noise (nb of e-)

        real_t qr;              # photometric Flux offset
        real_t *mr;             # guide star magnitude
        real_t bdw;             # bandwidth (A)

        real_t Tframe;
        real_t *lgsFlux;
        real_t *lgsExt;         #extension of lgs
        real_t *lgsTheta;       # angle of lgs extension
        real_t *noiseNGS;
        real_t *noiseLGSxx;
        real_t *noiseLGSyy;
        real_t *noiseLGSxy;

        #Ali
        long   *indexL0;
        real_t *L0diff;
        real_t *tabDPHI;
        real_t *u;
        real_t *v;
        real_t *sspSizeL;
        long    nsubaps_offaxis;

        char    sys_path[512];
        char    atm_path[512];


    int init_tomo_sys(tomo_struct *tomo);
    int init_tomo_atm(tomo_struct *tomo, int night_idx, int snapshots_per_night, int snapshot_idx, int obs_idx);
    int matcov_init_tomo_tiled(tomo_struct *tomo, char* sys_path, char* atm_path, int night_idx, int snapshots_per_night, int snapshot_idx, int obs_idx, real_t alphaX, real_t alphaY);
    void matcov_update_tomo_tiled(tomo_struct *tomo,int t);
    void matcov_free_tomo_tiled(tomo_struct *tomo);

    int matcov_getNumMeasurements(tomo_struct *tomo);
    int matcov_getNumMeasurementsTS(tomo_struct *tomo);

    int matcov_update_atm_params(tomo_struct *tomo, int night_idx, 
        int snapshots_per_night, int snapshot_idx, int obs_idx);    

    void cblas_Xtrsm(CBLAS_ORDER Order, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, CBLAS_TRANSPOSE TransA,
                 CBLAS_DIAG Diag, int M, int N,
                 real_t alpha, real_t *A, int lda,
                 real_t *B, int ldb);

    void cblas_Xsyr2k(CBLAS_ORDER Order, CBLAS_UPLO Uplo,
                  CBLAS_TRANSPOSE Trans, int N, int K,
                  real_t alpha, real_t *A, int lda,
                  real_t *B, int ldb, real_t beta,
                  real_t *C, int ldc);

    void cblas_Xsymm(CBLAS_ORDER Order, CBLAS_SIDE Side,
                 CBLAS_UPLO Uplo, int M, int N,
                 real_t alpha, real_t *A, int lda,
                 real_t *B, int ldb, real_t beta,
                 real_t *C, int ldc);

    void cblas_Xgemm(CBLAS_ORDER Order, CBLAS_TRANSPOSE TransA,
                 CBLAS_TRANSPOSE TransB, int M, int N,
                 int K, real_t alpha, real_t *A,
                 int lda, real_t *B, int ldb,
                 real_t beta, real_t *C, int ldc);

    int LAPACKE_Xpotrf( int matrix_layout, char uplo, int n, real_t* a,int lda )
    int LAPACKE_Xlacpy( int matrix_layout, char uplo, int m, int n, const real_t* a, int lda, real_t* b, int ldb )

cdef extern from "matcov.h":

    void matcov_comp_tile(real_t* data, int nrows, int ncols, int xoffset, 
        int yoffset, int lda, tomo_struct *tomo, int part);

cdef extern from "intersample.h":
    cdef struct isample_struct:
        long N;             # size of fft support
        long np;            # number of subaps
        long nidx;          # size of index map
        #double lambda;      # waelength #unused here
        double dactupix;    # interactuator distance

        long *idx;
        double *abs2fi;
        double *otf_tel;
        double *dphi;
    void  intersample_prepare(isample_struct *isample, long nidx, long nsubap, float Dtel, char* path);
    int intersample_init( isample_struct *isample);
    int intersample_process( isample_struct *isample, real_t *data);
    void intersample_free( isample_struct *isample);


cdef class Tomo_tiled:
    cdef tomo_struct tomo

cdef class Intersample:
    cdef isample_struct isample

cdef class FourierSpace:
  cdef readonly long N               # size of the support
  cdef readonly double uk            # size of pixels of the phase spectrum space
  cdef readonly double ud            # size of pixels of the Dphi space
  cdef readonly double champ_Dphi    # 
  cdef readonly double dactu         # 
  cdef readonly long dactupix        # 
  cdef readonly double lambdaIR      # 
  cdef readonly double uz            # 
  cdef readonly double uld           # 

  cdef readonly double e_fit         # 
  cdef readonly double e_bp          # 
  cdef readonly double e_others      # 
  
  cdef readonly np.ndarray k         # 
  cdef readonly np.ndarray kx        # 
  cdef readonly np.ndarray ky        # 

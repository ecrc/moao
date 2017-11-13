/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 *
 *MOAO general defines, typedefs and declarations
 **/


#ifndef __MOAO_DEFINES__
#define __MOAO_DEFINES__

#define RASC 206265.0 //!< convert radian to arcsec

#ifdef USE_SINGLE

typedef float    real_t;
#define REAL_IMG FLOAT_IMG
#define TREAL    TFLOAT
#define LAPACKE_Xpotrf LAPACKE_spotrf
#define LAPACKE_Xlacpy LAPACKE_slacpy
#define LAPACKE_Xlaset LAPACKE_slaset
#define cblas_Xtrsm    cblas_strsm
#define cblas_Xsyr2k   cblas_ssyr2k
#define cblas_Xsymm    cblas_ssymm
#define cblas_Xgemm    cblas_sgemm
#define cblas_Xscal    cblas_sscal

#define pXpotrf_    pspotrf_
#define pXtrsm_     pstrsm_
#define pXsyr2k_    pssyr2k_
#define pXsymm_     pssymm_
#define pXlaset_    pslaset_
#define pXlacpy_    pslacpy_
#define pXlascl_    pslascl_
#define pXgeadd_    psgeadd_
#define pXgemm_     psgemm_
#define pXgemr2d_   psgemr2d_
#define MPI_PREC    MPI_FLOAT


#else


typedef double   real_t;
#define REAL_IMG DOUBLE_IMG
#define TREAL    TDOUBLE
#define LAPACKE_Xpotrf LAPACKE_dpotrf
#define LAPACKE_Xlacpy LAPACKE_dlacpy
#define LAPACKE_Xlaset LAPACKE_dlaset
#define cblas_Xtrsm    cblas_dtrsm
#define cblas_Xsyr2k   cblas_dsyr2k
#define cblas_Xsymm    cblas_dsymm
#define cblas_Xgemm    cblas_dgemm
#define cblas_Xscal    cblas_dscal

#define pXpotrf_    pdpotrf_
#define pXtrsm_     pdtrsm_
#define pXsyr2k_    pdsyr2k_
#define pXsymm_     pdsymm_
#define pXlaset_    pdlaset_
#define pXlacpy_    pdlacpy_
#define pXlascl_    pdlascl_
#define pXgeadd_    pdgeadd_
#define pXgemm_     pdgemm_
#define pXgemr2d_   pdgemr2d_
#define MPI_PREC    MPI_DOUBLE

#endif

#define ROWMAJOR 1
#define COLMAJOR 0

#endif

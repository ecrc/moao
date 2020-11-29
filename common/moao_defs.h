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

namespace MOAO_CST{
const double RASC=206265.0; //!< convert radian to arcsec

const double Ga[11] = { 0, 12.067619015983075, 5.17183672113560444,
      0.795667187867016068, 0.0628158306210802181, 0.00301515986981185091,
      9.72632216068338833e-05, 2.25320204494595251e-06, 3.93000356676612095e-08,
      5.34694362825451923e-10, 5.83302941264329804e-12 };

const double Gma[11] = { -3.74878707653729304, -2.04479295083852408,
      -0.360845814853857083, -0.0313778969438136685, -0.001622994669507603,
      -5.56455315259749673e-05, -1.35720808599938951e-06,
      -2.47515152461894642e-08, -3.50257291219662472e-10,
      -3.95770950530691961e-12, -3.65327031259100284e-14 };

const double pi =  3.1415926535897932384626433; //!< PI

const double k1 = 0.1716613621245709486; //!<   2*gamma_R(11./6)*2^(-5./6)*pi^(-8./3)*(24*gamma_R(6./5)/5.)^(5./6);
const double k2 = 1.00563491799858928388289314170833; //!<  gamma_R(5./6)*2^(-1./6) (unsused)
const double k3 = 1.25331413731550012081;   //!< sqrt(pi/2)
const double a1 = 0.22222222222222222222;   //!< 2/9
const double a2 =-0.08641975308641974829;  //!< -7/89
const double a3 = 0.08001828989483310284;   //!< 175/2187

const double Xlim = 0.75*2*pi;
}





#ifdef USE_SINGLE

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

#endif

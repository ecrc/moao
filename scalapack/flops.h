/**
 *
 * @file flops.h
 *
 *  File provided by Univ. of Tennessee,
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @date 2010-12-20
 *
 **/
/*
 * This file provide the flops formula for all Level 3 BLAS and some
 * Lapack routines.  Each macro uses the same size parameters as the
 * function associated and provide one formula for additions and one
 * for multiplications. Example to use these macros:
 *
 *    FLOPS_ZGEMM( m, n, k )
 *
 * All the formula are reported in the LAPACK Lawn 41:
 *     http://www.netlib.org/lapack/lawns/lawn41.ps
 */
#ifndef _FLOPS_H_
#define _FLOPS_H_

#define FMULS_GEMM(__m, __n, __k) ((double)(__m) * (double)(__n) * (double)(__k))
#define FADDS_GEMM(__m, __n, __k) ((double)(__m) * (double)(__n) * (double)(__k))

#define FMULS_SYMM(__side, __m, __n) ( ( (__side) == 'L' ) ? FMULS_GEMM((__m), (__m), (__n)) : FMULS_GEMM((__m), (__n), (__n)) )
#define FADDS_SYMM(__side, __m, __n) ( ( (__side) == 'L' ) ? FADDS_GEMM((__m), (__m), (__n)) : FADDS_GEMM((__m), (__n), (__n)) )

#define FMULS_SYR2K(__n, __k) ((double)(__k) * (double)(__n) * (double)(__n)                )
#define FADDS_SYR2K(__n, __k) ((double)(__k) * (double)(__n) * (double)(__n) + (double)(__n))

#define FMULS_TRMM_2(__m, __n) (0.5 * (double)(__n) * (double)(__m) * ((double)(__m)+1.))
#define FADDS_TRMM_2(__m, __n) (0.5 * (double)(__n) * (double)(__m) * ((double)(__m)-1.))


#define FMULS_TRMM(__side, __m, __n) ( ( (__side) == 'L' ) ? FMULS_TRMM_2((__m), (__n)) : FMULS_TRMM_2((__n), (__m)) )
#define FADDS_TRMM(__side, __m, __n) ( ( (__side) == 'L' ) ? FADDS_TRMM_2((__m), (__n)) : FADDS_TRMM_2((__n), (__m)) )

#define FMULS_TRSM FMULS_TRMM
#define FADDS_TRSM FMULS_TRMM

#define FMULS_POTRF(__n) ((double)(__n) * (((1. / 6.) * (double)(__n) + 0.5) * (double)(__n) + (1. / 3.)))
#define FADDS_POTRF(__n) ((double)(__n) * (((1. / 6.) * (double)(__n)      ) * (double)(__n) - (1. / 6.)))

//Level 3 BLAS

//#define FLOPS_ZGEMM(__m, __n, __k) (6. * FMULS_GEMM((__m), (__n), (__k)) + 2.0 * FADDS_GEMM((__m), (__n), (__k)) )
//#define FLOPS_CGEMM(__m, __n, __k) (6. * FMULS_GEMM((__m), (__n), (__k)) + 2.0 * FADDS_GEMM((__m), (__n), (__k)) )
//#define FLOPS_DGEMM(__m, __n, __k) (     FMULS_GEMM((__m), (__n), (__k)) +       FADDS_GEMM((__m), (__n), (__k)) )
#define FLOPS_XGEMM(__m, __n, __k) (     FMULS_GEMM((__m), (__n), (__k)) +       FADDS_GEMM((__m), (__n), (__k)) )

//#define FLOPS_ZSYMM(__side, __m, __n) (6. * FMULS_SYMM(__side, (__m), (__n)) + 2.0 * FADDS_SYMM(__side, (__m), (__n)) )
//#define FLOPS_CSYMM(__side, __m, __n) (6. * FMULS_SYMM(__side, (__m), (__n)) + 2.0 * FADDS_SYMM(__side, (__m), (__n)) )
//#define FLOPS_DSYMM(__side, __m, __n) (     FMULS_SYMM(__side, (__m), (__n)) +       FADDS_SYMM(__side, (__m), (__n)) )
#define FLOPS_XSYMM(__side, __m, __n) (     FMULS_SYMM(__side, (__m), (__n)) +       FADDS_SYMM(__side, (__m), (__n)) )

//#define FLOPS_ZSYR2K(__n, __k) (6. * FMULS_SYR2K((__n), (__k)) + 2.0 * FADDS_SYR2K((__n), (__k)) )
//#define FLOPS_CSYR2K(__n, __k) (6. * FMULS_SYR2K((__n), (__k)) + 2.0 * FADDS_SYR2K((__n), (__k)) )
//#define FLOPS_DSYR2K(__n, __k) (     FMULS_SYR2K((__n), (__k)) +       FADDS_SYR2K((__n), (__k)) )
#define FLOPS_XSYR2K(__n, __k) (     FMULS_SYR2K((__n), (__k)) +       FADDS_SYR2K((__n), (__k)) )


//#define FLOPS_ZTRSM(__side, __m, __n) (6. * FMULS_TRSM(__side, (__m), (__n)) + 2.0 * FADDS_TRSM(__side, (__m), (__n)) )
//#define FLOPS_CTRSM(__side, __m, __n) (6. * FMULS_TRSM(__side, (__m), (__n)) + 2.0 * FADDS_TRSM(__side, (__m), (__n)) )
#define FLOPS_XTRSM(__side, __m, __n) (     FMULS_TRSM(__side, (__m), (__n)) +       FADDS_TRSM(__side, (__m), (__n)) )
//#define FLOPS_STRSM(__side, __m, __n) (     FMULS_TRSM(__side, (__m), (__n)) +       FADDS_TRSM(__side, (__m), (__n)) )

//#define FLOPS_ZPOTRF(__n) (6. * FMULS_POTRF((__n)) + 2.0 * FADDS_POTRF((__n)) )
//#define FLOPS_CPOTRF(__n) (6. * FMULS_POTRF((__n)) + 2.0 * FADDS_POTRF((__n)) )
#define FLOPS_XPOTRF(__n) (     FMULS_POTRF((__n)) +       FADDS_POTRF((__n)) )
//#define FLOPS_SPOTRF(__n) (     FMULS_POTRF((__n)) +       FADDS_POTRF((__n)) )

#define FLOPS_XCovMat(__m, __n, _l, _type_mat) ( (_type_mat) == 1 ? ((30. + (_l) * 115) * ( (double)(__m) * ((double)(__m) + 1.) / 2.0)) : ( ( _type_mat) == 3 ? ((30. + (_l) * 103) * (double)(__m) * (double)(__n) ) : ((30. + (_l) * 90) * (double)(__m) * (double)(__n) ) ) )

#define FMULS_LASCL(__m, __n) ((double)(__m) * (double)(__n))
#define FADDS_LASCL(__m, __n) (0)
#define FLOPS_XLASCL(__m, __n) (FMULS_LASCL(__m,__n) + FADDS_LASCL(__m,__n))

#define FMULS_GEADD(__m, __n) (2.0*(double)(__m) * (double)(__n))
#define FADDS_GEADD(__m, __n) ((double)(__m) * (double)(__n))
#define FLOPS_XGEADD(__m, __n) (FMULS_LASCL(__m,__n) + FADDS_LASCL(__m,__n))

#ifdef USE_INTERSAMPLE
  #define FLOPS_Intersample(__m, __n) ( (double)(__m) * (double)(__n) * 20.0)
  //TODO generate a correct formula
#endif

#endif /* _FLOPS_H_ */

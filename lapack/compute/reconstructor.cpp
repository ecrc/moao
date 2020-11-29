/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include <stdio.h>

#include "reconstructor.hpp"
#include "lapack_templates.hpp"

/**
* Compute the tomographic reconstructor
*
* input:
*   long nmeas   : number of measurements
*   long nmeasts : number of measurements of the truth sensor
*   real_t *Cmm  : covariance matrix between the sensors
*                  Cmm(nmeas, nmeas)
*   long ldCmm   : leading dimension of Cmm
*   real_t *Ctm  : covariance matrix between the truth sensors and the other sensors
*                  Ctm(nmeasts,nmeas)
*   long ldCtm   : leading dimension of Ctm
*
* output:
*   real_t *Ctm  : tomographic reconstructor
*
*/
template<typename T>
void MOAO_LAPACK::reconstructor(long nmeas, long nmeasts, T *Cmm, long ldCmm, T *Ctm, long ldCtm, long ngal){
    //int info = LAPACKE_Xpotrf(LAPACK_COL_MAJOR ,'L',nmeas,Cmm,ldCmm);
    int info = LAPACK<T>::potrf(LAPACK_COL_MAJOR ,'L',nmeas,Cmm,ldCmm);
    if (info != 0){
        printf("An error occured in dpotrf: %d\n", info);
        return;
    }

    long g;
    T *Ctm_g;
    for(g=0;g<ngal;g++){
        Ctm_g=&Ctm[nmeas*nmeasts*g];
        LAPACK<T>::trsm(CblasColMajor, CblasRight, CblasLower,
                    CblasTrans, CblasNonUnit,
                    nmeasts, nmeas,
                    1.,
                    Cmm, ldCmm,
                    Ctm_g, ldCtm);

        LAPACK<T>::trsm(CblasColMajor, CblasRight, CblasLower,
                    CblasNoTrans, CblasNonUnit,
                    nmeasts, nmeas,
                    1.,
                    Cmm, ldCmm,
                    Ctm_g, ldCtm);
    }
}

template void MOAO_LAPACK::reconstructor(long nmeas, long nmeasts, float *Cmm, long ldCmm, float *Ctm, long ldCtm, long ngal);
template void MOAO_LAPACK::reconstructor(long nmeas, long nmeasts, double *Cmm, long ldCmm, double *Ctm, long ldCtm, long ngal);

//for cython use
//void _dtrsm(int M, int N, double *A, int lda, double *B, int ldb){
//    //op( A )*X = alpha*B
//    cblas_dtrsm(CblasRowMajor,CblasLeft,CblasLower,CblasNoTrans,CblasNonUnit,M,N,1.,A,lda,B,ldb);
//}
//
//void _dsyr2k(int N, int M, double alpha, double *A, int lda, double *B, int ldb, double beta,double *C, int ldc){
//    //alpha*A*B' + alpha*B*A' + beta*C
//    cblas_dsyr2k(CblasRowMajor,CblasLower,CblasNoTrans,M,N,alpha,A,lda,B,ldb,beta,C,ldc);
//}
//
//void _dsymm(int M, int N, double alpha, double *A, int lda, double *B, int ldb, double beta, double *C, int ldc){
//    //C := alpha*B*A + beta*C
//    cblas_dsymm(CblasRowMajor,CblasRight,CblasLower,M,N,alpha,A,lda,B,ldb,beta,C,ldc);
//}
//
//void _dlacpy(int M,int N, double *A, int lda, double *B, int ldb){
//    //copy A in B
//    //cblas_dlacpy(CblasRowMajor,CblasLower,M,N,A,lda,B,ldb);
//    LAPACKE_dlacpy(CblasRowMajor,CblasLower,M,N,A,lda,B,ldb);
//}
//
//void _dgemm(int M,int N,int K,double alpha, double *A, int lda, double *B,int ldb,double beta, double *C,int ldc){
//
//    //C := alpha* A * B.T  + beta*C
//    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasTrans,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);
//
//}


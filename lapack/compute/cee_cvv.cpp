/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "cee_cvv.hpp"
#include "lapack_templates.hpp"

/**
* Compute the matrix Cee Cvv
*
* input:
*   long nmeas   : number of measurements
*   long nmeasts : number of measurements of the truth sensor
*   long nact    : number of actuators
*   real_t *Cmm  : covariance matrix between the sensors
*                  Cmm(nmeas, nmeas)
*   long ldCmm   : leading dimension of Cmm
*   real_t *Cpp  : covariance matrix of the truth sensor
*                  Cpp(nmeasts,nmeasts)
*   long ldCpp   : leading dimension of Cpp
*   real_t *Cpm  : covariance matrix between the truth sensors and the other sensors
*                  Cpm(nmeasts,nmeas)
*   long ldCpm   : leading dimension of Cpm
*   real_t *R    : tomographic reconstructor
*                  R(nmeasts,nmeas))
*   long ldR     : leading dimension of R
*   real_t *Dx   : special matrix
*                  Dx(nact,nmeasts)
*   long ldDx    : leading dimension of Dx
*   real_t *Cee  : covariance matrix
*                  Cee(nmeasts,nmeasts)
*   long ldCee   : leading dimension of Cee
*   real_t *Cvv  : covariance matrix
*                  Cvv(nact,nact)
*   long ldCvv   : leading dimension of Cvv
*   real_t *Tmp  : temporary buffer
*                  Tmp(nact,nmeasts)
*   long ldTmp   : leading dimension of Tmp
*
* output:
*   real_t *Cee : tomographic error
*   real_t *Cvv : covariance matrix
*
*/
template<typename T>
void MOAO_LAPACK::Cee_Cvv(long nmeas, long nmeasts, long nact, T *Cmm, long ldCmm, T *Cpp, long ldCpp,
                     T *Cpm, long ldCpm, T *R, long ldR, T *Dx, long ldDx, T *Cee, long ldCee, 
                     T *Cvv, long ldCvv, T *Tmp, long ldTmp){
    T alpha =-1.;
    T beta =1.;
    LAPACK<T>::syr2k(CblasColMajor, CblasLower, CblasNoTrans,
                 nmeasts, nmeas,
                 alpha,
                 Cpm, nmeasts,
                 R, nmeasts,
                 beta,
                 Cpp, nmeasts);

    alpha =1.;
    beta =0.;
    LAPACK<T>::symm(CblasColMajor, CblasRight, CblasLower,
                nmeasts, nmeas,
                alpha,
                Cmm, nmeas,
                R, nmeasts,
                beta,
                Cpm, nmeasts);

    int info = LAPACK<T>::laset(LAPACK_COL_MAJOR,'U',
                              nmeasts, nmeasts,
                              0., 0.,
                              Cee, nmeasts);

    LAPACK<T>::lacpy(LAPACK_COL_MAJOR,'L',
                   nmeasts, nmeasts,
                   Cpp, nmeasts,
                   Cee, nmeasts);

    alpha = 1.0;
    beta = 1.0;
    MOAO_LAPACK::geadd('T', nmeasts, nmeasts,
           alpha,
           Cee, ldCee,
           beta,
           Cee, ldCee);

    alpha = 0.5;
    LAPACK<T>::scal(nmeasts, alpha, Cee, nmeasts+1);


    alpha=1.;
    beta =1.;
    LAPACK<T>::gemm(CblasColMajor, CblasNoTrans, CblasTrans,
                nmeasts, nmeasts, nmeas,
                alpha,
                Cpm, nmeasts,
                R, nmeasts,
                beta,
                Cee, nmeasts);

    alpha=1.;
    beta =0.;
    LAPACK<T>::symm(CblasColMajor, CblasRight, CblasLower,
                nact, nmeasts,
                alpha,
                Cee, nmeasts,
                Dx, nact,
                beta,
                Tmp, nact);

    alpha=1.;
    beta =0.;
    LAPACK<T>::gemm(CblasColMajor, CblasNoTrans, CblasTrans,
                nact, nact, nmeasts,
                alpha,
                Tmp, nact,
                Dx, nact,
                beta,
                Cvv, nact);
}

template void MOAO_LAPACK::Cee_Cvv(long nmeas, long nmeasts, long nact, float *Cmm, long ldCmm, float *Cpp, long ldCpp,
                     float *Cpm, long ldCpm, float *R, long ldR, float *Dx, long ldDx, float *Cee, long ldCee, 
                     float *Cvv, long ldCvv, float *Tmp, long ldTmp);
template void MOAO_LAPACK::Cee_Cvv(long nmeas, long nmeasts, long nact, double *Cmm, long ldCmm, double *Cpp, long ldCpp,
                     double *Cpm, long ldCpm, double *R, long ldR, double *Dx, long ldDx, double *Cee, long ldCee, 
                     double *Cvv, long ldCvv, double *Tmp, long ldTmp);
/**
* Addition of 2 matrices
* Return C = beta*C + alpha*op(A)
*
* A    :(M,N)
* op(C):(M,N)
*
* op(C)=C            if trans='N'
*       transpose(C) if Trans='T'
*
*/
template<typename T>
void MOAO_LAPACK::geadd(char trans, long M, long N, T alpha, T *A, long ldA, T beta, T *B, long ldB){

    int i,j;

    if(trans=='N'){
        for(i=0; i<M;i++){
            for(j=0; j<N;j++){
                B[i+j*ldB] = beta*B[i+j*ldB] + alpha*A[i+j*ldA];
            }
        }
    }
    else if(trans=='T'){
        for(i=0; i<M;i++){
            for(j=0; j<N;j++){
                B[i*ldB+j] = beta*B[i+j*ldB] + alpha*A[i*ldA+j];
            }
        }
    }
}
template void MOAO_LAPACK::geadd(char trans, long M, long N, float alpha, float *A, long ldA, float beta, float *B, long ldB);
template void MOAO_LAPACK::geadd(char trans, long M, long N, double alpha, double *A, long ldA, double beta, double *B, long ldB);

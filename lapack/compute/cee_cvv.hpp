/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/
#ifndef CEE_CVV_H
#define CEE_CVV_H

#ifdef USE_MKL
#include <mkl_lapacke.h>
#include <mkl_cblas.h>
#else
#include <lapacke.h>
#include <cblas.h>
#endif

namespace MOAO_LAPACK{
template<typename T>
void Cee_Cvv(long nmeas, long nmeasts, long nact, T *Cmm, long ldCmm, T *Cpp, long ldCpp, T *Cpm, long ldCpm, T *R, long ldR, T *Dx, long ldDx, T *Cee, long ldCee, T *Cvv, long ldCvv, T *Tmp, long ldTmp);
template<typename T>
void geadd(char trans, long M, long  N, T alpha, T *A, long ldA, T beta, T *C, long ldC);
}
#endif // CEE_CVV_H

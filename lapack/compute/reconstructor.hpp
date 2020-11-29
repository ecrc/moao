/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/
#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#ifdef USE_MKL
#include <mkl_lapacke.h>
#include <mkl_cblas.h>
#else
#include <lapacke.h>
#include <cblas.h>
#endif

namespace MOAO_LAPACK{
template<typename T>
void reconstructor(long nmeas, long nmeasts, T *Cmm, long ldCmm, T *Ctm, long ldCtm, long ngal);
}

#endif // RECONSTRUCTOR_H

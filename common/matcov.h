/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef MATCOV_H
#define MATCOV_H

#include "tomo_struct.h"
#include "matcov_kernels.h"
#ifdef USE_CHAMELEON
#ifdef USE_GPU
#include "matcov_kernels_gpu.h"
#endif
#include <starpu.h>
#include <morse.h>
#endif

#ifdef __cplusplus
extern "C"{
#endif
void matcov_comp_tile(
  double* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  struct tomo_struct *tomo, int part);

#ifdef USE_CHAMELEON
int APPNAME_matcov_tile( MORSE_desc_t *descData, MORSE_sequence_t *sequence, MORSE_request_t *request, struct tomo_struct *tomo, int type_mat);
#endif

#ifdef __cplusplus
}
#endif
#endif //MATCOV_H

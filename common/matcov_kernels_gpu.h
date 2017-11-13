/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef MATCOV_KERNEL_GPU_H
#define MATCOV_KERNEL_GPU_H

#include <cuda.h>
#include <cuda_runtime.h>
#include "moao_defs.h"

#ifdef __cplusplus
extern "C" {
#endif
void matcov_comp_tile_gpu_4(  char uplo, char copy, double* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  double *sspSizeL, long *Nssp, double *u, double *v,
  double *indexL0, double *cn2, int Nw, int Nlayer,
  long *Nsubap, long Nx, double lgs_cst, double *noiseNGS, double *noiseLGSxx,
  double *noiseLGSyy, double *noiseLGSxy, int type_mat, int nlgs, double DiamTel,cudaStream_t stream);
#ifdef __cplusplus
}
#endif
#endif //MATCOV_KERNEL_GPU_H

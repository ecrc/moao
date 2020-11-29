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

template<typename T>
void matcov_tile_gpu_4(  char uplo, char copy, T* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  T *sspSizeL, long *Nssp, T *u, T *v,
  T *indexL0, T *cn2, int Nw, int Nlayer,
  long *Nsubap, long Nx, T lgs_cst, T *noiseNGS, T *noiseLGSxx,
  T *noiseLGSyy, T *noiseLGSxy, int type_mat, int nlgs, T DiamTel,cudaStream_t stream);
#endif //MATCOV_KERNEL_GPU_H

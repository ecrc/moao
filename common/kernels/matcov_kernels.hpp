/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef MATCOV_KERNEL_H
#define MATCOV_KERNEL_H


template<typename T>
void matcov_kernel_4(
  char uplo, char copy, T* data, int xoffset, int yoffset, int lda,
  T *sspSizeL, long *Nssp, T *u, T *v,
  T *indexL0, T *cn2, int Nw, int Nlayer,
  long *Nsubap, long Nx, T lgs_cst, T *noiseNGS, T *noiseLGSxx, 
  T *noiseLGSyy, T *noiseLGSxy, int type_mat, int nlgs, T teldiam, int lx, int ly);

template<typename T>
T compute_element_tiled_4(
    int ipos, int jpos, T *sspSizeL, long *Nssp, T *u, T *v,
    T *indexL0, T *cn2, int Nw, int Nlayer,
    long * Nsubap_wfs, long Nx, T lgs_cst, T *noiseNGS, T *noiseLGSxx,
    T *noiseLGSyy, T *noiseLGSxy, int type_mat, int nlgs, T teldiam);

template<typename T>
void matcov_ts_kernel_tile(
  T* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  T *X, T *Y, long *Nssp, 
  T *indexL0, T *cn2, int Nw, int Nlayer, int Nsubap, T teldiam,
  int lx, int ly);


#define ROWMAJOR 1
#define COLMAJOR 0

#endif //MATCOV_KERNEL_H

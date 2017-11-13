/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef MATCOV_KERNEL_H
#define MATCOV_KERNEL_H

#ifdef __cplusplus
extern "C" {
#endif

void matcov_kernel_4(
  char uplo, char copy, double* data, int xoffset, int yoffset, int lda,
  double *sspSizeL, long *Nssp, double *u, double *v,
  double *indexL0, double *cn2, int Nw, int Nlayer,
  long *Nsubap, long Nx, double lgs_cst, double *noiseNGS, double *noiseLGSxx, 
  double *noiseLGSyy, double *noiseLGSxy, int type_mat, int nlgs, double teldiam, int lx, int ly);

double compute_element_tiled_4(
    int ipos, int jpos, double *sspSizeL, long *Nssp, double *u, double *v,
    double *indexL0, double *cn2, int Nw, int Nlayer,
    long * Nsubap_wfs, long Nx, double lgs_cst, double *noiseNGS, double *noiseLGSxx,
    double *noiseLGSyy, double *noiseLGSxy, int type_mat, int nlgs, double teldiam);

void matcov_ts_kernel_tile(
  double* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  double *X, double *Y, long *Nssp, 
  double *indexL0, double *cn2, int Nw, int Nlayer, int Nsubap, double teldiam,
  int lx, int ly);

#ifdef __cplusplus
}
#endif

#define ROWMAJOR 1
#define COLMAJOR 0

#endif //MATCOV_KERNEL_H

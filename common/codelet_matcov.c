/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "codelet_matcov.h"
#include <math.h>
#include "matcov_kernels.h"

#ifdef USE_GPU
#include "matcov_kernels_gpu.h"
//static void cl_matcov_cpu_func(void *descr[], void *cl_arg);
static void cl_matcov_cuda_func(void *descr[], void *cl_arg);
static struct starpu_codelet cl_CovMat =
{
    .where = /*STARPU_CPU|*/STARPU_CUDA,
    //.cpu_funcs = {cl_matcov_cpu_func},
    .cuda_funcs = {cl_matcov_cuda_func},
    .nbuffers = 14
};
#else
static void cl_matcov_cpu_func(void *descr[], void *cl_arg);
static struct starpu_codelet cl_CovMat =
{
    .where = STARPU_CPU,
    .cpu_funcs = {cl_matcov_cpu_func},
    .nbuffers = 14
};
#endif

void APPNAME_TASK_matcov(APPNAME_desc *data,int nrows,int ncols, int m, int n,
                    int lda, struct tomo_struct tomo,
MORSE_desc_t *descSspSizeL,
MORSE_desc_t *descNssp,
MORSE_desc_t *descU,
MORSE_desc_t *descV,
MORSE_desc_t *descX,
MORSE_desc_t *descY,
MORSE_desc_t *descL0diff,
MORSE_desc_t *descCn2,
MORSE_desc_t *descNsubap,
MORSE_desc_t *descNoiseNGS,
MORSE_desc_t *descNoiseLGSxx,
MORSE_desc_t *descNoiseLGSyy,
MORSE_desc_t *descNoiseLGSxy,
                    int type_mat)
{
  char uplo, copy;
  int xoffset=m*data->mb;
  int yoffset=n*data->nb;

  uplo = 'f';   // full generation is enabled by default
  copy = 'n';

  if(type_mat == 1) // Caa matrix
  {
    // check if a square diagonal tile is generated then we set uplo to 'l' or 'u'
    // and then enable the copy
    // This also applies if the entire matrix will be generated
    // otherwise (off diagonal tile or non square submatrix) - full generation is assumed
    if((xoffset == yoffset) && (nrows == ncols))    // if sqaure & diagonal
    {
      uplo = 'l';
      copy = 'n';
    }
    else    // full generation, copy is ignored
    {
      uplo = 'l';
      copy = 'n';
    }
  }
  //else if(type_mat == 2) //
  else
  if(type_mat == 2 || type_mat == 3) // Cmaa matrix
  {
    uplo = 'f';             // full generation, copy is ignored
  }
  else
  if(type_mat != 4)
  {
    fprintf(stderr, "ERROR: unrecognized type_mat %d \n", type_mat);
  }

  int Nw=(int)tomo.Nw;

  struct starpu_codelet *codelet=&cl_CovMat;

  starpu_insert_task(codelet,
      STARPU_VALUE,&uplo            ,sizeof(char),
      STARPU_VALUE,&copy            ,sizeof(char),
      STARPU_W   ,RTBLKADDR(data,m,n),
      STARPU_VALUE,&nrows           ,sizeof(int),
      STARPU_VALUE,&ncols           ,sizeof(int),
      STARPU_VALUE,&xoffset         ,sizeof(int),
      STARPU_VALUE,&yoffset         ,sizeof(int),
      STARPU_VALUE,&lda             ,sizeof(int),
      STARPU_R    ,RTBLKADDR(descSspSizeL,0,0),
      STARPU_R    ,RTBLKADDR(descNssp,0,0),
      STARPU_R    ,RTBLKADDR(descU,0,0),
      STARPU_R    ,RTBLKADDR(descV,0,0),
      STARPU_R    ,RTBLKADDR(descX,0,0),
      STARPU_R    ,RTBLKADDR(descY,0,0),
      STARPU_R    ,RTBLKADDR(descL0diff,0,0),
      STARPU_R    ,RTBLKADDR(descCn2,0,0),
      STARPU_VALUE,&Nw              ,sizeof(int),
      STARPU_VALUE,&(tomo.Nlayer)   ,sizeof(long),
      STARPU_R    ,RTBLKADDR(descNsubap,0,0),
      STARPU_VALUE,&(tomo.Nx)       ,sizeof(int),
      STARPU_VALUE,&(tomo.lgs_cst)  ,sizeof(double),
      STARPU_R    ,RTBLKADDR(descNoiseNGS,0,0),
      STARPU_R    ,RTBLKADDR(descNoiseLGSxx,0,0),
      STARPU_R    ,RTBLKADDR(descNoiseLGSyy,0,0),
      STARPU_R    ,RTBLKADDR(descNoiseLGSxy,0,0),
      STARPU_VALUE,      &type_mat  ,sizeof(int),
      STARPU_VALUE,&(tomo.nlgs)     ,sizeof(int),
      STARPU_VALUE,&(tomo.DiamTel)  ,sizeof(double),
      0);
}






static void cl_matcov_cpu_func(void *descr[], void *cl_arg){

  double *data,*sspSizeL,*u,*v,*X,*Y,*L0diff,*cn2;
  double *noiseNGS, *noiseLGSxx, *noiseLGSyy, *noiseLGSxy;
  long *Nssp,*Nsubap;

  data      = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
  sspSizeL  = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
  Nssp      = (long   *)STARPU_MATRIX_GET_PTR(descr[2]);
  u         = (double *)STARPU_MATRIX_GET_PTR(descr[3]);
  v         = (double *)STARPU_MATRIX_GET_PTR(descr[4]);
  X         = (double *)STARPU_MATRIX_GET_PTR(descr[5]);
  Y         = (double *)STARPU_MATRIX_GET_PTR(descr[6]);
  L0diff    = (double *)STARPU_MATRIX_GET_PTR(descr[7]);
  cn2       = (double *)STARPU_MATRIX_GET_PTR(descr[8]);
  Nsubap    = (long   *)STARPU_MATRIX_GET_PTR(descr[9]);
  noiseNGS  = (double *)STARPU_MATRIX_GET_PTR(descr[10]);
  noiseLGSxx= (double *)STARPU_MATRIX_GET_PTR(descr[11]);
  noiseLGSyy= (double *)STARPU_MATRIX_GET_PTR(descr[12]);
  noiseLGSxy= (double *)STARPU_MATRIX_GET_PTR(descr[13]);

  double lgs_cst,DiamTel; 
  int Nw,nrows,ncols,xoffset,yoffset,lda,Nx,type_mat,nlgs;
  long Nlayer;
  char uplo,copy;

  starpu_codelet_unpack_args(cl_arg, &uplo,
                                      &copy,
                                      &nrows,
                                      &ncols,
                                      &xoffset,
                                      &yoffset,
                                      &lda,
                                      &Nw,
                                      &Nlayer,
                                      &Nx,
                                      &lgs_cst,
                                      &type_mat,
                                      &nlgs,
                                      &DiamTel
                                      );

  if(type_mat == 4){

    int lx,ly;
    for(lx = 0; lx < nrows; lx++){
      for(ly = 0; ly < ncols; ly++){
        matcov_ts_kernel_tile(
            data, nrows, ncols, xoffset, yoffset, lda,
            X, Y, Nssp,
            L0diff, cn2,
            Nw, Nlayer, Nsubap[Nw-1], DiamTel,
            lx, ly);
      }
    }
  }else{
    int lx,ly;
    for(lx = 0; lx < nrows; lx++){
      for(ly = 0; ly < ncols; ly++){
        matcov_kernel_4(
            uplo, copy, data, xoffset, yoffset, lda,  sspSizeL,
            Nssp, u, v, L0diff, cn2, Nw, Nlayer,
            Nsubap, Nx, lgs_cst, noiseNGS, noiseLGSxx,
            noiseLGSyy, noiseLGSxy, type_mat, nlgs, DiamTel,
            lx, ly);
      }
    }
  }

}

#ifdef USE_GPU
static void cl_matcov_cuda_func(void *descr[], void *cl_arg){

  double *data,*sspSizeL,*u,*v,*X,*Y,*L0diff,*cn2;
  double *noiseNGS, *noiseLGSxx, *noiseLGSyy, *noiseLGSxy;
  long *Nssp,*Nsubap;

  data      = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
  sspSizeL  = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
  Nssp      = (long   *)STARPU_MATRIX_GET_PTR(descr[2]);
  u         = (double *)STARPU_MATRIX_GET_PTR(descr[3]);
  v         = (double *)STARPU_MATRIX_GET_PTR(descr[4]);
  X         = (double *)STARPU_MATRIX_GET_PTR(descr[5]);
  Y         = (double *)STARPU_MATRIX_GET_PTR(descr[6]);
  L0diff    = (double *)STARPU_MATRIX_GET_PTR(descr[7]);
  cn2       = (double *)STARPU_MATRIX_GET_PTR(descr[8]);
  Nsubap    = (long   *)STARPU_MATRIX_GET_PTR(descr[9]);
  noiseNGS  = (double *)STARPU_MATRIX_GET_PTR(descr[10]);
  noiseLGSxx= (double *)STARPU_MATRIX_GET_PTR(descr[11]);
  noiseLGSyy= (double *)STARPU_MATRIX_GET_PTR(descr[12]);
  noiseLGSxy= (double *)STARPU_MATRIX_GET_PTR(descr[13]);

  double lgs_cst,DiamTel; 
  int Nw,nrows,ncols,xoffset,yoffset,lda,Nx,type_mat,nlgs;
  long Nlayer;
  char uplo,copy;

  starpu_codelet_unpack_args(cl_arg, &uplo,
                                      &copy,
                                      &nrows,
                                      &ncols,
                                      &xoffset,
                                      &yoffset,
                                      &lda,
                                      &Nw,
                                      &Nlayer,
                                      &Nx,
                                      &lgs_cst,
                                      &type_mat,
                                      &nlgs,
                                      &DiamTel
                                      );

  cudaStream_t stream = starpu_cuda_get_local_stream();

  matcov_comp_tile_gpu_4(
            uplo, copy, data, nrows, ncols, xoffset, yoffset, lda, sspSizeL,
            Nssp, u, v,
            L0diff, cn2, Nw, Nlayer,
            Nsubap, Nx, lgs_cst,noiseNGS, noiseLGSxx,
            noiseLGSyy, noiseLGSxy, type_mat, nlgs, DiamTel, stream);


#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif
}
#endif



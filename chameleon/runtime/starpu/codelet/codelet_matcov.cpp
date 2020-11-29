/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "codelet_matcov.hpp"
#include <math.h>

#ifdef USE_GPU
//#include "matcov_kernels_gpu.h"
static struct starpu_codelet cl_sCovMat =
{
    where : STARPU_CUDA,
    can_execute : NULL,     //int (*can_execute)(unsigned workerid, struct starpu_task *task, unsigned nimpl);
    type : STARPU_SEQ,      //enum starpu_codelet_type type;
    max_parallelism : 0,    //int max_parallelism;

    cpu_func : {},          //starpu_cpu_func_t cpu_func STARPU_DEPRECATED;
    cuda_func:{},           //starpu_cuda_func_t cuda_func STARPU_DEPRECATED;
    opencl_func:{},         //starpu_opencl_func_t opencl_func STARPU_DEPRECATED;

    cpu_funcs : {}, //starpu_cpu_func_t cpu_funcs[STARPU_MAXIMPLEMENTATIONS];
    cuda_funcs: {cl_matcov_cuda_func<float> },         //starpu_cuda_func_t cuda_funcs[STARPU_MAXIMPLEMENTATIONS];
    cuda_flags: {},         //char cuda_flags[STARPU_MAXIMPLEMENTATIONS];
    opencl_funcs : {},      //starpu_opencl_func_t opencl_funcs[STARPU_MAXIMPLEMENTATIONS];
    opencl_flags : {},      //char opencl_flags[STARPU_MAXIMPLEMENTATIONS];
    mic_funcs : {},         //starpu_mic_func_t mic_funcs[STARPU_MAXIMPLEMENTATIONS];
    //mpi_ms_funcs :{},       //starpu_mpi_ms_func_t mpi_ms_funcs[STARPU_MAXIMPLEMENTATIONS];
    scc_funcs : {},         //starpu_scc_func_t scc_funcs[STARPU_MAXIMPLEMENTATIONS];
    cpu_funcs_name :{},     //const char *cpu_funcs_name[STARPU_MAXIMPLEMENTATIONS];
    nbuffers : 14           //int nbuffers;
};
static struct starpu_codelet cl_dCovMat =
{
    where : STARPU_CUDA,
    can_execute : NULL,     //int (*can_execute)(unsigned workerid, struct starpu_task *task, unsigned nimpl);
    type : STARPU_SEQ,      //enum starpu_codelet_type type;
    max_parallelism : 0,    //int max_parallelism;

    cpu_func : {},          //starpu_cpu_func_t cpu_func STARPU_DEPRECATED;
    cuda_func:{},           //starpu_cuda_func_t cuda_func STARPU_DEPRECATED;
    opencl_func:{},         //starpu_opencl_func_t opencl_func STARPU_DEPRECATED;

    cpu_funcs : {}, //starpu_cpu_func_t cpu_funcs[STARPU_MAXIMPLEMENTATIONS];
    cuda_funcs: {cl_matcov_cuda_func<double> },         //starpu_cuda_func_t cuda_funcs[STARPU_MAXIMPLEMENTATIONS];
    cuda_flags: {},         //char cuda_flags[STARPU_MAXIMPLEMENTATIONS];
    opencl_funcs : {},      //starpu_opencl_func_t opencl_funcs[STARPU_MAXIMPLEMENTATIONS];
    opencl_flags : {},      //char opencl_flags[STARPU_MAXIMPLEMENTATIONS];
    mic_funcs : {},         //starpu_mic_func_t mic_funcs[STARPU_MAXIMPLEMENTATIONS];
    //mpi_ms_funcs :{},       //starpu_mpi_ms_func_t mpi_ms_funcs[STARPU_MAXIMPLEMENTATIONS];
    scc_funcs : {},         //starpu_scc_func_t scc_funcs[STARPU_MAXIMPLEMENTATIONS];
    cpu_funcs_name :{},     //const char *cpu_funcs_name[STARPU_MAXIMPLEMENTATIONS];
    nbuffers : 14           //int nbuffers;
};
#else
static struct starpu_codelet cl_sCovMat =
{
    where : STARPU_CPU,
    can_execute : NULL,     //int (*can_execute)(unsigned workerid, struct starpu_task *task, unsigned nimpl);
    type : STARPU_SEQ,      //enum starpu_codelet_type type;
    max_parallelism : 0,    //int max_parallelism;

    cpu_func : {},          //starpu_cpu_func_t cpu_func STARPU_DEPRECATED;
    cuda_func:{},           //starpu_cuda_func_t cuda_func STARPU_DEPRECATED;
    opencl_func:{},         //starpu_opencl_func_t opencl_func STARPU_DEPRECATED;

    cpu_funcs : {cl_matcov_cpu_func<float>}, //starpu_cpu_func_t cpu_funcs[STARPU_MAXIMPLEMENTATIONS];
    cuda_funcs: {},         //starpu_cuda_func_t cuda_funcs[STARPU_MAXIMPLEMENTATIONS];
    cuda_flags: {},         //char cuda_flags[STARPU_MAXIMPLEMENTATIONS];
    opencl_funcs : {},      //starpu_opencl_func_t opencl_funcs[STARPU_MAXIMPLEMENTATIONS];
    opencl_flags : {},      //char opencl_flags[STARPU_MAXIMPLEMENTATIONS];
    mic_funcs : {},         //starpu_mic_func_t mic_funcs[STARPU_MAXIMPLEMENTATIONS];
    //mpi_ms_funcs :{},       //starpu_mpi_ms_func_t mpi_ms_funcs[STARPU_MAXIMPLEMENTATIONS];
    scc_funcs : {},         //starpu_scc_func_t scc_funcs[STARPU_MAXIMPLEMENTATIONS];
    cpu_funcs_name :{},     //const char *cpu_funcs_name[STARPU_MAXIMPLEMENTATIONS];
    nbuffers : 14//int nbuffers;
};
static struct starpu_codelet cl_dCovMat =
{
    where : STARPU_CPU,
    can_execute : NULL,     //int (*can_execute)(unsigned workerid, struct starpu_task *task, unsigned nimpl);
    type : STARPU_SEQ,      //enum starpu_codelet_type type;
    max_parallelism : 0,    //int max_parallelism;

    cpu_func : {},          //starpu_cpu_func_t cpu_func STARPU_DEPRECATED;
    cuda_func:{},           //starpu_cuda_func_t cuda_func STARPU_DEPRECATED;
    opencl_func:{},         //starpu_opencl_func_t opencl_func STARPU_DEPRECATED;

    cpu_funcs : {cl_matcov_cpu_func<double>}, //starpu_cpu_func_t cpu_funcs[STARPU_MAXIMPLEMENTATIONS];
    cuda_funcs: {},         //starpu_cuda_func_t cuda_funcs[STARPU_MAXIMPLEMENTATIONS];
    cuda_flags: {},         //char cuda_flags[STARPU_MAXIMPLEMENTATIONS];
    opencl_funcs : {},      //starpu_opencl_func_t opencl_funcs[STARPU_MAXIMPLEMENTATIONS];
    opencl_flags : {},      //char opencl_flags[STARPU_MAXIMPLEMENTATIONS];
    mic_funcs : {},         //starpu_mic_func_t mic_funcs[STARPU_MAXIMPLEMENTATIONS];
    //mpi_ms_funcs :{},       //starpu_mpi_ms_func_t mpi_ms_funcs[STARPU_MAXIMPLEMENTATIONS];
    scc_funcs : {},         //starpu_scc_func_t scc_funcs[STARPU_MAXIMPLEMENTATIONS];
    cpu_funcs_name :{},     //const char *cpu_funcs_name[STARPU_MAXIMPLEMENTATIONS];
    nbuffers : 14//int nbuffers;
};
#endif


template<typename T>
void MOAO_STARPU::TASK_matcov(APPNAME_desc *data,int nrows,int ncols, int m, int n,
                    int lda, Tomo_struct<T> tomo,
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

  int Nw=(int)tomo.sys.nW;

  struct starpu_codelet *codelet ;//=&cl_CovMat;
  if(std::is_same<T,float>::value){
    codelet=&cl_sCovMat;
  }
  else{
    codelet=&cl_dCovMat;
  }

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
      STARPU_VALUE,&(tomo.atm.nLayer)   ,sizeof(long),
      STARPU_R    ,RTBLKADDR(descNsubap,0,0),
      STARPU_VALUE,&(tomo.Nx)       ,sizeof(int),
      STARPU_VALUE,&(tomo.sys.lgsCst)  ,sizeof(T),
      STARPU_R    ,RTBLKADDR(descNoiseNGS,0,0),
      STARPU_R    ,RTBLKADDR(descNoiseLGSxx,0,0),
      STARPU_R    ,RTBLKADDR(descNoiseLGSyy,0,0),
      STARPU_R    ,RTBLKADDR(descNoiseLGSxy,0,0),
      STARPU_VALUE,      &type_mat  ,sizeof(int),
      STARPU_VALUE,&(tomo.sys.nLgs)     ,sizeof(int),
      STARPU_VALUE,&(tomo.sys.diam)  ,sizeof(T),
      0);
}
template void MOAO_STARPU::TASK_matcov(APPNAME_desc *data,int nrows,int ncols, int m, int n,
                    int lda, Tomo_struct<float> tomo,
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
int type_mat);
template void MOAO_STARPU::TASK_matcov(APPNAME_desc *data,int nrows,int ncols, int m, int n,
                    int lda, Tomo_struct<double> tomo,
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
int type_mat);






template<typename T>
static void cl_matcov_cpu_func(void *descr[], void *cl_arg){

  T *data,*sspSizeL,*u,*v,*X,*Y,*L0diff,*cn2;
  T *noiseNGS, *noiseLGSxx, *noiseLGSyy, *noiseLGSxy;
  long *Nssp,*Nsubap;

  data      = (T *)STARPU_MATRIX_GET_PTR(descr[0]);
  sspSizeL  = (T *)STARPU_MATRIX_GET_PTR(descr[1]);
  Nssp      = (long   *)STARPU_MATRIX_GET_PTR(descr[2]);
  u         = (T *)STARPU_MATRIX_GET_PTR(descr[3]);
  v         = (T *)STARPU_MATRIX_GET_PTR(descr[4]);
  X         = (T *)STARPU_MATRIX_GET_PTR(descr[5]);
  Y         = (T *)STARPU_MATRIX_GET_PTR(descr[6]);
  L0diff    = (T *)STARPU_MATRIX_GET_PTR(descr[7]);
  cn2       = (T *)STARPU_MATRIX_GET_PTR(descr[8]);
  Nsubap    = (long   *)STARPU_MATRIX_GET_PTR(descr[9]);
  noiseNGS  = (T *)STARPU_MATRIX_GET_PTR(descr[10]);
  noiseLGSxx= (T *)STARPU_MATRIX_GET_PTR(descr[11]);
  noiseLGSyy= (T *)STARPU_MATRIX_GET_PTR(descr[12]);
  noiseLGSxy= (T *)STARPU_MATRIX_GET_PTR(descr[13]);

  T lgs_cst,DiamTel; 
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
template static void cl_matcov_cpu_func<float>(void *descr[], void *cl_arg);
template static void cl_matcov_cpu_func<double>(void *descr[], void *cl_arg);

#ifdef USE_GPU
template<typename T>
static void cl_matcov_cuda_func(void *descr[], void *cl_arg){

  T *data,*sspSizeL,*u,*v,*X,*Y,*L0diff,*cn2;
  T *noiseNGS, *noiseLGSxx, *noiseLGSyy, *noiseLGSxy;
  long *Nssp,*Nsubap;

  data      = (T *)STARPU_MATRIX_GET_PTR(descr[0]);
  sspSizeL  = (T *)STARPU_MATRIX_GET_PTR(descr[1]);
  Nssp      = (long   *)STARPU_MATRIX_GET_PTR(descr[2]);
  u         = (T *)STARPU_MATRIX_GET_PTR(descr[3]);
  v         = (T *)STARPU_MATRIX_GET_PTR(descr[4]);
  X         = (T *)STARPU_MATRIX_GET_PTR(descr[5]);
  Y         = (T *)STARPU_MATRIX_GET_PTR(descr[6]);
  L0diff    = (T *)STARPU_MATRIX_GET_PTR(descr[7]);
  cn2       = (T *)STARPU_MATRIX_GET_PTR(descr[8]);
  Nsubap    = (long   *)STARPU_MATRIX_GET_PTR(descr[9]);
  noiseNGS  = (T *)STARPU_MATRIX_GET_PTR(descr[10]);
  noiseLGSxx= (T *)STARPU_MATRIX_GET_PTR(descr[11]);
  noiseLGSyy= (T *)STARPU_MATRIX_GET_PTR(descr[12]);
  noiseLGSxy= (T *)STARPU_MATRIX_GET_PTR(descr[13]);

  T lgs_cst,DiamTel; 
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

  matcov_tile_gpu_4(
            uplo, copy, data, nrows, ncols, xoffset, yoffset, lda, sspSizeL,
            Nssp, u, v,
            L0diff, cn2, Nw, Nlayer,
            Nsubap, Nx, lgs_cst,noiseNGS, noiseLGSxx,
            noiseLGSyy, noiseLGSxy, type_mat, nlgs, DiamTel, stream);


#ifndef STARPU_CUDA_ASYNC
    cudaStreamSynchronize( stream );
#endif
}
template static void cl_matcov_cuda_func<float>(void *descr[], void *cl_arg);
template static void cl_matcov_cuda_func<double>(void *descr[], void *cl_arg);
#endif



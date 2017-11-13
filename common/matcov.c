/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "matcov.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

//#include <cuda.h>
//#include <cuda_runtime.h>
//#include <cuda_profiler_api.h>

#include "utils.h"

#ifdef USE_CHAMELEON
#include "codelet_matcov.h"

/*! matcov tile driver
*  Arguments
*  ==========
*  data         double pointer: A pointer to the matrix/submatrix to be generated. It
*                       should always point to the first element in a matrix/submatrix
*
*  nrows        integer: The number of rows of the matrix/submatrix to be generated
*
*  ncols        integer: The number of columns of the matrix/submatrix to be generated
*
*  xoffset      integer: The x-offset of the submatrix, must be zero if the entire matrix
*                       is generated. Its the x-coordinate of the first element in the matrix/submatrix
*
*  yoffset  integer: The y-offset of the submatrix, must be zero if the entire matrix
*                       is generated. Its the y-coordinate of the first element in the matrix/submatrix
*
*  lda          integer: The leading dimension of the matrix/submatrix
*
*  rowMajor     integer: specify if the tile is generated in row major or column major layout
*                        0:column major, else row major (default=1)
*/
 int APPNAME_matcov_tile(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t *request, struct tomo_struct *tomo, int type_mat)
{
MORSE_desc_t *descSspSizeL=NULL;
MORSE_desc_t *descNssp=NULL;
MORSE_desc_t *descU=NULL;
MORSE_desc_t *descV=NULL;
MORSE_desc_t *descX=NULL;
MORSE_desc_t *descY=NULL;
MORSE_desc_t *descL0diff=NULL;
MORSE_desc_t *descCn2=NULL;
MORSE_desc_t *descNsubap=NULL;
MORSE_desc_t *descAlphaX=NULL;
MORSE_desc_t *descAlphaY=NULL;
MORSE_desc_t *descNoiseNGS;
MORSE_desc_t *descNoiseLGSxx;
MORSE_desc_t *descNoiseLGSyy;
MORSE_desc_t *descNoiseLGSxy;


int Nw=(int)tomo->Nw, Nl=(int)tomo->Nlayer,Nx=(int)tomo->Nx;
  int i;
  int Nsubap=0;
  for(i=0;i<tomo->nlgs;i++){
  Nsubap+=tomo->Nsubap[i];
  }
  if(!Nsubap)
    Nsubap=1;

MORSE_Desc_Create(&(descSspSizeL) , tomo->sspSizeL, MorseRealDouble, Nw*Nl, 1, Nw*Nl, Nw*Nl, 1, 0, 0, Nw*Nl, 1,1,1);
MORSE_Desc_Create(&(descNssp), tomo->Nssp, MorseRealDouble, Nw, 1, Nw,Nw , 1, 0, 0, Nw, 1,1,1);
MORSE_Desc_Create(&(descU), tomo->u , MorseRealDouble, Nl*Nx, 1, Nl*Nx, Nl*Nx, 1, 0, 0, Nl*Nx, 1,1,1);
MORSE_Desc_Create(&(descV), tomo->v , MorseRealDouble, Nl*Nx, 1, Nl*Nx, Nl*Nx, 1, 0, 0, Nl*Nx, 1,1,1);
MORSE_Desc_Create(&(descX), tomo->X,MorseRealDouble,Nx,1,Nx,Nx,1,0,0,Nx,1,1,1);
MORSE_Desc_Create(&(descY), tomo->Y,MorseRealDouble,Nx,1,Nx,Nx,1,0,0,Nx,1,1,1);
MORSE_Desc_Create(&(descL0diff),tomo->L0diff,MorseRealDouble,Nl, 1,Nl,Nl, 1, 0, 0,Nl, 1,1,1);
MORSE_Desc_Create(&(descCn2), tomo->cn2,     MorseRealDouble, Nl, 1, Nl , Nl , 1, 0, 0, Nl,1,1,1);
MORSE_Desc_Create(&(descNsubap), tomo->Nsubap, MorseRealDouble, Nw, 1, Nw,Nw , 1, 0, 0, Nw, 1,1,1);
MORSE_Desc_Create(&(descAlphaX), tomo->alphaX, MorseRealDouble, Nw, 1, Nw,Nw , 1, 0, 0, Nw, 1,1,1);
MORSE_Desc_Create(&(descAlphaY), tomo->alphaY, MorseRealDouble, Nw, 1, Nw,Nw , 1, 0, 0, Nw, 1,1,1);
MORSE_Desc_Create(&(descNoiseNGS)  , tomo->noiseNGS, MorseRealDouble, Nw, 1, Nw,Nw , 1, 0, 0, Nw, 1,1,1);
MORSE_Desc_Create(&(descNoiseLGSxx), tomo->noiseLGSxx, MorseRealDouble, Nsubap, 1, Nsubap, Nsubap, 1, 0, 0, Nsubap, 1,1,1);
MORSE_Desc_Create(&(descNoiseLGSyy), tomo->noiseLGSyy, MorseRealDouble, Nsubap, 1, Nsubap, Nsubap, 1, 0, 0, Nsubap, 1,1,1);
MORSE_Desc_Create(&(descNoiseLGSxy), tomo->noiseLGSxy, MorseRealDouble, Nsubap, 1, Nsubap, Nsubap, 1, 0, 0, Nsubap, 1,1,1);

    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return sequence->status;
    RUNTIME_options_init(&options, morse, sequence, request);

    int m,n,ldam;
    int tempmm,tempnn;
    int ret=starpu_init(NULL);

    for(m=0;m<descA->mt;m++){
        tempmm = m == descA->mt-1 ? descA->m-m*descA->mb : descA->mb;
        ldam=descA->get_blkldd(descA, m);
        for(n=0;n<descA->nt;n++){
            tempnn = n == descA->nt-1 ? descA->n-n*descA->nb : descA->nb;
            APPNAME_TASK_matcov(descA,tempmm,tempnn,m,n,ldam,*tomo,
            descSspSizeL,
            descNssp,
            descU,
            descV,
            descX,
            descY,
            descL0diff,
            descCn2,
            descNsubap,
            descNoiseNGS,
            descNoiseLGSxx,
            descNoiseLGSyy,
            descNoiseLGSxy,
            type_mat);
        }
    }

    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();


MORSE_Desc_Destroy(&descSspSizeL);
MORSE_Desc_Destroy(&descNssp);
MORSE_Desc_Destroy(&descU);
MORSE_Desc_Destroy(&descV);
MORSE_Desc_Destroy(&descX);
MORSE_Desc_Destroy(&descY);
MORSE_Desc_Destroy(&descL0diff);
MORSE_Desc_Destroy(&descCn2);
MORSE_Desc_Destroy(&descNsubap);
MORSE_Desc_Destroy(&descAlphaX);
MORSE_Desc_Destroy(&descAlphaY);
MORSE_Desc_Destroy(&descNoiseNGS);
MORSE_Desc_Destroy(&descNoiseLGSxx);
MORSE_Desc_Destroy(&descNoiseLGSyy);
MORSE_Desc_Destroy(&descNoiseLGSxy);

    return 0;
}
#endif //USE_CHAMELEON

void matcov_comp_tile(
  double* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  struct tomo_struct *tomo, int type_mat)
{

  char uplo, copy;

  uplo = 'f';   // full generation is enabled by default
  copy = 'n';

  //int type_mat = tomo->part;

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
    exit(1);
  }

  if(type_mat == 4){
    const long Nw = tomo->Nw;
    const long Nsubap = tomo->Nsubap[Nw-1];

    int lx,ly;
    #ifdef USE_OPENMP
    #pragma omp parallel private(ly)
    #pragma omp for nowait
    #endif
    for(lx = 0; lx < nrows; lx++){
      for(ly = 0; ly < ncols; ly++){
        matcov_ts_kernel_tile(
            data, nrows, ncols, xoffset, yoffset, lda,
            tomo->X, tomo->Y, tomo->Nssp,
            tomo->L0diff, tomo->cn2,
            tomo->Nw, tomo->Nlayer, Nsubap, tomo->DiamTel,
            lx, ly);
      }
    }
  }else{
    //const long Nsubap = tomo->Nsubap[0];

    int lx,ly;
    #ifdef USE_OPENMP
    #pragma omp parallel private(ly)
    #pragma omp for nowait
    #endif
    for(lx = 0; lx < nrows; lx++){
      for(ly = 0; ly < ncols; ly++){
        matcov_kernel_4(
            uplo, copy, data, xoffset, yoffset, lda, tomo->sspSizeL,
            tomo->Nssp, tomo->u, tomo->v,
            tomo->L0diff, tomo->cn2, tomo->Nw, tomo->Nlayer,
            tomo->Nsubap, tomo->Nx, tomo->lgs_cst, 
            tomo->noiseNGS, tomo->noiseLGSxx, tomo->noiseLGSyy, tomo-> noiseLGSxy, 
            type_mat, tomo->nlgs, tomo->DiamTel, lx, ly);
      }
    }
  }
}

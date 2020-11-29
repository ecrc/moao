/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "matcov.hpp"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "kernels/matcov_kernels.hpp"


template<typename T>
void MOAO_COMMON::matcov_tile(
  T* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  Tomo_struct<T> *tomo, int type_mat)
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
    const long Nw = tomo->sys.nW;
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
            tomo->X.data(), tomo->Y.data(), tomo->Nssp.data(),
            tomo->L0diff.data(), tomo->atm.cn2.data(),
            tomo->sys.nW, tomo->atm.nLayer, Nsubap, tomo->sys.diam,
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
            uplo, copy, data, xoffset, yoffset, lda, tomo->sspSizeL.data(),
            tomo->Nssp.data(), tomo->u.data(), tomo->v.data(),
            tomo->L0diff.data(), tomo->atm.cn2.data(), tomo->sys.nW, tomo->atm.nLayer,
            tomo->Nsubap.data(), tomo->Nx, tomo->sys.lgsCst, 
            tomo->noiseNGS.data(), tomo->noiseLGSxx.data(), tomo->noiseLGSyy.data(), tomo-> noiseLGSxy.data(), 
            type_mat, tomo->sys.nLgs, tomo->sys.diam, lx, ly);
      }
    }
  }
}

template void MOAO_COMMON::matcov_tile( float* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  Tomo_struct<float> *tomo, int type_mat);

template void MOAO_COMMON::matcov_tile( double* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  Tomo_struct<double> *tomo, int type_mat);

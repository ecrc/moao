/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef CHECK_UTILS_H_
#define CHECK_UTILS_H_

//#define SAFEMALLOC 1
#include <stdio.h>
#include <string.h>
//#include "../common/matcov_tiled.h"

//#ifdef USE_PLASMA
//#include <plasma.h>
//#include <descriptor.h>
//#endif

//#include "../../common/moao_defs.h"


template<typename T, typename U>
void compareMatrices(T *A, U *B, long nrows, long ncols, char uplo, char trans, T &maxErr, T &frobeniusN );
//double compareMatrices2(real_t *A, double *B, long nrows, long ncols, char uplo, char trans, real_t *err, real_t *vA, real_t *vB);

//struct TileInfo{
//    //coordinates of the tile (in the tile grid
//    int m;
//    int n;
//    //coordinate of the first element of the tile in the global matrix
//    int i;
//    int j;
//    //dimensions of the tile
//    int nrows;
//    int ncols;
//    //tile address in tile representation
//    void *addr;
//};
//
//#ifdef USE_PLASMA
//#define HANDLE_ERROR(_result, _func)               \
//  if(PLASMA_SUCCESS != _result) {                  \
//    fprintf(stderr,"An error occurred (%d), when calling %s at line: %d... \n\n", _result, _func, __LINE__ -1);   \
//    break;                                         \
//  }
//#define HANDLE_ERROR_RET(_result, _func)           \
//  if(PLASMA_SUCCESS != _result) {                  \
//    fprintf(stderr,"An error occurred (%d), when calling %s at line: %d... \n\n", _result, _func, __LINE__ -1);   \
//    exit(-1);                                      \
//  }
//#endif



#endif

/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef CHECK_UTILS_H_
#define CHECK_UTILS_H_

#define SAFEMALLOC 1
#include <stdio.h>
#include <string.h>
//#include "../common/matcov_tiled.h"

#ifdef USE_PLASMA
#include <plasma.h>
#include <descriptor.h>
#endif

#include "moao_defs.h"


void *myMalloc(size_t size);
void *myCalloc(size_t num, size_t size);

struct TileInfo{
    //coordinates of the tile (in the tile grid
    int m;
    int n;
    //coordinate of the first element of the tile in the global matrix
    int i;
    int j;
    //dimensions of the tile
    int nrows;
    int ncols;
    //tile address in tile representation
    void *addr;
};

#ifdef USE_PLASMA
#define HANDLE_ERROR(_result, _func)               \
  if(PLASMA_SUCCESS != _result) {                  \
    fprintf(stderr,"An error occurred (%d), when calling %s at line: %d... \n\n", _result, _func, __LINE__ -1);   \
    break;                                         \
  }
#define HANDLE_ERROR_RET(_result, _func)           \
  if(PLASMA_SUCCESS != _result) {                  \
    fprintf(stderr,"An error occurred (%d), when calling %s at line: %d... \n\n", _result, _func, __LINE__ -1);   \
    exit(-1);                                      \
  }
real_t p_compareMatrices(real_t *A, PLASMA_desc descB, /*real_t *B,*/ char uplo, char trans);

//void p_get_TileInfo(PLASMA_desc descA, struct TileInfo ti, int m, int n);
void p_tile_to_lapack(PLASMA_desc descA, real_t *B, int rowMajor);
void p_lapack_to_tile(PLASMA_desc descA, real_t *B, int rowMajor);
#endif

//void get_thread_count(int *thrdnbr);

double compareMatrices(double *A, double *B, long nrows, long ncols, char uplo, char trans, double *err, double *vA, double *vB);
double compareMatrices2(real_t *A, double *B, long nrows, long ncols, char uplo, char trans, real_t *err, real_t *vA, real_t *vB);
void minMax(real_t *A,long nbElem,real_t *min,real_t*max);
void copy(double *src, real_t *dst,long nbElem, long ncols, char uplo,char trans);
void saveMat(char *fileName, real_t *mat, int nrows, int ncols);
//void getFitsDims(char *fileName,long *naxes);
//void concatenateFileName(char *fileName, char *filePath, char *typeMat, char *suffix);
//void readFits(char *fileName,real_t *Cmm);
//void writeFits(char *file_name, int nx, int ny, real_t* data);

//void printTomoStruct(struct tomo_struct tomo);

#endif

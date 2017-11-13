/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>

#include "moao_defs.h"

#ifdef __cplusplus
extern "C" {
#endif
void read_paraml(FILE *file, long *par);
void read_parami(FILE *file, int *par);
void read_paramd(FILE *file, real_t *par);
void read_arrayl(FILE *file, long *arr, int nmax);
void read_arrayi(FILE *file, int *arr, int nmax);
void read_arrayd(FILE *file, real_t *arr, int nmax);
double cWtime(void);
#ifdef __cplusplus
}
#endif

#ifdef DEBUG_MSG
#define FPRINTF fprintf
#else
#define FPRINTF 
#endif

#ifdef USE_CHAMELEON
#define HANDLE_ERROR(_result, _func)               \
  if(MORSE_SUCCESS != _result) {                  \
    fprintf(stderr,"An error occurred (%d), when calling %s at line: %d... \n\n", _result, _func, __LINE__ -1);   \
  }  //break;                                         
  
#define HANDLE_ERROR_RET(_result, _func)               \
  if(MORSE_SUCCESS != _result) {                  \
    fprintf(stderr,"An error occurred (%d), when calling %s at line: %d... \n\n", _result, _func, __LINE__ -1);   \
    exit(-1);                                         \
  }
#endif

#ifdef USE_PLASMA
#define HANDLE_ERROR(_result, _func)               \
  if(PLASMA_SUCCESS != _result) {                  \
    fprintf(stderr,"An error occurred (%d), when calling %s at line: %d... \n\n", _result, _func, __LINE__ -1);   \
  }  //break;                                         
  
#define HANDLE_ERROR_RET(_result, _func)               \
  if(PLASMA_SUCCESS != _result) {                  \
    fprintf(stderr,"An error occurred (%d), when calling %s at line: %d... \n\n", _result, _func, __LINE__ -1);   \
    exit(-1);                                         \
  }
#endif

#if defined(USE_TRACE)
#else 
#define START_TRACING() if( 0 ) {};
#define STOP_TRACING() if( 0 ) {};
#endif

#define START_TIMING(_t)               \
START_TRACING();                    \
_t = -cWtime();

#define STOP_TIMING(_t)                \
_t += cWtime();                      \
STOP_TRACING();

#endif //UTILS_H

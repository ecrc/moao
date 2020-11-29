/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef INTERSAMPLE_H
#define INTERSAMPLE_H

#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>


#ifdef __cplusplus
extern "C"{
#endif
struct isample_struct {
  long N;   // size of fft support
  long np;  // number of subaps
  long nidx;  // size of index map
  double lambda;  // waelength
  double dactupix;  // interactuator distance

  long *idx;
  double *abs2fi;
  double *otf_tel;

  double *Map;
  double *div;

  double *map;
  fftw_complex *tmp;
  double *dphi;

  fftw_plan plan_backward;
  fftw_plan plan_forward;

  char filename_idx[256], filename_abs2fi[256], filename_otftel[256];
};

void 
_eclat_float(float *ar, int nx, int ny);
void 
_eclat_double(double *ar, int nx, int ny);

void 
intersample_prepare(struct isample_struct *isample, long nidx, long nsubap, float Dtel, const char* path);
int 
intersample_init(struct isample_struct *isample);
void 
intersample_free(struct isample_struct *isample);

#ifdef __cplusplus
}
#endif

template<typename T>
int intersample_process(struct isample_struct *isample, T *data);

#endif //INTERSAMPLE_H

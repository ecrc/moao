#include<string.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<fftw3.h>

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
};

void 
_eclat_float(float *ar, int nx, int ny);
void 
_eclat_double(double *ar, int nx, int ny);

void
intersample(struct isample_struct *isample, double *data);
void 
init_intersample(struct isample_struct *isample, long np, long N, long nidx, 
		 double lambda, double dactupix, long *idx, double *abs2fi, 
		 double *otf_tel);
void 
free_intersample(struct isample_struct *isample);

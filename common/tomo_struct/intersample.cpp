/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include"intersample.hpp"
#include <fitsio.h>

#include <math.h>

#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

#ifdef DEBUG_MSG
#define FPRINTF(args...) fprintf(##args)
#else
#define FPRINTF(args...)
#endif

void _eclat_float(float *ar, int nx, int ny)
{
  int i,j,k1,k2;
  float a;

  for ( i=0 ; i<(nx/2) ; ++i ) {
    for ( j=0 ; j<(ny/2) ; ++j ) {
      k1 = i+j*nx;
      k2 = (i+nx/2)+(j+ny/2)*nx;
      a = ar[k1];
      ar[k1] = ar[k2];
      ar[k2] = a;
    }
  }
  for ( i=(nx/2) ; i<nx ; ++i ) {
    for ( j=0 ; j<(ny/2) ; ++j ) {
      k1 = i+j*nx;
      k2 = (i-nx/2)+(j+ny/2)*nx;
      a = ar[k1];
      ar[k1] = ar[k2];
      ar[k2] = a;
    }
  }
}

void _eclat_double(double *ar, int nx, int ny)
{
  int i,j,k1,k2;
  double a;

  for ( i=0 ; i<(nx/2) ; ++i ) {
    for ( j=0 ; j<(ny/2) ; ++j ) {
      k1 = i+j*nx;
      k2 = (i+nx/2)+(j+ny/2)*nx;
      a = ar[k1];
      ar[k1] = ar[k2];
      ar[k2] = a;
    }
  }
  for ( i=(nx/2) ; i<nx ; ++i ) {
    for ( j=0 ; j<(ny/2) ; ++j ) {
      k1 = i+j*nx;
      k2 = (i-nx/2)+(j+ny/2)*nx;
      a = ar[k1];
      ar[k1] = ar[k2];
      ar[k2] = a;
    }
  }
}

template<typename T>
int
intersample_process(struct isample_struct *isample, T *data) {
  int i,j;
  FPRINTF(stdout, "processing intersample...");fflush(stdout);

  // transformation d'un couple de decalages (dx,dy) en un indice du tableau 'Map'
  
  int nmap = (2*isample->np-1)*(2*isample->np-1);
  memset(isample->Map,0,nmap*sizeof(double));
  memset(isample->map,0,isample->N * isample->N * sizeof(double));
  for(i=0;i<isample->nidx;i++) {
    int indx = isample->idx[i]-1;
    if ((indx > nmap-1) || (indx < 0)){
        printf("idx not godd\n");
    }
    isample->Map[indx] += data[i];
  }

  
  for(i=0;i<nmap;i++) {
    if (isample->div[i] > 0)
      isample->Map[i] /= isample->div[i];
  }
  
  
  // size of the side of Cvvmap (always odd number)
  int ncmap = isample->N/2+1;
  int ncov  = isample->np*2-1;
  int nelem = (ncov-1)/2;

  int cpt=0;
  for (i=ncmap-nelem*isample->dactupix-1;i<ncmap+nelem*isample->dactupix;
       i+=isample->dactupix) {
    for (j=ncmap-nelem*isample->dactupix-1;j<ncmap+nelem*isample->dactupix;
	 j+=isample->dactupix) {
      isample->map[j+i*isample->N] = isample->Map[cpt];
      cpt++;
    }
  }
  
  fftw_execute(isample->plan_forward);
  
  for (i=0;i<isample->N;i++){
    for (j=0;j<isample->N/2+1; j++ )
    {
      (isample->tmp)[i*(isample->N/2+1)+j][0] *= isample->abs2fi[j*(isample->N)+i];
      (isample->tmp)[i*(isample->N/2+1)+j][1] = 0.0f;
    }
  }
  
  fftw_execute(isample->plan_backward);

  // From correlation to Dphi
  // Dphi(r) = 2*C(0) - 2*C(r)
  // We take advantage we need to do a multiplication to multiply by another factor
  // in the same line. This is to translate dphi from m^2 into rd^2
  // fact2 is here to normalize fft
  
  double fact = 2 * pow(2*3.14159/isample->lambda,2);
  double fact2 = isample->N*isample->N*isample->dactupix*isample->dactupix;
  //double fact2 = isample->dactupix*isample->dactupix;
  int idx_norm = isample->N*(ncmap-1)+ncmap-1;
  double norm = isample->dphi[idx_norm] / fact2 * fact;
  
//for(i=0;i<isample->N*isample->N;i++) {
//isample->dphi[i] = norm - isample->dphi[i] / fact2 * fact;
//}

  for(i=0;i<isample->N*isample->N;i++) {
    isample->dphi[i] = exp(-0.5*(norm - isample->dphi[i] / fact2 * fact));
  }
  
  _eclat_double(isample->dphi,isample->N,isample->N);

  /* compute PSF : re-use arrays */

  //memset(isample->tmp,0,isample->N*(isample->N/2+1)*sizeof(fftw_complex));
  //memset(isample->dphi,0,isample->N*isample->N*sizeof(double));

  for (i=0;i<isample->N;i++){
    for (j=0;j<isample->N/2+1; j++ )
    {
      (isample->tmp)[i*(isample->N/2+1)+j][0] = 
	isample->otf_tel[j*(isample->N)+i]*isample->dphi[j*(isample->N)+i];
      (isample->tmp)[i*(isample->N/2+1)+j][1] = 0.0f;
    }
  }
  
  fftw_execute(isample->plan_backward);
  /*
  for(i=0;i<isample->N*isample->N;i++) {
    isample->dphi[i] = isample->otf_tel[i];
  }
  */
  _eclat_double(isample->dphi,isample->N,isample->N);

  FPRINTF(stdout, "done...\n");fflush(stdout);
  return 1;
}
template int intersample_process(struct isample_struct *isample, float *data);
template int intersample_process(struct isample_struct *isample, double *data);

void 
intersample_prepare(struct isample_struct *isample, long nidx, long np, float Dtel, const char* path){
  /*
   > fourier.N
   512
   > fourier.ud
   0.15873
   > fourier.dactupix
   3
   > numberof(ind)
   28472896
   > dm.Ndiam
   85
   > fourier.lambdaIR
   1.65e-06
   * */
  isample->np       = np; //nb subap along the diamter +1 (Nssp[truth Sensot]+1 (nactu)
  isample->N        = MAX(512,8*pow(2,(long)(log(np)+0.5)));//512;
  isample->nidx     = nidx;
  isample->lambda   = 1.65e-06;
  //dactu=Dtel/(np-1);
  //champ_dphi=3*Dtel
  //ud=champ_dphi/(float)N
  //dactupix=dactu/ud+1
  isample->dactupix=isample->N/(3*(np-1))+1;

  strcpy(isample->filename_idx, path);
  strcat(isample->filename_idx, "/idx.fits");
  strcpy(isample->filename_abs2fi, path);
  strcat(isample->filename_abs2fi, "/abs2fi.fits");
  strcpy(isample->filename_otftel, path);
  strcat(isample->filename_otftel, "/otftel.fits");
}

int
fits_load_file(long rows, long cols, double* buf, char* file_name){
  
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int status,  nfound, anynull;
  long naxes[2], fpixel;//, npixels;

  status = 0;
  fpixel   = 1;
  long nullval  = 0;                /* don't check for null values in the image */

  if ( fits_open_file(&fptr, file_name, READONLY, &status) )
    fits_report_error(stderr, status);
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
    fits_report_error(stderr, status);
  if( naxes[0] != rows || naxes[1] != cols ) return 0;
  if ( fits_read_img(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1], &nullval,
                     buf, &anynull, &status) )
    fits_report_error(stderr, status);
  if ( fits_close_file(fptr, &status) )
    fits_report_error(stderr, status);
  
  return 1;
}

int 
intersample_init(struct isample_struct *isample){
  
  FPRINTF(stdout, "initializing intersample...");fflush(stdout);
  
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int status,  nfound, anynull;
  long naxes[2], fpixel, npixels;

  status = 0;
  fpixel   = 1;
  long nullval  = 0;                /* don't check for null values in the image */

  isample->idx     = (long*)malloc(isample->nidx * sizeof(long));
  if(!isample->idx) return 0;
  //printf("reading idx ... \n");
  if ( fits_open_file(&fptr, isample->filename_idx, READONLY, &status) )
    fits_report_error(stderr, status);
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
    fits_report_error(stderr, status);
  if( naxes[0] != isample->nidx) return 0;
  if ( fits_read_img(fptr, TLONG, fpixel, isample->nidx, &nullval,
                     isample->idx, &anynull, &status) )
    fits_report_error(stderr, status);
  if ( fits_close_file(fptr, &status) )
    fits_report_error(stderr, status);

  isample->abs2fi  = (double*)malloc(isample->N * isample->N * sizeof(double));
  if(!isample->abs2fi) return 0;
  //printf("reading abs2fi ... \n");
  if ( fits_open_file(&fptr, isample->filename_abs2fi, READONLY, &status) )
    fits_report_error(stderr, status);
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
    fits_report_error(stderr, status);
  if( naxes[0] != isample->N) return 0;
  npixels  = naxes[0] * naxes[1];    
  if ( fits_read_img(fptr, TDOUBLE, fpixel, npixels, &nullval,
                     isample->abs2fi, &anynull, &status) )
    fits_report_error(stderr, status);
  if ( fits_close_file(fptr, &status) )
    fits_report_error(stderr, status);
  
  isample->otf_tel = (double*)malloc(isample->N * isample->N * sizeof(double));
  if(!isample->otf_tel) return 0;
  //printf("reading otftel ... \n");
  if ( fits_open_file(&fptr, isample->filename_otftel, READONLY, &status) )
    fits_report_error(stderr, status);
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
    fits_report_error(stderr, status);
  if( naxes[0] != isample->N) return 0;
  npixels  = naxes[0] * naxes[1]; 
  if ( fits_read_img(fptr, TDOUBLE, fpixel, npixels, &nullval,
                      isample->otf_tel, &anynull, &status) )
    fits_report_error(stderr, status);
  if ( fits_close_file(fptr, &status) )
    fits_report_error(stderr, status);

  isample->Map     = (double*)malloc((2 * isample->np - 1) * (2 * isample->np - 1) * sizeof(double));
  isample->div     = (double*)malloc((2 * isample->np - 1) * (2 * isample->np - 1) * sizeof(double));

  memset(isample->Map, 0, (2 * isample->np - 1) * (2 * isample->np - 1) * sizeof(double));
  memset(isample->div, 0, (2 * isample->np - 1) * (2 * isample->np - 1) * sizeof(double));
  
  //printf("line %d... \n", __LINE__);
  long i;
  for(i = 0; i < isample->nidx; i++) {
    //printf("line %d... isample->idx[%d]=%d\n", __LINE__, i, isample->idx[i] - 1);
    isample->div[ isample->idx[i] - 1] += 1;
  }
  
  //printf("line %d... \n", __LINE__);
  isample->map  = (double*)fftw_malloc(isample->N * isample->N * sizeof(double));
  isample->tmp  = (fftw_complex*)fftw_malloc(isample->N * (isample->N / 2 + 1) * sizeof(fftw_complex));
  isample->dphi = (double*)fftw_malloc(isample->N * isample->N * sizeof(double));
  
  //printf("line %d... \n", __LINE__);
  /* init fftw plans */
  isample->plan_forward = fftw_plan_dft_r2c_2d (isample->N, isample->N,  isample->map,
                                                isample->tmp, FFTW_ESTIMATE );
  //printf("line %d... \n", __LINE__);
  isample->plan_backward = fftw_plan_dft_c2r_2d (isample->N, isample->N, isample->tmp, 
						  isample->dphi, FFTW_ESTIMATE );

  FPRINTF(stdout, "done...\n");fflush(stdout);
  return 1;
}

void 
intersample_free(struct isample_struct *isample){
  free(isample->idx);
  free(isample->abs2fi);
  free(isample->otf_tel);
  fftw_free(isample->map);
  fftw_free(isample->tmp);
  fftw_free(isample->dphi);
  fftw_destroy_plan(isample->plan_forward);
  fftw_destroy_plan(isample->plan_backward);
  free(isample->Map);
  free(isample->div);
}

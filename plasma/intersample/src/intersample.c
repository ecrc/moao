#include"intersample.h"

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

void
intersample(struct isample_struct *isample, double *data) {
  int i,j;

  // transformation d'un couple de decalages (dx,dy) en un indice du tableau 'Map'
  
  int nmap = (2*isample->np-1)*(2*isample->np-1);
  for(i=0;i<isample->nidx;i++) {
    int indx = isample->idx[i]-1;
    if ((indx > nmap-1) || (indx < 0)) printf("idx not godd\n");
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
      isample->map[i+j*isample->N] = isample->Map[cpt];
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
  int idx_norm = isample->N*ncmap+ncmap-1;
  double norm = isample->dphi[idx_norm] / fact2 * fact;
  
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
}

void 
init_intersample(struct isample_struct *isample, long np, long N, long nidx, 
		 double lambda, double dactupix, long *idx, double *abs2fi, 
		 double *otf_tel){

  isample->np       = np;
  isample->N        = N;
  isample->nidx     = nidx;
  isample->lambda   = lambda;
  isample->dactupix = dactupix;

  isample->idx     = (long*)malloc(nidx*sizeof(long));
  isample->abs2fi  = (double*)malloc(N*N*sizeof(double));
  isample->otf_tel = (double*)malloc(N*N*sizeof(double));

  memcpy(isample->idx, idx, nidx*sizeof(long));
  memcpy(isample->abs2fi, abs2fi, N*N*sizeof(double));
  memcpy(isample->otf_tel, otf_tel, N*N*sizeof(double));

  isample->Map     = (double*)malloc((2*np-1)*(2*np-1)*sizeof(double));
  isample->div     = (double*)malloc((2*np-1)*(2*np-1)*sizeof(double));

  memset(isample->Map,0,(2*np-1)*(2*np-1)*sizeof(double));
  memset(isample->div,0,(2*np-1)*(2*np-1)*sizeof(double));
  
  int i;
  for(i=0;i<isample->nidx;i++) {
    isample->div[isample->idx[i]-1] += 1;
  }

  isample->map  = (double*)fftw_malloc(N*N*sizeof(double));
  isample->tmp  = (fftw_complex*)fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
  isample->dphi = (double*)fftw_malloc(N*N*sizeof(double));

  /* init fftw plans */
  isample->plan_forward = fftw_plan_dft_r2c_2d (N, N,  isample->map,
						 isample->tmp, FFTW_ESTIMATE );
  isample->plan_backward = fftw_plan_dft_c2r_2d (N, N, isample->tmp, 
						  isample->dphi, FFTW_ESTIMATE );

}

void 
free_intersample(struct isample_struct *isample){
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

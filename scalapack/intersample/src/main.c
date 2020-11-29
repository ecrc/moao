#include"intersample.h"
#include <sys/time.h>
#include <fitsio.h>

double timeval_subtract(struct timeval *t2, struct timeval *t1)
{
  return ((t2->tv_usec + 1000000 * t2->tv_sec) - (t1->tv_usec + 1000000 * t1->tv_sec)) * 1e-6;
}

int
main(int argc, char *argv[]) {
  /*
> info,cvv
 array(double,204,204)
> fourier.N
512
> fourier.ud
0.233333
> fourier.dactupix
12
> numberof(ind)
41616
> dm.Ndiam;
16
> fourier.lambdaIR
1.65e-06
   */
  struct timeval start, end;

  struct isample_struct isample;

  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int status,  nfound, anynull;
  long naxes[2], fpixel, npixels;

  status = 0;
  fpixel   = 1;
  long nullval  = 0;                /* don't check for null values in the image */

  long np = 16;

  printf("reading cvv ... \n");
  if ( fits_open_file(&fptr, "cvv.fits", READONLY, &status) )
    fits_report_error(stderr, status);
  /* read the NAXIS1 and NAXIS2 keyword to get image size */
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
    fits_report_error(stderr, status);
  npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
  double *data = (double*)malloc(npixels*sizeof(double));
  if ( fits_read_img(fptr, TDOUBLE, fpixel, npixels, &nullval,
		     data, &anynull, &status) )
    fits_report_error(stderr, status);
  if ( fits_close_file(fptr, &status) )
    fits_report_error(stderr, status);

  printf("reading idx ... \n");
  if ( fits_open_file(&fptr, "idx.fits", READONLY, &status) )
    fits_report_error(stderr, status);
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
    fits_report_error(stderr, status);
  long nidx  = naxes[0];      
  long *idx = (long*)malloc(npixels*sizeof(long));
  if ( fits_read_img(fptr, TLONG, fpixel, npixels, &nullval,
		     idx, &anynull, &status) )
    fits_report_error(stderr, status);
  if ( fits_close_file(fptr, &status) )
    fits_report_error(stderr, status);

  printf("reading abs2fi ... \n");
  if ( fits_open_file(&fptr, "abs2fi.fits", READONLY, &status) )
    fits_report_error(stderr, status);
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
    fits_report_error(stderr, status);
  npixels  = naxes[0] * naxes[1];    
  double *abs2fi = (double*)malloc(npixels*sizeof(double));
  if ( fits_read_img(fptr, TDOUBLE, fpixel, npixels, &nullval,
		     abs2fi, &anynull, &status) )
    fits_report_error(stderr, status);
  if ( fits_close_file(fptr, &status) )
    fits_report_error(stderr, status);

  printf("reading otftel ... \n");
  if ( fits_open_file(&fptr, "otftel.fits", READONLY, &status) )
    fits_report_error(stderr, status);
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
    fits_report_error(stderr, status);
  npixels  = naxes[0] * naxes[1]; 
  double *otftel = (double*)malloc(npixels*sizeof(double));
  if ( fits_read_img(fptr, TDOUBLE, fpixel, npixels, &nullval,
		     otftel, &anynull, &status) )
    fits_report_error(stderr, status);
  if ( fits_close_file(fptr, &status) )
    fits_report_error(stderr, status);

  printf("reading psf ... \n");
  if ( fits_open_file(&fptr, "psf.fits", READONLY, &status) )
    fits_report_error(stderr, status);
  if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) )
    fits_report_error(stderr, status);
  npixels  = naxes[0] * naxes[1];       
  double *refdata = (double*)malloc(npixels*sizeof(double));
  if ( fits_read_img(fptr, TDOUBLE, fpixel, npixels, &nullval,
		     refdata, &anynull, &status) )
    fits_report_error(stderr, status);
  if ( fits_close_file(fptr, &status) )
    fits_report_error(stderr, status);

  printf("Init : %ld  ...", np);
  gettimeofday(&start, NULL);
  init_intersample(&isample,np,512,nidx,0.00000165,12,idx,abs2fi,otftel);
  gettimeofday(&end, NULL);
  printf(" in %0.3f seconds\n", timeval_subtract(&end, &start));

  printf("compute intersample ... ");
  gettimeofday(&start, NULL);
  intersample(&isample, data);
  gettimeofday(&end, NULL);
  printf(" in %0.3f seconds\n", timeval_subtract(&end, &start));

  printf("Error = ");
  double err = 0.0;
  int i;
  for (i=0;i<npixels;i++) err += fabs(refdata[i]>1.e6) ?  (refdata[i]-isample.dphi[i]) / refdata[i] : 0;
  //err /= npixels;
  printf("%0.6f\n",err);
  /*
  if (fits_create_file(&fptr, "result.fits", &status))
    fits_report_error(stderr, status);
  naxes[0] = naxes[1] = 512;
  if (fits_create_img(fptr, DOUBLE_IMG, 2, naxes, &status))
    fits_report_error(stderr, status);
  long fpix[] = { 1, 1 };
  if (fits_write_pix(fptr, TDOUBLE, fpix, 512*512, isample.dphi, &status))
    fits_report_error(stderr, status);
  if (fits_close_file(fptr, &status))
    fits_report_error(stderr, status);
  */
  printf("free struct ... \n");
  free_intersample(&isample);

  printf("free temp data ... \n");
  free(data);
  free(idx);
  free(abs2fi);
  free(otftel);

  return EXIT_SUCCESS;
}

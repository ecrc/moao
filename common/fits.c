/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "fits.h"
#include <fitsio.h>
#include <string.h>


void concatenateFileName(char *fileName, char *filePath, char *typeMat, char *suffix){
    strcpy(fileName,filePath);
    strcat(fileName,typeMat);
    strcat(fileName,suffix);
    strcat(fileName,".fits");
}

int fitsExists(char *fileName){
    fitsfile *fptr;
    int status;
    if(fits_open_file(&fptr,fileName,READONLY,&status)){
        fits_report_error(stderr,status);
        return 0;
    }
    return 1;
}

void getFitsDims(char *fileName,long *naxes){
    fitsfile *fptr;
    int status,nfound;

    status = 0;

    if(fits_open_file(&fptr,fileName,READONLY,&status)){
        naxes[0]=-1;
        fits_report_error(stderr,status);
        return;
    }
    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ){
        naxes[0]=-2;
        fits_report_error(stderr, status);
        return;
    }
}

void readFits(char *fileName,double *data){
    fitsfile *fptr;
    int status,nfound,nullval, anynull;
    long naxes[2], fpixel;

    status = 0;
    fpixel   = 1;
    nullval  = 0;                /* don't check for null values in the image */

    if(fits_open_file(&fptr,fileName,READONLY,&status)){
        fits_report_error(stderr,status);
        return;
    }
    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ){
        fits_report_error(stderr, status);
        return;
    }
    
    if ( fits_read_img(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1], &nullval,
                     data, &anynull, &status) )
    fits_report_error(stderr, status);


    fits_close_file(fptr,&status);
    fits_report_error(stderr,status);
}

void writeFits(char *file_name, int nx, int ny, real_t* data)
{
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int status;
  int bitpix    =  REAL_IMG;
  long naxis    =   2;  /* 2-dimensional image                            */
  long naxes[2] = { nx , ny};
  long fpixel, nelements;
  
  
  remove(file_name);               /* Delete old file if it already exists */
  status = 0;         /* initialize status before calling fitsio routines */
  
  if ( fits_create_file(&fptr, file_name, &status) )
    fits_report_error(stderr, status);
  
  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    fits_report_error(stderr, status);
  
  fpixel = 1;                               /* first pixel to write      */
  nelements = naxes[0] * naxes[1];          /* number of pixels to write */
  
  if ( fits_write_img(fptr, TREAL, fpixel, nelements, data, &status) )
    fits_report_error(stderr, status);
  
  if ( fits_close_file(fptr, &status) )
    fits_report_error(stderr, status);
}


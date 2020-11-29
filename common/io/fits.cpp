/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "fits.hpp"
#include <string.h>
#include <fitsio.h>
#include <type_traits>
#include <algorithm>

template<typename> struct fitsTypes;

template<>
struct fitsTypes <float>
{
    static const int TREAL =TFLOAT ;
    static const int REAL_IMG = FLOAT_IMG;
};
template<>
struct fitsTypes<double>
{
    static const int TREAL =TDOUBLE ;
    static const int REAL_IMG = DOUBLE_IMG;
};


//void concatenateFileName(char *fileName, char *filePath, char *typeMat, char *suffix){
//    strcpy(fileName,filePath);
//    strcat(fileName,typeMat);
//    strcat(fileName,suffix);
//    strcat(fileName,".fits");
//}

int fitsExists(const char *fileName){
    fitsfile *fptr;
    int status;
    if(fits_open_file(&fptr,fileName,READONLY,&status)){
        fits_report_error(stderr,status);
        return 0;
    }
    return 1;
}

void getFitsDims(const char *fileName,long *naxes){
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

template<typename T>
void readFits(const char *fileName,T *data){
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
    
    double * data_double;
    if(std::is_same<T,float>::value){
        data_double= (double *) calloc(naxes[0]*naxes[1],sizeof(double));
    }
    else{
        data_double=(double*)data;
    }

    if ( fits_read_img(fptr, TDOUBLE, fpixel, naxes[0]*naxes[1], &nullval,
                     data_double, &anynull, &status) )
    fits_report_error(stderr, status);

    if(std::is_same<T,float>::value){
        std::copy(data_double,data_double+(naxes[0]*naxes[1]),data);
        free(data_double);
    }

    fits_close_file(fptr,&status);
    fits_report_error(stderr,status);
}
template void readFits(const char *fileName,float *data);
template void readFits(const char *fileName,double *data);

template<typename T>
void writeFits(const char *file_name, int nx, int ny, T* data)
{
  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
  int status;
  int bitpix=DOUBLE_IMG;
  long naxis    =   2;  /* 2-dimensional image                            */
  long naxes[2] = { nx , ny};
  long fpixel, nelements;
  
  double *data_double;
  if(std::is_same<T,float>::value){
    data_double = (double*)calloc(nx*ny,sizeof(double));
    std::copy(data,data+nx*ny,data_double);
  }
  else{
    data_double=(double*)data;
  }
  
  remove(file_name);               /* Delete old file if it already exists */
  status = 0;         /* initialize status before calling fitsio routines */
  
  if ( fits_create_file(&fptr, file_name, &status) )
    fits_report_error(stderr, status);
  
  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
    fits_report_error(stderr, status);
  
  fpixel = 1;                               /* first pixel to write      */
  nelements = naxes[0] * naxes[1];          /* number of pixels to write */
  
  if ( fits_write_img(fptr, TDOUBLE, fpixel, nelements, data_double, &status) )
    fits_report_error(stderr, status);

  if(std::is_same<T,float>::value){
      free(data_double);
  }
  
  if ( fits_close_file(fptr, &status) )
    fits_report_error(stderr, status);
}
template void writeFits(const char *file_name, int nx, int ny, float* data);
template void writeFits(const char *file_name, int nx, int ny, double* data);

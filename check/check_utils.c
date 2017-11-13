/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "check_utils.h"

#include <stdlib.h>
#include <unistd.h>
#include <math.h>
//#include <lapacke.h>
//#include <fitsio.h>

void* myMalloc(size_t size){
    void *ptr=malloc(size);
#ifdef SAFEMALLOC 
    if(!ptr){
        fprintf(stderr, "Out of Memory\n");
        exit(1);
    }
#endif
    return ptr;
}

void* myCalloc(size_t num,size_t size){

    void *ptr=calloc(num, size);
#ifdef SAFEMALLOC
    if(!ptr){
        fprintf(stderr, "Out of Memory\n");
        exit(1);
    }   
#endif
    return ptr;
}


double compareMatrices(double *A, double *B, long nrows, long ncols, char uplo, char trans, double *err, double *vA, double *vB){
    //trans==T means B is to be transpose and B shape is (ncols,nrows)
    double diff=0, tmp;
    long i,j,/*i0,*/nA,nB;
    long start=0;
    long end=ncols;
    double ma=A[0],mb=B[0];

    double norm=0;

    *err=0;
    *vA=0;
    *vB=0;

    //i0=0;
    for(i=0;i<nrows;i++){
        if(uplo=='U')
            start=i;
        else if(uplo=='L')
            end=i;
        
        for(j=start;j<end;j++){
            //nA=ncols*i+j;
            nA=i+j*nrows;
            if(trans=='T')
                nB=nrows*j+i;
            else
                nB=ncols*i+j;
            tmp=fabs(A[nA]-B[nB]);
            norm+=tmp*tmp;
            if(tmp>diff){
                diff=tmp;
                //i0=nA;
                ma=A[nA];
                mb=B[nB];
            }
        }
    }
    //fprintf(stderr," max err at %ld: %e - %e\n%e\n",i0,ma,mb,diff);
    *err=diff;
    *vA=ma;
    *vB=mb;
    return norm;
}

// second parameter forced to real_t to fit with yorick data format
double compareMatrices2(real_t *A, double *B, long nrows, long ncols, char uplo, char trans, real_t *err, real_t *vA, real_t *vB){
    //trans==T means B is to be transpose and B shape is (ncols,nrows)
    real_t diff=0, tmp;
    long i,j,/*i0,*/nA,nB;
    long start=0;
    long end=ncols;
    real_t ma=A[0],mb=B[0];

    double norm=0;

    *err=0;
    *vA=0;
    *vB=0;

    //i0=0;
    for(i=0;i<nrows;i++){
        if(uplo=='U')
            start=i;
        else if(uplo=='L')
            end=i;
        
        for(j=start;j<end;j++){
            //nA=ncols*i+j;
            nA=i+j*nrows;
            if(trans=='T')
                nB=nrows*j+i;
            else
                nB=ncols*i+j;
            tmp=fabs(A[nA]-B[nB]);
            norm+=tmp*tmp;
            if(tmp>diff){
                diff=tmp;
                //i0=nA;
                ma=A[nA];
                mb=B[nB];
            }
        }
    }
    //fprintf(stderr," max err at %ld: %e - %e\n%e\n",i0,ma,mb,diff);
    *err=diff;
    *vA=ma;
    *vB=mb;
    return norm;

}

void minMax(real_t *A,long nbElem,real_t *min,real_t*max){
    *min=A[0];
    *max=A[0];
    int i;
    for(i=1;i<nbElem;i++){
        if(A[i]<*min){
            *min=A[i];
        }
        if(A[i]>*max){
            *max=A[i];
        }
    }   
}

void copy(double *src, real_t *dst,long nrows, long ncols, char uplo, char trans){
    long i,j;
    long nSrc,nDst;
    long start=0;
    long end=ncols;

    for(i=0;i<nrows;i++){

        if(uplo=='U')
            start=i;
        else if(uplo=='L')
            end=i;

        for(j=start;j<end;j++){
            //nDst=ncols*i+j;
            nDst=i+j*nrows;
            if(trans=='T')
                nSrc=nrows*j+i;
            else
                nSrc=ncols*i+j;
            dst[nDst]=src[nSrc];
        }
    }
}

void saveMat(char *fileName, real_t *mat, int nrows, int ncols){
    FILE *file;
    int i,j;

    file=fopen(fileName,"w+");
    for(i=0;i<nrows;i++){
        for(j=0;j<ncols;j++){
            fprintf(file,"%f ",mat[i*ncols+j]);
        }
        fprintf(file,"\n");
    }
    fclose(file);
}

#ifdef USE_PLASMA
#if 0
never used, need to provide proper compareMatrices function
real_t p_compareMatrices(real_t *A, PLASMA_desc descB, /*real_t *B,*/ char uplo, char trans){
    real_t *tmp=(real_t*)malloc(descB.lm*descB.ln*sizeof(real_t));
    p_tile_to_lapack(descB, tmp,0);
    real_t maxerr,errv1,errv2;
    compareMatrices(tmp,A,descB.lm,descB.ln,uplo,trans,&maxerr,&errv1,&errv2);
    free(tmp);
    return maxerr;
}
#endif


void p_get_TileInfo(PLASMA_desc descA, struct TileInfo *ti, int m, int n){
/** from descriptor.h (plasma_desc)
 *  Tile matrix descriptor
 *
 *  Matrices are stored in a contiguous data chunk containing in order
 *  A11, A21, A12, A22 with :
 *
 *           n1      n2
 *      +----------+---+
 *      |          |   |    With m1 = lm - (lm%mb)
 *      |          |   |         m2 = lm%mb
 *  m1  |    A11   |A12|         n1 = ln - (ln%nb)
 *      |          |   |         n2 = ln%nb
 *      |          |   |
 *      +----------+---+
 *  m2  |    A21   |A22|
 *      +----------+---+
 *
 */

    ti->m=m;
    ti->n=n;
    ti->addr=plasma_getaddr(descA,m,n);
    
    if(m<descA.lm1){
        //tile belongs to A11 or A12
        ti->nrows=descA.mb;
        ti->i    =descA.mb*m;
    }
    else{
        //tile belongs to A21 or A22
        ti->nrows=descA.lm%descA.mb;
        ti->i    =descA.mb*descA.lm1;
    }
    if(n<descA.ln1){
        //tile belongs to A11 or A21
        ti->ncols=descA.nb;
        ti->j    =descA.nb*n;
    }
    else{
        //tile belongs to A12 or A22
        ti->ncols=descA.ln%descA.nb;
        ti->j    =descA.nb*descA.ln1;
    }
}

void p_tile_to_lapack(PLASMA_desc descA, real_t *B, int rowMajor){

    int m,n,i,j;
    struct TileInfo ti;
    long start;
    real_t *A;
    for(m=0;m<descA.mt;m++){
        for(n=0;n<descA.nt;n++){
        //for each tile of descA
            //get info on tile
            p_get_TileInfo(descA,&ti,m,n);
            //copy tile
            A=(real_t*)ti.addr;
            for(i=0;i<ti.nrows;i++){
                for(j=0;j<ti.ncols;j++){
                    start=descA.ln*ti.i+ti.j;
                    if(rowMajor){
                        B[start+descA.ln*i+j]=A[ti.ncols*i+j];
                    }
                    else{
                        B[start+descA.ln*i+j]=A[ti.nrows*j+i];
                    }
                }
            }
        }
    }

}

void p_lapack_to_tile(PLASMA_desc descA, real_t *B, int rowMajor){

    int m,n,i,j;
    struct TileInfo ti;
    long start;
    real_t *A;
    for(m=0;m<descA.mt;m++){
        for(n=0;n<descA.nt;n++){
        //for each tile of descA
            //get info on tile
            p_get_TileInfo(descA,&ti,m,n);
            //copy tile
            A=(real_t*)ti.addr;
            for(i=0;i<ti.nrows;i++){
                for(j=0;j<ti.ncols;j++){
                    start=descA.ln*ti.i+ti.j;
                    if(rowMajor){
                        A[ti.ncols*i+j]=B[start+descA.ln*i+j];
                    }
                    else{
                        A[ti.nrows*j+i]=B[start+descA.ln*i+j];
                    }
                }
            }
        }
    }
}
#endif
//void concatenateFileName(char *fileName, char *filePath, char *typeMat, char *suffix){
//    strcpy(fileName,filePath);
//    strcat(fileName,typeMat);
//    strcat(fileName,suffix);
//    strcat(fileName,".fits");
//}
//
//void getFitsDims(char *fileName,long *naxes){
//    fitsfile *fptr;
//    int status,nfound,nullval, anynull;
//    long fpixel, npixels;
//
//    status = 0;
//    fpixel   = 1;
//    nullval  = 0;                /* don't check for null values in the image */
//
//    if(fits_open_file(&fptr,fileName,READONLY,&status)){
//        naxes[0]=-1;
//        fits_report_error(stderr,status);
//        return;
//    }
//    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ){
//        naxes[0]=-2;
//        fits_report_error(stderr, status);
//        return;
//    }
//}
//
//void readFits(char *fileName,real_t *data){
//    fitsfile *fptr;
//    int status,nfound,nullval, anynull;
//    long naxes[2], fpixel, npixels;
//
//    status = 0;
//    fpixel   = 1;
//    nullval  = 0;                /* don't check for null values in the image */
//
//    if(fits_open_file(&fptr,fileName,READONLY,&status)){
//        fits_report_error(stderr,status);
//        return;
//    }
//    if ( fits_read_keys_lng(fptr, "NAXIS", 1, 2, naxes, &nfound, &status) ){
//        fits_report_error(stderr, status);
//        return;
//    }
//    
//    if ( fits_read_img(fptr, Treal_t, fpixel, naxes[0]*naxes[1], &nullval,
//                     data, &anynull, &status) )
//    fits_report_error(stderr, status);
//
//
//    fits_close_file(fptr,&status);
//    fits_report_error(stderr,status);
//}
//
//void writeFits(char *file_name, int nx, int ny, real_t* data)
//{
//  fitsfile *fptr;       /* pointer to the FITS file, defined in fitsio.h */
//  int status;
//  int bitpix   =  real_t_IMG;
//  long naxis    =   2;  /* 2-dimensional image                            */
//  long naxes[2] = { nx , ny};
//  long fpixel, nelements;
//  
//  
//  remove(file_name);               /* Delete old file if it already exists */
//  status = 0;         /* initialize status before calling fitsio routines */
//  
//  if ( fits_create_file(&fptr, file_name, &status) )
//    fits_report_error(stderr, status);
//  
//  if ( fits_create_img(fptr,  bitpix, naxis, naxes, &status) )
//    fits_report_error(stderr, status);
//  
//  fpixel = 1;                               /* first pixel to write      */
//  nelements = naxes[0] * naxes[1];          /* number of pixels to write */
//  
//  if ( fits_write_img(fptr, Treal_t, fpixel, nelements, data, &status) )
//    fits_report_error(stderr, status);
//  
//  if ( fits_close_file(fptr, &status) )
//    fits_report_error(stderr, status);
//}

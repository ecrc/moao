#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include "moao_scalapack.h"
#include "flops.h"

void Xgeadd(char trans, int M, int  N, real_t alpha, real_t *A, int ldA, real_t beta, real_t *C, int ldC);
/**
* Compute the tomographic reconstructor
*
* input:
*   int  nmeas   : number of measurements
*   int  nmeasts : number of measurements of the truth sensor
*   real_t *Cmm  : covariance matrix between the sensors
*                  Cmm(mloc_nmeas, nloc_nmeas)
*   int *descCmm : descriptor for the covariance matrix Cmm
*   real_t *Ctm  : covariance matrix between the truth sensors and the other sensors
*                  Ctm(mloc_nmeasts,nloc_nmeas)
*   int *descCtm : descriptor for the covariance matrix Ctm
*
* output:
*   real_t *Ctm  : tomographic reconstructor
*
*/
int reconstructor(int  nmeas, int  nmeasts, real_t *Cmm, int *descCmm, real_t *Ctm, int *descCtm, long maxgal, long offset_gal) {
    int info;
    int i1 = 1, nbgal;
    real_t alpha = 1.0;

    pXpotrf_( "L", &nmeas, Cmm, &i1, &i1, descCmm, &info );
    HANDLE_ERROR(info, "pdpotrf");

    for(nbgal=0;nbgal<maxgal;nbgal++) {
        pXtrsm_( "R", "L", "T", "N", &nmeasts, &nmeas, &alpha, 
                 Cmm, &i1, &i1, descCmm,
                 Ctm+nbgal*offset_gal, &i1, &i1, descCtm );

        pXtrsm_( "R", "L", "N", "N", &nmeasts, &nmeas, &alpha, 
                 Cmm, &i1, &i1, descCmm,
                 Ctm+nbgal*offset_gal, &i1, &i1, descCtm );
    }

    return info;
}



/**
* Compute the matrix Cee Cvv
*
* input:
*   int  nmeas   : number of measurements
*   int  nmeasts : number of measurements of the truth sensor
*   int  nact    : number of actuators
*   real_t *Cmm  : covariance matrix between the sensors
*                  Cmm(nmeas, nmeas)
*   int *descCmm : descriptor for the covariance matrix Cmm
*   real_t *Ctt  : covariance matrix of the truth sensor
*                  Ctt(nmeasts,nmeasts)
*   int *descCtt : descriptor for the covariance matrix Ctt
*   real_t *Ctm  : covariance matrix between the truth sensors and the other sensors
*                  Ctm(nmeasts,nmeas)
*   int *descCtm : descriptor for the covariance matrix Ctm
*   real_t *R    : tomographic reconstructor
*                  R(nmeasts,nmeas))
*   int *descR   : descriptor for the tomographic reconstructor R
*   real_t *Dx   : interaction matrix
*                  Dx(nact,nmeasts)
*   int *descDx  : descriptor for the interaction matrix Dx
*   real_t *Cee  : covariance matrix
*                  Cee(nmeasts,nmeasts)
*   int *descCee : descriptor for the covariance matrix Cee
*   real_t *Cvv  : covariance matrix
*                  Cvv(nact,nact)
*   int *descCvv : descriptor for the covariance matrix Cvv
*   real_t *Tmp  : temporary buffer
*                  Tmp(nact,nmeasts)
*   int *descTmp : descriptor for the temporary buffer Tmp
*
* output:
*   real_t *Ctm : Ctm is overwritten!
*   real_t *Ctt : Ctt is overwritten!
*   real_t *Cee : tomographic error
*   real_t *Cvv : covariance matrix
*
*/
int compute_Cee_Cvv(int  nmeas, int  nmeasts, int  nact, real_t *Cmm, int *descCmm, real_t *Ctt, int *descCtt, real_t *Ctm, int *descCtm, real_t *R, int *descR, real_t *Dx, int *descDx, real_t *Cee, int *descCee, real_t *Cvv, int *descCvv, real_t *Tmp, int *descTmp, double* flops){
    real_t alpha, beta;
    int i1 = 1;
    int info, ioffd;


    alpha = -1.;
    beta = 1.;
    pXsyr2k_( "L", "N", &nmeasts, &nmeas, &alpha,
              Ctm, &i1, &i1, descCtm,
              R, &i1, &i1, descR,
              &beta,
              Ctt, &i1, &i1, descCtt );
    *flops += FLOPS_XSYR2K(nmeasts, nmeas);

    alpha = 1.;
    beta = 0.;
    pXsymm_( "R", "L", &nmeasts, &nmeas, &alpha,
             Cmm, &i1, &i1, descCmm,
             R, &i1, &i1, descR,
             &beta,
             Ctm, &i1, &i1, descCtm );
    *flops += FLOPS_XSYMM('R', nmeasts, nmeas);

    alpha = 0.;
    beta = 0.;
    pXlaset_( "U", &nmeasts, &nmeasts,
              &alpha, &beta,
              Cee, &i1, &i1, descCee );

    pXlacpy_( "L", &nmeasts, &nmeasts,
              Ctt, &i1, &i1, descCtt,
              Cee, &i1, &i1, descCee );


    //printf("descCtm %d\n", descCtm[0]);
    //printf("descR %d\n", descR[0]);
    //printf("descCee %d\n", descCee[0]);
    //printf("descCtt %d\n", descCtt[0]);
    //printf("descCtm %d\n", descCtm[0]);
    //printf("descCvv %d\n", descCvv[0]);
    //printf("descDx %d\n", descDx[0]);
    //printf("descTmp %d\n", descTmp[0]);
    //printf("nmeas %d, nmeasts %d nact %d\n", nmeas, nmeasts, nact);

    // This step could further improved by scaling directly the diagonal instead of scaling the whole upper part.
/*
    alpha = 1.0;
    beta = 1.0;

    pdgeadd_( "T", &nmeasts, &nmeasts, 
              &alpha, 
              Cee, &i1, &i1, descCee,
              &beta, 
              Cee, &i1, &i1, descCee );

    alpha = 0.5;
    cblas_dscal(nmeasts, alpha, Cee, nmeasts+1);

*/
    alpha = 2.0;
    beta = 1.0;
    pXlascl_( "U", &alpha, &beta,
              &nmeasts, &nmeasts,
              Cee, &i1, &i1, descCee, &info );
    HANDLE_ERROR(info, "pdlascl");
    *flops += FLOPS_XLASCL(nmeasts, nmeasts);

    pXlacpy_( "A", &nmeasts, &nmeasts,
              Cee, &i1, &i1, descCee,
              Ctt, &i1, &i1, descCtt );

    alpha = 1.0;
    beta = 1.0;
    // pdgeadd does not work properly if A and C are same matrices. We have to use Ctt for out of place addition.
    pXgeadd_( "T", &nmeasts, &nmeasts, 
              &alpha, 
              Ctt, &i1, &i1, descCtt,
              &beta, 
              Cee, &i1, &i1, descCee );
    *flops += FLOPS_XGEADD(nmeasts, nmeasts);

    alpha = 1.;
    beta = 1.;
    pXgemm_( "N", "T", &nmeasts, &nmeasts, &nmeas, &alpha,
             Ctm, &i1, &i1, descCtm,
             R, &i1, &i1, descR,
             &beta,
             Cee, &i1, &i1, descCee );
    *flops += FLOPS_XGEMM(nmeasts, nmeasts, nmeas);

    alpha = 1.;
    beta = 0.;
    pXsymm_( "R", "L", &nact, &nmeasts, &alpha,
             Cee, &i1, &i1, descCee,
             Dx, &i1, &i1, descDx,
             &beta,
             Tmp, &i1, &i1, descTmp );
    *flops += FLOPS_XSYMM('R', nact, nmeasts);

    alpha = 1.;
    beta = 0.;
    pXgemm_( "N", "T", &nact, &nact, &nmeasts, &alpha,
             Tmp, &i1, &i1, descTmp,
             Dx, &i1, &i1, descDx,
             &beta,
             Cvv, &i1, &i1, descCvv );
    *flops += FLOPS_XGEMM(nact, nact, nmeasts);

    return 0;
}
/**
* Addition of 2 matrices
* Return C = beta*C + alpha*op(A)
*
* A    :(M,N)
* op(C):(M,N)
*
* op(C)=C            if trans='N'
*       transpose(C) if Trans='T'
*
*/
void Xgeadd(char trans, int M, int N, real_t alpha, real_t *A, int ldA, real_t beta, real_t *B, int ldB){

    int i,j;

    if(trans=='N'){
        for(i=0; i<M;i++){
            for(j=0; j<N;j++){
                B[i+j*ldB] = beta*B[i+j*ldB] + alpha*A[i+j*ldA];
            }
        }
    }
    else if(trans=='T'){
        for(i=0; i<M;i++){
            for(j=0; j<N;j++){
                B[i*ldB+j] = beta*B[i+j*ldB] + alpha*A[i*ldA+j];
            }
        }
    }
}

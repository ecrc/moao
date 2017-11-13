/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include <math.h>
#include <stdlib.h>
#include <starpu.h>
#include <morse.h>
#include "dscaldiag.h"

#define RTBLKADDR( desc, m, n ) ( (starpu_data_handle_t)RUNTIME_desc_getaddr( desc, m, n ) )

enum blas_order_type {
            blas_rowmajor = 101,
            blas_colmajor = 102 };

enum blas_cmach_type {
            blas_base      = 151,
            blas_t         = 152,
            blas_rnd       = 153,
            blas_ieee      = 154,
            blas_emin      = 155,
            blas_emax      = 156,
            blas_eps       = 157,
            blas_prec      = 158,
            blas_underflow = 159,
            blas_overflow  = 160,
            blas_sfmin     = 161};

enum blas_norm_type {
            blas_one_norm       = 171,
            blas_real_one_norm  = 172,
            blas_two_norm       = 173,
            blas_frobenius_norm = 174,
            blas_inf_norm       = 175,
            blas_real_inf_norm  = 176,
            blas_max_norm       = 177,
            blas_real_max_norm  = 178 };

static
double
BLAS_dpow_di(double x, int n) {
  double rv = 1.0;

  if (n < 0) {
    n = -n;
    x = 1.0 / x;
  }

  for (; n; n >>= 1, x *= x) {
    if (n & 1)
      rv *= x;
  }

  return rv;
}



static void
BLAS_error(char *rname, int err, int val, int x) {
  fprintf( stderr, "%s %d %d %d\n", rname, err, val, x );
  //abort();
}

static
double
BLAS_dfpinfo(enum blas_cmach_type cmach) {
  double eps = 1.0, r = 1.0, o = 1.0, b = 2.0;
  int t = 53, l = 1024, m = -1021;
  char rname[] = "BLAS_dfpinfo";

  if ((sizeof eps) == sizeof(float)) {
    t = 24;
    l = 128;
    m = -125;
  } else {
    t = 53; 
    l = 1024;
    m = -1021;
  }

  /* for (i = 0; i < t; ++i) eps *= half; */
  eps = BLAS_dpow_di( b, -t );
  /* for (i = 0; i >= m; --i) r *= half; */
  r = BLAS_dpow_di( b, m-1 );

  o -= eps;
  /* for (i = 0; i < l; ++i) o *= b; */
  o = (o * BLAS_dpow_di( b, l-1 )) * b;

  switch (cmach) {
    case blas_eps: return eps;
    case blas_sfmin: return r;
    default:
      BLAS_error( rname, -1, cmach, 0 );
      break;
  }
  return 0.0;
}


static
void
BLAS_zge_norm(enum blas_order_type order, enum blas_norm_type norm,
  int m, int n, const double *a, int lda, double *res) {
  int i, j; float anorm, v;
  char rname[] = "BLAS_zge_norm";

  if (order != blas_colmajor) BLAS_error( rname, -1, order, 0 );

  if (norm == blas_frobenius_norm) {
    anorm = 0.0f;
    for (j = n; j; --j) {
      for (i = m; i; --i) {
        v = a[0];
        anorm += v * v;
        a++;
      }
      a += lda - m;
    }
    anorm = sqrt( anorm );
  } else if (norm == blas_inf_norm) {
    anorm = 0.0f;
    for (i = 0; i < m; ++i) {
      v = 0.0f;
      for (j = 0; j < n; ++j) {
        v += fabs( a[i + j * lda] );
      }
      if (v > anorm)
        anorm = v;
    }
  } else {
    BLAS_error( rname, -2, norm, 0 );
    return;
  }

  if (res) *res = anorm;
}

static void BLAS_zge_norm_func(void *descr[], void *cl_arg){
  double *data, *res;
  double v;
  int m,n,lda;
  int i,j;

  data=(double *)STARPU_MATRIX_GET_PTR(descr[0]);
  res =(double *)STARPU_MATRIX_GET_PTR(descr[1]);
  starpu_codelet_unpack_args(cl_arg,&m,
                                    &n,
                                    &lda
                            );

  res[0]=0.0f;
  for (i = 0; i < m; ++i) {
    v = 0.0f;
    for (j = 0; j < n; ++j) {
      v += fabs( data[i + j * lda] );
    }
    if (v > res[0])
      res[0] = v;
  }


}

static struct starpu_codelet cl_BLAS_zge_norm=
{
    .where=STARPU_CPU,
    .cpu_funcs={BLAS_zge_norm_func},
    .nbuffers=2
};

BLAS_zge_norm_Tile(enum blas_order_type order, enum blas_norm_type norm,
  MORSE_desc_t *descA, double *res, MORSE_sequence_t *sequence){

  int i, j,m,n,lda; float v;
  char rname[] = "BLAS_zge_norm";

  MORSE_desc_t *descR=NULL;
  MORSE_Desc_Create(&descR,NULL,MorseRealDouble, 1, 1, 1, descA->mt*descA->nt, 1, 0, 0, descA->mt*descA->nt, 1,1,1);


  if (order != blas_colmajor) BLAS_error( rname, -1, order, 0 );

    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();

    struct starpu_codelet *codelet=&cl_BLAS_zge_norm;    
    for(i=0;i<descA->mt;i++){
      m = i == descA->mt-1 ? descA->m-i*descA->mb : descA->mb;
      lda=descA->get_blkldd(descA, i);
      for(j=0;j<descA->nt;j++){
        n = j == descA->nt-1 ? descA->n-j*descA->nb : descA->nb;
        starpu_insert_task(codelet,
                            STARPU_R    ,RTBLKADDR(descA,i,j),
                            STARPU_W    ,RTBLKADDR(descR,i*descA->mt+j,0),
                            STARPU_VALUE,&m,sizeof(int),
                            STARPU_VALUE,&n,sizeof(int),
                            STARPU_VALUE,&lda,sizeof(int),
                            0);
        
      }
    }

  MORSE_Sequence_Wait(sequence);
  res[0]=0.0f;
  for(i=0;i<descA->mt*descA->nt;i++){
    res[0]+=((double*)(descR->mat))[i];
  }
  //if (res) *res = anorm;
  MORSE_Desc_Destroy(&descR);

}



/*------------------------------------------------------------------------
 *  Check the factorization of the matrix A2
 */
//static int check_factorization(int N, MORSE_Complex64_t *A1, MORSE_Complex64_t *A2, int LDA, int uplo, double eps)
static int check_factorization(MORSE_desc_t *A1, MORSE_desc_t *A2, int uplo, double eps, MORSE_sequence_t *sequence)
{
    printf("============\n");
    fflush(stdout);
    double Anorm, Rnorm;
    double alpha;
    int info_factorization;

    int N=A1->lm;
    int ts=A1->mb;

    MORSE_desc_t *Residual=NULL;
    MORSE_desc_t *L1=NULL;
    MORSE_desc_t *L2=NULL;
    MORSE_Desc_Create(&Residual,NULL,MorseRealDouble, ts, ts, ts*ts, N, N, 0, 0, N, N,1,1);
    MORSE_Desc_Create(&L1,NULL,MorseRealDouble, ts, ts, ts*ts, N, N, 0, 0, N, N,1,1);
    MORSE_Desc_Create(&L2,NULL,MorseRealDouble, ts, ts, ts*ts, N, N, 0, 0, N, N,1,1);


    MORSE_dlaset_Tile(MorseUpperLower, 0.,0.,Residual);
    MORSE_dlaset_Tile(MorseUpperLower, 0.,0.,L1);
    MORSE_dlaset_Tile(MorseUpperLower, 0.,0.,L2);

    alpha= 1.0;

    MORSE_dlacpy_Tile( uplo, A1, Residual);


    BLAS_zge_norm( blas_colmajor, blas_inf_norm, N, N, Residual->mat, N, &Anorm );
    //BLAS_zge_norm_Tile( blas_colmajor, blas_inf_norm,Residual, &Anorm ,sequence);

    /* Dealing with L'L or U'U  */
    if (uplo == MorseUpper){
        MORSE_dlacpy_Tile( MorseUpper, A2, L1);
        MORSE_dlacpy_Tile( MorseUpper, A2, L2);
        MORSE_dtrmm_Tile(MorseLeft, MorseUpper, MorseTrans,MorseNonUnit,alpha,L1,L2);
    }
    else{
        MORSE_dlacpy_Tile( MorseLower, A2, L1);
        MORSE_dlacpy_Tile( MorseLower, A2, L2);

        MORSE_dtrmm_Tile(MorseRight, MorseLower, MorseTrans,MorseNonUnit,alpha,L1,L2);

    }
    /* Compute the Residual || A -L'L|| */
    MORSE_dgeadd_Tile(MorseNoTrans, -1,L2, 1,Residual);

    BLAS_zge_norm( blas_colmajor, blas_inf_norm, N, N, Residual->mat, N, &Rnorm );
    //BLAS_zge_norm_Tile( blas_colmajor, blas_inf_norm,Residual, &Rnorm ,sequence);

    printf("Checking the Cholesky Factorization \n");
    printf("N=%d  ||A||=%e -- ||L'L-A||_oo/(||A||_oo.N.eps) = %e \n",N,Anorm,Rnorm/(Anorm*N*eps));

    if ( isnan(Rnorm/(Anorm*N*eps)) || isinf(Rnorm/(Anorm*N*eps)) || (Rnorm/(Anorm*N*eps) > 60.0) ){
        printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else{
        printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    free(Residual); free(L1); free(L2);

    return info_factorization;
}


/*------------------------------------------------------------------------
 *  Check the accuracy of the solution of the linear system
 */

static int check_solution(MORSE_desc_t *A1, MORSE_desc_t *B1, MORSE_desc_t *B2, double eps, MORSE_sequence_t *sequence )
{
    int info_solution;
    double Rnorm, Anorm, Xnorm, Bnorm, result;
    double alpha, beta;

    alpha = 1.0;
    beta  = -1.0;

    int N=A1->lm;
    int NRHS=B1->lm;
    int LDA=A1->ln;
    int LDB=B1->lm;
    int ts=B2->mb;

    MORSE_desc_t *B3=NULL;
    MORSE_Desc_Create(&B3,NULL,MorseRealDouble, ts, ts, ts*ts, B2->lm, B2->ln, 0, 0,  B2->lm, B2->ln,1,1);
    MORSE_dlacpy_Tile( MorseUpperLower, B2, B3);

    BLAS_zge_norm( blas_colmajor, blas_inf_norm, NRHS, N, B3->mat, LDB, &Xnorm );
    BLAS_zge_norm( blas_colmajor, blas_inf_norm, N, N,    A1->mat, LDA, &Anorm );
    BLAS_zge_norm( blas_colmajor, blas_inf_norm, NRHS, N, B1->mat, LDB, &Bnorm );
    //BLAS_zge_norm_Tile( blas_colmajor, blas_inf_norm,B3, &Xnorm ,sequence);
    //BLAS_zge_norm_Tile( blas_colmajor, blas_inf_norm,A1, &Anorm ,sequence);
    //BLAS_zge_norm_Tile( blas_colmajor, blas_inf_norm,B1, &Bnorm ,sequence);

    //cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, NRHS, N, CBLAS_SADDR(alpha), A1, LDA, B3, LDB, CBLAS_SADDR(beta), B1, LDB);
    MORSE_dtrmm_Tile(MorseRight,MorseLower,MorseNoTrans,MorseNonUnit,alpha,A1,B3);
    MORSE_dgeadd_Tile(MorseNoTrans, alpha,B3, beta,B1);
    BLAS_zge_norm( blas_colmajor, blas_inf_norm, N, NRHS, B1->mat, LDB, &Rnorm );
    //BLAS_zge_norm_Tile( blas_colmajor, blas_inf_norm,B1, &Rnorm ,sequence);

    if (getenv("MORSE_TESTING_VERBOSE"))
      printf( "||A||_oo=%f\n||X||_oo=%f\n||B||_oo=%f\n||X A - B||_oo=%e\n", Anorm, Xnorm, Bnorm, Rnorm );

    result = Rnorm / ( (Anorm*Xnorm+Bnorm)*N*eps ) ;
    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||xA-B||_oo/((||x||_oo||A||_oo+||B||_oo).N.eps) = %e \n", result);

    if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
     }
    else{
        printf("-- The solution is CORRECT ! \n");
        info_solution = 0;
    }


    return info_solution;
}


#ifndef LAPACK_TEMPLATES_H
#define LAPACK_TEMPLATES_H

#ifdef USE_MKL
#include <mkl_lapacke.h>
#include <mkl_cblas.h>
#else
#include <lapacke.h>
#include <cblas.h>
#endif

template<typename> struct LAPACK;

template<>
struct LAPACK <float>
{
    static inline lapack_int potrf(int matrix_order, char uplo, lapack_int n, float* a, lapack_int lda){
        return  LAPACKE_spotrf( matrix_order, uplo, n, a, lda );
    }

    static inline void lacpy( int matrix_order, char uplo, lapack_int m,
                           lapack_int n, const float* a, lapack_int lda, float* b,
                           lapack_int ldb ){
        LAPACKE_slacpy( matrix_order, uplo, m, n, a, lda, b, ldb );
    }

    static inline lapack_int laset( int matrix_order, char uplo, lapack_int m,
                           lapack_int n, float alpha, float beta, float* a,
                           lapack_int lda){
        return LAPACKE_slaset( matrix_order, uplo, m, n, alpha, beta, a, lda );
    }

    static inline void trsm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_SIDE Side, 
                        OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                        OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
                        OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda, 
                        float *B, OPENBLAS_CONST blasint ldb){
        return  cblas_strsm(Order, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
    }

    static inline void syr2k(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                            OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                            OPENBLAS_CONST blasint K, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A,
                            OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *B, OPENBLAS_CONST blasint ldb, 
                            OPENBLAS_CONST float beta, float *C, OPENBLAS_CONST blasint ldc){
        cblas_ssyr2k( Order,  Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C,  ldc);
    }

    static inline void symm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_SIDE Side,
                                    OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
                                    OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A, OPENBLAS_CONST blasint lda,
                                    OPENBLAS_CONST float *B, OPENBLAS_CONST blasint ldb, OPENBLAS_CONST float beta,
                                    float *C, OPENBLAS_CONST blasint ldc){
        cblas_ssymm( Order, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
    }

    static inline void gemm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                                OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
                                OPENBLAS_CONST blasint K, OPENBLAS_CONST float alpha, OPENBLAS_CONST float *A, 
                                OPENBLAS_CONST blasint lda, OPENBLAS_CONST float *B, OPENBLAS_CONST blasint ldb, 
                                OPENBLAS_CONST float beta, float *C, OPENBLAS_CONST blasint ldc){
         cblas_sgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
    }

    static inline void scal(OPENBLAS_CONST blasint N, OPENBLAS_CONST float alpha, float *X, OPENBLAS_CONST blasint incX){
        cblas_sscal(N, alpha, X, incX);
    }

};
template<>
struct LAPACK<double>
{
    static inline lapack_int potrf(int matrix_order, char uplo, lapack_int n, double* a, lapack_int lda){
        return  LAPACKE_dpotrf( matrix_order, uplo, n, a, lda );
    }

    static inline void lacpy( int matrix_order, char uplo, lapack_int m,
                           lapack_int n, const double* a, lapack_int lda, double* b,
                           lapack_int ldb){
        LAPACKE_dlacpy( matrix_order, uplo, m, n, a, lda, b, ldb );
    }

    static inline lapack_int laset( int matrix_order, char uplo, lapack_int m,
                           lapack_int n, double alpha, double beta, double* a,
                           lapack_int lda ){
        return LAPACKE_dlaset( matrix_order,uplo, m, n, alpha, beta, a, lda );
    }

    static inline void trsm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_SIDE Side, 
                        OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                        OPENBLAS_CONST enum CBLAS_DIAG Diag, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
                        OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda, 
                        double *B, OPENBLAS_CONST blasint ldb){
        return  cblas_dtrsm(Order, Side, Uplo, TransA, Diag, M, N, alpha, A, lda, B, ldb);
    }

    static inline void syr2k(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_UPLO Uplo,
                            OPENBLAS_CONST enum CBLAS_TRANSPOSE Trans, OPENBLAS_CONST blasint N,
                            OPENBLAS_CONST blasint K, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A,
                            OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *B, OPENBLAS_CONST blasint ldb, 
                            OPENBLAS_CONST double beta, double *C, OPENBLAS_CONST blasint ldc){
        cblas_dsyr2k( Order,  Uplo, Trans, N, K, alpha, A, lda, B, ldb, beta, C,  ldc);
    }

    static inline void symm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_SIDE Side,
                                    OPENBLAS_CONST enum CBLAS_UPLO Uplo, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
                                    OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A, OPENBLAS_CONST blasint lda,
                                    OPENBLAS_CONST double *B, OPENBLAS_CONST blasint ldb, OPENBLAS_CONST double beta,
                                    double *C, OPENBLAS_CONST blasint ldc){
        cblas_dsymm( Order, Side, Uplo, M, N, alpha, A, lda, B, ldb, beta, C, ldc);
    }

    static inline void gemm(OPENBLAS_CONST enum CBLAS_ORDER Order, OPENBLAS_CONST enum CBLAS_TRANSPOSE TransA,
                                OPENBLAS_CONST enum CBLAS_TRANSPOSE TransB, OPENBLAS_CONST blasint M, OPENBLAS_CONST blasint N,
                                OPENBLAS_CONST blasint K, OPENBLAS_CONST double alpha, OPENBLAS_CONST double *A, 
                                OPENBLAS_CONST blasint lda, OPENBLAS_CONST double *B, OPENBLAS_CONST blasint ldb, 
                                OPENBLAS_CONST double beta, double *C, OPENBLAS_CONST blasint ldc){
         cblas_dgemm(Order, TransA, TransB, M, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
    }

    static inline void scal(OPENBLAS_CONST blasint N, OPENBLAS_CONST double alpha, double *X, OPENBLAS_CONST blasint incX){
        cblas_dscal(N, alpha, X, incX);
    }

};

#endif // LAPACK_TEMPLATES_H

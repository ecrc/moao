#ifndef CHAMELEON_TEMPLATES_H
#define CHAMELEON_TEMPLATES_H

#include <morse.h>
#include "coreblas/lapacke.h"
#include "coreblas/cblas.h"

template<typename> struct MORSE;

template<>
struct MORSE <float>{
    static const int real = MorseRealFloat;
    static const int int64 = MorseRealDouble;
    static inline void cblas_Xscal(const int N, const float alpha, float *X, const int incX){
        cblas_sscal(N, alpha, X, incX);
    }

    static inline int potrf_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence,
                            MORSE_request_t *request){
    MORSE_spotrf_Tile_Async( uplo, A, sequence, request);
    }

    static inline int trsm_Tile_Async(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, 
                            float alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, 
                            MORSE_request_t *request){
         MORSE_strsm_Tile_Async(side, uplo, transA, diag, alpha, A, B, sequence, request);
    }

    static inline int syr2k_Tile_Async(MORSE_enum uplo, MORSE_enum trans, float alpha, MORSE_desc_t *A, MORSE_desc_t *B, float beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request){
        MORSE_ssyr2k_Tile_Async(uplo, trans, alpha, A, B, beta, C, sequence, request);
    }

    static inline int symm_Tile_Async(MORSE_enum side, MORSE_enum uplo, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request){
        MORSE_ssymm_Tile_Async(side, uplo, alpha, A, B, beta, C, sequence, request);
    }

    static inline int lacpy_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request){
        MORSE_slacpy_Tile_Async(uplo, A, B, sequence, request);
    }


    static inline int tradd_Tile_Async(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request){
        MORSE_stradd_Tile_Async(uplo, trans, alpha, A,  beta, B, sequence, request);
    }

    static inline int gemm_Tile_Async(MORSE_enum transA, MORSE_enum transB, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request){
        MORSE_sgemm_Tile_Async(transA, transB, alpha, A, B,  beta, C, sequence, request);
    }
};
template<>
struct MORSE <double>{
    static const int real = MorseRealDouble;
    static const int int64 = MorseRealDouble;

    static inline void cblas_Xscal(const int N, const double alpha, double *X, const int incX){
        cblas_dscal( N,  alpha, X, incX);
    }
    static inline int potrf_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence,
                            MORSE_request_t *request){
    MORSE_dpotrf_Tile_Async( uplo, A, sequence, request);
    }

    static inline int trsm_Tile_Async(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, 
                            double alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence,
                            MORSE_request_t *request){
         MORSE_dtrsm_Tile_Async(side, uplo, transA, diag, alpha, A, B, sequence, request);
    }

    static inline int syr2k_Tile_Async(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request){
        MORSE_dsyr2k_Tile_Async(uplo, trans, alpha, A, B, beta, C, sequence, request);
    }

    static inline int symm_Tile_Async(MORSE_enum side, MORSE_enum uplo, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request){
        MORSE_dsymm_Tile_Async(side, uplo, alpha, A, B, beta, C, sequence, request);
    }

    static inline int lacpy_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request){
        MORSE_dlacpy_Tile_Async(uplo, A, B, sequence, request);
    }


    static inline int tradd_Tile_Async(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request){
        MORSE_dtradd_Tile_Async(uplo, trans, alpha, A,  beta, B, sequence, request);
    }

    static inline int gemm_Tile_Async(MORSE_enum transA, MORSE_enum transB, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request){
        MORSE_dgemm_Tile_Async(transA, transB, alpha, A, B,  beta, C, sequence, request);
    }
};

#endif // CHAMELEON_TEMPLATES_H

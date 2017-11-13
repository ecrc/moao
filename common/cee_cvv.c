/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "cee_cvv.h"
#include "utils.h"
#include "dscaldiag.h"

void MORSE_Cee_Cvv(MORSE_desc_t *Cmm, MORSE_desc_t *Cpp, MORSE_desc_t *Cpm, MORSE_desc_t *R, MORSE_desc_t *Dx,
                   MORSE_desc_t *Cee, MORSE_desc_t *Cvv, MORSE_desc_t *Tmp, MORSE_sequence_t *sequence){
    MORSE_request_t request[8];
    double alpha,beta;
    alpha = -1.0;
    beta = 1.0;
    HANDLE_ERROR( MORSE_dsyr2k_Tile_Async(MorseLower, MorseNoTrans,alpha, Cpm, R,beta, Cpp,sequence, &request[0]),
                "MORSE_dsyr2k_Tile_Async");
    MORSE_Sequence_Wait(sequence);

    alpha = 1.0;
    beta = 0.0;
    // Cpm gets overwritten at this point and it is used as a buffer
    HANDLE_ERROR(MORSE_dsymm_Tile_Async(MorseRight, MorseLower,alpha, Cmm, R,beta,
                    Cpm,sequence, &request[1]),
                 "MORSE_dsymm_Tile_Async");
    MORSE_Sequence_Wait(sequence);

    double *yTmp=(double*)malloc(Cee->m*Cee->n*sizeof(double));
        // Rebuild the symmetry for Cpp
    HANDLE_ERROR(MORSE_dlacpy_Tile_Async( MorseLower, Cpp, Cee,sequence, &request[2]),
                     "PLASMA_dlacpy_Tile_Async");
    MORSE_Sequence_Wait(sequence);

    HANDLE_ERROR(MORSE_dtradd_Tile_Async( MorseUpper, MorseTrans, alpha, Cee, alpha, Cee, sequence, &request[3] ),
                "PLASMA_dgeadd_Tile_Async");
    MORSE_Sequence_Wait(sequence);
    
    alpha = 0.5;
    HANDLE_ERROR(APPNAME_dscaldiag(alpha, Cee, sequence, &request[4]),
                "PLASMA_dscaldiag_Tile_Async");
    MORSE_Sequence_Wait(sequence);
MORSE_Tile_to_Lapack(Cee,yTmp,Cee->m);

free(yTmp);


    alpha = 1.0;
    beta = 1.0;
    HANDLE_ERROR(MORSE_dgemm_Tile_Async(MorseNoTrans, MorseTrans, alpha, Cpm, R, beta, Cee, sequence, &request[5]),
                "PLASMA_dgemm_Tile_Async");
    MORSE_Sequence_Wait(sequence);

    alpha = 1.0;
    beta = 0.0;
    HANDLE_ERROR(MORSE_dsymm_Tile_Async(MorseRight, MorseLower, alpha, Cee, Dx,beta,  Tmp, sequence, &request[6]),
                "PLASMA_dsymm_Tile_Async");
    MORSE_Sequence_Wait(sequence);

    alpha = 1.0;
    beta = 0.0;
    HANDLE_ERROR(MORSE_dgemm_Tile_Async(MorseNoTrans,MorseTrans, alpha, Tmp, Dx, beta, Cvv, sequence, &request[7]),
                "PLASMA_dgemm_Tile_Async");
    MORSE_Sequence_Wait(sequence);

}


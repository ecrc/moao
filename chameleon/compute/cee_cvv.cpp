/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "cee_cvv.hpp"
//#include "utils.h"
#include "dscaldiag.hpp"
#include "chameleon_templates.hpp"

#ifdef USE_CHAMELEON
#include "stdio.h"
#define HANDLE_ERROR(_result, _func)               \
  if(MORSE_SUCCESS != _result) {                  \
    fprintf(stderr,"An error occurred (%d), when calling %s at line: %d... \n\n", _result, _func, __LINE__ -1);   \
  }  //break;

#define HANDLE_ERROR_RET(_result, _func)               \
  if(MORSE_SUCCESS != _result) {                  \
    fprintf(stderr,"An error occurred (%d), when calling %s at line: %d... \n\n", _result, _func, __LINE__ -1);   \
    exit(-1);                                         \
  }
#endif

#include <iostream>
template<typename T>
void MOAO_CHAMELEON::Cee_Cvv(MORSE_desc_t *Cmm, MORSE_desc_t *Cpp, MORSE_desc_t *Cpm, MORSE_desc_t *R, MORSE_desc_t *Dx,
                   MORSE_desc_t *Cee, MORSE_desc_t *Cvv, MORSE_desc_t *Tmp, MORSE_sequence_t *sequence,int sync){
    MORSE_request_t request[8];
    T alpha,beta;
    alpha = -1.0;
    beta = 1.0;
    HANDLE_ERROR( MORSE<T>::syr2k_Tile_Async(MorseLower, MorseNoTrans,alpha, Cpm, R,beta, Cpp,sequence, &request[0]),
                "MORSE_Xsyr2k_Tile_Async");
    if(sync)MORSE_Sequence_Wait(sequence);

    alpha = 1.0;
    beta = 0.0;
    // Cpm gets overwritten at this point and it is used as a buffer
    HANDLE_ERROR(MORSE<T>::symm_Tile_Async(MorseRight, MorseLower,alpha, Cmm, R,beta,
                    Cpm,sequence, &request[1]),
                 "MORSE_Xsymm_Tile_Async");
    if(sync)MORSE_Sequence_Wait(sequence);

    // Rebuild the symmetry for Cpp
    HANDLE_ERROR(MORSE<T>::lacpy_Tile_Async( MorseLower, Cpp, Cee,sequence, &request[2]),
                     "MORSE_Xlacpy_Tile_Async");
    if(sync)MORSE_Sequence_Wait(sequence);

    //necessary barrier (do not remove)
    MORSE_Sequence_Wait(sequence);
    HANDLE_ERROR(MORSE<T>::tradd_Tile_Async( MorseUpper, MorseTrans, alpha, Cee, alpha, Cee, sequence, &request[3] ),
                "MORSE_Xtradd_Tile_Async");
    if(sync)MORSE_Sequence_Wait(sequence);
    
    MORSE_Sequence_Wait(sequence);
    alpha = 0.5;
    HANDLE_ERROR(MOAO_CHAMELEON::Xscaldiag<T>(alpha, Cee, sequence, &request[4]),
                "MORSE_Xscaldiag_Tile_Async");
    MORSE_Sequence_Wait(sequence);


    alpha = 1.0;
    beta = 1.0;
    HANDLE_ERROR(MORSE<T>::gemm_Tile_Async(MorseNoTrans, MorseTrans, alpha, Cpm, R, beta, Cee, sequence, &request[5]),
                "MORSE_Xgemm_Tile_Async");
    if(sync)MORSE_Sequence_Wait(sequence);

    alpha = 1.0;
    beta = 0.0;
    HANDLE_ERROR(MORSE<T>::symm_Tile_Async(MorseRight, MorseLower, alpha, Cee, Dx,beta,  Tmp, sequence, &request[6]),
                "MORSE_Xsymm_Tile_Async");

    alpha = 1.0;
    beta = 0.0;
    HANDLE_ERROR(MORSE<T>::gemm_Tile_Async(MorseNoTrans,MorseTrans, alpha, Tmp, Dx, beta, Cvv, sequence, &request[7]),
                "MORSE_Xgemm_Tile_Async");
    if(sync)MORSE_Sequence_Wait(sequence);

}
template void MOAO_CHAMELEON::Cee_Cvv<float>(MORSE_desc_t *Cmm, MORSE_desc_t *Cpp, MORSE_desc_t *Cpm, MORSE_desc_t *R, MORSE_desc_t *Dx,
                   MORSE_desc_t *Cee, MORSE_desc_t *Cvv, MORSE_desc_t *Tmp, MORSE_sequence_t *sequence,int sync);
template void MOAO_CHAMELEON::Cee_Cvv<double>(MORSE_desc_t *Cmm, MORSE_desc_t *Cpp, MORSE_desc_t *Cpm, MORSE_desc_t *R, MORSE_desc_t *Dx,
                   MORSE_desc_t *Cee, MORSE_desc_t *Cvv, MORSE_desc_t *Tmp, MORSE_sequence_t *sequence, int sync);


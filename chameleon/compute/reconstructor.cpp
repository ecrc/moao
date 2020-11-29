/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "reconstructor.hpp"
#include "chameleon_templates.hpp"
//#include "utils.h"

//#include "checking_solve.c"

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


template<typename T>
void MOAO_CHAMELEON::reconstructor(MORSE_desc_t *Cmm, MORSE_desc_t *Ctm, MORSE_sequence_t *sequence,int check, int sync){


    MORSE_desc_t *CmmCpy=NULL;
    MORSE_desc_t *CtmCpy=NULL;
    //TODO
    //if(check>0){
    //    MORSE_Desc_Create(&CmmCpy,NULL,MorseRealDouble, Cmm->mb, Cmm->nb, Cmm->mb*Cmm->nb, Cmm->lm, Cmm->ln, 0, 0, Cmm->lm, Cmm->ln,1,1);
    //    MORSE_Desc_Create(&CtmCpy,NULL,MorseRealDouble, Ctm->mb, Ctm->nb, Ctm->mb*Ctm->nb, Ctm->lm, Ctm->ln, 0, 0, Ctm->lm, Ctm->ln,1,1);
    //    MORSE_dlacpy_Tile( MorseUpperLower, Cmm, CmmCpy);
    //    MORSE_dlacpy_Tile( MorseUpperLower, Ctm, CmmCpy);
    //}

    T alpha = 1.0;
    MORSE_request_t request[3] = { MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER};
    HANDLE_ERROR(MORSE<T>::potrf_Tile_Async(MorseLower,Cmm,sequence,&request[0]),"MORSE_Xpotrf_Async");
    if(sync)MORSE_Sequence_Wait(sequence);

    HANDLE_ERROR(MORSE<T>::trsm_Tile_Async(MorseRight, MorseLower, MorseTrans, MorseNonUnit,alpha, Cmm, Ctm, sequence, &request[1]),
            "MORSE_Xtrsm_Async");
    if(sync)MORSE_Sequence_Wait(sequence);

    HANDLE_ERROR(MORSE<T>::trsm_Tile_Async(MorseRight, MorseLower, MorseNoTrans, MorseNonUnit, alpha, Cmm, Ctm, sequence, &request[2]),
            "MORSE_Xtrsm_Async");
    if(sync)MORSE_Sequence_Wait(sequence);

    //TODO
    //if(check>0){
    //    double eps=BLAS_dfpinfo( blas_eps );
    //    printf("eps: %e\ncall checking the Cholesky Factorization \n",eps);fflush(stdout);
    //    check_factorization(CmmCpy,Cmm,MorseLower,eps,sequence);
    //    printf("call checking the Cholesky solve \n");fflush(stdout);
    //    check_solution(CmmCpy,CtmCpy,Ctm,eps,sequence);
    //}
}
template void MOAO_CHAMELEON::reconstructor<float>(MORSE_desc_t *Cmm, MORSE_desc_t *Ctm, MORSE_sequence_t *sequence,int check, int sync);
template void MOAO_CHAMELEON::reconstructor<double>(MORSE_desc_t *Cmm, MORSE_desc_t *Ctm, MORSE_sequence_t *sequence,int check, int sync);

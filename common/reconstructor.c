/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "reconstructor.h"
#include "utils.h"

#include "checking_solve.c"

void APPNAME_reconstructor(MORSE_desc_t *Cmm, MORSE_desc_t *Ctm, MORSE_sequence_t *sequence,int check){


    MORSE_desc_t *CmmCpy=NULL;
    MORSE_desc_t *CtmCpy=NULL;
    if(check>0){
        MORSE_Desc_Create(&CmmCpy,NULL,MorseRealDouble, Cmm->mb, Cmm->nb, Cmm->mb*Cmm->nb, Cmm->lm, Cmm->ln, 0, 0, Cmm->lm, Cmm->ln,1,1);
        MORSE_Desc_Create(&CtmCpy,NULL,MorseRealDouble, Ctm->mb, Ctm->nb, Ctm->mb*Ctm->nb, Ctm->lm, Ctm->ln, 0, 0, Ctm->lm, Ctm->ln,1,1);
        MORSE_dlacpy_Tile( MorseUpperLower, Cmm, CmmCpy);
        MORSE_dlacpy_Tile( MorseUpperLower, Ctm, CmmCpy);
    }

    double alpha = 1.0;
    MORSE_request_t request[3] = { MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER};
    HANDLE_ERROR(MORSE_dpotrf_Tile_Async(MorseLower,Cmm,sequence,&request[0]),"MORSE_dpotrf_Async");
    MORSE_Sequence_Wait(sequence);

    HANDLE_ERROR(MORSE_dtrsm_Tile_Async(MorseRight, MorseLower, MorseTrans, MorseNonUnit,alpha, Cmm, Ctm, sequence, &request[1]),
            "MORSE_dtrsm_Async");
    MORSE_Sequence_Wait(sequence);

    HANDLE_ERROR(MORSE_dtrsm_Tile_Async(MorseRight, MorseLower, MorseNoTrans, MorseNonUnit, alpha, Cmm, Ctm, sequence, &request[2]),
            "MORSE_dtrsm_Async");
    MORSE_Sequence_Wait(sequence);

    if(check>0){
        double eps=BLAS_dfpinfo( blas_eps );
        printf("eps: %e\ncall checking the Cholesky Factorization \n",eps);fflush(stdout);
        check_factorization(CmmCpy,Cmm,MorseLower,eps,sequence);
        printf("call checking the Cholesky solve \n");fflush(stdout);
        check_solution(CmmCpy,CtmCpy,Ctm,eps,sequence);
    }
}

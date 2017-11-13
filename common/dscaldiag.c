/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "dscaldiag.h"
#include "codelet_dscaldiag.h"
int APPNAME_dscaldiag(double alpha, MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request) {


    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    int k;
    int ldak;
    int tempkk;
    MORSE_desc_t A = *descA;
    //double *data;

   // struct starpu_codelet *cl=&cl_dscaldiag;
    for (k = 0; k < A.nt; k++) {
        tempkk = k == A.nt-1 ? A.m-k*A.mb : A.mb;
        ldak = descA->get_blkldd(descA, k);
        APPNAME_TASK_dscaldiag(tempkk,alpha,descA,k,ldak);
        //double *data=(double*)morse_getaddr_ccrb(descA,k,k);
        //starpu_insert_task(cl,
        //        STARPU_VALUE, &alpha,  sizeof(double),
        //        STARPU_VALUE, &tempkk, sizeof(int),
        //        STARPU_VALUE, &data,   sizeof(double*),
        //        STARPU_VALUE, &ldak,   sizeof(int),
        //        0);
    }    
    RUNTIME_options_finalize(&options, morse);
    MORSE_TASK_dataflush_all();
    return MORSE_SUCCESS;
}


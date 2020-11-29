/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "dscaldiag.hpp"
#include "codelet_dscaldiag.hpp"
#include "context.h"
template<typename T>
int MOAO_CHAMELEON::Xscaldiag(T alpha, MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request) {


    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return sequence->status;
    RUNTIME_options_init(&options, morse, sequence, request);

    int k;
    int ldak;
    int tempkk;
    MORSE_desc_t A = *descA;

   // struct starpu_codelet *cl=&cl_dscaldiag;
    for (k = 0; k < A.nt; k++) {
        tempkk = k == A.nt-1 ? A.m-k*A.mb : A.mb;
        ldak = descA->get_blkldd(descA, k);
        MOAO_STARPU::TASK_Xscaldiag(tempkk,alpha,descA,k,ldak);
    }    
    RUNTIME_options_finalize(&options, morse);
    MORSE_Desc_Flush(descA,sequence);
    return MORSE_SUCCESS;
}
template int MOAO_CHAMELEON::Xscaldiag<float>(float alpha, MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request);
template int MOAO_CHAMELEON::Xscaldiag<double>(double alpha, MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request);

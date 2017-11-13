/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "codelet_dscaldiag.h"
#include <starpu.h>

static struct starpu_codelet cl_dscaldiag =
{
    .where = STARPU_CPU,
    .cpu_funcs = {cl_dscaldiag_cpu_func},
    .nbuffers = 0
};

int APPNAME_TASK_dscaldiag(int M,double alpha,MORSE_desc_t *descA, int k, int lda){

    struct starpu_codelet *cl=&cl_dscaldiag;

    double *data=(double*)descA->get_blkaddr(descA,k,k);
    starpu_insert_task(cl,
            STARPU_VALUE, &alpha,  sizeof(double),
            STARPU_VALUE, &M, sizeof(int),
            STARPU_VALUE, &data,   sizeof(double*),
            STARPU_VALUE, &lda,   sizeof(int),
            0);
}


void cl_dscaldiag_cpu_func(void *buffers[],void *cl_arg){
    int M, LDA; 
    double *A, alpha;
      
    starpu_codelet_unpack_args(cl_arg, &alpha, &M, &A, &LDA);
    cblas_dscal( M, alpha, A, LDA+1 );
}

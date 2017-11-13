/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef CODELET_DSCALDIAG_H
#define CODELET_DSCALDIAG_H

#include <morse.h>

#ifdef __cplusplus
extern "C"{
#endif
int APPNAME_TASK_dscaldiag(int M, double alpha, MORSE_desc_t *descA, int k, int lda);
void cl_dscaldiag_cpu_func(void *buffers[],void *cl_arg);
#ifdef __cplusplus
}
#endif

#endif //CODELET_DSCALDIAG_H

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

namespace MOAO_STARPU
{
template<typename T>
int TASK_Xscaldiag(int M, T alpha, MORSE_desc_t *descA, int k, int lda);
}
template<typename T>
void cl_scaldiag_cpu_func(void *buffers[],void *cl_arg);
#endif //CODELET_DSCALDIAG_H

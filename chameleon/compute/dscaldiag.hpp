/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef DSCALDIAG_H
#define DSCALDIAG_H

#include <morse.h>
#include "codelet_dscaldiag.hpp"

namespace MOAO_CHAMELEON
{
template<typename T>
int Xscaldiag(T alpha, MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t  *request) ;
}

#endif //DSCALDIAG_H

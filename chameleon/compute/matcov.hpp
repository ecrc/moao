/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef MATCOV_H
#define MATCOV_H

#include "common.hpp"

#include <starpu.h>
#include <morse.h>

namespace MOAO_CHAMELEON
{
template<typename T>
int matcov_tile( MORSE_desc_t *descData, MORSE_sequence_t *sequence, MORSE_request_t *request, struct Tomo_struct<T> *tomo, int type_mat);
}
#endif //MATCOV_H

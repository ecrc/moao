/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef CEE_CVV_H
#define CEE_CVV_H

#include <morse.h>


namespace MOAO_CHAMELEON
{
template<typename T>
void Cee_Cvv(MORSE_desc_t *Cmm, MORSE_desc_t *Cpp, MORSE_desc_t *Cpm, MORSE_desc_t *R, MORSE_desc_t *Dx,
                   MORSE_desc_t *Cee, MORSE_desc_t *Cvv, MORSE_desc_t *Tmp, MORSE_sequence_t *sequence, int sync=0);
}

#endif // CEE_CVV_H

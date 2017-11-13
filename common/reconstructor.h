/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef RECONSTRUCTOR_H
#define RECONSTRUCTOR_H

#include <morse.h>

#ifdef __cplusplus
extern "C"{
#endif
void APPNAME_reconstructor(MORSE_desc_t *Cmm, MORSE_desc_t *Ctm, MORSE_sequence_t *sequence,int check);
#ifdef __cplusplus
}
#endif
#endif //RECONSTRUCTOR_H

/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef CODELET_MATCOV_H
#define CODELET_MATCOV_H

#include <starpu.h>
#include <morse.h>
#include "utils.h"
#include "tomo_struct.h"

#define APPNAME_desc MORSE_desc_t
#define RTBLKADDR( desc, m, n ) ( (starpu_data_handle_t)RUNTIME_desc_getaddr( desc, m, n ) )

#ifdef __cplusplus
extern "C"{
#endif
void APPNAME_TASK_matcov(APPNAME_desc *data,int nrows,int ncols, int m, int n,
                    int lda, struct tomo_struct tomo,
MORSE_desc_t *descSspSizeL,
MORSE_desc_t *descNssp,
MORSE_desc_t *descU,
MORSE_desc_t *descV,
MORSE_desc_t *descX,
MORSE_desc_t *descY,
MORSE_desc_t *descL0diff,
MORSE_desc_t *descCn2,
MORSE_desc_t *descNsubap,
MORSE_desc_t *descNoiseNGS,
MORSE_desc_t *descNoiseLGSxx,
MORSE_desc_t *descNoiseLGSyy,
MORSE_desc_t *descNoiseLGSxy,
                    int type_mat);
#ifdef __cplusplus
}
#endif

#endif //CODELET_MATCOV_H

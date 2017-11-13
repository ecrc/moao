/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef NOISE_H
#define NOISE_H

#include "moao_defs.h"
#include "tomo_struct.h"
#include<stdio.h>

#ifdef __cplusplus
extern "C"{
#endif

real_t flux(real_t qr, real_t mr, real_t throughput, real_t d, real_t t, real_t bdw);
real_t magnitudeFromFlux(real_t qr, real_t F, real_t throughput, real_t d, real_t t, real_t bdw);
void updateSeeing(struct tomo_struct *tomo);
void updateNoise(struct tomo_struct *tomo);
void generateElongation(real_t *x, real_t *y, real_t *alphaX, real_t *alphaY,int nlgs, real_t diam, real_t *lgsExt, real_t *lgsTheta, real_t dH_lgs, real_t alt_lgs, long *Nsubap );

#ifdef __cplusplus
}
#endif

#endif //NOISE_H

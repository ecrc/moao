/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef NOISE_H
#define NOISE_H

#define RASC 206265.

#include "tomo_struct.hpp"
#include <vector>
#include <stdio.h>


template<typename T>
T flux(T qr, T mr, T throughput, T d, T t, T bdw);
template<typename T>
T magnitudeFromFlux(T qr, T F, T throughput, T d, T t, T bdw);
//template<typename T>
//void generateElongation(std::vector<T> &x, std::vector<T> &y, std::vector<T> &alphaX, std::vector<T> &alphaY,int nlgs, T diam, std::vector<T> &lgsExt, std::vector<T> &lgsTheta, T dH_lgs, T alt_lgs, std::vector<long> &Nsubap );

//#ifdef __cplusplus
//}
//#endif

#include "noise.cpp"

#endif //NOISE_H

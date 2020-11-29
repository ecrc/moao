/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef  NOISE_IMPL_
#define  NOISE_IMPL_

#include "noise.hpp"
#include <math.h>

template<typename T>
T flux(T qr, T mr, T throughput, T d, T t, T bdw){
/** compute the flux of e- per subapertur and per (time) frame
*
* qr            : Photometric Flux offset
* mr            : R stars magnitude
* throughput    : transmission through atmospheric layers
* d             : subaperture diameter (m)
* t             : frame rate (Hz)
* bdw		: bandwidth in A
*
* return:
* flux: /ph/subap/frame
**/
    //flux: /ph/subap/frame
    return pow(10,qr-0.4*mr)*throughput*d*d*t*bdw*1.e4;
}
template float flux<float>(float qr, float mr, float throughput, float d, float t, float bdw);
template double flux<double>(double qr, double mr, double throughput, double d, double t, double bdw);

template<typename T>
T magnitudeFromFlux(T qr, T F, T throughput, T d, T t, T bdw){
/** compute the magnitude from the flux
*
* qr            : Photometric Flux offset
* flux: /ph/subap/frame
* throughput    : transmission through atmospheric layers
* d             : subaperture diameter (m)
* t             : frame rate (Hz)
* bdw		: bandwidth in A
*
* return:
* mr            : R stars magnitude
**/

    return -(log10(F/(throughput*d*d*bdw*1.e4*t))-qr)/0.4;

}
template float magnitudeFromFlux<float>(float qr, float F, float throughput, float d, float t, float bdw);
template double magnitudeFromFlux<double>(double qr, double F, double throughput, double d, double t, double bdw);



/**
*generate an array of elongation for the LGS
*
* x        : x coordinates of the subapertures
* y        : y coordinates of the subapertures
* alphaX   : direction(in x) of the lgs
* alphaY   : direction(in y) of the lgs
* lgsExt   : length of the lgs extension
* lgsTheta : angle  of the lgs extension
* dH_lgs 
* alt_lgs 
*
**/
template<typename T>
void generateElongation(std::vector<T> &x, std::vector<T> &y, std::vector<T> &alphaX, std::vector<T> &alphaY,int nlgs,
            T diam, std::vector<T> &lgsExt, std::vector<T> &lgsTheta, T dH_lgs, T alt_lgs, std::vector<long> &Nsubap ){

    int w,i,k=0;
    const T toFwhm=RASC*dH_lgs/(alt_lgs*alt_lgs);
    for(w=0;w<nlgs;w++){
        const T lltx=alphaX[w]*(diam/2./sqrt(alphaX[w]*alphaX[w]+alphaY[w]*alphaY[w]));
        const T llty=alphaY[w]*(diam/2./sqrt(alphaX[w]*alphaX[w]+alphaY[w]*alphaY[w]));
        for(i=0;i<Nsubap[w];i++){
            lgsExt[k]=sqrt((pow(x[k]-lltx,2)+pow(y[k]-llty,2)))*toFwhm;
            lgsTheta[k]=atan2((y[k]-llty)*toFwhm,(x[k]-lltx)*toFwhm);
            k++;
        }
    }
}


#endif //NOISE_IMPL_

/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "noise.h"
#include <math.h>

real_t flux(real_t qr, real_t mr, real_t throughput, real_t d, real_t t, real_t bdw){
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

real_t magnitudeFromFlux(real_t qr, real_t F, real_t throughput, real_t d, real_t t, real_t bdw){
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


void updateSeeing(struct tomo_struct *tomo){
    int i;
    const real_t lambda_m=500.*1.e-9; //(meters)
    real_t r0_500nm=0;
    for(i=0;i<tomo->Nlayer;i++){
        r0_500nm+=tomo->cn2[i];
    }
    //r0 in meters
    r0_500nm=pow(r0_500nm,-3./5.);

    // r0 = (lambda/lambda_500nm)^(6/5)*r0_500nm
    // s = lambda/r
    // s = lambda_500nm^(6/5) / (lambda^(1/5)*r0_500nm)
    // CONVERT IN ARCSEC
    tomo->sNGS2=lambda_m/r0_500nm*pow(lambda_m/tomo->lambdaNGS,0.2)*RASC;
    //get square of seeing
    tomo->sNGS2*=tomo->sNGS2;

    //same for lgs
    tomo->sLGS2=lambda_m/r0_500nm*pow(lambda_m/tomo->lambdaLGS,0.2)*RASC;
    tomo->sLGS2*=tomo->sLGS2;

}

void updateNoise(struct tomo_struct *tomo){
/** update the noise arrays of the tomo struct
*
**/
    //update square of seeing
    updateSeeing(tomo);

    int i,j,k=0;
    for(i=0;i<tomo->Nw;i++){
        const real_t d=(tomo->DiamTel/tomo->Nssp[i]); //already taken into account in F ?
        const real_t F=flux(tomo->qr, tomo->mr[i],tomo->throughput[i], d, tomo->Tframe,tomo->bdw);

        const real_t fact1=0.32/F;
        const real_t fact2=0.82*pow(tomo->RON/(F*tomo->pixSize[i]),2);

        real_t noiseLongAxis;
        if(i<tomo->nlgs){
            const real_t partialExt=tomo->sLGS2 + (tomo->spot_width*tomo->spot_width);

            const real_t noiseShortAxis=(fact1*partialExt+
                                         fact2* pow((partialExt),1.5)* pow(partialExt,0.5));

            for(j=0;j<tomo->Nsubap[i];j++){

                noiseLongAxis=(fact1*(partialExt+tomo->lgsExt[k]*tomo->lgsExt[k])+
                              fact2* pow((partialExt+tomo->lgsExt[k]*tomo->lgsExt[k]),1.5)* pow(partialExt,0.5));
                tomo->noiseLGSxx[k]=noiseLongAxis* pow(cos(tomo->lgsTheta[k]),2) + noiseShortAxis*pow(sin(tomo->lgsTheta[k]),2);
                tomo->noiseLGSyy[k]=noiseLongAxis* pow(sin(tomo->lgsTheta[k]),2) + noiseShortAxis*pow(cos(tomo->lgsTheta[k]),2);
                tomo->noiseLGSxy[k]=(noiseLongAxis-noiseShortAxis)*sin(tomo->lgsTheta[k])*cos(tomo->lgsTheta[k]);

		        k++;
            }
        }
        else{
	        tomo->noiseNGS[i-tomo->nlgs]=(fact1*tomo->sNGS2+fact2*tomo->sNGS2*tomo->sNGS2);
        }
    }
}

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
void generateElongation(real_t *x, real_t *y, real_t *alphaX, real_t *alphaY,int nlgs, real_t diam, real_t *lgsExt, real_t *lgsTheta, real_t dH_lgs, real_t alt_lgs, long *Nsubap ){

    int w,i,k=0;
    const real_t toFwhm=RASC*dH_lgs/(alt_lgs*alt_lgs);
    for(w=0;w<nlgs;w++){
        const real_t lltx=alphaX[w]*(diam/2./sqrt(alphaX[w]*alphaX[w]+alphaY[w]*alphaY[w]));
        const real_t llty=alphaY[w]*(diam/2./sqrt(alphaX[w]*alphaX[w]+alphaY[w]*alphaY[w]));
        for(i=0;i<Nsubap[w];i++){
            lgsExt[k]=sqrt((pow(x[k]-lltx,2)+pow(y[k]-llty,2)))*toFwhm;
            lgsTheta[k]=atan2((y[k]-llty)*toFwhm,(x[k]-lltx)*toFwhm);
            k++;
        }
    }
}

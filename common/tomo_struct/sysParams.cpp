/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/
#ifndef SYSPARAMS_IMPL_
#define SYSPARAMS_IMPL_


#include "sysParams.hpp"

#include <stdexcept>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "utils/utils_vector.hpp"
#include "utils/utils_IO.hpp"


/*!\brief Default constructor
 *
 */
template<typename T>
SysParams<T>::SysParams(){
}

/*!\brief Constructor
 *
 * Initialize parameters from file
 * @param[in] sf : string : parameter file to read
 */
template<typename T>
SysParams<T>::SysParams(std::string sf){
    read(sf);
}

/*!\brief Constructor
 *
 * create a SysParam from a limited number of parameters
 * other parameter are default ones
 *
 * @param[in] d      : meter    : telescope diameter
 * @param[in] nNGS   :          : number of NGS
 * @param[in] nLGS   :          : number of LGS
 * @param[in] NGSmag :          : magnitude of NGS (common for all)
 * @param[in] Flux   :(ph/m2/s) :  LGS photon return at M1
 * @param[in] seed   :          : seed for the NGS position (random)
 */
template<typename T>
SysParams<T>::SysParams(T d, int nNGS, int nLGS, T NGSmag, T Flux, int seed){
    diam        =d;
    obs         =0.2;
    tFrame      =0.004;
    nW          =nNGS+nLGS+1;
    nNgs        =nNGS;
    nLgs        =nLGS;
    nTarget     =1;
    nssp        =diam*2.;
    gsAlt       =std::vector<T>(nW);
    gsType      =std::vector<int>(nW);
    alphaX_as   =std::vector<T>(nW);
    alphaY_as   =std::vector<T>(nW);
    XPup        =std::vector<T>(nW);
    YPup        =std::vector<T>(nW);
    thetaML     =std::vector<T>(nW);
    thetaCam    =std::vector<T>(nW);
    sensibilite =std::vector<T>(nW,1.);
    tracking    =std::vector<T>(3,1.);
    pasDphi     =0.0001;
    ncpu        =1;
    mrNGS       =std::vector<T>(nNGS+1);
    lgsFlux     =Flux;
    ngsPixSize  =0.3;
    lgsPixSize  =1.;
    lambdaNGS   =6.5e-07;
    lambdaLGS   =5.89e-07;
    bdw_m       =3.3e-07;
    throughNGS  =0.425;
    throughLGS  =0.382;
    throughAtm  =0.84;
    RON         =3;
    lgsCst      =0.1;
    spotWidth   =1;
    lgsAlt      =100000;
    lgsDepth    =50000;
    targetX_as  =std::vector<T>(1);
    targetY_as  =std::vector<T>(1);

    std::fill(gsAlt.begin(), gsAlt.begin()+nLGS,1.e-5);
    std::fill(gsType.begin(),gsType.begin()+nLGS,2);    std::fill(gsType.begin()+nLGS,gsType.end(),1); 
    // launch angle diameter for lgs (arcmin)
    T radius=7.4; //for EELT
    // convert to arcsec
    radius=radius/2.*60;
    //scale to telescope diameter
    radius*=(diam/38.542);
    fillVectorCircle(alphaX_as ,alphaY_as ,0,nLGS,radius);
    fillVectorRand(alphaX_as ,nLGS,nNGS,radius/(T)2.,radius);
    fillVectorRand(alphaY_as ,nLGS,nNGS,radius/(T)2.,radius);

}

/*!\bief Check the size of the SysParam vectors
 *
 * each vector of AtmParam must have at least nLayer elements
 * throw an exception if this is not the case
 * return a warning if the vector ahas more elements than the number of layers
 */
template<typename T>
void SysParams<T>::vectorSize(){
    int err=0;
    err=std::max(err,checkVectorSize<T  >(gsAlt       ,nW         ,"gsAlt"        ));
    err=std::max(err,checkVectorSize<int>(gsType      ,nW         ,"type"         ));
    err=std::max(err,checkVectorSize<T  >(alphaX_as   ,nW         ,"alphaX_as "   ));
    err=std::max(err,checkVectorSize<T  >(alphaY_as   ,nW         ,"alphaY_as "   ));
    err=std::max(err,checkVectorSize<T  >(XPup        ,nW         ,"XPup"         ));
    err=std::max(err,checkVectorSize<T  >(YPup        ,nW         ,"YPup"         ));
    err=std::max(err,checkVectorSize<T  >(thetaML     ,nW         ,"thetaML"      ));
    err=std::max(err,checkVectorSize<T  >(thetaCam    ,nW         ,"thetaCam"     ));
    err=std::max(err,checkVectorSize<T  >(sensibilite ,nW         ,"sensibilite"  ));
    err=std::max(err,checkVectorSize<T  >(tracking    ,3          ,"tracking"     ));
    err=std::max(err,checkVectorSize<T  >(mrNGS       ,nW-nLgs    ,"mrNGS"        ));
    err=std::max(err,checkVectorSize<T  >(targetX_as  ,nTarget    ,"targetX_as"   ));
    err=std::max(err,checkVectorSize<T  >(targetY_as  ,nTarget    ,"targetY_as"   ));
    if(err>0)
        std::cerr<<std::endl;
    if(err>1)
        throw std::runtime_error("System Parameter vector size does not matches");
}

/*!\brief Read parameters to file
 *
 * each entry of the parameter file is preceed by a comment line
 *
 *@param[in] sf             : string    : input file name 
 *@param[in] checkComment   : bool      : check comment line preceeding data entry
 */
template<typename T>
void SysParams<T>::read(std::string sf, bool checkComment){
    std::cout<<"    Reading system parameters: "<<sf<<std::endl;
    std::string err("");
    readIO(sf, *this);
    nNgs=nW-nLgs-1;
    
    if(!err.empty()){
        std::cerr<<"WARNING:\n"<< err<<std::endl;;
    }
    vectorSize();
}

template<typename T>
void SysParams<T>::write(std::string sf){
    std::cout<<"Writing system parameters: "<<sf<<std::endl;
    writeIO(sf,*this);
}

/*!\brief return attribute list
 *
 * @param[inout] e : string : string on which errors will be appened
 * @return : list<IO> : list of IO with all SysParams attributes
 */
template<typename T>
std::list<IO> SysParams<T>::getParams(std::string & e)
{
    std::list<IO> l;
    appendNamedIO(l,diam       , e, "diam       : meter     : Telescope diameter    "                               );
    appendNamedIO(l,obs        , e, "obs        : percent   : Central obscuration"                                  );
    appendNamedIO(l,tFrame     , e, "tFrame     : second    : frame rate"                                           );
    appendNamedIO(l,nW         , e, "nW         :           : number of WFS"                                        );
    appendNamedIO(l,nLgs       , e, "nLgs       :           : number of LGS"                                        );
    appendNamedIO(l,nTarget    , e, "nTarget    :           : number of Target"                                     );
    appendNamedIO(l,nssp       , e, "Nssp       :           : number of subaperture per wfs along the diameter"     );
    appendNamedIO(l,gsAlt      , e, "gsAlt      : meter^-1  : inverse of lazer altitude"                            );
    appendNamedIO(l,gsType     , e, "type       :           : guide star type (1:NGS, 2:LGS)"                       );
    appendNamedIO(l,alphaX_as  , e, "alphaX_as  : arcsec    : pointing direction of the wfs on x axis"              );
    appendNamedIO(l,alphaY_as  , e, "alphaY_as  : arcsec    : pointing direction of the wfs on y axis"              );
    appendNamedIO(l,XPup       , e, "XPup       : meter     : pupil shift of the WFS"                               );
    appendNamedIO(l,YPup       , e, "YPup       : meter     : pupil shift of the WFS"                               );
    appendNamedIO(l,thetaML    , e, "thetaML    :           : rotation of the microlenses"                          );
    appendNamedIO(l,thetaCam   , e, "thetaCam   :           : rotation of the camera"                               );
    appendNamedIO(l,sensibilite, e, "sensibility:           : sensitivity coeff of this WFS"                        );
    appendNamedIO(l,tracking   , e, "tracking   : arcsec^2  : telescope tracking error parameters (x^2, y^2 and xy)");
    appendNamedIO(l,pasDphi    , e, "pasDPHI    :           : Precision of DPHI precomputation. //deprecated"       );
    appendNamedIO(l,ncpu       , e, "ncpu       :           : Number of CPU used (only with openMP)"                );
    appendNamedIO(l,mrNGS      , e, "mrNGS      :           : magnitude of NGS"                                     );
    appendNamedIO(l,lgsFlux    , e, "lgsFlux    : (ph/m2/s) : LGS photon return at M1"                              );
    appendNamedIO(l,ngsPixSize , e, "ngsPixSize : arcsec    : NGS pixel size"                                       );
    appendNamedIO(l,lgsPixSize , e, "lgsPixSize : arcsec    : LGS pixel size"                                       );
    appendNamedIO(l,lambdaNGS  , e, "lambdaNGS  : meter     : wave length for NGS"                                  );
    appendNamedIO(l,lambdaLGS  , e, "lambdaLGS  : meter     : wave length for LGS"                                  );
    appendNamedIO(l,bdw_m      , e, "bdw_m      : meter     : bandwidth"                                            );
    appendNamedIO(l,throughNGS , e, "throughNGS : percent   : transmission for NGS"                                 );
    appendNamedIO(l,throughLGS , e, "throughLGS : percent   : transmission for LGS"                                 );
    appendNamedIO(l,throughAtm , e, "throughAtm : percent   : atmosphere transmission"                              );
    appendNamedIO(l,RON        , e, "RON        : nb of e-  : Read Out Noise "                                      );
    appendNamedIO(l,lgsCst     , e, "lgsCst     :           : constant on lgs (simulate that LGS cannot measure tip-tilt and focus)");
    appendNamedIO(l,spotWidth  , e, "spotWidth  : arcsec    : lazer width"                                          );
    appendNamedIO(l,lgsAlt     , e, "lgsAlt     : meter     : sodium layer altitude"                                );
    appendNamedIO(l,lgsDepth   , e, "lgsDepth   : meter     : depth of the sodium layer"                            );
    appendNamedIO(l,targetX_as , e, "targetX_as : arcsec    :  taget direction on x axis"                           );
    appendNamedIO(l,targetY_as , e, "targetY_as : arcsec    :  taget direction on y axis"                           );
    return l;
   
}

#endif // SYSPARAMS_IMPL_

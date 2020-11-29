#ifndef TOMO_STRUCT_H
#define TOMO_STRUCT_H

#include <string>
#include "sysParams.hpp"
#include "atmParams.hpp"

template<typename T>
class Tomo_struct  {

    public:
    std::string sys_path;       //!< telescope parameter file
    std::string atm_path;       //!< atmosphere parameter file

    SysParams<T> sys;
    AtmParams<T> atm;

    std::vector<T> X;           //!< subapertures coordinates (in X) for all subapertures and all WFS
    std::vector<T> Y;           //!< subapertures coordinates (in Y) for all subapertures and all WFS
    std::vector<long> Nssp;      //!<           : number of subaperture for all wfs along the diameter
    std::vector<long> Nsubap;   //!< array of the number of subap of each WFS, contains Nw elements 
    std::vector<T> diamPup;     //!< WFS magnification factor
    std::vector<T> sspSize;     //!< WFS subaperture size
    T rmax;                     //!< maximum distance between subapertures

    T qr;                       //!<           : photometric Flux offset
    T bdw;                      //!< angstrom  : bandwidth
    std::vector<T>  mr;         //!<           : guide stars magnitude
     std::vector<T> throughput; //!<           :
    std::vector<T>  alphaX    ; //!< radian    : pointing direction of the wfs on x axis
    std::vector<T>  alphaY    ; //!< radian    : pointing direction of the wfs on y axis
    std::vector<T>  targetX   ; //!< radian    : taget direction on x axis 
    std::vector<T>  targetY   ; //!< radian    : taget direction on y axis 

    int Nx;                     //!<           : total number of valid subapertures
    int Nslopes;                //!< unused
    int part;                   // Computed part of the cov. mat. 0: complete 1: cmm 2: cpp 3: cpm ??
    T sNGS2;                    //!< square of the seeing at NGS lambda
    T sLGS2;                    //!< square of the seeing at LGS lambda
    std::vector<T> pixSize;     //!< arcsec    : pixel size for each WFS
    std::vector<T> lgsExt;      //!< extension of lgs
    std::vector<T> lgsTheta;    //!< angle of lgs extension
    std::vector<T> noiseNGS;    //!< WFS noise of the NGS 
    std::vector<T> noiseLGSxx;  //!< WFS noise of LGS for the xx covariances
    std::vector<T> noiseLGSyy;  //!< WFS noise of LGS for the yy covariances
    std::vector<T> noiseLGSxy;  //!< WFS noise of LGS for the xy covariances

    //Ali
    std::vector<long> indexL0;
    std::vector<T> L0diff;
    std::vector<T> u;           //!< subapertures coordinates (in X) for all projected subapertures and all WFS
    std::vector<T> v;           //!< subapertures coordinates (in Y) for all projected subapertures and all WFS
    std::vector<T> sspSizeL;
    long    nsubaps_offaxis;

    public:
    Tomo_struct();
    Tomo_struct(std::string sf, std::string af);
    void write_sys(std::string sf);
    void write_atm(std::string af);
    void write(std::string sf, std::string af);

    void read_sys(std::string sf="");
    void read_atm(std::string af="");
    void read(std::string sf="", std::string af="");

    void setSysParams();
    void getSysParams();
    void initSys();

    void setAtmParams();
    void getAtmParams();
    void initAtm();
    void updateAtm(std::string af="");
    void allocAtm();

    void setTarget(int tarInd);
    void setTargetCoord(T aX, T aY, int toRad=0);

    //TODO updateCn2(std::vector<T> newCn2);
    //TODO updateH  (std::vector<T> newH  );
    //TODO updateL0 (std::vector<T> newL0 );

    void updateSeeing();
    void updateNoise();

    void validSubap();
    void generateElongation();
    void projectedSubap();

    long getNMeas();
    long getNMeasTS();

};
template class Tomo_struct<float>;
template class Tomo_struct<double>;

#endif // TOMO_STRUCT_H

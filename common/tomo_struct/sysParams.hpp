#ifndef SYSPARAMS_H
#define SYSPARAMS_H
#include <vector>
#include <list>
#include <string>
#include "io/IO.hpp"

template<typename T>
class SysParams{
    public:
    T               diam;       //!< meter     : Telescope diameter    
    T               obs;        //!< percent   : Central obscuration
    T               tFrame;     //!< second    : frame rate
    int             nW;         //!<           : number of WFS
    int             nLgs;       //!<           : number of LGS
    int             nNgs;       //!<           : number of NGS
    int             nTarget;    //!<           : number of Target
    int             nssp;       //!<           : number of subaperture per wfs along the diameter
    std::vector<T>  gsAlt;      //!< meter^-1  : inverse of lazer altitude // replaced with 1/lgsAlt
    std::vector<int>gsType;     //!<           : guide star type (1:NGS, 2:LGS , 3= (not available yet)TipTilt-guide star)
    std::vector<T>  alphaX_as ; //!< arcsec    : pointing direction of the wfs on x axis
    std::vector<T>  alphaY_as ; //!< arcsec    : pointing direction of the wfs on y axis
    std::vector<T>  XPup;       //!< meter     : pupil shift of the WFS (on axis X) 
    std::vector<T>  YPup;       //!< meter     : pupil shift of the WFS (on axis Y)
    std::vector<T>  thetaML;    //!<           : rotation of the microlenses
    std::vector<T>  thetaCam;   //!<           : rotation of the camera
    std::vector<T>  sensibilite;//!<           : sensitivity coeff of this WFS
    std::vector<T>  tracking;   //!< arcsec^2  : telescope tracking error parameters (x^2, y^2 and xy)
    T               pasDphi;    //!<           : Precision of DPHI precomputation. //deprecated
    int             ncpu;       //!<           : Number of CPU used (only with openMP)
    std::vector<T>  mrNGS;      //!<           : magnitude of NGS
    T               lgsFlux;    //!< (ph/m2/s) : LGS photon return at M1
    T               ngsPixSize; //!< arcsec    : NGS pixel size
    T               lgsPixSize; //!< arcsec    : LGS pixel size
    T               lambdaNGS;  //!< meter     : wave length for NGS
    T               lambdaLGS;  //!< meter     : wave length for LGS
    T               bdw_m;      //!< meter     : bandwidth 
    T               throughNGS; //!< percent   : transmission for NGS
    T               throughLGS; //!< percent   : transmission for LGS
    T               throughAtm; //!< percent   : atmosphere transmission
    int             RON;        //!< nb of e-  : Read Out Noise 
    T               lgsCst;     //!<           : constant on lgs (simulate that LGS cannot measure tip-tilt and focus)
    T               spotWidth;  //!< arcsec    : lazer width
    T               lgsAlt;     //!< meter     : sodium layer altitude
    T               lgsDepth;   //!< meter     : depth of the sodium layer
    std::vector<T>  targetX_as; //!< arcsec    :  taget direction on x axis 
    std::vector<T>  targetY_as; //!< arcsec    :  taget direction on y axis 


    public:
        SysParams();
        SysParams(std::string sf);
        SysParams(T d, int nNGS, int nLGS, T NGSmag=13., T lgsFlux=1.e7, int seed=1234);
        void vectorSize();
        void read(std::string sf, bool checkComment=false);
        void write(std::string sf);

        std::list<IO> getParams(std::string & e);
};
template class SysParams<float>;
template class SysParams<double>;

#endif //SYSPARAMS_H

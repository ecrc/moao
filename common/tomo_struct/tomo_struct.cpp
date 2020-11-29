#include "tomo_struct.hpp"
#include "noise.hpp"

#include <algorithm>
#include <numeric>

template<typename T>
Tomo_struct<T>::Tomo_struct(){};

template<typename T>
Tomo_struct<T>::Tomo_struct(std::string sf, std::string af){
    sys_path=sf;
    initSys();

    atm_path=af;
    initAtm();
}

/*! write optical system parameters to file
 * 
 *@param[in] sf : string : output file name
 * */
template<typename T>
void Tomo_struct<T>::write_sys(std::string sf){
    sys.write(sf);
    setSysParams();
}

/*! write atmospheric parameters to file
 * 
 *@param[in] af : string : output file name
 * */
template<typename T>
void Tomo_struct<T>::write_atm(std::string af){
    atm.write(af);
}

/*! write optical system and atmospheric parameters to file
 * 
 *@param[in] sf : string : output file name for the telescope 
 *@param[in] af : string : output file name for the atmosphere
 * */
template<typename T>
void Tomo_struct<T>::write(std::string sf,std::string af){
    write_sys(sf);
    write_atm(af);
}

/*! read optical system parameters from file
 * 
 *@param[in] sf : string : input file name
 * */
template<typename T>
void Tomo_struct<T>::read_sys(std::string sf){
    if(sf!=""){
        sys_path=sf;
    }
    sys.read(sys_path);
    getSysParams();
}

/*! read optical atmospheric parameters from file
 * 
 *@param[in] af : string : input file name
 * */
template<typename T>
void Tomo_struct<T>::read_atm(std::string af){
    if(af!=""){
        atm_path=af;
    }
    atm.read(atm_path);
}

/*! read optical system and atmospheric parameters from file
 * 
 *@param[in] sf : string : input file name for the telescope 
 *@param[in] af : string : input file name for the atmosphere
 * */
template<typename T>
void Tomo_struct<T>::read(std::string sf, std::string af){
    read_sys(sf);
    read_atm(af);
}

/*!\brief Set up the optical systems attribute
 *
 * Based on the data from the SysParam 
 * */
template<typename T>
void Tomo_struct<T>::getSysParams(){
    qr=2.9; //TODO as param?
    bdw=sys.bdw_m*1.e10;
    
    Nssp = std::vector<long>(sys.nW, sys.nssp);

    pixSize = std::vector<T>(sys.nW);
    std::fill(pixSize.begin(),pixSize.begin()+sys.nLgs,sys.lgsPixSize);
    std::fill(pixSize.begin()+sys.nLgs,pixSize.end(),sys.ngsPixSize);

    throughput = std::vector<T>(sys.nW);
    std::fill(throughput.begin(),throughput.begin()+sys.nLgs,sys.throughLGS);
    std::fill(throughput.begin()+sys.nLgs,throughput.end(),sys.throughNGS);

    diamPup=std::vector<T>(Nssp.begin(),Nssp.end());

    //WFS pointing direction to radian
    alphaX=std::vector<T>(sys.nW);
    std::transform(sys.alphaX_as.begin(), sys.alphaX_as.end(), 
                    alphaX.begin(),std::bind1st(std::multiplies<T>(),1/RASC));
    alphaY=std::vector<T>(sys.nW);
    std::transform(sys.alphaY_as.begin(), sys.alphaY_as.end(),
                    alphaY.begin(),std::bind1st(std::multiplies<T>(),1/RASC));

    //target pointing direction to radian
    targetX=std::vector<T>(sys.nW);
    std::transform(sys.targetX_as.begin(), sys.targetX_as.end(),
                    targetX.begin(),std::bind1st(std::multiplies<T>(),1/RASC));
    targetY=std::vector<T>(sys.nW);
    std::transform(sys.targetY_as.begin(), sys.targetY_as.end(),
                    targetY.begin(),std::bind1st(std::multiplies<T>(),1/RASC));

    //transform: divide diam by Nssp (for all values)
    sspSize=std::vector<T>(sys.nW);
    std::transform(Nssp.begin(),Nssp.end(),sspSize.begin(),[&](int n){return sys.diam/n;} );

    mr=std::vector<T>(sys.nW);
    std::copy(sys.mrNGS.begin(),sys.mrNGS.end(),mr.begin()+sys.nLgs);
    //transfor: compute corresponding magnitude for LGS 
    T fluxPerSubapPerFrame = sys.lgsFlux*sys.tFrame*pow(sys.diam/Nssp[0],2.)*throughput[0]/sys.throughAtm;  //e-/subs/frame
    std::transform(throughput.begin(),std::next(throughput.begin(),sys.nLgs), sspSize.begin(), mr.begin(),
            [&] (T through, T ss){return magnitudeFromFlux(qr,fluxPerSubapPerFrame,through,ss, sys.tFrame, bdw);} ); 

}

/*!\brief Set up the optical systems parameter
 *
 * Based on the data, convert/extract datas describing the system
 * */
template<typename T>
void Tomo_struct<T>::setSysParams(){
    sys.bdw_m=bdw/1.e10;
    
    sys.nssp=Nssp[0];

    pixSize = std::vector<T>(sys.nW);
    if(sys.nLgs>0)
    {
        sys.lgsPixSize=pixSize[0];
        sys.throughLGS=throughput[0];
        sys.lgsFlux=flux(qr,mr[0],throughput[0],sspSize[0], sys.tFrame, bdw);
    }
    sys.ngsPixSize=pixSize[sys.nLgs];
    sys.throughNGS=throughput[sys.nLgs];

    diamPup=std::vector<T>(Nssp.begin(),Nssp.end());


    //WFS pointing direction to arcsec
    std::transform(alphaX.begin(), alphaX.end(), sys.alphaX_as.begin(),std::bind1st(std::multiplies<T>(),RASC));
    std::transform(alphaY.begin(), alphaY.end(), sys.alphaY_as.begin(),std::bind1st(std::multiplies<T>(),RASC));

    //target pointing direction to arcsec
    std::transform(targetX.begin(), targetX.end(), sys.targetX_as.begin(),std::bind1st(std::multiplies<T>(),RASC));
    std::transform(targetY.begin(), targetY.end(), sys.targetY_as.begin(),std::bind1st(std::multiplies<T>(),RASC));


}

template<typename T>
void Tomo_struct<T>::getAtmParams(){
    
}

template<typename T>
void Tomo_struct<T>::setAtmParams(){
}


/*!\brief Initialize optical system variables
 *
 * define the valid subapertures
 * compute the laser elongation
 *
 * */
template<typename T>
void Tomo_struct<T>::initSys(){
    std::cout<<"Initialize system\n"<<std::flush;
    sys=SysParams<T>(sys_path);
    getSysParams();
    validSubap();
    generateElongation();
    std::cout<<"done"<<std::endl;
    
}

/*!\brief Initialize atmospheric variables
 *
 * */
template<typename T>
void Tomo_struct<T>::initAtm(){
    atm=AtmParams<T>();
    updateAtm(atm_path);
}

/*!\brief Update atmospheric variables
 *
 * in particular cn2,h,L0.
 * update the WFS noise and subapertures projection
 *
 *@params[in] af : string : new parameter file
 * */
template<typename T>
void Tomo_struct<T>::updateAtm(std::string af){
    std::cout<<"Updating atmospheric parameters...\n"<<std::flush;
    if(!af.empty()){
        int nLayer_old=atm.nLayer;
        atm_path=af;
        atm.read(atm_path);
        if(nLayer_old!=atm.nLayer){
            allocAtm();
        }
    }

    int cc;
    for (cc = 0; cc < sys.nW * atm.nLayer; cc++) {
      int n = cc / atm.nLayer;
      int l = cc - n * atm.nLayer;
      if(n >= sys.nW) n-=1;
      sspSizeL[cc] = sspSize[n] * (1. - sys.gsAlt[n] * atm.h[l]);
    }
    //Search the different L0 and build indexL0
    const long cNlayer = atm.nLayer;
    long i, j;
    int cpt = 1;
    T tmp[cNlayer];

    tmp[0] = atm.L0[0];
    indexL0[0] = 0;
    
    for (i = 1; i < cNlayer; i++) {
      j = 0;
      const T l0 = atm.L0[i];
      
      while ((j < cpt) && (tmp[j] != l0)){
          j++;
      }
      indexL0[i] = j;      
      if (j == cpt) {
        tmp[j] = l0;
        cpt++;
      }
    }
    for (i = 0; i < cNlayer; i++)  {
      L0diff[i] = atm.L0[i];
    }
    //Computes  u and v
    projectedSubap();
    //Compute noise
    updateNoise();

    std::cout<<"done"<<std::endl;
}

template<typename T>
void Tomo_struct<T>::allocAtm(){
    indexL0  =std::vector<long>(atm.nLayer);
    L0diff   = std::vector<T>(atm.nLayer*Nx);
    u        = std::vector<T>(atm.nLayer*Nx);
    v        = std::vector<T>(atm.nLayer*Nx);
    sspSizeL = std::vector<T>(sys.nW*atm.nLayer);
}

/*!\brief change target using target list
 *
 *@param[in] tarInd : int : target index
 * */
template<typename T>
void Tomo_struct<T>::setTarget(int tarInd){
    setTargetCoord(targetX[tarInd],targetY[tarInd],0);
}

/*!\brief change target coordinates
 *
 * The pointing direction is expected in radian.
 * If given in arcsec, conversion flag: toRad must be set to 1
 *
 *@param[in] aX    : T   : target direction (X axis)
 *@param[in] aY    : T   : target direction (Y axis)
 *@param[in] toRad : int : convert pointing direction to radian
 * */
template<typename T>
void Tomo_struct<T>::setTargetCoord(T aX, T aY, int toRad){
    if(toRad){
        aX/=RASC;
        aY/=RASC;
    }
    alphaX[sys.nW-1]=aX;
    alphaY[sys.nW-1]=aY;
    //Computes  u and v
    projectedSubap();
    //Compute noise
    updateNoise();
}




/*!\brief Compute the seeing
 *
 * update the seeing for the ngs and lgs wavelength
 *
 * */
template<typename T>
void Tomo_struct<T>::updateSeeing(){
    int i;
    const T lambda_m=500.*1.e-9; //(meters) //TODO as param?
    T r0_500nm=std::accumulate(atm.cn2.begin(), atm.cn2.end(), 0.);
    //r0 in meters
    r0_500nm=pow(r0_500nm,-3./5.);

    // r0 = (lambda/lambda_500nm)^(6/5)*r0_500nm
    // s = lambda/r
    // s = lambda_500nm^(6/5) / (lambda^(1/5)*r0_500nm)
    // CONVERT IN ARCSEC
    sNGS2=lambda_m/r0_500nm*pow(lambda_m/sys.lambdaNGS,0.2)*RASC;
    //get square of seeing
    sNGS2*=sNGS2;

    //same for lgs
    sLGS2=lambda_m/r0_500nm*pow(lambda_m/sys.lambdaLGS,0.2)*RASC;
    sLGS2*=sLGS2;
}


/*!\brief Compute the WFS noise for each WFS
 *
 * */
template<typename T>
void Tomo_struct<T>::updateNoise(){
    noiseNGS    =std::vector<T>(sys.nNgs+1);
    noiseLGSxx  =std::vector<T>(Nx);
    noiseLGSyy  =std::vector<T>(Nx);
    noiseLGSxy  =std::vector<T>(Nx);
    //update square of seeing
    updateSeeing();

    int i,j,k=0;
    for(i=0;i<sys.nW;i++){
        const T d=(sys.diam/Nssp[i]); //already taken into account in F ?
        const T F=flux(qr, mr[i],throughput[i], d, sys.tFrame,bdw);

        const T fact1=0.32/F;
        const T fact2=0.82*pow(sys.RON/(F*pixSize[i]),2);


        T noiseLongAxis;
        if(i<sys.nLgs){
            const T partialExt=sLGS2 + (sys.spotWidth*sys.spotWidth);
            const T noiseShortAxis=(fact1*partialExt+
                                         fact2* pow((partialExt),1.5)* pow(partialExt,0.5));
            for(j=0;j<Nsubap[i];j++){
                noiseLongAxis=(fact1*(partialExt+lgsExt[k]*lgsExt[k])+
                              fact2* pow((partialExt+lgsExt[k]*lgsExt[k]),1.5)* pow(partialExt,0.5));
                noiseLGSxx[k]=noiseLongAxis* pow(cos(lgsTheta[k]),2) + noiseShortAxis*pow(sin(lgsTheta[k]),2);
                noiseLGSyy[k]=noiseLongAxis* pow(sin(lgsTheta[k]),2) + noiseShortAxis*pow(cos(lgsTheta[k]),2);
                noiseLGSxy[k]=(noiseLongAxis-noiseShortAxis)*sin(lgsTheta[k])*cos(lgsTheta[k]);
		        k++;
            }
        }
        else{
	        noiseNGS[i-sys.nLgs]=(fact1*sNGS2+fact2*sNGS2*sNGS2);
        }
    }
}

/*!\brief Generate the position (X,Y) of each subapertures of each WFS on the telescope pupil and the number of subapertures
 *
 * */
template<typename T>
void Tomo_struct<T>::validSubap()
{
    std::cout<<"    generating valid subapertures...\n"<< std::flush;
    Nsubap=std::vector<long>(sys.nW);

    const T bornemin = -sys.diam / 2.;
    const T Rtel2 = (sys.diam * sys.diam) / 4.;
    long NsubapTot = 0;
    long n;

    //Total number of subapertures (without obstruction)
    for (n = 0; n < sys.nW; n++) {
        NsubapTot += Nssp[n] * Nssp[n];
    }

    const long cNsubapTot = NsubapTot;
    T x[cNsubapTot], y[cNsubapTot];
    int index[cNsubapTot];

    int cpt = 0;
    int ioff = 0;

    //Computation of all the subapertures' positions
    for (n = 0; n < sys.nW; n++) {
        long Nsap = 0;
        T pas = sys.diam / (1. * Nssp[n]);
        int i;
        T Robs2;

        // to avoid some bug that eliminates useful central subapertures when obs=0.286
        if (Nssp[n] != 7 || (sys.obs <= 0.285 || sys.obs >= 0.29)) {
            Robs2 = sys.diam * sys.obs / 2. * sys.diam * sys.obs / 2.;
        } else {
            Robs2 = sys.diam * 0.285 / 2. * sys.diam * 0.285 / 2.;
        }

        if (Nssp[n] != 1) {
            for (i = 0; i < Nssp[n]; i++) {
                T tp = bornemin + pas / 2. * (2. * i + 1.); // y-coord of current subap
                int j;

                for (j = 0; j < Nssp[n]; j++) {
                    x[ioff + j] = bornemin + pas / 2. * (2. * j + 1.); // x-coord of current subap
                    y[ioff + j] = tp;

                    T r2 = x[ioff + j] * x[ioff + j] + y[ioff + j] * y[ioff + j];

                    //Search the non-valid subapertures
                    if (r2 < Robs2 || r2 >= Rtel2) {
                        index[cpt] = j + ioff; //list of the useless subapertures index
                        cpt++;
                    }
                    else {
                        Nsap++;
                    }
                }
                ioff += Nssp[n];
            }
            Nsubap[n] = Nsap;
        } else { //Special case (Nssp = 1)
            x[ioff] = 0.; // x-coord of current subap
            y[ioff] = 0.;
            ioff += Nssp[n];
            Nsubap[n] = 1;
        }
    }

    Nx = cNsubapTot-cpt;
    X = std::vector<T>(Nx);
    Y = std::vector<T>(Nx);
    std::cout<<"        "<<Nx<<" valid subapertures out of "<<cNsubapTot<<"\n"<< std::flush;
    
    int a = 0;
    int off = 0;
    int borne = 0;
    int i;
    //Suppress the non-valid subapertures
    while (a <= cpt) {
        if (a == cpt) {
            borne = cNsubapTot;
        } else {
            borne = index[a];
        }
        for (i = off; i < borne; i++) {
            X[i - a] = x[i];
            Y[i - a] = y[i];
        }
        off = index[a] + 1;
        a++;
    }
    std::cout<<"    done"<<std::endl;
}

/*!\brief generate the elongation of the laser spots for each subaperture
 *
 *
 */
template<typename T>
void Tomo_struct<T>::generateElongation(){
    std::cout<<"    generate LGS elongation... "<< std::flush;
    lgsExt      =std::vector<T>(Nx);
    lgsTheta    =std::vector<T>(Nx);

    int w,i,k=0;
    const T toFwhm=RASC*sys.lgsDepth/(sys.lgsAlt*sys.lgsAlt);
    for(w=0;w<sys.nLgs;w++){
        const T lltx=alphaX[w]*(sys.diam/2./sqrt(alphaX[w]*alphaX[w]+alphaY[w]*alphaY[w]));
        const T llty=alphaY[w]*(sys.diam/2./sqrt(alphaX[w]*alphaX[w]+alphaY[w]*alphaY[w]));
        for(i=0;i<Nsubap[w];i++){
            lgsExt[k]=sqrt((pow(X[k]-lltx,2)+pow(Y[k]-llty,2)))*toFwhm;
            lgsTheta[k]=atan2((Y[k]-llty)*toFwhm,(X[k]-lltx)*toFwhm);
            k++;
        }
    }
    std::cout<<"    done"<<std::endl;
}


/*!\brief Computes the projected coordinates of all subapertures  projected onto all the layer
 *
 */
template<typename T>
void Tomo_struct<T>::projectedSubap() {

  std::cout<<"    computing position of projected subapertures..."<< std::flush;
  long i, tid;
  long n = 0;
  long l;
  const T rad = 3.14159265358979323846 / 180.;
  
  //const long int cNsubap0 = Nsubap[0];
  //const long int cNw = sys.nW;

  long ioff[sys.nW];
  ioff[0] = 0;
  for (i = 1; i < sys.nW; i++) {
    ioff[i] = ioff[i-1] + Nsubap[i-1];
  }
  
  for (tid = 0; tid < atm.nLayer * Nx; tid++) {

    n = 0;
    l = tid / Nx;
    const int pos = tid - l * Nx;
    long Nsubapx = Nsubap[0];

    while(pos >= Nsubapx){
      n++;
      Nsubapx += Nsubap[n];
    }
    Nsubapx -= Nsubap[n];

    i = pos - Nsubapx;
    
    const T dX = alphaX[n] * atm.h[l];
    const T dY = alphaY[n] * atm.h[l];

    const T rr = 1. - atm.h[l] * sys.gsAlt[n];

    //const long Nsap = Nsubap[n];
    const long nssp = Nssp[n];

    //magnification factor
    const T G = diamPup[n] / (T) (nssp);

    //rotation angle
    const T th = sys.thetaML[n] * rad;

    //taking magnification factor into account
    const T xtp = X[ioff[n] + i] * G;
    const T ytp = Y[ioff[n] + i] * G;

    //taking rotation into account
    T uu = xtp * cos(th) - ytp * sin(th);
    T vv = xtp * sin(th) + ytp * cos(th);

    //taking pupil offset into account
    uu += sys.XPup[n];
    vv += sys.YPup[n];

    //Projection onto  the layer
    u[tid] = uu * rr + dX;
    v[tid] = vv * rr + dY;
  }
  std::cout<<"    done"<<std::endl;
}

/*!\brief return the number of measurements of the actual WFS
 */
template<typename T> int 
long Tomo_struct<T>::getNMeas(){
    return 2*(Nx-Nsubap[sys.nW-1]);
}

/*!\brief return the number of measurements of the truth sensor
 */
template<typename T> 
long Tomo_struct<T>::getNMeasTS(){
    return 2*Nsubap[sys.nLgs+sys.nNgs];
}

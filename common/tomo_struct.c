/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "tomo_struct.h"
#include "utils.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include<stdio.h>
#include "noise.h"

#ifdef __cplusplus
extern "C" {
#endif

  int matcov_init_tomo_tiled(struct tomo_struct *tomo, char* sys_path, char* atm_path, int night_idx, int snapshots_per_night, int snapshot_idx, int obs_idx, real_t alphaX, real_t alphaY){
    
      tomo->Nlayer=NAN;
    FPRINTF(stdout, "Initializing matcov...");fflush(stdout);

    strcpy( tomo->sys_path, sys_path);
    strcpy( tomo->atm_path, atm_path);
    //initialize tomo struct on cpu
    if(!init_tomo_sys(tomo)) return 0;
    tomo->alphaX[tomo->Nw-1] = alphaX / 206265.0; // convert to radian;
    tomo->alphaY[tomo->Nw-1] = alphaY / 206265.0; // convert to radian;
    if(!matcov_update_atm_params(tomo, night_idx, snapshots_per_night, snapshot_idx, obs_idx)) return 0;
                  
    FPRINTF(stdout, "done...\n");fflush(stdout);
    
    return 1;
  }

//------------------------------------------------------------------------------------
/*!
 * DOCUMENT         subap_position(tomo, u, v)
   <tomo>                : structure with all the needed information.
   <u> and <v>           : 3d arrays containing the sub-apertures projected coordinates onto all the layers. u[0][2][1] is the X-coordinate of the subap 2 of the WFS 0 on the layer 1.

   Computes the projected coordinates of all subapertures  projected onto all the layer
  */
void subap_position(struct tomo_struct *tomo) {

  FPRINTF(stdout, "Starting subap_position...");fflush(stdout);
  long i, tid;
  long n = 0;
  long l;
  const real_t rad = 3.14159265358979323846 / 180.;
  
  //const long int cNsubap0 = tomo->Nsubap[0];
  //const long int cNw = tomo->Nw;

  long ioff[tomo->Nw];
  ioff[0] = 0;
  for (i = 1; i < tomo->Nw; i++) {
    ioff[i] = ioff[i-1] + tomo->Nsubap[i-1];
  }
  
  for (tid = 0; tid < tomo->Nlayer * tomo->Nx; tid++) {

    n = 0;
    l = tid / tomo->Nx;
    const int pos = tid - l * tomo->Nx;
    long Nsubapx = tomo->Nsubap[0];

    while(pos >= Nsubapx){
      n++;
      Nsubapx += tomo->Nsubap[n];
    }
    Nsubapx -= tomo->Nsubap[n];

    i = pos - Nsubapx;
    
    const real_t dX = tomo->alphaX[n] * tomo->h[l];
    const real_t dY = tomo->alphaY[n] * tomo->h[l];

    const real_t rr = 1. - tomo->h[l] * tomo->GsAlt[n];

    //const long Nsap = tomo->Nsubap[n];
    const long nssp = tomo->Nssp[n];

    //magnification factor
    const real_t G = tomo->diamPup[n] / (real_t) (nssp);

    //rotation angle
    const real_t th = tomo->thetaML[n] * rad;

    //taking magnification factor into account
    const real_t xtp = tomo->X[ioff[n] + i] * G;
    const real_t ytp = tomo->Y[ioff[n] + i] * G;

    //taking rotation into account
    real_t uu = xtp * cos(th) - ytp * sin(th);
    real_t vv = xtp * sin(th) + ytp * cos(th);

    //taking pupil offset into account
    uu += tomo->XPup[n];
    vv += tomo->YPup[n];

    //Projection onto  the layer
    tomo->u[tid] = uu * rr + dX;
    tomo->v[tid] = vv * rr + dY;
  }
  FPRINTF(stdout, "done...\n");fflush(stdout);
}

//------------------------------------------------------------------------------------
void generateXY(struct tomo_struct *tomo)
/* DOCUMENT  generateXY(struct tomo_struct tomo, double *Nsubap)
 <tomo>               :  structure with all the needed information
 <tomo.X> & <tomo.Y>            :   arrays containing all the sub-apertures coordinates of all WFS, one after the other
<tomo.Nsubap>              :  number of subaperture of ezach WFS
 Generate the position (X,Y) of each subapertures of each WFS on the telescope pupil and the number of subaperture of ezach WFS (Nsubap) 
 */
{
  const real_t bornemin = -tomo->DiamTel / 2.;
  const real_t Rtel2 = (tomo->DiamTel * tomo->DiamTel) / 4.;
  long NsubapTot = 0;
  long n;

  //Total number of subapertures (without obstruction)
  for (n = 0; n < tomo->Nw; n++) {
    NsubapTot += tomo->Nssp[n] * tomo->Nssp[n];
  }

  const long cNsubapTot = NsubapTot;
  real_t x[cNsubapTot], y[cNsubapTot];
  int index[cNsubapTot];

  int cpt = 0;
  int ioff = 0;

  //Computation of all the subapertures' positions
  for (n = 0; n < tomo->Nw; n++) {
    long Nsap = 0;
    real_t pas = tomo->DiamTel / (1. * tomo->Nssp[n]);
    int i;
    real_t Robs2;

    // to avoid some bug that eliminates useful central subapertures when obs=0.286
    if (tomo->Nssp[n] != 7 || (tomo->obs <= 0.285 || tomo->obs >= 0.29)) {
      Robs2 = tomo->DiamTel * tomo->obs / 2. * tomo->DiamTel * tomo->obs / 2.;
    } else {
      Robs2 = tomo->DiamTel * 0.285 / 2. * tomo->DiamTel * 0.285 / 2.;
    }

    if (tomo->Nssp[n] != 1) {
      for (i = 0; i < tomo->Nssp[n]; i++) {
        real_t tp = bornemin + pas / 2. * (2. * i + 1.); // y-coord of current subap
        int j;

        for (j = 0; j < tomo->Nssp[n]; j++) {
          x[ioff + j] = bornemin + pas / 2. * (2. * j + 1.); // x-coord of current subap
          y[ioff + j] = tp;

          real_t r2 = x[ioff + j] * x[ioff + j] + y[ioff + j] * y[ioff + j];

          //Search the non-valid subapertures
          if (r2 < Robs2 || r2 >= Rtel2) {
            index[cpt] = j + ioff; //list of the useless subapertures index
            cpt++;
          }
          else {
            Nsap++;
          }
        }
        ioff += tomo->Nssp[n];
      }
      tomo->Nsubap[n] = Nsap;
   } else { //Special case (Nssp = 1)
      x[ioff] = 0.; // x-coord of current subap
      y[ioff] = 0.;
      ioff += tomo->Nssp[n];
      tomo->Nsubap[n] = 1;
    }
  }

  tomo->X=(real_t*)malloc((cNsubapTot-cpt)*sizeof(real_t));
  tomo->Y=(real_t*)malloc((cNsubapTot-cpt)*sizeof(real_t));
  tomo->Nx = cNsubapTot-cpt;
  fprintf(stdout,"\ncNsubapTot %d  cpt %d, NX: %d\n",cNsubapTot,cpt,tomo->Nx);
  
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
      tomo->X[i - a] = x[i];
      tomo->Y[i - a] = y[i];
    }

    off = index[a] + 1;
    a++;
  }
}

//------------------------------------------------------------------------------------
int init_tomo_sys(struct tomo_struct *tomo){
  int i;
  long nssp_tmp;
  double pasDPHI;

  char sys_filename[512];
  sprintf(sys_filename, "%ssys-params.txt", tomo->sys_path);
  FPRINTF(stdout, "Opening file %s ", sys_filename);
  
  FILE *file = fopen(sys_filename, "r");
  if(!file){
    fprintf(stderr, "ERROR: not able to open file %s\n", sys_filename);
    return 0;
  }
  read_paramd(file,&(tomo->DiamTel));   // telescope diameter
  read_paramd(file,&(tomo->obs));       // central obscuration
  read_paramd(file,&(tomo->Tframe));    // 
  read_paraml(file,&(tomo->Nw));        // number of wavefront sensors
  read_parami(file,&(tomo->nlgs));      // number of laser guide stars
  read_paraml(file,&(tomo->nTarget));   // number of targets

  tomo->Nssp        = (long*)   malloc(tomo->Nw*sizeof(long));
  tomo->Nsubap      = (long*)   malloc(tomo->Nw*sizeof(long));// array of the number of subap of each WFS, contains Nw elements
  tomo->GsAlt       = (real_t*) malloc(tomo->Nw*sizeof(real_t));
  tomo->type        = (int*)    malloc(tomo->Nw*sizeof(int));
  tomo->alphaX      = (real_t*) malloc(tomo->Nw*sizeof(real_t));
  tomo->alphaY      = (real_t*) malloc(tomo->Nw*sizeof(real_t));
  tomo->XPup        = (real_t*) malloc(tomo->Nw*sizeof(real_t)); // pupil shift of the WFS, in meters
  tomo->YPup        = (real_t*) malloc(tomo->Nw*sizeof(real_t)); // pupil shift of the WFS, in meters
  tomo->thetaML     = (real_t*) malloc(tomo->Nw*sizeof(real_t));  // rotation of microlenses
  tomo->thetaCam    = (real_t*) malloc(tomo->Nw*sizeof(real_t)); // rotation of microlenses
  tomo->sensibilite = (real_t*) malloc(tomo->Nw*sizeof(real_t)); // sensitivity coeff of this WFS
  tomo->diamPup     = (real_t*) malloc(tomo->Nw*sizeof(real_t)); 
  tomo->sspSize     = (real_t*) malloc(tomo->Nw*sizeof(real_t)); // subaperture size of this WFS
  tomo->tracking    = (real_t*) malloc(3*sizeof(real_t));
  tomo->mr          = (real_t*) calloc(tomo->Nw,sizeof(real_t));
  tomo->lgsFlux     = (real_t*) calloc(tomo->nlgs,sizeof(real_t));
  tomo->pixSize     = (real_t*) calloc(tomo->Nw,sizeof(real_t));
  tomo->throughput  = (real_t*) calloc(tomo->Nw,sizeof(real_t));
  tomo->targetX= (real_t*)malloc(tomo->nTarget*sizeof(real_t));
  tomo->targetY= (real_t*)malloc(tomo->nTarget*sizeof(real_t));

  // array of the number of subap of each WFS along the telescop diameter, contains Nw elements
  read_paraml(file,&nssp_tmp); // number of wavefront sensors
  //read_arrayl(file,tomo->Nssp,tomo->Nw);
  // array of the inverse of the guide star altitude (in 1/meters), contains Nw elements
  read_arrayd(file,tomo->GsAlt,tomo->Nw);
  // type of WFS, 0, 1, 2 or 3. 0 is unused, 1=NGS, 2=LGS, 3=TipTilt-guide star
  read_arrayi(file,tomo->type,tomo->Nw);
  // Pointing directions of WFS
  read_arrayd(file,tomo->alphaX,tomo->Nw);
  read_arrayd(file,tomo->alphaY,tomo->Nw);
  // Deviations of WFSs
  read_arrayd(file,tomo->XPup,tomo->Nw);
  read_arrayd(file,tomo->YPup,tomo->Nw);
  read_arrayd(file,tomo->thetaML,tomo->Nw);
  read_arrayd(file,tomo->thetaCam,tomo->Nw);
  read_arrayd(file,tomo->sensibilite,tomo->Nw);
  // telescope tracking error parameters (x^2, y^2 and xy), units : arcsec^2
  read_arrayd(file,tomo->tracking,3);
  read_paramd(file,&(pasDPHI));
  read_parami(file,&(tomo->ncpu));
  //magnitude parameters
  read_arrayd(file,&(tomo->mr[tomo->nlgs]),tomo->Nw-tomo->nlgs);
  read_paramd(file,&(tomo->lgsFlux[0]));//Photon return at M1 (ph/m2/s)
  //pixSize
  read_paramd(file,&(tomo->pixSize[tomo->nlgs]));//pixsize NGS
  read_paramd(file,&(tomo->pixSize[0]));         //pixsize LGS
  //wavelength parameters
  read_paramd(file,&(tomo->lambdaNGS));
  read_paramd(file,&(tomo->lambdaLGS));
  read_paramd(file,&(tomo->bdw));
  //throughput
  read_paramd(file,&(tomo->throughput[tomo->nlgs])); //throughput NGS
  read_paramd(file,&(tomo->throughput[0]));          //throughput LGS
  read_paramd(file,&(tomo->Tatmo));
  // noise stuff variables
  read_parami(file,&(tomo->RON));
  read_paramd(file,&(tomo->lgs_cst));
  read_paramd(file,&(tomo->spot_width)); 
  read_paramd(file,&(tomo->lgs_alt));
  read_paramd(file,&(tomo->lgs_depth));
  // Pointing directions of targets
  read_arrayd(file,tomo->targetX,tomo->nTarget);
  read_arrayd(file,tomo->targetY,tomo->nTarget);
  fclose(file);


  for(i=0; i<tomo->Nw; i++){
    tomo->Nssp[i]    = nssp_tmp; // all wfs have same number of subaps from the parameter file
    tomo->alphaX[i] /= 206265.0; // convert to radian
    tomo->alphaY[i] /= 206265.0; 
    tomo->diamPup[i] =(real_t)tomo->Nssp[i];
    tomo->sspSize[i] =(real_t)tomo->DiamTel/tomo->diamPup[i];
  }

  //computation of flux 
  tomo->qr=2.9;
  tomo->bdw=tomo->bdw*1.e10; //convert to Angstrom
  tomo->lgsFlux[0]=tomo->lgsFlux[0]*tomo->Tframe*pow(tomo->DiamTel/tomo->Nssp[0],2)*tomo->throughput[0]/tomo->Tatmo;  //e-/subs/frame
  for(i=1;i<tomo->nlgs;i++){
      tomo->lgsFlux[i]=tomo->lgsFlux[0];
      tomo->pixSize[i]=tomo->pixSize[0];
      tomo->throughput[i]=tomo->throughput[0];
      tomo->mr[i]=magnitudeFromFlux(tomo->qr,tomo->lgsFlux[i],tomo->throughput[i],tomo->DiamTel/tomo->Nssp[i], tomo->Tframe, tomo->bdw);
  }
  for(i=tomo->nlgs+1;i<tomo->Nw;i++){
      tomo->pixSize[i]=tomo->pixSize[tomo->nlgs];
      tomo->throughput[i]=tomo->throughput[tomo->nlgs];
  }
      tomo->mr[0]=magnitudeFromFlux(tomo->qr,tomo->lgsFlux[0],tomo->throughput[0],tomo->DiamTel/tomo->Nssp[0], tomo->Tframe, tomo->bdw);
  

  //Generate the subapertures positions and fill tomo.Nsubap
  generateXY(tomo);

  for(i=0;i<tomo->nlgs;i++){
      tomo->type[i]=2;
      tomo->GsAlt[i]=1/tomo->lgs_alt;
  }
  for(i=tomo->nlgs;i<tomo->Nw;i++){
      tomo->type[i]=1;
      tomo->GsAlt[i]=0;
  }


  //computation of laser elongation for the noise
  tomo->noiseNGS=(real_t*)calloc(tomo->Nw-tomo->nlgs,sizeof(real_t));
  int Nsubap=0;
  for(i=0;i<tomo->nlgs;i++){
  Nsubap+=tomo->Nsubap[i];
  }
  tomo->noiseLGSxx  =(real_t*)calloc(Nsubap,sizeof(real_t));
  tomo->noiseLGSyy  =(real_t*)calloc(Nsubap,sizeof(real_t));
  tomo->noiseLGSxy  =(real_t*)calloc(Nsubap,sizeof(real_t));
  tomo->lgsExt      =(real_t*)calloc(Nsubap,sizeof(real_t));
  tomo->lgsTheta    =(real_t*)calloc(Nsubap,sizeof(real_t));
  
  generateElongation(tomo->X,tomo->Y,tomo->alphaX,tomo->alphaY,tomo->nlgs,tomo->DiamTel,tomo->lgsExt,tomo->lgsTheta,tomo->lgs_depth,tomo->lgs_alt,tomo->Nsubap);

    tomo->cn2=NULL;
    tomo->h=NULL;
    tomo->L0=NULL;
    tomo->indexL0=NULL;
    tomo->u=NULL;
    tomo->v=NULL;
    tomo->sspSizeL=NULL;
    tomo->L0diff=NULL;

  return 1;
}

//------------------------------------------------------------------------------------
int init_tomo_atm(struct tomo_struct *tomo){
  tomo->cn2         = (real_t*)malloc(tomo->Nlayer*sizeof(real_t)); // profile strengh, units TBD ..
  tomo->h           = (real_t*)malloc(tomo->Nlayer*sizeof(real_t)); // altitude of layers (meters)
  tomo->L0          = (real_t*)malloc(tomo->Nlayer*sizeof(real_t)); // outer scale (meters)
  tomo->indexL0     = (long*)malloc(tomo->Nlayer * sizeof(long));
  tomo->u           = (real_t*)malloc(tomo->Nlayer * tomo->Nx * sizeof(real_t));
  tomo->v           = (real_t*)malloc(tomo->Nlayer * tomo->Nx * sizeof(real_t));
  tomo->sspSizeL    = (real_t*)malloc(tomo->Nw * tomo->Nlayer * sizeof(real_t));
  tomo->L0diff      = (real_t*)malloc(tomo->Nlayer*sizeof(real_t));
  return 1;
}

  //------------------------------------------------------------------------------------
  int matcov_update_atm_params(struct tomo_struct *tomo, int night_idx, int snapshots_per_night, int snapshot_idx, int obs_idx){
    FPRINTF(stdout, "Updating matcov atm-params...");fflush(stdout);

    long int Nlayer;
    char atm_filename[512];
    sprintf(atm_filename, "%sprof%d-atmos-night%d.txt", tomo->atm_path, snapshot_idx * snapshots_per_night + obs_idx, night_idx);

    FPRINTF(stdout, "Opening file %s ", atm_filename);
    
    FILE *file = fopen(atm_filename, "r");
    if(!file){
      fprintf(stderr, "ERROR: not able to open file %s!\n", atm_filename);
      return 0;
    }
      
    read_paraml(file,&Nlayer); // number of layers in the profile
    if(Nlayer != tomo->Nlayer){
        tomo->Nlayer=Nlayer;
        free_tomo_atm(tomo);
        init_tomo_atm(tomo);
    }
    read_paramd(file,&(tomo->r0)); 
    read_arrayd(file,tomo->cn2,tomo->Nlayer);
    read_arrayd(file,tomo->h,tomo->Nlayer);
    read_arrayd(file,tomo->L0,tomo->Nlayer);
    
    fclose(file);
    FPRINTF(stdout, "done...\n");fflush(stdout);

    matcov_update_tomo_tiled(tomo,0); // target the first element of the target list


    return 1;
  }
//  //------------------------------------------------------------------------------------
//  //equivalent to update_tomo_atm_gpu_gb()
  void matcov_update_tomo_tiled(struct tomo_struct *tomo, int t){

    if(t>-1){
        tomo->alphaX[tomo->Nw-1] = tomo->targetX[t];
        tomo->alphaY[tomo->Nw-1] = tomo->targetY[t];
    }

    FPRINTF(stdout, "Starting matcov_update_tomo_tiled...");fflush(stdout);
    FPRINTF(stdout, "Starting cc loop...");fflush(stdout);
    int cc;
    for (cc = 0; cc < tomo->Nw * tomo->Nlayer; cc++) {
      int n = cc / tomo->Nlayer;
      int l = cc - n * tomo->Nlayer;
      if(n >= tomo->Nw) n-=1;
      tomo->sspSizeL[cc] = tomo->sspSize[n] * (1. - tomo->GsAlt[n] * tomo->h[l]);
    }
    FPRINTF(stdout, "done...\n");fflush(stdout);
    FPRINTF(stdout, "Starting tmp loop...");fflush(stdout);
    //Search the different L0 and build indexL0
    const long cNlayer = tomo->Nlayer;
    long i, j;
    int cpt = 1;
    real_t tmp[cNlayer];

    tmp[0] = tomo->L0[0];
    tomo->indexL0[0] = 0;
    
    FPRINTF(stdout, "cpt: ...");fflush(stdout);
    for (i = 1; i < cNlayer; i++) {
      j = 0;
      const real_t l0 = tomo->L0[i];
      
      while ((j < cpt) && (tmp[j] != l0)){
          j++;
      }
      tomo->indexL0[i] = j;      
      if (j == cpt) {
        tmp[j] = l0;
        cpt++;
      }
    }
    FPRINTF(stdout, "%d done...\n",cpt);fflush(stdout);
    for (i = 0; i < cNlayer; i++)  {
      tomo->L0diff[i] = tomo->L0[i];
    }
    FPRINTF(stdout, "done...\n");fflush(stdout);
    //Computes  u and v
    subap_position(tomo);
    //Compute noise
    updateNoise(tomo);
    FPRINTF(stdout, "done...\n");fflush(stdout);
  }
  //------------------------------------------------------------------------------------
  int matcov_getNumMeasurements(struct tomo_struct *tomo){
    return tomo->Nx * 2;
    //return tomo->nsubaps_offaxis * 2;
  }
  //------------------------------------------------------------------------------------
  int matcov_getNumMeasurementsTS(struct tomo_struct *tomo){
    return tomo->Nsubap[tomo->Nw-1] * 2;
  }


  //------------------------------------------------------------------------------------
  void matcov_free_tomo_tiled(struct tomo_struct *tomo){

    //free_tomo_tile(tomo);
    free_tomo_sys(tomo);
    free_tomo_atm(tomo);
  }


void free_tomo_sys(struct tomo_struct *tomo){
    if(tomo->Nssp!=NULL)        {free(tomo->Nssp);          tomo->Nssp=NULL;} 
    if(tomo->Nsubap!=NULL)      {free(tomo->Nsubap);        tomo->Nsubap=NULL;} 
    if(tomo->GsAlt!=NULL)       {free(tomo->GsAlt);         tomo->GsAlt=NULL;} 
    if(tomo->type!=NULL)        {free(tomo->type);          tomo->type=NULL;} 
    if(tomo->alphaX!=NULL)      {free(tomo->alphaX);        tomo->alphaX=NULL;} 
    if(tomo->alphaY!=NULL)      {free(tomo->alphaY);        tomo->alphaY=NULL;} 
    if(tomo->XPup!=NULL)        {free(tomo->XPup);          tomo->XPup=NULL;} 
    if(tomo->YPup!=NULL)        {free(tomo->YPup);          tomo->YPup=NULL;} 
    if(tomo->thetaML!=NULL)     {free(tomo->thetaML);       tomo->thetaML=NULL;} 
    if(tomo->thetaCam!=NULL)    {free(tomo->thetaCam);      tomo->thetaCam=NULL;} 
    if(tomo->sensibilite!=NULL) {free(tomo->sensibilite);   tomo->sensibilite=NULL;} 
    if(tomo->diamPup!=NULL)     {free(tomo->diamPup);       tomo->diamPup=NULL;} 
    if(tomo->sspSize!=NULL)     {free(tomo->sspSize);       tomo->sspSize=NULL;} 
    if(tomo->tracking!=NULL)    {free(tomo->tracking);      tomo->tracking=NULL;} 
    if(tomo->mr!=NULL)          {free(tomo->mr);            tomo->mr=NULL;} 
    if(tomo->lgsFlux!=NULL)     {free(tomo->lgsFlux);       tomo->lgsFlux=NULL;} 
    if(tomo->pixSize!=NULL)     {free(tomo->pixSize);       tomo->pixSize=NULL;} 
    if(tomo->throughput!=NULL)  {free(tomo->throughput);    tomo->throughput=NULL;} 
    if(tomo->targetX!=NULL)     {free(tomo->targetX);       tomo->targetX=NULL;} 
    if(tomo->targetY!=NULL)     {free(tomo->targetY);       tomo->targetY=NULL;} 
    if(tomo->noiseNGS!=NULL)    {free(tomo->noiseNGS);      tomo->noiseNGS=NULL;} 
    if(tomo->noiseLGSxx!=NULL)  {free(tomo->noiseLGSxx);    tomo->noiseLGSxx=NULL;} 
    if(tomo->noiseLGSyy!=NULL)  {free(tomo->noiseLGSyy);    tomo->noiseLGSyy=NULL;} 
    if(tomo->noiseLGSxy!=NULL)  {free(tomo->noiseLGSxy);    tomo->noiseLGSxy=NULL;} 
    if(tomo->lgsExt!=NULL)      {free(tomo->lgsExt);        tomo->lgsExt=NULL;} 
    if(tomo->lgsTheta!=NULL)    {free(tomo->lgsTheta);      tomo->lgsTheta=NULL;} 
    if(tomo->X!=NULL)           {free(tomo->X);             tomo->X=NULL;} 
    if(tomo->Y!=NULL)           {free(tomo->Y);             tomo->Y=NULL;} 
}
void free_tomo_atm(struct tomo_struct *tomo){
    if(tomo->cn2 !=NULL)        {free(tomo->cn2);       tomo->cn2=NULL;}
    if(tomo->h !=NULL)          {free(tomo->h);         tomo->h=NULL;}
    if(tomo->L0 !=NULL)         {free(tomo->L0);        tomo->L0=NULL;}
    if(tomo->indexL0 !=NULL)    {free(tomo->indexL0);   tomo->indexL0=NULL;}
    if(tomo->u !=NULL)          {free(tomo->u);         tomo->u=NULL;}
    if(tomo->v !=NULL)          {free(tomo->v);         tomo->v=NULL;}
    if(tomo->sspSizeL!=NULL)    {free(tomo->sspSizeL);  tomo->sspSizeL=NULL;}
    if(tomo->L0diff != NULL)    {free(tomo->L0diff);    tomo->L0diff=NULL;}

}
//------------------------------------------------------------------------------------
void free_tomo(struct tomo_struct *tomo){
    free_tomo_sys(tomo);
    free_tomo_atm(tomo);
}

  //------------------------------------------------------------------------------------
  void matcov_set_gal_coords(struct tomo_struct *tomo, real_t alphaX, real_t alphaY){
    FPRINTF(stdout, "Setting matcov galaxy params...");fflush(stdout);
    tomo->alphaX[tomo->Nw-1] = alphaX/ 206265.0; // convert to radian;
    tomo->alphaY[tomo->Nw-1] = alphaY/ 206265.0; // convert to radian;
    
    real_t dmax = 0.0;
    real_t maxalt = tomo->h[tomo->Nlayer-1];
    int cc;
    for (cc=0;cc<tomo->Nw;cc++) {
      real_t tmp = sqrtf(tomo->alphaX[cc]*tomo->alphaX[cc] + tomo->alphaY[cc]*tomo->alphaY[cc]);
      if (tmp > dmax) dmax = tmp;
    }
    tomo->rmax = dmax * 2 * maxalt + tomo->DiamTel;
    
    //subap_position(tomo);
    matcov_update_tomo_tiled(tomo,-1); // use the current target
    FPRINTF(stdout, "done...\n");fflush(stdout);
  }

#ifdef __cplusplus
}
#endif

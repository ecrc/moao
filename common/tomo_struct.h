/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef TOMO_STRUCT_H
#define TOMO_STRUCT_H

#ifdef USE_CHAMELEON
#include <morse.h>
#endif

#include "moao_defs.h"

struct tomo_struct {
  long Nw;  // number of wavefront sensors
  long nTarget; //number of targets

  // pointers on arrays containing corrdinates of sub-apertures
  // X and Y are biiiiig arrays from yorick containing all the subap
  // coordinates of all the WFSs, one after the other.
  real_t *X;
  real_t *Y;

  real_t DiamTel; // telescope Diameter

  real_t obs;//telescope obscuration

  // array of the number of subap of each WFS, contains Nw elements
  long *Nsubap;

  // array of the number of subap of each WFS along the telescop diameter, contains Nw elements
  long *Nssp;

  // array of the inverse of the guide star altitude (in 1/meters), contains Nw elements
  real_t *GsAlt;

  // type of WFS, 0, 1, 2 or 3. 0 is unused, 1=NGS, 2=LGS, 3=TipTilt-guide star
  int *type;

  // Pointing directions of WFS
  real_t *alphaX;           // pointing direction in X, radian
  real_t *alphaY;           // pointing direction in Y, radian
  // Pointing directions of targets
  real_t *targetX;          // pointing direction in X, radian
  real_t *targetY;          // pointing direction in Y, radian

  // Deviations of WFSs
  real_t *XPup;             // pupil shift of the WFS, in meters
  real_t *YPup;             // pupil shift of the WFS, in meters
  real_t *thetaML;          // rotation of microlenses
  real_t *thetaCam;         // rotation of camera
  real_t *sensibilite;      // sensitivity coeff of this WFS
  real_t *diamPup;          // magnification factor of this WFS
  real_t *sspSize;          // subaperture size of this WFS

  // PROFILE
  long Nlayer;              // number of layers in the profile
  real_t r0;                // r0 at wfs lambda
  real_t *cn2;              // profile strengh, units as in  E-SPE-ESO-276-0206_atmosphericparameters
  real_t *h;                // altitude of layers (meters)
  real_t *L0;               // outer scale (meters)

  real_t rmax;              // maximum distance between subapertures (computed with yorick)
  real_t *tracking;         // telescope tracking error parameters (x^2, y^2 and xy), units : arcsec^2

  //real_t pasDPHI;           //Precision of DPHI precomputation.
  int ncpu;                 //Number of CPU used (only with openMP)
  int part;                 //Computed part of the cov. matrix. 0: complete 1: cmm 2: cpp 3: cpm ??
  int Nx;
  int Nslopes;
  int nlgs;

  real_t lgs_cst;
  real_t spot_width;
  real_t lgs_alt;
  real_t lgs_depth;

  real_t lambdaNGS;	    // lambda for the NGS in meter
  real_t lambdaLGS;         // lambda for the LGS in meter
  real_t sNGS2;             // square of the seeing at NGS lambda
  real_t sLGS2;             // square of the seeing at LGS lambda

  real_t *pixSize;

  real_t Tatmo;             // atmosphere transmission (in visible at z=30°
  real_t *throughput;       // throughput of NGS

  int RON;                  // read out noise (nb of e-)

  real_t qr;                // photometric Flux offset
  real_t *mr;               // guide star magnitude
  real_t bdw;               // bandwidth (A)

  real_t Tframe;
  real_t *lgsFlux;
  real_t *lgsExt;           //extension of lgs
  real_t *lgsTheta;         // angle of lgs extension
  real_t *noiseNGS;
  real_t *noiseLGSxx;
  real_t *noiseLGSyy;
  real_t *noiseLGSxy;

  //Ali
  long   *indexL0;
  real_t *L0diff;
  //real_t *tabDPHI;
  real_t *u;
  real_t *v;
  real_t *sspSizeL;
  long    nsubaps_offaxis;

  char    sys_path[512];
  char    atm_path[512];

};

#ifdef __cplusplus
extern "C" {
#endif
int matcov_init_tomo_tiled(struct tomo_struct *tomo, char* sys_path, char* atm_path, int night_idx, int snapshots_per_night, int snapshot_idx, int obs_idx, real_t alphaX, real_t alphaY);
void matcov_update_tomo_tiled(struct tomo_struct *tomo, int t);
int matcov_update_atm_params(struct tomo_struct *tomo, int night_idx, int snapshots_per_night, int snapshot_idx, int obs_idx);

void matcov_free_tomo_tiled(struct tomo_struct *tomo);
void free_tomo(struct tomo_struct *tomo);
void free_tomo_tile(struct tomo_struct *tomo);

int matcov_getNumMeasurements(struct tomo_struct *tomo);
int matcov_getNumMeasurementsTS(struct tomo_struct *tomo);

void matcov_set_gal_coords(struct tomo_struct *tomo, real_t alphaX, real_t alphaY);

int init_tomo_sys(struct tomo_struct *tomo);
//int init_tomo_atm(struct tomo_struct *tomo, int night_idx, int snapshots_per_night, int snapshot_idx, int obs_idx);
int init_tomo_atm(struct tomo_struct *tomo);
void free_tomo_sys(struct tomo_struct *tomo);
void free_tomo_atm(struct tomo_struct *tomo);
int init_tomo_tile(struct tomo_struct *tomo);
void generateXY(struct tomo_struct *tomo);

#ifdef __cplusplus
}
#endif

#endif //TOMO_STRUCT_H

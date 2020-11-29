
#include<stdio.h>
#include<ctype.h>
#include<stdlib.h>
#include<math.h>
#include<mkl_blas.h> 

#ifndef MATCOV
#define MATCOV

#define USE_OPENMP 1

struct tomo_struct {
  long Nw;  // number of wavefront sensors

  // pointers on arrays containing corrdinates of sub-apertures
  // X and Y are biiiiig arrays from yorick containing all the subap
  // coordinates of all the WFSs, one after the other.
  double *X;
  double *Y;

  double DiamTel; // telescope Diameter

  double obs;//telescope obscuration

  // array of the number of subap of each WFS, contains Nw elements
  long *Nsubap;

  // array of the number of subap of each WFS along the telescop diameter, contains Nw elements
  long *Nssp;

  // array of the inverse of the guide star altitude (in 1/meters), contains Nw elements
  double *GsAlt;

  // type of WFS, 0, 1, 2 or 3. 0 is unused, 1=NGS, 2=LGS, 3=TipTilt-guide star
  int *type;

  // Pointing directions of WFS
  double *alphaX;           // pointing direction in X, radian
  double *alphaY;           // pointing direction in Y, radian

  // Deviations of WFSs
  double *XPup;             // pupil shift of the WFS, in meters
  double *YPup;             // pupil shift of the WFS, in meters
  double *thetaML;          // rotation of microlenses
  double *thetaCam;         // rotation of camera
  double *sensibilite;      // sensitivity coeff of this WFS
  double *diamPup;          // magnification factor of this WFS
  double *sspSize;          // subaperture size of this WFS

  // PROFILE
  long Nlayer;              // number of layers in the profile
  double r0;                // r0 at wfs lambda
  double *cn2;              // profile strengh, units as in  E-SPE-ESO-276-0206_atmosphericparameters
  double *h;                // altitude of layers (meters)
  double *L0;               // outer scale (meters)

  double rmax;              // maximum distance between subapertures (computed with yorick)
  double *tracking;         // telescope tracking error parameters (x^2, y^2 and xy), units : arcsec^2

  double pasDPHI;           //Precision of DPHI precomputation.
  int ncpu;                 //Number of CPU used (only with openMP)
  int part;                 //Computed part of the cov. matrix. 0: complete 1: cmm 2: cpp 3: cpm
  int Nx;
  int Nslopes;
  int nlgs;

  double lgs_cst;
  double noise_var;
  double spot_width;
  double lgs_alt;
  double lgs_depth;
};

void initialize_tomo2(struct tomo_struct *tomo, int nssp);
void init_tomo_sys(struct tomo_struct *tomo);
void init_tomo_atm(struct tomo_struct *tomo);
void print_tomo(struct tomo_struct tomo);
void fill_Cmaa(struct tomo_struct tomo, double *Cmaa);
void fill_Caa(struct tomo_struct tomo, double *Caa);
void fill_Caa_dam(struct tomo_struct tomo, double *Caa);
void fill_Cmaa_dam(struct tomo_struct tomo, double *Caa);
void free_tomo(struct tomo_struct *tomo);

void matcov(struct tomo_struct tomo, double *data);
void matcov_dam(struct tomo_struct tomo, double *data);
void matcov_cpp(struct tomo_struct tomo, double *data);

void generateXY(struct tomo_struct *tomo);
void subap_position(struct tomo_struct tomo,  double ***u, double ***v );
void subap_position_dam(struct tomo_struct tomo,  double ***u, double ***v );

double cov_XX(double du, double dv, double ac, double ad, double bc, double bd, double *rr, double **tabDPHI, long indexL0, double convert);
double cov_YY(double du, double dv, double ac, double ad, double bc, double bd, double *rr, double **tabDPHI, long indexL0, double convert);
double cov_XY(double du, double dv, double s0, double *rr, double **tabDPHI, long indexL0, double convert);
double* compute_cov(double du, double dv,double ac, double ad, double bc, double bd,double s1, double s2, double *rr, double **tabDPHI, long indexL0, double convert, double pasDPHI, double units);

double DPHI(double x, double y, long indexL0, double *rr, double **tabDPHI, double convert);
double **tabulateDPHI(double* rr,struct tomo_struct tomo, long Ndphi, long *indexL0);

double **arr2dAlloc(long nbLin, long nbCol);
void arr2dFree(double **tableau);
double*** arr3dAlloc(long Nw, long *Nsubap, long Nlayer);
void arr3dFree(double ***array, long Nw, long *Nsubap);

#endif

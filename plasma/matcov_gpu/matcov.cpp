#include"matcov.h"
#include <algorithm> 

#ifdef USE_OPENMP
#include<omp.h>
#endif

void print_tomo(struct tomo_struct tomo){
  int i;
  printf("tomo->DiamTel = %f\n",tomo.DiamTel);
  printf("tomo->obs = %f\n",tomo.obs);

  printf("tomo->Nw = %lu\n",tomo.Nw);
  for(i=0; i<tomo.Nw; i++){
    printf("tomo->Nssp[%d] = %lu\n",i, tomo.Nssp[i]);
    printf("tomo->Nsubap[%d] = %lu\n",i, tomo.Nsubap[i]);
    printf("tomo->GsAlt[%d] = %f\n", i, tomo.GsAlt[i]);
    printf("tomo->type[%d] = %d\n", i, tomo.type[i]);
    printf("tomo->alphaX[%d] = %f\n", i, tomo.alphaX[i]);
    printf("tomo->alphaY[%d] = %f\n", i, tomo.alphaY[i]);
    printf("tomo->XPup[%d] = %f\n", i, tomo.XPup[i]);
    printf("tomo->YPup[%d] = %f\n", i, tomo.YPup[i]);
    printf("tomo->thetaML[%d] = %f\n", i, tomo.thetaML[i]);
    printf("tomo->thetaCam[%d] = %f\n", i, tomo.thetaCam[i]);
    printf("tomo->sensibilite[%d] L0= %f\n", i, tomo.sensibilite[i]);
    printf("tomo->diamPup[%d] = %f\n", i, tomo.diamPup[i]);
    printf("tomo->sspSize[%d] = %f\n", i, tomo.sspSize[i]);

  }
  printf("tomo->Nlayer = %lu\n", tomo.Nlayer);
  for(i=0; i<tomo.Nlayer; i++){
    printf("tomo->cn2[%d] = %f\n", i, tomo.cn2[i]);
    printf("tomo->h[%d] = %f\n", i, tomo.h[i]);
    printf("tomo->L0[%d] = %f\n", i, tomo.L0[i]);
  }
  printf("tomo->rmax = %f\n", tomo.rmax);

  printf("tomo->tracking[0] = %f\n", tomo.tracking[0]);
  printf("tomo->tracking[1] = %f\n", tomo.tracking[1]);
  printf("tomo->tracking[2] = %f\n", tomo.tracking[2]);

  printf("tomo->ncpu = %d\n", tomo.ncpu);

  printf("tomo->pasDPHI = %f\n", tomo.pasDPHI);
}

void read_paraml(FILE *file, long *par) {
  char line[128];
  fgets(line, sizeof(line), file);
  //printf("%s", line);
  fgets(line, sizeof(line), file);
  long tmpi;
  if( sscanf(line, "%ld", &tmpi) == 1 ) *par = tmpi;
  //printf("%ld\n", tmpi);
}

void read_parami(FILE *file, int *par) {
  char line[128];
  fgets(line, sizeof(line), file);
  //printf("%s", line);
  fgets(line, sizeof(line), file);
  int tmpi;
  if( sscanf(line, "%d", &tmpi) == 1 ) *par = tmpi;
  //printf("%d\n", tmpi);
}

void read_arrayl(FILE *file, long *arr, int nmax) {
  char line[128];
  char *tail = line;
  int i = 0;
  fgets(line, sizeof(line), file);
  //printf("%s", line);
  fgets(line, sizeof(line), file);
  while((*tail) && (i < nmax)){
    while (isspace (*tail)) tail++;
    arr[i++] = strtol(tail,&tail,10);
  }
  //for(i=0; i<nmax; i++) printf("%ld ",arr[i]);
  //printf("\n");
}

void read_arrayi(FILE *file, int *arr, int nmax) {
  char line[128];
  char *tail = line;
  int i = 0;
  fgets(line, sizeof(line), file);
  //printf("%s", line);
  fgets(line, sizeof(line), file);
  while((*tail) && (i < nmax)){
    while (isspace (*tail)) tail++;
    arr[i++] = (int)strtol(tail,&tail,10);
  }
  //for(i=0; i<nmax; i++) printf("%d ",arr[i]);
  //printf("\n");
}

void read_paramd(FILE *file, double *par) {
  char line[128];
  fgets(line, sizeof(line), file);
  //printf("%s", line);
  fgets(line, sizeof(line), file);
  float tmpf;
  if( sscanf(line, "%f", &tmpf) == 1 ) *par = tmpf;
  //printf("%f\n", tmpf);
}

void read_arrayd(FILE *file, double *arr, int nmax) {
  char line[128];
  char *tail = line;
  int i = 0;
  fgets(line, sizeof(line), file);
  //printf("%s", line);
  fgets(line, sizeof(line), file);
  while((*tail) && (i < nmax)){
    while (isspace (*tail)) tail++;
    arr[i++] = strtod(tail,&tail);
  }
  //for(i=0; i<nmax; i++) printf("%f ",arr[i]);
  //printf("\n");
}

void initialize_tomo2(struct tomo_struct *tomo, int nssp){
  printf("initialize 10 wfs %dx%d with 9 layers\n", nssp, nssp);
  int i;

  tomo->DiamTel = 40;
  tomo->obs = 0.2;

  tomo->Nw = 10; // number of wavefront sensors

  // array of the number of subap of each WFS along the telescop diameter, contains Nw elements
  tomo->Nssp = (long*)malloc(tomo->Nw*sizeof(long));

  // array of the number of subap of each WFS, contains Nw elements
  tomo->Nsubap = (long*)malloc(tomo->Nw*sizeof(long));

  // array of the inverse of the guide star altitude (in 1/meters), contains Nw elements
  tomo->GsAlt = (double*)malloc(tomo->Nw*sizeof(double));

  // type of WFS, 0, 1, 2 or 3. 0 is unused, 1=NGS, 2=LGS, 3=TipTilt-guide star
  tomo->type = (int*)malloc(tomo->Nw*sizeof(int));

  // Pointing directions of WFS
  double alphaX[] = { 0.000494509, 0.000247255, -0.000247255, -0.000494509,
      -0.000247255, 0.000247255, -7.75701e-05, -0.00015514, 0.000252103, 0 }; // pointing direction in X, arcseconds
  tomo->alphaX = (double*)malloc(tomo->Nw*sizeof(double));
  double alphaY[] = { 0, 0.000428258, 0.000428258, 6.05599e-20, -0.000428258,
      -0.000428258, -0.000174533, 0.000184229, -2.42407e-05, 0 }; // pointing direction in Y, arcseconds
  tomo->alphaY = (double*)malloc(tomo->Nw*sizeof(double));

  // Deviations of WFSs
  tomo->XPup = (double*)malloc(tomo->Nw*sizeof(double)); // pupil shift of the WFS, in meters
  tomo->YPup = (double*)malloc(tomo->Nw*sizeof(double)); // pupil shift of the WFS, in meters
  tomo->thetaML = (double*)malloc(tomo->Nw*sizeof(double));  // rotation of microlenses
  tomo->thetaCam = (double*)malloc(tomo->Nw*sizeof(double)); // rotation of microlenses
  tomo->sensibilite = (double*)malloc(tomo->Nw*sizeof(double)); // sensitivity coeff of this WFS

  tomo->diamPup = (double*)malloc(tomo->Nw*sizeof(double)); 
  tomo->sspSize = (double*)malloc(tomo->Nw*sizeof(double)); // subaperture size of this WFS

  for(i=0; i<tomo->Nw; i++){
    tomo->Nssp[i] = nssp;
    tomo->Nsubap[i] = -1;
    tomo->GsAlt[i] = (i<6? 1e-05:0);
    tomo->type[i] = (i<6? 2:1);

    tomo->alphaX[i] = alphaX[i];
    tomo->alphaY[i] = alphaY[i];
    tomo->XPup[i]=0;
    tomo->YPup[i]=0;
    tomo->thetaML[i]=0;
    tomo->thetaCam[i]=0;
    tomo->sensibilite[i]=1;

    tomo->diamPup[i]=nssp;
    tomo->sspSize[i]=tomo->DiamTel/tomo->diamPup[i];
  }

  tomo->Nlayer = 9; // number of layers in the profile

  // PROFILE
  double cn2[] = { 15.8791, 0.79031, 1.34961, 3.526, 3.00622, 0.896698, 1.81771,
      1.30705, 1.82379 }; // profile strengh, units TBD ..
  tomo->cn2 = (double*)malloc(tomo->Nlayer*sizeof(double));
  double h[] = { 47, 140, 281, 562, 1125, 2250, 4500, 9000, 18000 }; // altitude of layers (meters)
  tomo->h = (double*)malloc(tomo->Nlayer*sizeof(double));
  double L0[] = { 25, 25, 25, 25, 25, 25, 25, 25, 25 }; // outer scale (meters)
  tomo->L0 = (double*)malloc(tomo->Nlayer*sizeof(double));

  for(i=0; i<tomo->Nlayer; i++){
    tomo->cn2[i] = cn2[i];
    tomo->h[i] = h[i];
    tomo->L0[i] = L0[i];
  }

  tomo->rmax = 57.8023; // maximum distance between subapertures (computed with yorick)
 
  // telescope tracking error parameters (x^2, y^2 and xy), units : arcsec^2
  tomo->tracking = (double*)malloc(3*sizeof(double));
  tomo->tracking[0] = 0;
  tomo->tracking[1] = 0;
  tomo->tracking[2] = 0;

  tomo->ncpu = 12; //Number of CPU used (only with openMP)

  tomo->pasDPHI = 0.0001;

  //Generate the subapertures positions and fill tomo.Nsubap
  generateXY(tomo);
}

void init_tomo_sys(struct tomo_struct *tomo){
  int i;

  FILE *file = fopen("sys-params.txt", "r");
  read_paramd(file,&(tomo->DiamTel));
 
  read_paramd(file,&(tomo->obs));

  read_paraml(file,&(tomo->Nw)); // number of wavefront sensors

  // array of the number of subap of each WFS along the telescop diameter, contains Nw elements
  long nssp;
  read_paraml(file,&nssp); // number of wavefront sensors
  tomo->Nssp = (long*)malloc(tomo->Nw*sizeof(long));
  //read_arrayl(file,tomo->Nssp,tomo->Nw);

  // array of the number of subap of each WFS, contains Nw elements
  tomo->Nsubap = (long*)malloc(tomo->Nw*sizeof(long));

  // array of the inverse of the guide star altitude (in 1/meters), contains Nw elements
  tomo->GsAlt = (double*)malloc(tomo->Nw*sizeof(double));
  read_arrayd(file,tomo->GsAlt,tomo->Nw);

  // type of WFS, 0, 1, 2 or 3. 0 is unused, 1=NGS, 2=LGS, 3=TipTilt-guide star
  tomo->type = (int*)malloc(tomo->Nw*sizeof(int));
  read_arrayi(file,tomo->type,tomo->Nw);

  read_parami(file,&(tomo->nlgs));

  // Pointing directions of WFS
  tomo->alphaX = (double*)malloc(tomo->Nw*sizeof(double));
  read_arrayd(file,tomo->alphaX,tomo->Nw);
  tomo->alphaY = (double*)malloc(tomo->Nw*sizeof(double));
  read_arrayd(file,tomo->alphaY,tomo->Nw);

  // Deviations of WFSs
  tomo->XPup = (double*)malloc(tomo->Nw*sizeof(double)); // pupil shift of the WFS, in meters
  read_arrayd(file,tomo->XPup,tomo->Nw);
  tomo->YPup = (double*)malloc(tomo->Nw*sizeof(double)); // pupil shift of the WFS, in meters
  read_arrayd(file,tomo->YPup,tomo->Nw);
  tomo->thetaML = (double*)malloc(tomo->Nw*sizeof(double));  // rotation of microlenses
  read_arrayd(file,tomo->thetaML,tomo->Nw);
  tomo->thetaCam = (double*)malloc(tomo->Nw*sizeof(double)); // rotation of microlenses
  read_arrayd(file,tomo->thetaCam,tomo->Nw);
  tomo->sensibilite = (double*)malloc(tomo->Nw*sizeof(double)); // sensitivity coeff of this WFS
  read_arrayd(file,tomo->sensibilite,tomo->Nw);

  tomo->diamPup = (double*)malloc(tomo->Nw*sizeof(double)); 
  tomo->sspSize = (double*)malloc(tomo->Nw*sizeof(double)); // subaperture size of this WFS

  for(i=0; i<tomo->Nw; i++){
    tomo->Nssp[i]    = nssp; // all wfs have same number of subaps
    tomo->alphaX[i] /= 206265.0; // convert to radian
    tomo->alphaY[i] /= 206265.0; 
    tomo->diamPup[i] =(double)tomo->Nssp[i];
    tomo->sspSize[i] =(double)tomo->DiamTel/tomo->diamPup[i];
  }

  // telescope tracking error parameters (x^2, y^2 and xy), units : arcsec^2
  tomo->tracking = (double*)malloc(3*sizeof(double));
  read_arrayd(file,tomo->tracking,3);

  read_paramd(file,&(tomo->pasDPHI));

  read_parami(file,&(tomo->ncpu));
  /*
      noise stuff
   */

  read_paramd(file,&(tomo->lgs_cst));

  read_paramd(file,&(tomo->noise_var)); 

  read_paramd(file,&(tomo->spot_width)); 

  read_paramd(file,&(tomo->lgs_alt));

  read_paramd(file,&(tomo->lgs_depth));

  //Generate the subapertures positions and fill tomo.Nsubap
  generateXY(tomo);

  int Nslopes = 0;
  for (int i = 0; i < tomo->Nw; i++)
    Nslopes += tomo->Nsubap[i]*2;
  tomo->Nslopes = Nslopes;

  fclose(file);
  //printf("initialized %lu wfs \n", tomo->Nw);

}

void init_tomo_atm(struct tomo_struct *tomo){
  int i;

  FILE *file = fopen("atm-params.txt", "r");

  read_paraml(file,&(tomo->Nlayer)); // number of layers in the profile
  read_paramd(file,&(tomo->r0)); // number of layers in the profile
 
  // PROFILE
  tomo->cn2 = (double*)malloc(tomo->Nlayer*sizeof(double)); // profile strengh, units TBD ..
  read_arrayd(file,tomo->cn2,tomo->Nlayer);
  tomo->h = (double*)malloc(tomo->Nlayer*sizeof(double)); // altitude of layers (meters)
  read_arrayd(file,tomo->h,tomo->Nlayer);
  tomo->L0 = (double*)malloc(tomo->Nlayer*sizeof(double)); // outer scale (meters)
  read_arrayd(file,tomo->L0,tomo->Nlayer);
  /*
  double scn2 = 0.0d;
  for (int cc=0;cc<tomo->Nw;cc++) {
    scn2 += tomo->cn2[cc];
  }
  for (int cc=0;cc<tomo->Nw;cc++) {
    tomo->cn2[cc] /= scn2;
    tomo->cn2[cc] /= pow(tomo->r0,5.0/3.0);
  }
  */
  //rmax = max(abs(tomo.wfs.x,tomo.wfs.y))*2*max(tomo.learn.altitude)/206265.+tomo.tel.diam;
  double dmax = 0.0d;
  double maxalt = tomo->h[tomo->Nlayer-1];
  for (int cc=0;cc<tomo->Nw;cc++) {
    double tmp = sqrtf(tomo->alphaX[cc]*tomo->alphaX[cc] + tomo->alphaY[cc]*tomo->alphaY[cc]);
    if (tmp > dmax) dmax = tmp;
  }
  tomo->rmax = dmax * 2 * maxalt + tomo->DiamTel;

  fclose(file);
  //printf("initialized %lu layers \n", tomo->Nlayer);
}

void fill_Caa(struct tomo_struct tomo, double *Caa){
  tomo.part = 1; //compute the Caa matrix
  //print_tomo(tomo);
  matcov(tomo, Caa);
}

void fill_Caa_dam(struct tomo_struct tomo, double *Caa){
  tomo.part = 1; //compute the Caa matrix
  //print_tomo(tomo);
  matcov_dam(tomo, Caa);
}

void fill_Cmaa(struct tomo_struct tomo, double *Cmaa){
  tomo.part = 3; //compute the Cmaa matrix
  matcov(tomo, Cmaa);
}

void fill_Cmaa_dam(struct tomo_struct tomo, double *Cmaa){
  tomo.part = 3; //compute the Cmaa matrix
  matcov_dam(tomo, Cmaa);
}

void free_tomo(struct tomo_struct *tomo){
  free(tomo->Nssp);

  // array of the number of subap of each WFS, contains Nw elements
  free(tomo->Nsubap);

  // array of the inverse of the guide star altitude (in 1/meters), contains Nw elements
  free(tomo->GsAlt);

  // type of WFS, 0, 1, 2 or 3. 0 is unused, 1=NGS, 2=LGS, 3=TipTilt-guide star
  free(tomo->type);

  free(tomo->alphaX);
  free(tomo->alphaY);

  // Deviations of WFSs
  free(tomo->XPup); // pupil shift of the WFS, in meters
  free(tomo->YPup); // pupil shift of the WFS, in meters
  free(tomo->thetaML);  // rotation of microlenses
  free( tomo->thetaCam); // rotation of microlenses
  free(tomo->sensibilite); // sensitivity coeff of this WFS

  free(tomo->diamPup); 
  free(tomo->sspSize);

  free(tomo->cn2);
  free(tomo->h);
  free(tomo->L0);
  free(tomo->tracking);

  free(tomo->X);
  free(tomo->Y);
}

/* %%%%%%%%%%%%%%%%%%%%%% Matcov %%%%%%%%%%%%%%%%%%%%%%%% */
void
matcov(struct tomo_struct tomo, double *data) {
  
  // %%%%%%% Pre-computation of DPHI %%%%%%%%%%
  //Computes an array of DPHI (tabDPHI) for an array of subaperture distance rr for each DIFFERENT L0
  const long cNw = tomo.Nw;
  const long cNlayer = tomo.Nlayer;
  const double crmax = tomo.rmax;
  const double pasDPHI = 1./tomo.pasDPHI; //inverse du pas de rr
  const long cNdphi = floor(crmax*pasDPHI)+1;
  const double convert = (double)(cNdphi-1)/(crmax+1./pasDPHI);
  double rr[cNdphi];
  long i; 

 //rr varie de 0 à rmax+pasDphi (évite de traiter à part les cas où rr =rmax dans la fonction de DPHI)
  for (i=0;i<cNdphi;i++) {
    rr[i] = (double)(i)/ convert;
  }
  
  long indexL0[cNlayer]; //link between index in L0 and index in L0diff
  double **tabDPHI;
  tabDPHI = tabulateDPHI(rr, tomo, cNdphi, indexL0);

  

  // %%%%%%% Computation of the sub-apertures positions and sizes %%%%%%%%%%%
 // u, v :arrays containing all the sub-apertures coordinates of all WFS, one after the other
  // u[0][1][3] is the X-coordinate of subap number 3 of wfs number 0 at altitude 3
  double*** u = arr3dAlloc( cNw, tomo.Nsubap, cNlayer);
  double*** v = arr3dAlloc( cNw, tomo.Nsubap, cNlayer);

  //Computes  u and v
   subap_position(tomo, u, v);
  
  //Rescale the projected size of all subapertures at the different altitude layer
  double sspSizeL[cNw][cNlayer];
  long l, n;

  for (n = 0; n < cNw; n++) {
    for (l = 0; l < cNlayer; l++) {
      sspSizeL[n][l] = tomo.sspSize[n] * (1. - tomo.GsAlt[n] * tomo.h[l]);
    }
  }

  
  // %%%%%%% Computation of the covariance matrix %%%%%%%%%%%
  const double lambda2 = 0.00026942094446267851; //lambda2 = pow(206265.*0.5e-6/2./3.1415926535,2);
  long m, j;
  long ioff = 0, joff = 0;
  double units[cNlayer];
  long Nslopes=0; // very total number of slopes


 for (i = 0; i < cNw; i++)
    Nslopes += tomo.Nsubap[i]*2;
  
  //  Initialization and limits of the loops on WFS, function of the computed parts of the matrix
  long m0=0;
  long mf=0;
  long n0=0;
  long *nf=0;
  MKL_INT NL;
  long ts = cNw - 1;//Truth sensor : ts

  if (tomo.part == 0) { //Complete matrix
    m0 = 0;
    mf = cNw;
    n0 = 0.;
    nf = &m;

    NL=0; // very total number of slopes
    for (i = 0; i < cNw; i++)
      NL += tomo.Nsubap[i]*2;

  } else if (tomo.part == 1) { //Cmm
    m0 = 0;
    mf = cNw-1;
    n0 = 0; 
    nf = &m;

    Nslopes=0; // very total number of slopes
    for (i = 0; i < cNw; i++)
      Nslopes += tomo.Nsubap[i];
    Nslopes *= 2;
    NL = Nslopes-2.*tomo.Nsubap[ts];

  } else if  (tomo.part == 3) { //Cpm
    m0 = 0;
    mf = cNw - 1;
    n0 = cNw - 1;
    nf = &ts;
    NL = 2.*tomo.Nsubap[cNw - 1];
  }

  
  //WFS m
  for (m = m0; m < mf; m++) {

    const long Ni = tomo.Nsubap[m] + ioff;
    
    //WFS n
    for (n = n0; n < *nf+1; n++) {

      const long off_XY = tomo.Nsubap[n];
      const long off_YX = tomo.Nsubap[m] * NL;
      const long off_YY = off_XY+ off_YX;

      const long Nj = tomo.Nsubap[n] + joff;
      const double kk = 1. / (tomo.sspSize[m] * tomo.sspSize[n]);

      for (l = 0; l < cNlayer; l++) {
        units[l] = kk * lambda2 * tomo.cn2[l];
      }

#ifdef USE_OPENMP
#pragma omp parallel private(j,l) num_threads(tomo.ncpu)
#pragma omp for nowait
#endif
      //Subaperture i
      for (i = ioff; i < Ni; i++) {

        //Subaperture j
        for (j = joff; j < Nj; j++) {
	  double caa_xx = 0;
	  double caa_yy = 0;
	  double caa_xy = 0;

	  //Layer l
	   for (l = 0; l < cNlayer; l++) {

	    //test if the altitude layers is not higher than the LGS altitude
	    if ((sspSizeL[m][l] > 0) && (sspSizeL[n][l] > 0)) {
	      //Distances in x and y between the subapertures i and j
	      const double du = u[m][i-ioff][l] - u[n][j-joff][l];	      
	      const double dv =  v[m][i-ioff][l] - v[n][j-joff][l];

	      const double s1 = sspSizeL[m][l] * 0.5;
	      const double s2 = sspSizeL[n][l] * 0.5;

	      const double ac = s1 - s2;
	      const double ad = s1 + s2;
	      const double bc = -ad;   // initially -s1-s2;
	      const double bd = -ac;   // initially -s1+s2;

	      //Computation of the covariance on each layer
	      double *cov;

	      cov = compute_cov(du, dv, ac, ad, bc, bd, s1, s2,rr,tabDPHI,indexL0[l],convert,pasDPHI, units[l]);

	      caa_xx += cov[0];	    
	      caa_yy += cov[1];
	      caa_xy += cov[2];

	      free(cov);
	      }
	  }

	  const long i0 = i * NL + j;
	  data[i0] = caa_xx;   //xx
	  data[i0 + off_XY] = caa_xy;   //xy
	  data[i0 + off_YX] = caa_xy;   //yx
	  data[i0 + off_YY] = caa_yy; //yy

	}
      }
      joff = joff + 2 * tomo.Nsubap[n];
    }
    ioff = ioff + 2 * tomo.Nsubap[m];
    joff = 0;
  }
  
  //Recopie de la symétrie
  if (tomo.part == 0 || tomo.part == 1) { //Complete matrix
    MKL_INT size=NL-1;
    MKL_INT one=1;
    double *matL = (double*)data+1;
    double *matU = (double*)data+NL;
    do{
      //dcopy(&size, matU, &NL, matL, &one);
      //-> if there are issues with dcopy use this:
      for(int j=0; j<size; j++) matL[j]=matU[j*NL];
      size--;
      matL+=NL+1;
      matU+=NL+1;
    } while (size>0);
   }
  
  arr3dFree(u, cNw, tomo.Nsubap);
  arr3dFree(v, cNw, tomo.Nsubap);
  arr2dFree(tabDPHI);
}


/* %%%%%%%%%%%%%%%%%%%%%% Matcov %%%%%%%%%%%%%%%%%%%%%%%% */
void
matcov_dam(struct tomo_struct tomo, double *data) {
  
  // %%%%%%% Pre-computation of DPHI %%%%%%%%%%
  //Computes an array of DPHI (tabDPHI) for an array of subaperture distance rr for each DIFFERENT L0
  const long cNw = tomo.Nw;
  const long cNlayer = tomo.Nlayer;
  const double crmax = tomo.rmax;
  const double pasDPHI = 1./tomo.pasDPHI; //inverse du pas de rr
  const long cNdphi = floor(crmax*pasDPHI)+1;
  const double convert = (double)(cNdphi-1)/(crmax+1./pasDPHI);
  double rr[cNdphi];
  long i; 

 //rr varie de 0 à rmax+pasDphi (évite de traiter à part les cas où rr =rmax dans la fonction de DPHI)
  for (i=0;i<cNdphi;i++) {
    rr[i] = (double)(i)/ convert;
  }
  
  long indexL0[cNlayer]; //link between index in L0 and index in L0diff
  double **tabDPHI;
  tabDPHI = tabulateDPHI(rr, tomo, cNdphi, indexL0);

  // %%%%%%% Computation of the sub-apertures positions and sizes %%%%%%%%%%%
 // u, v :arrays containing all the sub-apertures coordinates of all WFS, one after the other
  // u[0][1][3] is the X-coordinate of subap number 3 of wfs number 0 at altitude 3
  double*** u = arr3dAlloc( cNw, tomo.Nsubap, cNlayer);
  double*** v = arr3dAlloc( cNw, tomo.Nsubap, cNlayer);

  //Computes  u and v
   subap_position_dam(tomo, u, v);
  
  //Rescale the projected size of all subapertures at the different altitude layer
  double sspSizeL[cNw][cNlayer];
  long l, n;

  for (n = 0; n < cNw; n++) {
    for (l = 0; l < cNlayer; l++) {
      sspSizeL[n][l] = tomo.sspSize[n] * (1. - tomo.GsAlt[n] * tomo.h[l]);
    }
  }

  /*
  double *sspSizeL2 = (double *)malloc(sizeof(double)*cNw*cNlayer);
  for (int cc = 0; cc < cNw * cNlayer; cc++) {
    n = cc / cNw;
    l = cc - n * cNw;
    sspSizeL2[cc] = tomo.sspSize[n] * (1. - tomo.GsAlt[n] * tomo.h[l]);
    printf("%f %f\n",sspSizeL2[cc],sspSizeL[n][l]);
  }
  free(sspSizeL2);
  */

  // %%%%%%% Computation of the covariance matrix %%%%%%%%%%%
  const double lambda2 = 0.00026942094446267851; //lambda2 = pow(206265.*0.5e-6/2./3.1415926535,2);
  long m, j;
  long ioff = 0, joff = 0;
  double units[cNlayer];
  long Nslopes=0; // very total number of slopes


 for (i = 0; i < cNw; i++)
    Nslopes += tomo.Nsubap[i]*2;
  
  //  Initialization and limits of the loops on WFS, function of the computed parts of the matrix
  long m0=0;
  long mf=0;
  long n0=0;
  long *nf=0;
  MKL_INT NL;
  long ts = cNw - 1;//Truth sensor : ts

  if (tomo.part == 0) { //Complete matrix
    m0 = 0;
    mf = cNw;
    n0 = 0.;
    nf = &m;

    NL=0; // very total number of slopes
    for (i = 0; i < cNw; i++)
      NL += tomo.Nsubap[i]*2;
    Nslopes = NL;

  } else if (tomo.part == 1) { //Cmm
    m0 = 0;
    mf = cNw-1;
    n0 = 0; 
    nf = &m;

    Nslopes=0; // very total number of slopes
    for (i = 0; i < cNw; i++)
      Nslopes += tomo.Nsubap[i];
    Nslopes *= 2;
    NL = Nslopes-2.*tomo.Nsubap[ts];
    Nslopes = NL;

  } else if  (tomo.part == 3) { //Cpm
    m0 = 0;
    mf = cNw - 1;
    n0 = cNw - 1;
    nf = &ts;
    Nslopes=0; // very total number of slopes
    for (i = 0; i < cNw; i++)
      Nslopes += tomo.Nsubap[i];
    Nslopes -= tomo.Nsubap[ts];
    Nslopes *= 2;
    NL = 2.*tomo.Nsubap[ts];
  }

  int *tab_wfs;
  tab_wfs = (int*)malloc(Nslopes*sizeof(int));
  int *tab_subap;
  tab_subap = (int*)malloc(Nslopes*sizeof(int));
  int *tab_xy;
  tab_xy = (int*)malloc(Nslopes*sizeof(int));
 
  int cpt = 0;
  for (int cc=0;cc<cNw;cc++) {
    if (cc != ts) {
      int nslps = tomo.Nsubap[cc]*2;
      for (int ccc=0;ccc<nslps;ccc++) {
	if (cc > ts) tab_wfs[ccc+cpt] = cc - 1;
	else tab_wfs[ccc+cpt] = cc;
	if (ccc < nslps/2) {
	  tab_subap[ccc+cpt] = ccc;
	  tab_xy[ccc+cpt] = 0;
	} else {
	  tab_subap[ccc+cpt] = ccc - nslps/2;
	  tab_xy[ccc+cpt] = 1;
	}
      }
      cpt += nslps;
    }
  }
  
  long mcc=0;
  for (mcc = 0;mcc<NL*Nslopes;mcc++) {
    int jpos = mcc / Nslopes;
    int ipos = mcc - jpos * Nslopes;

  //WFS m
    if  (tomo.part == 3) m = ts;
    else m = tab_wfs[jpos];
  //WFS n
    n = tab_wfs[ipos];
  //subap i
    if  (tomo.part == 3) i = (jpos < tomo.Nsubap[ts]) ? jpos : jpos - tomo.Nsubap[ts];
    else i = tab_subap[jpos];
  //subap j
    j = tab_subap[ipos];
  //xy i
    int xy_i;
    if  (tomo.part == 3) xy_i = (jpos < tomo.Nsubap[ts]) ? 0 : 1;
    else xy_i = tab_xy[jpos];
  //xy j
    int xy_j = tab_xy[ipos];

    const double kk = 1. / (tomo.sspSize[m] * tomo.sspSize[n]);

    for (l = 0; l < cNlayer; l++) {
      units[l] = kk * lambda2 * tomo.cn2[l];
    }

   int type;
    if ((xy_i == 0) && (xy_j == 0)) type = 0;
    if ((xy_i == 0) && (xy_j == 1)) type = 1;
    if ((xy_i == 1) && (xy_j == 0)) type = 2;
    if ((xy_i == 1) && (xy_j == 1)) type = 3;

    //Layer l
    for (l = 0; l < cNlayer; l++) {
      //test if the altitude layers is not higher than the LGS altitude
      if ((sspSizeL[m][l] > 0) && (sspSizeL[n][l] > 0)) {
	const double du = u[m][i][l] - u[n][j][l];	      
	const double dv =  v[m][i][l] - v[n][j][l];
	
	const double s1 = sspSizeL[m][l] * 0.5;
	const double s2 = sspSizeL[n][l] * 0.5;
	
	const double ac = s1 - s2;
	const double ad = s1 + s2;
	const double bc = -ad;   // initially -s1-s2;
	const double bd = -ac;   // initially -s1+s2;
	if (type == 0) data[mcc] += 0.5 * pasDPHI * cov_XX(du,dv,ac,ad,bc,bd,rr,tabDPHI,indexL0[l],convert) * units[l];
	if (type == 3) data[mcc] += 0.5 * pasDPHI * cov_YY(du,dv,ac,ad,bc,bd,rr,tabDPHI,indexL0[l],convert) * units[l];
	if ((type == 1) || (type == 2)) {
	  const double s0 = sqrt(s1 * s1 + s2 * s2); //half size of the subaperture equivalent to a convolution by s1 and s2
	  const double dd = (s1 > s2) ? 1. - s2 / s1 : 1. - s1 / s2; // Nono's style ....
	  data[mcc] += 0.25 * pasDPHI * cov_XY(du,dv,s0,rr,tabDPHI,indexL0[l],convert) * units[l] * (1. - dd * dd);
	}
      }
    }
  }
  free(tab_wfs);
  free(tab_subap);
  free(tab_xy);

  arr3dFree(u, cNw, tomo.Nsubap);
  arr3dFree(v, cNw, tomo.Nsubap);
  arr2dFree(tabDPHI);
}




/* %%%%%%%%%%%%%%%%%%%%%% CPP %%%%%%%%%%%%%%%%%%%%%%%% */


void
matcov_cpp(struct tomo_struct tomo, double *data)
/* DOCUMENT matcov_cpp(struct tomo_struct tomo, double *data) 
Compute the covariance matrix of the TS. The TS must be at the last WFS.

   <tomo>                :  structure with all the needed information ( see matcov.h). Can contain only the TS or all the WFS. In last case, the TS must be the last WFS.

 SEE ALSO:
*/
 {
  // %%%%%%% Pre-computation of DPHI %%%%%%%%%%
  //Computes an array of DPHI (tabDPHI) for an array of subaperture distance rr for each DIFFERENT L0
  const long cNw = tomo.Nw;
  const long cNlayer = tomo.Nlayer;
  const double crmax = tomo.rmax;
  const double pasDPHI = 1./tomo.pasDPHI; //inverse du pas de rr
  const long cNdphi = floor(crmax*pasDPHI)+1;
  const double convert = (double)(cNdphi-1)/(crmax+1./pasDPHI);
  double rr[cNdphi];
  long i; 

 //rr varie de 0 à rmax+pasDphi (évite de traiter à part les cas où rr =rmax dans la fonction de DPHI)
  for (i=0;i<cNdphi;i++) {
    rr[i] = (double)(i)/ convert;
  }

  long indexL0[cNlayer]; //link between index in L0 and index in L0diff
  double **tabDPHI;
  tabDPHI = tabulateDPHI(rr, tomo, cNdphi, indexL0);



  // %%%%%%% Computation of the covariance matrix %%%%%%%%%%%
  long l,j;
  const long ts = cNw - 1;
  const long cNsubapTS = tomo.Nsubap[ts];
  const long offXY = cNsubapTS* 2.*cNsubapTS;
  const long offYY = offXY+ cNsubapTS;
  const double lambda2 = 0.00026942094446267851; //lambda2 = pow(206265.*0.5e-6/2./3.1415926535,2);
  const double kk = 1. / (tomo.sspSize[ts] * tomo.sspSize[ts]);


  double units[cNlayer];
  for (l = 0; l < cNlayer; l++) {
    units[l] = kk * lambda2 * tomo.cn2[l];
  }


  long Nslopes=0; // very total number of slopes
  for (i = 0; i < cNw; i++)
    Nslopes += tomo.Nsubap[i];
  const long off = Nslopes - cNsubapTS;
  Nslopes *= 2;


  const double s = tomo.sspSize[ts] * 0.5;
  const double ac = 0.;
  const double ad = 2. * s;
  const double bc = -ad;   
  const double bd = 0.;  


#ifdef USE_OPENMP
#pragma omp parallel private(j,l) num_threads(tomo.ncpu)
#pragma omp for nowait
#endif
  //Subaperture i
  for (i = 0; i < cNsubapTS; i++) {

    //Subaperture j
    for (j = 0; j < cNsubapTS; j++) {
      double caa_xx = 0.;
      double caa_yy = 0.;
      double caa_xy = 0.;

      //Distances in x and y between the subapertures i and j
      const double du = tomo.X[i+off] - tomo.X[j+off];
      const double dv = tomo.Y[i+off] - tomo.Y[j+off];

      //Layer l
      for (l = 0; l < cNlayer; l++) {
        //test if the altitude layers is not higher than the LGS altitude
        if ((s > 0)) {
	  double *cov;
          //Computation of the covariance on each layer
	  cov = compute_cov(du, dv, ac, ad, bc, bd, s, s,rr,tabDPHI,indexL0[l],convert,pasDPHI, units[l]);
	  caa_xx += cov[0];	    
	  caa_yy += cov[1];
	  caa_xy += cov[2];
	  free(cov);
        }
      }

      const long i0 = i * 2.*cNsubapTS + j;

      data[i0] = caa_xx;          //xx

      data[i0 + cNsubapTS] =  caa_xy;          //xy

      data[i0 +  offXY] = caa_xy;          //yx
      
      data[i0 + offYY] = caa_yy; //yy
    }
  }
  arr2dFree(tabDPHI);
}

void generateXY(struct tomo_struct *tomo)
/* DOCUMENT  generateXY(struct tomo_struct tomo, double *Nsubap)
 <tomo>               :  structure with all the needed information
 <tomo.X> & <tomo.Y>            :   arrays containing all the sub-apertures coordinates of all WFS, one after the other
<tomo.Nsubap>              :  number of subaperture of ezach WFS
 Generate the position (X,Y) of each subapertures of each WFS on the telescope pupil and the number of subaperture of ezach WFS (Nsubap) 
 */
{
  const double bornemin = -tomo->DiamTel / 2.;
  const double Rtel2 = (tomo->DiamTel * tomo->DiamTel) / 4.;
  long NsubapTot = 0;
  long n;

  //Total number of subapertures (without obstruction)
  for (n = 0; n < tomo->Nw; n++) {
    NsubapTot += tomo->Nssp[n] * tomo->Nssp[n];
  }

  const long cNsubapTot = NsubapTot;
  double x[cNsubapTot], y[cNsubapTot];
  int index[cNsubapTot];

  int cpt = 0;
  int ioff = 0;

  //Computation of all the subapertures' positions
  for (n = 0; n < tomo->Nw; n++) {
    long Nsap = 0;
    double pas = tomo->DiamTel / (1. * tomo->Nssp[n]);
    int i;
    double Robs2;

    // to avoid some bug that eliminates useful central subapertures when obs=0.286
    if (tomo->Nssp[n] != 7 || (tomo->obs <= 0.285 || tomo->obs >= 0.29)) {
      Robs2 = tomo->DiamTel * tomo->obs / 2. * tomo->DiamTel * tomo->obs / 2.;
    } else {
      Robs2 = tomo->DiamTel * 0.285 / 2. * tomo->DiamTel * 0.285 / 2.;
    }

    if (tomo->Nssp[n] != 1) {
      for (i = 0; i < tomo->Nssp[n]; i++) {
        double tp = bornemin + pas / 2. * (2. * i + 1.); // y-coord of current subap
        int j;

        for (j = 0; j < tomo->Nssp[n]; j++) {
          x[ioff + j] = bornemin + pas / 2. * (2. * j + 1.); // x-coord of current subap
          y[ioff + j] = tp;

          double r2 = x[ioff + j] * x[ioff + j] + y[ioff + j] * y[ioff + j];

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

  tomo->X=(double*)malloc((cNsubapTot-cpt)*sizeof(double));
  tomo->Y=(double*)malloc((cNsubapTot-cpt)*sizeof(double));
  tomo->Nx = cNsubapTot-cpt;
  
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

void
subap_position(struct tomo_struct tomo, double ***u, double ***v) {
  /* DOCUMENT         subap_position(tomo, u, v)
   <tomo>                : structure with all the needed information.
   <u> and <v>           : 3d arrays containing the sub-apertures projected coordinates onto all the layers. u[0][2][1] is the X-coordinate of the subap 2 of the WFS 0 on the layer 1.

   Computes the projected coordinates of all subapertures  projected onto all the layer
   */
  long i;
  long n;
  long l;
  const double rad = 3.14159265358979323846 / 180.;

  for (l = 0; l < tomo.Nlayer; l++) {
    long ioff = 0;

    for (n = 0; n < tomo.Nw; n++) {

      const double dX = tomo.alphaX[n] * tomo.h[l];
      const double dY = tomo.alphaY[n] * tomo.h[l];

      const double rr = 1. - tomo.h[l] * tomo.GsAlt[n];

      const long Nsap = tomo.Nsubap[n];
      const long nssp = tomo.Nssp[n];

      //magnification factor
      const double G = tomo.diamPup[n] / (double) (nssp);

      //rotation angle
      const double th = tomo.thetaML[n] * rad;

      for (i = 0; i < Nsap; i++) {
        //taking magnification factor into account
        const double xtp = tomo.X[ioff + i] * G;
        const double ytp = tomo.Y[ioff + i] * G;

        //taking rotation into account
        double uu = xtp * cos(th) - ytp * sin(th);
        double vv = xtp * sin(th) + ytp * cos(th);

        //taking pupil offset into account
        uu += tomo.XPup[n];
        vv += tomo.YPup[n];

        //Projection onto  the layer
        u[n][i][l] = uu * rr + dX;
        v[n][i][l] = vv * rr + dY;
      }
      //index offset
      ioff += Nsap;
    }
  }
}

void
subap_position_dam(struct tomo_struct tomo, double ***u, double ***v) {
  /* DOCUMENT         subap_position(tomo, u, v)
   <tomo>                : structure with all the needed information.
   <u> and <v>           : 3d arrays containing the sub-apertures projected coordinates onto all the layers. u[0][2][1] is the X-coordinate of the subap 2 of the WFS 0 on the layer 1.

   Computes the projected coordinates of all subapertures  projected onto all the layer
   */
  long i;
  long n;
  long l;
  const double rad = 3.14159265358979323846 / 180.;
  long ioff[tomo.Nw];
  ioff[0] = 0;
  for (int i=1;i<tomo.Nw;i++) ioff[i] = ioff[i-1] + tomo.Nsubap[i-1];

  for (int cc = 0;cc < tomo.Nlayer * tomo.Nw * tomo.Nsubap[0]; cc++) {
    l = cc / (tomo.Nw * tomo.Nsubap[0]);

    const int pos = cc - l * (tomo.Nsubap[0] * tomo.Nw);

    i = pos / tomo.Nw;
    n = pos - i * tomo.Nw;

    const double dX = tomo.alphaX[n] * tomo.h[l];
    const double dY = tomo.alphaY[n] * tomo.h[l];
    
    const double rr = 1. - tomo.h[l] * tomo.GsAlt[n];
    
    const long nssp = tomo.Nssp[n];
    
    //magnification factor
    const double G = tomo.diamPup[n] / (double) (nssp);
    
    //rotation angle
    const double th = tomo.thetaML[n] * rad;

    //taking magnification factor into account
    const double xtp = tomo.X[ioff[n] + i] * G;
    const double ytp = tomo.Y[ioff[n] + i] * G;
    
    //taking rotation into account
    double uu = xtp * cos(th) - ytp * sin(th);
    double vv = xtp * sin(th) + ytp * cos(th);
    
    //taking pupil offset into account
    uu += tomo.XPup[n];
    vv += tomo.YPup[n];
    
    //Projection onto  the layer
    u[n][i][l] = uu * rr + dX;
    v[n][i][l] = vv * rr + dY;
  }

}

/*
 *  _____           _   _                       _        ____
 * |  ___|__  _ __ | |_(_) ___  _ __  ___    __| | ___  | __ )  __ _ ___  ___
 * | |_ / _ \| '_ \| __| |/ _ \| '_ \/ __|  / _` |/ _ \ |  _ \ / _` / __|/ _ \
 * |  _| (_) | | | | |_| | (_) | | | \__ \ | (_| |  __/ | |_) | (_| \__ \  __/
 * |_|  \___/|_| |_|\__|_|\___/|_| |_|___/  \__,_|\___| |____/ \__,_|___/\___|
 *
 */

double *compute_cov(double du, double dv, double ac, double ad, double bc, double bd,double s1, double s2, double *rr, double **tabDPHI, long indexL0, double convert, double pasDPHI, double units)
 /* DOCUMENT
   <du> & <dv>                 : X et Y coordinates of the distance between the deux considered subapertures.
   <ac> & <ad> & <bc> & <bd>  : precomputed values
   <s1> & <s2>                : half size of each subapertures
   <rr> & <tabDPHI>           : tabDPHI contains values of DPHI precomputed on rr
   <indexL0>                  : 
   <convert>                  :
   <pasDPHI>                  :
   <units>                    :
   Computes the XX, XY and YY covariance values for two subapertures.
   */
{        
  double *cov = (double*)malloc(3*sizeof(double));
  double cov_xx=0;
  double cov_yy=0;
  double cov_xy=0;


  //Computation of the covariance on each layer
  cov_xx = cov_XX(du,dv,ac,ad,bc,bd,rr,tabDPHI,indexL0,convert);
  cov_xx *= 0.5 * pasDPHI;
  
  cov_yy = cov_YY(du,dv,ac,ad,bc,bd,rr,tabDPHI,indexL0,convert);
  cov_yy *= 0.5 * pasDPHI;
  
  const double s0 = sqrt(s1 * s1 + s2 * s2); //half size of the subaperture equivalent to a convolution by s1 and s2
  cov_xy = cov_XY(du,dv,s0,rr,tabDPHI,indexL0,convert);
  cov_xy *= 0.25 * pasDPHI;
  
  // double cc = 1.-fmin(s1,s2)/fmax(s1,s2);
  // when s1==s2, then cc=0 and the computation is just void,
  // so we perform it only when s1!=s2
  const double cc = (s1 > s2) ? 1. - s2 / s1 : 1. - s1 / s2; // Nono's style ....
  cov_xy *= (1. - cc * cc);
  
  //units
  cov_xx *= units;
  cov_yy *= units;
  cov_xy *= units;

  cov[0] = cov_xx;
  cov[1] = cov_yy;
  cov[2] = cov_xy;
  return cov;
}



double cov_XX(double du, double dv, double ac, double ad, double bc, double bd, double *rr, double **tabDPHI, long indexL0, double convert)
 /* DOCUMENT
   Compute the XX-covariance with the distance sqrt(du2+dv2). DPHI is precomputed on tabDPHI.
 */
{
  return -DPHI(du + ac, dv, indexL0, rr, tabDPHI, convert)
    + DPHI(du + ad, dv, indexL0, rr, tabDPHI, convert)
    + DPHI(du + bc, dv, indexL0, rr, tabDPHI, convert)
    - DPHI(du + bd, dv, indexL0, rr, tabDPHI, convert);
}

double cov_YY(double du, double dv, double ac, double ad, double bc, double bd, double *rr, double **tabDPHI, long indexL0, double convert)
/* DOCUMENT
   Compute the YY-covariance with the distance sqrt(du2+dv2). DPHI is precomputed on tabDPHI.
 */
{ 
  return  -DPHI(du, dv + ac, indexL0, rr, tabDPHI, convert)
    + DPHI(du, dv + ad, indexL0, rr, tabDPHI, convert)
    + DPHI(du, dv + bc, indexL0, rr, tabDPHI, convert)
    - DPHI(du, dv + bd, indexL0, rr, tabDPHI, convert);
}


double cov_XY(double du, double dv, double s0, double *rr, double **tabDPHI, long indexL0, double convert)
/* DOCUMENT
   Compute the XY-covariance with the distance sqrt(du2+dv2). DPHI is precomputed on tabDPHI.
 */
{
  return -DPHI(du + s0, dv - s0, indexL0, rr, tabDPHI, convert)
    + DPHI(du + s0, dv + s0, indexL0, rr, tabDPHI, convert)
    + DPHI(du - s0, dv - s0, indexL0, rr, tabDPHI, convert)
    - DPHI(du - s0, dv + s0, indexL0, rr, tabDPHI, convert);
}


double
macdo_x56(double x, int k)
/* DOCUMENT  macdo_x56(x)

 Computation of the function
 f(x) = x^(5/6)*K_{5/6}(x)
 using a series for the esimation of K_{5/6}, taken from Rod Conan thesis :
 K_a(x)=1/2 \sum_{n=0}^\infty \frac{(-1)^n}{n!}
 \left(\Gamma(-n-a) (x/2)^{2n+a} + \Gamma(-n+a) (x/2)^{2n-a} \right) ,
 with a = 5/6.

 Setting x22 = (x/2)^2, setting uda = (1/2)^a, and multiplying by x^a,
 this becomes :
 x^a * Ka(x) = 0.5 $ -1^n / n! [ G(-n-a).uda x22^(n+a) + G(-n+a)/uda x22^n ]
 Then we use the following recurrence formulae on the following quantities :
 G(-(n+1)-a) = G(-n-a) / -a-n-1
 G(-(n+1)+a) = G(-n+a) /  a-n-1
 (n+1)! = n! * (n+1)
 x22^(n+1) = x22^n * x22
 and at each iteration on n, one will use the values already computed at step (n-1).
 The values of G(a) and G(-a) are hardcoded instead of being computed.

 The first term of the series has also been skipped, as it
 vanishes with another term in the expression of Dphi.

 SEE ALSO:
 */
{
  const double a = 5. / 6.;
  const double x2a = pow(x, 2. * a), x22 = x * x / 4.;
  double x2n;               // x^2.a, etc
  double s = 0.0;
  int n;

  const double Ga[11] = { 0, 12.067619015983075, 5.17183672113560444,
      0.795667187867016068, 0.0628158306210802181, 0.00301515986981185091,
      9.72632216068338833e-05, 2.25320204494595251e-06, 3.93000356676612095e-08,
      5.34694362825451923e-10, 5.83302941264329804e-12 };

  const double Gma[11] = { -3.74878707653729304, -2.04479295083852408,
      -0.360845814853857083, -0.0313778969438136685, -0.001622994669507603,
      -5.56455315259749673e-05, -1.35720808599938951e-06,
      -2.47515152461894642e-08, -3.50257291219662472e-10,
      -3.95770950530691961e-12, -3.65327031259100284e-14 };

  x2n = 0.5;                           // init (1/2) * x^0

  s = Gma[0] * x2a;
  s *= x2n;

  // prepare recurrence iteration for next step
  x2n *= x22;    // x^n

  for (n = 1; n <= 10; n++) {

    s += (Gma[n] * x2a + Ga[n]) * x2n;
    // prepare recurrence iteration for next step
    x2n *= x22;    // x^n
  }
  return s;
}

double
asymp_macdo(double x)
/* DOCUMENT asymp_macdo(x)

 Computes a term involved in the computation of the phase struct
 function with a finite outer scale according to the Von-Karman
 model. The term involves the MacDonald function (modified bessel
 function of second kind) K_{5/6}(x), and the algorithm uses the
 asymptotic form for x ~ infinity.
 Warnings :
 - This function makes a doubleing point interrupt for x=0
 and should not be used in this case.
 - Works only for x>0.

 SEE ALSO:
 */
{
  // k2 is the value for
  // gamma_R(5./6)*2^(-1./6)
  const double k2 = 1.00563491799858928388289314170833;
  const double k3 = 1.25331413731550012081;   //  sqrt(pi/2)
  const double a1 = 0.22222222222222222222;   //  2/9
  const double a2 = -0.08641975308641974829;  //  -7/89
  const double a3 = 0.08001828989483310284;   // 175/2187
  double res;
  double x_1;

  x_1 = 1. / x;
  res = k2
      - k3 * exp(-x) * pow(x, 1 / 3.)
          * (1.0 + x_1 * (a1 + x_1 * (a2 + x_1 * a3)));
  return res;
}

double
rodconan(double r, double L0, int k)
/* DOCUMENT rodconan(r,L0,k=)
 The phase structure function is computed from the expression
 Dphi(r) = k1  * L0^(5./3) * (k2 - (2.pi.r/L0)^5/6 K_{5/6}(2.pi.r/L0))

 For small r, the expression is computed from a development of
 K_5/6 near 0. The value of k2 is not used, as this same value
 appears in the series and cancels with k2.
 For large r, the expression is taken from an asymptotic form.

 SEE ALSO:
 */
{

  const double pi = 3.1415926535897932384626433;
  double res = 0;

  // k1 is the value of :
  // 2*gamma_R(11./6)*2^(-5./6)*pi^(-8./3)*(24*gamma_R(6./5)/5.)^(5./6);
  const double k1 = 0.1716613621245709486;
  const double dprf0 = (2 * pi / L0) * r;
  // k2 is the value for gamma_R(5./6)*2^(-1./6),
  // but is now unused
  // k2 = 1.0056349179985892838;

  // Xlim = 0.75*2*pi;   // = 4.71239
  if (dprf0 > 4.71239) {
    res = asymp_macdo(dprf0);
  } else {
    res = -macdo_x56(dprf0, k);
  }
  res *= k1 * pow(L0, 5. / 3);
  return res;
}

double
DPHI(double x, double y, long indexL0, double *rr, double **tabDPHI,
    double convert)
/* DOCUMENT dphi = DPHI(x,y,indexL0,rr,tabDPHI,convert) * r0^(-5./3)
 <x> & <y>         :  separation between apertures
 <indexL0>         :  index for the L0 taken into account
 <rr>              :  array of distance between apertures
 <tabDPHI>         :  array of precomputed DPHI
 <convert>         :  relation between the index on tabDPHI and (x,y)

 Computes the phase structure function for a separation (x,y).
 The r0 is not taken into account : the final result of DPHI(x,y,L0)
 has to be scaled with r0^-5/3, with r0 expressed in meters, to get
 the right value.

 SEE ALSO:
 */
{
  double r = sqrt(x * x + y * y);
  long i0 = (long) (r * convert);
  long i1 = i0 + 1;

  return ((r - rr[i0]) * tabDPHI[indexL0][i1]
      + (rr[i1] - r) * tabDPHI[indexL0][i0]);
}

double**
tabulateDPHI(double* rr,struct tomo_struct tomo, long Ndphi, long *indexL0)
/* DOCUMENT tabDPHI = tabulateDPHI(rr,tomo,Ndphi, indexL0)
 <rr>              :  array of distance between apertures
 <tomo>            :  structure with all the needed information
 <Ndphi>           :  size of rr
 <indexL0>         :  link between the index of the studied layer and the index of the precomputed one. 

 Computes the phase structure function for a separation rr(x,y).
 The r0 is not taken into account : the final result of DPHI(x,y,L0)
 has to be scaled with r0^-5/3, with r0 expressed in meters, to get
 the right value.

 Computes the phase structure for each different L0 and give a array (indexL0) to link the index of the layer i and the index of tabDPHI : for the layer l, DPHI = DPHI( du, dv, indexL0[l],rr,tabDPHI, convert).
 SEE ALSO: DPHI
 */
{
  //Search the different L0 and build indexL0
  const long cNlayer = tomo.Nlayer;
  long i, j, l;
  int cpt = 1;
  double tmp[cNlayer];

  tmp[0] = tomo.L0[0];
  indexL0[0] = 0;

  for (i = 1; i < cNlayer; i++) {
    j = 0;
    const double l0 = tomo.L0[i];

    while ((j < cpt) && (tmp[j] != l0)) {
      j++;
    }
    indexL0[i] = j;

    if (j == cpt) {
      tmp[j] = l0;
      cpt++;
    }
  }

  const int Nl0 = cpt;
  double L0diff[Nl0];
  for (i = 0; i < Nl0; i++) {
    L0diff[i] = tmp[i];
  }

  //précalcul de DPHI : que pour chaque différent L0
  double ** tabDPHI = arr2dAlloc(Nl0, Ndphi);

  for (l = 0; l < Nl0; l++) {
#ifdef USE_OPENMP
#pragma omp parallel num_threads(tomo.ncpu)
#pragma omp for nowait
#endif
    for (j = 0; j < Ndphi; j++) {
      tabDPHI[l][j] = rodconan(rr[j], L0diff[l], 10);
    }
  }
  return tabDPHI;
}

/*
 *
 *  _____           _   _                    ___  _ _                     
 * |  ___|__  _ __ | |_(_) ___  _ __  ___   / _ \| | | ___   __ 
 * | |_ / _ \| '_ \| __| |/ _ \| '_ \/ __| | /_\ | | |/ _ \ / _|
 * |  _| (_) | | | | |_| | (_) | | | \__ \ | __  | | | (_) | (_| 
 * |_|  \___/|_| |_|\__|_|\___/|_| |_|___/ |_| |_|_|_|\___/ \__|
 *
 */

double***
arr3dAlloc(long Nw, long *Nsubap, long Nlayer) {
  /* DOCUMENT  array = arr3dAlloc(Nw,Nsubap,Nlayer)
   <Nw>                  :  number of WFS
   <Nsubap>              :  array of the number of subap of each WFS
   <Nlayer>              :  number of layer

   Allocates a 3d array in one big array with a non-constant size for the 2rd dimension.
   */

  double ***array = (double***) malloc(Nw * sizeof(double**));
  long n, i;

  for (n = 0; n < Nw; n++) {
    array[n] = (double**) malloc( Nsubap[n] * sizeof(double*));
    for (i = 0; i < Nsubap[n]; i++) {
      array[n][i] = (double*) malloc(Nlayer * sizeof(double));
    }
  }

  return array;

}

void
arr3dFree(double ***array, long Nw, long *Nsubap) {

  /* DOCUMENT  array = arr3dFree(array,Nw,Nsubap)
   <Nw>                  :  number of WFS
   <Nsubap>              :  array of the number of subap of each WFS

   Free a 3d array with a non-constant size for the  2nd dimension.
   */
  long n, i;
  for (n = 0; n < Nw; n++) {
    for (i = 0; i < Nsubap[n]; i++) {
      free(array[n][i]);
    }
    free(array[n]);
  }
  free(array);
}

double **
arr2dAlloc(long nbLin, long nbCol)
/* DOCUMENT  array = arr2dAlloc(nblin,nbcol)

 Allocates a 2d array (double).
 */
{
  double **tableau = (double **) malloc(sizeof(double*) * nbLin);
  double *tableau2 = (double *) malloc(sizeof(double) * nbCol * nbLin);
  long i;

  for (i = 0; i < nbLin; i++) {
    tableau[i] = &tableau2[i * nbCol];
  }
  return tableau;
}

void
arr2dFree(double **tableau)
/* DOCUMENT  arr2dFree(array)

 Free a 2d array (double).
 */
{
  free(tableau[0]);
  free(tableau);
}

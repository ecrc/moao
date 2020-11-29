#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>
#include <plasma.h>
#include <quark.h>

#include <cblas.h>
#include <lapacke.h>

//#include <common.h>
#include "flops.h"
#ifdef USE_INTERSAMPLE
#include "intersample.h"
#endif

//TODO 
#include "check_utils.h"
#ifdef CHECK_PLASMA
#include "check_utils.h"
#include "matcov_tiled.h"
//#include "moao_lapack.h"
#endif

#ifdef USE_MATCOV_TILED
#include "matcov_tiled.h"
#endif

//#ifdef USE_TRACE
//#endif

//helper functions && globals ========================================
static int
startswith(const char *s, const char *prefix) {
  size_t n = strlen( prefix );
  if (strncmp( s, prefix, n ))
    return 0;
  return 1;
}

static void
get_thread_count(int *thrdnbr) {
#if defined WIN32 || defined WIN64
  sscanf( getenv( "NUMBER_OF_PROCESSORS" ), "%d", thrdnbr );
#else
  *thrdnbr = sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

double cWtime(void)
{
  struct timeval tp;
  gettimeofday( &tp, NULL );
  return tp.tv_sec + 1e-6 * tp.tv_usec;
}

#ifdef USE_MATCOV_TILED
int PLASMA_GenCovMat_Tile_Async(PLASMA_desc     *descA,
                                PLASMA_sequence *sequence,
                                PLASMA_request  *request,
                                struct tomo_struct *tomo,
                                int part);
#endif

#ifdef USE_INTERSAMPLE
int PLASMA_Intersample_genPSF(double     *A,
                              struct isample_struct *isample,
                              int nact, int iloop, int oloop, int gal,
                              PLASMA_sequence *sequence,
                              PLASMA_request  *request
);
void ao_pdtile_to_lapack_quark(PLASMA_desc* A, double *Af77, int lda,
                               PLASMA_sequence *sequence, PLASMA_request *request);
#endif

int PLASMA_dscaldiag_Tile_Async(double alpha, PLASMA_desc *descA, PLASMA_sequence *sequence, PLASMA_request  *request);



#if defined(USE_TRACE)
#else 
#define START_TRACING() if( 0 ) {};
#define STOP_TRACING() if( 0 ) {};
#endif

#define START_TIMING(_t)               \
START_TRACING();                    \
_t = -cWtime();

#define STOP_TIMING(_t)                \
_t += cWtime();                      \
STOP_TRACING();

#ifndef USE_MATCOV_TILED
#define USAGE fprintf(stderr,"usage: --n_cores=cores --tile=tile_size --nmeas=nmeas --nmeasts=nbmeasts --nact=nbact --maxrefine=maxrefine --maxobs=maxobs --quark --ompss [--v] [--s] [--h]\n"\
"--n_cores: number of cores\n" \
"--tile: tunable tile size\n" \
"--nmeas: number of measurements\n" \
"--nmeasts: number of measurements in the true sensor\n" \
"--nact: number of actuators\n" \
"--maxrefine: max number of refinements\n" \
"--maxobs: max number of observations\n" \
"--quark: select QUARK runtime\n" \
"--ompss: select OmpSs runtime\n" \
"--v: verbose\n"\
"--s: synchronize between stages, default off\n"\
"--h: print this help\n"\
);
#else
#define USAGE fprintf(stderr,"usage: --n_cores=cores --tile=tile_size --nssp=nssp --nact=nbact --maxrefine=maxrefine --maxobs=maxobs --quark --ompss [--v] [--s] [--h]\n"\
"--n_cores: number of cores\n" \
"--tile: tunable tile size\n" \
"--nssp: number of measurements\n" \
"--nact: number of actuators\n" \
"--maxrefine: max number of refinements\n" \
"--maxobs: max number of observations\n" \
"--quark: select QUARK runtime\n" \
"--ompss: select OmpSs runtime\n" \
"--v: verbose\n"\
"--s: synchronize between stages, default off\n"\
"--h: print this help\n"\
);
#endif

#define HANDLE_ERROR(_result, _func)               \
  if(PLASMA_SUCCESS != _result) {                  \
    fprintf(stderr,"An error occured (%d), when calling %s at line: %d... \n\n", _result, _func, __LINE__ -1);   \
    break;                                         \
  }
#define HANDLE_ERROR_RET(_result, _func)               \
  if(PLASMA_SUCCESS != _result) {                  \
    fprintf(stderr,"An error occured (%d), when calling %s at line: %d... \n\n", _result, _func, __LINE__ -1);   \
    exit(-1);                                         \
  }

//---------------------------------------------------------------------------
 
//#ifdef USE_YDATA
//#define DATA_FILES_PATH "../datafile/check/"
//#else
//#define DATA_FILES_PATH "../datafile/"
//#endif
  
//---------------------------------------------------------------------------
int isSPD(double* A, int m, int n, int lda){
  int c, sum, r;
  if(m != n) return 1;
  for(c = 0; c < m; c++){
    if(A[c + c* lda] <= 0)
      return 2;
    sum = 0;
    for(r = 0; r < m; r++){
      if(A[r + c*lda] != A[c + r*lda])
        return 3;
      if(r != c)
        sum += fabs(A[r + c*lda]);
    }
    if(sum >= A[c + c*lda])
      return 4;
  }
  return 0;
}
//---------------------------------------------------------------------------
int makeSPD(double* A, int m, int n, int lda, int bigM){
  int c;
  if(m != n) return 1;
  for(c = 0; c < m; c++){
    A[c + c* lda] += (double)bigM;
  }
  return 0;
}
//---------------------------------------------------------------------------
//#define BLKADDR(A, type, m, n)  (type *)plasma_getaddr(A, m, n)
//#define A(m,n) BLKADDR(A, double, m, n)
#define BLKLDD(A, k) ( ( (k) + (A).i/(A).mb) < (A).lm1 ? (A).mb : (A).lm%(A).mb )
#define BLKLDD2(A, k) ( ( (k) + (A).j/(A).nb) < (A).ln1 ? (A).nb : (A).ln%(A).nb )
#define A(m,n) (double *)plasma_getaddr(A, m, n)
#define my_imax(m,n) ( (m) > (n) ? (m) : (n) )

//main ===========================================================
int main( int argc, char** argv)
{
char DATA_FILES_PATH[100]="../datafile/";
char suffix[40]="";
char dataFile[150]="";

  int ncores = -1, result=0;
  int maxobs = 1, maxrefine = 1;
#ifdef USE_MATCOV_TILED
  fprintf(stderr,"use_matcov_tiled\n");
  struct tomo_struct tomo;
#endif
#ifdef USE_INTERSAMPLE
  struct isample_struct isample;
//  struct isample_struct isampleL;
#endif
  if(argc < 2)
  {
    USAGE;
    return 0;
  }
  double alpha, beta;
  int nact=512;
#ifdef USE_MATCOV_TILED
  int nssp = 10;
#endif
int i,j;
  int nmeas=0, nmeasts=0, ts=0, verbose=0, sync=0, printhelp=0, unknown=0, runtime=0;
  double flops=0.0;

  //parse options
  int argi;
  for(argi = 1; argi < argc; argi++)
  {
    if(startswith( argv[argi], "--n_cores=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &ncores );
    }
    else
    if(startswith( argv[argi], "--nact=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &nact );
      if(nact <= 0){
        fprintf(stderr,"Invalid number of actuators\n"); return 0;
      }
    }
#ifndef USE_MATCOV_TILED
    else
    if(startswith( argv[argi], "--nmeas=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &nmeas );
      if(nmeas <= 0){
        fprintf(stderr,"Invalid number of total measurements\n"); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--nmeasts=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &nmeasts );
      if(nmeasts <= 0){
        fprintf(stderr,"Invalid number of measurements of the true sensor\n"); return 0;
      }
    }
#else
    else
    if(startswith( argv[argi], "--nssp=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &nssp );
      if(nssp <= 0){
        fprintf(stderr,"Invalid number of total measurements\n"); return 0;
      }
    }
#endif
    else
    if(startswith( argv[argi], "--maxobs=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &maxobs );
      if(maxobs < 0){
        fprintf(stderr,"Invalid number of max obs\n"); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--maxrefine=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &maxrefine );
      if(maxrefine < 0){
        fprintf(stderr,"Invalid number of max refinement\n"); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--tile=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &ts );
      if(ts <= 0){
        fprintf(stderr,"Invalid tile size\n"); return 0;
      }
    }

    else if(startswith( argv[argi], "--filePath=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%s", &DATA_FILES_PATH );
    }
    else if(startswith( argv[argi], "--suffix=")){
        //--suffix start with --s :this condition must be tested before --s
        sscanf( strchr( argv[argi], '=' ) + 1, "%s", &suffix);
    }
    else
    if(startswith( argv[argi], "--quark"))
      runtime = 0;
    else
    if(startswith( argv[argi], "--ompss"))
      runtime = 1;
    else
    if(startswith( argv[argi], "--v"))
      verbose = 1;
    else
    if(startswith( argv[argi], "--s"))
      sync = 1;
    else
    if(startswith( argv[argi], "--h"))
      printhelp = 1;
    else
    {
      fprintf(stderr,"Unknown option %s, aborting...\n", argv[argi]); unknown=1;
    }
  }

  
#if USE_INTERSAMPLE
    long naxes[2];
  concatenateFileName(dataFile,DATA_FILES_PATH,"Dx",suffix);
  getFitsDims(dataFile,naxes);
  if(naxes[0]>0){
      nact=naxes[0];
  }
#endif

#ifdef USE_MATCOV_TILED
  result = !matcov_init_tomo_tiled(&tomo, nssp, DATA_FILES_PATH, 0, 0, 0, -1, 0., 0.);
  HANDLE_ERROR_RET(result, "matcov_init_tomo_tiled");

  nmeas = matcov_getNumMeasurements(&tomo);
  nmeasts = matcov_getNumMeasurementsTS(&tomo);
  nmeas-=nmeasts;
  //nact = nmeasts;
  matcov_set_gal_coords(&tomo, 0, 0);
#endif

#ifdef USE_INTERSAMPLE
  intersample_prepare(&isample, nact*nact, tomo.Nssp[tomo.Nw-1]+1, tomo.DiamTel, DATA_FILES_PATH);
  intersample_init(&isample);

  double *cumulatedPSF=(double *)calloc(isample.N*isample.N,sizeof(double));

#endif
  if(printhelp||unknown)
  {
    USAGE;
    return 0;
  }
  if(nmeas <= 0)
  {
    fprintf(stderr,"Please provide number of measurements\n");
    USAGE;
    return 0;
  }
  if(nmeasts <= 0)
  {
    fprintf(stderr,"Please provide number of measurements in the true sensor\n");
    USAGE;
    return 0;
  }
  if(nact <= 0)
  {
    fprintf(stderr,"Please provide number of actuators\n");
    USAGE;
    return 0;
  }
  if(maxrefine < 0)
  {
    fprintf(stderr,"Please provide number of max refine\n");
    USAGE;
    return 0;
  }
  if(maxobs < 0)
  {
    fprintf(stderr,"Please provide number of max obs\n");
    USAGE;
    return 0;
  }
  if(ts <= 0)
  {
    fprintf(stderr,"Please provide tile size\n");
    USAGE;
    return 0;
  }

  //defaults
  if(ncores < 0)
  {
    get_thread_count(&ncores);
  }

  if(verbose) fprintf(stderr, "\nProcessing %d number of measurements and %d number of actuators\n", nmeas, nact);
  if(verbose) fprintf(stderr, "Tile size (tunable): %d\n", ts);
  if(verbose) fprintf(stderr, "Working on: %d CPU cores\n\n", ncores);

  //PLASMA_Init(ncores);

#ifdef PLASMA_HAS_OMPSS
  if(runtime) {
    PLASMA_Runtime_Init(ncores, PLASMA_OMPSS);
    PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_OMPSS);
    if(verbose) fprintf(stderr, "PLASMA initializing with OmpSs... done\n\n");
  }
  else {
    PLASMA_Runtime_Init(ncores, PLASMA_QUARK);
    PLASMA_Set(PLASMA_RUNTIME_MODE, PLASMA_QUARK);
    if(verbose) fprintf(stderr, "PLASMA initializing with QUARK... done\n\n");
  }
#else
  PLASMA_Init(ncores);
  if(verbose) fprintf(stderr, "PLASMA initializing with QUARK... done\n\n");
#endif 
  /*
  plasma_context_t *plasma;
  plasma = plasma_context_self();
  plasma_setlapack_sequential(plasma);
  */
 
  PLASMA_Disable(PLASMA_AUTOTUNING);
  PLASMA_Set(PLASMA_TILE_SIZE, ts );
  PLASMA_Set(PLASMA_TRANSLATION_MODE, PLASMA_INPLACE);
  PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING);
  PLASMA_Enable(PLASMA_WARNINGS);
  PLASMA_Enable(PLASMA_ERRORS);


  double total_time = 0.0; 
  int nbobs, nbrefine;
  // time break down
  double alloc_cpu_time = 0.0;
  double preprocess_cpu_time = 0.0;
  double compute_time = 0.0;
  double total_refine_time = 0.0, refine_time = 0.0;
  double total_obs_time = 0.0, obs_time = 0.0;
  double total_intersample_time = 0.0, intersample_time = 0.0;
  double total_matcov_time = 0.0, matcov_time = 0.0;
  
  double *Cmm = NULL, *Cpm = NULL, *Cpp = NULL;
  PLASMA_desc *descCmm = NULL, *descCpm = NULL, *descCpp = NULL;

  double *R = NULL, *Cee = NULL, *Cvv = NULL, *Cvv_lap = NULL, *Dx = NULL, *Tmp = NULL;
  PLASMA_desc *descR = NULL, *descCee = NULL, *descCvv = NULL, *descDx = NULL, *descTmp = NULL;

#ifdef CHECK_PLASMA
    double *lCmm=NULL, *lCpm=NULL, *lCpp=NULL, *lR=NULL, *lCee=NULL, *lCvv=NULL, *lDx=NULL;
    double *lTmp;
#endif

  PLASMA_sequence *sequence;
  PLASMA_request request[19] = { PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,
                                 PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,
                                 PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,
                                 PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,
                                 PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER};
  PLASMA_Sequence_Create(&sequence);
  
#ifdef USE_INTERSAMPLE
  //PLASMA_sequence *sequence_intersample;
  //PLASMA_Sequence_Create(&sequence_intersample);
#endif
  // Allocate matrix memory
  {
    if(verbose) fprintf(stderr, "Allocating matrix memory...");
    alloc_cpu_time -= cWtime();
    Cmm = (double*)calloc( (size_t)nmeas * nmeas , sizeof(double) );
    if ( ! Cmm ) {
      fprintf(stderr, "Out of Memory for Covariance Matrix (Cmm)\n");
      return -1;
    }
#ifdef PLASMA_HAS_OMPSS
    if(runtime) {
       #pragma omp register ([nmeas*nmeas]Cmm)
    }
#endif
    Cpm = (double*)calloc( (size_t)nmeasts * nmeas , sizeof(double) );
    if ( ! Cpm ) {
      fprintf(stderr, "Out of Memory for Cpm Matrix\n");
      return -1;
    }
#ifdef PLASMA_HAS_OMPSS
    if(runtime) {
       #pragma omp register ([nmeasts*nmeas]Cpm)
    }
#endif
    Cpp = (double*)calloc( (size_t)nmeasts * nmeasts , sizeof(double) );
    if ( ! Cpp ) {
      fprintf(stderr, "Out of Memory for Cpp Matrix\n");
      return -1;
    }
#ifdef PLASMA_HAS_OMPSS
    if(runtime) {
       #pragma omp register ([nmeasts*nmeasts]Cpp)
    }
#endif
    R = (double*)calloc( (size_t)nmeasts * nmeas , sizeof(double) );
    if ( ! R ) {
      fprintf(stderr, "Out of Memory for ToR (R)\n");
      return -1;
    }
#ifdef PLASMA_HAS_OMPSS
    if(runtime) {
       #pragma omp register ([nmeasts*nmeas]R)
    }
#endif
    Cee = (double*)calloc( (size_t)nmeasts * nmeasts , sizeof(double) );
    if ( ! Cee ) {
      fprintf(stderr, "Out of Memory for Cee Matrix\n");
      return -1;
    }
#ifdef PLASMA_HAS_OMPSS
    if(runtime) {
       #pragma omp register ([nmeasts*nmeasts]Cee)
    }
#endif
    Cvv = (double*)calloc( (size_t)nact * nact , sizeof(double) );
    if ( ! Cvv ) {
      fprintf(stderr, "Out of Memory for Cvv Matrix\n");
      return -1;
    }
#ifdef PLASMA_HAS_OMPSS
    if(runtime) {
       #pragma omp register ([nact*nact]Cvv)
    }
#endif
#ifdef USE_INTERSAMPLE
    Cvv_lap = (double*)calloc( (size_t)nact * nact , sizeof(double) );
    if ( ! Cvv_lap ) {
      fprintf(stderr, "Out of Memory for Cvv_lap Matrix\n");
      return -1;
    }
#endif
    Dx = (double*)calloc( (size_t)nact * nmeasts , sizeof(double) );
    if ( ! Dx ) {
      fprintf(stderr, "Out of Memory for Dx Matrix\n");
      return -1;
    }
#ifdef PLASMA_HAS_OMPSS
    if(runtime) {
       #pragma omp register ([nact*nmeasts]Dx)
    }
#endif
    Tmp = (double*)calloc( (size_t)nact * nmeasts , sizeof(double) );
    if ( ! Tmp ) {
      fprintf(stderr, "Out of Memory for Tmp Matrix\n");
      return -1;
    }
#ifdef PLASMA_HAS_OMPSS
    if(runtime) {
       #pragma omp register ([nact*nmeasts]Tmp)
    }
#endif
    alloc_cpu_time += cWtime();
    total_time += alloc_cpu_time;
    if(verbose) fprintf(stderr, "Done in %f(s)\n\n", alloc_cpu_time);
  }


#ifdef CHECK_PLASMA
    lCmm = (double*)calloc( (size_t)nmeas * nmeas , sizeof(double) );
    if(!lCmm){fprintf(stderr,"Out of Memory for lCmm");return -1;}
    lCpm = (double*)calloc( (size_t)nmeasts * nmeas , sizeof(double) );
    if(!lCpm){fprintf(stderr,"Out of Memory for lCpm");return -1;}
    lCpp = (double*)calloc( (size_t)nmeasts * nmeasts , sizeof(double) );
    if(!lCpp){fprintf(stderr,"Out of Memory for lCpp");return -1;}
    lR   = (double*)calloc( (size_t)nmeasts * nmeas , sizeof(double) );
    if(!lR){fprintf(stderr,"Out of Memory for lR");return -1;}
    lCee = (double*)calloc( (size_t)nmeasts * nmeasts , sizeof(double) );
    if(!lCee){fprintf(stderr,"Out of Memory for lCee");return -1;}
    lCvv = (double*)calloc( (size_t)nact * nact , sizeof(double) );
    if(!lCvv){fprintf(stderr,"Out of Memory for lCvv");return -1;}
    lDx  = (double*)calloc( (size_t)nact * nmeasts , sizeof(double) );
    if(!lDx){fprintf(stderr,"Out of Memory for lDx");return -1;}
    lTmp = (double*)calloc( (size_t)nmeas * nmeas , sizeof(double) );
    if(!lTmp){fprintf(stderr,"Out of Memory for lTmp");return -1;}

double t_Cmm,t_Ctt,t_Ctm,t_R,t_CeeCvv;
double maxerr,errv1,errv2;
char logName[50];
strcpy(logName,"log");
strcat(logName,suffix);
strcat(logName,".txt");
FILE *logFile=fopen(logName,"w");
if(logFile==NULL){
  fprintf(stderr,"could not open log file\n");
}
fprintf(logFile,"nmeas: %d\nnmeasts: %d\nnatc: %d\nnLayers: %d\ntile size: %d\nncores: %d\npsf:0\n",nmeas,nmeasts,nact,tomo.Nlayer,ts,ncores);
fprintf(logFile,"Matrices times(P)   max_err(P/Y) err_val(P) err_val(Y)\n");
#else
    double *lDx,*lTmp;
    lDx  = (double*)calloc( (size_t)nact * nmeasts , sizeof(double) );
    if(!lDx){fprintf(stderr,"Out of Memory for lDx");return -1;}
    lTmp  = (double*)calloc( (size_t)nact * nmeasts , sizeof(double) );
    if(!lDx){fprintf(stderr,"Out of Memory for lDx");return -1;}

#endif

  PLASMA_Desc_Create(&descCmm, Cmm, PlasmaRealDouble, ts, ts, ts*ts, nmeas, nmeas, 0, 0, nmeas, nmeas);
  PLASMA_Desc_Create(&descCpm, Cpm, PlasmaRealDouble, ts, ts, ts*ts, nmeasts, nmeas, 0, 0, nmeasts, nmeas);
  PLASMA_Desc_Create(&descCpp, Cpp, PlasmaRealDouble, ts, ts, ts*ts, nmeasts, nmeasts, 0, 0, nmeasts, nmeasts);
  PLASMA_Desc_Create(&descR,   R,   PlasmaRealDouble, ts, ts, ts*ts, nmeasts, nmeas, 0, 0, nmeasts, nmeas);
  PLASMA_Desc_Create(&descCee, Cee, PlasmaRealDouble, ts, ts, ts*ts, nmeasts, nmeasts, 0, 0, nmeasts, nmeasts);
  PLASMA_Desc_Create(&descCvv, Cvv, PlasmaRealDouble, ts, ts, ts*ts, nact, nact, 0, 0, nact, nact);
  PLASMA_Desc_Create(&descDx,  Dx,  PlasmaRealDouble, ts, ts, ts*ts, nact, nmeasts, 0, 0, nact, nmeasts);
  PLASMA_Desc_Create(&descTmp, Tmp, PlasmaRealDouble, ts, ts, ts*ts, nact, nmeasts, 0, 0, nact, nmeasts);
  
  if(verbose) fprintf(stderr, "PLASMA: MOAO started\n\n");
  if(verbose) fprintf(stderr, "PLASMA: MOAO Loading the interaction matrix Dx\n");

  preprocess_cpu_time -= cWtime();
//read Dx from file
#ifdef USE_INTERSAMPLE
    concatenateFileName(dataFile,DATA_FILES_PATH,"Dx",suffix);
    readFits(dataFile,lTmp);
    copy(lTmp,lDx,nact,nmeasts,'A','T');
    p_lapack_to_tile(*descDx,lDx,0);
#else
  PLASMA_dplrnt_Tile_Async( descDx, rand()%100, sequence, &request[0] );
#endif
  if(sync){
#ifdef PLASMA_HAS_OMPSS
    if(runtime)
      #pragma omp taskwait
    else
#endif
      PLASMA_Sequence_Wait(sequence);
  }
  preprocess_cpu_time += cWtime();
  if(verbose) fprintf(stderr, "Done in %f(s)\n\n", preprocess_cpu_time);

  total_time += preprocess_cpu_time;

    
  if(verbose) fprintf(stderr, "PLASMA: MOAO main computation started with\n");
  if(verbose) fprintf(stderr, "Maximum number of refinements %d\n", maxrefine);
  if(verbose) fprintf(stderr, "Maximum number of observations %d\n", maxobs);

  START_TIMING(compute_time);

  for (nbrefine = 0; nbrefine < my_imax(maxrefine, 1); nbrefine++) {

    if(maxrefine > 0){
       START_TIMING(refine_time);

       if(verbose) fprintf(stderr, "\nnbrefine %d\n", nbrefine+1);
       // plasma_setlapack_multithreads(plasma->world_size);
       // MatCov generation
       // plasma_setlapack_sequential(plasma);
       
       if(sync){
#ifdef PLASMA_HAS_OMPSS
         if(runtime)
           #pragma omp taskwait
         else
#endif
           PLASMA_Sequence_Wait(sequence);
           START_TIMING(matcov_time);
       }
#ifdef USE_MATCOV_TILED
      //need to get new atm-params
      int night_idx = 0,
          snapshots_per_night = 8,
          snapshot_idx = nbrefine,
          obs_idx = -1;
#ifdef CHECK_PLASMA
START_TIMING(t_Cmm);
#endif
      result = !matcov_update_atm_params(&tomo, night_idx, snapshots_per_night, snapshot_idx, obs_idx);
        HANDLE_ERROR(result, "matcov_update_atm_params");
       PLASMA_GenCovMat_Tile_Async( descCmm, sequence, &request[1], &tomo, 1 );
       HANDLE_ERROR(result, "PLASMA_GenCovMat_Tile_Async");
       flops = flops + FLOPS_DCovMat(nmeas, nmeas, tomo.Nlayer, 1);
#else
       result = PLASMA_dplgsy_Tile_Async( (double)nmeas, descCmm, rand()%100, sequence, &request[1] );
       HANDLE_ERROR(result, "PLASMA_dplgsy_Tile_Async");
       flops = flops + FLOPS_DCovMat(nmeas, nmeas);
#endif
       if(sync){
#ifdef PLASMA_HAS_OMPSS
         if(runtime)
           #pragma omp taskwait
         else
#endif
           PLASMA_Sequence_Wait(sequence);
#ifdef CHECK_PLASMA
STOP_TIMING(t_Cmm);
#endif
	
       }

#ifdef USE_MATCOV_TILED
#ifdef CHECK_PLASMA
START_TIMING(t_Ctm);
#endif
       result = PLASMA_GenCovMat_Tile_Async( descR, sequence, &request[2], &tomo, 3 );
       HANDLE_ERROR(result, "PLASMA_GenCovMat_Tile_Async");
       flops = flops + FLOPS_DCovMat(nmeas, nmeasts, tomo.Nlayer, 3);
#else
       result = PLASMA_dplrnt_Tile_Async( descR, rand()%100, sequence, &request[2] );
       HANDLE_ERROR(result, "PLASMA_dplrnt_Tile_Async");
       flops = flops + FLOPS_DCovMat(nmeas, nmeasts);
#endif
    
       if(sync){
#ifdef PLASMA_HAS_OMPSS
         if(runtime)
           #pragma omp taskwait
         else
#endif
           PLASMA_Sequence_Wait(sequence);
         STOP_TIMING(matcov_time);
         total_matcov_time += matcov_time;
#ifdef CHECK_PLASMA
STOP_TIMING(t_Ctm);
#endif
       }


#ifdef CHECK_PLASMA

//    if(naxes[0]>0){
//    concatenateFileName(dataFile,DATA_FILES_PATH,"Dx",suffix);
//    readFits(dataFile,lTmp);
//    copy(lTmp,lDx,nact,nmeasts,'A','T');
//    p_lapack_to_tile(*descDx,lDx,0);
//    }
//    else{
//        p_tile_to_lapack(*descDx,lDx,0);
//    }

    fprintf(stderr,"CHECK: COMPARE MATRIX GENERATION\n");  
    p_tile_to_lapack(*descR, lR,0);
    concatenateFileName(dataFile,DATA_FILES_PATH,"ctm",suffix);
    readFits(dataFile,lTmp);
    compareMatrices2(lR,lTmp,nmeasts,nmeas,'A','T',&maxerr,&errv1,&errv2);
    fprintf(logFile,"Ctm      %e   %e     %e   %e\n",t_Ctm,maxerr,errv1,errv2);
    
    p_tile_to_lapack(*descCmm, lCmm,0);
    concatenateFileName(dataFile,DATA_FILES_PATH,"cmm",suffix);
    readFits(dataFile,lTmp);
    compareMatrices2(lCmm,lTmp,nmeas,nmeas,'L','T',&maxerr,&errv1,&errv2);
    fprintf(logFile,"Cmm      %e   %e     %e   %e\n",t_Cmm,maxerr,errv1,errv2);

    fprintf(stderr,"END CHECK: COMPARE MATRIX GENERATION\n\n");  
#endif

#ifdef CHECK_PLASMA
START_TIMING(t_R);
#endif
       //
       result = PLASMA_dpotrf_Tile_Async(PlasmaLower, descCmm, sequence, &request[3]);
       HANDLE_ERROR(result, "PLASMA_dpotrf_Tile_Async");
       flops = flops + FLOPS_DPOTRF(nmeas);

       if(sync){
#ifdef PLASMA_HAS_OMPSS
         if(runtime)
           #pragma omp taskwait
         else
#endif
           PLASMA_Sequence_Wait(sequence);
       }
       
       alpha = 1.0;
       result = PLASMA_dtrsm_Tile_Async(PlasmaRight, PlasmaLower, PlasmaTrans, PlasmaNonUnit,
                                        alpha, descCmm, descR, sequence, &request[4]);
       HANDLE_ERROR(result, "PLASMA_dtrsm_Tile_Async");
       flops = flops + FLOPS_DTRSM(PlasmaRight, nmeasts, nmeas);

       if(sync){
#ifdef PLASMA_HAS_OMPSS
         if(runtime)
           #pragma omp taskwait
         else
#endif
           PLASMA_Sequence_Wait(sequence);
       }

       alpha = 1.0;
       result = PLASMA_dtrsm_Tile_Async(PlasmaRight, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, 
                                        alpha, descCmm, descR, sequence, &request[5]);
       HANDLE_ERROR(result, "PLASMA_dtrsm_Tile_Async");
       flops = flops + FLOPS_DTRSM(PlasmaRight, nmeasts, nmeas);

       if(sync){
#ifdef PLASMA_HAS_OMPSS
         if(runtime)
           #pragma omp taskwait
         else
#endif
           PLASMA_Sequence_Wait(sequence);
#ifdef CHECK_PLASMA
STOP_TIMING(t_R);
#endif
       }
  

#ifdef CHECK_PLASMA
    fprintf(stderr,"CHECK RECONSTRUCTOR\n");

    p_tile_to_lapack(*descR, lR,0);
    concatenateFileName(dataFile,DATA_FILES_PATH,"R",suffix);
    readFits(dataFile,lTmp);
    compareMatrices2(lR,lTmp,nmeasts,nmeas,'A','T',&maxerr,&errv1,&errv2);
    fprintf(logFile,"R        %e   %e     %e   %e\n",t_R,maxerr,errv1,errv2);

    fprintf(stderr,"END CHECK RECONSTRUCTOR\n\n");
#endif

  /*
       result = PLASMA_dpotri_Tile_Async(PlasmaLower, descCmm, sequence, &request[3]);
       if(PLASMA_SUCCESS != result) break;
       flops = flops + FLOPS_DPOTRI(nmeas);

       if(sync){
         if(runtime)
           #pragma omp taskwait
         else
           PLASMA_Sequence_Wait(sequence);
       }
           
       alpha = 1.0;
       beta = 0.0;
       result = PLASMA_dsymm_Tile_Async(PlasmaRight, PlasmaLower,
                                        alpha, descCmm, descCpm,
                                        beta,  descR,
                                        sequence, &request[4]);
       if(PLASMA_SUCCESS != result) break;
       flops = flops + FLOPS_DSYMM(PlasmaRight, nmeasts, nmeas);

       if(sync){
         if(runtime)
           #pragma omp taskwait
         else
           PLASMA_Sequence_Wait(sequence);
       }
  */
       STOP_TIMING(refine_time);
       total_refine_time += refine_time;
    }
       for (nbobs = 0; nbobs < maxobs; nbobs++) {
            START_TIMING(obs_time);
            if(verbose) fprintf(stderr, " nbobs %d\n", nbobs+1);
            // MatCov generation
  /*
   */
            if(sync){
#ifdef PLASMA_HAS_OMPSS
              if(runtime)
                #pragma omp taskwait
              else
#endif
                PLASMA_Sequence_Wait(sequence);
                START_TIMING(matcov_time);
            }
#ifdef USE_MATCOV_TILED
            //need to get new atm-params
          int night_idx = 0,
              snapshots_per_night = 8,
              snapshot_idx = nbrefine,
              obs_idx = nbobs;
            result = !matcov_update_atm_params(&tomo, night_idx, snapshots_per_night, snapshot_idx, obs_idx);
            HANDLE_ERROR(result, "matcov_update_atm_params");
#ifdef CHECK_PLASMA
START_TIMING(t_Cmm);
#endif
            result = PLASMA_GenCovMat_Tile_Async( descCmm, sequence, &request[6], &tomo, 1 );
            HANDLE_ERROR(result, "PLASMA_GenCovMat_Tile_Async");
            flops = flops + FLOPS_DCovMat(nmeas, nmeas, tomo.Nlayer, 1);
#else
            result = PLASMA_dplgsy_Tile_Async( (double)nmeas, descCmm, rand()%100, sequence, &request[6] );
            HANDLE_ERROR(result, "PLASMA_dplgsy_Tile_Async");
            flops = flops + FLOPS_DCovMat(nmeas, nmeas);
#endif
            if(sync){
#ifdef PLASMA_HAS_OMPSS
              if(runtime)
                #pragma omp taskwait
              else
#endif
                PLASMA_Sequence_Wait(sequence);
#ifdef CHECK_PLASMA
STOP_TIMING(t_Cmm);
#endif
            }
            
#ifdef USE_MATCOV_TILED
#ifdef CHECK_PLASMA
START_TIMING(t_Ctm);
#endif
            result = PLASMA_GenCovMat_Tile_Async( descCpm, sequence, &request[7], &tomo, 3 );
            HANDLE_ERROR(result, "PLASMA_GenCovMat_Tile_Async");
            flops = flops + FLOPS_DCovMat(nmeasts, nmeas, tomo.Nlayer, 3);
#else
            result = PLASMA_dplrnt_Tile_Async( descCpm, rand()%100, sequence, &request[7] );
            HANDLE_ERROR(result, "PLASMA_dplrnt_Tile_Async");
            flops = flops + FLOPS_DCovMat(nmeasts, nmeas);
#endif
            if(sync){
#ifdef PLASMA_HAS_OMPSS
              if(runtime)
                #pragma omp taskwait
              else
#endif
                PLASMA_Sequence_Wait(sequence);
#ifdef CHECK_PLASMA
STOP_TIMING(t_Ctm);
#endif
            }
            
#ifdef USE_MATCOV_TILED
#ifdef CHECK_PLASMA
START_TIMING(t_Ctt);
#endif
            result = PLASMA_GenCovMat_Tile_Async( descCpp, sequence, &request[8], &tomo, 4 );
            HANDLE_ERROR(result, "PLASMA_GenCovMat_Tile_Async");
            flops = flops + FLOPS_DCovMat(nmeasts, nmeasts, tomo.Nlayer, 4);
#else
            result = PLASMA_dplgsy_Tile_Async( (double)nmeasts, descCpp, rand()%100, sequence, &request[8] );
            HANDLE_ERROR(result, "PLASMA_dplgsy_Tile_Async");
            flops = flops + FLOPS_DCovMat(nmeasts, nmeasts);
#endif

            if(sync){
#ifdef PLASMA_HAS_OMPSS
              if(runtime)
                #pragma omp taskwait
              else
#endif
                PLASMA_Sequence_Wait(sequence);
              STOP_TIMING(matcov_time);
              total_matcov_time += matcov_time;
#ifdef CHECK_PLASMA
STOP_TIMING(t_Ctt);
#endif
            }
//      break;
#ifdef CHECK_PLASMA
fprintf(logFile,"Cmm      %e   x               x              x \n",t_Cmm);
fprintf(logFile,"Ctm      %e   x               x              x \n",t_Ctm);
    concatenateFileName(dataFile,DATA_FILES_PATH,"ctt",suffix);
    readFits(dataFile,lTmp);
    p_tile_to_lapack(*descCpp, lCpp,0);
    compareMatrices2(lCpp,lTmp,nmeasts,nmeasts,'A','T',&maxerr,&errv1,&errv2);
fprintf(logFile,"Ctt      %e   %e     %e   %e\n",t_Ctt,maxerr,errv1,errv2);
#endif

#ifdef CHECK_PLASMA
START_TIMING(t_CeeCvv);
#endif
            alpha = -1.0;
            beta = 1.0;
            result = PLASMA_dsyr2k_Tile_Async(PlasmaLower, PlasmaNoTrans, 
                                              alpha, descCpm, descR, 
                                              beta, descCpp, 
                                              sequence, &request[9]);
            HANDLE_ERROR(result, "PLASMA_dsyr2k_Tile_Async");
            flops = flops + FLOPS_DSYR2K(nmeasts, nmeas);
     
            if(sync){
#ifdef PLASMA_HAS_OMPSS
              if(runtime)
                #pragma omp taskwait
              else
#endif
                PLASMA_Sequence_Wait(sequence);
            }

            alpha = 1.0;
            beta = 0.0;
            result = PLASMA_dsymm_Tile_Async(PlasmaRight, PlasmaLower,
                                             alpha, descCmm, descR,
                                             beta,  descCpm,
                                             sequence, &request[10]);
            // Cpm gets overwritten at this point and it is used as a buffer
            HANDLE_ERROR(result, "PLASMA_dsymm_Tile_Async");
            flops = flops + FLOPS_DSYMM(PlasmaRight, nmeasts, nmeas);
 
            if(sync){
#ifdef PLASMA_HAS_OMPSS
              if(runtime)
                #pragma omp taskwait
              else
#endif
                PLASMA_Sequence_Wait(sequence);
            }

            // Rebuild the symmetry for Cpp
            result = PLASMA_dlacpy_Tile_Async( PlasmaLower, descCpp, descCee, 
                                      sequence, &request[11] );
            HANDLE_ERROR(result, "PLASMA_dlacpy_Tile_Async");

            alpha = 1.0;
            result = PLASMA_dtradd_Tile_Async( PlasmaUpper, PlasmaTrans, alpha, descCee, alpha, descCee, 
                                      sequence, &request[12] );
            HANDLE_ERROR(result, "PLASMA_dgeadd_Tile_Async");

            alpha = 0.5;
            result = PLASMA_dscaldiag_Tile_Async(alpha, descCee, sequence, &request[13]);
            HANDLE_ERROR(result, "PLASMA_dscaldiag_Tile_Async");


            alpha = 1.0;
            beta = 1.0;
            result = PLASMA_dgemm_Tile_Async(PlasmaNoTrans, PlasmaTrans, 
                                             alpha, descCpm, descR, 
                                             beta, descCee, 
                                             sequence, &request[14]);
            HANDLE_ERROR(result, "PLASMA_dgemm_Tile_Async");
            flops = flops + FLOPS_DGEMM(nmeasts, nmeasts, nmeas);
     
            if(sync){
#ifdef PLASMA_HAS_OMPSS
              if(runtime)
                #pragma omp taskwait
              else
#endif
                PLASMA_Sequence_Wait(sequence);
            }
     
            alpha = 1.0;
            beta = 0.0;
            result = PLASMA_dsymm_Tile_Async(PlasmaRight, PlasmaLower,
                                             alpha, descCee, descDx,
                                             beta,  descTmp,
                                             sequence, &request[15]);
            HANDLE_ERROR(result, "PLASMA_dsymm_Tile_Async");
            flops = flops + FLOPS_DSYMM(PlasmaRight, nact, nmeasts);
     
            if(sync){
#ifdef PLASMA_HAS_OMPSS
              if(runtime)
                #pragma omp taskwait
              else
#endif
                PLASMA_Sequence_Wait(sequence);
            }
     
            alpha = 1.0;
            beta = 0.0;
            result = PLASMA_dgemm_Tile_Async(PlasmaNoTrans, PlasmaTrans, 
                                             alpha, descTmp, descDx, 
                                             beta, descCvv, 
                                             sequence, &request[16]);
            HANDLE_ERROR(result, "PLASMA_dgemm_Tile_Async");
            flops = flops + FLOPS_DGEMM(nact, nact, nmeasts);
     
            if(sync){
#ifdef PLASMA_HAS_OMPSS
              if(runtime)
                #pragma omp taskwait
              else
#endif
                PLASMA_Sequence_Wait(sequence);
#ifdef CHECK_PLASMA
STOP_TIMING(t_CeeCvv);
#endif
            }



#ifdef CHECK_PLASMA
    fprintf(stderr,"CHECK MATRICES CEE & CVV against yorick\n");
    double *Ytmp=(double*)malloc(nmeasts*nmeasts*sizeof(double));
    concatenateFileName(dataFile,DATA_FILES_PATH,"cee",suffix);
    readFits(dataFile,Ytmp);
    p_tile_to_lapack(*descCee, lTmp,0);
    compareMatrices2(lTmp,Ytmp,nmeasts,nmeasts,'A','T',&maxerr,&errv1,&errv2);
fprintf(logFile,"Cee      %e   %e     %e   %e\n",t_CeeCvv,maxerr,errv1,errv2);

    concatenateFileName(dataFile,DATA_FILES_PATH,"cvv",suffix);
    readFits(dataFile,Ytmp);
    p_tile_to_lapack(*descCvv, lTmp,0);
    compareMatrices2(lTmp,Ytmp,nact,nact,'A','T',&maxerr,&errv1,&errv2);
fprintf(logFile,"Cvv      %e   %e     %e   %e\n",t_CeeCvv,maxerr,errv1,errv2);
    fprintf(stderr,"END CHECK MATRICES CEE & CVV against yorick\n\n");

    free(lTmp);
    free(Ytmp);

#endif
     
            // Call to intersample
#ifdef USE_INTERSAMPLE
//START_TIMING();
            ao_pdtile_to_lapack_quark(descCvv, Cvv_lap, nact,
                                      sequence, &request[17]);
            /*PLASMA_Sequence_Wait(sequence);
            PLASMA_dgecfi_Async( nact, nact, Cvv,
                                PlasmaCCRB, ts, ts,
                                PlasmaCM, ts, ts,
                                sequence, &request[17] );*/
            PLASMA_Sequence_Wait(sequence);
            //PLASMA_Sequence_Wait(sequence_intersample);
            //char file_name[100] = "Cvv.fits";
            //fits_save_double_2d_square(nact, Cvv_lap, file_name);
            START_TIMING(intersample_time);
            result = PLASMA_Intersample_genPSF( Cvv_lap,
                                                &isample, nact, nbobs, 1, 1,
                                                sequence, &request[18]);


            flops = flops + FLOPS_Intersample(nact, nact);
            if(sync){
#ifdef PLASMA_HAS_OMPSS
              if(runtime)
                #pragma omp taskwait
              else
#endif
                PLASMA_Sequence_Wait(sequence);
            }
//cumulate all the psf
            for(i=0;i<isample.N*isample.N;i++){
                cumulatedPSF[i]+=isample.dphi[i];
            }
            STOP_TIMING(intersample_time);
            total_intersample_time += intersample_time;
//STOP_TIMING();
#endif
            //TODO ALi says: this is fishy, since we need to sync before reporting time (esp. for the async version)
            STOP_TIMING(obs_time);
            total_obs_time += obs_time;
#if defined USE_INTERSAMPLE && defined CHECK_PLASMA
    double *Ypsf=(double*)calloc((size_t)isample.N*isample.N, sizeof(double));
    concatenateFileName(dataFile,DATA_FILES_PATH,"psf",suffix);
    FILE *FileExists=fopen(dataFile,"r");
    if(FileExists!=NULL){
        fclose(FileExists);
        readFits(dataFile,Ypsf);
    }

    compareMatrices2(isample.dphi,Ypsf,isample.N,isample.N,'A','T',&maxerr,&errv1,&errv2);
fprintf(logFile,"psf      %e   %e     %e   %e\n",intersample_time,maxerr,errv1,errv2);
    free(Ypsf);

#endif

       }// End of OBS loop
  }// End of REFINE loop

#ifdef PLASMA_HAS_OMPSS
  if(runtime)
    #pragma omp taskwait
  else//Ali: need to wait to report proper time
#endif

#ifdef USE_INTERSAMPLE
//normalize the cumulated psf
    for(i=0;i<isample.N*isample.N;i++){
       cumulatedPSF[i]/=(maxobs*nbrefine);
    }
    char psf_file[256];
    sprintf(psf_file, "psf_nbrefine%d_maxobs%d.fits", nbrefine,maxobs);
    writeFits(psf_file,isample.N,isample.N,cumulatedPSF);
#endif



    PLASMA_Sequence_Wait(sequence);

  STOP_TIMING(compute_time);
    
  total_time += compute_time;
  
  PLASMA_Sequence_Destroy(sequence);
  
  PLASMA_Desc_Destroy( &descCmm );
  PLASMA_Desc_Destroy( &descCpm );
  PLASMA_Desc_Destroy( &descCpp );
  PLASMA_Desc_Destroy( &descR );
  PLASMA_Desc_Destroy( &descCee );
  PLASMA_Desc_Destroy( &descCvv );
  PLASMA_Desc_Destroy( &descDx );
  PLASMA_Desc_Destroy( &descTmp );
  
#ifdef USE_MATCOV_TILED
  matcov_free_tomo_tiled(&tomo);
#endif

#ifdef USE_INTERSAMPLE
  intersample_free(&isample);
#endif

  if(result != PLASMA_SUCCESS) {
    fprintf(stderr,"An error occured! sorry... %d\n\n", result);
  }
  else { 
     if(verbose) {
        fprintf(stderr, "Done in %f(s)\n", compute_time);
        fprintf(stderr, "Computation successful :-)\n\n");
     }
     fprintf(stderr," Cores  Tile  Refine  Obs  Total   Total  Total MOAO alloc MOAO preprocess  Matcov   MOAO refine    MOAO obs    MOAO PSF     MOAO compute |   Total    |  Gflops/s\n");
     fprintf(stderr," #cores size  #iter  #iter #Meas #Meas-ts #Actu  time (s)   time (s)        time (s)  time (s)      time(s)     time(s)        time (s)   |   time(s)  |  Perf\n");
     fprintf(stderr," ====== ====  ====== ===== ===== ======== ===== ========== =============== ========= =========== ============  ==========   ============= | ========== | ========\n");
     fprintf(stderr,"    %d   %d    %d     %d    %-8d%-8d%-4d    %-11.3f %-10.3f   %-10.3f%-10.3f   %-10.3f   %-10.3f      %-4.3f        %-10.3f   %-4.3f\n\n", ncores, ts, maxrefine, maxobs, nmeas, nmeasts, nact, alloc_cpu_time, preprocess_cpu_time, total_matcov_time, total_refine_time, total_obs_time, total_intersample_time, compute_time, total_time, flops / 1e9 / total_time);
  }



  PLASMA_Finalize();
  if(verbose) fprintf(stderr, "PLASMA: done...\n\n" );
  return 0;
}

#ifdef USE_MATCOV_TILED
//===========================================================================================
void CORE_dGenCovMat_quark(Quark *quark)
{
  struct tomo_struct *tomo;
  int m;
  int n;
  double *A;
  int lda;
  int m0;
  int n0;
  int bigM;
  int part;

  quark_unpack_args_9( quark, m, n, A, lda, m0, n0, bigM, tomo, part);
  matcov_comp_tile(
    A, m, n, m0, n0, lda,
    tomo, part,COLMAJOR);
  /*if(m0 == n0 && part == 1)
  {
    int spd;
    makeSPD(A, m, n, lda, bigM);
    / *
    spd = isSPD(A, m, n, lda);
    if(spd != 0){
      fprintf(stderr,"tile %d,%d is not SPD, returned %d!\n", m0/ts, n0/ts, spd);
    }* /
  }*/
}
//#define DAG_CORE_GenCovMat  DAG_SET_PROPERTIES( "MatCov" , "white"   )
void QUARK_CORE_dGenCovMat( Quark *quark, Quark_Task_Flags *task_flags,
                        int m, int n, double *A, int lda, int m0, int n0, int bigM,
                        struct tomo_struct *tomo, int part)
{
  //DAG_CORE_GenCovMat;
  QUARK_Insert_Task(quark, CORE_dGenCovMat_quark, task_flags,
                    sizeof(int),                      &m,    VALUE,
                    sizeof(int),                      &n,    VALUE,
                    sizeof(double)*lda*n,             A,     OUTPUT | LOCALITY,
                    sizeof(int),                      &lda,  VALUE,
                    sizeof(int),                      &m0,   VALUE,
                    sizeof(int),                      &n0,   VALUE,
                    sizeof(int),                      &bigM, VALUE,
                    sizeof(struct tomo_struct*),      &tomo, VALUE,
                    sizeof(int),                      &part, VALUE,
                    0);
}

int PLASMA_GenCovMat_Tile_Async(PLASMA_desc     *descA,
                                PLASMA_sequence *sequence,
                                PLASMA_request  *request,
                                struct tomo_struct *tomo,
                                int part) {

    Quark *quark;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    int m, n;
    int ldam,ldan;
    int tempnn, tempmm;
    PLASMA_desc A = *descA;

    PLASMA_Get_Quark(&quark);
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
    
    for (m = 0; m < A.mt; m++) {
      tempmm = m == A.mt-1 ? A.m-m*A.mb : A.mb;
      ldam = BLKLDD(A, m);
      
      for (n = 0; n < A.nt; n++) {
        //generate the Lower and diagonal tiles if symmetric
        if(part == 1 && n > m) break;
        ldan=BLKLDD2(A, n);

        tempnn = n == A.nt-1 ? A.n-n*A.nb : A.nb;
        
        QUARK_CORE_dGenCovMat(quark, &task_flags,
                              tempmm, tempnn, A(m, n), ldam,
                              m*A.mb, n*A.nb, A.lm,
                              tomo, part);
      }
    }
    return PLASMA_SUCCESS;
}
//===========================================================================================
void CORE_dUpdateMatCov_quark(Quark *quark)
{
  struct tomo_struct *tomo;

  quark_unpack_args_1( quark, tomo);
  matcov_update_tomo_tiled(tomo);
}
int PLASMA_UpdateMatCov_Async(PLASMA_desc     *descA,
                              PLASMA_sequence *sequence,
                              PLASMA_request  *request,
                                struct tomo_struct *tomo) {

    Quark *quark;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    PLASMA_desc A = *descA;

    PLASMA_Get_Quark(&quark);
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    //DAG_CORE_GenCovMat;
    QUARK_Insert_Task(quark, CORE_dUpdateMatCov_quark, &task_flags,
                      sizeof(struct tomo_struct*),         &tomo, INOUT,
                      0);
    return PLASMA_SUCCESS;
}
#endif
//===========================================================================================
#ifdef USE_INTERSAMPLE
void CORE_dIntersample_genPSF(Quark *quark)
{
  char isample_output_filename[256];
  struct isample_struct *isample;
  double *A;
  int iloop, oloop, gal;

  quark_unpack_args_5( quark, isample, A, iloop, oloop, gal);
  
if(
  intersample_process(isample, A) ==1){
  sprintf(isample_output_filename, "psf_oloop%d_iloop%d_gal%d.fits", oloop, iloop, gal);
  intersample_save(isample, isample_output_filename);
}
}
int PLASMA_Intersample_genPSF(double     *A,
                              struct isample_struct *isample,
                              int nact, int iloop, int oloop, int gal,
                              PLASMA_sequence *sequence,
                              PLASMA_request  *request
                             ) {

    Quark *quark;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    //PLASMA_desc A = *descA;

    PLASMA_Get_Quark(&quark);
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    //DAG_CORE_GenCovMat;
    QUARK_Insert_Task(quark, CORE_dIntersample_genPSF, &task_flags,
                      sizeof(struct isample_struct),  isample, INPUT,
                      sizeof(double)*nact*nact,        A,        INPUT,
                      sizeof(int),                     &iloop,   VALUE,
                      sizeof(int),                     &oloop,   VALUE,
                      sizeof(int),                     &gal,     VALUE,
                      0);
    return PLASMA_SUCCESS;
}

void AO_dlacpy_quark(Quark *quark)
{
  PLASMA_enum uplo;
  int M;
  int N;
  const double *A;
  int LDA;
  double *B;
  int LDB;
  
  quark_unpack_args_7(quark, uplo, M, N, A, LDA, B, LDB);
  LAPACKE_dlacpy_work(
    LAPACK_COL_MAJOR,
    'A',
                      M, N, A, LDA, B, LDB);
}

void AO_CORE_dlacpy(Quark *quark, Quark_Task_Flags *task_flags,
                        PLASMA_enum uplo, int m, int n, int nb,
                        const double *A, int lda,
                        double *B, int ldb)
{
  QUARK_Insert_Task(quark, AO_dlacpy_quark, task_flags,
                    sizeof(PLASMA_enum),                &uplo,  VALUE,
                    sizeof(int),                        &m,     VALUE,
                    sizeof(int),                        &n,     VALUE,
                    sizeof(double)*nb*nb,    A,             INPUT,
                    sizeof(int),                        &lda,   VALUE,
                    sizeof(double)*nb*nb,    B,             OUTPUT,
                    sizeof(int),                        &ldb,   VALUE,
                    0);
}
  
#define AF77(m, n) &(Af77[ ((int64_t)A.nb*(int64_t)lda*(int64_t)(n)) + (int64_t)(A.mb*(m)) ])
#define ABDL(m, n) BLKADDR(A, double, m, n)
void ao_pdtile_to_lapack_quark(PLASMA_desc *descA, double *Af77, int lda,
                                   PLASMA_sequence *sequence, PLASMA_request *request)
{
    double *f77;
    double *bdl;
    int X1, Y1;
    int X2, Y2;
    int n, m, ldt;
    Quark *quark;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    PLASMA_desc A = *descA;

    PLASMA_Get_Quark(&quark);
    if (sequence->status != PLASMA_SUCCESS)
        return;
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for (m = 0; m < A.mt; m++)
    {
        ldt = BLKLDD(A, m);
        for (n = 0; n < A.nt; n++)
        {
            X1 = n == 0 ? A.j%A.nb : 0;
            Y1 = m == 0 ? A.i%A.mb : 0;
            X2 = n == A.nt-1 ? (A.j+A.n-1)%A.nb+1 : A.nb;
            Y2 = m == A.mt-1 ? (A.i+A.m-1)%A.mb+1 : A.mb;

            f77 = AF77(m, n);
            bdl = (double *)plasma_getaddr(A, m, n);
            //bdl = ABDL(m, n);
            AO_CORE_dlacpy(
                quark, &task_flags,
                PlasmaUpperLower, (Y2-Y1), (X2-X1), A.mb,
                &(bdl[X1*lda+Y1]), ldt,
                &(f77[X1*lda+Y1]), lda);
        }
    }
}
#endif


void CORE_dscaldiag_quark(Quark *quark)
{
  int M, LDA;
  double *A, alpha;
  
  quark_unpack_args_4(quark, alpha, M, A, LDA);
  cblas_dscal( M, (alpha), A, LDA+1 );
}


int PLASMA_dscaldiag_Tile_Async(double alpha, PLASMA_desc *descA, PLASMA_sequence *sequence, PLASMA_request  *request) {

    Quark *quark;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    int k;
    int ldak;
    int tempkk;
    PLASMA_desc A = *descA;

    PLASMA_Get_Quark(&quark);
    QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);

    for (k = 0; k < A.nt; k++) {
       tempkk = k == A.nt-1 ? A.m-k*A.mb : A.mb;
       ldak = BLKLDD(A, k);
       QUARK_Insert_Task(quark, CORE_dscaldiag_quark, &task_flags,
                         sizeof(double),                  &alpha,   VALUE,
                         sizeof(int),                     &tempkk,  VALUE,
                         sizeof(double)*tempkk*tempkk,    A(k,k),   INOUT,
                         sizeof(int),                     &ldak,    VALUE,
                         0);
    }
    return PLASMA_SUCCESS;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <mpi.h>
#include <plasma.h>
#include <lapacke.h>
#include "flops.h"
#include "myscalapack.h"
#ifdef USE_MATCOV_TILED
#include "matcov_tiled.h"
#endif
#ifdef USE_INTERSAMPLE
#include "intersample.h"
#endif

//helper functions && globals ========================================
static int
startswith(const char *s, const char *prefix) {
  size_t n = strlen( prefix );
  if (strncmp( s, prefix, n ))
    return 0;
  return 1;
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

void ao_pdlapack_to_tile_quark(double *Af77, int lda, PLASMA_desc *descA,
                                   PLASMA_sequence *sequence, PLASMA_request *request);

#define START_TIMING(_t)               \
_t =- MPI_Wtime();

#define STOP_TIMING(_t)                \
_t += MPI_Wtime();                      \

#ifndef USE_MATCOV_TILED
#define USAGE if (myrank_mpi == 0) fprintf(stderr, "usage: --p=nprow --q=npcol --nb=nb_size --tile=tile_size --nmeas=nmeas --nmeasts=nbmeasts --nact=nbact --maxrefine=maxrefine --maxobs=maxobs [--v] [--s] [--h]\n"\
"--nprow: number of grid row MPI Processes\n" \
"--npcol: number of grid col MPI Processes\n" \
"--nb: Block size\n" \
"--tile: tunable tile size\n" \
"--nmeas: number of measurements\n" \
"--nmeasts: number of measurements in the true sensor\n" \
"--nact: number of actuators\n" \
"--ngal: number of galaxies\n" \
"--maxrefine: max number of refinements\n" \
"--maxobs: max number of observations\n" \
"--v: verbose\n"\
"--s: synchronize between stages, default off\n"\
"--h: print this help\n"\
);
#else
#define USAGE if (myrank_mpi == 0) fprintf(stderr, "usage: --p=nprow --q=npcol --nb=nb_size --tile=tile_size --nssp=nssp --nact=nbact --maxrefine=maxrefine --maxobs=maxobs [--v] [--s] [--h]\n"\
"--nprow: number of grid row MPI Processes\n" \
"--npcol: number of grid col MPI Processes\n" \
"--nb: Block size\n" \
"--tile: tunable tile size\n" \
"--nssp: number of measurements\n" \
"--nact: number of actuators\n" \
"--ngal: number of galaxies\n" \
"--maxrefine: max number of refinements\n" \
"--maxobs: max number of observations\n" \
"--v: verbose\n"\
"--s: synchronize between stages, default off\n"\
"--h: print this help\n"\
);
#endif

#define HANDLE_ERROR(_result, _func)               \
  if(0 != _result) {                  \
    fprintf(stderr,"An error occured (%d), when calling %s at line: %d... \n\n", _result, _func, __LINE__ -1);   \
    MPI_Barrier(MPI_COMM_WORLD);      \
    MPI_Finalize();     \
    return _result;                                         \
  }
  
#define BLKLDD(A, k) ( ( (k) + (A).i/(A).mb) < (A).lm1 ? (A).mb : (A).lm%(A).mb )
#define A(m,n) (double *)plasma_getaddr(A, m, n)
#define my_imax(m,n) ( (m) > (n) ? (m) : (n) )

#define DATA_FILES_PATH "../datafiles/"

//main ===========================================================
int main( int argc, char** argv)
{
  int ncores = 32;
  int nprow = -1, npcol = -1;
  int myrank_mpi, nprocs_mpi;
  int maxobs = 1, maxrefine = 1;
#ifdef USE_MATCOV_TILED
  //struct tomo_struct *gal_tomos;
  struct tomo_struct plasma_tomo;
  struct tomo_struct tomo;
#endif
#ifdef USE_INTERSAMPLE
  struct isample_struct isample;
#endif
  double alpha, beta;
  int nact=0, nmeas=0, nmeasts=0, ts=0, nb=0, verbose=0, sync=0, printhelp=0, unknown=0, ngal = 0;
#ifdef USE_MATCOV_TILED
  int nssp = 30;
#endif
  double flops=0.0; 
  double mc_flops=0.0;
  double sl_flops=0.0;
  double pl_flops=0.0;
/**/
  int nights_obs = 1;
  int hrs_per_night = 1;
  int snapshot_per_hr = 1;
  int snapshot_per_night = hrs_per_night * snapshot_per_hr;
  int total_snapshots = snapshot_per_night * nights_obs;
  MPI_Comm newcomm_snapshot;
  int myrank_s, nprocs_s;
  int nprocs_per_snapshot;
  int color;
  int result;
/**/

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);

  if (verbose & (myrank_mpi == 0)) fprintf(stderr, "ScaLAPACK & PLASMA: MOAO starts... \n");

  if(argc < 2)
  {
    USAGE;
    return 0;
  }
 
  //parse options
  if (verbose & (myrank_mpi == 0)) fprintf(stderr, "Checking arguments...");
  int argi;
  for(argi = 1; argi < argc; argi++)
  {
    if(startswith( argv[argi], "--nprow=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &nprow );
    }
    else
    if(startswith( argv[argi], "--npcol=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &npcol );
    }
    else
    if(startswith( argv[argi], "--nact=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &nact );
      if(nact <= 0){
        if (myrank_mpi == 0) fprintf(stderr, "Invalid number of actuators\n"); MPI_Finalize(); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--ngal=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &ngal );
      if(ngal <= 0){
        if (myrank_mpi == 0) fprintf(stderr, "Invalid number of galaxies\n"); MPI_Finalize(); return 0;
      }
    }
#ifndef USE_MATCOV_TILED
    else
    if(startswith( argv[argi], "--nmeas=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &nmeas );
      if(nmeas <= 0){
        if (myrank_mpi == 0) fprintf(stderr, "Invalid number of total measurements\n"); MPI_Finalize(); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--nmeasts=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &nmeasts );
      if(nmeasts <= 0){
        if (myrank_mpi == 0) fprintf(stderr, "Invalid number of measurements of the true sensor\n"); MPI_Finalize(); return 0;
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
        if (myrank_mpi == 0) fprintf(stderr, "Invalid number of max obs\n"); MPI_Finalize(); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--maxrefine=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &maxrefine );
      if(maxrefine <= 0){
        if (myrank_mpi == 0) fprintf(stderr, "Invalid number of max refinement\n"); MPI_Finalize(); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--nb=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &nb );
      if(nb <= 0){
        if (myrank_mpi == 0) fprintf(stderr, "Invalid tile size\n"); MPI_Finalize(); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--tile=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &ts );
      if(ts <= 0){
        if (myrank_mpi == 0) fprintf(stderr,"Invalid tile size\n"); MPI_Finalize(); return 0;
      }
    }
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
      if (myrank_mpi == 0) fprintf(stderr, "Unknown option %s, aborting...\n", argv[argi]); unknown=1;
    }
  }

  if(printhelp||unknown)
  {
    USAGE;
    MPI_Finalize(); return 0;
  }
  
  /*
  if(nmeas <= 0)
  {
    if (myrank_mpi == 0) fprintf(stderr, "Please provide number of measurements\n");
    USAGE;
    MPI_Finalize(); return 0;
  }
  if(nmeasts <= 0)
  {
    if (myrank_mpi == 0) fprintf(stderr, "Please provide number of measurements in the true sensor\n");
    USAGE;
    MPI_Finalize(); return 0;
  }
  */
  if(nact <= 0)
  {
    if (myrank_mpi == 0) fprintf(stderr, "Please provide number of actuators\n");
    USAGE;
    MPI_Finalize(); return 0;
  }
  if(ngal <= 0)
  {
    if (myrank_mpi == 0) fprintf(stderr, "Please provide number of galaxies/directions\n");
    USAGE;
    MPI_Finalize(); return 0;
  }
  if(maxrefine <= 0)
  {
    if (myrank_mpi == 0) fprintf(stderr, "Please provide number of max refine\n");
    USAGE;
    MPI_Finalize(); return 0;
  }
  if(maxobs < 0)
  {
    if (myrank_mpi == 0) fprintf(stderr, "Please provide number of max obs\n");
    USAGE;
    MPI_Finalize(); return 0;
  }
  if(nb <= 0)
  {
    if (myrank_mpi == 0) fprintf(stderr, "Please provide block size\n");
    USAGE;
    MPI_Finalize(); return 0;
  }
  if(ts <= 0)
  {
    if (myrank_mpi == 0) fprintf(stderr,"Please provide tile size\n");
    USAGE;
    MPI_Finalize(); return 0;
  }

  //defaults
  if(nprow < 0 || npcol < 0)
  {
    nprow = 1;
    npcol = 1;
  }
  if (verbose & (myrank_mpi == 0)) fprintf(stderr, " Done\n");

  int ictxt, ictxtsys;
  int myrow, mycol;
  int i0 = 0, i1 = 1;
  int i, j, k;
  int iseed, info;
  int mloc_nmeas, nloc_nmeas;
  int mloc_nmeasts, nloc_nmeasts;
  int mloc_nact, nloc_nact;
  int mloc_n_totnmeasts;
  double total_time = 0.0; 
  int nbobs, nbrefine;
  // time break down
  double alloc_cpu_time = 0.0, alloc_cpu =0.0;
  double preprocess_cpu_time = 0.0, preprocess_cpu = 0.0;
  double compute_cpu_time = 0.0, compute_cpu =0.0;
  double total_refine_time = 0.0, refine_cpu_time = 0.0, refine_cpu = 0.0;
  double total_comm_time = 0.0, comm_cpu_time = 0.0, comm_cpu = 0.0;
  double total_obs_time = 0.0, obs_time = 0.0, obs_cpu_time = 0.0, obs_cpu = 0.0;
  double total_intersample_time = 0.0, intersample_time = 0.0;
  double total_matcov_time = 0.0, matcov_time = 0.0;
  
  double *Cmm = NULL, *Cpm = NULL, *Cpp = NULL;
  double *Rglob = NULL, *R = NULL, *Cee = NULL, *Cvv = NULL, *Dx = NULL, *Tmp = NULL;
  int descCmm[9], descCpm[9], descCpp[9], descTmp[9];
  int descRglob[9], descR[9], descCee[9], descCvv[9], descDx[9];

  nprocs_per_snapshot = nprocs_mpi/total_snapshots;
  color               = myrank_mpi/nprocs_per_snapshot;
  MPI_Comm_split(MPI_COMM_WORLD, color, myrank_mpi, &newcomm_snapshot);

  MPI_Comm_size(newcomm_snapshot, &nprocs_s);
  MPI_Comm_rank(newcomm_snapshot, &myrank_s);
  int *newcomm_snapshot_procs = (int *)malloc(nprocs_s*sizeof(int)) ;
  MPI_Allgather(&myrank_mpi, 1, MPI_INT, newcomm_snapshot_procs, 1, MPI_INT, newcomm_snapshot);

  //fprintf(stderr, "my global rank is %d and my local rank is %d\n", myrank_mpi, myrank_s);
  //if (myrank_s == 0)
     //for (i = 0; i < nprocs_s; i++)
         //fprintf(stderr, "my global rank is %d and my local rank is %d: %d\n", myrank_mpi, myrank_s, newcomm_snapshot_procs[i]);

  // Initialize BLACS business
  //nprow = sqrt(nprocs_s);
  npcol = nprocs_s/nprow;
  int *imap = (int *)malloc(nprow*npcol*sizeof(int));
  k = 0;
  for (i = 0; i < nprow; i++){
     for (j = 0; j < npcol; j++){
	 *(imap + i + j * nprow) = newcomm_snapshot_procs[k];
	 k = k + 1;
     }
  }
  //MPI_Barrier(MPI_COMM_WORLD);
  //if (myrank_mpi==0) fprintf(stderr, "==========================> nprow %d npcol %d nprocs_s %d\n", nprow, npcol, nprocs_s);
  //MPI_Finalize(); return 0;
  //if (myrank_s == 0)
     //for (i = 0; i < nprocs_per_snapshot; i++)
         //fprintf(stderr, "my global rank is %d and my local rank is %d: %d %d %d\n", myrank_mpi, myrank_s, newcomm_snapshot_procs[i], nprow, npcol);


#ifdef USE_MATCOV_TILED
  /*gal_tomos = (struct tomo_struct*)malloc( ngal * sizeof(struct tomo_struct));
  if(!gal_tomos){
    if (myrank_s == 0) fprintf(stderr, "Out of Memory for gal_tomos\n");
    MPI_Finalize(); return -1;
  }*/
  //double alphaX[16] = {0.0,3.0,41.0,-12.0,71.0,4.0,66.0,-61.0,-1.0,-46.0,-8,75.0,12.0,-37.0,68.0,-74.0};
  //double alphaY[16] = {0.0,69.0,-43.0,-71.0,-38.0,-57.0,-68.0,-53.0,62.0,75.0,-46.0,44.0,32.0,2.0,55.0,-14.0};
  double alphaX[32] ={5,-5,0,0,0,0,0,0,35,0,-35,0,35,0,-35,0,53,75,53,0,-53,-75,-53,0,53,75,53,0,-53,-75,-53,0};
  double alphaY[32] = {0,0,-5,5,-15,15,-15,15,0,-35,0,35,0,-35,0,35,53,0,-53,-75,-53,0,53,75,53,0,-53,-75,-53,0,53,75};

  int gal_index = myrank_s / ncores;
  int night_idx = myrank_mpi / (nprocs_s * snapshot_per_night),
        snapshot_idx = (myrank_mpi / nprocs_s) % snapshot_per_night,
        obs_idx = nbobs, t;
  /*for(t = 0; t < ngal; t++){
    result = !matcov_init_tomo_tiled(&gal_tomos[t], nssp, DATA_FILES_PATH, night_idx, snapshot_per_night, snapshot_idx, -1, alphaX[t], alphaY[t]);
    HANDLE_ERROR(result, "matcov_init_tomo_tiled");
  }*/
  
  result = !matcov_init_tomo_tiled(&tomo, nssp, DATA_FILES_PATH, night_idx, snapshot_per_night, snapshot_idx, -1, 0., 0.);
  HANDLE_ERROR(result, "matcov_init_tomo_tiled");
  nmeas = matcov_getNumMeasurements(&tomo);
  nmeasts = matcov_getNumMeasurementsTS(&tomo);

  
  if ( (myrank_s % ncores) == 0 ) {
    
    //need to set galaxy coordinates on tomo
    //matcov_set_gal_coords(&plasma_tomo, alphaX[gal_index], alphaY[gal_index]);

    result = !matcov_init_tomo_tiled(&plasma_tomo, nssp, DATA_FILES_PATH, night_idx, snapshot_per_night, snapshot_idx, -1, alphaX[gal_index], alphaY[gal_index]);
    HANDLE_ERROR(result, "matcov_init_tomo_tiled");
  }
  //nact = nmeasts;
#endif


#ifdef USE_INTERSAMPLE
  if ( (myrank_s % ncores) == 0 ) {
    intersample_prepare(&isample, nact*nact, DATA_FILES_PATH);
    intersample_init(&isample);
  }
#endif



  if (verbose & (myrank_s == 0)) fprintf(stderr, "BLACS Init...");
  Cblacs_get(0, 0, &ictxt);
  Cblacs_gridmap(&ictxt, imap, nprow, nprow, npcol);
  Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
  if (verbose & (myrank_s == 0)) fprintf(stderr, " Done\n");
  if (verbose & (myrank_s == 0)) fprintf(stderr, "GLOBAL RANK: %d nprow %d npcol %d myrow %d mycol %d\n", myrank_mpi, nprow, npcol, myrow, mycol);

  //MPI_Barrier(MPI_COMM_WORLD);
  //fprintf(stderr, "R %d r %d nprow %d npcol %d\n", myrank_mpi, myrank_s, nprow, npcol);
  //return 0;


  if (myrank_s == 0) {
     fprintf(stderr, "# \n");
     fprintf(stderr, "# nmeas %d nmeasts %d nact %d nb %d\n", nmeas, nmeasts, nact, nb);
     fprintf(stderr, "# nprocs %d P %d Q %d\n", nprocs_s, nprow, npcol);
     //fprintf(stderr, "# niter %d\n", niter);
     //fprintf(stderr, "# n_range %d:%d:%d mode: %d cond: %2.4e \n", start, stop, step, mode, cond);
     fprintf(stderr, "# \n");
  }

  int n_totnmeasts  = ngal * nmeasts;
  mloc_nmeas        = numroc_( &nmeas,        &nb, &myrow, &i0, &nprow );
  nloc_nmeas        = numroc_( &nmeas,        &nb, &mycol, &i0, &npcol );
  mloc_nmeasts      = numroc_( &nmeasts,      &nb, &myrow, &i0, &nprow );
  nloc_nmeasts      = numroc_( &nmeasts,      &nb, &mycol, &i0, &npcol );
  mloc_nact         = numroc_( &nact,         &nb, &myrow, &i0, &nprow );
  nloc_nact         = numroc_( &nact,         &nb, &mycol, &i0, &npcol );
  mloc_n_totnmeasts = numroc_( &n_totnmeasts, &nb, &myrow, &i0, &nprow );

  /*/initialize galaxy indices
  double *alphaX = NULL, *alphaY = NULL;
  alphaX = (double*)malloc(ngal * sizeof(double));
  alphaY = (double*)malloc(ngal * sizeof(double));
  if ( ! alphaX || !alphaY ) {
    if (myrank_s == 0) fprintf(stderr, "Out of Memory for alphaY or alphaX\n");
    MPI_Finalize(); return -1;
  }
  srand(time(NULL));
  int g;
  for(g = 0; g < ngal; g++){
    //TODO make sure this value is the same for all procs participating in one field of view
    alphaX[g] = rand() % 100;
    alphaY[g] = rand() % 100;
  }*/
  

  if (verbose & (myrank_s == 0)) fprintf(stderr, "mloc_nmeas %d nloc_nmeas %d mloc_nmeasts %d nloc_nmeasts %d mloc_nact %d nloc_nact %d\n", mloc_nmeas, nloc_nmeas, mloc_nmeasts, nloc_nmeasts, mloc_nact, nloc_nact);


  if (verbose & (myrank_s == 0)) fprintf(stderr, "Desc Init...");
  descinit_( descCmm,   &nmeas,        &nmeas, &nb, &nb, &i0, &i0, &ictxt, &mloc_nmeas,        &info );
  descinit_( descRglob, &n_totnmeasts, &nmeas, &nb, &nb, &i0, &i0, &ictxt, &mloc_n_totnmeasts, &info );
  if (verbose & (myrank_s == 0)) fprintf(stderr, " Done\n");


  // Allocate matrix memory
  {
    if(verbose & (myrank_s == 0)) fprintf(stderr, "Allocating matrix memory...");
    START_TIMING(alloc_cpu_time);
    Cmm = (double*)calloc( (size_t)mloc_nmeas * nloc_nmeas, sizeof(double) );
    if ( ! Cmm ) {
      if (myrank_s == 0) fprintf(stderr, "Out of Memory for Covariance Matrix (Cmm)\n");
      MPI_Finalize(); return -1;
    }
    Rglob = (double*)calloc( (size_t)mloc_n_totnmeasts * nloc_nmeas, sizeof(double) );
    if ( ! Rglob ) {
      if (myrank_s == 0) fprintf(stderr, "Out of Memory for Global ToR (R)\n");
      MPI_Finalize(); return -1;
    }
    STOP_TIMING(alloc_cpu_time);
    MPI_Allreduce( &alloc_cpu_time, &alloc_cpu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    total_time += alloc_cpu;
    if(verbose & (myrank_s == 0)) fprintf(stderr, " Done in %f(s)\n\n", alloc_cpu);
  }
  
    
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                       ScaLAPACK                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if(verbose & (myrank_s == 0)) fprintf(stderr, "Maximum number of refinements %d\n", maxrefine);
  if(verbose & (myrank_s == 0)) fprintf(stderr, "Maximum number of observations %d\n", maxobs);
  if(verbose & (myrank_s == 0)) fprintf(stderr, "ScaLAPACK: MOAO main computation...");

  START_TIMING(compute_cpu_time);
    
  START_TIMING(refine_cpu_time);
  //if(verbose & (myrank_s == 0)) fprintf(stderr, "\nnbrefine %d\n", nbrefine+1);

  //PLASMA_dplgsy_Tile_Async( (double)nmeas, descCmm, rand()%100, sequence, &request[1] );
  //Cmm = (double*)malloc( (size_t)mloc_nmeas * nloc_nmeas * sizeof(double) );
  /*
  */

  START_TIMING(matcov_time);
#ifdef USE_MATCOV_TILED
  //matcov_comp_tile_scalapack(Cmm, nmeas, nmeas, 0, 0, mloc_nmeas, myrow, mycol, nb, nprow, npcol, &tomo, 1);
  matcov_comp_tile_scalapack(Cmm, mloc_nmeas, nloc_nmeas, 0, 0, mloc_nmeas, myrow, mycol, nb, nprow, npcol, &tomo, 1);
  flops = flops + FLOPS_DCovMat(nmeas, nmeas, tomo.Nlayer, 1);
  mc_flops +=FLOPS_DCovMat(nmeas, nmeas, tomo.Nlayer, 1);
#else
  iseed = myrank_mpi * mloc_nmeas * nloc_nmeas;
  pdmatgen_( &ictxt, "S","D", &nmeas, &nmeas, &nb, &nb, Cmm, &mloc_nmeas,
             &i0, &i0, &iseed, 
             &i0, &mloc_nmeas, &i0, &nloc_nmeas, 
             &myrow, &mycol, &nprow, &npcol );
#endif
  STOP_TIMING(matcov_time);
  total_matcov_time += matcov_time;

  //PLASMA_dplrnt_Tile_Async( descR, rand()%100, sequence, &request[2] );
  //R = (double*)malloc( (size_t)mloc_nmeasts * nloc_nmeas * sizeof(double) );
  /*
  */
  START_TIMING(matcov_time);
#ifdef USE_MATCOV_TILED
  //for each galaxy
  int g;
  for(g = 0; g < ngal; g++){
      //set the galaxy coords
      matcov_set_gal_coords(&tomo, alphaX[g], alphaY[g]);
      matcov_comp_tile_scalapack( Rglob + (size_t)mloc_nmeasts * g, 
                                  mloc_nmeasts, nloc_nmeas,
      //matcov_comp_tile_scalapack( Rglob + (size_t)mloc_nmeasts * nloc_nmeas * g, 
                                  //nmeasts, nmeas, 
                                  0, 0, 
                                  mloc_n_totnmeasts, myrow, mycol,
                                  nb, nprow, npcol, 
                                  &tomo, 3);
                                  //&gal_tomos[g], 3);
      //flops = flops + FLOPS_DCovMat(nmeasts, nmeas, gal_tomos[g].Nlayer, 3);
      flops = flops + FLOPS_DCovMat(nmeasts, nmeas, tomo.Nlayer, 3);
      mc_flops += FLOPS_DCovMat(nmeasts, nmeas, tomo.Nlayer, 3);
      //matcov_comp_tile_scalapack(R, nmeasts, nmeas, 0, 0, mloc_nmeasts, myrow, mycol, nb, nprow, npcol, &tomo, 3);
  }
#else
  iseed = myrank_mpi * mloc_n_totnmeasts * nloc_nmeas;
  pdmatgen_( &ictxt, "R", "ND", &n_totnmeasts, &nmeas, &nb, &nb, Rglob, &mloc_n_totnmeasts,
             &i0, &i0, &iseed, 
             &i0, &mloc_nmeasts, &i0, &nloc_nmeas, 
             &myrow, &mycol, &nprow, &npcol );
  flops = flops + FLOPS_DCovMat(nmeasts, nmeas);
  mc_flops += FLOPS_DCovMat(nmeasts, nmeas);
#endif
  STOP_TIMING(matcov_time);
  total_matcov_time += matcov_time;
  if(verbose & (myrank_s == 0)) fprintf(stderr, " Done Matcov flops=%f \n", mc_flops/total_matcov_time);
  //result = PLASMA_dpotrf_Tile_Async(PlasmaLower, descCmm, sequence, &request[3]);
  pdpotrf_( "L", &nmeas, Cmm, &i1, &i1, descCmm, &info );
  HANDLE_ERROR(info, "pdpotrf");
  flops = flops + FLOPS_DPOTRF(nmeas);
  sl_flops += FLOPS_DPOTRF(nmeas);
  alpha = 1.0;
  //result = PLASMA_dtrsm_Tile_Async(PlasmaRight, PlasmaLower, PlasmaTrans, PlasmaNonUnit, 
                                        //alpha, descCmm, descR, sequence, &request[4]);
  pdtrsm_( "R", "L", "T", "N", &n_totnmeasts, &nmeas, &alpha, 
           Cmm, &i1, &i1, descCmm,
           Rglob, &i1, &i1, descRglob );
  flops = flops + FLOPS_DTRSM(PlasmaRight, n_totnmeasts, nmeas);
  sl_flops += FLOPS_DTRSM(PlasmaRight, n_totnmeasts, nmeas);
  alpha = 1.0;
  //result = PLASMA_dtrsm_Tile_Async(PlasmaRight, PlasmaLower, PlasmaNoTrans, PlasmaNonUnit, 
                                     //alpha, descCmm, descR, sequence, &request[5]);
  pdtrsm_( "R", "L", "N", "N", &n_totnmeasts, &nmeas, &alpha,
           Cmm, &i1, &i1, descCmm,
           Rglob, &i1, &i1, descRglob );
  flops = flops + FLOPS_DTRSM(PlasmaRight, n_totnmeasts, nmeas);
  sl_flops += FLOPS_DTRSM(PlasmaRight, n_totnmeasts, nmeas);
  STOP_TIMING(refine_cpu_time);
  MPI_Allreduce( &refine_cpu_time, &refine_cpu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  total_refine_time += refine_cpu;
  if(verbose & (myrank_s == 0)) fprintf(stderr, " Done flops=%f \n", sl_flops/total_refine_time);
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                       GATHER                                                              //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  START_TIMING(comm_cpu_time);
  if(verbose & (myrank_s == 0)) fprintf(stderr, "Gather each single ToR on each single node...");
  MPI_Datatype sltype;
  int gsizes[2], distribs[2], dargs[2], psizes[2];

  gsizes[0] = nmeasts; /* no. of rows in global array */
  gsizes[1] = nmeas; /* no. of columns in global array*/

  distribs[0] = MPI_DISTRIBUTE_CYCLIC;
  distribs[1] = MPI_DISTRIBUTE_CYCLIC;

  dargs[0] = nb; // no of rows in block
  dargs[1] = nb; // no of cols in block

  psizes[0] = nprow; /* no. of processes in vertical dimension
         of process grid */
  psizes[1] = npcol; /* no. of processes in horizontal dimension
         of process grid */

  MPI_Type_create_darray( nprocs_s, myrank_s, 2, gsizes, distribs, dargs, psizes,
                          MPI_ORDER_C, MPI_DOUBLE, &sltype );

  int locsize;
  MPI_Type_commit(&sltype);
  MPI_Type_size(sltype, &locsize);
  int nelements = locsize / sizeof(double);

  //fprintf(stderr, "I AM HERE %d %d %s \n", mloc_nmeasts * nloc_nmeas, nelements, mloc_nmeasts * nloc_nmeas == nelements ? "OK" : "NOT OK");
  R = (double *)malloc(nmeasts*nmeas*sizeof(double)); 
 
  int iroot, root;
  for(iroot=0;iroot<nprocs_s/ncores;iroot++){
      root = iroot *ncores;
      if(myrank_s==root) R = (double *)malloc(gsizes[0] * gsizes[1] * sizeof(double)); 
      int rank;
      int tag = 0;
      MPI_Status status;
    
      MPI_Request sreq;
      int reqcount= nprocs_s -1, count=0;
      MPI_Request rreq[reqcount];
      MPI_Status sstatus;
      MPI_Status rstatus[reqcount];
    
      for(i=0; i<nprow; i++){
          for(j=0; j<npcol; j++){
              rank=i*npcol+j;
              if(( rank== myrank_s)&& (rank != root)) {
                   MPI_Isend(Rglob+iroot*mloc_nmeasts * nloc_nmeas,  mloc_nmeasts * nloc_nmeas, MPI_DOUBLE, root, tag, newcomm_snapshot, &sreq);
              }
              if(myrank_s == root){
                 if (rank == root) {
                    MPI_Sendrecv(Rglob+iroot*mloc_nmeasts * nloc_nmeas, mloc_nmeasts * nloc_nmeas, MPI_DOUBLE, root , tag, R, 1, sltype, root, tag,  newcomm_snapshot, &status);
                 }
                 else {
                    MPI_Irecv(R+nb*j+i*(npcol*nloc_nmeas)*nb -(root%npcol)*nb- (root/npcol)*(npcol*nloc_nmeas)*nb, 1, sltype, rank, tag, newcomm_snapshot, &rreq[count]);
                    count ++;
                 }
              }
          }
      }
      if(myrank_s == root)
         MPI_Waitall(reqcount, rreq, rstatus);
      else 
         MPI_Wait(&sreq, &sstatus);
  }
  MPI_Type_free(&sltype); 

  STOP_TIMING(comm_cpu_time);
  MPI_Allreduce( &comm_cpu_time, &comm_cpu, 1, MPI_DOUBLE, MPI_MAX, newcomm_snapshot);
  total_comm_time += comm_cpu;


  // Let's free some memory for PLASMA
  free( Rglob );
  free( Cmm );

#ifdef USE_MATCOV_TILED
  matcov_free_tomo_tiled(&tomo);
  //for each galaxy
  /*for(g = 0; g < ngal; g++){
    matcov_free_tomo_tiled(&gal_tomos[g]);
  }*/
#endif

  if(verbose & (myrank_s == 0)) fprintf(stderr, " Done\n");
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                       PLASMA                                                              //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
//                                                                                                                           //
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if ( (myrank_s % ncores) == 0 ) {
    
#define AFFINITY 1
#ifndef  AFFINITY
    PLASMA_Init(ncores);
#else
    int c;
    int* core_ids = (int*) malloc((ncores)*sizeof(int));
    for(c = 0; c < ncores; c++)
      core_ids[c] = c;
     PLASMA_Init_Affinity(ncores, core_ids);
#endif
    
     PLASMA_Disable(PLASMA_AUTOTUNING);
     PLASMA_Set(PLASMA_TILE_SIZE, ts );
     //PLASMA_Set(PLASMA_TRANSLATION_MODE, PLASMA_INPLACE);
     PLASMA_Set(PLASMA_SCHEDULING_MODE, PLASMA_DYNAMIC_SCHEDULING);
     PLASMA_Enable(PLASMA_WARNINGS);
     PLASMA_Enable(PLASMA_ERRORS);
 
     int result;

     double *Cmm = NULL, *Cpm = NULL, *Cpp = NULL;
     PLASMA_desc *descCmm = NULL, *descCpm = NULL, *descCpp = NULL;

     double *Cee = NULL, *Cvv = NULL, *Cvv_lap = NULL, *Dx = NULL, *Tmp = NULL;
     PLASMA_desc *descR = NULL, *descCee = NULL, *descCvv = NULL, *descDx = NULL, *descTmp = NULL;
 
     PLASMA_sequence *sequence;
     PLASMA_request request[13] = { PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,
                                    PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,
                                    PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER,PLASMA_REQUEST_INITIALIZER};
     PLASMA_Sequence_Create(&sequence);

     // Allocate matrix memory
     //TODO why not allocating as a preprocess step
     {
       if(verbose & (myrank_s == 0)) fprintf(stderr, "Allocating matrix memory...");
       START_TIMING(alloc_cpu_time);
       Cmm = (double*)calloc( (size_t)nmeas * nmeas, sizeof(double) );
       if ( ! Cmm ) {
         fprintf(stderr, "Out of Memory for Covariance Matrix (Cmm)\n");
         return -1;
       }
#ifdef PLASMA_HAS_OMPSS
       if(runtime) {
          #pragma omp register ([nmeas*nmeas]Cmm)
       }
#endif
       Cpm = (double*)calloc( (size_t)nmeasts * nmeas, sizeof(double) );
       if ( ! Cpm ) {
         fprintf(stderr, "Out of Memory for Cpm Matrix\n");
         return -1;
       }
#ifdef PLASMA_HAS_OMPSS
       if(runtime) {
          #pragma omp register ([nmeasts*nmeas]Cpm)
       }
#endif
       Cpp = (double*)calloc( (size_t)nmeasts * nmeasts, sizeof(double) );
       if ( ! Cpp ) {
         fprintf(stderr, "Out of Memory for Cpp Matrix\n");
         return -1;
       }
#ifdef PLASMA_HAS_OMPSS
       if(runtime) {
          #pragma omp register ([nmeasts*nmeasts]Cpp)
       }
#endif
       Cee = (double*)calloc( (size_t)nmeasts * nmeasts, sizeof(double) );
       if ( ! Cee ) {
         fprintf(stderr, "Out of Memory for Cee Matrix\n");
         return -1;
       }
#ifdef PLASMA_HAS_OMPSS
       if(runtime) {
          #pragma omp register ([nmeasts*nmeasts]Cee)
       }
#endif
       Cvv = (double*)calloc( (size_t)nact * nact, sizeof(double) );
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
       Dx = (double*)calloc( (size_t)nact * nmeasts, sizeof(double) );
       if ( ! Dx ) {
         fprintf(stderr, "Out of Memory for Dx Matrix\n");
         return -1;
       }
#ifdef PLASMA_HAS_OMPSS
       if(runtime) {
          #pragma omp register ([nact*nmeasts]Dx)
       }
#endif
       Tmp = (double*)calloc( (size_t)nact * nmeasts, sizeof(double) );
       if ( ! Tmp ) {
         fprintf(stderr, "Out of Memory for Tmp Matrix\n");
         return -1;
       }
#ifdef PLASMA_HAS_OMPSS
       if(runtime) {
          #pragma omp register ([nact*nmeasts]Tmp)
       }
#endif
       STOP_TIMING(alloc_cpu_time);
       total_time += alloc_cpu_time;
       if(verbose & (myrank_s == 0)) fprintf(stderr, "Done in %f(s)\n\n", alloc_cpu_time);
     }

     //TODO we should use the Async version
     //PLASMA_dgecfi( nmeasts, nmeas, R,
                    //PlasmaCM, ts, ts,
                    //PlasmaCCRB, ts, ts);
     PLASMA_Desc_Create(&descR, R, PlasmaRealDouble, ts, ts, ts*ts, nmeasts, nmeas, 0, 0, nmeasts, nmeas);

     ao_pdlapack_to_tile_quark(R, nmeasts, descR, sequence, &request[0]);

     PLASMA_Desc_Create(&descCmm, Cmm, PlasmaRealDouble, ts, ts, ts*ts, nmeas,   nmeas,   0, 0, nmeas,   nmeas);
     PLASMA_Desc_Create(&descCpm, Cpm, PlasmaRealDouble, ts, ts, ts*ts, nmeasts, nmeas,   0, 0, nmeasts, nmeas);
     PLASMA_Desc_Create(&descCpp, Cpp, PlasmaRealDouble, ts, ts, ts*ts, nmeasts, nmeasts, 0, 0, nmeasts, nmeasts);
     PLASMA_Desc_Create(&descCee, Cee, PlasmaRealDouble, ts, ts, ts*ts, nmeasts, nmeasts, 0, 0, nmeasts, nmeasts);
     PLASMA_Desc_Create(&descCvv, Cvv, PlasmaRealDouble, ts, ts, ts*ts, nact,    nact,    0, 0, nact,    nact);
     PLASMA_Desc_Create(&descDx,  Dx,  PlasmaRealDouble, ts, ts, ts*ts, nact,    nmeasts, 0, 0, nact,    nmeasts);
     PLASMA_Desc_Create(&descTmp, Tmp, PlasmaRealDouble, ts, ts, ts*ts, nact,    nmeasts, 0, 0, nact,    nmeasts);
  
     if(verbose & (myrank_s == 0)) fprintf(stderr, "ts %d nact %d nmeasts %d\n", ts, nact, nmeasts);
     if(verbose & (myrank_s == 0)) fprintf(stderr, "Loading the interaction matrix Dx...");
     START_TIMING(preprocess_cpu_time);
     PLASMA_dplrnt_Tile_Async( descDx, rand()%100, sequence, &request[0] );
     if(sync){
#ifdef PLASMA_HAS_OMPSS
       if(runtime)
         #pragma omp taskwait
       else
#endif
         PLASMA_Sequence_Wait(sequence);
     }
     STOP_TIMING(preprocess_cpu_time);
     if(verbose & (myrank_s == 0)) fprintf(stderr, " Done in %f(s)\n\n", preprocess_cpu_time);

     total_time += preprocess_cpu_time;

     if(verbose & (myrank_s % nprocs_s == 0)) fprintf(stderr, "PLASMA ID %d: MOAO main computation on %d observing sequence using %d cores\n", myrank_s, maxobs, ncores);


     for (nbobs = 0; nbobs < maxobs; nbobs++) {
          START_TIMING(obs_time);
          if(verbose & (myrank_s == 0)) fprintf(stderr, " nbobs %d\n", nbobs+1);
 
          // MatCov generation
          START_TIMING(matcov_time);
#ifdef USE_MATCOV_TILED
          srand(time(NULL));
          //need to get new atm-params
          int night_idx = myrank_mpi / (nprocs_s * snapshot_per_night),
              snapshot_idx = rand() % snapshot_per_night,//(myrank_mpi / nprocs_s) % snapshot_per_night,
              obs_idx = nbobs;
          result = !matcov_update_atm_params(&plasma_tomo, night_idx, snapshot_per_night, snapshot_idx, obs_idx);
          HANDLE_ERROR(result, "matcov_update_atm_params");
          
          result = PLASMA_GenCovMat_Tile_Async( descCmm, sequence, &request[1], &plasma_tomo, 1 );
          HANDLE_ERROR(result, "PLASMA_GenCovMat_Tile_Async");
          flops = flops + 16 * FLOPS_DCovMat(nmeas, nmeas, plasma_tomo.Nlayer, 1);
          mc_flops = mc_flops + 16 * FLOPS_DCovMat(nmeas, nmeas, plasma_tomo.Nlayer, 1);
#else
          result = PLASMA_dplgsy_Tile_Async( (double)nmeas, descCmm, rand()%100, sequence, &request[1] );
          HANDLE_ERROR(result, "PLASMA_dplgsy_Tile_Async");
          flops = flops + 16 * FLOPS_DCovMat(nmeas, nmeas);
          mc_flops = mc_flops + 16 * FLOPS_DCovMat(nmeas, nmeas);

#endif
          if(sync){
#ifdef PLASMA_HAS_OMPSS
            if(runtime)
              #pragma omp taskwait
            else
#endif
              PLASMA_Sequence_Wait(sequence);
          }
          STOP_TIMING(matcov_time);
          total_matcov_time += matcov_time;

          START_TIMING(matcov_time);
#ifdef USE_MATCOV_TILED
          result = PLASMA_GenCovMat_Tile_Async( descCpm, sequence, &request[2], &plasma_tomo, 3 );
          HANDLE_ERROR(result, "PLASMA_GenCovMat_Tile_Async");
          flops = flops + 16 * FLOPS_DCovMat(nmeasts, nmeas, plasma_tomo.Nlayer, 3);
          mc_flops = mc_flops + 16 * FLOPS_DCovMat(nmeasts, nmeas, plasma_tomo.Nlayer, 3);
#else
          result = PLASMA_dplrnt_Tile_Async( descCpm, rand()%100, sequence, &request[2] );
          HANDLE_ERROR(result, "PLASMA_dplrnt_Tile_Async");
          flops = flops + 16 * FLOPS_DCovMat(nmeasts, nmeas);
          mc_flops = mc_flops + 16 * FLOPS_DCovMat(nmeasts, nmeas);
#endif
          if(sync){
#ifdef PLASMA_HAS_OMPSS
            if(runtime)
              #pragma omp taskwait
            else
#endif
            PLASMA_Sequence_Wait(sequence);
          }
          STOP_TIMING(matcov_time);
          total_matcov_time += matcov_time;
        
          //PLASMA_dplgsy_Tile_Async( (double)nmeas, descCpp, rand()%100, sequence, &request[8] );
          START_TIMING(matcov_time);
#ifdef USE_MATCOV_TILED
          result = PLASMA_GenCovMat_Tile_Async( descCpp, sequence, &request[3], &plasma_tomo, 4 );
          HANDLE_ERROR(result, "PLASMA_GenCovMat_Tile_Async");
          mc_flops = mc_flops + 16 * FLOPS_DCovMat(nmeasts, nmeasts, plasma_tomo.Nlayer, 4);
#else
          result = PLASMA_dplgsy_Tile_Async( (double)nmeasts, descCpp, rand()%100, sequence, &request[3] );
          HANDLE_ERROR(result, "PLASMA_dplgsy_Tile_Async");
          flops = flops + 16 * FLOPS_DCovMat(nmeasts, nmeasts);
          mc_flops = mc_flops + 16 * FLOPS_DCovMat(nmeasts, nmeasts, plasma_tomo.Nlayer, 4);
#endif
          if(sync){
#ifdef PLASMA_HAS_OMPSS
            if(runtime)
              #pragma omp taskwait
            else
#endif
            PLASMA_Sequence_Wait(sequence);
          }
          STOP_TIMING(matcov_time);
          total_matcov_time += matcov_time;
        
          alpha = -1.0;
          beta = 1.0;
          result = PLASMA_dsyr2k_Tile_Async( PlasmaLower, PlasmaNoTrans, 
                                             alpha, descCpm, descR, 
                                             beta, descCpp, 
                                             sequence, &request[4] );
          HANDLE_ERROR(result, "PLASMA_dsyr2k_Tile_Async");
          flops = flops + 16 * FLOPS_DSYR2K(nmeasts, nmeas);
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
          result = PLASMA_dsymm_Tile_Async( PlasmaRight, PlasmaLower,
                                            alpha, descCmm, descR,
                                            beta,  descCpm,
                                            sequence, &request[5] );
          // Cpm gets overwritten at this point and it is used as a buffer
          HANDLE_ERROR(result, "PLASMA_dsymm_Tile_Async");
          flops = flops + 16 * FLOPS_DSYMM(PlasmaRight, nmeasts, nmeas);
          if(sync){
#ifdef PLASMA_HAS_OMPSS
            if(runtime)
              #pragma omp taskwait
            else
#endif
            PLASMA_Sequence_Wait(sequence);
          }
             

          result = PLASMA_dlacpy_Tile_Async( PlasmaLower, descCpp, descCee, 
                                             sequence, &request[6] );
          HANDLE_ERROR(result, "PLASMA_dlacpy_Tile_Async");
          if(sync){
#ifdef PLASMA_HAS_OMPSS
            if(runtime)
              #pragma omp taskwait
            else
#endif
            PLASMA_Sequence_Wait(sequence);
          }

          alpha = 0.5;
          result = PLASMA_dgeadd_Tile_Async( PlasmaTrans, alpha, descCee, alpha, descCpp, 
                                             sequence, &request[7] );
          HANDLE_ERROR(result, "PLASMA_dgeadd_Tile_Async");
          if(sync){
#ifdef PLASMA_HAS_OMPSS
            if(runtime)
              #pragma omp taskwait
            else
#endif
            PLASMA_Sequence_Wait(sequence);
          }

          alpha = 1.0;
          beta = 1.0;
          result = PLASMA_dgemm_Tile_Async( PlasmaNoTrans, PlasmaTrans, 
                                            alpha, descCpm, descR, 
                                            beta, descCpp, 
                                            sequence, &request[8] );
          HANDLE_ERROR(result, "PLASMA_dgemm_Tile_Async");
          flops = flops + 16 * FLOPS_DGEMM(nmeasts, nmeasts, nmeas);
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
          result = PLASMA_dsymm_Tile_Async( PlasmaRight, PlasmaLower,
                                            alpha, descCpp, descDx,
                                            beta,  descTmp,
                                            sequence, &request[9] );
          HANDLE_ERROR(result, "PLASMA_dsymm_Tile_Async");
          flops = flops + 16 * FLOPS_DSYMM(PlasmaRight, nact, nmeasts);
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
          result = PLASMA_dgemm_Tile_Async( PlasmaNoTrans, PlasmaTrans, 
                                            alpha, descTmp, descDx, 
                                            beta, descCvv, 
                                            sequence, &request[10] );
          HANDLE_ERROR(result, "PLASMA_dgemm_Tile_Async");
          flops = flops + 16 * FLOPS_DGEMM(nact, nact, nmeasts);
          if(sync){
#ifdef PLASMA_HAS_OMPSS
            if(runtime)
              #pragma omp taskwait
            else
#endif
            PLASMA_Sequence_Wait(sequence);
          }
     
          // Call to intersample
#ifdef USE_INTERSAMPLE
          ao_pdtile_to_lapack_quark(descCvv, Cvv_lap, nact,
                                    sequence, &request[11]);
          PLASMA_Sequence_Wait(sequence);
          START_TIMING(intersample_time);
          int oloop_index = myrank_mpi/nprocs_per_snapshot;
          result = PLASMA_Intersample_genPSF( Cvv_lap,
                                              &isample, nact, nbobs, oloop_index, gal_index,
                                              sequence, &request[12]);
          HANDLE_ERROR(result, "PLASMA_Intersample_genPSF");
          flops = flops + 16 * FLOPS_Intersample(nact, nact);
          if(sync){
#ifdef PLASMA_HAS_OMPSS
            if(runtime)
              #pragma omp taskwait
            else
#endif
            PLASMA_Sequence_Wait(sequence);
          }
          STOP_TIMING(intersample_time);
          total_intersample_time += intersample_time;
#endif
          STOP_TIMING(obs_cpu_time);
          //MPI_Allreduce( &obs_cpu_time, &obs_cpu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
          total_obs_time += obs_cpu;

     }
    PLASMA_Sequence_Wait(sequence);
     if(verbose & (myrank_s == 0)) fprintf(stderr, " Done\n");

     PLASMA_Sequence_Destroy(sequence);
  
     PLASMA_Desc_Destroy( &descCmm );
     PLASMA_Desc_Destroy( &descCpm );
     PLASMA_Desc_Destroy( &descCpp );
     PLASMA_Desc_Destroy( &descR );
     PLASMA_Desc_Destroy( &descCee );
     PLASMA_Desc_Destroy( &descCvv );
     PLASMA_Desc_Destroy( &descDx );
     PLASMA_Desc_Destroy( &descTmp );
  
     free( Cmm );
     free( Cpm );
     free( Cpp );
     free( R );
     free( Cee );
     free( Cvv );
     free( Cvv_lap );
     free( Dx );
     free( Tmp );
#ifdef USE_MATCOV_TILED
     matcov_free_tomo_tiled(&plasma_tomo);
#endif

#ifdef  AFFINITY
    free(core_ids);
#endif

#ifdef USE_INTERSAMPLE
    intersample_free(&isample);
#endif
  }




  STOP_TIMING(compute_cpu_time);
  //MPI_Allreduce( &compute_cpu_time, &compute_cpu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
  total_time += compute_cpu_time;
  
  //MPI_Barrier(MPI_COMM_WORLD); 

  if(info != 0) {
    if(myrank_s % ncores == 0) fprintf(stderr, "An error occured! sorry... %d\n\n", info);
  }
  else { 
     if(verbose & (myrank_s == 0)) fprintf(stderr, "Done in %f(s)\n", compute_cpu);
     if(verbose & (myrank_s == 0)) fprintf(stderr, "Computation successful :-)\n\n");

     if(myrank_s == 0) fprintf(stderr, "  Procs   P    Q    Block  Refine  Obs  Total   Total  Total MOAO alloc MOAO preprocess Matcov  MOAO refine MOAO obs     MOAO compute |   Total    |  Gflops/s\n");
     if(myrank_s == 0) fprintf(stderr, "    #                size  #iter  #iter #Meas #Meas-ts #Actu  time (s)   time (s)       time (s)    time (s)   time(s)   time (s)     |   time(s)  |  Perf\n");
     if(myrank_s == 0) fprintf(stderr, " ======= ==== ===== ====== ====== ===== ===== ======== ===== ========== =============== ========= =========== ======== ============== | ========== | ========\n");
     if(myrank_s == 0) fprintf(stderr, "    %d    %d    %d    %d     %d     %d   %-8d%-8d%-4d  %-11.3f    %-10.3f   %-10.3f%-10.3f%-10.3f %-4.3f            %-10.3f %-4.3f\n\n", nprocs_mpi, nprow, npcol, nb, maxrefine, maxobs, nmeas, nmeasts, nact, alloc_cpu, preprocess_cpu, total_matcov_time, total_refine_time, total_obs_time, compute_cpu_time, total_time, flops / 1e9 / total_time);
     if(myrank_s == 0) fprintf(stderr, "Comm time    %f\n", total_comm_time );
  }

  sleep(20000);
  MPI_Barrier(MPI_COMM_WORLD); 
  Cblacs_gridexit(0);
  
  MPI_Finalize(); return 0;
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
    tomo, part);
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
                                int part)
{
    Quark *quark;
    Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
    int m, n;
    int ldam;
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
  
  intersample_process(isample, A);
  sprintf(isample_output_filename, "psf_oloop%d_iloop%d_gal%d.fits", oloop, iloop, gal);
  intersample_save(isample, isample_output_filename);
}

int PLASMA_Intersample_genPSF(double     *A,
                              struct isample_struct *isample,
                              int nact, int iloop, int oloop, int gal,
                              PLASMA_sequence *sequence,
                              PLASMA_request  *request) {
  
  Quark *quark;
  Quark_Task_Flags task_flags = Quark_Task_Flags_Initializer;
  //PLASMA_desc A = *descA;
  
  PLASMA_Get_Quark(&quark);
  QUARK_Task_Flag_Set(&task_flags, TASK_SEQUENCE, (intptr_t)sequence->quark_sequence);
  
  //DAG_CORE_GenCovMat;
  QUARK_Insert_Task(quark, CORE_dIntersample_genPSF, &task_flags,
                    sizeof(struct isample_struct*),  &isample, VALUE,
                    sizeof(double)*nact*nact,        A,        INPUT,
                    sizeof(int),                     &iloop,   VALUE,
                    sizeof(int),                     &oloop,   VALUE,
                    sizeof(int),                     &gal,     VALUE,
                    0);
  return PLASMA_SUCCESS;
}
#endif

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
    //lapack_const(uplo),
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
/*
*/
void ao_pdlapack_to_tile_quark(double *Af77, int lda, PLASMA_desc *descA,
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
            //bdl = ABDL(m, n);
            bdl = (double *)plasma_getaddr(A, m, n);
            AO_CORE_dlacpy(
                quark, &task_flags,
                PlasmaUpperLower, (Y2-Y1), (X2-X1), A.mb,
                &(f77[X1*lda+Y1]), lda,
                &(bdl[X1*lda+Y1]), ldt);
        }
    }
}

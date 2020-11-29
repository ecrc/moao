#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <mpi.h>

#include "moao_defs.h"
#include "moao_scalapack.h"
#include "matcov_tiled.h"
#include "flops.h"

#ifdef USE_INTERSAMPLE
#include "intersample.h"
#endif

#include "check_utils.h"
#include "fits.h"


//helper functions && globals ========================================
static int
startswith(const char *s, const char *prefix) {
  size_t n = strlen( prefix );
  if (strncmp( s, prefix, n ))
    return 0;
  return 1;
}

#define START_TIMING(_t)               \
_t =- MPI_Wtime();

#define STOP_TIMING(_t)                \
_t += MPI_Wtime();                      \

#define USAGE if (myrank_mpi == 0) fprintf(stderr, "usage: --nprow=p --npcol=q --nb=nb_size --maxrefine=maxrefine --maxobs=maxobs --snap_per_night=spnight --output [--v] [--h]\n"\
"--nprow: number of grid row MPI Processes\n" \
"--npcol: number of grid col MPI Processes\n" \
"--nb: block size\n" \
"--maxrefine: max number of refinements\n" \
"--maxobs: max number of observations\n" \
"--filePath: path of the input files\n"\
"--suffix: filename suffix\n"\
"--snap_per_night: number of snapshots per night\n"\
"--output: dump PSFs to disk\n"\
"--v: verbose\n"\
"--h: print this help\n"\
);

#define CHECK_PSF(uplo, trans, M, N, label_time, label_matrix, matrix1, matrix2) \
{ \
        double maxerr=0, errv1=0, errv2=0, norm=0; \
        concatenateFileName(dataFile, files_path, #matrix1, suffix); \
        readFits(dataFile, matrix2); \
        norm=compareMatrices(matrix1, matrix2, M, N, uplo, trans, &maxerr, &errv1, &errv2); \
        fprintf(logFile,"%s\t      %e   %e     %e   %e    %e\n", #label_matrix, label_time, maxerr, errv1, errv2, norm); \
}

#define CHECK(uplo, trans, M, N, label_time, label_matrix, matrix1, matrix2) \
{ \
        real_t maxerr=0, errv1=0, errv2=0; \
        double norm=0; \
        real_t *Stmp = NULL; \
        if(myrank_mpi == 0){ \
            descinit_( descStmp, &M, &N, &nb, &nb, &i0, &i0, &ctxt_0, &M, &info ); \
            Stmp = (real_t*)malloc(M*N*sizeof(real_t)); \
        } \
        pXgemr2d_(&M, &N, matrix1, &i1, &i1, desc##label_matrix, Stmp, &i1, &i1, descStmp, &ictxt); \
        if (myrank_mpi == 0) { \
        concatenateFileName(dataFile, files_path, #label_matrix, suffix); \
        readFits(dataFile,matrix2); \
        norm=compareMatrices2(Stmp, matrix2, M, N, uplo, trans, &maxerr, &errv1, &errv2); \
        fprintf(logFile,"%s\t      %e   %e     %e   %e    %e\n", #label_matrix, label_time, maxerr, errv1, errv2, norm); \
        } \
        free( Stmp ); \
}





//main ===========================================================
int main( int argc, char** argv)
{
  int nprow = -1, npcol = -1;
  int myrank_mpi, nprocs_mpi;
  int maxobs = 1, maxrefine = 1;
  struct tomo_struct tomo;

  int nact = 0, nmeas = 0, nmeasts = 0, nb = 0, verbose = 0, printhelp = 0, unknown = 0, output = 0;
  int nssp = 10;
  double flops=0.0; 

  char files_path[100];
  char suffix[40]="";
  char dataFile[150]="";
#ifdef USE_INTERSAMPLE
  struct isample_struct isample;
#endif
#ifdef CHECK_SCALAPACK
  strcpy(files_path,"../datafile/check/");
#else
  strcpy(files_path,"../datafile/");
#endif
  int night_idx=0;
  int snapshots_per_night=8;
  int snapshot_idx=0;
  int obs_idx=-1;
  real_t alphaX=0.0;
  real_t alphaY=0.0;
  double *psf;
  double t_Cmm,t_Ctt,t_Ctm,t_R,t_CeeCvv,t_psf;
  double t_Cmm_time,t_Ctt_time,t_Ctm_time,t_R_time,t_CeeCvv_time,t_psf_time;
  char logName[50];
  FILE *logFile;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank_mpi);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs_mpi);

  if(argc < 2)
  {
    USAGE;
    return 0;
  }
 
  //parse options
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
    if(startswith( argv[argi], "--maxobs=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &maxobs );
      if(maxobs < 0){
        if (myrank_mpi == 0) fprintf(stderr, "Invalid number of long exposures\n"); MPI_Finalize(); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--snap_per_night=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &snapshots_per_night );
    }
    else
    if(startswith( argv[argi], "--maxrefine=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &maxrefine );
      if(maxrefine <= 0){
        if (myrank_mpi == 0) fprintf(stderr, "Invalid number of short exposures\n"); MPI_Finalize(); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--nb=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &nb );
      if(nb <= 0){
        if (myrank_mpi == 0) fprintf(stderr, "Invalid block size\n"); MPI_Finalize(); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--filePath=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%s", files_path );
    }
    else
    if(startswith( argv[argi], "--suffix=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%s", suffix );
    }
    else
    if(startswith( argv[argi], "--output"))
      output = 1;
    else
    if(startswith( argv[argi], "--v"))
      verbose = 1;
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

#ifdef USE_INTERSAMPLE
    long naxes[2];
    concatenateFileName(dataFile,files_path,"Dx",suffix);
    getFitsDims(dataFile,naxes);
    nact=naxes[0];
#endif

  matcov_init_tomo_tiled(&tomo, nssp, files_path, night_idx, snapshots_per_night, snapshot_idx, obs_idx, alphaX, alphaY);
  
  nmeas = matcov_getNumMeasurements(&tomo);
  nmeasts = matcov_getNumMeasurementsTS(&tomo);
  nmeas = nmeas - nmeasts;
  //nact = nmeasts;
  
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

  //defaults
  if(nprow < 0 || npcol < 0)
  {
    nprow = 1;
    npcol = 1;
  }

  if (verbose & (myrank_mpi == 0)) fprintf(stderr, "ScaLAPACK: MOAO starts... \n");
  if (verbose & (myrank_mpi == 0)) fprintf(stderr, "Checking arguments done\n");
  int ictxt, myrow, mycol, ioffd;
  int i0 = 0, i1 = 1;
  int iseed, info;
  int mloc_nmeas, nloc_nmeas;
  int mloc_nmeasts, nloc_nmeasts;
  int mloc_nact, nloc_nact;
  double total_time = 0.0; 
  int nbobs, nbrefine, nbgal;
  // time break down
  double alloc_cpu_time = 0.0, alloc_cpu =0.0;
  double preprocess_cpu_time = 0.0, preprocess_cpu = 0.0;
  double compute_cpu_time = 0.0, compute_cpu =0.0;
  double total_refine_time = 0.0, refine_cpu_time = 0.0, refine_cpu = 0.0;
  double total_obs_time = 0.0, obs_cpu_time = 0.0, obs_cpu = 0.0;
  double total_matcov_time = 0.0, matcov_time = 0.0;
  
  real_t *Cmm = NULL, *Ctm = NULL, *Ctt = NULL;
  real_t *R = NULL, *Cee = NULL, *Cvv = NULL, *Dx = NULL, *Tmp = NULL;
  double *Ytmp = NULL;
  real_t *Dxglo = NULL, *Cvvglo = NULL;
  int descCmm[9], descCtm[9], descCtt[9], descTmp[9];
  int descR[9], descCee[9], descCvv[9], descDx[9];
  int descDxglo[9]; descDxglo[1] = -1;
  int descCvvglo[9]; descCvvglo[1] = -1;
  int descStmp[9]; descStmp[1] = -1;

#ifdef CHECK_SCALAPACK
    if (myrank_mpi == 0) {
    strcpy(logName,"log");
    strcat(logName,suffix);
    strcat(logName,".txt");
    logFile=fopen(logName,"w");
    if(logFile==NULL){
      fprintf(stderr,"could not open log file\n");
    }
    else {
    fprintf(logFile,"maxrefine: %d\nmaxobs: %d\nnmeas: %d\nnmeasts: %d\nnact: %d\nnLayers: %ld\ngal: %ld\npsf: %ld\n", 
                            maxrefine, maxobs, nmeas, nmeasts, nact, tomo.Nlayer, tomo.nTarget,
#ifdef USE_INTERSAMPLE
    maxrefine*maxobs*tomo.nTarget
#else
    0
#endif
    );
    fprintf(logFile,"Matrices\t times(P)       max_err(L-Y)     err_val(L)     err_val(Y)      norm(L-Y)\n");
    fflush(logFile);
    }
    }
#endif

  matcov_set_gal_coords(&tomo, 0, 0);    

#ifdef USE_INTERSAMPLE
    intersample_prepare(&isample, nact*nact, tomo.Nssp[tomo.Nw-1]+1,tomo.DiamTel, files_path);
    intersample_init(&isample);
#endif

  Cblacs_get( -1, 0, &ictxt );
  Cblacs_gridinit( &ictxt, "C", nprow, npcol );
  Cblacs_gridinfo( ictxt, &nprow, &npcol, &myrow, &mycol );
  if (verbose & (myrank_mpi == 0)) fprintf(stderr, "BLACS Init done\n");

  if (myrank_mpi == 0) {
     fprintf(stderr, "# \n");
     fprintf(stderr, "# maxgal %ld nmeas %d nmeasts %d nact %d nb %d\n", tomo.nTarget, nmeas, nmeasts, nact, nb);
     fprintf(stderr, "# nprocs %d P %d Q %d\n", nprocs_mpi, nprow, npcol);
     //fprintf(stderr, "# niter %d\n", niter);
     //fprintf(stderr, "# n_range %d:%d:%d mode: %d cond: %2.4e \n", start, stop, step, mode, cond);
     fprintf(stderr, "# \n");
  }

  mloc_nmeas    = numroc_( &nmeas,   &nb, &myrow, &i0, &nprow );
  nloc_nmeas    = numroc_( &nmeas,   &nb, &mycol, &i0, &npcol );
  mloc_nmeasts  = numroc_( &nmeasts, &nb, &myrow, &i0, &nprow );
  nloc_nmeasts  = numroc_( &nmeasts, &nb, &mycol, &i0, &npcol );
  mloc_nact     = numroc_( &nact,    &nb, &myrow, &i0, &nprow );
  nloc_nact     = numroc_( &nact,    &nb, &mycol, &i0, &npcol );

  if (verbose & (myrank_mpi == 0)) fprintf(stderr, "Desc Init starts\n");
  if (verbose & (myrank_mpi == 0)) fprintf(stderr, "mloc_nmeas %d nloc_nmeas %d mloc_nmeasts %d nloc_nmeasts %d mloc_nact %d nloc_nact %d\n", 
                                                   mloc_nmeas, nloc_nmeas, mloc_nmeasts, nloc_nmeasts, mloc_nact, nloc_nact);
  descinit_( descCmm, &nmeas,  &nmeas,    &nb, &nb, &i0, &i0, &ictxt, &mloc_nmeas,   &info );
  descinit_( descCtm, &nmeasts, &nmeas,   &nb, &nb, &i0, &i0, &ictxt, &mloc_nmeasts, &info );
  descinit_( descCtt, &nmeasts, &nmeasts, &nb, &nb, &i0, &i0, &ictxt, &mloc_nmeasts, &info );
  descinit_( descR,   &nmeasts, &nmeas,   &nb, &nb, &i0, &i0, &ictxt, &mloc_nmeasts, &info );
  descinit_( descCee, &nmeasts, &nmeasts, &nb, &nb, &i0, &i0, &ictxt, &mloc_nmeasts, &info );
  descinit_( descCvv, &nact,    &nact,    &nb, &nb, &i0, &i0, &ictxt, &mloc_nact,    &info );
  descinit_( descDx,  &nact,    &nmeasts, &nb, &nb, &i0, &i0, &ictxt, &mloc_nact,    &info );
  descinit_( descTmp, &nact,    &nmeasts, &nb, &nb, &i0, &i0, &ictxt, &mloc_nact,    &info );

  int ctxt_0 = ictxt;
  Cblacs_gridinit( &ctxt_0, "C", 1, 1 );

  if (myrank_mpi == 0) {
     // Set up a process grid ctxt_0 of size 1 
     descinit_( descDxglo, &nact, &nmeasts, &nb, &nb, &i0, &i0, &ctxt_0, &nact, &info );
     Dxglo = (real_t*)calloc( (size_t)nact * nmeasts , sizeof(real_t) );
     if ( !Dxglo ) {
       fprintf(stderr, "Out of Memory for Dxglo Matrix\n");
       MPI_Finalize(); return -1;
     }
     descinit_( descCvvglo, &nact, &nact, &nb, &nb, &i0, &i0, &ctxt_0, &nact, &info );
     Cvvglo = (real_t*)malloc(nact*nact*sizeof(real_t));
     if ( !Cvvglo ) {
        fprintf(stderr, "Out of Memory for Cvvglo Matrix\n");
        MPI_Finalize(); return -1;
     }
  }

  if (verbose & (myrank_mpi == 0)) fprintf(stderr, "Desc Init ends\n");

  // Allocate matrix memory
  {
    if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Allocating matrix memory...");
    START_TIMING(alloc_cpu_time);
    Cmm = (real_t*)calloc( (size_t)mloc_nmeas * nloc_nmeas, sizeof(real_t) );
    if ( !Cmm ) {
      if (myrank_mpi == 0) fprintf(stderr, "Out of Memory for Covariance Matrix (Cmm)\n");
      MPI_Finalize(); return -1;
    }
    Ctm = (real_t*)calloc( (size_t)mloc_nmeasts * nloc_nmeas, sizeof(real_t) );
    if ( !Ctm ) {
      if (myrank_mpi == 0) fprintf(stderr, "Out of Memory for Ctm Matrix\n");
      MPI_Finalize(); return -1;
    }
    Ctt = (real_t*)calloc( (size_t)mloc_nmeasts * nloc_nmeasts, sizeof(real_t) );
    if ( !Ctt ) {
      if (myrank_mpi == 0) fprintf(stderr, "Out of Memory for Ctt Matrix\n");
      MPI_Finalize(); return -1;
    }
    R = (real_t*)calloc( (size_t)mloc_nmeasts * nloc_nmeas * tomo.nTarget, sizeof(real_t) );
    if ( !R ) {
      if (myrank_mpi == 0) fprintf(stderr, "Out of Memory for ToR (R)\n");
      MPI_Finalize(); return -1;
    }
    Cee = (real_t*)calloc( (size_t)mloc_nmeasts * nloc_nmeasts, sizeof(real_t) );
    if ( !Cee ) {
      if (myrank_mpi == 0) fprintf(stderr, "Out of Memory for Cee Matrix\n");
      MPI_Finalize(); return -1;
    }
    Cvv = (real_t*)calloc( (size_t)mloc_nact * nloc_nact, sizeof(real_t) );
    if ( !Cvv ) {
      if (myrank_mpi == 0) fprintf(stderr, "Out of Memory for Cvv Matrix\n");
      MPI_Finalize(); return -1;
    }
    Dx = (real_t*)calloc( (size_t)mloc_nact * nloc_nmeasts, sizeof(real_t) );
    if ( !Dx ) {
      if (myrank_mpi == 0) fprintf(stderr, "Out of Memory for Dx Matrix\n");
      MPI_Finalize(); return -1;
    }
    Tmp = (real_t*)calloc( (size_t)mloc_nact * nloc_nmeasts, sizeof(real_t) );
    if ( !Tmp ) {
      if (myrank_mpi == 0) fprintf(stderr, "Out of Memory for Tmp Matrix\n");
      MPI_Finalize(); return -1;
    }
    STOP_TIMING(alloc_cpu_time);
    MPI_Allreduce( &alloc_cpu_time, &alloc_cpu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    total_time += alloc_cpu;
    if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Done in %f(s)\n\n", alloc_cpu);
  }

if (myrank_mpi == 0) {
#ifdef CHECK_SCALAPACK
    Ytmp = (double*)malloc(nmeas*nmeas*sizeof(double));
#else
    Ytmp = (double*)malloc(nact*nmeasts*sizeof(double));
#endif
}

  START_TIMING(preprocess_cpu_time);
#ifdef USE_INTERSAMPLE
  //read Dx from file
  if (myrank_mpi == 0) {
     concatenateFileName(dataFile, files_path, "Dx", suffix);
     readFits(dataFile, Ytmp);
     copy(Ytmp, Dxglo, nact, nmeasts, 'A', 'T');
  }
  pXgemr2d_(&nact, &nmeasts, Dxglo, &i1, &i1, descDxglo, Dx, &i1, &i1, descDx, &ictxt);
#endif
  if (myrank_mpi == 0)
     free(Dxglo);
  STOP_TIMING(preprocess_cpu_time);
  MPI_Allreduce( &preprocess_cpu_time, &preprocess_cpu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if(verbose & (myrank_mpi == 0)) fprintf(stderr, "ScaLAPACK: MOAO Loading the interaction matrix Dx\n");
  if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Done in %f(s)\n\n", preprocess_cpu);

  total_time += preprocess_cpu;

  if(verbose & (myrank_mpi == 0)) fprintf(stderr, "ScaLAPACK: MOAO main computation started with\n");
  if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Maximum number of refinements %d\n", maxrefine);
  if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Maximum number of observations %d\n", maxobs);


  START_TIMING(compute_cpu_time);
    
  for (nbrefine = 0; nbrefine<maxrefine; nbrefine++) {
       if (verbose & (myrank_mpi == 0)) {
          fprintf(stdout,"\n\n=====>begin refine %d\n", nbrefine+1);
          fflush(stdout);
       }
       //need to get new atm-params
       night_idx = 0;
       snapshot_idx = nbrefine;
       obs_idx = -1;
       int result = !matcov_update_atm_params(&tomo, night_idx, snapshots_per_night, snapshot_idx, obs_idx);
       if(result){
          return 1;
       }
       matcov_update_tomo_tiled(&tomo);

       START_TIMING(refine_cpu_time);

       if(verbose & (myrank_mpi == 0)) fprintf(stderr, "\n Generating Cmm... ");
       START_TIMING(t_Cmm_time);
       matcov_comp_tile_scalapack(1, mloc_nmeas, nloc_nmeas, Cmm, mloc_nmeas, 0, 0, myrow, mycol, nb, nprow, npcol, &tomo);
       flops = flops + FLOPS_XCovMat(nmeas, nmeas, tomo.Nlayer, 1);
       STOP_TIMING(t_Cmm_time);
       MPI_Allreduce( &t_Cmm_time, &t_Cmm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
       total_time += t_Cmm;
       if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Done in %f(s)", t_Cmm);
       total_matcov_time += t_Cmm;

#ifdef CHECK_SCALAPACK
       //check Cmm
       CHECK('L', 'T', nmeas, nmeas, t_Cmm, Cmm, Cmm, Ytmp);
#endif

       for(nbgal=0; nbgal<tomo.nTarget; nbgal++){
	   matcov_set_gal_coords(&tomo, tomo.targetX[nbgal], tomo.targetY[nbgal]);

           if(verbose & (myrank_mpi == 0)) fprintf(stderr, "\n Galaxy %d Generating Ctm... ", nbgal+1);
           START_TIMING(t_Ctm_time);
           matcov_comp_tile_scalapack(3, mloc_nmeasts, nloc_nmeas, R+nbgal*mloc_nmeasts*nloc_nmeas, mloc_nmeasts, 0, 0, myrow, mycol, nb, nprow, npcol, &tomo);
           flops = flops + FLOPS_XCovMat(nmeas, nmeasts, tomo.Nlayer, 3);
           STOP_TIMING(t_Ctm_time);
           MPI_Allreduce( &t_Ctm_time, &t_Ctm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
           if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Done in %f(s)", t_Ctm);
           total_matcov_time += t_Ctm;
#ifdef CHECK_SCALAPACK
           //check Ctm
           CHECK('A', 'T', nmeasts, nmeas, t_Ctm, Ctm, R+nbgal*mloc_nmeasts*nloc_nmeas, Ytmp);
#endif
       }
    
       if(verbose & (myrank_mpi == 0)) fprintf(stderr, "\n Computing R... ");
       START_TIMING(t_R_time);
       info = reconstructor(nmeas, nmeasts, Cmm, descCmm, R, descR, tomo.nTarget, mloc_nmeasts*nloc_nmeas);
       STOP_TIMING(t_R_time);
       if ( info != 0 ) break;
       MPI_Allreduce( &t_R_time, &t_R, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
       if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Done in %f(s)\n\n", t_R);
       flops = flops + FLOPS_XPOTRF(nmeas);
       flops = flops + FLOPS_XTRSM('R', nmeasts, nmeas);
       flops = flops + FLOPS_XTRSM('R', nmeasts, nmeas);

#ifdef CHECK_SCALAPACK
       for (nbgal=0; nbgal<tomo.nTarget; nbgal++){
           //check R: the reconstructor
           CHECK('A','T', nmeasts, nmeas, t_R, R, R+nbgal*mloc_nmeasts*nloc_nmeas, Ytmp);
       }
#endif

       STOP_TIMING(refine_cpu_time);
       MPI_Allreduce( &refine_cpu_time, &refine_cpu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
       total_refine_time += refine_cpu;
     
       for (nbobs = 0; nbobs < maxobs; nbobs++) {
            if (verbose & (myrank_mpi == 0)) {
               fprintf(stdout,"==========>begin obs %d\n", nbobs+1);
               fflush(stdout);
            }
            //need to get new atm-params
            night_idx = 0;
            snapshot_idx = nbrefine;
            obs_idx = nbobs;
            result = !matcov_update_atm_params(&tomo, night_idx, snapshots_per_night, snapshot_idx, obs_idx);
	    if(result){
                return 1;
	    }
	    matcov_update_tomo_tiled(&tomo);

            START_TIMING(obs_cpu_time);

            if(verbose & (myrank_mpi == 0)) fprintf(stderr, "\n         Generating Cmm... ");
            START_TIMING(t_Cmm_time);
            matcov_comp_tile_scalapack(1, mloc_nmeas, nloc_nmeas, Cmm, mloc_nmeas, 0, 0, myrow, mycol, nb, nprow, npcol, &tomo);
            flops = flops + FLOPS_XCovMat(nmeas, nmeas, tomo.Nlayer, 1);
            STOP_TIMING(t_Cmm_time);
            MPI_Allreduce( &t_Cmm_time, &t_Cmm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Done in %f(s)", t_Cmm);
            total_matcov_time += t_Cmm;
#ifdef CHECK_SCALAPACK
            //check Cmm (only time: error already done)
            if (myrank_mpi == 0) {
               fprintf(logFile,"Matrices\t times(P)       max_err(L-Y)     err_val(L)     err_val(Y)      norm(L-Y)\n");
               fprintf(logFile,"Cmm\t      %e   x                x              x \n",t_Cmm);
            }
#endif

            for(nbgal=0; nbgal<tomo.nTarget; nbgal++){
                matcov_set_gal_coords(&tomo, tomo.targetX[nbgal], tomo.targetY[nbgal]);
                if(verbose & (myrank_mpi == 0)) fprintf(stderr, "\n         Galaxy %d Generating Ctm... ", nbgal+1);
                START_TIMING(t_Ctm_time);
                matcov_comp_tile_scalapack(3, mloc_nmeasts, nloc_nmeas, Ctm, mloc_nmeasts, 0, 0, myrow, mycol, nb, nprow, npcol, &tomo);
                flops = flops + FLOPS_XCovMat(nmeasts, nmeas, tomo.Nlayer, 3);
                STOP_TIMING(t_Ctm_time);
                MPI_Allreduce( &t_Ctm_time, &t_Ctm, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Done in %f(s)", t_Ctm);
                total_matcov_time += t_Ctm;
#ifdef CHECK_SCALAPACK
                //check Ctm (only time: error already done)
                if (myrank_mpi == 0) fprintf(logFile,"Ctm\t      %e   x                x              x \n",t_Ctm);
#endif

                if(verbose & (myrank_mpi == 0)) fprintf(stderr, "\n         Galaxy %d Generating Ctt... ", nbgal+1);
                START_TIMING(t_Ctt_time);
                matcov_comp_tile_scalapack(4, mloc_nmeasts, nloc_nmeasts, Ctt, mloc_nmeasts, 0, 0, myrow, mycol, nb, nprow, npcol, &tomo);
                flops = flops + FLOPS_XCovMat(nmeasts, nmeasts, tomo.Nlayer, 4);
                STOP_TIMING(t_Ctt_time);
                MPI_Allreduce( &t_Ctt_time, &t_Ctt, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Done in %f(s)", t_Ctt);
                total_matcov_time += t_Ctt;
#ifdef CHECK_SCALAPACK
                CHECK('L','T', nmeasts, nmeasts, t_Ctt, Ctt, Ctt, Ytmp);
#endif
/*
*/
                if(verbose & (myrank_mpi == 0)) fprintf(stderr, "\n         Galaxy %d Computing Cee Cvv... ", nbgal+1);
                START_TIMING(t_CeeCvv_time);
                info = compute_Cee_Cvv( nmeas, nmeasts, nact, 
                                        Cmm, descCmm, 
                                        Ctt, descCtt, 
                                        Ctm, descCtm, 
                                        R+nbgal*mloc_nmeasts*nloc_nmeas, descR, 
                                        Dx, descDx, 
                                        Cee, descCee, 
                                        Cvv, descCvv, 
                                        Tmp, descTmp,
                                        &flops );
                STOP_TIMING(t_CeeCvv_time);
                if ( info != 0 ) break;
                MPI_Allreduce( &t_CeeCvv_time, &t_CeeCvv, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Done in %f(s)", t_CeeCvv);

#ifdef CHECK_SCALAPACK
                //check Cee
                CHECK('A','T', nmeasts, nmeasts, t_CeeCvv, Cee, Cee, Ytmp);
            
                //check Cvv
                CHECK('A','T', nact, nact, t_CeeCvv, Cvv, Cvv, Ytmp);
#endif
        
#ifdef USE_INTERSAMPLE
                if(verbose & (myrank_mpi == 0)) fprintf(stderr, "\n         Galaxy %d Computing psf... ", nbgal+1);
                START_TIMING(t_psf_time);
                pXgemr2d_(&nact, &nact, Cvv, &i1, &i1, descCvv, Cvvglo, &i1, &i1, descCvvglo, &ictxt);

                //all to intersample
                if (myrank_mpi == 0) {
	           //output file for psf
                   char isample_output_filename[256];
                   int oloop = nbrefine, iloop = nbobs, gal = nbgal;
	           //compute psf
                   intersample_process(&isample, Cvvglo);
                   if (output) {
                      sprintf(isample_output_filename, "psf_oloop%d_iloop%d_gal%d.fits", oloop, iloop, gal);
                      intersample_save(&isample, isample_output_filename);
                   }
                }
                flops = flops + FLOPS_Intersample(nact, nact);
#endif
                STOP_TIMING(t_psf_time);
                MPI_Allreduce( &t_psf_time, &t_psf, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Done in %f(s)\n\n", t_psf);

#if defined(CHECK_SCALAPACK) && defined(USE_INTERSAMPLE)
                if (myrank_mpi == 0) {
                   //check psf
                   psf=isample.dphi;
                   CHECK_PSF('A','N', isample.N, isample.N, t_psf, psf, psf, Ytmp);
                }
#endif
            }// End of Galaxy loop 

            STOP_TIMING(obs_cpu_time);
            MPI_Allreduce( &obs_cpu_time, &obs_cpu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            total_obs_time += obs_cpu;
            if (verbose & (myrank_mpi == 0)) {
               fprintf(stdout,"==========>end obs %d\n\n", nbobs+1);
               fflush(stdout);
            }
#ifdef CHECK_SCALAPACK
            if (myrank_mpi == 0) {
               fprintf(logFile,"-----end obs\n\n");
               fflush(logFile);
            }
#endif
       }// End of OBS loop
#ifdef CHECK_SCALAPACK
       if (myrank_mpi == 0) {
          fprintf(logFile,"=====end refine\n\n");
          fflush(logFile);
       }
#endif
       if (verbose & (myrank_mpi == 0)) {
          fprintf(stdout,"=====>end refine %d\n\n", nbrefine+1);
          fflush(stdout);
       }
  }// End of REFINE loop

  STOP_TIMING(compute_cpu_time);
  MPI_Allreduce( &compute_cpu_time, &compute_cpu, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    
  total_time += compute_cpu;
  
  MPI_Barrier(MPI_COMM_WORLD); 

  if(info != 0) {
    if(myrank_mpi == 0) fprintf(stderr, "An error occured in MOAO! sorry...\n\n");
    return info;
  }
  else { 
     if(verbose & (myrank_mpi == 0)) fprintf(stderr, "\n\n\nFull pipeline done in %f(s)\n", compute_cpu);
     if(verbose & (myrank_mpi == 0)) fprintf(stderr, "Computation successful :-)\n\n\n\n");

     if(myrank_mpi == 0) fprintf(stderr, "  Procs   P    Q    Block  Refine  Obs  Total   Total  Total MOAO alloc MOAO preprocess Matcov  MOAO refine MOAO obs     MOAO compute |   Total    |  Gflops/s\n");
     if(myrank_mpi == 0) fprintf(stderr, "    #                size  #iter  #iter #Meas #Meas-ts #Actu  time (s)   time (s)       time (s)    time (s)   time(s)   time (s)     |   time(s)  |  Perf\n");
     if(myrank_mpi == 0) fprintf(stderr, " ======= ==== ===== ====== ====== ===== ===== ======== ===== ========== =============== ========= =========== ======== ============== | ========== | ========\n");
     if(myrank_mpi == 0) fprintf(stderr, "    %d    %d    %d    %d     %d     %d   %-8d%-8d%-4d  %-11.3f    %-10.3f   %-10.3f%-10.3f%-10.3f %-4.3f            %-10.3f %-4.3f\n\n", nprocs_mpi, nprow, npcol, nb, maxrefine, maxobs, nmeas, nmeasts, nact, alloc_cpu, preprocess_cpu, total_matcov_time, total_refine_time, total_obs_time, compute_cpu, total_time, flops / 1e9 / total_time);
  }

  if(verbose & (myrank_mpi == 0)) fprintf(stderr, "ScaLAPACK: MOAO done...\n\n" );

  matcov_free_tomo_tiled(&tomo);

  free( Cmm );
  free( Ctm );
  free( Ctt );
  free( R );
  free( Cee );
  free( Cvv );
  free( Dx );
  free( Tmp );
  if (myrank_mpi == 0) {
     free(Ytmp);
     free(Cvvglo);
  }
  
  MPI_Barrier(MPI_COMM_WORLD); 
  MPI_Finalize(); 

  return 0;
}

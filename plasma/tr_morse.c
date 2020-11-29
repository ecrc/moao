//
//  morse_test.c
//
//
//  Created by Ali M Charara on ....
//
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#include <cublas.h>
#include <magma.h>

#include <magma_lapack.h>

#include <coreblas.h>

#include "morse.h"
#include "auxiliary.h"
#include "descriptor.h"
#include "tile.h"

#include "compute_d.h"
#include "context.h"
//#include "codelet_d.h"
//#include "morse_kernels.h"

#ifdef USE_TRACE
#include <starpu_fxt.h>
//#include <eztrace.h>
#endif


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
void magma_dmake_hpd( magma_int_t N, double* A, magma_int_t lda );
//void cl_dgemm_restrict_where(uint32_t where);

short verbose = 0, warmUp = 0, synchronize = 0, calibrate=1, merged=0,genprofile=0,restrict_gemm=0;

#if defined(USE_TRACE)
short tracing_on = 0;
#endif
std::string sched_policy;


int MORSE_GenTomRecTiled(int N, int TS, int cmaa_rows, short verbose, double& total_time, short sync);
int morse_pdpotrimm(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *C,
                    MORSE_sequence_t *sequence, MORSE_request_t *request);
double cWtime(void)
{
  struct timeval tp;
  gettimeofday( &tp, NULL );
  return tp.tv_sec + 1e-6 * tp.tv_usec;
}

#if defined(USE_TRACE)
#define START_TRACING()       \
starpu_fxt_start_profiling();//eztrace_start();
#define STOP_TRACING()        \
starpu_fxt_stop_profiling();//eztrace_stop();
#else /* defined(MORSE_TRACE) */
#define START_TRACING() if( 0 ) {};
#define STOP_TRACING() if( 0 ) {};
#endif /* defined(MORSE_EZTRACE) */

#define START_TIMING(_t)               \
START_TRACING();                    \
_t = -cWtime();

#define STOP_TIMING(_t)                \
_t += cWtime();                      \
STOP_TRACING();

#define USAGE printf("usage: --n_cores=cores --n_gpus=gpus --N=matrix_size --tile=tile_size --cmaa=cmaa_rows [-m][-v][-s][-p][-r]\n"\
"--n_cores: number of CPU cores to use, maximum value should be total_machine_cores - n_gpus, default: total_machine_cores - n_gpus\n" \
"--n_gpus: number of GPU's to use, default: 1\n"  \
"--N: matrix size\n" \
"--tile: tile size\n" \
"--cmaa: number of rows of Cmaa matrix (default 3500)\n" \
"-v: verbose\n"\
"-s: synchronize between stages, default off\n"\
"-m: force merged fours stages\n"\
"-r: restrict GEMM on GPU only, default off\n"\
);

//main ===========================================================
int main( int argc, char** argv)
{
  //#ifdef USE_TRACE
  //  if(tracing_on)
  //    starpu_fxt_stop_profiling();
  //#endif
  int ncores = -1, ngpus = -1;
  if(argc < 3)
  {
    USAGE;
    return 0;
  }
  int N=0, TS=0, cmaa_rows = -1;
  int workers_per_cuda=-1;
  
  //parse options
  int argi;
  for(argi = 1; argi < argc; argi++)
  {
    if(startswith( argv[argi], "--n_cores=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", &ncores );
    }
    else
      if(startswith( argv[argi], "--n_gpus=")){
        sscanf( strchr( argv[argi], '=' ) + 1, "%d", &ngpus );
        if(1 <= ngpus){
          int ndevices;
          cudaGetDeviceCount( &ndevices );
          if(MagmaMaxGPUs < ngpus || ndevices < ngpus){
            printf("asked for too many GPUs (%d), maximum is %d, detected %d\n", ngpus, min(MagmaMaxGPUs,ndevices), ndevices); return 0;
          }
        }
      }
      else
        if(startswith( argv[argi], "--N=")){
          sscanf( strchr( argv[argi], '=' ) + 1, "%d", &N );
          if(N <= 0){
            printf("matrix size is too small! please use valid size\n"); return 0;
          }
        }
        else
          if(startswith( argv[argi], "--tile=")){
            sscanf( strchr( argv[argi], '=' ) + 1, "%d", &TS );
            if(TS <= 0){
              printf("tile size is too small! please use valid size\n"); return 0;
            }
          }
          else
            if(startswith( argv[argi], "--cmaa=")){
              sscanf( strchr( argv[argi], '=' ) + 1, "%d", &cmaa_rows );
              if(cmaa_rows <= 0){
                printf("cmaa rows is too small! please use valid size\n"); return 0;
              }
            }
            else
              if(startswith( argv[argi], "-v"))
                verbose = 1;
              else
                if(startswith( argv[argi], "-w"))
                  warmUp = 1;
                else
                  if(startswith( argv[argi], "-s"))
                    synchronize = 1;
                  else
                    if(startswith( argv[argi], "-m"))
                      merged = 1;
                    else
                      if(startswith( argv[argi], "-p"))
                        genprofile = 1;
                    else
                      if(startswith( argv[argi], "-r"))
                        restrict_gemm = 1;
                      else
                      {
                        printf("unkown option %s, aborting...\n", argv[argi]); return 0;
                      }
  }
  //defaults
  if(ngpus < 0)
    ngpus = 1;
  if(ncores < 0)
  {
    get_thread_count(&ncores);
    ncores -= ngpus;
  }
  if(cmaa_rows < 0)
    cmaa_rows = 3500;
  if(N <= 0)
  {
    printf("please provide Matrix size\n");
    USAGE;
    return 0;
  }
  if(TS <= 0)
  {
    printf("please provide Tile size\n");
    USAGE;
    return 0;
  }
  if (!getenv("STARPU_SCHED"))
    sched_policy = "dmdas";
  else
  {
    sched_policy = getenv("STARPU_SCHED");
    //printf("env: %s, sp: %s\n", getenv("STARPU_SCHED"), sched_policy.c_str());
  }
  if (!getenv("STARPU_CALIBRATE"))
    calibrate = 1;
  else
    sscanf( getenv("STARPU_CALIBRATE"), "%d", &calibrate );
  //calibrate = getenv("STARPU_CALIBRATE");
  
  if(verbose) printf( "processing matrix of size: %d with tile size: %d\n", N, TS);
  if(verbose) printf( "input cmaa matrix rows: %d \n", cmaa_rows);
  if(verbose) printf( "working on: %d GPUs, %d CPU cores\n", ngpus, ncores);
  if(verbose) printf( "using scheduling policy: %s \n", sched_policy.c_str());
  /*
   if(getenv("MORSE_NWORKER_PER_CUDA"))
   {
   //sscanf( getenv("MORSE_NWORKER_PER_CUDA"), "%d", &workers_per_cuda );
   if(verbose) printf( "using %d workers per cuda\n", workers_per_cuda);
   MORSE_InitPar(ncores, ngpus, workers_per_cuda);
   }
   else*/
  MORSE_Init(ncores, ngpus);
  
  MORSE_Disable(MORSE_AUTOTUNING);
  MORSE_Set(MORSE_INNER_BLOCK_SIZE, 32 );
  MORSE_Set(MORSE_HOUSEHOLDER_MODE, MORSE_FLAT_HOUSEHOLDER);
  MORSE_Set(MORSE_TRANSLATION_MODE, MORSE_INPLACE);
  //printf("env: %s, sp: %s\n", morse->schedopt.starpu->sched_policy_name, sched_policy.c_str());
  if(genprofile)
    MORSE_Enable(MORSE_PROFILING_MODE);
  if(verbose) printf( "MORSE initialize... done\n");
#ifdef USE_TRACE
  STOP_TRACING();
#endif
  double total_time;
  int result = MORSE_GenTomRecTiled(N, TS, cmaa_rows, verbose, total_time, synchronize);
  if(result != MORSE_SUCCESS)
    printf("an error occured! sorry... %d\n", result);
  
  //save results to file
  {
    std::string fileName;
    if(synchronize)
      fileName ="timing_tr_morse_synch.txt";
    else
      fileName ="timing_tr_morse.txt";
    FILE* outfile = fopen(fileName.c_str(), "a");
    if(outfile == NULL)
    {
      printf( "could not open output file!\n" );
    }
    else
    {
      double flops = (long int)N * ((long int)N * ((long int)N + 1.0f) + 1.0f )/(double)1e9
      + 2.0 * (long int)cmaa_rows * (long int)N * (long int)N/(double)1e9;
      fprintf(outfile, "%d %d %d %d %f %lf %f %s %d\n", N, ngpus, TS, cmaa_rows, total_time, flops, flops / total_time, sched_policy.c_str(), calibrate);
      fclose(outfile);
      printf("%-8.3f\n",flops / total_time);
    }
  }
  //if(verbose) printf( "MORSE: done...\n" );
  
  
  MORSE_Finalize();
  if(verbose) printf( "MORSE: done...\n" );
  return 0;
}

// ===========================================================
int MORSE_GenTomRecTiled(int N, int TS, int cmaa_rows, short verbose, double& total_time, short sync = 0)
{
  int result = MORSE_SUCCESS;
  // time break down
  double alloc_cpu_time = 0.0;
  double init_cpu_time = 0.0;
  double r_time = 0.0;
  total_time = 0.0;
  
  size_t Cmaa_ld = cmaa_rows, caa_ld = N, caa_size = caa_ld * N, Cmaa_size = Cmaa_ld * N;
  
  
  double *Cmaa = NULL, *caa = NULL, *res = NULL;
  MORSE_desc_t *descCmaa = NULL, *descCaa = NULL, *descRes = NULL;
  
  MORSE_sequence_t *sequence;
  MORSE_request_t request[3] = { MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER };
  
  MORSE_Set(MORSE_TILE_SIZE,        TS );
  MORSE_Sequence_Create(&sequence);
  
  //allocate matrix memory
  {
    if(verbose) printf( "allocating matrix memory...");
    alloc_cpu_time -= cWtime();
    Cmaa = (double*)malloc( (size_t)Cmaa_size * sizeof(double) );
    if ( ! Cmaa ) {
      fprintf(stderr, "Out of Memory for Cmaa Matrix\n");
      return -1;
    }
    caa = (double*)malloc( (size_t)caa_size * sizeof(double) );
    if ( ! caa ) {
      fprintf(stderr, "Out of Memory for Covariance Matrix (caa)\n");
      return -1;
    }
    res = (double*)malloc( (size_t)Cmaa_size * sizeof(double) );
    if ( ! res ) {
      fprintf(stderr, "Out of Memory for result Matrix (res)\n");
      return -1;
    }
    alloc_cpu_time += cWtime();
    total_time += alloc_cpu_time;
    if(verbose) printf( "done in %f(s)\n", alloc_cpu_time);
  }
  MORSE_Desc_Create(&descCmaa, Cmaa, MorseRealDouble, TS, TS, TS*TS, Cmaa_ld, N, 0, 0, Cmaa_ld, N,1,1);
  MORSE_Desc_Create(&descCaa, caa, MorseRealDouble, TS, TS, TS*TS, caa_ld, N, 0, 0, caa_ld, N,1,1);
  MORSE_Desc_Create(&descRes, res, MorseRealDouble, TS, TS, TS*TS, Cmaa_ld, N, 0, 0, Cmaa_ld, N,1,1);
  
  //initiating random matrix
  {
    if(verbose) printf( "generating the random matrices...");
    init_cpu_time -= cWtime();
    
    MORSE_dplrnt_Tile(descCmaa, rand()%100);
    //lapackf77_dlarnv( &ione, ISEED, (magma_int_t*)&Cmaa_size, Cmaa );
    MORSE_dplgsy_Tile( (double)N, descCaa, rand()%100 );
    
    init_cpu_time += cWtime();
    total_time += init_cpu_time;
    if(verbose) printf( "done in %f\n", init_cpu_time);
  }
  
  
  
  //do the Cholesky Factorization and Inversion
  //then compute Tomograpich Reconstruction matrix
  //all in one sequance
  //if(0)
  do
  {
    if(verbose) printf( "MORSE: Matrix Factorization and Inversion started... in %d stage(s) ... ", sync ? 3 : 1);
    if(verbose) fflush(stdout);
    
    double alpha = 1.0, beta = 0.0;
    //r_time -= cWtime();
    START_TIMING(r_time);
    
    if(merged && !sync)
    {
      //result = MORSE_dpotrimm_Tile_Async(MorseLower, descCaa, descCmaa, descRes, sequence, &request[0]);
      result = morse_pdpotrimm(MorseLower, descCaa, descCmaa, descRes, sequence, &request[0]);
      if(MORSE_SUCCESS != result) break;
    }
    else
    {
      result = MORSE_dpotrf_Tile_Async(MorseLower, descCaa, sequence, &request[0]);
      if(MORSE_SUCCESS != result) break;
      if(sync)
        MORSE_Sequence_Wait(sequence);
      result = MORSE_dpotri_Tile_Async(MorseLower, descCaa, sequence, &request[1]);
      if(MORSE_SUCCESS != result) break;
      if(sync)
        MORSE_Sequence_Wait(sequence);
      
      result = MORSE_dsymm_Tile_Async(MorseRight, MorseLower,
                                      alpha, descCaa, descCmaa,
                                      beta,  descRes,
                                      sequence, &request[2]);
    }
    MORSE_Sequence_Wait(sequence);
    
    STOP_TIMING(r_time);
    
    MORSE_Desc_getoncpu( descRes );
    
    //STOP_TIMING(r_time);
    //r_time += cWtime();
    total_time += r_time;
    if(verbose) printf( "done in %f(s)\n", r_time);
  }while(0);
  
  MORSE_Desc_getoncpu( descCmaa );
  MORSE_Desc_getoncpu( descCaa );
  
  MORSE_Sequence_Destroy(sequence);
  
  MORSE_Desc_Destroy( &descCmaa );
  MORSE_Desc_Destroy( &descCaa );
  MORSE_Desc_Destroy( &descRes );
  
  if(verbose) printf( "computation successful :-)\n");
  printf(" Matrix   cpu alloc   cpu init   R=Cmaa * inv(Caa) | Total      |  GFLOP/s\n");
  printf(" Dim       time (s)   time (s)     time (s)        | time(s)    |  \n");
  printf(" ======   =========   ========   ================= | ========== | ========\n");
  printf(" %-8d %-11.5f %-10.3f %-20.3f %-12.3f", N, alloc_cpu_time, init_cpu_time, r_time, total_time);

  total_time = r_time;
  return result;
}

#define AA(i,j)  AA[i + j*lda]
void magma_dmake_hpd( magma_int_t N, double* AA, magma_int_t lda )
{
  magma_int_t i, j;
  for( i=0; i<N; ++i ) {
    AA(i,i) = MAGMA_D_MAKE( MAGMA_D_REAL( AA(i,i) ) + N, 0. );
    for( j=0; j<i; ++j ) {
      AA(j,i) = MAGMA_D_CNJG( AA(i,j) );
    }
  }
}

short enable_priority = 1;

#define A(m,n) A,  m,  n
#define B(m,n) B,  m,  n
#define C(m,n) C,  m,  n
#define SetPrio(op, p) if(enable_priority) op.priority += p;
#define ResetPrio(op) if(enable_priority) op.priority = def_prio;
int morse_pdpotrimm(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *C,
                    MORSE_sequence_t *sequence, MORSE_request_t *request)
{
  MORSE_context_t *morse;
  MORSE_option_t options;
  
  int k, m, n;
  int lda, ldab, ldb, ldc;
  int ldak, ldam, ldan;
  int tempkm, tempmm, tempnn, tempkn;
  
  double alpha  = (double) 1.0;
  double beta  = (double) 0.0;
  double zbeta;
  double zone  = (double) 1.0;
  double mzone = (double)-1.0;
  
  
  morse = morse_context_self();
  if (sequence->status != MORSE_SUCCESS)
    return MORSE_ERR_UNEXPECTED;
  RUNTIME_options_init(&options, morse, sequence, request);
  
#ifdef MAGMAMORSE_USE_MAGMA
  {
    int nb = A->nb;//magma_get_dpotrf_nb(A->n);
    RUNTIME_options_ws_alloc( &options, nb*nb, 0 );
  }
#endif
  //if(restrict_gemm) cl_dgemm_restrict_where( STARPU_CUDA );
  //RUNTIME_dlocality_onerestrict(MORSE_GEMM, STARPU_CUDA); 
  
  
  /*
   *  MorseLower
   */
  if (uplo == MorseLower) {
    /*
     *  DPOTRF
     */
    int def_prio = options.priority;
    //enum {prio_high = 10000, prio_med = 1000, prio_low = 100};
    //if(verbose) printf("default priority is: %d\n", def_prio);
    
    for (k = 0; k < A->mt; k++) {
	    tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
	    ldak = BLKLDD(A, k);
      SetPrio(options, 20000);
	    MORSE_TASK_dpotrf(
                        &options,
                        MorseLower, tempkm, A->mb,
                        A(k, k), ldak, A->nb*k);
      ResetPrio(options);
	    
      for (m = k+1; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
        ldam = BLKLDD(A, m);
        if(m == k+1) {SetPrio(options, 15000);} else {SetPrio(options, 10000);}
        MORSE_TASK_dtrsm(
                         &options,
                         MorseRight, MorseLower, MorseTrans, MorseNonUnit,
                         tempmm, A->mb, A->mb,
                         zone, A(k, k), ldak,
                         A(m, k), ldam);
        ResetPrio(options);
      }
      MORSE_TASK_dataflush( &options, A(k, k) );
      
	    for (n = k+1; n < A->nt; n++) {
        tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
        ldan = BLKLDD(A, n);
        SetPrio(options, 5000);
        MORSE_TASK_dsyrk(
                         &options,
                         MorseLower, MorseNoTrans,
                         tempnn, A->nb, A->mb,
                         -1.0, A(n, k), ldan,
                         1.0, A(n, n), ldan);
        ResetPrio(options);
        
        for (m = n+1; m < A->mt; m++) {
          tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;
          ldam = BLKLDD(A, m);
          SetPrio(options, 500);
          MORSE_TASK_dgemm(
                           &options,
                           MorseNoTrans, MorseTrans,
                           tempmm, tempnn, A->mb, A->mb,
                           mzone, A(m, k), ldam,
                           A(n, k), ldan,
                           zone,  A(m, n), ldam);
          ResetPrio(options);
        }
        MORSE_TASK_dataflush( &options, A(n, k) );
	    }
    }
    /*
     *  DTRTRI
     */
    for (n = 0; n < A->nt; n++) {
      tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
      ldan = BLKLDD(A, n);
      for (m = n+1; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
        ldam = BLKLDD(A, m);
        SetPrio(options, 15000);
        MORSE_TASK_dtrsm(
                         &options,
                         MorseRight, uplo, MorseNoTrans, MorseNonUnit,
                         tempmm, tempnn, A->mb,
                         mzone, A(n, n), ldan,
                         A(m, n), ldam);
        ResetPrio(options);
      }
      for (m = n+1; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
        ldam = BLKLDD(A, m);
        for (k = 0; k < n; k++) {
          tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
          SetPrio(options, 3000);
          MORSE_TASK_dgemm(
                           &options,
                           MorseNoTrans, MorseNoTrans,
                           tempmm, tempkn, tempnn, A->mb,
                           zone, A(m, n), ldam,
                           A(n, k), ldan,
                           zone, A(m, k), ldam);
          ResetPrio(options);
        }
      }
      for (m = 0; m < n; m++) {
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
        SetPrio(options, 2000);
        MORSE_TASK_dtrsm(
                         &options,
                         MorseLeft, uplo, MorseNoTrans, MorseNonUnit,
                         tempnn, tempmm, A->mb,
                         zone, A(n, n), ldan,
                         A(n, m), ldan);
        ResetPrio(options);
      }
      SetPrio(options, 20000);
      MORSE_TASK_dtrtri(
                        &options,
                        uplo, MorseNonUnit,
                        tempnn, A->mb,
                        A(n, n), ldan, A->nb*n);
      ResetPrio(options);
    }
    /*
     *  DLAUUM
     */
    for (m = 0; m < A->mt; m++) {
      tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
      ldam = BLKLDD(A, m);
      for(n = 0; n < m; n++) {
        tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
        SetPrio(options, 1000);
        MORSE_TASK_dsyrk(
                         &options,
                         uplo, MorseTrans,
                         tempnn, tempmm, A->mb,
                         1.0, A(m, n), ldam,
                         1.0, A(n, n), A->mb);
        ResetPrio(options);
        
        for(k = n+1; k < m; k++) {
          tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
          SetPrio(options, 200);
          MORSE_TASK_dgemm(
                           &options,
                           MorseTrans, MorseNoTrans,
                           tempkm, tempnn, tempmm, A->mb,
                           zone, A(m, k), ldam,
                           A(m, n), ldam,
                           zone, A(k, n), A->mb);
          ResetPrio(options);
        }
      }
      for (n = 0; n < m; n++) {
        tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
        SetPrio(options, 200);
        MORSE_TASK_dtrmm(
                         &options,
                         MorseLeft, uplo, MorseTrans, MorseNonUnit,
                         tempmm, tempnn, A->mb,
                         zone, A(m, m), ldam,
                         A(m, n), ldam);
        ResetPrio(options);
      }
      SetPrio(options, 20000);
      MORSE_TASK_dlauum(
                        &options,
                        uplo,
                        tempmm,
                        A->mb, A(m, m), ldam);
      ResetPrio(options);
    }
    /*
     *  DSYMM
     */
    for (m = 0; m < C->mt; m++) {
      tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
      ldc = BLKLDD(C, m);
      for (n = 0; n < C->nt; n++) {
        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
        lda = BLKLDD(A, n);
        ldb = BLKLDD(B, m);
        for (k = 0; k < C->nt; k++) {
          tempkn = k == C->nt-1 ? C->n-k*C->nb : C->nb;
          ldak = BLKLDD(A, k);
          zbeta = k == 0 ? beta : zone;
          if (k < n) {
            MORSE_TASK_dgemm(
                             &options,
                             MorseNoTrans, MorseTrans,
                             tempmm, tempnn, tempkn, A->mb,
                             alpha, B(m, k), ldb,  /* ldb * K */
                             A(n, k), lda,  /* lda * K */
                             zbeta, C(m, n), ldc); /* ldc * Y */
          }
          else {
            if (k == n) {
              MORSE_TASK_dsymm(
                               &options,
                               MorseRight, uplo,
                               tempmm, tempnn, A->mb,
                               alpha, A(k, k), ldak, /* ldak * Y */
                               B(m, k), ldb,  /* ldb  * Y */
                               zbeta, C(m, n), ldc); /* ldc  * Y */
            }
            else {
              MORSE_TASK_dgemm(
                               &options,
                               MorseNoTrans, MorseNoTrans,
                               tempmm, tempnn, tempkn, A->mb,
                               alpha, B(m, k), ldb,  /* ldb  * K */
                               A(k, n), ldak, /* ldak * Y */
                               zbeta, C(m, n), ldc); /* ldc  * Y */
            }
          }
        }
      }
    }
  }
  /*
   *  MorseUpper
   */
  else {
    /*
     *  DPOTRF
     */
    for (k = 0; k < A->nt; k++) {
      tempkm = k == A->nt-1 ? A->n-k*A->nb : A->nb;
      ldak = BLKLDD(A, k);
      MORSE_TASK_dpotrf(
                        &options,
                        MorseUpper,
                        tempkm, A->mb,
                        A(k, k), ldak, A->nb*k);
      
      for (n = k+1; n < A->nt; n++) {
        tempnn = n == A->nt-1 ? A->n - n*A->nb : A->nb;
        MORSE_TASK_dtrsm(
                         &options,
                         MorseLeft, MorseUpper, MorseTrans, MorseNonUnit,
                         A->mb, tempnn, A->mb,
                         zone, A(k, k), ldak,
                         A(k, n), ldak);
      }
      MORSE_TASK_dataflush( &options, A(k, k) );
      
      
	    for (m = k+1; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m - m*A->mb : A->mb;
        ldam = BLKLDD(A, m);
        
        MORSE_TASK_dsyrk(
                         &options,
                         MorseUpper, MorseTrans,
                         tempmm, A->mb, A->mb,
                         -1.0, A(k, m), ldak,
                         1.0, A(m, m), ldam);
        
        
        for (n = m+1; n < A->nt; n++) {
          tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
          
          MORSE_TASK_dgemm(
                           &options,
                           MorseTrans, MorseNoTrans,
                           tempmm, tempnn, A->mb, A->mb,
                           mzone, A(k, m), ldak,
                           A(k, n), ldak,
                           zone,  A(m, n), ldam);
        }
        MORSE_TASK_dataflush( &options, A(k, m) );
	    }
    }
    /*
     *  DTRTRI
     */
    for (m = 0; m < A->mt; m++) {
      tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
      ldam = BLKLDD(A, m);
      for (n = m+1; n < A->nt; n++) {
        tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
        MORSE_TASK_dtrsm(
                         &options,
                         MorseLeft, uplo, MorseNoTrans, MorseNonUnit,
                         tempmm, tempnn, A->mb,
                         mzone, A(m, m), ldam,
                         A(m, n), ldam);
      }
      for (n = 0; n < m; n++) {
        tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
        ldan = BLKLDD(A, n);
        for (k = m+1; k < A->nt; k++) {
          tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
          MORSE_TASK_dgemm(
                           &options,
                           MorseNoTrans, MorseNoTrans,
                           tempnn, tempkn, tempmm, A->mb,
                           zone, A(n, m), ldan,
                           A(m, k), ldam,
                           zone, A(n, k), ldan);
        }
        MORSE_TASK_dtrsm(
                         &options,
                         MorseRight, uplo, MorseNoTrans, MorseNonUnit,
                         tempnn, tempmm, A->mb,
                         zone, A(m, m), ldam,
                         A(n, m), ldan);
      }
      MORSE_TASK_dtrtri(
                        &options,
                        uplo, MorseNonUnit,
                        tempmm, A->mb,
                        A(m, m), ldam, A->mb*m);
    }
    /*
     *  DLAUUM
     */
    for (m = 0; m < A->mt; m++) {
      tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
      ldam = BLKLDD(A, m);
      for (n = 0; n < m; n++) {
        tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
        MORSE_TASK_dsyrk(
                         &options,
                         uplo, MorseNoTrans,
                         tempnn, tempmm, A->mb,
                         1.0, A(n, m), A->mb,
                         1.0, A(n, n), A->mb);
        
        for (k = n+1; k < m; k++){
          tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
          MORSE_TASK_dgemm(
                           &options,
                           MorseNoTrans, MorseTrans,
                           tempnn, tempkm, tempmm, A->mb,
                           zone, A(n, m), A->mb,
                           A(k, m), A->mb,
                           zone, A(n, k), A->mb);
        }
      }
      for (n = 0; n < m; n++) {
        tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
        MORSE_TASK_dtrmm(
                         &options,
                         MorseRight, uplo, MorseTrans, MorseNonUnit,
                         tempnn, tempmm, A->mb,
                         zone, A(m, m), ldam,
                         A(n, m), A->mb);
      }
      MORSE_TASK_dlauum(
                        &options,
                        uplo,
                        tempmm,
                        A->mb, A(m, m), ldam);
    }
    /*
     *  DSYMM
     */
    for (m = 0; m < C->mt; m++) {
      tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
      ldc = BLKLDD(C, m);
      for (n = 0; n < C->nt; n++) {
        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
        lda = BLKLDD(A, n);
        ldb = BLKLDD(B, m);
        for (k = 0; k < C->nt; k++) {
          tempkn = k == C->nt-1 ? C->n-k*C->nb : C->nb;
          ldak = BLKLDD(A, k);
          zbeta = k == 0 ? beta : zone;
          if (k < n) {
            MORSE_TASK_dgemm(
                             &options,
                             MorseNoTrans, MorseNoTrans,
                             tempmm, tempnn, tempkn, A->mb,
                             alpha, B(m, k), ldb,  /* ldb  * K */
                             A(k, n), ldak, /* ldak * Y */
                             zbeta, C(m, n), ldc); /* ldc  * Y */
          }
          else {
            if (k == n) {
              MORSE_TASK_dsymm(
                               &options,
                               MorseRight, uplo,
                               tempmm, tempnn, A->mb,
                               alpha, A(k, k), ldak, /* ldak * Y */
                               B(m, k), ldb,  /* ldb  * Y */
                               zbeta, C(m, n), ldc); /* ldc  * Y */
            }
            else {
              MORSE_TASK_dgemm(
                               &options,
                               MorseNoTrans, MorseTrans,
                               tempmm, tempnn, tempkn, A->mb,
                               alpha, B(m, k), ldb,  /* ldb * K */
                               A(n, k), lda,  /* lda * K */
                               zbeta, C(m, n), ldc); /* ldc * Y */
            }
          }
        }
      }
    }
  }
  
#ifdef MAGMAMORSE_USE_MAGMA
  RUNTIME_options_ws_free(&options);
#endif
  RUNTIME_options_finalize(&options, morse);
  
  return MORSE_SUCCESS;
}

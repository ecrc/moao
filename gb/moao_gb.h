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
void ao_pdlapack_to_tile_quark(double *Af77, int lda, PLASMA_desc *A,
                                   PLASMA_sequence *sequence, PLASMA_request *request);

#define START_TIMING(_t)               \
_t =- MPI_Wtime();

#define STOP_TIMING(_t)                \
_t += MPI_Wtime();                      \

#ifndef USE_MATCOV_TILED 
void PRINT_USAGE( int __nb, int __ts, int __nights_obs, int __nhrs_per_night, int __nsnapshot_per_hr, int __nmeas, int __nmeasts, int __nact, int __ngal){
  fprintf(stderr, "usage: --p=nprow --q=npcol --nb=nb_size --tile=tile_size  --nnights=nights_tosimulate --nhrs_night=num_hours_per_night --nsnap_hr=num_snaps_per_hour --nmeas=nmeas --nmeasts=nbmeasts --nact=nbact --maxrefine=maxrefine --maxobs=maxobs [--v] [--s] [--h]\n");
  fprintf(stderr, "--nprow: number of grid row MPI Processes\n" );
  fprintf(stderr, "--npcol: number of grid col MPI Processes\n" );
  fprintf(stderr, "--nb: Block size (%d)\n", __nb );
  fprintf(stderr, "--tile: tunable tile size (%d)\n",__ts );
  fprintf(stderr, "--nnights: number of nights to simulate (%d)\n",__nights_obs );
  fprintf(stderr, "--nhrs_night: number of hours per night to simulate (%d)\n",__nhrs_per_night );
  fprintf(stderr, "--nsnap_hr: number of snapshots per hours to simulate (%d)\n",__nsnapshot_per_hr );
  fprintf(stderr, "--nmeas: number of measurements (%d)\n", __nmeas);
  fprintf(stderr, "--nmeasts: number of measurements in the true sensor (%d)\n", __nmeasts);
  fprintf(stderr, "--nact: number of actuators (%d)\n", __nact);
  fprintf(stderr, "--ngal: number of galaxies (%d)\n", __ngal );
  fprintf(stderr, "--maxrefine: max number of refinements\n");
  fprintf(stderr, "--maxobs: max number of observations\n" );
  fprintf(stderr, "--ncores: number of cores per node\n" );
  fprintf(stderr, "--v: verbose\n");
  fprintf(stderr, "--s: synchronize between stages, default off\n");
  fprintf(stderr, "--h: print this help\n");
}
#define USAGE( __nb, __ts, __nights_obs, __nhrs_per_night, __nsnapshot_per_hr, __nssp, __nmeas, __nmeasts, __nact, __ngal) \
if (myrank_mpi == 0){ \
  PRINT_USAGE( __nb, __ts, __nights_obs, __nhrs_per_night, __nsnapshot_per_hr, __nmeas, nmeasts, __nact, __ngal); \
}
#else
void PRINT_USAGE( int __nb, int __ts, int __nights_obs, int __nhrs_per_night, int __nsnapshot_per_hr, int __nssp, int __nact, int __ngal) {
  fprintf(stderr, "usage: --p=nprow --q=npcol --nb=nb_size --tile=tile_size --nnights=nights_tosimulate --nhrs_night=num_hours_per_night --nsnap_hr=num_snaps_per_hour --nssp=nssp --nact=nbact --maxrefine=maxrefine --maxobs=maxobs [--v] [--s] [--h]\n");
  fprintf(stderr, "--nprow: number of grid row MPI Processes\n" );
  fprintf(stderr, "--npcol: number of grid col MPI Processes\n" );
  fprintf(stderr, "--nb: Block size (%d)\n", __nb );
  fprintf(stderr, "--tile: tunable tile size (%d)\n",__ts );
  fprintf(stderr, "--nnights: number of nights to simulate (%d)\n",__nights_obs );
  fprintf(stderr, "--nhrs_night: number of hours per night to simulate (%d)\n",__nhrs_per_night );
  fprintf(stderr, "--nsnap_hr: number of snapshots per hours to simulate (%d)\n",__nsnapshot_per_hr );
  fprintf(stderr, "--nssp: number of measurements (%d)\n", __nssp );
  fprintf(stderr, "--nact: number of actuators (%d)\n", __nact); 
  fprintf(stderr, "--ngal: number of galaxies (%d)\n", __ngal ); 
  fprintf(stderr, "--maxrefine: max number of refinements\n"); 
  fprintf(stderr, "--maxobs: max number of observations\n" );
  fprintf(stderr, "--ncores: number of cores per node\n" ); 
  fprintf(stderr, "--v: verbose\n"); 
  fprintf(stderr, "--s: synchronize between stages, default off\n"); 
  fprintf(stderr, "--h: print this help\n"); 
}
#define USAGE( __nb, __ts, __nights_obs, __nhrs_per_night, __nsnapshot_per_hr, __nssp, __nmeas, __nmeasts, __nact, __ngal) \
if (myrank_mpi == 0){ \
  PRINT_USAGE( __nb, __ts, __nights_obs, __nhrs_per_night, __nsnapshot_per_hr, __nssp, __nact, __ngal); \
  }
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

//#define DATA_FILES_PATH "../datafiles/"
#define DATA_FILES_PATH "./"


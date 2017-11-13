/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "stdio.h"
#include <string.h>
#include "arguments.h"

#include <unistd.h>

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

int read_args(int argc, char **argv,
    int *ncores,
    int *ngpus,
    int *nact,
    int *nmeas,
    int *nmeasts,
    int *maxobs,
    int *maxrefine,
    int *ts,
    char *sys_path,
    char *atm_path,
    char *suffix,
#ifdef USE_PLASMA
    int *runtime,
#endif
    int *verbose,
    int *sync
    ){
      //parse options
  int argi;
  int unknown=0;

  //use current directory by default
  strcpy(sys_path,"./");
  strcpy(atm_path,"./");

  if(argc < 2){
    unknown=1;
  }

  for(argi = 1; argi < argc; argi++)
  {
    if(startswith( argv[argi], "--n_cores=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", ncores );
    }
    else
    if(startswith( argv[argi], "--n_gpus=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", ngpus );
    }
    else
    if(startswith( argv[argi], "--nact=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", nact );
      if(nact <= 0){
        fprintf(stderr,"Invalid number of actuators\n"); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--nmeas=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", nmeas );
      if(nmeas <= 0){
        fprintf(stderr,"Invalid number of total measurements\n"); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--nmeasts=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", nmeasts );
      if(nmeasts <= 0){
        fprintf(stderr,"Invalid number of measurements of the true sensor\n"); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--maxobs=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", maxobs );
      if(maxobs < 0){
        fprintf(stderr,"Invalid number of max obs\n"); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--maxrefine=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", maxrefine );
      if(maxrefine < 0){
        fprintf(stderr,"Invalid number of max refinement\n"); return 0;
      }
    }
    else
    if(startswith( argv[argi], "--tile=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%d", ts );
      if(ts <= 0){
        fprintf(stderr,"Invalid tile size\n"); return 0;
      }
    }

    else if(startswith( argv[argi], "--sys_path=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%s", sys_path );
    }

    else if(startswith( argv[argi], "--atm_path=")){
      sscanf( strchr( argv[argi], '=' ) + 1, "%s", atm_path );
    }
    else if(startswith( argv[argi], "--suffix=")){
        //--suffix start with --s :this condition must be tested before --s
        sscanf( strchr( argv[argi], '=' ) + 1, "%s", suffix);
    }
#ifdef USE_PLASMA
    else
    if(startswith( argv[argi], "--quark"))
      runtime = 0;
    else
    if(startswith( argv[argi], "--ompss"))
      runtime = 1;
#endif
    else
    if(startswith( argv[argi], "--v"))
      *verbose = 1;
    else
    if(startswith( argv[argi], "--s"))
      *sync = 1;
    else
    if(startswith( argv[argi], "--h"))
     unknown = 1;
    else
    {
      fprintf(stderr,"Unknown option %s, aborting...\n", argv[argi]); unknown=1;
    }
  }
    //defaults
  if(*ncores < 0)
  {
    get_thread_count(ncores);
  }
#ifdef USE_GPU
  fprintf(stdout,"GPU support enabled\n");
  if(*ngpus<0)
  {
    *ngpus=0;
  }
#else
    *ngpus=0;
#endif

  if(unknown!=0){
      USAGE;
  }
  return unknown;
}


int arg_isMissing(int nmeas, int nmeasts, int nact, int maxobs, int maxrefine, int ts){
    int stop=0;
  if(nmeas <= 0)
  {
    fprintf(stderr,"Please provide number of measurements\n");
    stop=1;
  }
  if(nmeasts <= 0)
  {
    fprintf(stderr,"Please provide number of measurements in the true sensor\n");
    stop=1;
  }
  if(nact <= 0)
  {
    fprintf(stderr,"Please provide number of actuators\n");
    stop=1;
  }
  if(maxrefine < 0)
  {
    fprintf(stderr,"Please provide number of max refine\n");
    stop=1;
  }
  if(maxobs < 0)
  {
    fprintf(stderr,"Please provide number of max obs\n");
    stop=1;
  }
  if(ts <= 0)
  {
    fprintf(stderr,"Please provide tile size\n");
    stop=1;
  }
  if(stop){
    USAGE;
    return 1;
  }
return 0;      
}

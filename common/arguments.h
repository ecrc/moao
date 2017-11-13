/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#ifndef ARGUMENTS_H
#define ARGUMENTS_H

#ifdef USE_PLASMA
#define HELP_PLASMA "--quark --ompss "
#else
#define HELP_PLASMA ""
#endif

#ifdef USE_GPU
#define HELP_GPU "--n_gpus=ngpus "
#else
#define HELP_GPU ""
#endif

#define HELP_MATCOV "--nmeas=nmeas --nmeasts=nbmeasts "



#define USAGE fprintf(stderr,"usage: --n_cores=ncores %s--tile=tile_size %s%s--maxrefine=maxrefine --maxobs=maxobs [--v] [--s] [--h]\nor\n"\
"usage:  --n_cores=ncores %s--tile=tile_size --sys_path=\"path\" --atm_path=\"path\" --suffix=\"suffix\" --maxrefine=maxrefine --maxobs=maxobs [--v] [--s] [--h]\n"\
"--n_cores: number of cores\n" \
"--n_gpus: number of GPUs\n" \
"--tile: tunable tile size\n" \
"--nact: number of actuators\n" \
"--maxrefine: max number of refinements\n" \
"--maxobs: max number of observations\n" \
"--quark: select QUARK runtime\n" \
"--ompss: select OmpSs runtime\n" \
"--sys_path: path to input files (system parameter)\n"\
"--atm_path: path to input files (atmosphere parameters)\n"\
"--suffix: suffix of input files\n"\
"--v: verbose\n"\
"--s: synchronize between stages, default off\n"\
"--h: print this help\n",HELP_GPU,HELP_MATCOV,HELP_GPU);

#ifdef __cplusplus
extern "C"{
#endif
int read_args(int argc, char **argv,
    int *n_cores, int *n_gpus, int *nact,  int *nmeas, int *nmeasts,
    int *maxobs, int *maxrefine, int *tile, char *sys_path, char *atm_path, char *suffix,
#ifdef USE_PLASMA
    int *runtime,
#endif
    int *verbose, int *sync);

int arg_isMissing(int nmeas, int nmeasts, int nact, int maxobs, int maxrefine, int tile);

#ifdef __cplusplus
}
#endif

#endif //ARGUMENTS_H_

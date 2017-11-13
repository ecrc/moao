/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>

#include <morse.h>
#include <starpu.h>
//#include <quark.h>

#include "arguments.h"
#include "utils.h"
#include "tomo_struct.h"

#include "flops.h"
#ifdef USE_INTERSAMPLE
#include "intersample.h"
#endif

#ifdef CHECK_CHAMELEON
#include "check_utils.h"
#endif



//main ===========================================================
int main( int argc, char** argv)
{
    char sys_path[100];
    char atm_path[100];
    char suffix[40]="";
    char dataFile[150]="";
    char isample_output_filename[256];

    int ncores = -1,ngpus=0, result=0;
    int maxobs = 1, maxrefine = 1;
    struct tomo_struct tomo;
#ifdef USE_INTERSAMPLE
    struct isample_struct isample;
#endif

  double alpha, beta;
  int nact;
  int nmeas=0, nmeasts=0, ts=0, verbose=0, sync=0, printhelp=0, unknown=0, runtime=0;
  double flops=0.0;
  double matcov_flops=0.0;


    //parse options
    read_args(argc,argv,&ncores,&ngpus,&nact,&nmeas,&nmeasts,&maxobs,&maxrefine,&ts,sys_path,atm_path,suffix,&verbose,&sync);

    MORSE_Init(ncores,ngpus);
    //TODO see settings: autotunning, tile size etc

    //get number of actuators: define by input matrix Dx
    long naxes[2];
    concatenateFileName(dataFile,sys_path,"Dx",suffix);
    getFitsDims(dataFile,naxes);
    if(naxes[0]>0){
        nact=naxes[0];
    }

    result = !matcov_init_tomo_tiled(&tomo, sys_path,atm_path, 0, 0, 0, -1, 0., 0.);
    HANDLE_ERROR_RET(result, "matcov_init_tomo_tiled");

    nmeas = matcov_getNumMeasurements(&tomo);
    nmeasts = matcov_getNumMeasurementsTS(&tomo);
    nmeas-=nmeasts;

#ifdef USE_INTERSAMPLE
    intersample_prepare(&isample, nact*nact, tomo.Nssp[tomo.Nw-1]+1, tomo.DiamTel, sys_path);
    //isample.np       = 16;
    intersample_init(&isample);
    double *cumulatedPSF=(double *)calloc(isample.N*isample.N,sizeof(double));
#endif

    arg_isMissing(nmeas,nmeasts,nact,maxobs,maxrefine,ts);

    if(verbose) fprintf(stdout, "\nProcessing %d number of measurements and %d number of actuators\n", nmeas, nact);
    if(verbose) fprintf(stdout, "Tile size (tunable): %d\n", ts);
    if(verbose) fprintf(stdout, "Working on: %d CPU cores and %d GPUs\n\n", ncores, ngpus);

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
    MORSE_desc_t *descCmm = NULL, *descCpm = NULL, *descCpp = NULL;

    double *R = NULL, *Cee = NULL, *Cvv = NULL, *Cvv_lap = NULL, *Dx = NULL, *Tmp = NULL;
    MORSE_desc_t *descR = NULL, *descCee = NULL, *descCvv = NULL, *descDx = NULL, *descTmp = NULL;

    MORSE_sequence_t *sequence;
    MORSE_request_t request[19] = { MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,
                                   MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,
                                   MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,
                                   MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,
                                   MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER};
    MORSE_Sequence_Create(&sequence);


    // Allocate matrix memory
    Cmm = (double*)calloc( (size_t)nmeas * nmeas , sizeof(double) );
    if ( ! Cmm ) {
      fprintf(stderr, "Out of Memory for Covariance Matrix (Cmm)\n");
      return -1;
    }
    Cpm = (double*)calloc( (size_t)nmeasts * nmeas , sizeof(double) );
    if ( ! Cpm ) {
      fprintf(stderr, "Out of Memory for Cpm Matrix\n");
      return -1;
    }
    Cpp = (double*)calloc( (size_t)nmeasts * nmeasts , sizeof(double) );
    if ( ! Cpp ) {
      fprintf(stderr, "Out of Memory for Cpp Matrix\n");
      return -1;
    }
    R = (double*)calloc( (size_t)nmeasts * nmeas , sizeof(double) );
    if ( ! R ) {
      fprintf(stderr, "Out of Memory for ToR (R)\n");
      return -1;
    }
    Cee = (double*)calloc( (size_t)nmeasts * nmeasts , sizeof(double) );
    if ( ! Cee ) {
      fprintf(stderr, "Out of Memory for Cee Matrix\n");
      return -1;
    }
    Cvv = (double*)calloc( (size_t)nact * nact , sizeof(double) );
    if ( ! Cvv ) {
      fprintf(stderr, "Out of Memory for Cvv Matrix\n");
      return -1;
    }
    Cvv_lap = (double*)calloc( (size_t)nact * nact , sizeof(double) );
    if ( ! Cvv_lap ) {
      fprintf(stderr, "Out of Memory for Cvv_lap Matrix\n");
      return -1;
    }
    Dx = (double*)calloc( (size_t)nact * nmeasts , sizeof(double) );
    if ( ! Dx ) {
      fprintf(stderr, "Out of Memory for Dx Matrix\n");
      return -1;
    }
    Tmp = (double*)calloc( (size_t)nact * nmeasts , sizeof(double) );
    if ( ! Tmp ) {
      fprintf(stderr, "Out of Memory for Tmp Matrix\n");
      return -1;
    }

    double t_Cmm,t_Ctt,t_Ctm,t_R,t_CeeCvv;
    #ifdef CHECK_CHAMELEON
    double maxerr,errv1,errv2;
    char logName[50];
    strcpy(logName,"log");
    strcat(logName,suffix);
    strcat(logName,".txt");
    FILE *logFile=fopen(logName,"w");
    if(logFile==NULL){
      fprintf(stderr,"could not open log file\n");
    }
    fprintf(logFile,"nmeas: %d\nnmeasts: %d\nnatc: %d\nnLayers: %d\ntile size: %d\nncores: %d\nngpus: %d\npsf:0\n",nmeas,nmeasts,nact,tomo.Nlayer,ts,ncores,ngpus);
    fprintf(logFile,"Matrices times(P)   max_err(P/Y) err_val(P) err_val(Y)\n");
    double *lTmp=(double*)malloc((long)nmeas*(long)nmeas*sizeof(double));
    #endif

    MORSE_Desc_Create(&descCmm, Cmm, MorseRealDouble, ts, ts, ts*ts, nmeas, nmeas, 0, 0, nmeas, nmeas,1,1);
    MORSE_Desc_Create(&descCpm, Cpm, MorseRealDouble, ts, ts, ts*ts, nmeasts, nmeas, 0, 0, nmeasts, nmeas,1,1);
    MORSE_Desc_Create(&descCpp, Cpp, MorseRealDouble, ts, ts, ts*ts, nmeasts, nmeasts, 0, 0, nmeasts, nmeasts,1,1);
    MORSE_Desc_Create(&descR,   R,   MorseRealDouble, ts, ts, ts*ts, nmeasts, nmeas, 0, 0, nmeasts, nmeas,1,1);
    MORSE_Desc_Create(&descCee, Cee, MorseRealDouble, ts, ts, ts*ts, nmeasts, nmeasts, 0, 0, nmeasts, nmeasts,1,1);
    MORSE_Desc_Create(&descCvv, Cvv, MorseRealDouble, ts, ts, ts*ts, nact, nact, 0, 0, nact, nact,1,1);
    MORSE_Desc_Create(&descDx,  Dx,  MorseRealDouble, ts, ts, ts*ts, nact, nmeasts, 0, 0, nact, nmeasts,1,1);
    MORSE_Desc_Create(&descTmp, Tmp, MorseRealDouble, ts, ts, ts*ts, nact, nmeasts, 0, 0, nact, nmeasts,1,1);

    if(verbose) fprintf(stdout, "CHAMELEON: MOAO started\n\n");
    if(verbose) fprintf(stdout, "CHAMELEON: MOAO Loading the interaction matrix Dx\n");


    preprocess_cpu_time -= cWtime();
#ifdef USE_INTERSAMPLE
    double *yTmp=(double*)malloc((long)nmeas*(long)nmeas*sizeof(double));
    concatenateFileName(dataFile,sys_path,"Dx",suffix);
    readFits(dataFile,yTmp);
    MORSE_Lapack_to_Tile(yTmp,nact,descDx);
#endif
    if(sync){
        MORSE_Sequence_Wait(sequence);
    }
    preprocess_cpu_time += cWtime();
    if(verbose) fprintf(stdout, "Done in %f(s)\n\n", preprocess_cpu_time);
    total_time += preprocess_cpu_time;
    
    if(verbose) fprintf(stdout, "CHAMELEON: MOAO main computation started with\n");
    if(verbose) fprintf(stdout, "Maximum number of refinements %d\n", maxrefine);
    if(verbose) fprintf(stdout, "Maximum number of observations %d\n", maxobs);

    START_TIMING(compute_time);
// add loop over max_nbrefine
    for(nbrefine=0;nbrefine<maxrefine;nbrefine++){
        if(maxrefine>0){
            START_TIMING(refine_time);

            if(verbose) fprintf(stdout, "\nnbrefine %d\n", nbrefine+1);

            if(sync){
                MORSE_Sequence_Wait(sequence);
                START_TIMING(matcov_time);
            }

            //need to get new atm-params
            int night_idx = 0,
            snapshots_per_night = 8,
            snapshot_idx = nbrefine,
            obs_idx = -1;
            result = !matcov_update_atm_params(&tomo, night_idx, snapshots_per_night, snapshot_idx, obs_idx);
            HANDLE_ERROR(result, "matcov_update_atm_params");
    
            START_TIMING(t_Cmm);
            APPNAME_matcov_tile(descCmm, sequence,&request[0],&tomo, 1);
            flops = flops + FLOPS_DCovMat(nmeas, nmeas, tomo.Nlayer, 1);
            matcov_flops = matcov_flops + FLOPS_DCovMat(nmeas, nmeas, tomo.Nlayer, 1);
            if(sync){
                MORSE_Sequence_Wait(sequence);
                STOP_TIMING(t_Cmm);
            }

            START_TIMING(t_Ctm);
            APPNAME_matcov_tile(descR, sequence,&request[1],&tomo, 3);
            flops = flops + FLOPS_DCovMat(nmeasts, nmeas, tomo.Nlayer, 3);
            matcov_flops = matcov_flops + FLOPS_DCovMat(nmeasts, nmeas, tomo.Nlayer, 3);
            if(sync){
                MORSE_Sequence_Wait(sequence);
                STOP_TIMING(matcov_time);
                total_matcov_time += matcov_time;
                STOP_TIMING(t_Ctm);
            }

            /* Wait for the call to complete.  */
            #pragma starpu wait
            fprintf(stdout,"matcov (Cmm): in %fs: %fGFLOPS\n",t_Cmm,matcov_flops*1.e-9/t_Cmm);
            fprintf(stdout,"matcov (Ctm): in %fs: %fGFLOPS\n",t_Ctm,matcov_flops*1.e-9/t_Ctm);
            fflush(stdout);
    


#ifdef CHECK_CHAMELEON
///////////////////////////////////////////////////////////////
// import data from yorick files
            concatenateFileName(dataFile,sys_path,"cmm",suffix);
            readFits(dataFile,yTmp);
            MORSE_Tile_to_Lapack(descCmm,lTmp,nmeas);
            MORSE_Sequence_Wait(sequence);
            //fprintf(stdout,"\tError on Cmm (Yorick) %f\n\n",compareMatrices(yCmm,yTmp,nmeas,nmeas,'U','N'));
            compareMatrices2(lTmp,yTmp,nmeas,nmeas,'U','N',&maxerr,&errv1,&errv2);
            fprintf(logFile,"Cmm      %e   %e     %e   %e\n",t_Cmm,maxerr,errv1,errv2);
    
            concatenateFileName(dataFile,sys_path,"ctm",suffix);
            readFits(dataFile,yTmp);
            MORSE_Tile_to_Lapack(descR,lTmp,nmeasts);
            MORSE_Sequence_Wait(sequence);
            //fprintf(stdout,"\tError on Ctm (Yorick) %f\n\n",compareMatrices(yCmm,yTmp,nmeasts,nmeas,'A','N'));
            compareMatrices2(lTmp,yTmp,nmeasts,nmeas,'A','N',&maxerr,&errv1,&errv2);
            fprintf(logFile,"Ctm      %e   %e     %e   %e\n",t_Ctm,maxerr,errv1,errv2);
///////////////////////////////////////////////////////////////
#endif

            START_TIMING(t_R);
            APPNAME_reconstructor(descCmm,descR,sequence,0);
            flops = flops + FLOPS_DPOTRF(nmeas);
            flops = flops + FLOPS_DTRSM(MorseRight, nmeasts, nmeas);
            flops = flops + FLOPS_DTRSM(MorseRight, nmeasts, nmeas);
            if(sync){
                MORSE_Sequence_Wait(sequence);
                STOP_TIMING(t_R);
                long nFlops=FLOPS_DPOTRF(nmeas)+FLOPS_DTRSM(MorseRight, nmeasts, nmeas)+FLOPS_DTRSM(MorseRight, nmeasts, nmeas);
                fprintf(stdout,"reconstructor: %f, %ld perf: %f\n",t_R,nFlops, nFlops/1e9/t_R);
            }

#ifdef CHECK_CHAMELEON
            concatenateFileName(dataFile,sys_path,"R",suffix);
            readFits(dataFile,yTmp);
            MORSE_Tile_to_Lapack(descR,lTmp,nmeasts);
            MORSE_Sequence_Wait(sequence);
            compareMatrices2(lTmp,yTmp,nmeasts,nmeas,'A','N',&maxerr,&errv1,&errv2);
            fprintf(logFile,"R        %e   %e     %e   %e\n",t_R,maxerr,errv1,errv2);
            fprintf(stdout,"RECONSTRUCTOR DONE\n");
#endif
    
    
            STOP_TIMING(refine_time);
            total_refine_time += refine_time;
        }

        for (nbobs = 0; nbobs < maxobs; nbobs++) {
            START_TIMING(obs_time);
            if(verbose) fprintf(stdout, " nbobs %d\n", nbobs+1);
        
            if(sync){
                MORSE_Sequence_Wait(sequence);
                START_TIMING(matcov_time);
            }

            START_TIMING(t_Cmm);
            APPNAME_matcov_tile(descCmm, sequence,&request[2],&tomo, 1);
            flops = flops + FLOPS_DCovMat(nmeas, nmeas, tomo.Nlayer, 1);
            if(sync){
                MORSE_Sequence_Wait(sequence);
                STOP_TIMING(t_Cmm);
            }

            START_TIMING(t_Ctm);
            APPNAME_matcov_tile(descCpm, sequence,&request[3],&tomo, 3);
            flops = flops + FLOPS_DCovMat(nmeasts, nmeas, tomo.Nlayer, 3);
            if(sync){
                MORSE_Sequence_Wait(sequence);
                STOP_TIMING(t_Ctm);
            }

            START_TIMING(t_Ctt);
            APPNAME_matcov_tile(descCpp, sequence,&request[4],&tomo, 4);
            flops = flops + FLOPS_DCovMat(nmeasts, nmeasts, tomo.Nlayer, 4);
            if(sync){
                MORSE_Sequence_Wait(sequence);
                STOP_TIMING(matcov_time);
                total_matcov_time += matcov_time;

                STOP_TIMING(t_Ctt);
            }
        
#ifdef CHECK_CHAMELEON
///////////////////////////////////////////////////////////////
            fprintf(logFile,"Cmm      %e   x               x              x \n",t_Cmm);
            fprintf(logFile,"Ctm      %e   x               x              x \n",t_Ctm);
            concatenateFileName(dataFile,sys_path,"ctt",suffix);
            readFits(dataFile,yTmp);
            MORSE_Tile_to_Lapack(descCpp,lTmp,nmeasts);
            MORSE_Sequence_Wait(sequence);
            compareMatrices2(lTmp,yTmp,nmeasts,nmeasts,'A','N',&maxerr,&errv1,&errv2);
            fprintf(logFile,"Ctt      %e   %e     %e   %e\n",t_Ctt,maxerr,errv1,errv2);
///////////////////////////////////////////////////////////////
#endif
        
        
           START_TIMING(t_CeeCvv);
            MORSE_Cee_Cvv(descCmm, descCpp, descCpm, descR, descDx,
                          descCee, descCvv, descTmp, sequence);
            flops = flops + FLOPS_DSYR2K(nmeasts, nmeas);
            flops = flops + FLOPS_DSYMM(MorseRight, nmeasts, nmeas);
            flops = flops + FLOPS_DGEMM(nmeasts, nmeasts, nmeas);
            flops = flops + FLOPS_DSYMM(MorseRight, nact, nmeasts);
            flops = flops + FLOPS_DGEMM(nact, nact, nmeasts);
            if(sync){
                MORSE_Sequence_Wait(sequence);
                STOP_TIMING(t_CeeCvv);
            }
        
#ifdef CHECK_CHAMELEON
///////////////////////////////////////////////////////////////
// compare Cee and Cvv
            concatenateFileName(dataFile,sys_path,"cee",suffix);
            readFits(dataFile,yTmp);
            MORSE_Sequence_Wait(sequence);
            MORSE_Tile_to_Lapack(descCee,lTmp,nmeasts);
            compareMatrices2(lTmp,yTmp,nmeasts,nmeasts,'A','N',&maxerr,&errv1,&errv2);
            fprintf(logFile,"Cee      %e   %e     %e   %e\n",t_CeeCvv,maxerr,errv1,errv2);
        
            concatenateFileName(dataFile,sys_path,"cvv",suffix);
            readFits(dataFile,yTmp);
            MORSE_Sequence_Wait(sequence);
            MORSE_Tile_to_Lapack(descCvv,lTmp,nact);
            compareMatrices2(lTmp,yTmp,nact,nact,'A','N',&maxerr,&errv1,&errv2);
            fprintf(logFile,"Cvv      %e   %e     %e   %e\n",t_CeeCvv,maxerr,errv1,errv2);
///////////////////////////////////////////////////////////////
#endif
        

#ifdef USE_INTERSAMPLE
            MORSE_Sequence_Wait(sequence);
            MORSE_Tile_to_Lapack(descCvv,Cvv_lap,nact);
        
            MORSE_Sequence_Wait(sequence);
            START_TIMING(intersample_time);

            //save short exposure psf to fits file
            intersample_process(&isample, Cvv_lap);
            sprintf(isample_output_filename, "psf_oloop%d_iloop%d_gal%d.fits", nbobs, nbrefine, 1);
            intersample_save(&isample, isample_output_filename);

            flops = flops + FLOPS_Intersample(nact, nact);
            if(sync){
                MORSE_Sequence_Wait(sequence);
                STOP_TIMING(intersample_time);
                total_intersample_time += intersample_time;
                STOP_TIMING(obs_time);
                total_obs_time += obs_time;
            }
            //cumulate all the psf
            int i;
            for(i=0;i<isample.N*isample.N;i++){
                cumulatedPSF[i]+=isample.dphi[i];
            }

        
#ifdef CHECK_CHAMELEON
            double *Ypsf=(double*)calloc((size_t)isample.N*isample.N, sizeof(double));
            concatenateFileName(dataFile,sys_path,"psf",suffix);
            FILE *FileExists=fopen(dataFile,"r");
            if(FileExists!=NULL){
                fclose(FileExists);
                readFits(dataFile,Ypsf);
            }
            compareMatrices2(isample.dphi,Ypsf,isample.N,isample.N,'A','T',&maxerr,&errv1,&errv2);
            fprintf(logFile,"psf      %e   %e     %e   %e\n",intersample_time,maxerr,errv1,errv2);
	    free(Ypsf);
#endif
#endif
        }//end loop on nbobs
    } //end loop on nbrefine

#ifdef USE_INTERSAMPLE
//normalize the cumulated psf
    int i;
    for(i=0;i<isample.N*isample.N;i++){
        cumulatedPSF[i]/=(maxobs*nbrefine);
    }
    char psf_file[256];
    sprintf(psf_file, "psf_nbrefine%d_maxobs%d.fits", nbrefine,maxobs);
    writeFits(psf_file,isample.N,isample.N,cumulatedPSF);
#endif


    MORSE_Sequence_Wait(sequence);

    STOP_TIMING(compute_time);
    total_time += compute_time;

    MORSE_Sequence_Destroy(sequence);

    MORSE_Desc_Destroy( &descCmm );
    MORSE_Desc_Destroy( &descCpm );
    MORSE_Desc_Destroy( &descCpp );
    MORSE_Desc_Destroy( &descR );
    MORSE_Desc_Destroy( &descCee );
    MORSE_Desc_Destroy( &descCvv );
    MORSE_Desc_Destroy( &descDx );
    MORSE_Desc_Destroy( &descTmp );

    matcov_free_tomo_tiled(&tomo);
#ifdef CHECK_CHAMELEON
    free(lTmp);
#endif

#ifdef USE_INTERSAMPLE
    intersample_free(&isample);
#endif




    if(result != MORSE_SUCCESS) {
        fprintf(stdout,"An error occurred! sorry... %d\n\n", result);
    }
    else {
        if(verbose) {
            fprintf(stdout, "Done in %f(s)\n", compute_time);
            fprintf(stdout, "Computation successful :-)\n\n");
        }
        fprintf(stdout," Cores  Tile  Refine  Obs  Total   Total  Total MOAO alloc MOAO preprocess  Matcov   MOAO refine    MOAO obs    MOAO PSF     MOAO compute |   Total    |  Gflops/s\n");
        fprintf(stdout," #cores size  #iter  #iter #Meas #Meas-ts #Actu  time (s)   time (s)        time (s)  time (s)      time(s)     time(s)        time (s)   |   time(s)  |  Perf\n");
        fprintf(stdout," ====== ====  ====== ===== ===== ======== ===== ========== =============== ========= =========== ============  ==========   ============= | ========== | ========\n");
        fprintf(stdout,"    %d   %d    %d     %d    %-8d%-8d%-4d    %-11.3f %-10.3f   %-10.3f%-10.3f   %-10.3f   %-10.3f      %-4.3f        %-10.3f   %-4.3f\n\n", ncores, ts, maxrefine, maxobs, nmeas, nmeasts, nact, alloc_cpu_time, preprocess_cpu_time, total_matcov_time, total_refine_time, total_obs_time, total_intersample_time, compute_time, total_time, flops / 1e9 / total_time);
    }


    MORSE_Finalize();
    return 0;
}



//#ifdef USE_INTERSAMPLE
//void CORE_dIntersample_genPSF_starpu(void *buffers[],void *cl_arg){
//    char isample_output_filename[256];
//    struct isample_struct isample;
//    double *Cvv;
//    int oloop,iloop,gal;
//    starpu_codelet_unpack_args(cl_arg,&isample,&Cvv,&oloop,&iloop,&gal);
//  
//    intersample_process(&isample, Cvv);
//    sprintf(isample_output_filename, "psf_oloop%d_iloop%d_gal%d.fits", oloop, iloop, gal);
//    intersample_save(&isample, isample_output_filename);
//}
//static struct starpu_codelet cl_Intersample =
//{
//    .where = STARPU_CPU,
//    .cpu_funcs = {CORE_dIntersample_genPSF_starpu},
//    .nbuffers = 0
//};
//int MORSE_Intersample_genPSF(double *Cvv, struct isample_struct *isample,int oloop, int iloop, int gal, MORSE_sequence_t *sequence, MORSE_request_t *request){
//    MORSE_context_t *morse;
//    MORSE_option_t options;
//    morse = morse_context_self();
//    if (sequence->status != MORSE_SUCCESS)
//        return;
//    RUNTIME_options_init(&options, morse, sequence, request);
//
//    struct starpu_codelet *cl=&cl_Intersample;
//    starpu_insert_task(cl,
//        STARPU_VALUE, isample,  sizeof(struct isample_struct),
//        STARPU_VALUE, &Cvv,     sizeof(double*),
//        STARPU_VALUE, &oloop,    sizeof(int),
//        STARPU_VALUE, &iloop,    sizeof(int),
//        STARPU_VALUE, &gal,      sizeof(int),
//        0);
//
//    RUNTIME_options_finalize(&options, morse);
//    MORSE_TASK_dataflush_all();
//
//    return 0;
//}
//#endif



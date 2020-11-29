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


#include "flops.h"
#include "moao_chameleon.hpp"


std::list<Arg> args={Arg("sys_path"     ,"path to the system parameter file"            ,""),
                     Arg("atm_path"     ,"path to the atmospheric parameter file"       ,""),
                     Arg("nmeas"        ,"number of sensor measurements"                ,"0" ),
                     Arg("nmeasTS"      ,"number of measurements for the truth sensor"  ,"0" ),
                     Arg("nact"         ,"number of actuators"                          ,"0" ),
                     Arg("maxobs"       ,"number of observation"                        ,"1" ),
                     Arg("maxrefine"    ,"number of refinment per observation"          ,"1" ),
                     Arg("snapPerNight" ,"number of snapshot per night"                 ,"8" ),
                     Arg("ncores"       ,"number of cores for the execution"            ,"1" ),
                     Arg("ngpus"        ,"number of gpus for the execution"             ,"0" ),
                     Arg("ts"           ,"tile size"                                    ,"200" ),
                     Arg("sync"         ,"activate synchronization between steps"       ,"0" ),
                     Arg("verbose"      ,""                                             ,"0" )
                     };
                     //TODO tile size, ncores,ngpus...
CmdLine cmd(args);

using real_t=float;

//main ===========================================================
int main( int argc, char** argv)
{

    cmd.parse(argc,argv);
    std::string sys_path=cmd.getString("sys_path");
    std::string atm_path=cmd.getString("atm_path");
    int nmeas           = cmd.getInt("nmeas");
    int nmeasts         = cmd.getInt("nmeasTS");
    int nact            = cmd.getInt("nact");
    int maxobs          = cmd.getInt("maxobs");
    int maxrefine       = cmd.getInt("maxrefine");
    int snapPerNight    = cmd.getInt("snapPerNight");
    int ncores          = cmd.getInt("ncores");
    int ngpus           = cmd.getInt("ngpus");
    int ts              = cmd.getInt("ts");
    int sync            = cmd.getInt("sync"); 
    int verbose         = cmd.getInt("verbose");
    cmd.recap(true);

    int obsIdx=-1;
    int nightIdx=0;
    int snapIdx=0;

    int result;
    double matcov_flops=0., reconstructor_flops=0., CeeCvv_flops=0., psf_flops=0.;
    Timer t_matcov, t_R, t_CeeCvv, t_psf, t_preprocess, t_alloc, t_glob, t_obs, t_refine;

    t_glob.start();
    t_preprocess.start();
    std::string psfFileName;
    std::string atmFileName ="/prof"+std::to_string(obsIdx)+"-atmos-night"+std::to_string(nightIdx)+".txt";
    Tomo_struct<real_t> tomo(sys_path+"/sys-params.txt",atm_path+atmFileName);

#ifdef USE_INTERSAMPLE
    struct isample_struct isample;
    //get number of actuators: define by input matrix Dx
    long naxes[2];
    std::string DxFile=sys_path+"/Dx.fits";
    getFitsDims(DxFile.c_str(),naxes);
    //get number of actuators: define by input matrix Dx //TODO nact>0
    nact=naxes[0];
    intersample_prepare(&isample, nact*nact, tomo.Nssp[tomo.sys.nW-1]+1, tomo.sys.diam, sys_path.c_str());
    intersample_init(&isample);
#endif

    nmeas   =tomo.getNMeas();
    nmeasts =tomo.getNMeasTS();
    std::cout<<"using "<<nmeas<<" meas and "<<nmeasts<<" ts meas"<<std::endl;
    std::cout<<"Using "<<nact<<" actuators"<<std::endl;
    std::cout<<"Tile size:"<<ts<<std::endl;
    std::cout<<"Working on "<<ncores<<" cores and "<<ngpus<<" gpus"<<std::endl;
    t_preprocess.stop();


    MORSE_Init(ncores,ngpus);
    //TODO see settings: autotunning, tile size etc

    double total_time = 0.0;
    int nbobs, nbrefine;

    t_alloc.start();
    real_t *Cmm = NULL, *Cpm = NULL, *Cpp = NULL;
    MORSE_desc_t *descCmm = NULL, *descCpm = NULL, *descCpp = NULL;

    real_t *R = NULL, *Cee = NULL, *Cvv = NULL, *Cvv_lap = NULL, *Dx = NULL, *Tmp = NULL;
    MORSE_desc_t *descR = NULL, *descCee = NULL, *descCvv = NULL, *descDx = NULL, *descTmp = NULL;

    MORSE_sequence_t *sequence;
    MORSE_request_t request[19] = { MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,
                                   MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,
                                   MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,
                                   MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,
                                   MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER};
    MORSE_Sequence_Create(&sequence);


    // Allocate matrix memory
    Cmm = (real_t*)calloc( (size_t)nmeas * nmeas , sizeof(real_t) );
    if ( ! Cmm ) {
      fprintf(stderr, "Out of Memory for Covariance Matrix (Cmm)\n");
      return -1;
    }
    Cpm = (real_t*)calloc( (size_t)nmeasts * nmeas , sizeof(real_t) );
    if ( ! Cpm ) {
      fprintf(stderr, "Out of Memory for Cpm Matrix\n");
      return -1;
    }
    Cpp = (real_t*)calloc( (size_t)nmeasts * nmeasts , sizeof(real_t) );
    if ( ! Cpp ) {
      fprintf(stderr, "Out of Memory for Cpp Matrix\n");
      return -1;
    }
    R = (real_t*)calloc( (size_t)nmeasts * nmeas , sizeof(real_t) );
    if ( ! R ) {
      fprintf(stderr, "Out of Memory for ToR (R)\n");
      return -1;
    }
    Cee = (real_t*)calloc( (size_t)nmeasts * nmeasts , sizeof(real_t) );
    if ( ! Cee ) {
      fprintf(stderr, "Out of Memory for Cee Matrix\n");
      return -1;
    }
    Cvv = (real_t*)calloc( (size_t)nact * nact , sizeof(real_t) );
    if ( ! Cvv ) {
      fprintf(stderr, "Out of Memory for Cvv Matrix\n");
      return -1;
    }
    Cvv_lap = (real_t*)calloc( (size_t)nact * nact , sizeof(real_t) );
    if ( ! Cvv_lap ) {
      fprintf(stderr, "Out of Memory for Cvv_lap Matrix\n");
      return -1;
    }
    Dx = (real_t*)calloc( (size_t)nact * nmeasts , sizeof(real_t) );
    if ( ! Dx ) {
      fprintf(stderr, "Out of Memory for Dx Matrix\n");
      return -1;
    }
    Tmp = (real_t*)calloc( (size_t)nact * nmeasts , sizeof(real_t) );
    if ( ! Tmp ) {
      fprintf(stderr, "Out of Memory for Tmp Matrix\n");
      return -1;
    }

    MORSE_Desc_Create(&descCmm, Cmm, MORSE<real_t>::real, ts, ts, ts*ts, nmeas  , nmeas  , 0, 0, nmeas  , nmeas  ,1,1);
    MORSE_Desc_Create(&descCpm, Cpm, MORSE<real_t>::real, ts, ts, ts*ts, nmeasts, nmeas  , 0, 0, nmeasts, nmeas  ,1,1);
    MORSE_Desc_Create(&descCpp, Cpp, MORSE<real_t>::real, ts, ts, ts*ts, nmeasts, nmeasts, 0, 0, nmeasts, nmeasts,1,1);
    MORSE_Desc_Create(&descR  , R,   MORSE<real_t>::real, ts, ts, ts*ts, nmeasts, nmeas  , 0, 0, nmeasts, nmeas  ,1,1);
    MORSE_Desc_Create(&descCee, Cee, MORSE<real_t>::real, ts, ts, ts*ts, nmeasts, nmeasts, 0, 0, nmeasts, nmeasts,1,1);
    MORSE_Desc_Create(&descCvv, Cvv, MORSE<real_t>::real, ts, ts, ts*ts, nact   , nact   , 0, 0, nact   , nact   ,1,1);
    MORSE_Desc_Create(&descDx , Dx,  MORSE<real_t>::real, ts, ts, ts*ts, nact   , nmeasts, 0, 0, nact   , nmeasts,1,1);
    MORSE_Desc_Create(&descTmp, Tmp, MORSE<real_t>::real, ts, ts, ts*ts, nact   , nmeasts, 0, 0, nact   , nmeasts,1,1);
    t_alloc.stop();
    if(verbose){
        std::cout<<"CHAMELEON: MOAO started\n\n"<<
                "CHAMELEON: MOAO Loading the interaction matrix Dx"<<std::endl;
    }


#ifdef USE_INTERSAMPLE
    t_preprocess.start();
    //read Dx from file
    readFits(DxFile.c_str(), Tmp);
    MORSE_Lapack_to_Tile(Tmp,nact,descDx);
    t_preprocess.stop();
#endif
    if(sync){
        MORSE_Sequence_Wait(sequence);
    }

    if(verbose){
        std::cout<<"CHAMELEON: MOAO main computation started with\n"<<
                "Maximum number of refinements "<< maxrefine<<"\n"<<
                "Maximum number of observations " << maxobs<<std::endl;
    }

// add loop over max_nbrefine
    for(nbrefine=0;nbrefine<maxrefine;nbrefine++){
        if(maxrefine>0){
            t_refine.start();

            if(verbose) std::cout<<"\nnbrefine"<< nbrefine+1<<std::endl;;

            if(sync){
                MORSE_Sequence_Wait(sequence);
            }

            //need to get new atm-params
            nightIdx=0;
            snapIdx=nbrefine;
            obsIdx=0;
            //update atm params
            atmFileName ="/prof"+std::to_string(snapIdx * snapPerNight + obsIdx)+"-atmos-night"+std::to_string(nightIdx)+".txt";
            tomo.updateAtm(atm_path+atmFileName);
    
            //----------------------
            //generating Cmm
            if(sync){t_matcov.start();}
            MOAO_CHAMELEON::matcov_tile(descCmm, sequence,&request[0],&tomo, 1);
            matcov_flops = matcov_flops + FLOPS_DCovMat(nmeas, nmeas, tomo.atm.nLayer, 1);

            //----------------------
            //generating Ctm(s)
            MOAO_CHAMELEON::matcov_tile(descR, sequence,&request[1],&tomo, 3);
            matcov_flops = matcov_flops + FLOPS_DCovMat(nmeasts, nmeas, tomo.atm.nLayer, 3);
            if(sync){
                MORSE_Sequence_Wait(sequence);
                t_matcov.stop();
            }
            /* Wait for the call to complete.  */
            #pragma starpu wait

            //----------------------
            //computing reconstructor
            if(sync){t_R.start();}
            MOAO_CHAMELEON::reconstructor<real_t>(descCmm,descR,sequence,0);
            reconstructor_flops = FLOPS_reconstructor(nmeas,nmeasts);
            if(sync){
                MORSE_Sequence_Wait(sequence);
                t_R.stop();
            }

            t_refine.stop();
        }

        t_obs.start();
        for (nbobs = 0; nbobs < maxobs; nbobs++) {
            if(verbose) std::cout<<" nbobs "<<nbobs+1<<std::endl;
            //need to get new atm-params
            nightIdx=0;
            snapIdx=nbrefine;
            obsIdx=nbobs;
            //update atm params
            atmFileName ="/prof"+std::to_string(snapIdx * snapPerNight + obsIdx)+"-atmos-night"+std::to_string(nightIdx)+".txt";
            tomo.updateAtm(atm_path+atmFileName);
        
            if(sync){
                MORSE_Sequence_Wait(sequence);
                t_matcov.start();
            }
            //----------------------
            //generating Cmm (destroyed by reconstructor computation)
            MOAO_CHAMELEON::matcov_tile(descCmm, sequence,&request[2],&tomo, 1);
            matcov_flops = matcov_flops + FLOPS_DCovMat(nmeas, nmeas, tomo.atm.nLayer, 1);
            if(sync){
                MORSE_Sequence_Wait(sequence);
                t_matcov.stop();
            }

            // TODO for(int nbgal=0; nbgal<tomo.sys.nTarget; nbgal++){
            // TODO     if(verbose) std::cout<<"  target "<<nbgal<<std::endl;
            // TODO     tomo.setTarget(nbgal);
                if(sync){
                    MORSE_Sequence_Wait(sequence);
                    t_matcov.start();
                }
                //----------------------
                //generating Ctm (destroyed by reconstructor computation)
                MOAO_CHAMELEON::matcov_tile(descCpm, sequence,&request[3],&tomo, 3);
                matcov_flops = matcov_flops + FLOPS_DCovMat(nmeasts, nmeas, tomo.atm.nLayer, 3);

                //----------------------
                //generating Ctt
                MOAO_CHAMELEON::matcov_tile(descCpp, sequence,&request[4],&tomo, 4);
                matcov_flops = matcov_flops + FLOPS_DCovMat(nmeasts, nmeasts, tomo.atm.nLayer, 4);
                if(sync){
                    MORSE_Sequence_Wait(sequence);
                    t_matcov.stop();
                    t_CeeCvv.start();
                }
        
                //----------------------
                //generating error and voltage
                std::cerr<<"Cee_Cvv"<<std::endl;
                MOAO_CHAMELEON::Cee_Cvv<real_t>(descCmm, descCpp, descCpm, descR, descDx,
                              descCee, descCvv, descTmp, sequence);
                CeeCvv_flops += FLOPS_CeeCvv(nmeas,nmeasts,nact);
                if(sync){
                    MORSE_Sequence_Wait(sequence);
                    t_CeeCvv.stop();
                std::cerr<<"done"<<std::endl;
                }
        

#ifdef USE_INTERSAMPLE
                MORSE_Sequence_Wait(sequence);
                MORSE_Tile_to_Lapack(descCvv,Cvv_lap,nact);
        
                MORSE_Sequence_Wait(sequence);
                t_psf.start();

                //save short exposure psf to fits file
                intersample_process(&isample, Cvv_lap);
                psfFileName="psf_oloop"+std::to_string(nbrefine)+"_iloop"+std::to_string(nbobs)+"_gal"+std::to_string(0)+".fits";
                writeFits(psfFileName.c_str(), isample.N, isample.N, isample.dphi);

                psf_flops = psf_flops + FLOPS_Intersample(nact, nact);
                if(sync){
                    MORSE_Sequence_Wait(sequence);
                    t_psf.stop();
                }
                //cumulate all the psf
                int i;
                for(i=0;i<isample.N*isample.N;i++){
                //TODO    cumulatedPSF[i]+=isample.dphi[i];
                }
#endif
            //TODO} //End of galaxies loop
            if(verbose) std::cout<<"-----end obs\n"<<std::endl;
        }//End of OBS loop
        t_obs.stop();
        if(verbose) std::cout<<"=====end refine\n"<<std::endl;
    }//End of REFINE loop

#ifdef USE_INTERSAMPLE
//normalize the cumulated psf
    int i;
    for(i=0;i<isample.N*isample.N;i++){
        //TODO cumulatedPSF[i]/=(maxobs*nbrefine);
    }
    psfFileName="psf_nbrefine"+std::to_string(nbrefine)+"_maxobs"+std::to_string(maxobs)+".fits";
    //writeFits(psfFileName.c_str(), isample.N, isample.N, cumulatedPSF);
#endif


    MORSE_Sequence_Wait(sequence);

    MORSE_Sequence_Destroy(sequence);

    MORSE_Desc_Destroy( &descCmm );
    MORSE_Desc_Destroy( &descCpm );
    MORSE_Desc_Destroy( &descCpp );
    MORSE_Desc_Destroy( &descR );
    MORSE_Desc_Destroy( &descCee );
    MORSE_Desc_Destroy( &descCvv );
    MORSE_Desc_Destroy( &descDx );
    MORSE_Desc_Destroy( &descTmp );

#ifdef USE_INTERSAMPLE
    intersample_free(&isample);
#endif
    t_glob.stop();




    //if(result != MORSE_SUCCESS) {
    //    fprintf(stdout,"An error occurred! sorry... %d\n\n", result);
    //}
    //else {
    std::cout<<"========================================"
            <<"\nsimu characteristics"
            <<"\nncores   : "<<ncores
            <<"\nngpus    : "<<ngpus
            <<"\ntile size: "<<ts   
            <<"\nnmeas    : "<<nmeas
            <<"\nnmeasts  : "<<nmeasts
            <<"\nnact     : "<<nact
            <<"\nnLayers  : "<<tomo.atm.nLayer
            <<"\nnTargets : "<<tomo.sys.nTarget
            <<"\nnrefine  : "<<maxrefine
            <<"\nnObs     : "<<maxobs
            <<"\nnTarget  : "<<tomo.sys.nTarget
            <<"\n========================================"
            <<std::endl;

        double compute_time=t_R.elapsed()+t_CeeCvv.elapsed();
        double flops=matcov_flops+reconstructor_flops+CeeCvv_flops+psf_flops;
        std::cout<<"Performances:  time (s)   GFlop/s"
                <<"\nallocation     "<<t_alloc.elapsed()
                <<"\npreprocess     "<<t_preprocess.elapsed()
                <<std::endl;
        if(sync){
        std::cout<<"matcov         "<<t_matcov.elapsed()<<"  "<<matcov_flops*1.e-9/t_matcov.elapsed()
                <<"\nreconstructor  "<<t_R.elapsed()     <<"  "<<reconstructor_flops*1.e-9/t_R.elapsed()
                <<"\nCee_Cvv        "<<t_CeeCvv.elapsed()<<"  "<<CeeCvv_flops*1.e-9/t_CeeCvv.elapsed()
                <<"\nintersample    "<<t_psf.elapsed()   <<"  "<<psf_flops*1.e-9/t_psf.elapsed()
                <<"\nrefine         "<<t_refine.elapsed()
                <<"\nobservation    "<<t_obs.elapsed()
                <<"\ncomputation    "<<compute_time      <<"  "<<flops*1.e-9/compute_time   
                <<std::endl;
        }
        std::cout<<"total          "<<t_glob.elapsed()  <<"  "<<flops*1.e-9/t_glob.elapsed()
                <<std::endl;



    MORSE_Finalize();
    return 0;
}


/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/


#include "moao_lapack.hpp"



std::list<Arg> args={Arg("sys_path"     ,"path to the system parameter file"            ,""),
                     Arg("atm_path"     ,"path to the atmospheric parameter file"       ,""),
                     Arg("nmeas"        ,"number of sensor measurements"                ,"0" ),
                     Arg("nmeasTS"      ,"number of measurements for the truth sensor"  ,"0" ),
                     Arg("nact"         ,"number of actuators"                          ,"0" ),
                     Arg("maxobs"       ,"number of observation"                        ,"1" ),
                     Arg("maxrefine"    ,"number of refinment per observation"          ,"1" ),
                     Arg("snapPerNight" ,"number of snapshot per night"                 ,"8" ),
                     Arg("verbose"      ,""                                             ,"0" )
                     };
CmdLine cmd(args);

using real_t=double;

int main(int argc, char *argv[]){



    cmd.parse(argc,argv);
    std::string sys_path=cmd.getString("sys_path");
    std::string atm_path=cmd.getString("atm_path");
    int nmeas           = cmd.getInt("nmeas");
    int nmeasts         = cmd.getInt("nmeasTS");
    int nact            = cmd.getInt("nact");
    int maxobs          = cmd.getInt("maxobs");
    int maxrefine       = cmd.getInt("maxrefine");
    int snapPerNight    = cmd.getInt("snapPerNight");
    int verbose         = cmd.getInt("verbose");
    cmd.recap(true);

    int obsIdx=-1;
    int nightIdx=0;
    int snapIdx=0;

    //TODO double matcov_flops=0., reconstructor_flops=0., CeeCvv_flops=0., psf_flops=0.;
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
    nact=naxes[0];
    intersample_prepare(&isample, nact*nact, tomo.Nssp[tomo.sys.nW-1]+1, tomo.sys.diam, sys_path.c_str());
    intersample_init(&isample);
#endif


    nmeas   =tomo.getNMeas();
    nmeasts =tomo.getNMeasTS();
    std::cout<<"using "<<nmeas<<" meas and "<<nmeasts<<" ts meas"<<std::endl;
    std::cout<<"Using "<<nact<<" actuators"<<std::endl;
    t_preprocess.stop();

    //allocate matrices
    t_alloc.start();
    real_t *Cmm = (real_t*)calloc(nmeas  *nmeas                , sizeof(real_t));
    real_t *R   = (real_t*)calloc(nmeasts*nmeas*tomo.sys.nTarget, sizeof(real_t));
    real_t *Ctm = (real_t*)calloc(nmeasts*nmeas                , sizeof(real_t));
    real_t *Ctt = (real_t*)calloc(nmeasts*nmeasts              , sizeof(real_t));
    real_t *Cee = (real_t*)calloc(nmeasts*nmeasts              , sizeof(real_t));
    real_t *Cvv = (real_t*)calloc(nact   *nact                 , sizeof(real_t));
    real_t *Dx  = (real_t*)calloc(nact   *nmeasts              , sizeof(real_t));
    real_t *Tmp = (real_t*)calloc(nact   *nmeasts              , sizeof(real_t));
    real_t *Rg;
    t_alloc.stop();

#ifdef USE_INTERSAMPLE
    t_preprocess.start();
    //read Dx from file
    readFits(DxFile.c_str(), Dx);
    t_preprocess.stop();
#endif

    for(int nbrefine=0; nbrefine<maxrefine; nbrefine++){
        if(maxrefine>0){
            t_refine.start();
            if(verbose) std::cout<<"\nnbrefine"<< nbrefine+1<<std::endl;;
            nightIdx=0;
            snapIdx=nbrefine;
            obsIdx=-1;
            atmFileName ="/prof"+std::to_string(snapIdx * snapPerNight + obsIdx)+"-atmos-night"+std::to_string(nightIdx)+".txt";
            tomo.updateAtm(atm_path+atmFileName);

            //----------------------
            //generating Cmm
            t_matcov.start();
            MOAO_COMMON::matcov_tile(Cmm, nmeas, nmeas, 0, 0, nmeas, &tomo, 1);
            t_matcov.stop();

            //----------------------
            //generating Ctm(s)
	        for(int nbgal=0; nbgal<tomo.sys.nTarget; nbgal++){

                if(verbose)std::cout<<"target "<<nbgal<<std::endl;
                Rg = &R[nmeas*nmeasts*nbgal];
                tomo.setTarget(nbgal);

                t_matcov.start();
                MOAO_COMMON::matcov_tile(Rg, nmeasts, nmeas, 0, 0, nmeasts, &tomo, 3);
                t_matcov.stop();
            }

            //----------------------
            //computing reconstructor
            t_R.start();
            MOAO_LAPACK::reconstructor(nmeas, nmeasts, Cmm, nmeas, R, nmeasts, tomo.sys.nTarget);
            t_R.stop();
            t_refine.stop();
        }

        t_obs.start();
        for(int nbobs=0; nbobs<maxobs; nbobs++){
            if(verbose) std::cout<<" nbobs "<<nbobs+1<<std::endl;
            //need to get new atm-params
            nightIdx=0;
            snapIdx=nbrefine;
            obsIdx=nbobs;
            //update atm params
            atmFileName ="/prof"+std::to_string(snapIdx * snapPerNight + obsIdx)+"-atmos-night"+std::to_string(nightIdx)+".txt";
            tomo.updateAtm(atm_path+atmFileName);

            //----------------------
            //generating Cmm (destroyed by reconstructor computation)
            t_matcov.start();
            MOAO_COMMON::matcov_tile(Cmm, nmeas, nmeas, 0, 0, nmeas, &tomo, 1);
            t_matcov.stop();

            for(int nbgal=0; nbgal<tomo.sys.nTarget; nbgal++){
                if(verbose) std::cout<<"  target "<<nbgal<<std::endl;
                Rg = &R[nmeas*nmeasts*nbgal];
                tomo.setTarget(nbgal);
                //----------------------
                //generating Ctm (destroyed by reconstructor computation)
                t_matcov.start();
                MOAO_COMMON::matcov_tile(Ctm, nmeasts, nmeas, 0, 0, nmeasts, &tomo, 3);
                t_matcov.stop();

                //----------------------
                //generating Ctt
                t_matcov.start();
                MOAO_COMMON::matcov_tile(Ctt, nmeasts, nmeasts, 0, 0, nmeasts, &tomo, 4);
                t_matcov.stop();

                //----------------------
                //generating error and voltage
                t_CeeCvv.start();
                MOAO_LAPACK::Cee_Cvv(nmeas, nmeasts, nact, Cmm, nmeas, Ctt, nmeasts, Ctm, nmeasts, Rg, nmeasts, Dx, nact, Cee, nmeasts, Cvv, nact, Tmp, nact);
                t_CeeCvv.stop();

#ifdef USE_INTERSAMPLE
                //----------------------
                //output file for psf
                char isample_output_filename[256];
                //compute psf
                t_psf.start();
                intersample_process(&isample, Cvv);
                t_psf.stop();

                psfFileName="psf_oloop"+std::to_string(nbrefine)+"_iloop"+std::to_string(nbobs)+"_gal"+std::to_string(nbgal)+".fits";
                writeFits(psfFileName.c_str(), isample.N, isample.N, isample.dphi);

#endif
            }//End of galaxies loop
            if(verbose) std::cout<<"-----end obs\n"<<std::endl;
        }//End of OBS loop
        t_obs.stop();
        if(verbose) std::cout<<"=====end refine\n"<<std::endl;
    }// End of REFINE loop
    t_glob.stop();

    std::cout<<"========================================"
            <<"\nsimu characteristics"
            <<"\nncores   : "<<1
            <<"\nngpus    : "<<0
            <<"\ntile size: "<<"n/a"   
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
        //double flops=matcov_flops+reconstructor_flops+CeeCvv_flops+psf_flops;
        std::cout<<"Performances:  time (s)   GFlop/s"
                <<"\nallocation     "<<t_alloc.elapsed()
                <<"\npreprocess     "<<t_preprocess.elapsed()
                <<"\nmatcov         "<<t_matcov.elapsed()<<"  "//<<matcov_flops*1.e-9/t_matcov.elapsed()
                <<"\nreconstructor  "<<t_R.elapsed()     <<"  "//<<reconstructor_flops*1.e-9/t_R.elapsed()
                <<"\nCee_Cvv        "<<t_CeeCvv.elapsed()<<"  "//<<CeeCvv_flops*1.e-9/t_CeeCvv.elapsed()
                <<"\nintersample    "<<t_psf.elapsed()   <<"  "//<<psf_flops*1.e-9/t_psf.elapsed()
                <<"\nrefine         "<<t_refine.elapsed()
                <<"\nobservation    "<<t_obs.elapsed()
                <<"\ncomputation    "<<compute_time      <<"  "//<<flops*1.e-9/compute_time   
                <<"\ntotal          "<<t_glob.elapsed()  <<"  "//<<flops*1.e-9/t_glob.elapsed()
                <<std::endl;


    free(Cmm);
    free(Ctm);
    free(Ctt);
    free(R);
    free(Cee);
    free(Cvv);
    free(Dx);
    free(Tmp);

    return 0;
}

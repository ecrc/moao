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
#include <algorithm>

#include "moao_chameleon.hpp"

#include "flops.h"

std::list<Arg> args={Arg("sys_path"     ,"path to the system parameter file"                ,""),
                     Arg("sys_check"    ,"path to the system parameter file"                ,""),
                     Arg("atm_path"     ,"path to the atmospheric parameter file"           ,""),
                     Arg("out_path"     ,"output folder for the generated matrices"         ,""),
                     Arg("compare_all"  ,"compare all matrices if true, only psf otherwise" ,"0"),
                     Arg("write_all"    ,"write all matrices if true, only psf otherwise"   ,"0"),
                     Arg("accuracy"     ,"expected accuracy   "                             ,"1.e-7"),
                     Arg("ncores"       ,"number of cpu cores for computation"              ,"1"),
                     Arg("ngpus"        ,"number of gpu devices for computation"            ,"0"),
                     Arg("ts"           ,"tile size"                                        ,"200")

                     };
CmdLine cmd(args);

using real_t=double;

int main(int argc, char *argv[]){

    //get command line args
    cmd.parse(argc,argv);
    std::string sys_path=cmd.getString("sys_path");
    std::string sys_check=cmd.getString("sys_check");
    std::string atm_path=cmd.getString("atm_path");
    std::string out_path=cmd.getString("out_path");
    int compareAll = cmd.getInt("compare_all");
    int writeAll = cmd.getInt("write_all");
    double accuracy = cmd.getReal("accuracy");
    int ncores = cmd.getInt("ncores");
    int ngpus = cmd.getInt("ngpus");
    int ts = cmd.getInt("ts");

    cmd.recap(true);

    Timer t_Cmm, t_Ctt, t_Ctm, t_R, t_CeeCvv, t_psf, t_glob;
    real_t err_fro_Cmm=0., err_fro_Ctt=0., err_fro_Ctm=0., err_fro_R=0., err_fro_Cee=0., err_fro_Cvv=0.;
    real_t err_max_Cmm=0., err_max_Ctt=0., err_max_Ctm=0., err_max_R=0., err_max_Cee=0., err_max_Cvv=0.;
    double  err_fro_psf=0.,err_max_psf=0.;
    double flops_Cmm=0., flops_Ctm=0., flops_Ctt=0.,flops_R=0., flops_CeeCvv=0.,flops_psf=0.;
    std::string fname;

    t_glob.start();
    //init tomo_struct
    std::string atmFileName = "/prof-1-atmos-night0.txt";
    Tomo_struct<real_t> tomo(sys_path+"/sys-params.txt",atm_path+atmFileName);
    long int nmeas   =tomo.getNMeas();
    long int nmeasts =tomo.getNMeasTS();

    std::cout<<"Working on "<<ncores<<" cores and "<<ngpus<<" gpus"<<std::endl;
    std::cout<<"Tile size:"<<ts<<std::endl;
    std::cout<<"using "<<nmeas<<" meas and "<<nmeasts<<" ts meas"<<std::endl;

#ifdef USE_INTERSAMPLE
    //init isample if needed
    //get number of actuators: define by input matrix Dx
    struct isample_struct isample;
    long naxes[2];
    std::string DxFile=sys_path+"/Dx.fits";
    getFitsDims(DxFile.c_str(),naxes);
    int nact=naxes[0];
    std::cout<<"Using "<<nact<<" actuators"<<std::endl;
    intersample_prepare(&isample, nact*nact, tomo.Nssp[tomo.sys.nW-1]+1, tomo.sys.diam, sys_path.c_str());
    intersample_init(&isample);
#endif

    MORSE_Init(ncores,ngpus);


    real_t *Cmm = (real_t*)calloc(nmeas  *nmeas     , sizeof(real_t));
    real_t *R   = (real_t*)calloc(nmeasts*nmeas     , sizeof(real_t));
    real_t *Ctm = (real_t*)calloc(nmeasts*nmeas     , sizeof(real_t));
    real_t *Ctt = (real_t*)calloc(nmeasts*nmeasts   , sizeof(real_t));
    real_t *Cee = (real_t*)calloc(nmeasts*nmeasts   , sizeof(real_t));
    real_t *Cvv = (real_t*)calloc(nact   *nact      , sizeof(real_t));
    real_t *Dx  = (real_t*)calloc(nact   *nmeasts   , sizeof(real_t));
    real_t *Tmp = (real_t*)calloc(nmeas  *nmeas     , sizeof(real_t));
    real_t *REF = (real_t*)calloc(nmeas  *nmeas     , sizeof(real_t));

    MORSE_desc_t *descCmm = NULL, *descCtm = NULL, *descCtt = NULL;
    MORSE_desc_t *descR = NULL, *descCee = NULL, *descCvv = NULL, *descDx = NULL, *descTmp = NULL;

    MORSE_sequence_t *sequence;
    MORSE_request_t request[5]={MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER,
                                MORSE_REQUEST_INITIALIZER,MORSE_REQUEST_INITIALIZER};
    MORSE_Sequence_Create(&sequence);

    MORSE_Desc_Create(&descCmm, Cmm, MORSE<real_t>::real, ts, ts, ts*ts, nmeas  , nmeas  , 0, 0, nmeas  , nmeas  ,1,1);
    MORSE_Desc_Create(&descCtm, Ctm, MORSE<real_t>::real, ts, ts, ts*ts, nmeasts, nmeas  , 0, 0, nmeasts, nmeas  ,1,1);
    MORSE_Desc_Create(&descCtt, Ctt, MORSE<real_t>::real, ts, ts, ts*ts, nmeasts, nmeasts, 0, 0, nmeasts, nmeasts,1,1);
    MORSE_Desc_Create(&descR  , R,   MORSE<real_t>::real, ts, ts, ts*ts, nmeasts, nmeas  , 0, 0, nmeasts, nmeas  ,1,1);
    MORSE_Desc_Create(&descCee, Cee, MORSE<real_t>::real, ts, ts, ts*ts, nmeasts, nmeasts, 0, 0, nmeasts, nmeasts,1,1);
    MORSE_Desc_Create(&descCvv, Cvv, MORSE<real_t>::real, ts, ts, ts*ts, nact   , nact   , 0, 0, nact   , nact   ,1,1);
    MORSE_Desc_Create(&descDx , Dx,  MORSE<real_t>::real, ts, ts, ts*ts, nact   , nmeasts, 0, 0, nact   , nmeasts,1,1);
    MORSE_Desc_Create(&descTmp, Tmp, MORSE<real_t>::real, ts, ts, ts*ts, nact   , nmeasts, 0, 0, nact   , nmeasts,1,1);


    tomo.updateAtm(atm_path+atmFileName);

    //----------------------
    //generating Cmm
    std::cout<<"Cmm     "<<std::flush;
    t_Cmm.start();
    MOAO_CHAMELEON::matcov_tile(descCmm, sequence, &request[0], &tomo, 1);
    MORSE_Sequence_Wait(sequence);
    t_Cmm.stop();
    flops_Cmm+=FLOPS_DCovMat(nmeas, nmeas, tomo.atm.nLayer, 1);
    std::cout<<"done "<<t_Cmm.elapsed()<<std::endl;
    t_glob.stop();
    if(!out_path.empty() && writeAll || compareAll){
        MORSE_Tile_to_Lapack(descCmm,Tmp,nmeas);
    }
    if(!out_path.empty() && writeAll){
        writeFits("Cmm.fits", nmeas, nmeas,Tmp);
    }
    if(compareAll){
        fname=sys_check+"Cmm.fits";
        readFits(fname.c_str(),REF);
        compareMatrices(Tmp,REF, nmeas  , nmeas  ,'N','T',err_max_Cmm,err_fro_Cmm);
    }
    t_glob.start();

    //----------------------
    //generating Ctm
    std::cout<<"Ctm     "<<std::flush;
    t_Ctm.start();
    MOAO_CHAMELEON::matcov_tile(descR, sequence, &request[1], &tomo, 3);
    MORSE_Sequence_Wait(sequence);
    t_Ctm.stop();
    flops_Ctm+=FLOPS_DCovMat(nmeasts, nmeas, tomo.atm.nLayer, 3);
    std::cout<<"done "<<t_Ctm.elapsed()<<std::endl;
    t_glob.stop();
    if(!out_path.empty() && writeAll || compareAll){
        MORSE_Tile_to_Lapack(descR,Tmp,nmeasts);
    }
    if(!out_path.empty() && writeAll){
        writeFits("Ctm.fits", nmeasts, nmeas,Tmp);
    }
    if(compareAll){
        fname=sys_check+"Ctm.fits";
        readFits(fname.c_str(),REF);
        compareMatrices(Tmp,REF, nmeasts, nmeas  ,'N','T',err_max_Ctm,err_fro_Ctm);
    }
    t_glob.start();

    //----------------------
    //computing reconstructor
    std::cout<<"R       "<<std::flush;
    t_R.start();
    MOAO_CHAMELEON::reconstructor<real_t>(descCmm,descR,sequence,0);
    MORSE_Sequence_Wait(sequence);
    t_R.stop();
    flops_R+=FLOPS_reconstructor(nmeas,nmeasts);
    std::cout<<"done "<<t_R.elapsed()<<std::endl;
    t_glob.stop();
    if(!out_path.empty() && writeAll || compareAll){
        MORSE_Tile_to_Lapack(descR,Tmp,nmeasts);
    }
    if(!out_path.empty() && writeAll){
        writeFits("R.fits", nmeasts, nmeas,Tmp);
    }
    if(compareAll){
        fname=sys_check+"R.fits";
        readFits(fname.c_str(),REF);
        compareMatrices(Tmp,REF, nmeasts, nmeas  ,'N','T',err_max_R  ,err_fro_R  );
    }
    t_glob.start();

    //----------------------
    //generating Cmm (destroyed by reconstructor computation)
    std::cout<<"Cmm     "<<std::flush;
    t_Cmm.start();
    MOAO_CHAMELEON::matcov_tile(descCmm, sequence, &request[2], &tomo, 1);
    MORSE_Sequence_Wait(sequence);
    t_Cmm.stop();
    flops_Cmm+=FLOPS_DCovMat(nmeas, nmeas, tomo.atm.nLayer, 1);
    std::cout<<"done"<<std::endl;
    
    //----------------------
    //generating Ctm (destroyed by reconstructor computation)
    std::cout<<"Ctm     "<<std::flush;
    t_Ctm.start();
    MOAO_CHAMELEON::matcov_tile(descCtm, sequence, &request[3], &tomo, 3);
    MORSE_Sequence_Wait(sequence);
    t_Ctm.stop();
    flops_Ctm+=FLOPS_DCovMat(nmeasts, nmeas, tomo.atm.nLayer, 3);
    std::cout<<"done"<<std::endl;

    //----------------------
    //generating Ctt
    std::cout<<"Ctt     "<<std::flush;
    t_Ctt.start();
    MOAO_CHAMELEON::matcov_tile(descCtt, sequence, &request[4], &tomo, 4);
    MORSE_Sequence_Wait(sequence);
    t_Ctt.stop();
    flops_Ctt+=FLOPS_DCovMat(nmeasts, nmeasts, tomo.atm.nLayer, 4);
    std::cout<<"done"<<std::endl;
    t_glob.stop();
    if(!out_path.empty() && writeAll || compareAll){
        MORSE_Tile_to_Lapack(descCtt,Tmp,nmeasts);
    }
    if(!out_path.empty() && writeAll){
        writeFits("Ctt.fits", nmeasts, nmeasts,Tmp);
    }
    if(compareAll){
        fname=sys_check+"Ctt.fits";
        readFits(fname.c_str(),REF);
        compareMatrices(Tmp,REF, nmeasts, nmeasts,'N','T',err_max_Ctt,err_fro_Ctt);
    }
    t_glob.start();

#ifdef USE_INTERSAMPLE
    //----------------------
    //read Dx from file
    fname=sys_path+"Dx.fits";
    readFits(fname.c_str(), Tmp);
    MORSE_Lapack_to_Tile(Tmp,nact,descDx);

    //----------------------
    //computing error and voltage
    std::cout<<"Cee,Cvv "<<std::flush;
    t_CeeCvv.start();
    MOAO_CHAMELEON::Cee_Cvv<real_t>(descCmm, descCtt, descCtm, descR, descDx,
                          descCee, descCvv, descTmp, sequence);
    MORSE_Sequence_Wait(sequence);
    t_CeeCvv.stop();
    flops_CeeCvv+=FLOPS_CeeCvv(nmeas,nmeasts,nact);
    std::cout<<"done"<<std::endl;
    t_glob.stop();
    if(!out_path.empty() && writeAll || compareAll){
        MORSE_Tile_to_Lapack(descCee,Tmp,nmeasts);
    }
    if(!out_path.empty() && writeAll){
        writeFits("Cee.fits", nmeasts, nmeasts,Tmp);
    }
    if(compareAll){
        fname=sys_check+"Cee.fits";
        readFits(fname.c_str(),REF);
        compareMatrices(Tmp,REF, nmeasts, nmeasts,'N','T',err_max_Cee,err_fro_Cee);
    }

    if(!out_path.empty() && writeAll || compareAll){
        MORSE_Tile_to_Lapack(descCvv,Tmp,nact);
    }
    if(!out_path.empty() && writeAll){
        writeFits("Cvv.fits", nact, nact,Tmp);
    }
    if(compareAll){
        fname=sys_check+"Cvv.fits";
        readFits(fname.c_str(),REF);
        compareMatrices(Tmp,REF, nact   , nact   ,'N','T',err_max_Cvv,err_fro_Cvv);
    }
    t_glob.start();

    //----------------------
    //compute psf
    std::cout<<"psf     "<<std::flush;
    t_psf.start();
    MORSE_Tile_to_Lapack(descCvv,Tmp,nact);
    intersample_process(&isample, Tmp);
    t_psf.stop();
    flops_psf+=FLOPS_Intersample(nact,nact);
    std::cout<<"done"<<std::endl;
    t_glob.stop();
    if(!out_path.empty() ){
        writeFits("psf.fits", isample.N, isample.N, isample.dphi);
    }
    fname=sys_check+"psf.fits";
    readFits(fname.c_str(),REF);
    compareMatrices(isample.dphi,REF, isample.N, isample.N, 'N','T',err_max_psf,err_fro_psf);
#endif

    //converting nb of flop to GFlops
    double flops_glob=flops_Cmm+flops_Ctm+flops_Ctt+flops_R+flops_CeeCvv+flops_psf;
    flops_Cmm= t_Cmm.elapsed()>0? flops_Cmm/t_Cmm.elapsed()*1.e-9: 0;
    flops_Ctm= t_Ctm.elapsed()>0? flops_Ctm/t_Ctm.elapsed()*1.e-9: 0;
    flops_Ctt= t_Ctt.elapsed()>0? flops_Ctt/t_Ctt.elapsed()*1.e-9: 0;
    flops_R  = t_R  .elapsed()>0? flops_R  /t_R  .elapsed()*1.e-9: 0;
    flops_CeeCvv= t_CeeCvv.elapsed()>0? flops_CeeCvv/t_CeeCvv.elapsed()*1.e-9: 0;
    flops_psf= t_psf.elapsed()>0? flops_psf/t_psf.elapsed()*1.e-9: 0;
    flops_glob= t_psf.elapsed()>0? flops_psf/t_psf.elapsed()*1.e-9: 0;
    
    if(compareAll<1){
    std::cout<<"WARNING: error only for psf!"<<std::endl;
    }
    //TODO 
    std::cout<<"=============================="<<std::endl;
    fprintf(stdout,"operation: Cmm           Ctm           Ctt           R             Cee           Cvv           psf          total\n");
    fprintf(stdout,"time (s) : %e  %e  %e  %e  %e                %e %e\n", t_Cmm.elapsed(),t_Ctm.elapsed(),t_Ctt.elapsed(),t_R.elapsed(),t_CeeCvv.elapsed(),t_psf.elapsed(),t_glob.elapsed());
    fprintf(stdout,"GFlops   : %e  %e  %e  %e  %e                %e %e\n", flops_Cmm,flops_Ctm, flops_Ctt,flops_R, flops_CeeCvv, flops_psf, flops_glob);
    
    fprintf(stdout,"F norm   : %e  %e  %e  %e  %e  %e  %e    \n", err_fro_Cmm, err_fro_Ctm, err_fro_Ctt,err_fro_R ,err_fro_Cee,err_fro_Cvv,err_fro_psf);
    fprintf(stdout,"max err  : %e  %e  %e  %e  %e  %e  %e    \n", err_max_Cmm, err_max_Ctm, err_max_Ctt,err_max_R ,err_max_Cee,err_max_Cvv,err_max_psf);
    
    std::cout<<"=============================="<<std::endl;

    free(Cmm);
    free(Ctm);
    free(Ctt);
    free(R);
    free(Cee);
    free(Cvv);
    free(Dx);
    free(Tmp);

if(err_max_Cmm > accuracy ||
    err_max_Ctm > accuracy ||
    err_max_Ctt > accuracy ||
    err_max_R   > accuracy ||
    err_max_Cee > accuracy ||
    err_max_Cvv > accuracy ||
    err_max_psf > accuracy 
){
std::cerr<<"test failed: bad accuracy\n\t expected "<<accuracy<<std::endl;
    MORSE_Sequence_Wait(sequence);
    MORSE_Sequence_Destroy(sequence);
    MORSE_Finalize();
    return 1;
}

    return 0;
}

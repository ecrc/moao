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

#include "moao_lapack.hpp"

//#include "flops.h"

std::list<Arg> args={Arg("sys_path"     ,"path to the system parameter file"                ,""),
                     Arg("sys_check"    ,"path to the system parameter file"                ,""),
                     Arg("atm_path"     ,"path to the atmospheric parameter file"           ,""),
                     Arg("out_path"     ,"output folder for the generated matrices"         ,""),
                     Arg("compare_all"  ,"compare all matrices if true, only psf otherwise" ,"0"),
                     Arg("write_all"    ,"write all matrices if true, only psf otherwise"   ,"0"),
                     Arg("accuracy"     ,"expected accuracy   "                             ,"1.e-7")

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

    cmd.recap(true);

    Timer t_Cmm, t_Ctt, t_Ctm, t_R, t_CeeCvv, t_psf, t_glob;
    real_t err_fro_Cmm=0., err_fro_Ctt=0., err_fro_Ctm=0., err_fro_R=0., err_fro_Cee=0., err_fro_Cvv=0.;
    real_t err_max_Cmm=0., err_max_Ctt=0., err_max_Ctm=0., err_max_R=0., err_max_Cee=0., err_max_Cvv=0.;
    double  err_fro_psf=0.,err_max_psf=0.;
    //double flops_Cmm=0., flops_Ctm=0., flops_Ctt=0.,flops_R=0., flops_CeeCvv=0.,flops_psf=0.;
    std::string fname;

    t_glob.start();
    //init tomo_struct
    std::string atmFileName = "/prof-1-atmos-night0.txt";
    Tomo_struct<real_t> tomo(sys_path+"/sys-params.txt",atm_path+atmFileName);
    int nmeas   =tomo.getNMeas();
    int nmeasts =tomo.getNMeasTS();
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


    real_t *Cmm = (real_t*)calloc(nmeas  *nmeas     , sizeof(real_t));
    real_t *R   = (real_t*)calloc(nmeasts*nmeas     , sizeof(real_t));
    real_t *Ctm = (real_t*)calloc(nmeasts*nmeas     , sizeof(real_t));
    real_t *Ctt = (real_t*)calloc(nmeasts*nmeasts   , sizeof(real_t));
    real_t *Cee = (real_t*)calloc(nmeasts*nmeasts   , sizeof(real_t));
    real_t *Cvv = (real_t*)calloc(nact   *nact      , sizeof(real_t));
    real_t *Dx  = (real_t*)calloc(nact   *nmeasts   , sizeof(real_t));
    real_t *Tmp = (real_t*)calloc(nact   *nmeasts   , sizeof(real_t));
    real_t *REF = (real_t*)calloc(nmeas  *nmeas     , sizeof(real_t));

    tomo.updateAtm(atm_path+atmFileName);

    //----------------------
    //generating Cmm
    std::cout<<"Cmm     "<<std::flush;
    t_Cmm.start();
    MOAO_COMMON::matcov_tile(Cmm, nmeas, nmeas, 0, 0, nmeas, &tomo, 1);
    t_Cmm.stop();
    //flops_Cmm+=FLOPS_DCovMat(nmeas, nmeas, tomo.atm.nLayer, 1);
    std::cout<<"done"<<std::endl;
    t_glob.stop();
    if(!out_path.empty() && writeAll){
        writeFits("Cmm.fits", nmeas, nmeas,Cmm);
    }
    if(compareAll){
        fname=sys_check+"Cmm.fits";
        readFits(fname.c_str(),REF);
        compareMatrices(Cmm,REF, nmeas  , nmeas  ,'N','T',err_max_Cmm,err_fro_Cmm);
    }
    t_glob.start();

    //----------------------
    //generating Ctm
    std::cout<<"Ctm     "<<std::flush;
    t_Ctm.start();
    MOAO_COMMON::matcov_tile(R, nmeasts, nmeas, 0, 0, nmeasts, &tomo, 3);
    t_Ctm.stop();
    //flops_Ctm+=FLOPS_DCovMat(nmeasts, nmeas, tomo.atm.nLayer, 3);
    std::cout<<"done"<<std::endl;
    t_glob.stop();
    if(!out_path.empty() && writeAll){
        writeFits("Ctm.fits", nmeasts, nmeas,R);
    }
    if(compareAll){
        fname=sys_check+"Ctm.fits";
        readFits(fname.c_str(),REF);
        compareMatrices(R,REF, nmeasts, nmeas  ,'N','T',err_max_Ctm,err_fro_Ctm);
    }
    t_glob.start();

    //----------------------
    //computing reconstructor
    std::cout<<"R       "<<std::flush;
    t_R.start();
    MOAO_LAPACK::reconstructor(nmeas, nmeasts, Cmm, nmeas, R, nmeasts, 1);
    t_R.stop();
    //flops_R+=FLOPS_reconstructor(nmeasts,nmeas);
    std::cout<<"done"<<std::endl;
    t_glob.stop();
    if(!out_path.empty() && writeAll){
        writeFits("R.fits", nmeasts, nmeas,R);
    }
    if(compareAll){
        fname=sys_check+"R.fits";
        readFits(fname.c_str(),REF);
        compareMatrices(R,REF, nmeasts, nmeas  ,'N','T',err_max_R  ,err_fro_R  );
    }
    t_glob.start();

    //----------------------
    //generating Cmm (destroyed by reconstructor computation)
    std::cout<<"Cmm     "<<std::flush;
    t_Cmm.start();
    MOAO_COMMON::matcov_tile(Cmm, nmeas, nmeas, 0, 0, nmeas, &tomo, 1);
    t_Cmm.stop();
    //flops_Cmm+=FLOPS_DCovMat(nmeas, nmeas, tomo.atm.nLayer, 1);
    std::cout<<"done"<<std::endl;
    
    //----------------------
    //generating Ctm (destroyed by reconstructor computation)
    std::cout<<"Ctm     "<<std::flush;
    t_Ctm.start();
    MOAO_COMMON::matcov_tile(Ctm, nmeasts, nmeas, 0, 0, nmeasts, &tomo, 3);
    t_Ctm.stop();
    //flops_Ctm+=FLOPS_DCovMat(nmeasts, nmeas, tomo.atm.nLayer, 3);
    std::cout<<"done"<<std::endl;

    //----------------------
    //generating Ctt
    std::cout<<"Ctt     "<<std::flush;
    t_Ctt.start();
    MOAO_COMMON::matcov_tile(Ctt, nmeasts, nmeasts, 0, 0, nmeasts, &tomo, 4);
    t_Ctt.stop();
    //flops_Ctt+=FLOPS_DCovMat(nmeasts, nmeasts, tomo.atm.nLayer, 4);
    std::cout<<"done"<<std::endl;
    t_glob.stop();
    if(!out_path.empty() && writeAll){
        writeFits("Ctt.fits", nmeasts, nmeasts,Ctt);
    }
    if(compareAll){
        fname=sys_check+"Ctt.fits";
        readFits(fname.c_str(),REF);
        compareMatrices(Ctt,REF, nmeasts, nmeasts,'N','T',err_max_Ctt,err_fro_Ctt);
    }
    t_glob.start();

#ifdef USE_INTERSAMPLE
    //----------------------
    //read Dx from file
    fname=sys_path+"Dx.fits";
    readFits(fname.c_str(), Dx);

    //----------------------
    //computing error and voltage
    std::cout<<"Cee,Cvv "<<std::flush;
    t_CeeCvv.start();
    MOAO_LAPACK::Cee_Cvv(nmeas, nmeasts, nact, Cmm, nmeas, Ctt, nmeasts, Ctm, nmeasts, R, nmeasts, Dx, nact, Cee, nmeasts, Cvv, nact, Tmp, nact);
    t_CeeCvv.stop();
    //flops_CeeCvv+=FLOPS_CeeCvv(nmeas,nmeasts,nact);
    std::cout<<"done"<<std::endl;
    t_glob.stop();
    if(!out_path.empty() && writeAll){
        writeFits("Cee.fits", nmeasts, nmeasts,Cee);
    }
    if(compareAll){
        fname=sys_check+"Cee.fits";
        readFits(fname.c_str(),REF);
        compareMatrices(Cee,REF, nmeasts, nmeasts,'N','T',err_max_Cee,err_fro_Cee);
    }

    if(!out_path.empty() && writeAll){
        writeFits("Cvv.fits", nact, nact,Cvv);
    }
    if(compareAll){
        fname=sys_check+"Cvv.fits";
        readFits(fname.c_str(),REF);
        compareMatrices(Cvv,REF, nact   , nact   ,'N','T',err_max_Cvv,err_fro_Cvv);
    }
    t_glob.start();

    //----------------------
    //compute psf
    std::cout<<"psf     "<<std::flush;
    t_psf.start();
    intersample_process(&isample, Cvv);
    t_psf.stop();
    //flops_psf+=FLOPS_Intersample(nact,nact);
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
    //double flops_glob=flops_Cmm+flops_Ctm+flops_Ctt+flops_R+flops_CeeCvv+flops_psf;
    //flops_Cmm= t_Cmm.elapsed()>0? flops_Cmm/t_Cmm.elapsed()*1.e-9: 0;
    //flops_Ctm= t_Ctm.elapsed()>0? flops_Ctm/t_Ctm.elapsed()*1.e-9: 0;
    //flops_Ctt= t_Ctt.elapsed()>0? flops_Ctt/t_Ctt.elapsed()*1.e-9: 0;
    //flops_R  = t_R  .elapsed()>0? flops_R  /t_R  .elapsed()*1.e-9: 0;
    //flops_CeeCvv= t_CeeCvv.elapsed()>0? flops_CeeCvv/t_CeeCvv.elapsed()*1.e-9: 0;
    //flops_psf= t_psf.elapsed()>0? flops_psf/t_psf.elapsed()*1.e-9: 0;
    //flops_glob= t_psf.elapsed()>0? flops_psf/t_psf.elapsed()*1.e-9: 0;
    
    if(compareAll<1){
    std::cout<<"WARNING: error only for psf!"<<std::endl;
    }
    //TODO 
    std::cout<<"=============================="<<std::endl;
    fprintf(stdout,"operation: Cmm           Ctm           Ctt           R             Cee           Cvv           psf          total\n");
    fprintf(stdout,"time (s) : %e  %e  %e  %e  %e                %e %e\n", t_Cmm.elapsed(),t_Ctm.elapsed(),t_Ctt.elapsed(),t_R.elapsed(),t_CeeCvv.elapsed(),t_psf.elapsed(),t_glob.elapsed());
    //fprintf(stdout,"GFlops   : %e  %e  %e  %e  %e                %e %e\n", flops_Cmm,flops_Ctm, flops_Ctt,flops_R, flops_CeeCvv, flops_psf, flops_glob);
    
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
    return 1;
}

    return 0;
}

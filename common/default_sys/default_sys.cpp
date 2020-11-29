
#include <string>
#include <list>
#include "common.hpp"

#define real_t  double


std::list<Arg> args={Arg("sys_file" ,"output file name"                             ,"sys-params.txt"),
                     Arg("diameter" ,"telescope diameter"                           ,"38"            ),
                     Arg("seed"     ,"seed for random generation of the guide star" ,"1234"          ),
                     Arg("nngs"     ,"number of NGS: Natural Guide Stars"           ,"4"             ),
                     Arg("nlgs"     ,"number of LGS: Laser Guide Stars"             ,"6"             ),
                     Arg("mag"      ,"magnitude of the NGS"                         ,"13"            ),
                     Arg("flux"     ,"LGS photon return at M1"                      ,"7.e6"          )
};
CmdLine cmd(args);

int main(int argc, char ** argv){
    cmd.parse(argc,argv);
    cmd.recap();

    std::string  fileName=cmd.getString("sys_file");
    int     seed    =cmd.getInt("seed");
    real_t  diam    =cmd.getReal("diameter");
    int     nNGS    =cmd.getInt("nngs");
    int     nLGS    =cmd.getInt("nlgs");
    real_t  NGSmag  =cmd.getReal("mag");
    real_t  lgsFlux =cmd.getReal("flux");


    SysParams<real_t> def( diam, nNGS, nLGS, NGSmag, lgsFlux, seed);
    def.write(fileName);
}


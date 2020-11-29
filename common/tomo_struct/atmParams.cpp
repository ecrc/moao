/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/
#ifndef ATMPARAMS_IMPL_
#define ATMPARAMS_IMPL_

#include "atmParams.hpp"
#include "utils/utils_vector.hpp"
#include "utils/utils_IO.hpp"
#include <iostream>

/*!\brief Default constructor
 */
template<typename T>
AtmParams<T>::AtmParams(){
    nLayer=0;
}

/*!\brief Constructor
 *
 * Initialize parameters from file
 * @param[in] af : string : parameter file to read
 */
template<typename T>
AtmParams<T>::AtmParams(std::string af){
    nLayer=0;
    read(af);
}

/*!\bief Check the size of the AtmParam vectors
 *
 * each vector of AtmParam must have at least nLayer elements
 * throw an exception if this is not the case
 * return a warning if the vector ahas more elements than the number of layers
 */
template<typename T>
void AtmParams<T>::vectorSize(){
    int err=0;
    err=std::max(err,checkVectorSize<T  >(cn2    ,nLayer     ,"cn2"        ));
    err=std::max(err,checkVectorSize<T  >(h      ,nLayer     ,"h  "        ));
    err=std::max(err,checkVectorSize<T  >(L0     ,nLayer     ,"L0 "        ));
    if(err>0)
        std::cerr<<std::endl;
    if(err>1)
        throw std::runtime_error("System Parameter vector size does not matches");
    
}

/*!\brief Read parameters to file
 *
 * each entry of the parameter file is preceed by a comment line
 *
 *@param[in] af             : string    : input file name 
 *@param[in] checkComment   : bool      : check comment line preceeding data entry
 */
template<typename T>
void AtmParams<T>::read(std::string af,bool checkComment){
    std::cout<<"    Reading atmospheric parameters: "<<af<<std::endl;
    std::string line;
    std::string err("");

    if(!cn2.empty()){cn2.clear();}
    if(!L0.empty()){L0.clear();}
    if(!h.empty()){h.clear();}

    readIO(af,*this);
    if(!err.empty()){
        std::cerr<<"WARNING:\n"<< err<<std::endl;;
    }
    vectorSize();
}

/*!\brief Write parameters to file
 *
 *@param[in] af : string : output file name 
 */
template<typename T>
void AtmParams<T>::write(std::string af){
    std::cout<<"Writing atmospheric parameters: "<<af<<std::endl;
    writeIO(af,*this);
}

/*!\brief  return attribute list
 *
 * @param[inout] e : string : string on which errors will be appened
 * @return : list<IO> : list of IO with all AtmParams attributes
 */
template<typename T>
std::list<IO> AtmParams<T>::getParams(std::string & e)
{
    std::list<IO> l;
     appendNamedIO(l,nLayer , e,"nLayer :       : number of turbulent layers    ");
     appendNamedIO(l,r0     , e,"r0     :       : global fried parameter        ");
     appendNamedIO(l,cn2    , e,"cn2    :       : layer strength r0^(-5/3)      ");
     appendNamedIO(l,h      , e,"h      : meter : layer altitude                ");
     appendNamedIO(l,L0     , e,"L0     : meter : outer scale                   ");
    return l;
}
    
#endif //ATMPARAMS_IMPL_

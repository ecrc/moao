#ifndef ATMPARAMS_H
#define ATMPARAMS_H
#include <vector>
#include <string>
#include "io/IO.hpp"

template<typename T>
class AtmParams{
    public:
    static constexpr const char *name="Atmospheric";
    int             nLayer; //!<       : number of turbulent layers
    T               r0;     //!<       : global fried parameter         
    std::vector<T>  cn2;    //!<       : layer strength r0^(-5/3)      
    std::vector<T>  h;      //!< meter : layer altitude                
    std::vector<T>  L0;     //!< meter : outer scale                   
 
    public:
        AtmParams();
        AtmParams(std::string af);
        void vectorSize();
        void read(std::string af, bool checkComment=false);
        void write(std::string af);

        std::list<IO> getParams(std::string & e);

};
template class AtmParams<float>;
template class AtmParams<double>;


#endif //ATMPARAMS_H

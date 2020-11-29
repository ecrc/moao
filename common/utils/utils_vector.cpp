#ifndef UTILS_VECTOR_IMPL_
#define UTILS_VECTOR_IMPL_

#include "utils_vector.hpp"

#include <iostream>
#include <string.h>
#include <algorithm>
#include <cmath>


/*!\brief check vector size against required value
 *
 * print a message if the size does not matches
 *
 *
 * @param[in] v     : vector<T> : vector to check
 * @param[in] size  : int       : required value
 * @return  0 if size matches; 1 if the vector is too large; 2 if hte vector is too small
 */
template<typename T>
int checkVectorSize(std::vector<T> v, int size, std::string name){
    int err=0;
    if(v.size()>size){
        std::cerr<<"WARNING:\n\tvector "<<name<<" as to many elements: "<<v.size()
            <<".\n\trequire "<<size<<"\n\tTRUNCATING VECTOR.";
        err=1;
    }
    if(v.size()<size){
        std::cerr<<"ERROR:\n\tvector "<<name<<" as only "<<v.size()<<"elements.\n\trequire "<<size<<'\n';
        return 2;
    }
    if(err>0)
        std::cerr<<std::endl;
    return err;
}
template int checkVectorSize<int>(std::vector<int> v, int size, std::string name);
template int checkVectorSize<float>(std::vector<float> v, int size, std::string name);
template int checkVectorSize<double>(std::vector<double> v, int size, std::string name);

/*!\brief generate the coordinates of points on a circle
 *
 * The points are regularly spaced on the circle, the first point is (radius,0)
 *
 * @param[inout]    vX      : vector<T> : X coordinates
 * @param[inout]    vX      : vector<T> : Y coordinates
 * @param[in]       start   : int       : first element of the vectors to set
 * @param[in]       n       : int       : number of points
 * @param[in]       radius  : T         : circle radius
 */
template<typename T>
void fillVectorCircle(std::vector<T> &vX, std::vector<T> &vY,int start, int n, T radius){
    int s=std::min(vX.size(),vY.size());
    if(start+n>s){
        std::cerr<<"error in fillVectorCircle:\n\t out of bounds ("
                 << start+n<<" out of "<<s<<")"<<std::endl;
    }
    const T a=2.*M_PI/n;
    for(int i=start;i<start+n;i++){
        vX[i]=radius*std::cos(i*a);
        vY[i]=radius*std::sin(i*a);
    }
}
template void fillVectorCircle<float>(std::vector<float> &vX, std::vector<float> &vY,int start, int n, float radius);
template void fillVectorCircle<double>(std::vector<double> &vX, std::vector<double> &vY,int start, int n, double radius);

/*!\brief generate random entries in a vector
 *
 * @param[inout]    v       : vector<T> : 
 * @param[in]       start   : int       : first entry to set
 * @param[in]       n       : int       : number of entry to generate
 * @param[in]       min     : T         : minimum generated value
 * @param[in]       range   : T         : range of the generated values
 */
template<typename T>
void fillVectorRand(std::vector<T> &v, int start, int n, T min, T range){
    if(start+n>v.size()){
        std::cerr<<"error in fillVectorRand:\n\t out of bounds ("
                 << start+n<<" out of "<<v.size()<<")"<<std::endl;
    }
    for(auto it=v.begin()+start;it!=v.begin()+start+n;it++){
        *it=rand()/double(RAND_MAX)*range-min;
    }
}
template void fillVectorRand<float>(std::vector<float> &v, int start, int n, float min, float range);
template void fillVectorRand<double>(std::vector<double> &v, int start, int n, double min, double range);

#endif // UTILS_VECTOR_IMPL_

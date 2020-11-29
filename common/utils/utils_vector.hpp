#ifndef UTILs_VECTOR_H
#define UTILs_VECTOR_H

#include <string>
#include <vector>



template<typename T>
int checkVectorSize(std::vector<T> v, int size, std::string name="");

template<typename T>
void fillVectorCircle(std::vector<T> &vX, std::vector<T> &vY,int start, int n, T radius);

template<typename T>
void fillVectorRand(std::vector<T> &v, int start, int n, T min, T range);


#endif // UTILs_VECTOR_H

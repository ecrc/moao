#ifndef NAMED_H
#define NAMED_H

#include <string>
#include <list>
#include "IO.hpp"

template<typename T>
struct Named
{
    T & t_;
    std::string name_;
    std::string & err_;
    Named(T & t, std::string & err, std::string name="");
};

template<typename T>
Named<T> named(T & t, std::string & err, std::string name="");

template<typename T>
std::istream& operator>>(std::istream& in, Named<T> & n);

template<typename T>
std::ostream& operator<<(std::ostream& os, const Named<T> & n);


template<typename T>
IO namedIO(T & t, std::string & err, std::string name="");


#include "named.cpp"

#endif//  NAMED_H

#ifndef UTILS_IO_H
#define UTILS_IO_H

#include <string>

#include "io/IO.hpp"
#include "io/named.hpp"



template<typename T>
//inline
void appendNamedIO(std::list<IO> &l,T &value, std::string &err, std::string comment)
{
    l.push_back(namedIO(value,err,comment));
}


#endif // UTILs_IO_H

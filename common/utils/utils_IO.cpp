#include "utils_IO.hpp"
#include <vector>
template<typename T>
void appendNamedIO(std::list<IO> &l,T &value, std::string &err, std::string comment){
    l.push_back(namedIO(value,err,comment));
}
